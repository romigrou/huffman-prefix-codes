/*
 * Copyright (c) 2024 Romain BAILLY
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy of
 * this software and associated documentation files (the "Software"), to deal in
 * the Software without restriction, including without limitation the rights to
 * use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
 * of the Software, and to permit persons to whom the Software is furnished to do
 * so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

/*
 * This file contains a straightforward implementation of prefix coding using Huffman's algorithm.
 */

#include "utils.h"
#include <algorithm>
#include <cassert>
#include <cstring>
#include <vector>


namespace
{


constexpr size_t   maxSymbolCount = 256;
constexpr unsigned maxCodeLength  = 32;
using Code = std::conditional<maxCodeLength <= 32u, uint32_t, uint64_t>::type;


struct TreeItem
{
    size_t          weight;
    size_t          number; // For leaves, this is the symbol's number
    const TreeItem* children[2] = {};

    friend bool operator>(const TreeItem& a, const TreeItem& b) noexcept
    {
        return (a.weight > b.weight || (a.weight == b.weight && a.number > b.number));
    }
};


const TreeItem* build_huffman_tree(TreeItem tree[2*maxSymbolCount-1], const size_t weights[maxSymbolCount]) noexcept
{
    // Create tree leaves
    size_t usedSymbolCount = 0;
    for (size_t symbol=0; symbol<maxSymbolCount; ++symbol)
        if (weights[symbol] != 0u)
            tree[usedSymbolCount++] = {weights[symbol], symbol};

    // Special case if we have less than two items
    if (usedSymbolCount <= 1u)
        return (usedSymbolCount == 1u) ? tree : nullptr;

    // Make a min heap with the leaves
    using HeapPredicate = std::greater<TreeItem>; // Because the default heap is a max heap and we want a min heap
    std::make_heap(tree, tree+usedSymbolCount, HeapPredicate());
    size_t heapSize = usedSymbolCount;

    // To keep the tree flatter, when items have the same weight those created earlier must be picked first.
    // This is enforced by having all nodes have a number higher than those of leaves and monotonously increasing.
    size_t itemNumber = maxSymbolCount;

    // Create nodes
    size_t itemCount = usedSymbolCount;
    do
    {
        // Pick 1st child (won't be overwritten, can be simply referenced)
        std::pop_heap(tree, tree+heapSize, HeapPredicate());
        const TreeItem& leftChild = tree[--heapSize];

        // Pick 2nd child (must be copied as it will be overwritten by the next node)
        const TreeItem& rightChild = tree[itemCount++] = tree[0];
        std::pop_heap(tree, tree+heapSize, HeapPredicate());

        // Make a node out of them and push it back into the heap
        tree[heapSize-1] = {leftChild.weight+rightChild.weight, itemNumber++, {&leftChild, &rightChild}};
        std::push_heap(tree, tree+heapSize, HeapPredicate());
    }
    while (heapSize != 1u);

    return tree; // Tree root
}


void retrieve_codes_and_lengths(Code codes[maxSymbolCount], uint8_t lengths[maxSymbolCount], const TreeItem* item, Code path=0, unsigned depth=0) noexcept
{
    if (item->children[0] == nullptr) // Leaf
    {
        codes[item->number]   = path;
        lengths[item->number] = uint8_t(depth);
    }
    else                              // Node
    {
        retrieve_codes_and_lengths(codes, lengths, item->children[0], (path<<1)|0, depth+1);
        retrieve_codes_and_lengths(codes, lengths, item->children[1], (path<<1)|1, depth+1);
    }
 }


} // namespace


void encode(std::vector<uint8_t>& encodedData, const uint8_t* data, size_t dataSize)
{
    const Clock::time_point t0 = Clock::now();

    // Compute each symbol's weight
    size_t weights[maxSymbolCount];
    compute_histogram(weights, data, dataSize);
    const Clock::time_point tHisto = Clock::now();

    // Build the Huffman tree
    TreeItem tree[2*maxSymbolCount-1];
    const TreeItem* treeRoot = build_huffman_tree(tree, weights);
    const Clock::time_point tHuffman = Clock::now();

    // Retrieve symbols' codes and lengths
    Code codes[maxSymbolCount] = {};
    uint8_t lengths[maxSymbolCount] = {};
    retrieve_codes_and_lengths(codes, lengths, treeRoot);
    const Clock::time_point tCodes = Clock::now();

    // Compute the size of the encoded data
    size_t encodedDataSizeInBits = 0;
    for (size_t i=0; i<dataSize; ++i)
        encodedDataSizeInBits += lengths[data[i]];
    const size_t headerSize = sizeof(size_t) * maxSymbolCount; // The weight of each symbol
    const size_t encodedDataSize = headerSize + (encodedDataSizeInBits + 7) / 8;
    encodedData.resize(encodedDataSize);

    // Write header
    uint8_t* writePtr = encodedData.data();
    memcpy(writePtr, weights, headerSize);
    writePtr += headerSize;

    // Encode data
    uint8_t  bitBuffer     = 0;
    unsigned bitBufferSize = 0;
    for (size_t i=0; i<dataSize; ++i)
    {
        // Retrieve code
        size_t   symbol = data[i];
        Code     code   = codes[symbol];
        unsigned length = lengths[symbol];

        // Push code into the bit buffer
        constexpr unsigned codeBits = sizeof(code) * 8;
        code <<= codeBits - length; // Left align code for easier consumption by the bit buffer
        do
        {
            unsigned count = std::min(8u-bitBufferSize, length);
            bitBuffer      = uint8_t((bitBuffer << count) | (code >> (codeBits-count)));
            bitBufferSize += count;
            length        -= count;
            code         <<= count;
            if (bitBufferSize == 8u)
            {
                *writePtr++ = bitBuffer;
                bitBufferSize = 0;
            }
        }
        while (length != 0u);
    }

    // Flush bit buffer
    if (bitBufferSize != 0u)
    {
        bitBuffer <<= 8 - bitBufferSize;
        *writePtr++ = bitBuffer;
    }

    const Clock::time_point tDone = Clock::now();

    assert(size_t(writePtr - encodedData.data()) == encodedDataSize);

    printf("Encoding: %5.3f ms\n", to_ms(t0, tDone));
    printf("  Histogram:     %10.3f us\n", to_us(t0,       tHisto));
    printf("  Huffman:       %10.3f us\n", to_us(tHisto,   tHuffman));
    printf("  Compute codes: %10.3f us\n", to_us(tHuffman, tCodes));
    printf("  Encoding loop: %10.3f us\n", to_us(tCodes,   tDone));
    printf("\n");
}


void decode(std::vector<uint8_t>& decodedData, const uint8_t* encodedData, size_t encodedDataSize)
{
    (void)encodedDataSize; // Only used for assert

    const Clock::time_point t0 = Clock::now();

    // Read header
    const size_t*  weights = reinterpret_cast<const size_t*>(encodedData);
    const uint8_t* readPtr = encodedData + maxSymbolCount * sizeof(size_t);

    // Build Huffman tree
    TreeItem tree[2*maxSymbolCount-1];
    const TreeItem* treeRoot = build_huffman_tree(tree, weights);
    const Clock::time_point tHuffman = Clock::now();

    // Decode data
    uint8_t  bitBuffer     = 0;
    unsigned bitBufferSize = 0;
    uint8_t* writePtr      = decodedData.data();
    for (size_t i=0; i<decodedData.size(); ++i)
    {
        // Decode symbol by descending the tree
        const TreeItem* item = treeRoot;
        do
        {
            // Refill buffer if needed
            if (bitBufferSize == 0u)
            {
                assert(readPtr < encodedData+encodedDataSize);
                bitBuffer     = *readPtr++;
                bitBufferSize = 8;
            }

            // Read one bit and descend one level of the tree
            item = item->children[bitBuffer >> 7];
            bitBuffer <<= 1;
            bitBufferSize--;
        }
        while (item->children[0] != nullptr);

        *writePtr++ = uint8_t(item->number);
    }

    const Clock::time_point tDone = Clock::now();

    printf("Decoding: %5.3f ms\n", to_ms(t0, tDone));
    printf("  Huffman:       %10.3f us\n", to_us(t0, tHuffman));
    printf("  Decoding loop: %10.3f us\n", to_us(tHuffman, tDone));
    printf("\n");
}
