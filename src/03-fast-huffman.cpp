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
 * Changes from previous file: linear-time Huffman tree construction.
 * For this we need two things:
 *  - Linear-time sorting of symbols by weight (achieved via radix sort)
 *  - Linear-time tree building from sorted symbols
 */

#include "utils.h"
#include <algorithm>
#include <cassert>
#include <vector>


namespace
{


constexpr size_t   maxSymbolCount = 256;
constexpr unsigned maxCodeLength  = 32;
using Code = std::conditional<maxCodeLength <= 32u, uint32_t, uint64_t>::type;
using BitBuffer = Code;
static_assert(sizeof(BitBuffer) >= sizeof(Code), "");


void sort_symbols_by_weight(uint8_t symbols[maxSymbolCount], const size_t weights[]) noexcept
{
    // 10 bits per round is a good tradeoff
    constexpr unsigned bitsPerRound = 10;
    constexpr uint32_t mask = (1 << bitsPerRound) - 1;

    uint8_t histo[1 << bitsPerRound]; // 8 bits are enough in our case (256 symbols)
    uint8_t tmp[256];

    // Write to tmp in the first round: this requires one memcpy in single-round cases
    // but none in dual-round cases. This is more balanced (and we assume more than 2
    // rounds to be a rare occurrence)
    uint8_t* src = symbols;
    uint8_t* dst = tmp;

    // Compute histogram of weight bits (and compute maxWeight on the fly)
    size_t maxWeight = 0;
    memset(histo, 0, sizeof(histo));
    for (size_t symbol=0; symbol < maxSymbolCount; ++symbol)
    {
        size_t w = weights[symbol];
        maxWeight = std::max(maxWeight, w);
        histo[w & mask]++;
    }

    // Turn histogram into prefix sum
    unsigned sum = 0;
    for (size_t i=0; i < (1u << bitsPerRound); ++i)
    {
        unsigned h = histo[i];
        histo[i] = uint8_t(sum);
        sum += h;
    }

    // 1st round of LSD-first radix sort
    for (size_t symbol=0; symbol < maxSymbolCount; ++symbol)
    {
        size_t bits     = weights[symbol] & mask;
        size_t dstIndex = histo[bits]++;
        dst[dstIndex]   = uint8_t(symbol);
    }

    // Some more rounds required, keep doing the same stuff
    for (unsigned shift=bitsPerRound; (maxWeight >> shift) > 0u; shift+=bitsPerRound)
    {
        std::swap(src, dst);

        // Compute histogram
        memset(histo, 0, sizeof(histo));
        for (size_t symbol=0; symbol < maxSymbolCount; ++symbol)
            histo[(weights[symbol] >> shift) & mask]++;

        // Turn histogram into prefix sum
        sum = 0;
        for (size_t i=0; i < (1u << bitsPerRound); ++i)
        {
            unsigned h = histo[i];
            histo[i] = uint8_t(sum);
            sum += h;
        }

        // LSD-first radix sort
        for (size_t srcIndex=0; srcIndex < maxSymbolCount; ++srcIndex)
        {
            size_t symbol   = src[srcIndex];
            size_t bits     = (weights[symbol] >> shift) & mask;
            size_t dstIndex = histo[bits]++;
            dst[dstIndex]   = uint8_t(symbol);
        }
    }

    if (dst != symbols)
        memcpy(symbols, dst, maxSymbolCount);
}


struct TreeItem
{
    size_t          weight;
    uint8_t         number;
    const TreeItem* children[2] = {};
};


const TreeItem* build_huffman_tree(TreeItem tree[2*maxSymbolCount-1], const size_t weights[maxSymbolCount], const uint8_t symbolsSortedByWeight[maxSymbolCount]) noexcept
{
    // Skip unused symbols
    size_t unusedSymbolCount = 0;
    while (unusedSymbolCount < maxSymbolCount && weights[symbolsSortedByWeight[unusedSymbolCount]] == 0u)
        ++unusedSymbolCount;
    const size_t usedSymbolCount = maxSymbolCount - unusedSymbolCount;

    // Create tree leaves
    for (size_t i=0; i < usedSymbolCount; ++i)
    {
        const uint8_t symbol = symbolsSortedByWeight[unusedSymbolCount + i];
        tree[i] = {weights[symbol], symbol};
    }

    // Special case if we have less than two items
    if (usedSymbolCount <= 1u)
        return (usedSymbolCount == 1u) ? tree : nullptr;

    // Build initial node: its children are the first two leaves
    TreeItem* lastNode = tree + usedSymbolCount;
    *lastNode = {tree[0].weight + tree[1].weight, 0, {tree+0, tree+1}};

    // Now, build the rest of the tree
    TreeItem* oldestUnconsumedLeaf = tree + 2; // The first 2 leaves have just been consumed
    TreeItem* oldestUnconsumedNode = lastNode; // The first node has not yet been consumed
    TreeItem* endOfLeaves = tree + usedSymbolCount;
    for (size_t i=0; i < usedSymbolCount-2; ++i)
    {
        ++lastNode;
        lastNode->weight = 0;

        #define HUFFMAN_PICK_CHILD(num, child) \
        {                                      \
            lastNode->weight += child->weight; \
            lastNode->children[num] = child;   \
            ++child;                           \
        }

        // Select 1st child
        if (oldestUnconsumedLeaf < endOfLeaves && oldestUnconsumedLeaf->weight <= oldestUnconsumedNode->weight)
            HUFFMAN_PICK_CHILD(0, oldestUnconsumedLeaf)
        else
            HUFFMAN_PICK_CHILD(0, oldestUnconsumedNode)

        // Select 2nd child
        if ((oldestUnconsumedLeaf < endOfLeaves && oldestUnconsumedLeaf->weight <= oldestUnconsumedNode->weight) || oldestUnconsumedNode == lastNode)
            HUFFMAN_PICK_CHILD(1, oldestUnconsumedLeaf)
        else
            HUFFMAN_PICK_CHILD(1, oldestUnconsumedNode)

        #undef HUFFMAN_PICK_CHILD
    }

    return lastNode; // Tree root
}


void get_code_Lengths(uint8_t lengths[maxSymbolCount], const TreeItem* item, unsigned depth=0) noexcept
{
    if (item->children[0] == nullptr) // Leaf
    {
        lengths[item->number] = uint8_t(depth);
    }
    else                              // Node
    {
        get_code_Lengths(lengths, item->children[0], depth+1);
        get_code_Lengths(lengths, item->children[1], depth+1);
    }
}


uint8_t* sort_symbols_by_length(uint8_t symbolsSortedByLength[maxSymbolCount], const uint8_t lengths[maxSymbolCount]) noexcept
{
    for (size_t i=0; i<maxSymbolCount; ++i)
        symbolsSortedByLength[i] = uint8_t(i);
    std::stable_sort(symbolsSortedByLength, symbolsSortedByLength+maxSymbolCount, [&](size_t a, size_t b) noexcept->bool
    {
        return (lengths[a] < lengths[b] || (lengths[a] == lengths[b] && a < b));
    });

    // Skip unused symbols
    uint8_t* usedSymbols = symbolsSortedByLength;
    while (!lengths[*usedSymbols])
        ++usedSymbols;
    return usedSymbols;
}


void compute_canonical_codes(Code codes[], const uint8_t lengths[], const uint8_t symbolsSortedByLength[], size_t symbolCount) noexcept
{
    Code     nextCode   = 0;
    unsigned prevLength = 0;
    for (size_t i=0; i<symbolCount; ++i)
    {
        size_t   symbol = symbolsSortedByLength[i];
        unsigned length = lengths[symbol];
        assert(0u < length && length <= 64u);
        Code code     = nextCode << (length - prevLength);
        codes[symbol] = code;
        nextCode   = code + 1;
        prevLength = length;
    }
}


struct DecodingEntry
{
    BitBuffer firstCode;
    uint8_t   firstSymbolIndex;
};


void build_decoding_table(DecodingEntry decodingTable[maxCodeLength+2], const uint8_t usedSymbolsSortedByLength[], size_t usedSymbolCount, const uint8_t lengths[maxSymbolCount], const Code codes[maxSymbolCount]) noexcept
{
    unsigned previousLength = 0;
    for (size_t i=0; i<usedSymbolCount;)
    {
        size_t   symbol = usedSymbolsSortedByLength[i];
        unsigned length = lengths[symbol];
        assert(0u < length && length <= maxCodeLength);

        BitBuffer code = BitBuffer(codes[symbol]) << (8*sizeof(BitBuffer) - length); // Left align the code in the bit buffer
        for (unsigned l=previousLength+1; l<=length; ++l)
            decodingTable[l] = {code, uint8_t(i)};

        while (i<usedSymbolCount && lengths[usedSymbolsSortedByLength[i]]==length)
            ++i;

        previousLength = length;
    }
    decodingTable[previousLength+1] = {std::numeric_limits<Code>::max(), UINT8_MAX};
}


} // namespace



void encode(std::vector<uint8_t>& encodedData, const uint8_t* data, size_t dataSize)
{
    const Clock::time_point t0 = Clock::now();

    // Compute each symbol's weight
    size_t weights[maxSymbolCount];
    compute_histogram(weights, data, dataSize);
    const Clock::time_point tHisto = Clock::now();

    // Sort symbols by weight
    uint8_t symbolsSortedByWeight[maxSymbolCount];
    sort_symbols_by_weight(symbolsSortedByWeight, weights);
    const Clock::time_point tSortByWeight = Clock::now();

    // Build the Huffman tree
    TreeItem tree[2*maxSymbolCount-1];
    const TreeItem* treeRoot;
    for (int i=0; i<100; ++i)
    treeRoot = build_huffman_tree(tree, weights, symbolsSortedByWeight);
    const Clock::time_point tHuffman = Clock::now();

    // Retrieve symbols' lengths
    uint8_t lengths[maxSymbolCount] = {};
    get_code_Lengths(lengths, treeRoot);
    const Clock::time_point tLengths = Clock::now();

    // Sort symbols by length
    uint8_t symbolsSortedByLength[maxSymbolCount];
    const uint8_t* usedSymbols     = sort_symbols_by_length(symbolsSortedByLength, lengths);
    const size_t   usedSymbolCount = maxSymbolCount - (usedSymbols - symbolsSortedByLength);
    const Clock::time_point tSortByLength = Clock::now();

    // Compute canonical codes
    Code codes[maxSymbolCount] = {};
    compute_canonical_codes(codes, lengths, usedSymbols, usedSymbolCount);
    const Clock::time_point tCodes = Clock::now();

    // Compute the size of the encoded data
    size_t encodedDataSizeInBits = 0;
    for (size_t i=0; i<dataSize; ++i)
        encodedDataSizeInBits += lengths[data[i]];
    const size_t headerSize = maxSymbolCount; // The length of each symbol on 1 byte
    const size_t encodedDataSize = headerSize + (encodedDataSizeInBits + 7) / 8;
    encodedData.resize(encodedDataSize);

    // Write header
    uint8_t* writePtr = encodedData.data();
    memcpy(writePtr, lengths, headerSize);
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
            unsigned mask  = (1 << count) - 1;
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

    assert(writePtr - encodedData.data() == encodedDataSize);

    printf("Encoding: %5.3f ms\n", to_ms(t0, tDone));
    printf("  Histogram:       %10.3f us\n", to_us(t0,            tHisto));
    printf("  Sort by weight:  %10.3f us\n", to_us(tHisto,        tSortByWeight));
    printf("  Huffman:         %10.3f us\n", to_us(tSortByWeight, tHuffman));
    printf("  Compute lengths: %10.3f us\n", to_us(tHuffman,      tLengths));
    printf("  Sort by length:  %10.3f us\n", to_us(tLengths,      tSortByLength));
    printf("  Compute codes:   %10.3f us\n", to_us(tSortByLength, tCodes));
    printf("  Encoding itself: %10.3f us\n", to_us(tCodes,        tDone));
    printf("\n");
}


void decode(std::vector<uint8_t>& decodedData, const uint8_t* encodedData, size_t encodedDataSize)
{
    const Clock::time_point t0 = Clock::now();

    constexpr unsigned bitBufferCapacity = sizeof(BitBuffer) * 8;

    // Read header
    const uint8_t* lengths = encodedData;
    const uint8_t* readPtr = encodedData + maxSymbolCount;

    // Sort symbols by length
    uint8_t symbolsSortedByLength[maxSymbolCount];
    const uint8_t* usedSymbols     = sort_symbols_by_length(symbolsSortedByLength, lengths);
    const size_t   usedSymbolCount = maxSymbolCount - (usedSymbols - symbolsSortedByLength);
    const Clock::time_point tSortByLength = Clock::now();

    // Compute canonical codes
    Code codes[maxSymbolCount] = {};
    compute_canonical_codes(codes, lengths, usedSymbols, usedSymbolCount);
    const Clock::time_point tCodes = Clock::now();

    // Build decoding table
    DecodingEntry decodingTable[maxSymbolCount+2];
    build_decoding_table(decodingTable, usedSymbols, usedSymbolCount, lengths, codes);
    const Clock::time_point tTable = Clock::now();

    // Decode data
    const unsigned minLength = lengths[usedSymbols[0]];
    Code     bitBuffer     = 0;
    unsigned bitBufferSize = 0;
    uint8_t* writePtr      = decodedData.data();
    for (size_t i=0; i<decodedData.size(); ++i)
    {
        // Figure out symbol's length
        unsigned length = minLength + 1;
        for (;;)
        {
            while (length < bitBufferSize && bitBuffer >= decodingTable[length].firstCode)
                ++length;
            if (length < bitBufferSize)
                break;

            // Read some more bits and retry
            bitBuffer     |= *readPtr++ << (bitBufferCapacity - bitBufferSize - 8);
            bitBufferSize += 8;
            assert(bitBufferSize <= bitBufferCapacity);
        }
        --length;

        // Decode symbol
        const DecodingEntry& e = decodingTable[length];
        unsigned shift     = bitBufferCapacity - length; // To right align the code
        Code     code      = bitBuffer >> shift;
        Code     firstCode = decodingTable[length].firstCode >> shift;
        uint8_t  symbol    = usedSymbols[e.firstSymbolIndex + (code - firstCode)];
        *writePtr++ = symbol;
        bitBuffer    <<= length;
        bitBufferSize -= length;
    }

    const Clock::time_point tDone = Clock::now();

    printf("Decoding: %5.3f ms\n", to_ms(t0, tDone));
    printf("  Sort by length:  %10.3f us\n", to_us(t0,            tSortByLength));
    printf("  Compute codes:   %10.3f us\n", to_us(tSortByLength, tCodes));
    printf("  Build table:     %10.3f us\n", to_us(tCodes,        tTable));
    printf("  Decoding itself: %10.3f us\n", to_us(tTable,        tDone));
    printf("\n");
}
