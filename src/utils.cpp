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

#include "utils.h"
#include <algorithm>
#include <cstring>
#include <string>


size_t compute_histogram(size_t histogram[256], const uint8_t data[], size_t size) noexcept
{
    memset(histogram, 0, 256*sizeof(histogram[0]));
    for (size_t i=0; i<size; ++i)
        histogram[data[i]]++;

    size_t usedSymbolCount = 0;
    for (size_t i=0; i<256u; ++i)
        if (histogram[i] != 0u)
            ++usedSymbolCount;

    return usedSymbolCount;
}


static std::string as_binary(unsigned val, unsigned length) noexcept
{
    char buffer[65];
    for (unsigned i=0; i<length; ++i)
        buffer[i] = '0' + ((val >> (length-1-i)) & 1);
    buffer[length] = 0;
    return buffer;
}


void print_canonical_codes(const uint8_t lengths[], size_t symbolCount)
{
    // Sort symbols by length
    const size_t maxSymbolCount = 256;
    uint8_t symbolsSortedByLength[maxSymbolCount];
    for (size_t i=0; i<symbolCount; ++i)
        symbolsSortedByLength[i] = uint8_t(i);
    std::stable_sort(symbolsSortedByLength, symbolsSortedByLength+symbolCount, [&](size_t a, size_t b) noexcept->bool
    {
        return (lengths[a] < lengths[b]);
    });

    // Skip zero-length symbols
    const uint8_t* usedSymbols = symbolsSortedByLength;
    while (!lengths[*usedSymbols])
        ++usedSymbols;
    const size_t usedSymbolCount = symbolCount - (usedSymbols - symbolsSortedByLength);

    // Print codes for used symbols
    unsigned nextCode   = 0;
    unsigned prevLength = 0;
    for (unsigned i=0; i<usedSymbolCount; ++i)
    {
        unsigned symbol = usedSymbols[i];
        const char* symbolFormat = (isdigit(symbol) || isalpha(symbol)) ? "'%c':  " : "0x%02X: ";
        printf(symbolFormat, symbol);

        unsigned length = lengths[symbol];
        unsigned code   = nextCode << (length - prevLength);
        printf("%s\n", as_binary(code, length).c_str());

        nextCode   = code + 1;
        prevLength = length;
    }
}


static void print_edge(FILE* file, unsigned code, unsigned length)
{
    fprintf(file, "\tn%s -- n%s [label=\"%u\"];\n", as_binary(code>>1, length-1).c_str(), as_binary(code, length).c_str(), code&1);
}


static void print_leaf(FILE* file, unsigned symbol, unsigned code, unsigned length)
{
    print_edge(file, code, length);
}

void dump_canonical_tree(const char* filepath, const uint8_t lengths[], size_t symbolCount)
{
    // Sort symbols by length
    const size_t maxSymbolCount = 256;
    uint8_t symbolsSortedByLength[maxSymbolCount];
    for (size_t i=0; i<symbolCount; ++i)
        symbolsSortedByLength[i] = uint8_t(i);
    std::stable_sort(symbolsSortedByLength, symbolsSortedByLength+symbolCount, [&](size_t a, size_t b) noexcept->bool
    {
        return (lengths[a] < lengths[b]);
    });

    // Skip zero-length symbols
    const uint8_t* usedSymbols = symbolsSortedByLength;
    while (!lengths[*usedSymbols])
        ++usedSymbols;
    const size_t usedSymbolCount = symbolCount - (usedSymbols - symbolsSortedByLength);

    FILE* file = fopen(filepath, "w");
    if (file == nullptr)
    {
        fprintf(stderr, "Failed to open file \"%s\" for writing\n", filepath);
        return;
    }
    fprintf(file, "graph {\n");

    unsigned code = 0;
    fprintf(file, "\tn [label=\"\"];\n");
    for (unsigned length=1, usedSymbolNum=0; usedSymbolNum < usedSymbolCount; ++length)
    {
        for (;usedSymbolNum < usedSymbolCount && lengths[usedSymbols[usedSymbolNum]] == length; ++usedSymbolNum, ++code)
        {
            const unsigned symbol = usedSymbols[usedSymbolNum];
            fprintf(file, "\tn%s [label=<", as_binary(code,length).c_str());
            const char* labelFormat = (isdigit(symbol) || isalpha(symbol)) ? "'%c'" : "0x%02X";
            fprintf(file, labelFormat, symbol);
            fprintf(file, ">];\n");
            print_edge(file, code, length);
        }

        for (unsigned c=code; c < (1u << length); ++c)
        {
            fprintf(file, "\tn%s [label=\"\"];\n", as_binary(c,length).c_str());
            print_edge(file, c, length);
        }

        code <<= 1;
    }

    fprintf(file, "}\n");
    fclose(file);
}
