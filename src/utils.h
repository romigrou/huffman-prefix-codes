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

#ifndef PREFIX_CODES_H
#define PREFIX_CODES_H


#include <cctype>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <chrono>


using Clock = std::chrono::steady_clock;

inline double to_ms(Clock::time_point t1, Clock::time_point t2) noexcept
{
    return std::chrono::duration_cast<std::chrono::nanoseconds>(t2-t1).count() / 1000000.0;
}

inline double to_us(Clock::time_point t1, Clock::time_point t2) noexcept
{
    return std::chrono::duration_cast<std::chrono::nanoseconds>(t2-t1).count() / 1000.0;
}


size_t compute_histogram(size_t histogram[256], const uint8_t data[], size_t size) noexcept;


template<typename TreeItem>
void dump_tree_node(FILE* file, const TreeItem* item, unsigned& num) noexcept
{
    const unsigned weightFontSize = 8;
    if (item->children[0] == nullptr) // Leaf
    {
        fprintf(file, "\tn%u [label=<", num);
        const unsigned symbol = unsigned(item->number);
        const char* labelFormat = (isdigit(symbol) || isalpha(symbol)) ? "'%c'" : "0x%02X";
        fprintf(file, labelFormat, symbol);
        fprintf(file, "<br/><font point-size='%u'>weight=%u</font>>];\n", weightFontSize, unsigned(item->weight));
    }
    else
    {
        const unsigned myNum = num;
        fprintf(file, "\tn%u [label=<<font point-size='%u'>weight=%u</font>>];\n", num, weightFontSize, unsigned(item->weight));
        fprintf(file, "\tn%u -- n%u [label=\"\"];\n", myNum, num+1);
        dump_tree_node(file, item->children[0], ++num);
        fprintf(file, "\tn%u -- n%u [label=\"1\"];\n", myNum, num+1);
        dump_tree_node(file, item->children[1], ++num);
    }
}


template<typename TreeItem>
void dump_tree(const char* filepath, const TreeItem* tree) noexcept
{
    FILE* file = fopen(filepath, "w");
    if (file == nullptr)
    {
        fprintf(stderr, "Failed to open file \"%s\" for writing\n", filepath);
        return;
    }

    fprintf(file, "graph {\n");
    unsigned num = 0;
    dump_tree_node(file, tree, num);
    fprintf(file, "}\n");

    fclose(file);
}


void print_canonical_codes(const uint8_t lengths[], size_t symbolCount);
void dump_canonical_tree(const char* filepath, const uint8_t lengths[], size_t symbolCount);


#endif // PREFIX_CODES_H
