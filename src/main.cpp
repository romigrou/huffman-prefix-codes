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

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <vector>


void encode(std::vector<uint8_t>& encodedData, const uint8_t* data, size_t dataSize);
void decode(std::vector<uint8_t>& decodedData, const uint8_t* encodedData, size_t encodedDataSize);


const unsigned defaultWeights[] =
{
    1192, 290, 637, 475, 1566, 254, 347, 421, 1059, 28, 155, 770, 423, 934, 1005, 444
};


const uint64_t fibonacciWeights[] =
{
    // The first 16 numbers
    1, 1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144, 233, 377, 610, 987,
    // The next 16 ones
    1597,2584,4181,6765,10946,17711,28657,46368,75025,121393,196418,317811,514229,832040,1346269,2178309,
    // The next 32 ones
    3524578,5702887,9227465,14930352,24157817,39088169,63245986,102334155,
    165580141,267914296,433494437,701408733,1134903170,1836311903,2971215073,4807526976,
    7778742049,12586269025,20365011074,32951280099,53316291173,86267571272,139583862445,225851433717,
    365435296162,591286729879,956722026041,1548008755920,2504730781961,4052739537881,6557470319842,10610209857723
};


//#define defaultWeights fibonacciWeights



int main(int argc, char* argv[])
{
    std::vector<uint8_t> data;

    if (argc > 1u)
    {
        FILE* file = fopen(argv[1], "rb");
        if (file == nullptr)
        {
            fprintf(stderr, "Failed to open file \"%s\" for reading\n", argv[1]);
            return EXIT_FAILURE;
        }
        fseek(file, 0, SEEK_END);
        int size = ftell(file);
        data.resize(size);

        fseek(file, 0, SEEK_SET);
        fread(data.data(), 1, size, file);
        fclose(file);
    }
    else
    {
        const unsigned weightCount = sizeof(defaultWeights) / sizeof(defaultWeights[0]);
        size_t dataSize = 0;
        for (unsigned i=0; i<weightCount; ++i)
            dataSize += defaultWeights[i];

        data.resize(dataSize);
        uint8_t* ptr = data.data();
        for (unsigned i=0; i<weightCount; ++i)
        {
            memset(ptr, 'A'+i, defaultWeights[i]);
            ptr += defaultWeights[i];
        }

        std::shuffle(data.data(), data.data()+dataSize, std::mt19937_64());
    }

    std::vector<uint8_t> encodedData;
    encode(encodedData, data.data(), data.size());
    printf("Compressed size: %4.2f%% of original size\n\n", (100.0*encodedData.size()) / data.size());

    std::vector<uint8_t> decodedData;
    decodedData.resize(data.size());
    decode(decodedData, encodedData.data(), encodedData.size());

    if (memcmp(decodedData.data(), data.data(), data.size()) != 0)
    {
        fprintf(stderr, "Decoded data is not the same as initial data\n");
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
