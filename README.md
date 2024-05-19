Prefix-Codes
============

This repository contains both the source code for the page [The (Complete?) Guide to Huffman and Prefix Codes](https://romigrou.github.io/prefix-codes/)
and its companion example implementations.

Those example implementations start from the easiest one and each incrementally
improves upon the previous one:
  1. [Straightforward Huffman encoding and decoding](src/01-naive-huffman.cpp)
  2. [Use the canonical representation to reduce header size and decode faster](src/02-canonical-huffman.cpp)
  3. [Build the Huffman tree with linear time complexity](src/03-fast-huffman.cpp)
  4. [Compute code lengths faster by avoiding to recurse the Huffman tree](src/04-fast-lengths.cpp)

**Note**: the code is this repository was written so as to remain as legible as possible,
          it is not meant to be the most optimized implementation around

Licensing
---------

The code in this repository is licensed under the zlib license, as follows:

> Copyright © 2024, Romain Bailly
>
> This software is provided 'as-is', without any express or implied
> warranty.  In no event will the authors be held liable for any damages
> arising from the use of this software.
>
> Permission is granted to anyone to use this software for any purpose,
> including commercial applications, and to alter it and redistribute it
> freely, subject to the following restrictions:
>
> 1. The origin of this software must not be misrepresented; you must not
>    claim that you wrote the original software. If you use this software
>    in a product, an acknowledgment in the product documentation would be
>    appreciated but is not required.
> 2. Altered source versions must be plainly marked as such, and must not be
>    misrepresented as being the original software.
> 3. This notice may not be removed or altered from any source distribution.
