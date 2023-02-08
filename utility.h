/* The MIT License (MIT)

   Copyright (c) 2023 Anrui Liu <liuar6@gmail.com>

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   “Software”), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
 */

#include <htslib/sam.h>

#define getAuxInteger(b, tag) bam_aux2i(bam_aux_get((b), (tag)))
bool isProper(bam1_t* b, bool isPaired, bool unique){
    if (isPaired && !isProperPair(b)) return false;
    if (unique && getAuxInteger(b, "NH")!=1) return false;
    return true;
}
#define FRUNSTRANDED 0
#define FRFIRSTSTRAND 1
#define FRSECONDSTRAND 2
char getStrand(bam1_t *b, bool isPaired, int libraryType){
    if (libraryType==FRFIRSTSTRAND){
        if (isPaired){
            if ((isFirstMate(b) && isReverse(b)) || (isSecondMate(b) && isMateReverse(b)))
                return '+';
            else if ((isFirstMate(b) && isMateReverse(b)) || (isSecondMate(b) && isReverse(b)))
                return '-';
        }
        else ((isReverse(b)?'+':'-'));
    }
    else if (libraryType==FRSECONDSTRAND){
        if (isPaired){
            if ((isFirstMate(b) && isMateReverse(b)) || (isSecondMate(b) && isReverse(b)))
                return '+';
            else if ((isFirstMate(b) && isReverse(b)) || (isSecondMate(b) && isMateReverse(b)))
                return '-';
        }
        else ((isReverse(b)?'-':'+'));
    }
    else if (libraryType==FRUNSTRANDED) return '.';
    exit(1);
}
#define max(a, b) (((a)>(b))?(a):(b))
#define min(a, b) (((a)>(b))?(b):(a))