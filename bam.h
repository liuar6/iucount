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

#include <unordered_map>
#include <map>
#include <vector>
#include <string.h>
#include <htslib/bgzf.h>
#include <htslib/sam.h>
using namespace std;

#define isUnmapped(b) (((b)->core.flag & BAM_FUNMAP) != 0)
#define isReverse(b) (((b)->core.flag&BAM_FREVERSE) != 0)
#define isMateReverse(b) (((b)->core.flag&BAM_FMREVERSE) != 0)
#define isProperPair(b) (((b)->core.flag & BAM_FPROPER_PAIR) !=0u)
#define isFirstMate(b) (((b)->core.flag & BAM_FREAD1) !=0u)
#define isSecondMate(b) (((b)->core.flag & BAM_FREAD2) !=0u)

//#define getStrand(b) ((isReverse(b)?'-':'+'))
#define getName(b) ((char*)(b)->data)
#define getReferenceID(b) ((b)->core.tid)
#define getReferenceName(b, h) (h->target_name[(b)->core.tid])
#define getPosition(b) ((b)->core.pos)

#define getCigar(b) ((uint32_t*)((b)->data + (b)->core.l_qname))
#define getCigarNum(b) ((b)->core.n_cigar)
#define getCigarOp(c) (c & 0xF)
#define getCigarOpchr(c) (BAM_CIGAR_STR[c & 0xFu])
#define getCigarOplen(c) ((c & 0xFFFFFFF0u)>>4u)
#define getCigarType(c) (BAM_CIGAR_TYPE>>((c & 0xFu)<<1u)&3u)

#define getQueryLength(b) ((b)->core.l_qseq)
#define getQuery(b) ((b)->data + ((b)->core.n_cigar<<2) + (b)->core.l_qname)
#define getQuality(b)  ((b)->data + ((b)->core.n_cigar<<2) + (b)->core.l_qname + (((b)->core.l_qseq + 1)>>1))
#define getBase(s, i) (seq_nt16_str[(s)[(i)>>1] >> ((~(i)&1)<<2) & 0xf])


class bamReader{
public:
    BGZF *fp;
    bam_hdr_t *header;
    bam1_t *b;
    bam1_t *b2;
    unordered_map<string, uint32_t> chr2index;
    unordered_map<int, uint64_t> offset;

    explicit bamReader(const char *fn){
        fp=nullptr;
        header=nullptr;
        b=nullptr;
        b2=nullptr;
        open(fn);
    }
    explicit bamReader(){
        fp=nullptr;
        header=nullptr;
        b=nullptr;
        b2=nullptr;
    }
    int open(const char *fn){
        fp=bgzf_open(fn, "r");
        header=bam_hdr_read(fp);
        b=bam_init1();
        b2=bam_init1();
        for (int i=0; i<header->n_targets; ++i) chr2index[header->target_name[i]]=i;
        return 1;
    }
    int reopen(){
        if (header) bam_hdr_destroy(header);
        if (b) bam_destroy1(b);
        if (b2) bam_destroy1(b2);
        bgzf_seek(fp, 0, SEEK_SET);
        header=bam_hdr_read(fp);
        b=bam_init1();
        b2=bam_init1();
        return 1;
    }
    int jumpToChrom(const char* chrom){
        return bgzf_seek(fp, offset[chr2index[chrom]], SEEK_SET);
    }
    int seek(int id){
        return bgzf_seek(fp, offset[id], SEEK_SET);
    }
    bam1_t* next(){
        if (bam_read1(fp, b)>0) return b;
        else return nullptr;
    }
    bam1_t* next2(){
        if (bam_read1(fp, b2)>0) return b2;
        else return nullptr;
    }
    static int return_with_error(const char* message, int ret){
        if (ret){
            cerr<<message<<endl;
            return ret;
        }
        return 0;
    }
    int parse(){
        reopen();
        int last_coor=0;
        int last_tid=-1;
        u_int64_t last_offset=bgzf_tell(fp);
        int tid;
        bool no_coor=false;
        while (next()){
            tid=b->core.tid;
            if (isUnmapped(b)) no_coor=true;
            if (tid>=0 && last_tid!=tid){
                if (no_coor) return_with_error("unsorted bam", -1);
                if (offset.find(tid)!=offset.end()) return_with_error("chromosome not continuous", -1);
                else offset[tid]=last_offset;
            }
            else if (tid >= 0 && last_coor > (b->core.pos)) return_with_error("unsorted bam", 0);
            last_tid=tid;
            last_offset=bgzf_tell(fp);
            last_coor=b->core.pos;
        }
        reopen();
        return 1;
    }
    ~bamReader(){
        if (fp) bgzf_close(fp);
        if (header) bam_hdr_destroy(header);
        if (b) bam_destroy1(b);
        if (b2) bam_destroy1(b2);
    }
};