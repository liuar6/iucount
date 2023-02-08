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

#include <iostream>
#include <unordered_map>
#include <vector>
#include <fstream>
#include <string.h>
#include <getopt.h>
#include <algorithm>
#include "bam.h"
#include "utility.h"

struct Parameter{
    char* intronFile;
    char* bamFile;
    char* outFile;
    bool isPaired;
    int libraryType;
    int span;
    int readLen;
    bool calculate;
    bool unique;
};
struct Parameter *parameters;
char * newstr(const char *str){
    char *cstr=new char[strlen(str)+1];
    strcpy(cstr, str);
    return cstr;
}
char ** split(char * line, char ** results, int length,char c='\t'){
    char *start=line;
    char *end=nullptr;
    int i=0;
    while ((end=strchr(start, c))!=nullptr && i<length){
        end[0]='\0';
        results[i]=start;
        start=end+1;
        i=i+1;
    }
    if (i<length && start[0]!='\0') {
        results[i]=start;
        i=i+1;
    }
    for (;i<length;++i) results[i]=nullptr;
    return results;
}
using namespace std;
struct Intron{
    char* chrom;
    uint32_t chromStart;
    uint32_t chromEnd;
    char strand;
    int incCount=0;
    int cntCount=0;
    int skipCount=0;
};
void readIntrons(vector<struct Intron*> *introns){
    char line[100000];
    char *items[8];
    ifstream infile;
    infile.open(parameters->intronFile);
    infile.getline(line, 100000);
    while (!infile.eof()){
        if (line[0]!='#'){
            split(line, items, 8);
            auto intron=new struct Intron;
            intron->chrom=newstr(items[0]);
            intron->chromStart=strtol(items[1], nullptr, 10);
            intron->chromEnd=strtol(items[2], nullptr, 10);
            intron->strand=*items[3];
            introns->push_back(intron);
        }
        infile.getline(line, 1000000);
    }
    infile.close();
};
void deleteIntron(struct Intron* intron){
    delete [] intron->chrom;
    delete intron;
}
void deleteIntrons(vector<struct Intron*> *introns){
    for (auto intron: *introns) deleteIntron(intron);
}

bool compare(struct Intron* i, struct Intron* j){
    return i->chromStart<j->chromStart;
}
int countInc(int32_t chromStart, int32_t chromEnd, char strand, vector<struct Intron*>* index, int guide, unordered_map<struct Intron *,int>* incRecorder, unordered_map<struct Intron *,int>* cntRecorder){
    int i=guide;
    int32_t overlap, leftSpan, rightSpan;
    while(i<index->size() && (*index)[i]->chromEnd<=chromStart) ++i;
    guide=i;
    while(i<index->size() && (*index)[i]->chromStart<chromEnd){
        struct Intron* intron=(*index)[i];
        if (strand == '.' || strand == (*index)[i]->strand){
            //an overlap less than zero mean there is no overlap
            overlap=min((int32_t)intron->chromEnd, chromEnd)-max((int32_t)intron->chromStart, chromStart);
            if (overlap>=parameters->span){
                leftSpan=(int32_t)intron->chromStart-chromStart;
                rightSpan=chromEnd-(int32_t)intron->chromEnd;
                if (leftSpan>=parameters->span || rightSpan>=parameters->span) (*incRecorder)[intron]=1;
                if (leftSpan<=0 && rightSpan<=0) (*cntRecorder)[intron]=1;
            }
        }
        ++i;
    }
    return guide;
}
void countSkip(int32_t chromStart, int32_t chromEnd, int32_t lastChromStart, int32_t lastChromEnd, char strand, unordered_map<uint64_t, struct Intron*>* index){
    if (lastChromEnd-lastChromStart<=parameters->span || chromEnd-chromStart<=parameters->span) return;
    unordered_map<uint64_t, struct Intron*>::iterator it;
    if (strand=='+' || strand=='.'){
        uint64_t i=((uint64_t)(lastChromEnd)<<32u)+chromStart;
        it=index->find(i);
        if (it!=index->end()) it->second->skipCount++;
    }
    if (strand=='-' || strand=='.'){
        uint64_t i=((uint64_t)(chromStart)<<32u)+lastChromEnd;
        it=index->find(i);
        if (it!=index->end()) it->second->skipCount++;
    }
}
void parseArgs(int, char *[]);
void calculateEffectiveLength(vector<struct Intron*>*);
int main(int argc, char *argv[]){
    parseArgs(argc, argv);
    //read introns
    cerr<<"loading intron file"<<endl;
    auto introns= new vector<struct Intron*>;
    readIntrons(introns);

    //check if calculate is true
    if (parameters->calculate) calculateEffectiveLength(introns);

    //read bam file
    cerr<<"loading bam file"<<endl;
    bamReader bam;
    if (!bam.open(parameters->bamFile) || !bam.parse()) exit(1);

    //prepare index for introns
    auto skipIndices=new unordered_map<uint64_t, struct Intron*> [bam.chr2index.size()];
    for (auto i:* introns) {
        if (i->strand=='+') skipIndices[bam.chr2index[i->chrom]][((uint64_t)(i->chromStart)<<32u)+i->chromEnd]=i;
        else if (i->strand=='-') skipIndices[bam.chr2index[i->chrom]][((uint64_t)(i->chromEnd)<<32u)+i->chromStart]=i;
    }
    auto incIndices=new vector<struct Intron*>[bam.chr2index.size()];
    for (auto i:* introns)  incIndices[bam.chr2index[i->chrom]].push_back(i);
    for (auto i=0; i<bam.offset.size(); ++i) sort(incIndices[i].begin(), incIndices[i].end(), compare);
    auto incRecorder=new unordered_map<struct Intron *, int>;
    auto cntRecorder=new unordered_map<struct Intron *, int>;

    int readCount=0;
    for (auto it=bam.offset.begin(); it != bam.offset.end(); ++it)
    {
        bam1_t* b;
        int chromId=it->first;
        bam.seek(chromId);
        cerr<<"processing "<<bam.header->target_name[chromId]<<endl;

        int32_t lastPosition=0;
        auto incIndex=&incIndices[chromId];
        auto skipIndex=&skipIndices[chromId];
        int guide=0;

        while ((b=bam.next())!=nullptr && b->core.tid==chromId){
            incRecorder->clear();
            cntRecorder->clear();

            if (lastPosition>getPosition(b)) exit(12);
            lastPosition=getPosition(b);

            if (!isProper(b, parameters->isPaired, parameters->unique)) continue;

            char strand=getStrand(b, parameters->isPaired, parameters->libraryType);
            const uint32_t *cigar=getCigar(b);
            const int cigarNum=getCigarNum(b);
            int32_t chromStart, chromEnd, lastChromStart, lastChromEnd;
            chromStart=getPosition(b);
            chromEnd=chromStart;
            lastChromStart=lastChromEnd=0;
            int i=0;
            int subGuide=guide;
            while(i<cigarNum){
                while(i<cigarNum && getCigarOp(cigar[i])!=3){ //Until an 'N' is reached or the cigar is finished
                    if (getCigarType(cigar[i]) & 2u) chromEnd+=getCigarOplen(cigar[i]);
                    ++i;
                }
                subGuide=countInc(chromStart, chromEnd, strand, incIndex, subGuide, incRecorder, cntRecorder);
                /*if this is the first segment, record the index value returned, then for
                 * next read the index of first overlaping intron should be more or equal to this
                 * index                                                               */
                if (lastChromStart==0) guide=subGuide;
                // nothing will be count if there is no last positions
                countSkip(chromStart, chromEnd, lastChromStart, lastChromEnd, strand, skipIndex);
                if (i<cigarNum){
                    lastChromStart=chromStart;
                    lastChromEnd=chromEnd;
                    chromStart=chromEnd+getCigarOplen(cigar[i]);
                    chromEnd=chromStart;
                    ++i;
                }
            }
            for (auto intron:*incRecorder) intron.first->incCount++;
            for (auto intron:*cntRecorder) intron.first->cntCount++;
        }
    }
    if (parameters->outFile!=nullptr){
        ofstream outfile;
        outfile.open(parameters->outFile);
        for (auto i:*introns){
            outfile<<i->chrom<<'\t'<<i->chromStart<<'\t'<<i->chromEnd<<'\t'<<i->strand<<'\t'<<i->incCount<<'\t'<<i->cntCount<<'\t'<<i->skipCount<<endl;
        }
        outfile.close();
    }
    else for (auto i:*introns){
            cerr<<i->chrom<<'\t'<<i->chromStart<<'\t'<<i->chromEnd<<'\t'<<i->strand<<'\t'<<i->incCount<<'\t'<<i->cntCount<<'\t'<<i->skipCount<<endl;
        }
    delete []incIndices;
    delete []skipIndices;
    delete incRecorder;
    delete cntRecorder;
    deleteIntrons(introns);
    return 0;
}
void calculateEffectiveLength(vector<struct Intron*>* introns){
    cerr<<"the read length is "<<parameters->readLen<<endl;
    if (parameters->readLen==-1) exit(1);
    int skipFormLen, incFormLen, cntFormLen, intronLen;
    int span=parameters->span;
    int readLen=parameters->readLen;
    ofstream outfile;
    outfile.open(parameters->outFile);
    for (auto intron: *introns){
        intronLen=intron->chromEnd-intron->chromStart;
        skipFormLen=readLen-2*span+1;
        incFormLen=readLen-2*span+1+min(readLen-2*span+1, intronLen);
        cntFormLen=incFormLen+max(0, intronLen-readLen+1);
        outfile<<intron->chrom<<'\t'<<intron->chromStart<<'\t'<<intron->chromEnd<<'\t'<<intron->strand<<
        '\t'<<incFormLen<<'\t'<<cntFormLen<<'\t'<<skipFormLen<<endl;
    }
    outfile.close();
    exit(0);
}
void usage()
{
    fprintf(stderr, "%s", "iucount: count intron usage information from alignment file.\n\
Usage:  iucount [options] --bam <alignment file> --intron <intron file>\n \
[options]\n\
-i/--intron                    : intron file.[required]\n\
-b/--bam                       : sorted bam alignment file.[required]\n\
-v/--version                   : show version\n\
-h/--help                      : show help informations\n\
-t/--library-type              : library type, one of fr-firststrand, fr-secondstrand or fr-unstranded\n\
-p/--paried                    : whether the library is paired-end\n\
-u/--unique                    : use unique alignment only.\n\
-s/--span                      : minimal span for segments when counting \"skip\" and \"include\", default 6. \n\
-r/--read-length               : read length of the library, currently no need to provide except for -c. \n\
-c/--calculate                 : calculate the effective length for each intron, read length must be provided.\n\
-o/--output                    : output file\n\
");

    exit(1);
}
void parseArgs(int argc, char *argv[]){
    char c;
    int showHelp=0;
    int showVersion=0;
    parameters=new struct Parameter;
    const char* version="1.0.0";
    if (argc == 1) //if no arguments are provided
    {
        //usage();
    }
    const char *shortOptions = "vhcpuo:b:i:t:s:r:";
    const struct option longOptions[] =
            {
                    { "help" , no_argument , NULL, 'h' },
                    { "version" , no_argument , NULL, 'V' },
                    { "output" , required_argument , NULL, 'o' },
                    { "bam" , required_argument, NULL, 'b' },
                    { "intron" , required_argument, NULL, 'i' },
                    { "library-type" , required_argument, NULL, 't' },
                    { "unique" , no_argument, NULL, 'u' },
                    { "span" , required_argument, NULL, 's' },
                    { "calculate" , no_argument, NULL, 'c' },
                    { "read-length" , required_argument, NULL, 'r' },
                    { "paired" , no_argument, NULL, 'p' },
                    {NULL, 0, NULL, 0} ,  /* Required at end of array. */
            };

    parameters->intronFile=nullptr;
    parameters->bamFile=nullptr;
    parameters->outFile=nullptr;
    parameters->isPaired=false;
    parameters->unique=false;
    parameters->libraryType=FRUNSTRANDED;
    parameters->span=6;
    parameters->calculate=false;
    parameters->readLen=-1;

    while ((c = getopt_long(argc, argv, shortOptions, longOptions, NULL)) >= 0)
    {
        switch (c)
        {
            case 'h':
                showHelp = 1;
                break;
            case 'v':
                showVersion = 1;
                break;
            case 'o':
                parameters->outFile=optarg;
                break;
            case 'b':
                parameters->bamFile=optarg;
                break;
            case 'i':
                parameters->intronFile=optarg;
                break;
            case 't':
                if (strcmp(optarg,"fr-firststrand")==0) parameters->libraryType=FRFIRSTSTRAND;
                else if (strcmp(optarg,"fr-secondstrand")==0) parameters->libraryType=FRSECONDSTRAND;
                else if (strcmp(optarg,"fr-unstranded")==0) parameters->libraryType=FRUNSTRANDED;
                else exit(1);
                break;
            case 'p':
                parameters->isPaired=true;
                break;
            case 'u':
                parameters->unique=true;
                break;
            case 's':
                parameters->span=strtol(optarg, nullptr, 10);
                break;
            case 'c':
                parameters->calculate=true;
                break;
            case 'r':
                parameters->readLen=strtol(optarg, nullptr, 10);
                break;
            case '?':
                showHelp = 1;
                break;
            default:
                usage();
        }
    }
    //if (argc != optind) usage();

    if (showVersion)
    {
        fprintf(stderr, "iucount-%s\n", version);
        exit(1);
    }

    if (showHelp)
    {
        usage();
        exit(1);
    }
    if (!parameters->bamFile && !parameters->calculate) {
        cerr<<"[warning] bam file not provided, read from standard input"<<endl;
        parameters->bamFile = "/dev/stdin";
    }

    if (!parameters->outFile) {
        cerr<<"[warning] out file not provided, write result to standard output"<<endl;
        parameters->outFile = "/dev/stdout";
    }

    if (!parameters->intronFile) {
        cerr<<"[error] please provide intron file"<<endl;
        exit(1);
    }
}