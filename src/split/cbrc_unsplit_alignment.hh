// Copyright 2012 Risa Kawaguchi
// Copyright 2013, 2014 Martin C. Frith

//
// 1. UnsplitAlignment -> container of maf file information
//

#ifndef CBRC_UNSPLIT_ALIGNMENT_HH
#define CBRC_UNSPLIT_ALIGNMENT_HH

#include <stddef.h>

#include <string>
#include <vector>

namespace cbrc {

typedef std::vector<std::string>::iterator       StringIt;
typedef std::vector<std::string>::const_iterator StringCi;

class UnsplitAlignment {
public:
    StringIt linesBeg;
    StringIt linesEnd;
    const char *qname;
    unsigned qstart;
    unsigned qend;
    char qstrand;
    unsigned rstart;
    unsigned rend;
    const char *rname;
    const char *ralign;
    const char *qalign;
    const char *qQual;
    UnsplitAlignment(){}
    UnsplitAlignment(StringIt linesBegIn,
		     StringIt linesEndIn, bool isTopSeqQuery)
      : linesBeg(linesBegIn), linesEnd(linesEndIn) { init(isTopSeqQuery); }
    void init(bool isTopSeqQuery);
    bool isForwardStrand() const { return qstrand < 2; }
    bool isFlipped() const { return qstrand % 2; }
};

// Appends maf "s", "q", and "p" lines to outputText.
// Appends an extra "p" line for "probs".
// Returns the line length (including a newline).
size_t mafSlice(std::vector<char> &outputText, const UnsplitAlignment &aln,
		unsigned alnBeg, unsigned alnEnd, const double *probs);

void mafSliceBeg(const char* rAln, const char* qAln,
		 unsigned qBeg, unsigned& qSliceBeg, unsigned& alnBeg);

void mafSliceEnd(const char* rAln, const char* qAln,
		 unsigned qEnd, unsigned& qSliceEnd, unsigned& alnEnd);

double pLinesToErrorProb(const char *line1, const char *line2);

}

#endif
