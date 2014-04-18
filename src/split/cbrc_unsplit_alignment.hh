// Copyright 2012 Risa Kawaguchi
// Copyright 2013, 2014 Martin C. Frith

//
// 1. UnsplitAlignment -> container of maf file information
//

#ifndef CBRC_UNSPLIT_ALIGNMENT_HH
#define CBRC_UNSPLIT_ALIGNMENT_HH

#include <string>
#include <vector>

namespace cbrc {

typedef std::vector<std::string>::iterator       StringIt;
typedef std::vector<std::string>::const_iterator StringCi;

class UnsplitAlignment {
public:
    std::string qname;
    unsigned qstart;
    unsigned qend;
    char qstrand;
    unsigned rstart;
    unsigned rend;
    std::string rname;
    std::string ralign;
    std::string qalign;
    std::string qQual;
    std::vector<std::string> lines;
    UnsplitAlignment(){}
    UnsplitAlignment(const std::vector<std::string>& linesIn)
      : lines(linesIn) { init(); }
    void init();
};

void flipMafStrands(StringIt linesBeg, StringIt linesEnd);

std::vector<std::string> mafSlice(StringCi linesBeg, StringCi linesEnd,
				  unsigned alnBeg, unsigned alnEnd);

void mafSliceBeg(const std::string& rAln, const std::string& qAln,
		 unsigned qBeg, unsigned& qSliceBeg, unsigned& alnBeg);

void mafSliceEnd(const std::string& rAln, const std::string& qAln,
		 unsigned qEnd, unsigned& qSliceEnd, unsigned& alnEnd);

void printMaf(const std::vector<std::string>& maf);

std::string pLineFromProbs(const std::vector<double>& p);

}

#endif
