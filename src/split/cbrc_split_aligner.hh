// Copyright 2012 Risa Kawaguchi
// Copyright 2013 Martin C. Frith

#ifndef CBRC_SPLIT_ALIGNER_HH
#define CBRC_SPLIT_ALIGNER_HH

#include "cbrc_unsplit_alignment.hh"

#include "Alphabet.hh"
#include "MultiSequence.hh"

#include <string>
#include <vector>
#include <cmath>
#include <climits>
#include <map>

namespace cbrc {

class SplitAligner {
public:
    SplitAligner() {
      setParams(-7, -1, -40, INT_MIN/2, 1.0);  // xxx ???
      setSpliceParams(0.0, 7.0, 1.7);
    }

    // A gap of length k scores: gapExistenceScore + k * gapExtensionScore.

    // "jumpScore" is the (negative) score for a trans-splice.

    // "restartScore" is the (negative) score for re-starting an
    // alignment, in the repeated matches algorithm from chapter 2 of
    // Durbin, Eddy, et al.  In that book, it is called "-T".

    void setParams(int gapExistenceScoreIn, int gapExtensionScoreIn,
		   int jumpScoreIn, int restartScoreIn, double scaleIn);

    void setSpliceParams(double splicePriorIn,
			 double meanLogDistIn, double sdevLogDistIn);

    void setScoreMat(const std::vector< std::vector<int> >& matrix,
		     const std::string& rowNames,
		     const std::string& colNames);

    void readGenome(const std::string& baseName);

    // XXX this should allow us to specify scores for gt-ag, at-ac, etc.
    void setSpliceSignals();

    // Outputs some algorithm parameters on lines starting with "#"
    void printParameters() const;

    // Prepares to analyze some candidate alignments for one query sequence
    void initForOneQuery(std::vector<UnsplitAlignment>::const_iterator beg,
			 std::vector<UnsplitAlignment>::const_iterator end);

    long viterbi();  // returns the optimal split-alignment score

    // Gets the chunks of an optimal split alignment.
    // For each chunk, it gets:
    // 1. The index of the candidate alignment that the chunk comes from
    // 2. The chunk's start coordinate in the query sequence
    // 3. The chunk's end coordinate in the query sequence
    // It gets the chunks in reverse order, from query end to query start.
    void traceBack(long viterbiScore,
		   std::vector<unsigned>& alnNums,
		   std::vector<unsigned>& queryBegs,
		   std::vector<unsigned>& queryEnds) const;

    // Calculates the alignment score for a segment of an alignment
    int segmentScore(unsigned alnNum,
		     unsigned queryBeg, unsigned queryEnd) const;

    void forward();

    void backward();

    // Returns one probability per column, for a segment of an alignment
    std::vector<double> marginalProbs(unsigned queryBeg, unsigned alnNum,
				      unsigned alnBeg, unsigned alnEnd) const;

private:
    typedef std::vector< std::vector<int> > MatrixInt;
    typedef std::vector< std::vector<long> > MatrixLong;
    typedef std::vector< std::vector<unsigned> > MatrixUnsigned;
    typedef std::vector< std::vector<double> > MatrixDouble;

    static const int numQualCodes = 64;
    static int score_mat[64][64][numQualCodes];
    int gapExistenceScore;
    int gapExtensionScore;
    int jumpScore;
    int restartScore;
    double jumpProb;
    double restartProb;
    double scale;
    unsigned numAlns;  // the number of candidate alignments (for 1 query)
    std::vector<UnsplitAlignment>::const_iterator alns;  // the candidates
    unsigned minBeg;  // the minimum query start coordinate of any candidate
    unsigned maxEnd;  // the maximum query end coordinate of any candidate
    std::vector<unsigned> dpBegs;  // dynamic programming begin coords
    std::vector<unsigned> dpEnds;  // dynamic programming end coords
    MatrixInt Amat;  // scores at query bases, for each candidate
    MatrixInt Dmat;  // scores between query bases, for each candidate
    MatrixLong Vmat;  // DP matrix for Viterbi algorithm
    std::vector<long> Vvec;  // DP vector for Viterbi algorithm
    MatrixDouble Aexp;
    MatrixDouble Dexp;
    MatrixDouble Fmat;  // DP matrix for Forward algorithm
    MatrixDouble Bmat;  // DP matrix for Backward algorithm
    std::vector<double> rescales;  // the usual scaling for numerical stability
    double zF;  // sum of probabilities from the forward algorithm
    double zB;  // sum of probabilities from the backward algorithm

    std::vector<unsigned> sortedAlnIndices;
    std::vector<unsigned> oldInplayAlnIndices;
    std::vector<unsigned> newInplayAlnIndices;

    double splicePrior;
    double meanLogDist;
    double sdevLogDist;
    double spliceTerm1;
    double spliceTerm2;
    MatrixUnsigned spliceBegCoords;
    MatrixUnsigned spliceEndCoords;
    std::vector<unsigned> rnameAndStrandIds;
    MultiSequence genome;
    Alphabet alphabet;
    typedef std::map<std::string, unsigned> StringNumMap;
    StringNumMap chromosomeIndex;
    int spliceBegScores[4 * 4 + 1];  // donor score for any dinucleotide
    int spliceEndScores[4 * 4 + 1];  // acceptor score for any dinucleotide
    double spliceBegProbs[4 * 4 + 1];
    double spliceEndProbs[4 * 4 + 1];
    unsigned spliceBegSignal(unsigned coordinate, char strand) const;
    unsigned spliceEndSignal(unsigned coordinate, char strand) const;
    int spliceBegScore(unsigned i, unsigned j) const;
    int spliceEndScore(unsigned i, unsigned j) const;
    double spliceBegProb(unsigned i, unsigned j) const;
    double spliceEndProb(unsigned i, unsigned j) const;
    int spliceScore(double dist) const;
    double spliceProb(double dist) const
    { return std::exp(spliceScore(dist) / scale); }
    void initSpliceCoords();
    void initRnameAndStrandIds();

    void updateInplayAlnIndicesF(unsigned& sortedAlnPos,
				 unsigned& oldNumInplay,
				 unsigned& newNumInplay, unsigned j);

    void updateInplayAlnIndicesB(unsigned& sortedAlnPos,
				 unsigned& oldNumInplay,
				 unsigned& newNumInplay, unsigned j);

    unsigned findScore(unsigned j, long score) const;
    unsigned findSpliceScore(unsigned i, unsigned j, long score) const;
    long scoreFromSplice(unsigned i, unsigned j, unsigned oldNumInplay,
			 unsigned& oldInplayPos) const;
    long endScore() const;
    unsigned findEndScore(long score) const;

    // "dp" means "dynamic programming":
    unsigned dpBeg(unsigned i) const { return dpBegs[i]; }
    unsigned dpEnd(unsigned i) const { return dpEnds[i]; }

    template<typename T> T&
    cell(std::vector<T>& v, unsigned j) const
    { return v[j - minBeg]; }

    template<typename T> const T&
    cell(const std::vector<T>& v, unsigned j) const
    { return v[j - minBeg]; }

    template<typename T> T&
    cell(std::vector< std::vector<T> >& v, unsigned i, unsigned j) const
    { return v[i][j - dpBeg(i)]; }

    template<typename T> const T&
    cell(const std::vector< std::vector<T> >& v, unsigned i, unsigned j) const
    { return v[i][j - dpBeg(i)]; }

    template<typename T>
    void resizeVector(T& v, int extraCells) const
    { v.resize(maxEnd - minBeg + extraCells); }

    template<typename T>
    void resizeMatrix(T& m, int extraCells) const {
      m.resize(numAlns);
      for (unsigned i = 0; i < numAlns; ++i)
	m[i].resize(dpEnd(i) - dpBeg(i) + extraCells);
    }

    double probFromSpliceF(unsigned i, unsigned j, unsigned oldNumInplay,
			   unsigned& oldInplayPos) const;

    double probFromSpliceB(unsigned i, unsigned j, unsigned oldNumInplay,
			   unsigned& oldInplayPos) const;

    void calcBaseScores(unsigned i);
    void calcInsScores(unsigned i);
    void calcDelScores(unsigned i);
    void calcScoreMatrices();
    void initForwardBackward();
    void initDpBounds();

    double IB(unsigned i, unsigned j) const
    { return (alns[i].qstart == j) ? 1.0 : 0.0; }

    double IE(unsigned i, unsigned j) const
    { return (alns[i].qend == j) ? 1.0 : 0.0; }

    int JB(unsigned i, unsigned j) const
    { return (alns[i].qstart == j) ? 0 : INT_MIN/2; }

    long scoreIndel(unsigned i, unsigned j) const {
      return cell(Vmat, i, j) + cell(Dmat, i, j);
    }
};

}

#endif
