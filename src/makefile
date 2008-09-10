CXX = g++

CXXFLAGS = -Wall -W -Wshadow -Wcast-qual -Wconversion -Wwrite-strings	\
-Wno-long-long -ansi -pedantic -O3 -g
# use -W instead of -Wextra, so that it works on old versions of g++

DBSRC = Alphabet.cc MultiSequence.cc PeriodicSpacedSeed.cc	\
SuffixArray.cc LastdbArguments.cc io.cc lastdb.cc

DBINC = Alphabet.hh MultiSequence.hh PeriodicSpacedSeed.hh	\
SuffixArray.hh LastdbArguments.hh io.hh stringify.hh version.hh

ALSRC = Alphabet.cc MultiSequence.cc PeriodicSpacedSeed.cc	\
SuffixArray.cc LastalArguments.cc io.cc ScoreMatrix.cc		\
DiagonalTable.cc SegmentPair.cc Alignment.cc XdropAligner.cc	\
SegmentPairPot.cc AlignmentPot.cc GeneralizedAffineGapCosts.cc	\
lastal.cc

ALINC = Alphabet.hh MultiSequence.hh PeriodicSpacedSeed.hh	\
SuffixArray.hh LastalArguments.hh io.hh ScoreMatrix.hh		\
DiagonalTable.hh SegmentPair.hh Alignment.hh XdropAligner.hh	\
SegmentPairPot.hh AlignmentPot.hh GeneralizedAffineGapCosts.hh	\
stringify.hh version.hh

all: lastdb lastal

lastdb: $(DBSRC) $(DBINC) makefile
	$(CXX) $(CXXFLAGS) -o $@ $(DBSRC)

lastal: $(ALSRC) $(ALINC) makefile
	$(CXX) $(CXXFLAGS) -o $@ $(ALSRC)
