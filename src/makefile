CXX = g++
CC  = gcc

CXXFLAGS = -O3
# -fomit-frame-pointer ?
CFLAGS = -Wall

STRICT = -Wall -W -Wcast-qual -Wswitch-enum -Wundef	\
-Wcast-align -Wno-long-long -ansi -pedantic
# -Wconversion
# use -W instead of -Wextra, so that it works with old versions of g++

DBSRC = Alphabet.cc MultiSequence.cc CyclicSubsetSeed.cc	\
SubsetSuffixArray.cc LastdbArguments.cc io.cc fileMap.cc	\
SubsetSuffixArraySort.cc lastdb.cc

DBINC = Alphabet.hh MultiSequence.hh CyclicSubsetSeed.hh	\
SubsetSuffixArray.hh LastdbArguments.hh io.hh fileMap.hh	\
Mmap.hh VectorOrMmap.hh ScoreMatrixRow.hh stringify.hh version.hh

ALSRC = Alphabet.cc MultiSequence.cc CyclicSubsetSeed.cc	\
SubsetSuffixArray.cc LastalArguments.cc io.cc fileMap.cc	\
ScoreMatrix.cc DiagonalTable.cc SegmentPair.cc Alignment.cc	\
XdropAligner.cc SegmentPairPot.cc AlignmentPot.cc		\
GeneralizedAffineGapCosts.cc Centroid.cc LambdaCalculator.cc	\
OneQualityScoreMatrix.cc QualityPssmMaker.cc GeneticCode.cc	\
Xdrop3FrameAligner.cc AlignmentWrite.cc MultiSequenceQual.cc	\
SegmentPairPssm.cc XdropAlignerPssm.cc lastal.cc

ALINC = Alphabet.hh MultiSequence.hh CyclicSubsetSeed.hh	\
SubsetSuffixArray.hh LastalArguments.hh io.hh fileMap.hh	\
ScoreMatrix.hh DiagonalTable.hh SegmentPair.hh Alignment.hh	\
XdropAligner.hh SegmentPairPot.hh AlignmentPot.hh		\
GeneralizedAffineGapCosts.hh Centroid.hh LambdaCalculator.hh	\
OneQualityScoreMatrix.hh QualityPssmMaker.hh GeneticCode.hh	\
Xdrop3FrameAligner.hh Mmap.hh VectorOrMmap.hh ScoreMatrixRow.hh	\
SequenceFormat.hh qualityScoreUtil.hh stringify.hh version.hh

EXSRC = Alphabet.cc LastexArguments.cc io.cc ScoreMatrix.cc lastex.cc	\
gumbel_params/*.cpp

EXINC = Alphabet.hh LastexArguments.hh io.hh ScoreMatrix.hh	\
gumbel_params/*.hpp gumbel_params/*/*.h* version.hh

OBJ = lambda_calculator.o

all: lastdb lastal lastex

lastdb: $(DBSRC) $(DBINC) makefile
	$(CXX) $(STRICT) $(CXXFLAGS) -o $@ $(DBSRC)

lastal: $(ALSRC) $(ALINC) makefile $(OBJ)
	$(CXX) $(STRICT) $(CXXFLAGS) -o $@ $(ALSRC) $(OBJ)

lastex: $(EXSRC) $(EXINC) makefile
	$(CXX) -Igumbel_params $(CXXFLAGS) -o $@ $(EXSRC)

$(OBJ): CA_code/*.c CA_code/*.h makefile
	$(CC) $(CFLAGS) -c CA_code/lambda_calculator.c
