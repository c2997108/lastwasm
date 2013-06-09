CXX = g++
CC  = gcc

CXXFLAGS = -O3 -Wall -Wextra -Wcast-qual -Wswitch-enum -Wundef	\
-Wcast-align -Wno-long-long -ansi -pedantic
# -Wconversion
# -fomit-frame-pointer ?

CFLAGS = -Wall

DBSRC = Alphabet.cc MultiSequence.cc CyclicSubsetSeed.cc	\
SubsetSuffixArray.cc LastdbArguments.cc io.cc fileMap.cc	\
SubsetSuffixArraySort.cc MultiSequenceQual.cc lastdb.cc

DBINC = Alphabet.hh MultiSequence.hh CyclicSubsetSeed.hh		\
SubsetSuffixArray.hh LastdbArguments.hh io.hh fileMap.hh Mmap.hh	\
VectorOrMmap.hh ScoreMatrixRow.hh SequenceFormat.hh			\
qualityScoreUtil.hh stringify.hh version.hh

ALSRC = Alphabet.cc MultiSequence.cc CyclicSubsetSeed.cc		\
SubsetSuffixArray.cc LastalArguments.cc io.cc fileMap.cc		\
ScoreMatrix.cc DiagonalTable.cc SegmentPair.cc Alignment.cc		\
GappedXdropAligner.cc SegmentPairPot.cc AlignmentPot.cc			\
GeneralizedAffineGapCosts.cc Centroid.cc LambdaCalculator.cc		\
TwoQualityScoreMatrix.cc OneQualityScoreMatrix.cc QualityPssmMaker.cc	\
GeneticCode.cc gaplessXdrop.cc gaplessPssmXdrop.cc			\
gaplessTwoQualityXdrop.cc AlignmentWrite.cc MultiSequenceQual.cc	\
GappedXdropAlignerPssm.cc GappedXdropAligner2qual.cc			\
GappedXdropAligner3frame.cc lastal.cc

ALINC = Alphabet.hh MultiSequence.hh CyclicSubsetSeed.hh		\
SubsetSuffixArray.hh LastalArguments.hh io.hh fileMap.hh		\
ScoreMatrix.hh DiagonalTable.hh SegmentPair.hh Alignment.hh		\
GappedXdropAligner.hh SegmentPairPot.hh AlignmentPot.hh			\
GeneralizedAffineGapCosts.hh Centroid.hh LambdaCalculator.hh		\
TwoQualityScoreMatrix.hh OneQualityScoreMatrix.hh QualityPssmMaker.hh	\
GeneticCode.hh gaplessXdrop.hh gaplessPssmXdrop.hh			\
gaplessTwoQualityXdrop.hh Mmap.hh VectorOrMmap.hh ScoreMatrixRow.hh	\
SequenceFormat.hh GappedXdropAlignerInl.hh qualityScoreUtil.hh		\
stringify.hh version.hh

EXSRC = Alphabet.cc LastexArguments.cc io.cc ScoreMatrix.cc lastex.cc	\
gumbel_params/*.cpp

EXINC = Alphabet.hh LastexArguments.hh io.hh ScoreMatrix.hh	\
gumbel_params/*.hpp gumbel_params/*/*.h* version.hh

OBJ = lambda_calculator.o

all: lastdb lastal lastex

lastdb: $(DBSRC) $(DBINC) makefile
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ $(DBSRC)

lastal: $(ALSRC) $(ALINC) makefile $(OBJ)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ $(ALSRC) $(OBJ)

lastex: $(EXSRC) $(EXINC) makefile
	$(CXX) -Igumbel_params $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ $(EXSRC)

$(OBJ): CA_code/*.c CA_code/*.h makefile
	$(CC) $(CPPFLAGS) $(CFLAGS) -c CA_code/lambda_calculator.c

clean:
	rm -f lastdb lastal lastex $(OBJ)

VERSION = \"`hg id -n`\"
VERSION2 = \"`grep latesttagdistance ../.hg_archival.txt | cut -d' ' -f2`\"
VERSION3 = \"UNKNOWN\"

version.hh: FORCE
	if test -e ../.hg ; \
	then echo $(VERSION) | cmp -s $@ - || echo $(VERSION) > $@ ; \
	elif test -e ../.hg_archival.txt ; \
	then echo $(VERSION2) | cmp -s $@ - || echo $(VERSION2) > $@ ; \
	else test -e $@ || echo $(VERSION3) > $@ ; \
	fi

FORCE:
