CXXFLAGS = -msse4 -O3 -Wall -g -std=c++11 -pthread
# -fomit-frame-pointer ?

# The number of symbol types that we can keep in sequences (e.g. 20
# amino acids, plus ambiguous ones, in upper & lower case, plus one
# delimiter):
ALPHABET_CAPACITY = 66
CPPF = -DALPHABET_CAPACITY=$(ALPHABET_CAPACITY) -DHAS_CXX_THREADS $(CPPFLAGS)

CFLAGS = -Wall -O2

alpObj = alp/sls_alignment_evaluer.o alp/sls_pvalues.o		\
alp/sls_alp_sim.o alp/sls_alp_regression.o alp/sls_alp_data.o	\
alp/sls_alp.o alp/sls_basic.o alp/njn_localmaxstatmatrix.o	\
alp/njn_localmaxstat.o alp/njn_localmaxstatutil.o		\
alp/njn_dynprogprob.o alp/njn_dynprogprobproto.o		\
alp/njn_dynprogproblim.o alp/njn_ioutil.o alp/njn_random.o	\
alp/sls_falp_alignment_evaluer.o alp/sls_fsa1_pvalues.o		\
alp/sls_fsa1_utils.o alp/sls_fsa1.o alp/sls_fsa1_parameters.o

indexObj = Alphabet.o CyclicSubsetSeed.o LambdaCalculator.o		\
MultiSequence.o MultiSequenceQual.o ScoreMatrix.o			\
SubsetMinimizerFinder.o SubsetSuffixArray.o SubsetSuffixArraySort.o	\
TantanMasker.o dna_words_finder.o fileMap.o cbrc_linalg.o		\
mcf_substitution_matrix_stats.o tantan.o LastdbArguments.o lastdb.o

alignObj = Alphabet.o Centroid.o CyclicSubsetSeed.o			\
LambdaCalculator.o MultiSequence.o MultiSequenceQual.o ScoreMatrix.o	\
SubsetMinimizerFinder.o SubsetSuffixArray.o SubsetSuffixArraySearch.o	\
TantanMasker.o dna_words_finder.o fileMap.o tantan.o			\
LastalArguments.o lastal.o Alignment.o AlignmentPot.o			\
AlignmentWrite.o GappedXdropAligner.o GappedXdropAlignerDna.o		\
GappedXdropAlignerPssm.o GappedXdropAligner2qual.o			\
GappedXdropAligner3frame.o GappedXdropAlignerFrame.o			\
mcf_alignment_path_adder.o mcf_frameshift_xdrop_aligner.o		\
mcf_gap_costs.o GeneticCode.o GreedyXdropAligner.o LastEvaluer.o	\
OneQualityScoreMatrix.o QualityPssmMaker.o SegmentPair.o		\
SegmentPairPot.o TwoQualityScoreMatrix.o cbrc_linalg.o			\
mcf_substitution_matrix_stats.o split/cbrc_split_aligner.o		\
split/cbrc_unsplit_alignment.o split/last_split_options.o		\
split/mcf_last_splitter.o $(alpObj)

splitObj = Alphabet.o LambdaCalculator.o MultiSequence.o fileMap.o	\
cbrc_linalg.o mcf_substitution_matrix_stats.o				\
split/cbrc_unsplit_alignment.o split/last_split_options.o		\
split/last-split-main.o split/cbrc_split_aligner.o			\
split/mcf_last_splitter.o split/last-split.o

PPOBJ = last-pair-probs.o last-pair-probs-main.o

MBOBJ = last-merge-batches.o

ALL = ../bin/lastdb ../bin/lastal ../bin/last-split	\
../bin/last-merge-batches ../bin/last-pair-probs

all: $(ALL)

../bin/lastdb: $(indexObj)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $(indexObj) -lz

../bin/lastal: $(alignObj)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $(alignObj) -lz

../bin/last-split: $(splitObj)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $(splitObj)

../bin/last-pair-probs: $(PPOBJ)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $(PPOBJ) -lz

../bin/last-merge-batches: $(MBOBJ)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $(MBOBJ)

.SUFFIXES:
.SUFFIXES: .o .c .cc .cpp

.c.o:
	$(CC) $(CPPFLAGS) $(CFLAGS) -c -o $@ $<

.cc.o:
	$(CXX) $(CPPF) $(CXXFLAGS) -I. -c -o $@ $<

.cpp.o:
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c -o $@ $<

clean:
	rm -f $(ALL) *.o* */*.o*

CyclicSubsetSeedData.hh: ../data/*.seed
	../build/seed-inc.sh ../data/*.seed > $@

GeneticCodeData.hh: ../data/gc.prt
	../build/gc-inc.sh ../data/gc.prt > $@

ScoreMatrixData.hh: ../data/*.mat
	../build/mat-inc.sh ../data/*.mat > $@

VERSION1 = git describe --dirty
VERSION2 = echo '$Format:%d$ ' | sed -e 's/.*tag: *//' -e 's/[,) ].*//'

VERSION = \"`test -e ../.git && $(VERSION1) || $(VERSION2)`\"

version.hh: FORCE
	echo $(VERSION) | cmp -s $@ - || echo $(VERSION) > $@

FORCE:

depend:
	sed '/[m][v]/q' makefile > m
	$(CXX) -MM -I. -std=c++11 *.cc >> m
	$(CC) -MM *.c >> m
	$(CXX) -MM alp/*.cpp | sed 's|.*:|alp/&|' >> m
	$(CXX) -MM -I. split/*.cc | sed 's|.*:|split/&|' >> m
	mv m makefile
Alignment.o: Alignment.cc Alignment.hh Centroid.hh GappedXdropAligner.hh \
 mcf_big_seq.hh mcf_contiguous_queue.hh mcf_reverse_queue.hh \
 mcf_gap_costs.hh mcf_simd.hh ScoreMatrixRow.hh OneQualityScoreMatrix.hh \
 mcf_substitution_matrix_stats.hh GreedyXdropAligner.hh SegmentPair.hh \
 mcf_frameshift_xdrop_aligner.hh Alphabet.hh GeneticCode.hh \
 TwoQualityScoreMatrix.hh
AlignmentPot.o: AlignmentPot.cc AlignmentPot.hh Alignment.hh Centroid.hh \
 GappedXdropAligner.hh mcf_big_seq.hh mcf_contiguous_queue.hh \
 mcf_reverse_queue.hh mcf_gap_costs.hh mcf_simd.hh ScoreMatrixRow.hh \
 OneQualityScoreMatrix.hh mcf_substitution_matrix_stats.hh \
 GreedyXdropAligner.hh SegmentPair.hh mcf_frameshift_xdrop_aligner.hh
AlignmentWrite.o: AlignmentWrite.cc Alignment.hh Centroid.hh \
 GappedXdropAligner.hh mcf_big_seq.hh mcf_contiguous_queue.hh \
 mcf_reverse_queue.hh mcf_gap_costs.hh mcf_simd.hh ScoreMatrixRow.hh \
 OneQualityScoreMatrix.hh mcf_substitution_matrix_stats.hh \
 GreedyXdropAligner.hh SegmentPair.hh mcf_frameshift_xdrop_aligner.hh \
 GeneticCode.hh LastEvaluer.hh alp/sls_alignment_evaluer.hpp \
 alp/sls_pvalues.hpp alp/sls_basic.hpp MultiSequence.hh VectorOrMmap.hh \
 Mmap.hh fileMap.hh stringify.hh Alphabet.hh
Alphabet.o: Alphabet.cc Alphabet.hh mcf_big_seq.hh
cbrc_linalg.o: cbrc_linalg.cc cbrc_linalg.hh
Centroid.o: Centroid.cc Centroid.hh GappedXdropAligner.hh mcf_big_seq.hh \
 mcf_contiguous_queue.hh mcf_reverse_queue.hh mcf_gap_costs.hh \
 mcf_simd.hh ScoreMatrixRow.hh OneQualityScoreMatrix.hh \
 mcf_substitution_matrix_stats.hh GappedXdropAlignerInl.hh
CyclicSubsetSeed.o: CyclicSubsetSeed.cc CyclicSubsetSeed.hh \
 CyclicSubsetSeedData.hh zio.hh mcf_zstream.hh stringify.hh
dna_words_finder.o: dna_words_finder.cc dna_words_finder.hh
fileMap.o: fileMap.cc fileMap.hh stringify.hh
GappedXdropAligner2qual.o: GappedXdropAligner2qual.cc \
 GappedXdropAligner.hh mcf_big_seq.hh mcf_contiguous_queue.hh \
 mcf_reverse_queue.hh mcf_gap_costs.hh mcf_simd.hh ScoreMatrixRow.hh \
 GappedXdropAlignerInl.hh TwoQualityScoreMatrix.hh \
 mcf_substitution_matrix_stats.hh
GappedXdropAligner3frame.o: GappedXdropAligner3frame.cc \
 GappedXdropAligner.hh mcf_big_seq.hh mcf_contiguous_queue.hh \
 mcf_reverse_queue.hh mcf_gap_costs.hh mcf_simd.hh ScoreMatrixRow.hh \
 GappedXdropAlignerInl.hh
GappedXdropAligner3framePssm.o: GappedXdropAligner3framePssm.cc \
 GappedXdropAligner.hh mcf_big_seq.hh mcf_contiguous_queue.hh \
 mcf_reverse_queue.hh mcf_gap_costs.hh mcf_simd.hh ScoreMatrixRow.hh \
 GappedXdropAlignerInl.hh
GappedXdropAligner.o: GappedXdropAligner.cc GappedXdropAligner.hh \
 mcf_big_seq.hh mcf_contiguous_queue.hh mcf_reverse_queue.hh \
 mcf_gap_costs.hh mcf_simd.hh ScoreMatrixRow.hh GappedXdropAlignerInl.hh
GappedXdropAlignerDna.o: GappedXdropAlignerDna.cc GappedXdropAligner.hh \
 mcf_big_seq.hh mcf_contiguous_queue.hh mcf_reverse_queue.hh \
 mcf_gap_costs.hh mcf_simd.hh ScoreMatrixRow.hh GappedXdropAlignerInl.hh
GappedXdropAlignerFrame.o: GappedXdropAlignerFrame.cc \
 GappedXdropAligner.hh mcf_big_seq.hh mcf_contiguous_queue.hh \
 mcf_reverse_queue.hh mcf_gap_costs.hh mcf_simd.hh ScoreMatrixRow.hh \
 GappedXdropAlignerInl.hh
GappedXdropAlignerPssm.o: GappedXdropAlignerPssm.cc GappedXdropAligner.hh \
 mcf_big_seq.hh mcf_contiguous_queue.hh mcf_reverse_queue.hh \
 mcf_gap_costs.hh mcf_simd.hh ScoreMatrixRow.hh GappedXdropAlignerInl.hh
GeneticCode.o: GeneticCode.cc GeneticCode.hh GeneticCodeData.hh \
 Alphabet.hh mcf_big_seq.hh zio.hh mcf_zstream.hh
GreedyXdropAligner.o: GreedyXdropAligner.cc GreedyXdropAligner.hh \
 ScoreMatrixRow.hh
LambdaCalculator.o: LambdaCalculator.cc LambdaCalculator.hh \
 cbrc_linalg.hh
LastalArguments.o: LastalArguments.cc LastalArguments.hh \
 SequenceFormat.hh split/last_split_options.hh stringify.hh getoptUtil.hh \
 version.hh
lastal.o: lastal.cc last.hh Alphabet.hh mcf_big_seq.hh \
 CyclicSubsetSeed.hh MultiSequence.hh ScoreMatrixRow.hh VectorOrMmap.hh \
 Mmap.hh fileMap.hh stringify.hh SequenceFormat.hh \
 SubsetMinimizerFinder.hh SubsetSuffixArray.hh dna_words_finder.hh \
 mcf_packed_array.hh qualityScoreUtil.hh LastalArguments.hh \
 split/last_split_options.hh QualityPssmMaker.hh OneQualityScoreMatrix.hh \
 mcf_substitution_matrix_stats.hh TwoQualityScoreMatrix.hh LastEvaluer.hh \
 mcf_frameshift_xdrop_aligner.hh mcf_gap_costs.hh \
 alp/sls_alignment_evaluer.hpp alp/sls_pvalues.hpp alp/sls_basic.hpp \
 GeneticCode.hh AlignmentPot.hh Alignment.hh Centroid.hh \
 GappedXdropAligner.hh mcf_contiguous_queue.hh mcf_reverse_queue.hh \
 mcf_simd.hh GreedyXdropAligner.hh SegmentPair.hh SegmentPairPot.hh \
 ScoreMatrix.hh TantanMasker.hh tantan.hh DiagonalTable.hh \
 gaplessXdrop.hh gaplessPssmXdrop.hh gaplessTwoQualityXdrop.hh zio.hh \
 mcf_zstream.hh threadUtil.hh split/mcf_last_splitter.hh \
 split/cbrc_split_aligner.hh split/cbrc_unsplit_alignment.hh \
 split/cbrc_int_exponentiator.hh Alphabet.hh MultiSequence.hh \
 split/last_split_options.hh version.hh
LastdbArguments.o: LastdbArguments.cc LastdbArguments.hh \
 SequenceFormat.hh stringify.hh getoptUtil.hh version.hh
lastdb.o: lastdb.cc last.hh Alphabet.hh mcf_big_seq.hh \
 CyclicSubsetSeed.hh MultiSequence.hh ScoreMatrixRow.hh VectorOrMmap.hh \
 Mmap.hh fileMap.hh stringify.hh SequenceFormat.hh \
 SubsetMinimizerFinder.hh SubsetSuffixArray.hh dna_words_finder.hh \
 mcf_packed_array.hh qualityScoreUtil.hh LastdbArguments.hh \
 TantanMasker.hh tantan.hh zio.hh mcf_zstream.hh threadUtil.hh version.hh
LastEvaluer.o: LastEvaluer.cc LastEvaluer.hh ScoreMatrixRow.hh \
 mcf_frameshift_xdrop_aligner.hh mcf_gap_costs.hh \
 alp/sls_alignment_evaluer.hpp alp/sls_pvalues.hpp alp/sls_basic.hpp \
 GeneticCode.hh mcf_alignment_path_adder.hh \
 alp/sls_falp_alignment_evaluer.hpp alp/sls_fsa1_pvalues.hpp \
 LastEvaluerData.hh
last-pair-probs.o: last-pair-probs.cc last-pair-probs.hh zio.hh \
 mcf_zstream.hh stringify.hh
last-pair-probs-main.o: last-pair-probs-main.cc last-pair-probs.hh \
 stringify.hh version.hh
mcf_alignment_path_adder.o: mcf_alignment_path_adder.cc \
 mcf_alignment_path_adder.hh
mcf_frameshift_xdrop_aligner.o: mcf_frameshift_xdrop_aligner.cc \
 mcf_frameshift_xdrop_aligner.hh mcf_gap_costs.hh
mcf_gap_costs.o: mcf_gap_costs.cc mcf_gap_costs.hh
mcf_substitution_matrix_stats.o: mcf_substitution_matrix_stats.cc \
 mcf_substitution_matrix_stats.hh LambdaCalculator.hh cbrc_linalg.hh
MultiSequence.o: MultiSequence.cc MultiSequence.hh mcf_big_seq.hh \
 ScoreMatrixRow.hh VectorOrMmap.hh Mmap.hh fileMap.hh stringify.hh io.hh
MultiSequenceQual.o: MultiSequenceQual.cc MultiSequence.hh mcf_big_seq.hh \
 ScoreMatrixRow.hh VectorOrMmap.hh Mmap.hh fileMap.hh stringify.hh
OneQualityScoreMatrix.o: OneQualityScoreMatrix.cc \
 OneQualityScoreMatrix.hh mcf_substitution_matrix_stats.hh \
 ScoreMatrixRow.hh qualityScoreUtil.hh stringify.hh
QualityPssmMaker.o: QualityPssmMaker.cc QualityPssmMaker.hh \
 ScoreMatrixRow.hh qualityScoreUtil.hh stringify.hh
ScoreMatrix.o: ScoreMatrix.cc ScoreMatrix.hh ScoreMatrixRow.hh \
 ScoreMatrixData.hh qualityScoreUtil.hh stringify.hh zio.hh \
 mcf_zstream.hh
SegmentPair.o: SegmentPair.cc SegmentPair.hh mcf_big_seq.hh
SegmentPairPot.o: SegmentPairPot.cc SegmentPairPot.hh SegmentPair.hh \
 mcf_big_seq.hh
SubsetMinimizerFinder.o: SubsetMinimizerFinder.cc \
 SubsetMinimizerFinder.hh CyclicSubsetSeed.hh
SubsetSuffixArray.o: SubsetSuffixArray.cc SubsetSuffixArray.hh \
 CyclicSubsetSeed.hh dna_words_finder.hh mcf_big_seq.hh \
 mcf_packed_array.hh VectorOrMmap.hh Mmap.hh fileMap.hh stringify.hh \
 io.hh
SubsetSuffixArraySearch.o: SubsetSuffixArraySearch.cc \
 SubsetSuffixArray.hh CyclicSubsetSeed.hh dna_words_finder.hh \
 mcf_big_seq.hh mcf_packed_array.hh VectorOrMmap.hh Mmap.hh fileMap.hh \
 stringify.hh
SubsetSuffixArraySort.o: SubsetSuffixArraySort.cc SubsetSuffixArray.hh \
 CyclicSubsetSeed.hh dna_words_finder.hh mcf_big_seq.hh \
 mcf_packed_array.hh VectorOrMmap.hh Mmap.hh fileMap.hh stringify.hh
tantan.o: tantan.cc tantan.hh mcf_simd.hh
TantanMasker.o: TantanMasker.cc TantanMasker.hh ScoreMatrixRow.hh \
 tantan.hh ScoreMatrix.hh mcf_substitution_matrix_stats.hh
TwoQualityScoreMatrix.o: TwoQualityScoreMatrix.cc \
 TwoQualityScoreMatrix.hh mcf_substitution_matrix_stats.hh \
 ScoreMatrixRow.hh qualityScoreUtil.hh stringify.hh
last-merge-batches.o: last-merge-batches.c version.hh
alp/njn_dynprogprob.o: alp/njn_dynprogprob.cpp alp/njn_dynprogprob.hpp \
 alp/njn_dynprogprobproto.hpp alp/njn_memutil.hpp alp/njn_ioutil.hpp
alp/njn_dynprogproblim.o: alp/njn_dynprogproblim.cpp \
 alp/njn_dynprogproblim.hpp alp/njn_dynprogprob.hpp \
 alp/njn_dynprogprobproto.hpp alp/njn_memutil.hpp alp/njn_ioutil.hpp
alp/njn_dynprogprobproto.o: alp/njn_dynprogprobproto.cpp \
 alp/njn_dynprogprobproto.hpp
alp/njn_ioutil.o: alp/njn_ioutil.cpp alp/njn_ioutil.hpp
alp/njn_localmaxstat.o: alp/njn_localmaxstat.cpp alp/sls_basic.hpp \
 alp/njn_localmaxstat.hpp alp/njn_memutil.hpp alp/njn_ioutil.hpp \
 alp/njn_dynprogproblim.hpp alp/njn_dynprogprob.hpp \
 alp/njn_dynprogprobproto.hpp alp/njn_function.hpp alp/njn_doubletype.hpp \
 alp/njn_integer.hpp alp/njn_localmaxstatutil.hpp alp/njn_matrix.hpp \
 alp/njn_approx.hpp alp/njn_vector.hpp
alp/njn_localmaxstatmatrix.o: alp/njn_localmaxstatmatrix.cpp \
 alp/njn_localmaxstatmatrix.hpp alp/njn_localmaxstat.hpp \
 alp/njn_localmaxstatutil.hpp alp/njn_matrix.hpp alp/njn_approx.hpp \
 alp/njn_doubletype.hpp alp/njn_ioutil.hpp alp/njn_vector.hpp \
 alp/njn_memutil.hpp
alp/njn_localmaxstatutil.o: alp/njn_localmaxstatutil.cpp \
 alp/njn_localmaxstatutil.hpp alp/njn_matrix.hpp alp/njn_approx.hpp \
 alp/njn_doubletype.hpp alp/njn_ioutil.hpp alp/njn_vector.hpp \
 alp/njn_dynprogproblim.hpp alp/njn_dynprogprob.hpp \
 alp/njn_dynprogprobproto.hpp alp/njn_integer.hpp alp/njn_memutil.hpp \
 alp/njn_root.hpp alp/njn_function.hpp alp/sls_basic.hpp
alp/njn_random.o: alp/njn_random.cpp alp/njn_random.hpp
alp/sls_alignment_evaluer.o: alp/sls_alignment_evaluer.cpp \
 alp/sls_alignment_evaluer.hpp alp/sls_pvalues.hpp alp/sls_basic.hpp \
 alp/sls_alp.hpp alp/sls_alp_data.hpp alp/sls_alp_regression.hpp \
 alp/njn_random.hpp alp/njn_uniform.hpp alp/sls_alp_sim.hpp \
 alp/njn_localmaxstatmatrix.hpp alp/njn_localmaxstat.hpp \
 alp/njn_localmaxstatutil.hpp alp/njn_matrix.hpp alp/njn_approx.hpp \
 alp/njn_doubletype.hpp alp/njn_ioutil.hpp alp/njn_vector.hpp
alp/sls_alp.o: alp/sls_alp.cpp alp/sls_alp.hpp alp/sls_alp_data.hpp \
 alp/sls_basic.hpp alp/sls_alp_regression.hpp alp/njn_random.hpp \
 alp/njn_uniform.hpp
alp/sls_alp_data.o: alp/sls_alp_data.cpp alp/sls_alp_data.hpp \
 alp/sls_basic.hpp alp/sls_alp_regression.hpp alp/njn_random.hpp \
 alp/njn_uniform.hpp
alp/sls_alp_regression.o: alp/sls_alp_regression.cpp \
 alp/sls_alp_regression.hpp alp/sls_basic.hpp
alp/sls_alp_sim.o: alp/sls_alp_sim.cpp alp/sls_alp_sim.hpp \
 alp/sls_alp_data.hpp alp/sls_basic.hpp alp/sls_alp_regression.hpp \
 alp/njn_random.hpp alp/njn_uniform.hpp alp/sls_alp.hpp
alp/sls_basic.o: alp/sls_basic.cpp alp/sls_basic.hpp
alp/sls_falp_alignment_evaluer.o: alp/sls_falp_alignment_evaluer.cpp \
 alp/sls_falp_alignment_evaluer.hpp alp/sls_fsa1_pvalues.hpp \
 alp/sls_basic.hpp alp/sls_fsa1_parameters.hpp alp/sls_fsa1_utils.hpp \
 alp/njn_random.hpp alp/njn_uniform.hpp alp/sls_alp_regression.hpp \
 alp/njn_localmaxstatmatrix.hpp alp/njn_localmaxstat.hpp \
 alp/njn_localmaxstatutil.hpp alp/njn_matrix.hpp alp/njn_approx.hpp \
 alp/njn_doubletype.hpp alp/njn_ioutil.hpp alp/njn_vector.hpp \
 alp/sls_fsa1.hpp
alp/sls_fsa1.o: alp/sls_fsa1.cpp alp/sls_fsa1.hpp alp/sls_alp_regression.hpp \
 alp/sls_basic.hpp alp/sls_fsa1_utils.hpp alp/njn_random.hpp \
 alp/njn_uniform.hpp alp/sls_fsa1_parameters.hpp alp/sls_fsa1_pvalues.hpp \
 alp/njn_localmaxstatmatrix.hpp alp/njn_localmaxstat.hpp \
 alp/njn_localmaxstatutil.hpp alp/njn_matrix.hpp alp/njn_approx.hpp \
 alp/njn_doubletype.hpp alp/njn_ioutil.hpp alp/njn_vector.hpp
alp/sls_fsa1_parameters.o: alp/sls_fsa1_parameters.cpp \
 alp/sls_fsa1_parameters.hpp alp/sls_fsa1_utils.hpp alp/sls_basic.hpp \
 alp/njn_random.hpp alp/njn_uniform.hpp alp/sls_alp_regression.hpp \
 alp/sls_fsa1_pvalues.hpp
alp/sls_fsa1_pvalues.o: alp/sls_fsa1_pvalues.cpp alp/sls_fsa1_pvalues.hpp \
 alp/sls_basic.hpp alp/sls_fsa1_utils.hpp alp/njn_random.hpp \
 alp/njn_uniform.hpp alp/sls_normal_distr_array.hpp
alp/sls_fsa1_utils.o: alp/sls_fsa1_utils.cpp alp/sls_fsa1_utils.hpp \
 alp/sls_basic.hpp alp/njn_random.hpp alp/njn_uniform.hpp
alp/sls_pvalues.o: alp/sls_pvalues.cpp alp/sls_pvalues.hpp alp/sls_basic.hpp \
 alp/sls_alp_data.hpp alp/sls_alp_regression.hpp alp/njn_random.hpp \
 alp/njn_uniform.hpp alp/sls_normal_distr_array.hpp
split/cbrc_split_aligner.o: split/cbrc_split_aligner.cc \
 split/cbrc_split_aligner.hh split/cbrc_unsplit_alignment.hh \
 split/cbrc_int_exponentiator.hh Alphabet.hh mcf_big_seq.hh \
 MultiSequence.hh ScoreMatrixRow.hh VectorOrMmap.hh Mmap.hh fileMap.hh \
 stringify.hh mcf_substitution_matrix_stats.hh
split/cbrc_unsplit_alignment.o: split/cbrc_unsplit_alignment.cc \
 split/cbrc_unsplit_alignment.hh
split/last-split.o: split/last-split.cc split/last-split.hh \
 split/last_split_options.hh split/mcf_last_splitter.hh \
 split/cbrc_split_aligner.hh split/cbrc_unsplit_alignment.hh \
 split/cbrc_int_exponentiator.hh Alphabet.hh mcf_big_seq.hh \
 MultiSequence.hh ScoreMatrixRow.hh VectorOrMmap.hh Mmap.hh fileMap.hh \
 stringify.hh
split/last-split-main.o: split/last-split-main.cc split/last-split.hh \
 split/last_split_options.hh stringify.hh version.hh
split/last_split_options.o: split/last_split_options.cc \
 split/last_split_options.hh
split/mcf_last_splitter.o: split/mcf_last_splitter.cc \
 split/mcf_last_splitter.hh split/cbrc_split_aligner.hh \
 split/cbrc_unsplit_alignment.hh split/cbrc_int_exponentiator.hh \
 Alphabet.hh mcf_big_seq.hh MultiSequence.hh ScoreMatrixRow.hh \
 VectorOrMmap.hh Mmap.hh fileMap.hh stringify.hh \
 split/last_split_options.hh
