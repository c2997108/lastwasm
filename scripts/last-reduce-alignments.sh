#! /bin/sh

# This script reads MAF-format alignments with lastal header
# information, removes "uninteresting" alignments, and writes the
# remainder.  Specifically, it removes alignments that are "dominated"
# in both genomes.  See last-remove-dominated.py for the definition of
# "dominated".  This procedure is likely to remove many paralogs, and
# it is unlikely to remove one-to-one orthologs.

# If option "-d" is specified, then it removes alignments that are
# dominated in either genome.

opt=
if [ "$1" = "-d" ]
then
    opt=$1
    shift
fi

dir=`dirname $0`  # assume the other scripts are in the same directory

{
    # remove alignments that are dominated in the upper sequence
    $dir/maf-sort.sh "$@" |
    $dir/last-remove-dominated.py

    # remove alignments that are dominated in the lower sequence
    $dir/maf-swap.py "$@" |
    $dir/maf-sort.sh |
    $dir/last-remove-dominated.py |
    $dir/maf-swap.py

} | $dir/maf-sort.sh $opt  # merge the results
