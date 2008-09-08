#! /bin/sh

# Sort MAF-format alignments by sequence name, then start position,
# then end position, of the top sequence.  Also, merge identical
# alignments.  Comment lines starting with "#" are written at the top,
# in unchanged order.  If option "-d" is specified, then alignments
# that appear only once are omitted (like uniq -d).

opt=
if [ "$1" = "-d" ]
then
    opt=$1
    shift
fi

cat "$@" > /tmp/$$  # can we avoid using a temporary file?

grep '^#' /tmp/$$

grep -v '^#' /tmp/$$       |  # remove comment lines
sed '/^a/y/ /!/'           |  # change spaces to '!'s in 'a' lines
perl -pe 's/\n/#/ if /\S/' |  # join each alignment into one big line
sort -k2,2 -k3,3n -k4,4n   |  # sort the lines
uniq $opt                  |  # merge identical alignments (don't use sort -u)
perl -pe 's/#/\n/g'        |  # undo the line-joining
sed '/^a/y/!/ /'              # change '!'s back to spaces in 'a' lines

rm /tmp/$$
