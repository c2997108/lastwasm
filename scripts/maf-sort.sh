#! /bin/sh

# Sort MAF-format alignments by sequence name, then strand, then start
# position, then end position, of the top sequence.  Also, merge
# identical alignments.  Comment lines starting with "#" are written
# at the top, in unchanged order.  If option "-d" is specified, then
# alignments that appear only once are omitted (like uniq -d).

# XXX Preceding whitespace is considered part of the sequence name.  I
# want to use sort -b, but it seems to be broken in different ways for
# different versions of sort!

# XXX Alignments with differences in whitespace are considered
# non-identical.

# This script uses perl instead of specialized commands like uniq.
# The reason is that, on some systems (e.g. Mac OS X), uniq doesn't
# work with long lines.

# Make "sort" use a standard ordering:
LC_ALL=C
export LC_ALL

uniqOpt=1
while getopts hd opt
do
    case $opt in
	h)  cat <<EOF
Usage: $(basename $0) [options] my-alignments.maf

Options:
  -h  show this help message and exit
  -d  only print duplicate alignments
EOF
	    exit
	    ;;
	d)  uniqOpt=2
            ;;
    esac
done
shift $((OPTIND - 1))

tmpfile=${TMPDIR-/tmp}/maf-sort.$$

cat "$@" | tee $tmpfile | perl -ne 'print if /^#/'

perl -ne 'print unless /^#/' $tmpfile |  # remove comment lines
perl -pe 'y/ /!/  if /^a/'     |  # change spaces to '!'s in 'a' lines
perl -pe 's/\n/#/ if /\S/'     |  # join each alignment into one big line
sort -k2,2 -k5,5 -k3,3n -k4,4n |  # sort the lines

# print only the first (or second) of each run of identical lines:
perl -ne '$c = 0 if $x ne $_; $x = $_; print if ++$c == '$uniqOpt |

perl -pe 's/#/\n/g'            |  # undo the line-joining
perl -pe 'y/!/ / if /^a/'         # change '!'s back to spaces in 'a' lines

rm $tmpfile
