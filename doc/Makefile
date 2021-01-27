# Rules for converting the txt documents to html.

DOCS = FAQ.html bisulfite.html last-dotplot.html last-evalues.html	\
last-map-probs.html last-matrices.html last-pair-probs.html		\
last-papers.html last-parallel.html last-postmask.html			\
last-repeats.html last-seeds.html last-split.html last-train.html	\
last-tuning.html last-tutorial.html last.html lastal.html lastdb.html	\
maf-convert.html maf-cut.html

all: ${DOCS}

${DOCS}: last-doc.css Makefile

.SUFFIXES: .html .txt

RST2HTML = rst2html  # xxx .py ???

# Ugh!  Is there a better way?
RST_CSS = `locate html4css1.css | tail -n1`
#RST_CSS = html4css1.css

RSTFLAGS = --initial-header-level=2 --no-compact-lists	\
--no-compact-field-lists --option-limit=0 --no-doc-info

.txt.html:
	${RST2HTML} --stylesheet-path=${RST_CSS},last-doc.css ${RSTFLAGS} $< $@

last-matrices.txt: ../data/*.mat
	../build/mat-doc.sh ../data/*.mat > $@

last-seeds.txt: ../data/*.seed
	cd ../data && ../build/seed-doc.sh [!R]*d RY?-*d RY??-*d > ../doc/$@
