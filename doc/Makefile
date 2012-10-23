# Rules for converting the txt documents to html.

DOCS = bisulfite.html last-map-probs.html last-pair-probs.html	\
last-split-probs.html last-tutorial.html last.html

all: ${DOCS}

${DOCS}: last-doc.css Makefile

.SUFFIXES: .html .txt

RST_CSS = `locate html4css1.css | head -n1`
# Ugh!  Is there a better way?

.txt.html:
	rst2html --stylesheet-path ${RST_CSS},last-doc.css --initial-header-level=2 $< $@
