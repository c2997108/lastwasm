# Rules for converting the txt documents to html.

DOCS = bisulfite.html last-map-probs.html last-pair-probs.html	\
last-split-probs.html last-tutorial.html

all: ${DOCS}

${DOCS}: last-doc.css Makefile

.SUFFIXES: .html .txt

RST_CSS = /usr/share/docutils/writers/html4css1/html4css1.css

.txt.html:
	rst2html --stylesheet-path ${RST_CSS},last-doc.css --initial-header-level=2 $< $@
