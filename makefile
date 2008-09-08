all: dist
	@cd src && $(MAKE)

dist: version log
	mkdir -p last-`svnversion`/src
	mkdir -p last-`svnversion`/doc
	mkdir -p last-`svnversion`/scripts
	ln src/*.hh src/*.cc src/makefile last-`svnversion`/src
	ln doc/*.txt doc/HOXD70 last-`svnversion`/doc
	ln scripts/*.sh scripts/*.py last-`svnversion`/scripts
	ln *.txt last-`svnversion`
	zip -qrm archive/last-`svnversion` last-`svnversion`

version: FORCE
	echo \"`svnversion`\" > src/version.hh

log: FORCE
	svn log > ChangeLog.txt

FORCE:
