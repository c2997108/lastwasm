all: version log
	@cd src && $(MAKE)

dist: version log
	rsync -rC --exclude 'last??' doc examples s* *.txt last-`svnversion .`
	zip -qrm archive/last-`svnversion .` last-`svnversion .`

version: FORCE
	echo \"`svnversion .`\" > src/version.hh

log: FORCE
	svn log > ChangeLog.txt

FORCE:
