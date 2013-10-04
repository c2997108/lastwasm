CXXFLAGS = -O3
all:
	@cd src && $(MAKE) CXXFLAGS="$(CXXFLAGS)"

prefix = /usr/local
exec_prefix = $(prefix)
bindir = $(exec_prefix)/bin
install: all
	mkdir -p $(bindir)
	cp src/last?? src/last-split scripts/*.?? $(bindir)

clean:
	@cd src && $(MAKE) clean

html:
	@cd doc && $(MAKE)

distdir = last-`hg id -n`

dist: log html
	@cd src && $(MAKE) version.hh
	rsync -rC --exclude 'last??' doc examples makefile s* *.txt $(distdir)
	zip -qrm $(distdir) $(distdir)

log:
	hg log --style changelog > ChangeLog.txt
