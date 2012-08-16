CXXFLAGS = -O3
all:
	@cd src && $(MAKE) CXXFLAGS="$(CXXFLAGS)"

prefix = /usr/local
exec_prefix = $(prefix)
bindir = $(exec_prefix)/bin
install: all
	mkdir -p $(bindir)
	cp src/last?? scripts/*.?? $(bindir)

clean:
	@cd src && $(MAKE) clean

distdir = last-`svnversion .`

dist: log
	@cd src && $(MAKE) version.hh
	rsync -rC --exclude 'last??' doc examples makefile s* *.txt $(distdir)
	zip -qrm archive/$(distdir) $(distdir)

log:
	svn log > ChangeLog.txt
