all:
	@cd src && $(MAKE)

prefix = /usr/local
exec_prefix = $(prefix)
bindir = $(exec_prefix)/bin
install: all
	mkdir -p $(bindir)
	cp src/last?? scripts/*.?? $(bindir)

clean:
	@cd src && $(MAKE) clean

VERSION = `svnversion .`

dist: log
	@cd src && $(MAKE) version.hh
	rsync -rC --exclude 'last??' d* e* makefile s* *.txt last-$(VERSION)
	zip -qrm archive/last-$(VERSION) last-$(VERSION)

log:
	svn log > ChangeLog.txt
