# Makefile

PREFIX=/usr/local

all : compairr

compairr:
	make -C src compairr

test: compairr
	make -C test

install: compairr test
	/usr/bin/install -c src/compairr $(PREFIX)/bin

clean:
	make -C src clean
