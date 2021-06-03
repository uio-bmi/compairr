# Makefile for CompAIRR

ifndef PREFIX
	PREFIX=/usr/local
endif

all : compairr

compairr:
	make -C src compairr

test: compairr
	make -C test

install: compairr test
	/usr/bin/install -c src/compairr $(PREFIX)/bin/compairr

clean:
	make -C src clean
