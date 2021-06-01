#!/bin/sh

if ! [ -e ../src/compairr ] ; then
    echo The compairr binary is missing
    echo Test failed.
    exit 1
fi

../src/compairr -m seta.tsv setb.tsv -d 1 -i -l compairr.log -o output.tsv

if diff -q output.tsv expected.tsv; then
    echo Test completed successfully.
else
    echo Test failed.
    exit 1
fi
