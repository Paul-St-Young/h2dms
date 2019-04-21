#!/bin/bash
files="h2hp.cpp h2hp.h common.h h2wf.h h2ham.h h2ham.cpp j2pade.h j2pade.cpp"
for file in $files; do
  ln -s ../purecpp/$file .
done
