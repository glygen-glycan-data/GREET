#!/bin/sh

for C in run/*.ini; do
  BASE=`basename "$C" .ini`
  if [ ! -f "run/$BASE.tsv" ]; then
    ./greet.py -c "$C" -v -o "run/$BASE.tsv" "$@" > "run/$BASE.log" 2>&1
  fi
done
