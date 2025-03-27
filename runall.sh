#!/bin/sh

for C in run/*.ini; do
  BASE=`basename "$C" .ini`
  if [ ! =f "$BASE.tsv" ]; then
    ./greet.py -c "$C" -v -o "$BASE.tsv" "$@" > "$BASE.log" 2>&1
  fi
done
