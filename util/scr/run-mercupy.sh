#!/bin/bash
for i in $(seq 1 1)
do
    make prepare run out ref continue
done
date +%s.%N > end.sig
