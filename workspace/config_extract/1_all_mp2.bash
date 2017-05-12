#!/bin/bash

rm result.log
rm all.xyz
rm *.txt

for x in *.log
do
  base=${x%.log}
  cp ${base}.log result.log
  python extrac_from_log.py
  mv result.txt ${base}.txt
  cat ${base}.txt >> all.xyz
done



