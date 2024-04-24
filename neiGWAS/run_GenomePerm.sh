#!/bin/bash

for i in {1..99}
do
Rscript ./neiGWAS/GenomePerm_i.R $i
done

