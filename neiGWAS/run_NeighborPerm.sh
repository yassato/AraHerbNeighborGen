#!/bin/bash

for i in {1..99}
do
Rscript ./neiGWAS/NeighborPerm_i.R $i
done

