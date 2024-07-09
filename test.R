require(tidyverse)
require(dynamicTreeCut)
require(sparseBC)
require(s4vd)
require(biclust)
require(SCBiclust)
require(future.apply)
require(BCBC)
require(BCEL)

args <- commandArgs(trailingOnly = TRUE)

R.version

R.home()

print(args)
