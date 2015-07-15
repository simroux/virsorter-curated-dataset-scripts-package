# We need the vegan library for the matrix randomization
library(vegan)
# lp_brim package originally developed by Tim Poissot : https://github.com/PoisotLab/lpbrim, slightly modified on the plot part
source("lp_brim.R")
# Reading the matrix of virus-host asociations
mat_2_whole<-as.matrix(read.csv("PhageSorter_all_vs_host.csv",header=T))
# Checking the modularity
result_whole<-bBRIM(mat_2_whole)
# Generating the network plot
plotModules(result_whole,0.4)
# Generating the matrix plot
plotMatrixModules(result_whole,0.4)
# Permuting the matrix 99 times (mat_2_whole_random is a collection of randomly permuted matrix)
mat_2_whole_random<-permatfull(mat_2_whole)
# Checking the modularity of all these randomly permuted matrices and printing the output in a file
capture.output(source("compute_stats.R"),file="Randomization_result.txt")
