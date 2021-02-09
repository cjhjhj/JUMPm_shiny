rm(list = ls())

file = "example_fully_aligned.feature"
data = read.table(file, header = T, sep = "\t", check.names = F)