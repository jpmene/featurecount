#!/usr/bin/env Rscript
library(readr)

args <- commandArgs(TRUE)
print(args)

data <- read_delim(args[1], 
                   "\t", escape_double = FALSE, trim_ws = TRUE, 
                   skip = 1)

data2 <- read_delim(args[2], 
                   "\t", escape_double = FALSE, trim_ws = TRUE, 
                   skip = 1)

colnames(data) <- c("Geneid","Chr","Start","End","Strand","Length","Count")
colnames(data2) <- c("Geneid","Chr","Start","End","Strand","Length","Count")

data <- data[,c(1,7)]

data2 <- data2[,c(1,7)]

data_merge <- merge(data,data2,by = "Geneid")

name_file <- paste(args[1],"plot",sep = "_")

jpeg(filename = name_file)

plot(data_merge$Count.x,
     data_merge$Count.y,
     log = "xy",
     xlab=args[1],
     ylab=args[2])

dev.off()