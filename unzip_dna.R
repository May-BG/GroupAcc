args <- commandArgs(T)

require("rphast")
require("ape")

fname_zipped <- args[1]

 unzip_dna <- function(x){ 
 read.dna(unzip(fname_zipped, x), format="fasta") -> out
 unlink(paste(x))
 return(out)
 }
 dna <- lapply(names, unzip_dna)



