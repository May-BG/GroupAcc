library(phangorn)
library(stringr)
library(rphast)
library(ape)
library(apex)
library(rlist)
args = commandArgs(trailingOnly=TRUE)

# USAGE:
# Rscript groupLRT.R IN_DIR_ZIP REF_FA(or REF_MOD) TREE_NEWICK FOREGROUND_STRING OUT_DIR

# load the input concatenated fasta files from zip
fname_zipped <- args[1] #"ZZZ3.zip"

# read the reference alignment (no acceleration)
##ref_alignment <- as.phyDat(read.dna(args[2], "fasta"))

# read phylogenetic tree
tree <- read.tree(args[3])

# foreground species (support one more species as foreground, separated by comma)
foreground <- args[4]

# define the output directory
out_dir <- args[5]

if (grepl(".fa",args[2], fixed=TRUE)) {
# read the reference alignment (no acceleration)
ref_alignment <- as.phyDat(read.dna(args[2], "fasta"))
# build the reference phylogenetic model
ncat <- as.integer(4)
model <- pml(tree, ref_alignment, model="GTR", k=ncat)
model <- optim.pml(model, optBf=T, optQ=T, optGamma=T,
                   control = pml.control(epsilon=1e-08, maxit=1000, trace = 0))

} else {
model <- list.load(args[2])
}
print(model)


# utility function to scale a phylogenetic tree
scale_tree <- function(tree, foreground, foreground_scale, background_scale){
    species <- str_split(foreground, ",")[[1]]

    # find all nodes in foreground
    if (length(species) != 1){
        tip <- which(tree$tip.label %in% species)
        mrca <- mrca.phylo(tree, tip)
        des <- Descendants(tree, mrca, type="all")
        all_node <- c(mrca, des)
    }
    else{
        all_node <- which(tree$tip.label %in% species)
    }

    # find all edges in foreground and background
    foreground_edge <- tree$edge[, 2] %in% all_node
    background_edge <- !foreground_edge

    # rescale tree
    tree$edge.length[foreground_edge] = tree$edge.length[foreground_edge] * foreground_scale
    tree$edge.length[background_edge] = tree$edge.length[background_edge] * background_scale

    return(tree)
}

# model H1 with two scaling factors
H1 <- function(scale, model, data, foreground){
    # rescale tree
    model$tree <- scale_tree(model$tree, foreground, scale[1], scale[2])

    # create a new phylogenetic model
    new_model <- pml(model$tree, data, model$bf, model$Q, model$inv, model$k, model$shape, model=model$model)

    return(-new_model$logLik)
}

# model H0 with one scaling factor. This is a simple wrapper function of model H1.
H0 <- function(scale, model, data, foreground){
    scale <- rep(scale, 2)
    return(H1(scale, model, data, foreground))
}

# likelihood ratio test
LRT <- function(model, data, foreground){
    # optimize model H0
    res_H0 <- optim(1, H0, model=model, data=data, foreground=foreground, method="L-BFGS-B", lower=c(1e-3, 1e-3))

    # optimize model H1, with initial scales from model H0
    res_H1 <- optim(rep(res_H0$par, 2), H1, model=model, data=data, foreground=foreground, method="L-BFGS-B", lower=c(1e-3, 1e-3))

    # likelihood ratio statistic
    likelihood_ratio <- 2 * (res_H0$value - res_H1$value)

    # p-value (chi-square distribution with df = 1)
    p_value <- 1 - pchisq(likelihood_ratio, df=1)

    # BIC of H0
    BIC_0 <- 2 * res_H0$value + 1 * log(N)

    # BIC of H1
    BIC_1 <- 2 * res_H1$value + 2 * log(N)

    # print result
    cat(sprintf("H0 lnl = %f\n", -res_H0$value))
    cat(sprintf("H1 lnl = %f\n", -res_H1$value))
    cat(sprintf("BIC of H0 = %f\n", BIC_0))
    cat(sprintf("BIC of H1 = %f\n", BIC_1))
    cat(sprintf("H0 scale = %f\n", res_H0$par))
    cat(sprintf("H1 scale = %f %f\n", res_H1$par[1], res_H1$par[2]))
    cat(sprintf("likelihood ratio statistic = %f\n", likelihood_ratio))
    cat(sprintf("gamma shape = %f\n", model$shape))
    cat(sprintf("p-value = %e\n", p_value))
}

# read the alignment for acceleration test (acceleration might occur)
#table <- read.csv("partition1.txt", header=FALSE)
#align <- read.dna("catphyml.fas", format="fasta")


#fasta_names <- list.files(input_dir, pattern="*.fasta", full.names=TRUE)
#data.list <- lapply(fasta_names, read.msa)
#total_align <- concat.msa(data.list)
total_name <- paste(out_dir, "/total_new.fa", sep="")
#write.msa(total_align, total_name, format="FASTA")

unzip_dna <- function(x){ 
read.dna(pipe(paste("unzip -p",fname_zipped, x)), format="fasta") -> out
# unlink(paste(x))
 return(out)
 }
names <- unzip(fname_zipped, list=TRUE)$Name
index_fasta <- grep("fasta", names,fixed=TRUE)
fasta_names <- names[index_fasta]

 dna <- lapply(fasta_names, unzip_dna)
x <- new("multidna", dna)
phy <- multidna2multiphyDat(x)
concatenate(x) -> totaldna
ncol(totaldna) -> N
print(N)
concatenate(phy) -> alignment
write.dna(totaldna, paste(total_name), format="fasta")
#alignment <- as.phyDat(read.dna(in_align, "fasta"))

#total_align <- read.msa(in_align, format="FASTA")
#N <- ncol.msa(total_align)

# perform likleihood ratio test
LRT(model, alignment, foreground)
