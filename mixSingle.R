library(phangorn)
library(stringr)
library(rphast)
library(rlist)
args = commandArgs(trailingOnly=TRUE)

# USAGE:
# Rscript mixSingle.R IN_DIR REF_FA TREE_NEWICK FOREGROUND_STRING OUT_DIR

# load the input directory containing fasta files
input_dir <- unzip(args[1], list=TRUE)$Name[1]
# read the reference alignment (no acceleration)
#ref_alignment <- as.phyDat(read.dna(args[2], "fasta"))

# read phylogenetic tree
tree <- read.tree(args[3])

# foreground species (support one more species as foreground, separated by comma)
foreground <- args[4]

# define the output directory
out_dir <- args[5]

# build the reference phylogenetic model
if (grepl(".fa",args[2], fixed=TRUE)) {
# read the reference alignment (no acceleration)
ref_alignment <- as.phyDat(read.dna(args[2], "fasta"))
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
#    cat(sprintf("H0 lnl = %f\n", -res_H0$value))
#    cat(sprintf("H1 lnl = %f\n", -res_H1$value))
#    cat(sprintf("H0 scale = %f\n", res_H0$par))
#    cat(sprintf("H1 scale = %f %f\n", res_H1$par[1], res_H1$par[2]))
#    cat(sprintf("likelihood ratio statistic = %f\n", likelihood_ratio))
#    cat(sprintf("gamma shape = %f\n", model$shape))
#    cat(sprintf("p-value = %e\n", p_value))
return(c(res_H0$par, likelihood_ratio))
}

# read the alignment for acceleration test (acceleration might occur)
fasta_names <- list.files(input_dir, pattern="*.fasta", full.names=TRUE)
data.list <- lapply(fasta_names, read.msa)
total_align <- concat.msa(data.list)
total_name <- paste(out_dir, "/total.fa", sep="")
write.msa(total_align, total_name, format="FASTA")
alignment <- as.phyDat(read.dna(total_name, "fasta"))
N <- ncol.msa(total_align)

# read the alignment individually and run LRT on each alignment
align.list <- lapply(fasta_names, read.dna, format="fasta")
alignment.list <- lapply(align.list, as.phyDat)
# perform likleihood ratio test
empirical_result <- data.frame(lik=c())
for (i in 1:length(alignment.list))
{
empirical_result[i,1] <- LRT(model,alignment.list[[i]], foreground)[[2]]
}
#write.csv(empirical_result, file=paste(out_dir,"/emp.csv", sep=""))

# sample 10000 times from the length distribution
length.list <- lapply(data.list, ncol.msa)
data_len <- data.frame(length.list)
write.csv(data_len, file=paste(out_dir,"/length_total.csv", sep=""))
len_sam <- sample(length.list, 10000, replace=TRUE)

# perform likleihood ratio test on the concatenated alignment and get the scaling factor of H0
scale <- LRT(model, alignment, foreground)[[1]]
cat(sprintf("scale = %e\n", scale))

# rescale the reference tree
model$tree$edge.length = model$tree$edge.length * scale

# simulate 10000 sequences and run LRT on each simulated sequence
sim_result <- data.frame(sim_value=c())
sim_data_tree <- list()
for (i in 1:length(len_sam))
{
sim_data_tree[[i]] <- simSeq(model$tree, l=len_sam[[i]], Q=model$Q, bf=model$bf)
#sim_data_tree[[i]]
#print(c(i, len_sam[[i]]))
# perform likleihood ratio test
sim_result[i,1] <- LRT(model, sim_data_tree[[i]], foreground)[[2]]
    }
write.csv(sim_result, file=paste(out_dir,"/sim.csv", sep=""))

# use mixture model to calculate p value
library(ClassComparison)
p2 <- data.frame()
simulation <- sim_result$V1
    for (i in 1:nrow(empirical_result)) {
        p2[i,1] <- length(simulation[simulation >= empirical_result[i,1]])/length(simulation)
        }
colnames(p2) <- "simP"
new <- cbind(empirical_result, p2)
bp <- Bum(new$simP)
empirical_new <- cbind(new, bp@pvals)
write.csv(empirical_new, file=paste(out_dir,"/emp.csv", sep=""))
acc_p <- 1-bp@pihat
cat(sprintf("proportion of acc = %e\n", acc_p))
print(acc_p)
