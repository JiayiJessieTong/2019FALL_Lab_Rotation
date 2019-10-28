# this code is written by Jessie
# Last update: 10/26/2019
# this code is the simulation code for batch effect 
# aim: to figure out if the order of the batch when applying MNN method matters
# Based on the discussion with Prof. Li on 10/21/2019
# The batched effects can be generated right after the true counts are generated
# which means, after use the first funtion to generate the true with 7 population
# then extract 3, then add batches
# extract 4th, then add batches
# etc.

# install package
library("devtools")
# devtools::install_github("YosefLab/SymSim")
library("SymSim")
library("ape")
library(dplyr)
library(Seurat)
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("batchelor")
# BiocManager::install("scater")
library(batchelor)
library(scater)
library(Rtsne)

########### Simulation Idea
# 1000 cells for each batch
# 20000 genes for each cell
### ---- Batch One: 3 population [0.20, 0.30, 0.50]
### ---- Batch Two: 4 population [0.05, 0.15, 0.35, 0.45]
### ---- Batch Three: 5 population [0.10, 0.15, 0.20, 0.25, 0.30]
### ---- Batch Four: 7 population [0.05, 0.08, 0.12, 0.15,  0.18, 0.20, 0.22 ]

###########
ngenes = 200
true_counts_res <- SimulateTrueCounts(ncells_total=1000,
                                      min_popsize=100, #number of cells in the rarest population
                                      randseed = 0, # random seed to reproduce the results
                                      ngenes=ngenes, # number of genes
                                      phyla=rtree(7),
                                      evf_type="discrete",
                                      scale_s = 0.9,
                                      nevf=30,
                                      Sigma=0.1,
                                      n_de_evf=10) 

tsne_true_counts <- PlotTsne(meta=true_counts_res[[3]], 
                             data=log2(true_counts_res[[1]]+1), 
                             evf_type="discrete", 
                             n_pc=20, label='pop', 
                             saving = F, 
                             plotname="Total - 7 pops")

# plot of true counts with 7 different cells population
tsne_true_counts[[2]]

# extract cell-meta from true_counts_res
cell_meta = true_counts_res$cell_meta

# Pooled data: check the number of cells in each population
print(table(cell_meta[,2]))

data(gene_len_pool)
# function
extract_batch <- function(id_list, num_batch){
  
  # extract the population with the id I want
  id = c(true_counts_res[[3]]$pop %in% id_list)
  
  # extract the meta_cell info with the id
  meta_cell = true_counts_res[[3]][id,]
  
  # extract the true_counts info with the id
  true_counts = true_counts_res[[1]][,id]
  
  # get the observed counts with function True2ObservedCounts
  # num_batch is the number of batches that the id_list (cell population) has
  gene_len <- sample(gene_len_pool, ngenes, replace = FALSE)
  
  observed_counts <- True2ObservedCounts(true_counts=true_counts, 
                                         meta_cell=meta_cell, 
                                         protocol="UMI", 
                                         alpha_mean=0.1, 
                                         alpha_sd=0.05, 
                                         gene_len=gene_len, 
                                         depth_mean=1e5, 
                                         depth_sd=3e3,
                                         nbatch = num_batch)
  
  return(observed_counts)
}


pop123 = extract_batch(c(1,2,3), 4)
pop4 = extract_batch(c(4), 3)
pop5 = extract_batch(c(5), 2)
pop67 = extract_batch(c(6,7), 1)


########################################
# Bacth one: with three cell populations
batch1_id = which(pop123$cell_meta$batch == 1)
batch1_meta = pop123[[2]][batch1_id,]
batch1_count = pop123[[1]][,batch1_id]
PlotTsne(meta=batch1_meta, 
         data=log2(batch1_count+1),
         evf_type="Batch One (3 pop)", 
         n_pc=20, 
         label='pop', 
         saving = F, 
         plotname="Batch One (3 pop)")[[2]]
# proportion of each population: number of cells in each pop / total number of cells
print(paste("The proportions of three cells pop are",
            round(table(batch1_meta$pop) / dim(batch1_meta)[1],2)))

########################################
# Bacth two: with four cell populations
batch2_id = which(pop123$cell_meta$batch == 2)
batch2_id_2 = which(pop4$cell_meta$batch == 1)
batch2_meta = rbind(pop123[[2]][batch2_id,], pop4[[2]][batch2_id_2,])
batch2_count = cbind(pop123[[1]][,batch2_id], pop4[[1]][,batch2_id_2])
PlotTsne(meta=batch2_meta, 
         data=log2(batch2_count+1),
         evf_type="Batch Two (4 pop)", 
         n_pc=20, 
         label='pop', 
         saving = F, 
         plotname="Batch Two (4 pop)")[[2]]
# proportion of each population: number of cells in each pop / total number of cells
print(paste("The proportions of four cells pop are",
            round(table(batch2_meta$pop) / dim(batch2_meta)[1],2)))

########################################
# Bacth three: with five cell populations
batch3_id = which(pop123$cell_meta$batch == 3)
batch3_id_2 = which(pop4$cell_meta$batch == 2)
batch3_id_3 = which(pop5$cell_meta$batch == 1)
batch3_meta = rbind(pop123[[2]][batch3_id,], 
                    pop4[[2]][batch3_id_2,],
                    pop5[[2]][batch3_id_3,])
batch3_count = cbind(pop123[[1]][,batch3_id], 
                     pop4[[1]][,batch3_id_2],
                     pop5[[1]][,batch3_id_3])
PlotTsne(meta=batch3_meta, 
         data=log2(batch3_count+1),
         evf_type="Batch Three (5 pop)", 
         n_pc=20, 
         label='pop', 
         saving = F, 
         plotname="Batch Three (5 pop)")[[2]]
# proportion of each population: number of cells in each pop / total number of cells
print(paste("The proportions of five cells pop are",
            round(table(batch3_meta$pop) / dim(batch3_meta)[1],2)))

########################################
# Bacth four: with seven cell populations
batch4_id = which(pop123$cell_meta$batch == 4)
batch4_id_2 = which(pop4$cell_meta$batch == 3)
batch4_id_3 = which(pop5$cell_meta$batch == 2)
batch4_id_4 = which(pop67$cell_meta$batch == 1)
batch4_meta = rbind(pop123[[2]][batch4_id,], 
                    pop4[[2]][batch4_id_2,],
                    pop5[[2]][batch4_id_3,],
                    pop67[[2]][batch4_id_4,])
batch4_count = cbind(pop123[[1]][,batch4_id], 
                     pop4[[1]][,batch4_id_2],
                     pop5[[1]][,batch4_id_3],
                     pop67[[1]][,batch4_id_4])
PlotTsne(meta=batch4_meta, 
         data=log2(batch4_count+1),
         evf_type="Batch Four (7 pop)", 
         n_pc=20, 
         label='pop', 
         saving = F, 
         plotname="Batch Four (7 pop)")[[2]]
# proportion of each population: number of cells in each pop / total number of cells
print(paste("The proportions of seven cells pop are",
            round(table(batch4_meta$pop) / dim(batch4_meta)[1],2)))

#########################################################
########## WITHOUT CORRECTION
#########################################################
logcounts_1 <- log2(batch1_count + 1)
logcounts_2 <- log2(batch2_count + 1)
logcounts_3 <- log2(batch3_count + 1)
logcounts_4 <- log2(batch4_count + 1)

########## combine batch 1 and batch 2 without correction
PlotTsne(meta=rbind(batch1_meta, batch2_meta), 
         data=log2(cbind(batch1_count, batch2_count)+1),
         evf_type="Combined Batch One and Two w/o Correction", 
         n_pc=20, 
         label='pop', 
         saving = F, 
         plotname="Combined Batch One and Two w/o Correction")[[2]]
batch12_sce <- SingleCellExperiment(
  assays = list(counts = cbind(batch1_count, batch2_count),
                logcounts = log2(cbind(batch1_count, batch2_count)+1))
) 
batch12_sce <- runTSNE(batch12_sce)
batch12_sce$batch <- factor(rep(c(1,2), 
                         c(ncol(logcounts_1), 
                           ncol(logcounts_2))))
plotTSNE(batch12_sce, run_args=list(perplexity = 10), colour_by = "batch")

########## combine batch 1, batch 2, batch3 without correction
PlotTsne(meta=rbind(batch1_meta, batch2_meta, batch3_meta), 
         data=log2(cbind(batch1_count, batch2_count, batch3_count)+1),
         evf_type="Combined Batch One, Two, and Three w/o Correction", 
         n_pc=20, 
         label='pop', 
         saving = F, 
         plotname="Combined Batch One, Two, and Three w/o Correction")[[2]]
batch123_sce <- SingleCellExperiment(
  assays = list(counts = cbind(batch1_count, batch2_count, batch3_count),
                logcounts = log2(cbind(batch1_count, batch2_count, batch3_count)+1))
) 
batch123_sce <- runTSNE(batch123_sce)
batch123_sce$batch <- factor(rep(c(1,2, 3), 
                                c(ncol(logcounts_1), 
                                  ncol(logcounts_2),
                                  ncol(logcounts_3))))
plotTSNE(batch123_sce, run_args=list(perplexity = 10), colour_by = "batch")

########## combine batch 1, batch 2, batch3, batch 4 without correction
PlotTsne(meta=rbind(batch1_meta, batch2_meta, batch3_meta, batch4_meta), 
         data=log2(cbind(batch1_count, batch2_count, batch3_count, batch4_count)+1),
         evf_type="Combined Batch One, Two, Three, Four w/o Correction", 
         n_pc=20, 
         label='pop', 
         saving = F, 
         plotname="Combined Batch One, Two, Three, Four w/o Correction")[[2]]
batch1234_sce <- SingleCellExperiment(
  assays = list(counts = cbind(batch1_count, batch2_count, batch3_count,
                               batch4_count),
                logcounts = log2(cbind(batch1_count, batch2_count, 
                                       batch3_count, batch4_count)+1))
) 
batch1234_sce <- runTSNE(batch1234_sce)
batch1234_sce$batch <- factor(rep(c(1,2,3,4), 
                                 c(ncol(logcounts_1), 
                                   ncol(logcounts_2),
                                   ncol(logcounts_3),
                                   ncol(logcounts_4))))
plotTSNE(batch1234_sce, run_args=list(perplexity = 10), colour_by = "batch")

#########################################################
########## WITH CORRECTION： MNN
#########################################################

###########
out1 <- mnnCorrect(logcounts_1, logcounts_2, logcounts_3, logcounts_4)
counts(out1) <- assay(out1)
logcounts(out1) <- assay(out1)
out1 <- runTSNE(out1)
out1$batch <- factor(rep(1:4, c(ncol(logcounts_1), 
                                ncol(logcounts_2),
                                ncol(logcounts_3),
                                ncol(logcounts_4))))
plotTSNE(out1, run_args=list(perplexity = 10), colour_by = "batch")

# out1@colData@listData[["batch"]]

##########
out2 <- mnnCorrect(logcounts_4, logcounts_1, logcounts_2, logcounts_3)
counts(out2) <- assay(out2)
logcounts(out2) <- assay(out2)
out2 <- runTSNE(out2)
out2$batch <- factor(rep(c(4,1,2,3), 
                         c(ncol(logcounts_4), 
                           ncol(logcounts_1),
                           ncol(logcounts_2),
                           ncol(logcounts_3))))
plotTSNE(out2, run_args=list(perplexity = 10), colour_by = "batch")

# out2@colData@listData[["batch"]]



# out2 <- fastMNN(logcounts_1, logcounts_2, logcounts_3, logcounts_4)
# dim(reducedDim(out2)) 
# assay(out2, "reconstructed")
# set.seed(607)
# out2 <- runTSNE(out2, use_dimred="corrected")
# out2$batch <- factor(out2$batch)
# plotTSNE(out2, colour_by="batch")
# 
# out3 <- fastMNN(logcounts_4, logcounts_1, logcounts_2, logcounts_3)
# dim(reducedDim(out3)) 
# assay(out3, "reconstructed")
# set.seed(607)
# out3 <- runTSNE(out3, use_dimred="corrected")
# out3$batch <- factor(out3$batch)
# plotTSNE(out3, colour_by="batch")
# 
# 
# logNormCounts(batch1_count)
# 
# 
# 
# example_sce <- mockSCE()
# example_sce <- logNormCounts(example_sce)
# example_sce <- runPCA(example_sce)
# 
# ## Examples plotting PC1 and PC2
# plotPCA(fastMNN(logcounts_1, logcounts_2, logcounts_3, logcounts_4))
# plotPCA(example_sce)
# plotPCA(example_sce, colour_by = "Cell_Cycle")
# plotPCA(example_sce, colour_by = "Cell_Cycle", shape_by = "Treatment")
# plotPCA(example_sce, colour_by = "Cell_Cycle", shape_by = "Treatment",
#         size_by = "Mutation_Status")
# 
# ## Force legend to appear for shape:
# example_subset <- example_sce[, example_sce$Treatment == "treat1"]
# plotPCA(example_subset, colour_by = "Cell_Cycle", shape_by = "Treatment", 
#         by_show_single = TRUE)
# 
# ## Examples plotting more than 2 PCs
# plotPCA(example_sce, ncomponents = 4, colour_by = "Treatment",
#         shape_by = "Mutation_Status")
# 
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("scater")
# library(scater)
# BiocManager::valid()
# install.packages("BiocManager")
# BiocManager::install("scran")
# 
# 
# 
# data("sc_example_counts")
# data("sc_example_cell_info")
# 
# pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
# example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd)
# drop_genes <- apply(exprs(example_sceset), 1, function(x) {var(x) == 0})
# example_sceset <- example_sceset[!drop_genes, ]
# 
# ## Examples plotting PC1 and PC2
# plotTSNE(example_sceset, perplexity = 10)
# plotTSNE(example_sceset, colour_by = "Cell_Cycle", perplexity = 10)
# plotTSNE(example_sceset, colour_by = "Cell_Cycle", shape_by = "Treatment",
#          perplexity = 10)
# plotTSNE(example_sceset, colour_by = "Cell_Cycle", shape_by = "Treatment",
#          size_by = "Mutation_Status", perplexity = 10)
# plotTSNE(example_sceset, shape_by = "Treatment", size_by = "Mutation_Status",
#          perplexity = 5)
# plotTSNE(example_sceset, feature_set = 1:100, colour_by = "Treatment",
#          shape_by = "Mutation_Status", perplexity = 5)
# 
# plotTSNE(example_sceset, shape_by = "Treatment", return_SCESet = TRUE,
#          perplexity = 10)
# source("https://bioconductor.org/biocLite.R")
# biocLite("scater")
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("scater")
# library(scater)
