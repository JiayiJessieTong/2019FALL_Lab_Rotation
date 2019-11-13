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
library(pracma)


#### seurat workflow function
seurat_workflow <- function(count_matrix){
  data_seurat <- CreateSeuratObject(counts = count_matrix,
                                    min.cells = 3, min.features = 200)
  data_seurat <- NormalizeData(data_seurat)
  data_seurat <- ScaleData(data_seurat)
  data_seurat <- FindVariableFeatures(data_seurat, selection.method = "vst", nfeatures = 2000)
  data_seurat <- RunPCA(data_seurat, features = VariableFeatures(object = data_seurat))
  data_seurat <- FindNeighbors(data_seurat, dims = 1:10)
  data_seurat <- FindClusters(data_seurat, resolution = 0.5)
  data_seurat <- RunUMAP(data_seurat, dims = 1:10)
  DimPlot(data_seurat, reduction = "umap")
}

# function to extract batch
extract_batch <- function(id_list, num_batch, 
                          alpha_mean, alpha_sd,
                          depth_mean, depth_sd){
  
  data(gene_len_pool)
  
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
                                         alpha_mean=alpha_mean,
                                         alpha_sd=alpha_sd,
                                         gene_len=gene_len,
                                         depth_mean=depth_mean,
                                         depth_sd=depth_sd,
                                         nbatch = num_batch)
  
  
  return(observed_counts)
}

# create data
ngenes = 500
ncells_total = 2000
true_counts_res <- SimulateTrueCounts(ncells_total=ncells_total,
                                      min_popsize=200, #number of cells in the rarest population
                                      randseed = 0, # random seed to reproduce the results
                                      ngenes=ngenes, # number of genes
                                      phyla=rtree(7),
                                      evf_type="discrete",
                                      scale_s = 0.9,
                                      nevf=100,
                                      Sigma=0.5,
                                      n_de_evf=50)

# plot cells
tsne_true_counts <- PlotTsne(meta=true_counts_res[[3]],
                             data=log2(true_counts_res[[1]]+1),
                             evf_type="discrete",
                             n_pc=20, label='pop',
                             saving = F,
                             plotname="Total - 7 pops")
# plot of true counts with 7 different cells population
tsne_true_counts[[2]]


# plot clustering cells
colnames(true_counts_res[[1]]) = true_counts_res[[3]][,1]
row.names(true_counts_res[[1]]) = c(1:ngenes)
seurat_workflow(true_counts_res[[1]])


## print the proprotion of zeroes in the true counts data
sum(true_counts_res[[1]] == 0)/(ngenes * ncells_total)

# extract cell-meta from true_counts_res
cell_meta = true_counts_res$cell_meta

# Pooled data: check the number of cells in each population
print(table(cell_meta[,2]))


#### create cells population
pop123 = extract_batch(c(1,2,3), 4, 
                       alpha_mean = 0.05, alpha_sd = 0.02,
                        depth_mean = 1e5, depth_sd = 3e3)
pop4 = extract_batch(c(4), 3, 
                     alpha_mean = 0.01, alpha_sd = 0.05,
                     depth_mean = 1e4, depth_sd = 3e3)
pop5 = extract_batch(c(5), 2, 
                     alpha_mean = 0.1, alpha_sd = 0.03,
                     depth_mean = 2e3, depth_sd = 3e3)
pop67 = extract_batch(c(6,7), 1, 
                      alpha_mean = 0.06, alpha_sd = 0.01,
                      depth_mean = 5e4, depth_sd = 3e3)


########################################
# Bacth one: with three cell populations
batch1_id = which(pop123$cell_meta$batch == 1)
batch1_meta = pop123[[2]][batch1_id,]
batch1_count = pop123[[1]][,batch1_id]

## print the proprotion of zeroes in the observed data
sum(batch1_count == 0)/(dim(batch1_count)[1] * dim(batch1_count)[2])

## plot
PlotTsne(meta=batch1_meta,
         data=log2(batch1_count+1),
         evf_type="Batch One (3 pop)",
         n_pc=20,
         label='pop',
         saving = F,
         plotname="Batch One (3 pop)")[[2]]

# # proportion of each population: number of cells in each pop / total number of cells
print(paste("The proportions of three cells pop are",
            round(table(batch1_meta$pop) / dim(batch1_meta)[1],2)))

### signle cell experinment to plot cells
batch1_sce <- SingleCellExperiment(
  assays = list(counts = batch1_count,
                logcounts = log2(batch1_count+1))
) 

batch1_tsne <- runTSNE(batch1_sce)
batch1_tsne$batch <- factor(rep(c(1), dim(batch1_sce)[2]))
batch1_tsne$cell_type <- factor(batch1_meta[,2])
plotTSNE(batch1_tsne, run_args=list(perplexity = 10), colour_by = "batch")
plotTSNE(batch1_tsne, run_args=list(perplexity = 10), colour_by = "cell_type")


########################################
# Bacth two: with four cell populations
batch2_id = which(pop123$cell_meta$batch == 2)
batch2_id_2 = which(pop4$cell_meta$batch == 1)
batch2_meta = rbind(pop123[[2]][batch2_id,], pop4[[2]][batch2_id_2,])
batch2_count = cbind(pop123[[1]][,batch2_id], pop4[[1]][,batch2_id_2])

## print the proprotion of zeroes in the observed data
sum(batch2_count == 0)/(dim(batch2_count)[1] * dim(batch2_count)[2])

## plot
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

### signle cell experinment to plot cells
batch2_sce <- SingleCellExperiment(
  assays = list(counts = batch2_count,
                logcounts = log2(batch2_count+1))
) 

batch2_tsne <- runTSNE(batch2_sce)
batch2_tsne$batch <- factor(rep(c(2), dim(batch2_sce)[2]))
batch2_tsne$cell_type <- factor(batch2_meta[,2])
plotTSNE(batch2_tsne, run_args=list(perplexity = 10), colour_by = "batch")
plotTSNE(batch2_tsne, run_args=list(perplexity = 10), colour_by = "cell_type")



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

## print the proprotion of zeroes in the observed data
sum(batch3_count == 0)/(dim(batch3_count)[1] * dim(batch3_count)[2])

## plot
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

### signle cell experinment to plot cells
batch3_sce <- SingleCellExperiment(
  assays = list(counts = batch3_count,
                logcounts = log2(batch3_count+1))
) 

batch3_tsne <- runTSNE(batch3_sce)
batch3_tsne$batch <- factor(rep(c(3), dim(batch3_sce)[2]))
batch3_tsne$cell_type <- factor(batch3_meta[,2])
plotTSNE(batch3_tsne, run_args=list(perplexity = 10), colour_by = "batch")
plotTSNE(batch3_tsne, run_args=list(perplexity = 10), colour_by = "cell_type")


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
## print the proprotion of zeroes in the observed data
sum(batch4_count == 0)/(dim(batch4_count)[1] * dim(batch4_count)[2])

## plot
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

### signle cell experinment to plot cells
batch4_sce <- SingleCellExperiment(
  assays = list(counts = batch4_count,
                logcounts = log2(batch4_count+1))
) 

batch4_tsne <- runTSNE(batch4_sce)
batch4_tsne$batch <- factor(rep(c(4), dim(batch4_sce)[2]))
batch4_tsne$cell_type <- factor(batch4_meta[,2])
plotTSNE(batch4_tsne, run_args=list(perplexity = 10), colour_by = "batch")
plotTSNE(batch4_tsne, run_args=list(perplexity = 10), colour_by = "cell_type")

########################################################
################ MNN ####################################
#########################################################
#########################################################
logcounts_1 <- log2(batch1_count + 1)
logcounts_2 <- log2(batch2_count + 1)
logcounts_3 <- log2(batch3_count + 1)
logcounts_4 <- log2(batch4_count + 1)

########## combine batch 1 and batch 2 without correction
#########################################################
########## WITHOUT CORRECTION
#########################################################
# PlotTsne(meta=rbind(batch1_meta, batch2_meta), 
#          data=log2(cbind(batch1_count, batch2_count)+1),
#          evf_type="Combined Batch One and Two w/o Correction", 
#          n_pc=20, 
#          label='pop', 
#          saving = F, 
#          plotname="Combined Batch One and Two w/o Correction")[[2]]
batch12_sce <- SingleCellExperiment(
  assays = list(counts = cbind(batch1_count, batch2_count),
                logcounts = log2(cbind(batch1_count, batch2_count)+1))
) 

# tsne
batch12_tsne <- runTSNE(batch12_sce)
batch12_tsne$batch <- factor(rep(c(1,2), 
                                c(ncol(logcounts_1), 
                                  ncol(logcounts_2))))
batch12_tsne$cell_type <- factor(c(batch1_meta[,2], batch2_meta[,2]))
plotTSNE(batch12_tsne, run_args=list(perplexity = 10), colour_by = "batch")
plotTSNE(batch12_tsne, run_args=list(perplexity = 10), colour_by = "cell_type")

# pca
batch12_pca <- runPCA(batch12_sce)
batch12_pca$batch <- factor(rep(c(1,2), 
                                 c(ncol(logcounts_1), 
                                   ncol(logcounts_2))))
batch12_pca$cell_type <- factor(c(batch1_meta[,2], batch2_meta[,2]))
plotTSNE(batch12_pca, run_args=list(perplexity = 10), colour_by = "batch")
plotTSNE(batch12_pca, run_args=list(perplexity = 10), colour_by = "cell_type")


# #########################################################
# ########## WITH CORRECTION： MNN
# #########################################################
out1 <- mnnCorrect(logcounts_1, logcounts_2)

names(assays(out1))

counts(out1) <- assay(out1)
logcounts(out1) <- assay(out1)

assays(out1)["corrected"]

# tsne
out1_tsne <- runTSNE(out1)
out1_tsne$batch <- factor(rep(1:2, c(ncol(logcounts_1), 
                                ncol(logcounts_2))))
out1_tsne$cell_type <- factor(c(batch1_meta[,2], batch2_meta[,2]))
plotTSNE(out1_tsne, run_args=list(perplexity = 10), colour_by = "batch")
plotTSNE(out1_tsne, run_args=list(perplexity = 10), colour_by = "cell_type")

# pca
out1_pca <- runPCA(out1)
out1_pca$batch <- factor(rep(1:2, c(ncol(logcounts_1), 
                                ncol(logcounts_2))))
out1_pca$cell_type <- factor(c(batch1_meta[,2], batch2_meta[,2]))
plotPCA(out1_pca, run_args=list(perplexity = 10), colour_by = "batch")
plotPCA(out1_pca, run_args=list(perplexity = 10), colour_by = "cell_type")



########## combine batch 1, batch 2, batch3 without correction
#########################################################
########## WITHOUT CORRECTION
#########################################################
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
# tsne
batch123_tsne <- runTSNE(batch123_sce)
batch123_tsne$batch <- factor(rep(c(1,2,3), 
                                 c(ncol(logcounts_1), 
                                   ncol(logcounts_2),
                                   ncol(logcounts_3))))
batch123_tsne$cell_type <- factor(c(batch1_meta[,2], 
                                    batch2_meta[,2],
                                    batch3_meta[,2]))
plotTSNE(batch123_tsne, run_args=list(perplexity = 10), colour_by = "batch")
plotTSNE(batch123_tsne, run_args=list(perplexity = 10), colour_by = "cell_type")

# pca
batch123_pca <- runPCA(batch123_sce)
batch123_pca$batch <- factor(rep(c(1,2,3), 
                                  c(ncol(logcounts_1), 
                                    ncol(logcounts_2),
                                    ncol(logcounts_3))))
batch123_pca$cell_type <- factor(c(batch1_meta[,2], 
                                    batch2_meta[,2],
                                    batch3_meta[,2]))
plotPCA(batch123_pca, run_args=list(perplexity = 10), colour_by = "batch")
plotPCA(batch123_pca, run_args=list(perplexity = 10), colour_by = "cell_type")



# #########################################################
# ########## WITH CORRECTION： MNN
# #########################################################
out2 <- mnnCorrect(logcounts_1, logcounts_2, logcounts_3)
out2 <- mnnCorrect(logcounts_2, logcounts_3, logcounts_1)
counts(out2) <- assay(out2)
logcounts(out2) <- assay(out2)

# tsne
out2_tsne <- runTSNE(out2)
out2_tsne$batch <- factor(rep(1:3, c(ncol(logcounts_1), 
                                     ncol(logcounts_2),
                                     ncol(logcounts_3))))
out2_tsne$cell_type <- factor(c(batch1_meta[,2], 
                                batch2_meta[,2],
                                batch3_meta[,2]))
plotTSNE(out2_tsne, run_args=list(perplexity = 10), colour_by = "batch")
plotTSNE(out2_tsne, run_args=list(perplexity = 10), colour_by = "cell_type")

# pca
out2_pca <- runPCA(out2)
out2_pca$batch <- factor(rep(1:3, c(ncol(logcounts_1), 
                                    ncol(logcounts_2),
                                    ncol(logcounts_3))))
out2_pca$cell_type <- factor(c(batch1_meta[,2], 
                               batch2_meta[,2],
                               batch3_meta[,2]))
plotPCA(out2_pca, run_args=list(perplexity = 10), colour_by = "batch")
plotPCA(out2_pca, run_args=list(perplexity = 10), colour_by = "cell_type")





########## combine batch 1, batch 2, batch3, batch 4 without correction
#########################################################
########## WITHOUT CORRECTION
#########################################################
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


# tsne
batch1234_tsne <- runTSNE(batch1234_sce)
batch1234_tsne$batch <- factor(rep(c(1,2,3,4), 
                                   c(ncol(logcounts_1), 
                                     ncol(logcounts_2),
                                     ncol(logcounts_3),
                                     ncol(logcounts_4))))
batch1234_tsne$cell_type <- factor(c(batch1_meta[,2], 
                                    batch2_meta[,2],
                                    batch3_meta[,2],
                                    batch4_meta[,2]))
plotTSNE(batch1234_tsne, run_args=list(perplexity = 10), colour_by = "batch")
plotTSNE(batch1234_tsne, run_args=list(perplexity = 10), colour_by = "cell_type")

# pca
batch1234_tsne <- runPCA(batch1234_sce)
batch1234_tsne$batch <- factor(rep(c(1,2,3,4), 
                                 c(ncol(logcounts_1), 
                                   ncol(logcounts_2),
                                   ncol(logcounts_3),
                                   ncol(logcounts_4))))
batch1234_tsne$cell_type <- factor(c(batch1_meta[,2], 
                                   batch2_meta[,2],
                                   batch3_meta[,2],
                                   batch4_meta[,2]))
plotPCA(batch1234_tsne, run_args=list(perplexity = 10), colour_by = "batch")
plotPCA(batch1234_tsne, run_args=list(perplexity = 10), colour_by = "cell_type")




# #########################################################
# ########## WITH CORRECTION： MNN
# #########################################################
out3 <- mnnCorrect(logcounts_1, logcounts_2, logcounts_3, logcounts_4)
counts(out3) <- assay(out3)
logcounts(out3) <- assay(out3)

# tsne
out3_tsne <- runTSNE(out3)
out3_tsne$batch <- factor(rep(1:4, c(ncol(logcounts_1), 
                                     ncol(logcounts_2),
                                     ncol(logcounts_3),
                                     ncol(logcounts_4))))
out3_tsne$cell_type <- factor(c(batch1_meta[,2], 
                                batch2_meta[,2],
                                batch3_meta[,2],
                                batch4_meta[,2]))
plotTSNE(out3_tsne, run_args=list(perplexity = 10), colour_by = "batch")
plotTSNE(out3_tsne, run_args=list(perplexity = 10), colour_by = "cell_type")

# pca
out3_pca <- runPCA(out3)
out3_pca$batch <- factor(rep(1:4,  c(ncol(logcounts_1), 
                                     ncol(logcounts_2),
                                     ncol(logcounts_3),
                                     ncol(logcounts_4))))
out3_pca$cell_type <- ffactor(c(batch1_meta[,2], 
                                batch2_meta[,2],
                                batch3_meta[,2],
                                batch4_meta[,2]))
plotPCA(out3_pca, run_args=list(perplexity = 10), colour_by = "batch")
plotPCA(out3_pca, run_args=list(perplexity = 10), colour_by = "cell_type")





# #########################################################
# ########## WITH CORRECTION： MNN
# #########################################################
# 
# ###########
# out1 <- mnnCorrect(logcounts_1, logcounts_2, logcounts_3, logcounts_4)
# counts(out1) <- assay(out1)
# logcounts(out1) <- assay(out1)
# out1 <- runTSNE(out1)
# out1$batch <- factor(rep(1:4, c(ncol(logcounts_1), 
#                                 ncol(logcounts_2),
#                                 ncol(logcounts_3),
#                                 ncol(logcounts_4))))
# plotTSNE(out1, run_args=list(perplexity = 10), colour_by = "cell")
# 
# 
# # out1@colData@listData[["batch"]]
# 
# ##########
# out2 <- mnnCorrect(logcounts_4, logcounts_1, logcounts_2, logcounts_3)
# counts(out2) <- assay(out2)
# logcounts(out2) <- assay(out2)
# out2 <- runTSNE(out2)
# out2$batch <- factor(rep(c(4,1,2,3), 
#                          c(ncol(logcounts_4), 
#                            ncol(logcounts_1),
#                            ncol(logcounts_2),
#                            ncol(logcounts_3))))
# plotTSNE(out2, run_args=list(perplexity = 10), colour_by = "batch")
# # out2@colData@listData[["batch"]]
