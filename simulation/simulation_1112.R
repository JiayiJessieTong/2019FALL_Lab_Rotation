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
library(gridExtra)
library(randomcoloR)

#### seurat workflow function
#### seurat workflow function
seurat_workflow <- function(count_matrix, cellType_info, batch_info){
  data_seurat <- CreateSeuratObject(counts = count_matrix,
                                    min.cells = 3, min.features = 10)
  data_seurat <- NormalizeData(data_seurat)
  data_seurat <- ScaleData(data_seurat)
  data_seurat <- FindVariableFeatures(data_seurat, selection.method = "vst", nfeatures = 2000)
  data_seurat <- RunPCA(data_seurat, features = VariableFeatures(object = data_seurat))
  data_seurat <- FindNeighbors(data_seurat, dims = 1:10)
  data_seurat <- FindClusters(data_seurat, resolution = 0.5)
  data_seurat <- RunUMAP(data_seurat, dims = 1:10)
  
  n1 = length(levels(data_seurat@meta.data[["RNA_snn_res.0.5"]]))
  c1 = distinctColorPalette(n1)
  # clustering
  p1 = DimPlot(data_seurat, reduction = "umap", cols = c1)
  
  # based on real cell type
  data_seurat$CellType = Idents(data_seurat)
  # Idents(data_seurat) <- c(batch1_meta[,2], batch2_meta[,2]) - 1
  Idents(data_seurat) <- cellType_info
  n2 = length(unique(Idents(data_seurat)))
  c2 = distinctColorPalette(n2)
  p2 = DimPlot(data_seurat, reduction = "umap", cols = c2)
  
  # based on batch
  data_seurat$Batch = Idents(data_seurat)
  # Idents(data_seurat) = c(batch1_meta$batch, batch2_meta$batch)
  Idents(data_seurat) = batch_info
  n3 = length(unique(Idents(data_seurat)))
  c3 = distinctColorPalette(n3)
  p3 = DimPlot(data_seurat, reduction = "umap", cols = c3)
  
  # plot = grid.arrange(p1, p2, p3, nrow = 1)
  plot <- arrangeGrob(p1, p2, p3, nrow = 1)
  
  return(list(plot = plot))
}
# function to extract batch
extract_batch <- function(true_counts_res, ngenes,
                          id_list, num_batch, 
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
## plotting
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
# colnames(true_counts_res[[1]]) = true_counts_res[[3]][,1]
# row.names(true_counts_res[[1]]) = c(1:ngenes)
# tmp_0 = seurat_workflow(true_counts_res[[1]])


## print the proprotion of zeroes in the true counts data
sum(true_counts_res[[1]] == 0)/(ngenes * ncells_total)

# extract cell-meta from true_counts_res
cell_meta = true_counts_res$cell_meta

# Pooled data: check the number of cells in each population
print(table(cell_meta[,2]))


#### create cells population
pop123 = extract_batch(true_counts_res, ngenes,
                       c(1,2,3), 4, 
                       alpha_mean = 0.05, alpha_sd = 0.02,
                       depth_mean = 1e5, depth_sd = 3e3)
pop4 = extract_batch(true_counts_res, ngenes,
                     c(4), 3, 
                     alpha_mean = 0.01, alpha_sd = 0.05,
                     depth_mean = 1e4, depth_sd = 3e3)
pop5 = extract_batch(true_counts_res, ngenes,
                     c(5), 2, 
                     alpha_mean = 0.1, alpha_sd = 0.03,
                     depth_mean = 2e3, depth_sd = 3e3)
pop67 = extract_batch(true_counts_res, ngenes,
                      c(6,7), 1, 
                      alpha_mean = 0.06, alpha_sd = 0.01,
                      depth_mean = 5e4, depth_sd = 3e3)


########################################
# Bacth one: with three cell populations
batch1_id = which(pop123$cell_meta$batch == 1)
batch1_meta = pop123[[2]][batch1_id,]
batch1_count = pop123[[1]][,batch1_id]

## print the proprotion of zeroes in the observed data
sum(batch1_count == 0)/(dim(batch1_count)[1] * dim(batch1_count)[2])

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

# ########## combine batch 1 and batch 2 without correction
# #########################################################
# ########## WITHOUT CORRECTION
# #########################################################
# batch12_sce <- SingleCellExperiment(
#   assays = list(counts = cbind(batch1_count, batch2_count),
#                 logcounts = log2(cbind(batch1_count, batch2_count)+1))
# ) 
# 
# # tsne
# batch12_tsne <- runTSNE(batch12_sce)
# batch12_tsne$batch <- factor(rep(c(1,2), 
#                                  c(ncol(logcounts_1), 
#                                    ncol(logcounts_2))))
# batch12_tsne$cell_type <- factor(c(batch1_meta[,2], batch2_meta[,2]))
# plotTSNE(batch12_tsne, run_args=list(perplexity = 10), colour_by = "batch")
# plotTSNE(batch12_tsne, run_args=list(perplexity = 10), colour_by = "cell_type")
# 
# # clustering analysis
# batch12_tsne_count_matrix = batch12_tsne@assays@data@listData[["counts"]]
# colnames(batch12_tsne_count_matrix) = as.character(c(1:dim(batch12_tsne_count_matrix)[2]))
# row.names(batch12_tsne_count_matrix) = as.character(c(1:ngenes))
# # seurat_workflow(count_matrix = batch12_tsne_count_matrix)
# r = seurat_workflow(count_matrix = batch12_tsne_count_matrix,
#                 cellType_info =
#                   c(batch1_meta[,2], batch2_meta[,2]) - 1,
#                 batch_info = c(batch1_meta$batch, batch2_meta$batch))
# grid.newpage()
# grid.draw(r$plot)
# 
# 
# # #########################################################
# # ########## WITH CORRECTION： MNN
# # #########################################################
# # option 1
# out1 <- mnnCorrect(logcounts_1, logcounts_2)
# counts(out1) <- out1@assays@data@listData[["corrected"]]
# logcounts(out1) <- log2(out1@assays@data@listData[["corrected"]] + 1)
# 
# # tsne
# out1_tsne <- runTSNE(out1)
# out1_tsne$batch <- factor(rep(1:2, c(ncol(logcounts_1), 
#                                      ncol(logcounts_2))))
# out1_tsne$cell_type <- factor(c(batch1_meta[,2], batch2_meta[,2]))
# plotTSNE(out1_tsne, run_args=list(perplexity = 10), colour_by = "batch")
# plotTSNE(out1_tsne, run_args=list(perplexity = 10), colour_by = "cell_type")
# 
# # clustering analysis
# out1_count_matrix = out1@assays@data@listData[["corrected"]]
# 
# # proportion of negative numbers in the corrected count matrix
# length(which(out1_count_matrix<0))/length(out1_count_matrix)
# 
# colnames(out1_count_matrix) = as.character(c(1:dim(out1_count_matrix)[2]))
# row.names(out1_count_matrix) = as.character(c(1:ngenes))
# # seurat_workflow(count_matrix = out1_count_matrix)
# r = seurat_workflow(count_matrix = out1_count_matrix,
#                 cellType_info =
#                   c(batch1_meta[,2], batch2_meta[,2]) - 1,
#                 batch_info = c(batch1_meta$batch, batch2_meta$batch))
# grid.newpage()
# grid.draw(r$plot)
# 
# # option 2
# out2 <- mnnCorrect(logcounts_2, logcounts_1)
# counts(out2) <- out2@assays@data@listData[["corrected"]]
# logcounts(out2) <- log2(out2@assays@data@listData[["corrected"]] + 1)
# 
# #### check similar
# par(mfrow = c(1,2))
# # batch1
# check_similar_1 = out1@assays@data@listData[["corrected"]][1:ncol(logcounts_1),]
# check_similar_2 = out2@assays@data@listData[["corrected"]][1:ncol(logcounts_1),]
# row_mean_1 = apply(check_similar_1, 1, mean)
# row_mean_2 = apply(check_similar_2, 1, mean)
# plot(row_mean_1, row_mean_2, 
#      xlab = "Order option 1",
#      ylab = "Order option 2",
#      main = "Batch 1")
# lines(x = c(0,100), y = c(0,100), col = "red", lwd = 2)
# 
# # batch2
# check_similar_1 = out1@assays@data@listData[["corrected"]][1:ncol(logcounts_2),]
# check_similar_2 = out2@assays@data@listData[["corrected"]][1:ncol(logcounts_2),]
# row_mean_1 = apply(check_similar_1, 1, mean)
# row_mean_2 = apply(check_similar_2, 1, mean)
# plot(row_mean_1, row_mean_2, 
#      xlab = "Order option 1",
#      ylab = "Order option 2",
#      main = "Batch 2")
# lines(x = c(0,100), y = c(0,100), col = "red", lwd = 2)
# par(mfrow = c(1,1))
# 
# 
# ########## combine batch 1, batch 2, batch3 without correction
# #########################################################
# ########## WITHOUT CORRECTION
# #########################################################
# batch123_sce <- SingleCellExperiment(
#   assays = list(counts = cbind(batch1_count, batch2_count, batch3_count),
#                 logcounts = log2(cbind(batch1_count, batch2_count, batch3_count)+1))
# ) 
# # tsne
# batch123_tsne <- runTSNE(batch123_sce)
# batch123_tsne$batch <- factor(rep(c(1,2,3), 
#                                   c(ncol(logcounts_1), 
#                                     ncol(logcounts_2),
#                                     ncol(logcounts_3))))
# batch123_tsne$cell_type <- factor(c(batch1_meta[,2], 
#                                     batch2_meta[,2],
#                                     batch3_meta[,2]))
# plotTSNE(batch123_tsne, run_args=list(perplexity = 10), colour_by = "batch")
# plotTSNE(batch123_tsne, run_args=list(perplexity = 10), colour_by = "cell_type")
# 
# # clustering analysis
# batch123_tsne_count_matrix = batch123_tsne@assays@data@listData[["counts"]]
# colnames(batch123_tsne_count_matrix) = as.character(c(1:dim(batch123_tsne_count_matrix)[2]))
# row.names(batch123_tsne_count_matrix) = as.character(c(1:ngenes))
# # seurat_workflow(count_matrix = batch123_tsne_count_matrix)
# 
# r = seurat_workflow(count_matrix = batch123_tsne_count_matrix,
#                 cellType_info =
#                   c(batch1_meta[,2], batch2_meta[,2],
#                     batch3_meta[,2]) - 1,
#                 batch_info = c(batch1_meta$batch, batch2_meta$batch,
#                                batch3_meta$batch))
# grid.newpage()
# grid.draw(r$plot)
# 
# 
# 
# 
# # #########################################################
# # ########## WITH CORRECTION： MNN
# # #########################################################
# # function
# three_func <- function(l1, l2, l3,
#                        batch_list,
#                        batch_num_list,
#                        batch_meta_list,
#                        batch_info){
#   out <- mnnCorrect(l1, l2, l3)
#   counts(out) <- out@assays@data@listData[["corrected"]]
#   logcounts(out) <- log2(out@assays@data@listData[["corrected"]] + 1)
#   
#   # tsne
#   out_tsne <- runTSNE(out)
#   out_tsne$batch <- factor(rep(batch_list, batch_num_list))
#   out_tsne$cell_type <- factor(batch_meta_list)
#   plot1 = plotTSNE(out_tsne, run_args=list(perplexity = 10), colour_by = "batch")
#   plot2 = plotTSNE(out_tsne, run_args=list(perplexity = 10), colour_by = "cell_type")
#   
#   
#   # clustering analysis
#   out_count_matrix = out@assays@data@listData[["corrected"]]
#   colnames(out_count_matrix) = as.character(c(1:dim(out_count_matrix)[2]))
#   row.names(out_count_matrix) = as.character(c(1:ngenes))
#   plot3 = seurat_workflow(count_matrix = out_count_matrix,
#                           cellType_info = batch_meta_list,
#                           batch_info = batch_info)
#   
#   return(list(tsne_batch = plot1, 
#               tsne_cell_type = plot2, 
#               clustering = plot3,
#               correct_count = out_count_matrix))
# }
# ############## RUN three different options
# ncol_list = c(ncol(logcounts_1), ncol(logcounts_2), ncol(logcounts_3))
# # option 1
# list_1 = c(1,2,3)
# batch_num_list_1 = ncol_list[list_1]
# meta_list_1 = c(batch1_meta[,2], 
#                 batch2_meta[,2],
#                 batch3_meta[,2])
# batch_info_1 = c(batch1_meta$batch,
#                batch2_meta$batch,
#                batch3_meta$batch)
# option1 <- three_func(logcounts_1,
#                       logcounts_2,
#                       logcounts_3,
#                       list_1,
#                       batch_num_list_1,
#                       meta_list_1,
#                       batch_info = batch_info_1)
# option1$tsne_batch
# option1$tsne_cell_type
# r = option1$clustering
# grid.newpage()
# grid.draw(r$plot)
# # m1 <- option1$correct_count
# 
# # option 2
# list_2 = c(2,3,1)
# batch_num_list_2 = ncol_list[list_2]
# meta_list_2 =  c(batch2_meta[,2], 
#                  batch3_meta[,2],
#                  batch1_meta[,2])
# batch_info_2 = c(batch2_meta$batch,
#                batch3_meta$batch,
#                batch1_meta$batch)
# option2 <- three_func(logcounts_2,
#                       logcounts_3,
#                       logcounts_1,
#                       list_2,
#                       batch_num_list_2,
#                       meta_list_2,
#                       batch_info = batch_info_2)
# option2$tsne_batch
# option2$tsne_cell_type
# r = option2$clustering
# grid.newpage()
# grid.draw(r$plot)
# # m2 <- option2$correct_count
# 
# # option 3
# list_3 = c(3,2,1)
# batch_num_list_3 = ncol_list[list_3]
# meta_list_3 =  c(batch3_meta[,2], 
#                  batch2_meta[,2],
#                  batch1_meta[,2])
# batch_info_3 = c(batch3_meta$batch,
#                  batch2_meta$batch,
#                  batch1_meta$batch)
# option3 <- three_func(logcounts_3,
#                       logcounts_2,
#                       logcounts_1,
#                       list_3,
#                       batch_num_list_3,
#                       meta_list_3,
#                       batch_info = batch_info_3)
# option3$tsne_batch
# option3$tsne_cell_type
# r = option3$clustering
# grid.newpage()
# grid.draw(r$plot)
# 
# #### check similar
# par(mfrow = c(1,3))
# # batch1
# check_similar_1 = option1$correct_count[1:ncol(logcounts_1),]
# check_similar_2 = option2$correct_count[1:ncol(logcounts_1),]
# check_similar_3 = option3$correct_count[1:ncol(logcounts_1),]
# row_mean_1 = apply(check_similar_1, 1, mean)
# row_mean_2 = apply(check_similar_2, 1, mean)
# row_mean_3 = apply(check_similar_1, 1, mean)
# plot(row_mean_1, row_mean_2, 
#      xlab = "Order option 1",
#      ylab = "Order option 2",
#      main = "Batch 1")
# lines(x = c(0,100), y = c(0,100), col = "red", lwd = 2)
# plot(row_mean_1, row_mean_3, 
#      xlab = "Order option 1",
#      ylab = "Order option 3",
#      main = "Batch 1")
# lines(x = c(0,100), y = c(0,100), col = "red", lwd = 2)
# plot(row_mean_2, row_mean_3, 
#      xlab = "Order option 2",
#      ylab = "Order option 3",
#      main = "Batch 1")
# lines(x = c(0,100), y = c(0,100), col = "red", lwd = 2)
# 
# # batch2
# check_similar_1 = option1$correct_count[1:ncol(logcounts_2),]
# check_similar_2 = option2$correct_count[1:ncol(logcounts_2),]
# check_similar_3 = option3$correct_count[1:ncol(logcounts_2),]
# row_mean_1 = apply(check_similar_1, 1, mean)
# row_mean_2 = apply(check_similar_2, 1, mean)
# row_mean_3 = apply(check_similar_1, 1, mean)
# plot(row_mean_1, row_mean_2, 
#      xlab = "Order option 1",
#      ylab = "Order option 2",
#      main = "Batch 1")
# lines(x = c(0,100), y = c(0,100), col = "red", lwd = 2)
# plot(row_mean_1, row_mean_3, 
#      xlab = "Order option 1",
#      ylab = "Order option 3",
#      main = "Batch 1")
# lines(x = c(0,100), y = c(0,100), col = "red", lwd = 2)
# plot(row_mean_2, row_mean_3, 
#      xlab = "Order option 2",
#      ylab = "Order option 3",
#      main = "Batch 1")
# lines(x = c(0,100), y = c(0,100), col = "red", lwd = 2)
# 
# # batch2
# check_similar_1 = option1$correct_count[1:ncol(logcounts_3),]
# check_similar_2 = option2$correct_count[1:ncol(logcounts_3),]
# check_similar_3 = option3$correct_count[1:ncol(logcounts_3),]
# row_mean_1 = apply(check_similar_1, 1, mean)
# row_mean_2 = apply(check_similar_2, 1, mean)
# row_mean_3 = apply(check_similar_1, 1, mean)
# plot(row_mean_1, row_mean_2, 
#      xlab = "Order option 1",
#      ylab = "Order option 2",
#      main = "Batch 1")
# lines(x = c(0,100), y = c(0,100), col = "red", lwd = 2)
# plot(row_mean_1, row_mean_3, 
#      xlab = "Order option 1",
#      ylab = "Order option 3",
#      main = "Batch 1")
# lines(x = c(0,100), y = c(0,100), col = "red", lwd = 2)
# plot(row_mean_2, row_mean_3, 
#      xlab = "Order option 2",
#      ylab = "Order option 3",
#      main = "Batch 1")
# lines(x = c(0,100), y = c(0,100), col = "red", lwd = 2)
# 


########## combine batch 1, batch 2, batch3, batch 4 without correction
#########################################################
########## WITHOUT CORRECTION
#########################################################
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

# clustering analysis
batch1234_tsne_count_matrix = batch1234_tsne@assays@data@listData[["counts"]]
colnames(batch1234_tsne_count_matrix) = as.character(c(1:dim(batch1234_tsne_count_matrix)[2]))
row.names(batch1234_tsne_count_matrix) = as.character(c(1:ngenes))
# seurat_workflow(count_matrix = batch1234_tsne_count_matrix)

r = seurat_workflow(count_matrix = batch1234_tsne_count_matrix,
                    cellType_info =
                      c(batch1_meta[,2], batch2_meta[,2],
                        batch3_meta[,2], batch4_meta[,2]) - 1,
                    batch_info = c(batch1_meta$batch, batch2_meta$batch,
                                   batch3_meta$batch, batch4_meta$batch))
grid.newpage()
grid.draw(r$plot)





# #########################################################
# ########## WITH CORRECTION： MNN
# #########################################################
# function
four_func <- function(l1, l2, l3, l4,
                      batch_list,
                      batch_num_list,
                      batch_meta_list,
                      batch_info){
  out <- mnnCorrect(l1, l2, l3, l4)
  counts(out) <- out@assays@data@listData[["corrected"]]
  logcounts(out) <- log2(out@assays@data@listData[["corrected"]] + 1)
  
  # tsne
  out_tsne <- runTSNE(out)
  out_tsne$batch <- factor(rep(batch_list, batch_num_list))
  out_tsne$cell_type <- factor(batch_meta_list)
  plot1 = plotTSNE(out_tsne, run_args=list(perplexity = 10), colour_by = "batch")
  plot2 = plotTSNE(out_tsne, run_args=list(perplexity = 10), colour_by = "cell_type")
  
  
  # clustering analysis
  out_count_matrix = out@assays@data@listData[["corrected"]]
  colnames(out_count_matrix) = as.character(c(1:dim(out_count_matrix)[2]))
  row.names(out_count_matrix) = as.character(c(1:ngenes))
  plot3 = seurat_workflow(count_matrix = out_count_matrix,
                          cellType_info = batch_meta_list,
                          batch_info = batch_info)
  
  return(list(tsne_batch = plot1, 
              tsne_cell_type = plot2, 
              clustering = plot3,
              correct_count = out_count_matrix, 
              out = out))
}

############## RUN three different options
ncol_list = c(ncol(logcounts_1), ncol(logcounts_2), 
              ncol(logcounts_3), ncol(logcounts_4))
# option 1
list_1 = c(1,2,3,4)
batch_num_list_1 = ncol_list[list_1]
meta_list_1 = c(batch1_meta[,2], 
                batch2_meta[,2],
                batch3_meta[,2],
                batch4_meta[,2])
batch_info_1 = c(batch1_meta$batch,
                 batch2_meta$batch,
                 batch3_meta$batch,
                 batch4_meta$batch)
option1 <- four_func(logcounts_1,
                     logcounts_2,
                     logcounts_3,
                     logcounts_4,
                     list_1,
                     batch_num_list_1,
                     meta_list_1,
                     batch_info = batch_info_1)
option1$tsne_batch
option1$tsne_cell_type
r = option1$clustering
grid.newpage()
grid.draw(r$plot)
out1 = option1$out
# out1 = out1@assays@data@listData[["corrected"]]



# option 2
list_2 = c(4,1,2,3)
batch_num_list_2 = ncol_list[list_2]
meta_list_2 = c(batch4_meta[,2], 
                batch1_meta[,2],
                batch2_meta[,2],
                batch3_meta[,2])
batch_info_2 = c(batch4_meta$batch,
                 batch1_meta$batch,
                 batch2_meta$batch,
                 batch3_meta$batch)
option2 <- four_func(logcounts_4,
                     logcounts_1,
                     logcounts_2,
                     logcounts_3,
                     list_2,
                     batch_num_list_2,
                     meta_list_2,
                     batch_info = batch_info_2)
option2$tsne_batch
option2$tsne_cell_type
r = option2$clustering
grid.newpage()
grid.draw(r$plot)
out2 = option2$out
# m2 <- option2$correct_count


# option 3
list_3 = c(3,4,2,1)
batch_num_list_3 = ncol_list[list_3]
meta_list_3 = c(batch3_meta[,2], 
                batch4_meta[,2],
                batch2_meta[,2],
                batch1_meta[,2])
batch_info_3 = c(batch3_meta$batch,
                 batch4_meta$batch,
                 batch2_meta$batch,
                 batch1_meta$batch)
option3 <- four_func(logcounts_3,
                     logcounts_4,
                     logcounts_2,
                     logcounts_1,
                     list_3,
                     batch_num_list_3,
                     meta_list_3,
                     batch_info = batch_info_3)
option3$tsne_batch
option3$tsne_cell_type
r = option3$clustering
grid.newpage()
grid.draw(r$plot)
out3 = option3$out
# m3 <- option3$correct_count

check_similar <- function(out1, out2){
  #### check similar
  par(mfrow = c(2,2))
  # batch1
  check_similar_1 = out1@assays@data@listData[["corrected"]][,1:ncol(logcounts_1)]
  check_similar_2 = out2@assays@data@listData[["corrected"]][,1:ncol(logcounts_1)]
  row_mean_1 = apply(check_similar_1, 1, mean)
  row_mean_2 = apply(check_similar_2, 1, mean)
  plot(row_mean_1, row_mean_2,
       xlab = "Order option 1",
       ylab = "Order option 2",
       main = "Batch 1")
  lines(x = c(0,100), y = c(0,100), col = "red", lwd = 2)
  
  # batch2
  check_similar_1 = out1@assays@data@listData[["corrected"]][,ncol(logcounts_1)+1:ncol(logcounts_2)]
  check_similar_2 = out2@assays@data@listData[["corrected"]][,ncol(logcounts_1)+1:ncol(logcounts_2)]
  row_mean_1 = apply(check_similar_1, 1, mean)
  row_mean_2 = apply(check_similar_2, 1, mean)
  plot(row_mean_1, row_mean_2,
       xlab = "Order option 1",
       ylab = "Order option 2",
       main = "Batch 2")
  lines(x = c(0,100), y = c(0,100), col = "red", lwd = 2)
  
  # batch3
  check_similar_1 = out1@assays@data@listData[["corrected"]][,ncol(logcounts_1)+1+ncol(logcounts_2):ncol(logcounts_3)]
  check_similar_2 = out2@assays@data@listData[["corrected"]][,ncol(logcounts_1)+1+ncol(logcounts_2):ncol(logcounts_3)]
  row_mean_1 = apply(check_similar_1, 1, mean)
  row_mean_2 = apply(check_similar_2, 1, mean)
  plot(row_mean_1, row_mean_2,
       xlab = "Order option 1",
       ylab = "Order option 2",
       main = "Batch 3")
  lines(x = c(0,100), y = c(0,100), col = "red", lwd = 2)
  
  # batch4
  check_similar_1 = out1@assays@data@listData[["corrected"]][,ncol(logcounts_1)+1+ncol(logcounts_2)+ncol(logcounts_3):ncol(logcounts_4)]
  check_similar_2 = out2@assays@data@listData[["corrected"]][,ncol(logcounts_1)+1+ncol(logcounts_2)+ncol(logcounts_3):ncol(logcounts_4)]
  row_mean_1 = apply(check_similar_1, 1, mean)
  row_mean_2 = apply(check_similar_2, 1, mean)
  plot(row_mean_1, row_mean_2,
       xlab = "Order option 1",
       ylab = "Order option 2",
       main = "Batch 4")
  lines(x = c(0,100), y = c(0,100), col = "red", lwd = 2)
}


# Frobeius norm a matrix
# norm(m1-m2, type = "F")
# norm(m1-m3, type = "F")
# norm(m2-m3, type = "F")
