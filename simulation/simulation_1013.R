# this code is written by Jessie
# Last update: 10/20/2019
# this code is the simulation code for batch effect 
# aim: to figure out if the order of the batch when applying MNN method matters

# install package
library("devtools")
devtools::install_github("YosefLab/SymSim")
library("SymSim")
library("ape")

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
table(cell_meta[,2])

data(gene_len_pool)
gene_len <- sample(gene_len_pool[which(gene_len_pool < 200)], ngenes, replace = FALSE)

# Start the clock! -----------------------------
ptm <- proc.time()

observed_counts <- True2ObservedCounts(true_counts=true_counts_res[[1]], 
                                       meta_cell=true_counts_res[[3]], 
                                       protocol="nonUMI", 
                                       alpha_mean=0.1, 
                                       alpha_sd=0.05, 
                                       gene_len=gene_len, 
                                       depth_mean=1e5, 
                                       depth_sd=3e3,
                                       nbatch = 4)

# Stop the clock ----------------------------
proc.time() - ptm

# plot the combined population from four different batches
PlotTsne(meta=observed_counts[[2]], 
         data=log2(observed_counts[[1]]+1),
         evf_type="7.population", 
         n_pc=20, 
         label='pop', 
         saving = F, 
         plotname="7.population")[[2]]

# funtion to divide the batch and population based on the combined data
batch_divider <- function(batch_id, batch_name, size_batch_num){
  
  # extract the batch with batch_id from the whole population
  batch_id = which(observed_counts[[2]]$batch == batch_id)
  
  # based on the id, extract the cell meta
  batch_observed_cell_meta = observed_counts[[2]][batch_id,]
  
  # print the cell meta information
  batch_cell_meta_dim = dim(batch_observed_cell_meta) 
  
  # based on the id, extract the observed counts
  batch_observed_counts = observed_counts[[1]][,batch_id]
  
  # print the count information
  batch_count_dim = dim(batch_observed_counts) 
  
  # extract the batch with size_batch_num: 
  # e.g., batch 1 has three pop, batch 2 has 4 pop, batch 3 has 5 pop, batch 4 has 7 pop
  batch_sub_id = batch_observed_cell_meta$pop %in% c(1:size_batch_num)
  
  # based on the batch sub id, extract the cell meta
  batch_observed_cell_meta_sub = batch_observed_cell_meta[batch_sub_id,]
  
  # print the cell meta sub information
  batch_cell_meta_dim_sub = dim(batch_observed_cell_meta_sub) 
  
  # based on the sub id, extract the observed counts
  batch_observed_counts_sub = batch_observed_counts[,batch_sub_id]
  
  # print the counbt information
  batch_count_dim_sub = dim(batch_observed_counts_sub) 
  
  # plot the cells population, subgroup 3 populations
  tsne_true_counts = PlotTsne(meta=batch_observed_cell_meta_sub, 
                              data=log2(batch_observed_counts_sub+1),
                              evf_type=paste("Batch", batch_name, size_batch_num,"population"), 
                              n_pc=20, 
                              label='pop', 
                              saving = F, 
                              plotname=paste("Batch", batch_name, size_batch_num,"population"))
  
  return(list(tsne = tsne_true_counts,
              batch_observed_cell_meta_sub = batch_observed_cell_meta_sub,
              batch_cell_meta_dim_sub = batch_cell_meta_dim_sub,
              batch_observed_counts_sub = batch_observed_counts_sub,
              batch_count_dim_sub = batch_count_dim_sub))
}

########################################
# Bacth one: with three cell populations
batch1 = batch_divider(batch_id = 1,
                       batch_name = "One",
                       size_batch_num = 3)
# plot
batch1$tsne[[2]]

# calculate the proportion of each population
round(table(batch1$batch_observed_cell_meta_sub[,2])/
  sum(table(batch1$batch_observed_cell_meta_sub[,2])),2)

########################################
# Bacth two: with four cell populations
batch2 = batch_divider(batch_id = 2,
                       batch_name = "Two",
                       size_batch_num = 4)
batch2$tsne[[2]]

# calculate the proportion of each population
round(table(batch2$batch_observed_cell_meta_sub[,2])/
        sum(table(batch2$batch_observed_cell_meta_sub[,2])),2)

########################################
# Bacth three: with five cell populations
batch3 = batch_divider(batch_id = 3,
                       batch_name = "Three",
                       size_batch_num = 5)
batch3$tsne[[2]]

# calculate the proportion of each population
round(table(batch3$batch_observed_cell_meta_sub[,2])/
        sum(table(batch3$batch_observed_cell_meta_sub[,2])),2)

########################################
# Bacth four: with seven cell populations
batch4 = batch_divider(batch_id = 4,
                       batch_name = "Four",
                       size_batch_num = 7)
batch4$tsne[[2]]

# calculate the proportion of each population
round(table(batch4$batch_observed_cell_meta_sub[,2])/
        sum(table(batch4$batch_observed_cell_meta_sub[,2])),2)





