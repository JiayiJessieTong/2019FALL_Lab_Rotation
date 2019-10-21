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
ngenes = 500
true_counts_res <- SimulateTrueCounts(ncells_total=500,
                                      min_popsize=50, #number of cells in the rarest population
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


PlotTsne(meta=observed_counts[[2]], 
         data=log2(observed_counts[[1]]+1),
         evf_type="7.population", 
         n_pc=20, 
         label='pop', 
         saving = F, 
         plotname="7.population")[[2]]


########################################
# Bacth one: with three cell populations
# number of cell - 1000
# sampling 

########################################
# Bacth two: with four cell populations
# number of cell - 1000

########################################
# Bacth three: with five cell populations
# number of cell - 1000

########################################
# Bacth four: with seven cell populations
# number of cell - 1000






