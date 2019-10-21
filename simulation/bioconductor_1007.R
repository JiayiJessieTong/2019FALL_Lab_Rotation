if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("splatter")


## ---- include = FALSE------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup-----------------------------------------------------------------
library("splatter")
library("scater")
library("ggplot2")

## ----splat-params----------------------------------------------------------
params <- newSplatParams()
params

## ----nGenes----------------------------------------------------------------
# Set the number of genes to 1000
params <- setParam(params, "nGenes", 1000)

sim <- splatSimulate(params, verbose = FALSE)
dim(sim)

## ----batches---------------------------------------------------------------
# Simulation with two batches of 100 cells
sim <- splatSimulate(params, batchCells = c(100, 100), verbose = FALSE)

# PCA plot using scater
sim <- normalize(sim)
plotPCA(sim, colour_by = "Batch")

## ----batch-factors---------------------------------------------------------
# Simulation with small batch effects
sim1 <- splatSimulate(params, batchCells = c(100, 100),
                      batch.facLoc = 0.001, batch.facScale = 0.001,
                      verbose = FALSE)
sim1 <- normalize(sim1)
plotPCA(sim1, colour_by = "Batch") + ggtitle("Small batch effects")

# Simulation with big batch effects
sim2 <- splatSimulate(params, batchCells = c(100, 100),
                      batch.facLoc = 0.5, batch.facScale = 0.5,
                      verbose = FALSE)
sim2 <- normalize(sim2)
plotPCA(sim2, colour_by = "Batch") + ggtitle("Big batch effects")

## ----outlier-prob----------------------------------------------------------
# Few outliers
sim1 <- splatSimulate(out.prob = 0.001, verbose = FALSE)
ggplot(as.data.frame(rowData(sim1)),
       aes(x = log10(GeneMean), fill = OutlierFactor != 1)) +
  geom_histogram(bins = 100) +
  ggtitle("Few outliers")

# Lots of outliers
sim2 <- splatSimulate(out.prob = 0.2, verbose = FALSE)
ggplot(as.data.frame(rowData(sim2)),
       aes(x = log10(GeneMean), fill = OutlierFactor != 1)) +
  geom_histogram(bins = 100) +
  ggtitle("Lots of outliers")

## ----groups----------------------------------------------------------------
params.groups <- newSplatParams(batchCells = 500, nGenes = 1000)

# One small group, one big group
sim1 <- splatSimulateGroups(params.groups, group.prob = c(0.9, 0.1),
                            verbose = FALSE)
sim1 <- normalize(sim1)
plotPCA(sim1, colour_by = "Group") + ggtitle("One small group, one big group")

# Five groups
sim2 <- splatSimulateGroups(params.groups,
                            group.prob = c(0.2, 0.2, 0.2, 0.2, 0.2),
                            verbose = FALSE)
sim2 <- normalize(sim2)
plotPCA(sim2, colour_by = "Group") + ggtitle("Five groups")

## ----de-prob---------------------------------------------------------------
# Few DE genes
sim1 <- splatSimulateGroups(params.groups, group.prob = c(0.5, 0.5),
                            de.prob = 0.01, verbose = FALSE)
sim1 <- normalize(sim1)
plotPCA(sim1, colour_by = "Group") + ggtitle("Few DE genes")

# Lots of DE genes
sim2 <- splatSimulateGroups(params.groups, group.prob = c(0.5, 0.5),
                            de.prob = 0.3, verbose = FALSE)
sim2 <- normalize(sim2)
plotPCA(sim2, colour_by = "Group") + ggtitle("Lots of DE genes")

## ----de-factors------------------------------------------------------------
# Small DE factors
sim1 <- splatSimulateGroups(params.groups, group.prob = c(0.5, 0.5),
                            de.facLoc = 0.01, verbose = FALSE)
sim1 <- normalize(sim1)
plotPCA(sim1, colour_by = "Group") + ggtitle("Small DE factors")

# Big DE factors
sim2 <- splatSimulateGroups(params.groups, group.prob = c(0.5, 0.5),
                            de.facLoc = 0.3, verbose = FALSE)
sim2 <- normalize(sim2)
plotPCA(sim2, colour_by = "Group") + ggtitle("Big DE factors")

## ----complex-de------------------------------------------------------------
# Different DE probs
sim1 <- splatSimulateGroups(params.groups,
                            group.prob = c(0.2, 0.2, 0.2, 0.2, 0.2),
                            de.prob = c(0.01, 0.01, 0.1, 0.1, 0.3),
                            verbose = FALSE)
sim1 <- normalize(sim1)
plotPCA(sim1, colour_by = "Group") +
  labs(title = "Different DE probabilities",
       caption = paste("Groups 1 and 2 have very few DE genes,",
                       "Groups 3 and 4 have the default number,",
                       "Group 5 has many DE genes"))


# Different DE factors
sim2 <- splatSimulateGroups(params.groups,
                            group.prob = c(0.2, 0.2, 0.2, 0.2, 0.2),
                            de.facLoc = c(0.01, 0.01, 0.1, 0.1, 0.2),
                            de.facScale = c(0.2, 0.5, 0.2, 0.5, 0.4),
                            verbose = FALSE)
sim2 <- normalize(sim2)
plotPCA(sim2, colour_by = "Group") +
  labs(title = "Different DE factors",
       caption = paste("Group 1 has factors with small location (value),",
                       "and scale (variability),",
                       "Group 2 has small location and greater scale.\n",
                       "Groups 3 and 4 have greater location with small,",
                       "and large scales",
                       "Group 5 has bigger factors with moderate variability"))

# Combination of everything
sim3 <- splatSimulateGroups(params.groups,
                            group.prob = c(0.05, 0.2, 0.2, 0.2, 0.35),
                            de.prob = c(0.3, 0.1, 0.2, 0.01, 0.1),
                            de.downProb = c(0.1, 0.4, 0.9, 0.6, 0.5),
                            de.facLoc = c(0.6, 0.1, 0.1, 0.01, 0.2),
                            de.facScale = c(0.1, 0.4, 0.2, 0.5, 0.4),
                            verbose = FALSE)
sim3 <- normalize(sim3)
plotPCA(sim3, colour_by = "Group") +
  labs(title = "Different DE factors",
       caption = paste(
         "Group 1 is small with many very up-regulated DE genes,",
         "Group 2 has the default DE parameters,\n",
         "Group 3 has many down-regulated DE genes,",
         "Group 4 has very few DE genes,",
         "Group 5 is large with moderate DE factors")
  )

## ----paths-----------------------------------------------------------------
# Linear paths
sim1 <- splatSimulatePaths(params.groups,
                           group.prob = c(0.25, 0.25, 0.25, 0.25),
                           de.prob = 0.5, de.facLoc = 0.2,
                           path.from = c(0, 1, 2, 3),
                           verbose = FALSE)
sim1 <- normalize(sim1)
plotPCA(sim1, colour_by = "Group") + ggtitle("Linear paths")

# Branching path
sim2 <- splatSimulatePaths(params.groups,
                           group.prob = c(0.25, 0.25, 0.25, 0.25),
                           de.prob = 0.5, de.facLoc = 0.2,
                           path.from = c(0, 1, 1, 3),
                           verbose = FALSE)
sim2 <- normalize(sim2)
plotPCA(sim2, colour_by = "Group") + ggtitle("Branching path")

## ----paths-steps-----------------------------------------------------------
# Few steps
sim1 <- splatSimulatePaths(params.groups, path.nSteps = 3,
                           de.prob = 0.5, de.facLoc = 0.2, verbose = FALSE)
sim1 <- normalize(sim1)
plotPCA(sim1, colour_by = "Step") + ggtitle("Few steps")

# Lots of steps
sim2 <- splatSimulatePaths(params.groups, path.nSteps = 1000,
                           de.prob = 0.5, de.facLoc = 0.2, verbose = FALSE)
sim2 <- normalize(sim2)
plotPCA(sim2, colour_by = "Step") + ggtitle("Lots of steps")

## ----paths-skew------------------------------------------------------------
# Skew towards the end
sim1 <- splatSimulatePaths(params.groups, path.skew = 0.1,
                           de.prob = 0.5, de.facLoc = 0.2, verbose = FALSE)
sim1 <- normalize(sim1)
plotPCA(sim1, colour_by = "Step") + ggtitle("Skewed towards the end")

