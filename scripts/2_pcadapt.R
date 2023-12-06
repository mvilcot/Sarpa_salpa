# ---- read data ----
source("scripts/0_wrapper.R")

# filters <- "callrateind0.50_callrateloci0.70_maf0.05"
filters <- "callrateind0.50_callrateloci0.70_maf0.05_RMindoutliers"

genlight <- 
  readRDS(paste0("intermediate/Genlight_Sarpa_salpa_", filters, ".RDS"))




# ---- PCADAPT ----
## ---- read Genlight into PCAdapt format ----

# transform to data frame
gdf <- as.data.frame(genlight)
genotype <- t(gdf)
gdf$pop <- genlight$pop # set population information

# transform to pcadapt format
genotype_pca <- read.pcadapt(genotype)



## ---- PCAdapt ----
# first run
PCADAPT <- pcadapt(genotype_pca, K = 15)

# check number of PCs (=K) to retain
plotSCREE <- plot(PCADAPT, option = "screeplot") 
ggsave(paste0("results/pcadapt_Screeplot_", filters, ".png"), 
       plotSCREE, height = 5, width = 8)

# rerun with good number of PCs
K <- 2
PCADAPT <- pcadapt(genotype_pca, K = K)


## ---- Graphic visualisation ----

# plot PC1-PC2
plotPCs <- plot(PCADAPT, option = "scores", i = 1, j = 2, pop = gdf$pop)
ggsave(paste0("results/pcadapt_PC1PC2_", filters, "_K", K, ".png"), 
       plotPCs, height = 5, width = 8)

# Manhattan plot
plotMH <- plot(PCADAPT, option = "manhattan") 
ggsave(paste0("results/pcadapt_Manhattan_", filters, "_K", K, ".png"), 
       plotMH, height = 5, width = 8)

# Check the expected uniform distribution of the p-values
plotQQ <- plot(PCADAPT, option = "qqplot", threshold = 0.1) 
ggsave(paste0("results/pcadapt_QQplot_", filters, "_K", K, ".png"), 
       plotQQ, height = 5, width = 8)

# The excess of small p-values indicates the presence of outliers
plotDIST <- plot(PCADAPT, option = "stat.distribution") 
ggsave(paste0("results/pcadapt_Distribution_", filters, "_K", K, ".png"), 
       plotDIST, height = 5, width = 8)

# presence of outliers is also visible when plotting a histogram of the test statistic Dj.
plotHIST <- ggplot() + xlab("p-values") +
  geom_histogram(aes(PCADAPT$pvalues), bins = 50, fill = "orange", colour = "black")
ggsave(paste0("results/pcadapt_Histogram_", filters, "_K", K, ".png"), 
       plotHIST, height = 5, width = 8)



## ---- Choose a cutoff for outlier detection ----
# q-values
qval <- qvalue(PCADAPT$pvalues)$qvalues
alpha <- 0.1
outliersQ <- which(qval < alpha)
length(outliersQ) ## N = 856 outliers

# # Bonferroni correction (conservative)
# padj <- p.adjust(PCADAPT$pvalues,method="bonferroni")
# alpha <- 0.1
# outliersB <- which(padj < alpha)
# length(outliersB) ## N = 153 outliers
# 
# # compare both methods
# setdiff(outliersB, outliersQ) # in outliersB but not in outliersQ
# setdiff(outliersQ, outliersB) # in outliersQ but not in outliersB

# keep outliers based on 10% Q-values
outliers <- outliersQ

# extract corresponding SNP names and positions
outlier_loci <- data.frame(loci = rownames(genotype)[outliers], 
                           position = outliers)
write.table(outlier_loci, 
            file = paste0('intermediate/pcadapt_outliersQ_loci_position_', filters, '.txt'), 
            sep = "\t", quote = F, row.names = F, col.names = F)

cat(length(outliers), "outlier loci")



## ---- Remove outlier loci from genlight ----
genlightNEUTRAL <- gl.drop.loc(genlight, loc.list = outlier_loci$loci)
genlightNEUTRAL 
# 23469 SNPs remaining

genlightADAPT <- gl.keep.loc(genlight, loc.list = outlier_loci$loci)
genlightADAPT 
# 856  SNPs remaining

genlightNEUTRAL %>% 
  saveRDS(paste0("intermediate/Genlight_Sarpa_salpa_", filters, "_RMpcadaptQ0.1.RDS"))

genlightADAPT %>% 
  saveRDS(paste0("intermediate/Genlight_Sarpa_salpa_", filters, "_pcadaptQ0.1.RDS"))





# ---- dbRDA - environment ----


