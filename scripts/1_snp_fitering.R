source("scripts/0_wrapper.R")


# ---- Read data ----
DART <- 
  read.csv("data/Report_DSarp22-7382_SNP_2.csv", header = F)

metadata <- readRDS('data/metadata_DSarp_complete.RDS')


# ---- Reorder individuals ----
# reorder Dart columns (individuals) according to metadata order
DARTt <- as.data.frame(t(DART))
DARTt <- DARTt[25:nrow(DARTt),c(3,4,5,7)]
colnames(DARTt) <- c('plate_id', 'row', 'column', 'id')
DARTt <- as_tibble(DARTt)
DARTt$plate_id <- as.double(DARTt$plate_id)
DARTt$column <- as.double(DARTt$column)

data_samples <- 
  metadata %>% 
  right_join(DARTt, by = join_by(plate_id, row, column)) %>% 
  relocate(id)
data_samples %>% write_csv('intermediate/metadata_samples_sequenced.csv')

samples_not_sequenced <-
  metadata %>% 
  full_join(DARTt, by = join_by(plate_id, row, column)) %>% 
  dplyr::filter(is.na(id)) %>% 
  relocate(id)
samples_not_sequenced %>%  write_csv('intermediate/metadata_samples_NOT_sequenced.csv')

cat("nb individuals sequenced =", nrow(data_samples), '/', nrow(metadata))


# get individuals order from metadata
ind_order <- 
  match(data_samples$id, DART[7, ]) # DART 7th row stores individuals names


# reorder Dart
DARTord <- DART[, c(1:24, ind_order)]


# save reordered dart SNP file
DARTord %>% write.table("intermediate/Report_DSarp22-7382_SNP_2_ordered.csv",
              sep = ",", row.names = F, col.names = F)


# ---- Format genlight ----

# read sorted DART file as genlight
gl <- 
  gl.read.dart(filename = "intermediate/Report_DSarp22-7382_SNP_2_ordered.csv",
               ind.metafile = "intermediate/metadata_samples_sequenced.csv") # takes the order of the DART columns for individuals

# add location as pop info
gl@pop <- 
  gl@other$ind.metrics$location

# save as RDS 
gl %>% 
  saveRDS(file = "intermediate/Genlight_Sarpa_salpa_ordered.RDS")

# number of SNPs and individuals: 105144 SNPs, 369 individuals
gl
gl@ind.names
gl@pop



# ---- Graphic visualisation ----
gl <- 
  readRDS(file = "intermediate/Genlight_Sarpa_salpa_ordered.RDS")

## bases ----
png("results/Genlight_raw_report_bases.png", height = 8, width = 10, units = 'in', res = 500)
gl.report.bases(gl)
dev.off()


## callrate by individual ----
png("results/Genlight_raw_report_callrate_ind.png", height = 8, width = 10, units = 'in', res = 500)
gl.report.callrate(gl, method = 'ind') # missing data by individuals
dev.off()

# save callrate by individual in dataframe
options(max.print=2000)
t1 <- capture.output(gl.report.callrate(gl, method = "ind"))
t2 <- t1[38:length(t1)-1]
t3 <- as.data.frame(t2)
t4 <- 
  read.table(text=sub("^(\\S+)\\s+.*\\s+(\\S+)$", "\\1 \\2", t3$t2),
             header=FALSE, stringsAsFactors= FALSE) %>% 
  as_tibble()
colnames(t4) <- t4[1,]
t4 <- t4[-1,]
t5 <- 
  t4 %>%
  dplyr::rename(id = "ind_name") %>% 
  dplyr::rename(site = "pop") %>% 
  dplyr::rename(callrate = "missing_data")

options(max.print=1000)

t5 %>% 
  write_csv("results/Genlight_raw_report_callrate_ind.csv")


## callrate by loci ----
png("results/Genlight_raw_report_callrate_loc.png", height = 8, width = 10, units = 'in', res = 500)
gl.report.callrate(gl, method = 'loc')
dev.off()

## maf ----
png("results/Genlight_raw_report_MAF.png", height = 8, width = 10, units = 'in', res = 500)
gl.report.maf(gl) # mininum allele frequency
dev.off()

## reproducibility ----
png("results/Genlight_raw_report_reproducibility.png", height = 8, width = 10, units = 'in', res = 500)
gl.report.reproducibility(gl)
dev.off()

## read depth ----
png("results/Genlight_raw_report_read_depth.png", height = 8, width = 10, units = 'in', res = 500)
gl.report.rdepth(gl)
dev.off()

## secondaries ----
png("results/Genlight_raw_report_secondaries.png", height = 8, width = 10, units = 'in', res = 500)
gl.report.secondaries(gl)
dev.off()

## heterozygosity ----
png("results/Genlight_raw_report_heterozygosity.png", height = 8, width = 10, units = 'in', res = 500)
gl.report.heterozygosity(gl)
dev.off()



# ---- Filtering ----
# 278 genotypes,  112001 SNPs

## 1 - Individual callrate > 0.50 ----
# Remove individual "SS38 (S_F224)"
# and delete monomorphic loci here, that's why there is a different number of SNPs
t5 <- read_csv("results/Genlight_raw_report_callrate_ind.csv")

ind_to_remove <-
  t5 %>% 
  dplyr::filter(callrate < 0.5) %>% 
  dplyr::pull(id)
  
gl1 <- 
  gl.drop.ind(gl, ind_to_remove, recalc = T, mono.rm = T) 
gl1


## 2 - Callrate by loci > 0.70 ----
gl2 <- 
  gl.filter.callrate(gl1, threshold = 0.70)
gl2

## 3 - Minor Allele Frequency > 0.05 ----
gl3 <- 
  gl.filter.maf(gl2, threshold = 0.05)
gl3


# ## 4 - reproducibility ----
# gl4 <-
#   gl.filter.reproducibility(gl3, threshold = 1)
# gl4
# 
# 
# ## 5 - secondaries ----
# gl5 <-
#   gl.filter.secondaries(gl4)
# gl5



# ---- Save ----
filters <- "callrateind0.50_callrateloci0.70_maf0.05"
# filters <- "callrateind0.50_callrateloci0.70_maf0.05_reprod1"
# filters <- "callrateind0.50_callrateloci0.70_maf0.05_reprod1_secondaries"

gl3 %>% 
  saveRDS(paste0("intermediate/Genlight_Sarpa_salpa_", filters, ".RDS"))


## rm ind outliers...
gl3sub <-
  gl.drop.ind(gl3, c("SS01", "SS02", "SS30", "SS124"), recalc = T, mono.rm = T)

filters <- "callrateind0.50_callrateloci0.70_maf0.05_RMindoutliers"

gl3sub %>% 
  saveRDS(paste0("intermediate/Genlight_Sarpa_salpa_", filters, ".RDS"))

