# ---- read data ----
source("scripts/0_wrapper.R")

filters <- "callrateind0.50_callrateloci0.70_maf0.05"
filters <- "callrateind0.50_callrateloci0.70_maf0.05_RMindoutliers"
filters <- "callrateind0.50_callrateloci0.70_maf0.05_RMindoutliers_pcadaptQ0.1"
filters <- "callrateind0.50_callrateloci0.70_maf0.05_RMindoutliers_RMpcadaptQ0.1"
# filters <- "callrateind0.50_callrateloci0.70_maf0.05_reprod1"
# filters <- "callrateind0.50_callrateloci0.70_maf0.05_reprod1_secondaries"

genlight <- 
  readRDS(paste0("intermediate/Genlight_Sarpa_salpa_", filters, ".RDS"))

data_samples <- 
  read_csv('intermediate/metadata_samples_sequenced.csv')

genlight


# ---- run PCA ----
npca <- length(genlight$ind.names) - 1

PCA <- adegenet::glPca(genlight, nf = npca)
PCA %>% saveRDS(paste0("intermediate/PCA_output_", filters, "_", npca, "PCA.RDS"))
PCA <- readRDS(paste0("intermediate/PCA_output_", filters, "_", npca, "PCA.RDS"))



# ---- plot perso ----
# get dapc values and metadata
pca_scores <- 
  as.data.frame(PCA$scores) %>% 
  rownames_to_column("id") %>% 
  dplyr::left_join(data_samples, by = "id")

# setup axis
percent_explained <- PCA$eig / sum(PCA$eig) * 100
pretty_pe <- format(round(percent_explained, digits =1), nsmall=1, trim=TRUE)
labels <- c(glue("PC1 ({pretty_pe[1]}%)"),
            glue("PC2 ({pretty_pe[2]}%)"),
            glue("PC3 ({pretty_pe[3]}%)"))

# save  
gg1 <- 
  ggplot(pca_scores, aes(x=PC1, y=PC2, label = id, color = location)) +
  geom_point(size = 1, alpha = 0.5) +
  labs(x=labels[1], y=labels[2]) +
  # scale_color_manual(values = color_perso) +
  theme_light()

gg2 <- 
  ggplot(pca_scores, aes(x=PC1, y=PC3, label = id, color = location)) +
  geom_point(size = 1, alpha = 0.5) +
  labs(x=labels[1], y=labels[3]) +
  # scale_color_manual(values = color_perso) +
  theme_light()

# plotly::ggplotly(gg1)
# plotly::ggplotly(gg2)

gg <- 
  gg1 + gg2 + patchwork::plot_layout(guides = 'collect')
gg

ggsave(paste0("results/PCA_", filters, ".png"),
       gg, height = 4, width = 9)





