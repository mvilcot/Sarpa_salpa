source("scripts/0_wrapper.R")


# ---- read data ----
filters <- "callrateind0.50_callrateloci0.70_maf0.05"

genlight <- 
  readRDS(paste0("intermediate/Genlight_Sarpa_salpa_", filters, ".RDS"))



# ---- convert to LEA data format ----

# SNP presence/absence lfmm (package LEA, SilicoDArT)
dartR::gl2geno(genlight,
        outfile = paste0("Geno_Sarpa_salpa_", filters),
        outpath='intermediate')


# ---- SNMF ----
Kmax = 8
nrun = 2

## run snmf ----
obj.snmf <-
  LEA::snmf(paste0("intermediate/Geno_Sarpa_salpa_", filters, ".geno"),
       K = 1:Kmax, alpha = 100, project = "new", repetitions = nrun, entropy = TRUE)

obj.snmf %>%
  saveRDS(paste0("intermediate/snmf_Sarpa_salpa_K1-", Kmax, "_nrun", nrun, ".RDS"))


## barplot ----
obj.snmf <-
  readRDS(paste0("intermediate/snmf_Sarpa_salpa_K1-", Kmax, "_nrun", nrun, ".RDS"))

pdf(paste0("results/snmf_barplot_Sarpa_salpa_K1-", Kmax, "_nrun", nrun, ".pdf"),
    height = 4, width = 10)
plot(obj.snmf, cex = 1.2, pch = 19)
for (K in 2:Kmax){
  # get the cross-entropy of the 10 runs
  ce = LEA::cross.entropy(obj.snmf, K = K)
  # select the run with the lowest cross-entropy for K = 4
  bestrun = which.min(ce)
  # get best qmatrix
  qmatrix <- LEA::Q(obj.snmf, K = K, run = bestrun)


  ## perso plot
  # create an object with membership probabilities
  probs <- as.data.frame(qmatrix)

  # put probabilities in a tibble with IDS and labels for sites
  probs <-
    as.data.frame(qmatrix) %>%
    cbind(id = genlight@ind.names) %>%
    dplyr::as_tibble() 

  # melt into long format
  probs_long <-
    probs %>%
    tidyr::pivot_longer(colnames(qmatrix), names_to = "cluster", values_to = "prob") %>% 
    left_join(data_samples, by = 'id')

  gg <-
    ggplot(probs_long, aes(factor(id), prob, fill = factor(cluster))) +
    geom_col() +
    facet_grid(~location, switch = "x", scales = "free", space = "free") +
    labs(x = "", y = "membership probability") +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_discrete(expand = expansion(add = 0.5)) +
    scale_fill_viridis_d() +
    theme(
      panel.spacing.x = unit(0.5, "lines"), #0.15
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      # panel.grid = element_blank(),
      # panel.background = element_rect(fill = 'white', color = 'white'),
      strip.background = element_rect(fill = "white")) +
    labs(fill = "cluster")

  print(gg)

}
dev.off()







# ---- STRUCTURE ----
# install STUCTURE (non GUI version) from here: https://web.stanford.edu/group/pritchardlab/structure_software/release_versions/v2.3.4/html/structure.html

# change levels to alphabetic, otherwise issue with gl.map.structure
genlight$pop <- factor(as.character(genlight$pop))


## run structure ----
kmin <- 1
kmax <- 8
nrep <- 2
nloci <- genlight$n.loc

structure <- 
  dartR::gl.run.structure(
    genlight,
    k.range = kmin:kmax,
    num.k.rep = nrep,
    exec = "C:/Users/vilcot/Documents/Structure2.3.4/console/structure.exe",
    plot.out = FALSE)

structure %>% saveRDS(paste0("intermediate/STRUCTURE_Sarpa_salpa_", nloci, "loci_K", kmin, "-", kmax, "_r", nrep, ".RDS"))


## output ---- 
structure <- readRDS(paste0("intermediate/STRUCTURE_Sarpa_salpa_", nloci, "loci_K", kmin, "-", kmax, "_r", nrep, ".RDS"))


# evanno plot
ev <- gl.evanno(structure)

evgg <- # arrange plots
  wrap_plots(ev$plots) +
  patchwork::plot_annotation(tag_levels = 'a', # add tags
                             tag_prefix = '(',
                             tag_suffix = ')',
                             tag_sep = "")

evgg <- # add common x axis
  patchwork::wrap_elements(panel = evgg) +
  labs(tag = "K") + 
  theme(plot.tag = element_text(size = rel(1)),
        plot.tag.position = "bottom")

evgg
ggsave(paste0("results/STRUCTURE_evanno_", nloci, "loci_", sites, "_K", kmin, "-", kmax, "_r", nrep, ".png"),
       evgg, 
       height = 6, width = 8)


# automatic ancestry plot
k <- 4
qmat <- gl.plot.structure(structure, K=k, colors_clusters = viridis(k))

# map barplot
gl.map.structure(qmat, genlightLL, K=k, scalex=0.5, scaley=0.2) 



## plot structure ----

### preparing data ----
qmat_plot_large <- 
  as.data.frame(qmat[[1]]) %>% 
  dplyr::select(c("orig.pop", "Label", contains("cluster"))) 
qmat_plot_long <- 
  qmat_plot_large %>% 
  pivot_longer(contains("cluster"), names_to = "cluster", values_to = "ancestry_prop")
qmat_plot_long$orig.pop <- factor(qmat_plot_long$orig.pop, levels = levels(data_sites$site))
qmat_plot_long <- 
  qmat_plot_long %>% 
  dplyr::arrange(orig.pop, Label, cluster, ancestry_prop)



### perso plot ----
gg <- 
  ggplot(qmat_plot_long, aes(Label, ancestry_prop, fill = cluster)) +
  geom_col(width = 1.1) +
  facet_grid(~orig.pop, switch = "x", scales = "free", space = "free") +
  labs(x = "", y = "Ancestry proportion") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expansion(add = 0.15)) +
  scale_fill_viridis_d() +
  theme(
    panel.spacing.x = unit(0.15, "lines"),
    # axis.text.x = element_text(size=2, angle = 90, vjust = 0.5, hjust=1),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    strip.text.x.bottom = element_text(angle = 90, vjust = 0.5, hjust=1),
    panel.background = element_rect(fill = 'white', color = 'white'),
    strip.background = element_rect(fill = 'white', color = "white"))

gg
ggsave(paste0("results/STRUCTURE_barplot_", nloci, "loci_", sites, "_K", k, ".png"),
       gg,
       height = 4, width = 10)




### perso map ----
# # Proportion of individual from each cluster
# qmat_plot_large <-
#   qmat_plot_large %>% 
#   mutate(CLUSTER = ifelse(cluster1 >= 0.8, "cluster1",
#                           ifelse(cluster2 >= 0.8, "cluster2", "mixed")))
# 
# prop_cluster <-
#   qmat_plot_large %>% 
#   dplyr::rename(site = orig.pop) %>% 
#   group_by(site, CLUSTER) %>% 
#   dplyr::summarise(value = n(), .groups = "keep") %>% 
#   pivot_wider(names_from = CLUSTER, values_from = value)
# prop_cluster$site <- ordered(prop_cluster$site, levels=levels(data_sites$site))
# 
# prop_cluster <-
#   prop_cluster %>% 
#   left_join(data_sites, by = "site")
# 
# prop_cluster <- prop_cluster %>% replace(is.na(.), 0)


# OR Proportion of ancestry in each site
prop_cluster <-
  qmat_plot_large %>% 
  dplyr::rename(site = orig.pop) %>% 
  group_by(site) %>% 
  dplyr::summarise(cluster1 = sum(cluster1), 
                   cluster2 = sum(cluster2), 
                   cluster3 = sum(cluster3), 
                   cluster4 = sum(cluster4),
                   cluster5 = sum(cluster5),
                   .groups = "keep") 
prop_cluster$site <- ordered(prop_cluster$site, levels=levels(data_sites$site))

prop_cluster <-
  prop_cluster %>% 
  left_join(data_sites, by = "site")

# load maps --
# full map
map1 <-
  fortify(maps::map(fill=TRUE, plot=FALSE)) %>% 
  as_tibble()

# add map +360
map2 <- 
  map1 %>% 
  mutate(long = long + 360,
         group = group + max(map1$group) + 1)

# crop lon & lat extent
map <- 
  rbind(map1, map2) %>% 
  filter(long > 20  & long < 210 & lat <48 & lat > -48)

library(scatterpie)
gg <- 
  ggplot() +
  ## Countries
  geom_polygon(data = map, aes(x = long, y = lat, group = group), fill="grey20") +
  ## Sites
  geom_scatterpie(data = shift.lon(prop_cluster),
             aes(x = longitude, y = latitude), #r = number_samples/10
             cols = c("cluster1", "cluster2", "cluster3", "cluster4", "cluster5"), # 
             pie_scale = 1.5, color=NA, alpha=.8) + 
  ggrepel::geom_text_repel(data = shift.lon(prop_cluster),
                           aes(x = longitude, y = latitude, color = site, label = site),
                           hjust=0.5, vjust=0, max.overlaps = 5, nudge_x = -10,
                           bg.color = "grey70", bg.r = 0.02) +
  ## Theme
  scale_color_viridis(discrete = T) +
  scale_fill_viridis(discrete = T, direction = 1) +
  theme_minimal() +
  theme(plot.background = element_rect(fill="white", color = "white")) +
  labs(x = "Longitude", y = "Latitude") +
  coord_equal()
  
gg
ggsave(paste0("results/STRUCTURE_mappiechart_", nloci, "loci_", sites, "_K", k, "_propancestry.png"),
       gg,
       height = 8, width = 15)


