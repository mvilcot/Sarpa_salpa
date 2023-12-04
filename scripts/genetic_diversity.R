
# ---- read SNPs dataset ----

# parameters
filters = "missind1_callrate0.70_maf0.05"
level = "site"
sites = "allsites"

# read genlight
genlight <- 
  read.genlight(filters, level,
                site2drop = NULL,
                site2keep = NULL,
                station2drop = NULL,
                station2keep = NULL)

genlight 

# 
# genlight <- 
#   genlight %>% 
#   gl.drop.ind("ECO0919")


# ---- convert to other formats ----

# genind 
genind <- dartR::gl2gi(genlight)

# hierfstat format
ghierfstat <- hierfstat::genind2hierfstat(genind)

# SNP presence/absence lfmm (package LEA, SilicoDArT)
geno <- dartR::gl2geno(genlight)

# matrix
gmatrix <- as.matrix(genlight) # one column by SNP
gmatrixPA <- gmatrix[ , colSums(is.na(gmatrix))==0]
gmatrixPA[gmatrixPA >= 1] <- 1
# gmatrix2 <- as.matrix(genind@tab) # one column by allele
# gmatrix2PA <- gmatrix2[ , colSums(is.na(gmatrix2))==0]
# gmatrix2PA[gmatrix2PA >= 1] <- 1

gmatrixPAloc <- as.data.frame(gmatrixPA)
gmatrixPAloc[[level]] <- genlight@pop
gmatrixPAloc <- 
  aggregate(gmatrixPAloc[,-ncol(gmatrixPAloc)], by=list(location=gmatrixPAloc[[level]]), FUN=max) %>% ##by site
  column_to_rownames("location")



# ---- mean genetic diversity ----

# basic stats
BS <- basic.stats(genind)
BSo <- BS$overall

# Jost additive framework Dst = Ht - Hs
Hs = BSo[["Hs"]]
Ht = BSo[["Ht"]]
Hst = (Ht-Hs)/(1-Hs)
Dst <- Ht - Hs

# Jost multiplicative framework Jst = Jt / Js
Js = 1/(1-Hs)
Jt = 1/(1-Ht)
Jst = Jt/Js # or 1/(1-Hst)

# Hedrick Gst" (Meirmans and Hedrick 2011)
GstPP.hed <- mmod::Gst_Hedrick(genind)

# Jost D (Jost 2008)
JostD <- mmod::D_Jost(genind)

# Fst population-specific (Weir and Goudet 2017)
popFst <- hierfstat::betas(ghierfstat, nboot=1000)

# # Confidence interval
# bs <- chao_bootstrap(genind)
# bs_D <- summarise_bootstrap(bs, JostD)
# bias <- bs.D$summary.global.het[1]- obs.D$global.het
# bs.D$summary.global.het- bias

# Jaccard
### either pairwise by individual, then average
# jac <- beta.multi(gmatrixPA, index.family="jaccard")
# jac2 <- beta.pair(gmatrix2PA, index.family="jaccard")
# 
# jacSIM <- generate_pw_jaccard(geno = gmatrix, 
#                            pop.label = genlight$pop,
#                            plot_it = FALSE)

### or PA by site first, then pairwise jaccard
jac_multi <- betapart::beta.multi(gmatrixPAloc, index.family="jaccard")

# create global table
gd_global <-
  cbind(Ho = BSo["Ho"], Hs, Ht, Hst, Dst,
        Js, Jt, Jst,
        Fst = BSo["Fst"], 
        Fis = BSo["Fis"],
        Dstp = BSo["Dstp"], 
        Dest = BSo["Dest"],
        Gstpp.hed = GstPP.hed$global,
        D.Jost = JostD$global.het,
        popFst.WG = popFst$betaW,
        jtu = jac_multi$beta.JTU,
        jac = jac_multi$beta.JAC,
        jne = jac_multi$beta.JNE) %>%
  dplyr::as_tibble()




# ---- alpha gd by location ----
## ---- by hand ----
gd_alpha <-
  data.frame(Hs = colMeans(BS$Hs, na.rm = T),
             Ho = colMeans(BS$Ho, na.rm = T)) %>% 
  tibble::rownames_to_column(level) %>% 
  dplyr::left_join(
    data.frame(popFst.WG = popFst$betaiovl) %>% 
      tibble::rownames_to_column(level), 
    by = level)


# ## ---- with dartR ----
# ## !!!!!! Similar values, except for Christmas Island... issue ECO019 Ho very high ####
# alpha_pop <-
#   gl.report.heterozygosity(genlight, method = "pop") %>% 
#   dplyr::rename(site = pop) %>% 
#   dplyr::left_join(
#     data.frame(popFst.WG = popFst$betaiovl) %>% 
#       tibble::rownames_to_column(level), 
#     by = level)
# 
# alpha_ind <-
#   gl.report.heterozygosity(genlight, method = "ind") 
# 
# alpha_pop %>% write_csv(paste0("results/1_genetic_diversity/gd_table_", level, "_dartR_pop.csv"))
# alpha_ind %>% write_csv(paste0("results/1_genetic_diversity/gd_table_dartR_ind.csv"))



# ---- beta gd pairwise ----
# Fst
Fst_pair <- hierfstat::genet.dist(genind, method = "WC84")

# Hedrick G"st
GstPP.hed_pair <- mmod::pairwise_Gst_Hedrick(genind)

# Jost D
JostD_pair <- mmod::pairwise_D(genind)

# Jaccard
### either PA by site first, then pairwise jaccard
jac_pair <- betapart::beta.pair(gmatrixPAloc, index.family="jaccard")

### or pairwise by individual, then average
jac_pair <- betapart::beta.pair(gmatrixPA, index.family="jaccard")


# put into list
list_gd_beta_pair <-
  list(Fst = Fst_pair, 
       GstPP.hed = GstPP.hed_pair,
       D.Jost = JostD_pair,
       jtu = jac_pair$beta.jtu,
       jac = jac_pair$beta.jac,
       jne = jac_pair$beta.jne)


# ---- export ----
BS %>% saveRDS(paste0("intermediate/1_genetic_diversity/basic_stats_", level, ".RDS"))
gd_global %>% write.csv(paste0("results/1_genetic_diversity/gd_table_global_", level, ".csv"), row.names = F, quote = F)
gd_alpha %>% write.csv(paste0("results/1_genetic_diversity/gd_table_", level, ".csv"), row.names = F, quote = F)
list_gd_beta_pair %>% saveRDS(paste0("results/1_genetic_diversity/gd_list_pairwise_", level, ".RDS"))



# ---- *** DRAFTS ----

## *** Pairwise (diveRsity) ----
# test <- diffCalc(infile = "Intermediate/GenepopRadiator_DartSeq_Etelis_coruscans_grouped_missind_callrate0.70_maf0.05_sites_noCocos.gen", 
#          outfile = "test", fst = T, pairwise = T)
# 
# plot(test$pairwise$D ~ test$pairwise$GGst)
# 
# 
# test2 <- fastDivPart(infile = "Intermediate/GenepopRadiator_DartSeq_Etelis_coruscans_grouped_missind_callrate0.70_maf0.05_sites_noCocos.gen", outfile = "myresults", gp = 3, 
#                  pairwise = T, plot = TRUE)



## *** Jaccard ---- 
# ## HierJd 
# GDbetaJAC <- HierJd("Intermediate/test.gen", ncode = 3, r = 1, nreg = 1)
# saveRDS(GDbetaJAC, "Results/test_HierJd_betaGD_stations_noinf2.RDS")
# # !!!!!!!!!!!!!!!!! Error in HierJd("Intermediate/test.gen", ncode = 3) : !!!!!!!!!!!!
# # argument "r" is missing, with no default
# # In addition: There were 27 warnings (use warnings() to see them)
# 
# 
# 
## TEST Jaccard Betapart
# PAalleles <- as.matrix(as.data.frame(genind@tab))
# PAallelesTEST <- PAalleles[ , colSums(is.na(PAalleles))==0]
# PAallelesTEST[PAallelesTEST > 0] <- 1
# test <- beta.pair(PAallelesTEST, index.family="jaccard")
# 
# PAallelesTEST$site <- genind@pop
# PAallelesTESTagg <- aggregate(PAallelesTEST[,-ncol(PAallelesTEST)], by=list(Site=PAallelesTEST$site), FUN=max) ##by site
# test <- beta.multi(PAallelesTEST[,-ncol(PAallelesTEST)], index.family="jaccard")
# 
# 
# Jaccard
# library(HierDpart)
# JacHier <- HierJd("_archive/-36_radiator_genomic_converter_20221115@1129/radiator_data_20221115@1129_genepop.gen",
#                   ncode=3, nreg=1, r=1)
# print(JacHier)
# 
# 
# 
# 
# 
## *** Compute mean alpha by site ---- 
# gd_alpha_site2 <-
#   as.data.frame(BS$Hs) %>%
#   rownames_to_column("loci") %>%
#   pivot_longer(-loci, names_to = "site") %>%
#   group_by(site) %>%
#   summarise(Hs=mean(value, na.rm=T), .groups = "keep")
# 
# gd_alpha_site3 <-
#   data.frame(Hs = apply(BS$Hs, MARGIN = 2, FUN = mean, na.rm = T)) # mean by population
# 
# 
# 
# 
## *** convert to other formats ---- 
# # Genpop class
# genpop <- adegenet::genind2genpop(genind)
# genpop
# 
# # Genepop (Rousset)
# genomic_converter(genlight, output = "genepop")
# genepop <- dartR::gl2genepop(genlight, outfile = paste0('Intermediate/Genepop_Etelis_coruscans_ordered_', filters, '.genepop.gen', outpath = "."))
# genepop
# 
# # data frame format
# gdf <- genind2df(genind)
# gdf2 <- 
#   gdf %>% 
#   select(-pop)
# 
# alleles <- 
#   read.table(text = gsub("/", "-", colnames(gdf2)), sep = "-", as.is = TRUE) %>%# get base variant for each SNP
#   select(3, 4) %>%
#   cbind(colnames(gdf2), .) %>%
#   dplyr::rename(name = 1, var1 = 2, var2 = 3)
#   
# for (i in 1:nrow(alleles)){
#   allele1 <- alleles[i, "var1"]
#   allele2 <- alleles[i, "var2"]
#   
#   gdf2[!is.na(gdf[,i]) & gdf2[,i] == paste0(allele1, allele1), i] <- 0
#   gdf2[!is.na(gdf[,i]) & gdf2[,i] == paste0(allele1, allele2), i] <- 1
#   gdf2[!is.na(gdf[,i]) & gdf2[,i] == paste0(allele2, allele1), i] <- 1
#   gdf2[!is.na(gdf[,i]) & gdf2[,i] == paste0(allele2, allele2), i] <- 2
#   
# }
# test <- as.data.frame(sapply(gdf2, as.numeric))
# rownames(test) <- rownames(gdf2)
# 
# 
# 
# # plink
# genlight$chromosome <- as.factor("1")
# gl2plink(genlight, bed_file = F, 
#          outpath = getwd(), outfile = "test.plink")



# ---- load ----
level = "site"

gd_alpha <- read.csv(paste0("results/1_genetic_diversity/gd_table_", level, ".csv"))
gd_beta <- readRDS(paste0("results/1_genetic_diversity/gd_list_pairwise_", level, ".RDS"))
BS <- readRDS(paste0("intermediate/1_genetic_diversity/basic_stats_", level, ".RDS"))


# ---- alpha GD by location ----

## Hs by loci ----
gd_alpha_loci <- 
  dplyr::as_tibble(BS[["Hs"]]) %>% 
  tidyr::pivot_longer(cols = everything(), 
                      names_to = level, 
                      values_to = "Hs")

# relevel sites
gd_alpha_loci[[level]] <- factor(gd_alpha_loci[[level]], levels = levels(data_samples[[level]]))

# plot
ggplot(gd_alpha_loci, aes(x=.data[[level]], y=.data[["Hs"]], fill = .data[[level]])) + 
  geom_boxplot() +
  scale_fill_manual(values = color_perso) +
  theme_light() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 30, vjust = 1, hjust=1))

# save
ggsave(paste0("results/1_genetic_diversity/boxplot_alpha_Hs_", level, "_loci_hierfstat.png"),
       width = 8, height = 6)



## Hs by ind dartR ----
alpha_ind <- read_csv(paste0("results/1_genetic_diversity/gd_table_ind_dartR.csv"))

# add metadata
alpha_ind <-
  alpha_ind %>% 
  dplyr::rename(id = ind.name) %>% 
  full_join(data_samples, by= "id")

# plot
ggA <- 
  ggplot(alpha_ind, aes(x=.data[[level]], y=.data[["Ho"]], fill = .data[[level]])) + 
  geom_boxplot() +
  scale_fill_manual(values = color_perso) +
  theme_light() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 30, vjust = 1, hjust=1))

ggsave(paste0("results/1_genetic_diversity/boxplot_alpha_Ho_", level, "_ind_dartR.png"),
       width = 8, height = 6)


# ---- Hs ~ Fst pop specific ----

# Pearson, lm
corP <- cor.test(gd_alpha$Hs, gd_alpha$popFst.WG)
# model <- summary(lm(Hs ~ popFst.WG, data = gd_alpha))

# relevel sites
gd_alpha[[level]] <- factor(gd_alpha[[level]], levels = levels(data_samples[[level]]))

# plot
ggplot(gd_alpha, aes(x=popFst.WG, y=Hs, color = .data[[level]])) + 
  geom_smooth(method='lm', formula=y~x, color = "grey35") +
  geom_point(size = 4, alpha = 0.8) +
  geom_text_repel(data = gd_alpha[gd_alpha$site %in% c("Seychelles","Christmas_Island","Hawaii","W_Australia"),],
                  aes(label=.data[[level]]), box.padding = 0.5, size = 3.5,
                  bg.color = "grey70", bg.r = 0.02) +
  scale_color_manual(values = color_perso) +
  annotate('text', x=max(gd_alpha$popFst.WG), y=max(gd_alpha$Hs), 
           hjust=1, vjust=0.5, size=3,
           label=paste0("r = ", signif(corP$estimate, 3), "\n p-value = ", signif(corP$p.value, 1))) + 
  theme_light()

# save
ggsave(paste0("results/1_genetic_diversity/plot_Hs_Fst_pop_specific_", level, "_light.png"),
       width = 7, height = 5)



# Hs ~ longitude ----
## load
df <-
  data_sites %>% 
  left_join(gd_alpha)

## relevel
df$site <- factor(df$site, levels = unique(df$site))
df <- shift.lon(df)# Pearson, lm

## correlation
corP <- cor.test(df$longitude, df$Ho)

## plot
ggplot(df, aes(x=longitude, y=Ho, color = .data[[level]])) + 
  # geom_smooth(method='lm', formula=y~x, colour = "grey35") +
  geom_point() +
  scale_color_manual(values = color_perso) +
  theme_light()+
  # geom_text_repel(aes(label=.data[[level]]), ) +
  annotate('text', x=max(df$longitude), y=max(df$Ho), 
           hjust=1, vjust=1, size=4.5,
           label=paste0("Pearson = ", round(corP$estimate, 4),
                        "\n p = ", round(corP$p.value, 4)))


# save
ggsave(paste0("results/1_genetic_diversity/alpha_gd_Ho_longitude_site.png"),
       width = 8, height = 6)




# ---- comparison beta GD metrics ----
list_GDbeta <- list()

for (metricGD in names(gd_beta)){
  
  # get distance matrix
  mat_GDbeta <- as.matrix(gd_beta[[metricGD]])
  
  # order rows alphabetically
  mat_GDbeta <- mat_GDbeta[order(rownames(mat_GDbeta)), order(colnames(mat_GDbeta))]
  
  # pivot longer distance matrix
  melt_GDbeta <- melt.dist(dist = mat_GDbeta, metric = metricGD)
  
  # put into list
  list_GDbeta[[metricGD]] <- melt_GDbeta
  
}

# merge two distance matrix into one df
merge_gd_beta <- 
  plyr::join_all(list_GDbeta,
                 by = c(paste0(level, "1"), paste0(level, "2")))

merge_gd_beta$site <- paste(merge_gd_beta$site1,
                            merge_gd_beta$site2,
                            sep = "-")

merge_gd_beta$Christmas <- grepl("Christmas", merge_gd_beta$site)


ggplot(merge_gd_beta, aes(Fst, GstPP.hed, color = Christmas)) +
  geom_point() 

ggplot(merge_gd_beta, aes(D.Jost, GstPP.hed, color = Christmas)) +
  geom_point()

ggplot(merge_gd_beta, aes(Fst, D.Jost, color = Christmas)) +
  geom_point()


# library(plotly)
# plot_ly(x=merge_gd_beta[["GstPP.hed"]], 
#         y=merge_gd_beta[["JostD"]], 
#         z=merge_gd_beta[["Fst"]], 
#         type="scatter3d", mode="markers", size = 0.5)
# 
#  
# 
# library(corrplot)
# temp <- cor(merge_gd_beta[, names(gd_beta)])
# corrplot(temp)
# 


# Hs ~ dist to CT ----
## !!!!!!!! TO DO !!!!!!!!!!!!!!!! ----



