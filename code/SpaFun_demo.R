
Rcpp::sourceCpp("./rkhscov2_gauss.cpp")
source("./rkhsprep.R")


load('../data/human_prostate_cancer_ffpe.RData')
annotations <- read.csv('../data/prostate_cancer_annotation_region_for_FPCA.csv')


library(ggplot2)
# plot results at single cell level
df <- data.frame(x = annotations$x, y = annotations$y, domain = rowSums(count))
ggplot(df, aes(x = x, y = y, color = domain)) +           
  geom_point() + coord_fixed(ratio = 1)


# plot results at single cell level
library(ggplot2)
df <- data.frame(x = annotations$x, y = annotations$y, domain = annotations$area)
ggplot(df, aes(x = x, y = y, color = domain)) +           
  geom_point() + coord_fixed(ratio = 1) + 
  scale_color_manual(values = c('domain_1' = 'orange', 'domain_2' = 'darkgreen', 'domain_3' = 'blue', 'domain_4' = 'purple', 
                                'domain_5' = 'purple', 'domain_6' = 'darkgreen', 'other' = 'white' ))+
  theme_classic() + labs(color = "") + 
  theme(axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        plot.title=element_text(hjust = 0.5),
        legend.position="none",
        text = element_text(size=15))


#count_all <- count
loc_all <- loc

colnames(loc_all) <- c('x', 'y')
count_all <- count
colnames(count_all) <- gene_name

layer_name <- setdiff(unique(annotations$area), c('other'))

for (l in layer_name[c(2)]){

  load('human_prostate_cancer_ffpe.RData')
  count <- count[annotations$area == l, ]
  loc <- loc_all[annotations$area == l, ]
  print(dim(loc))
  
  re_filter <- filter_count(count, loc, min_percentage = 0.3, min_total = 1)
  loc_f <- re_filter[[1]]
  count_f <- re_filter[[2]]
  
  loc <- loc_f
  
  gene_expression_levels <- normalization22(count_f)
  print(dim(gene_expression_levels))
  
  remove(count)
  remove(count_f)
  remove(re_filter)
  
  library(parallel)
  n = ncol(gene_expression_levels)
  m = nrow(gene_expression_levels)
  
  means = rowMeans(gene_expression_levels)
  
  #library(rkhscovfun)
  Xs <- as.vector(gene_expression_levels)
  Xs <- Xs - rep(means, n)
  
  Sub = rep(1:n, each = m)
  
  
  S1 <- (loc[,1] - min(loc[,1]))/(max(loc[,1]) - min(loc[,1]))
  S2 <- (loc[,2] - min(loc[,2]))/(max(loc[,2]) - min(loc[,2]))
  
  dis2 = get_distance2(S1, S2)
  gamma = median(sqrt(dis2))
  print(gamma)
  Kernel = "gauss"
  
  CV = rkhscov.cv(S1, S2, Xs, Sub, control.alg=list(preptol = 1e-05, gamma = gamma))
  res = CV$sobj
  
  #### fPCA ####
  fres <- fpca.rkhscov(res, Kernel = "gauss", gamma = gamma)
  efun <- compute.fpca(S1, S2, fres, Kernel = "gauss", gamma = gamma) # eigenfunctions
  
  save(efun, fres, gamma, file = paste0('prostate_annotation_fpca_result_',l, '.RData'))

  print(l)
}


########################################
############ DRG Generation ############
########################################

colnames(count) <- gene_name

source("/Users/yanghongguo/Desktop/Research/SpaFun/Figures/code/helper_function.R")

# domain_1 strome
# domain_2 tumor

l <- 'domain_2'
count_d <- count[annotations$area == l, ]
loc_d <- loc[annotations$area == l, ]
# filter genes
re_filter <- filter_count(count_d, loc_d, min_percentage = 0.3, min_total = 1)
loc_f <- re_filter[[1]]
count_f <- re_filter[[2]]

# normalization
gene_expression_levels <- normalization22(count_f)
print(dim(gene_expression_levels))

load(paste0('/Users/yanghongguo/Desktop/Research/SpaFun/From Xi Jiang/results/prostate_annotation_fpca_result_', l, '.RData'))

# count gene numbers
########################
# gene_list
n_gene <- ncol(count_f)

p_values <- matrix(0, nrow = n_gene, ncol = 10)
rownames(p_values) <- colnames(count_f)
beta <- matrix(0, nrow = n_gene, ncol = 10)
rownames(beta) <- colnames(count_f)

for (p in 1:10){
  p_temp <- numeric(n_gene)
  for (i in 1:n_gene){
    rr <- cor.test(gene_expression_levels[, i], efun[p, ])
    beta[i, p] <- rr$estimate
    p_temp[i] <- rr$p.value
  }
  p_values[, p] <- p.adjust(p_temp, 'BH')
}

# save for sending to Lin Xu
adj_p <- p_values[,1:5]
colnames(adj_p) <- c("PC1", "PC2", "PC3", "PC4", "PC5")
write.csv(adj_p, paste0('/Users/yanghongguo/Desktop/prostate_cancer_', l,'_region_adj_pvalue_SpaFun.csv'))


correlation <- beta[,1:5]
colnames(correlation) <- c("PC1", "PC2", "PC3", "PC4", "PC5")
write.csv(correlation, paste0('/Users/yanghongguo/Desktop/prostate_cancer_', l,'_region_corr_coef_SpaFun.csv'))




################# Plot the Results #################

library(ggplot2)

plot_eigen_function <- function(anno, i, domain){
  pc_tumor_domain <- data.frame(x = anno$x, y = anno$y, pc1 = NA)
  pc_tumor_domain$pc1[anno$area == domain] <- -efun[i, ]
  df_clean <- na.omit(pc_tumor_domain)
  filename = paste0("../figures/prostate/PC_", i,  ".png")
  g <- ggplot(df_clean) + geom_point(mapping = aes(x = x, y = y, color = pc1), size = 1.5) + 
    coord_fixed(ratio = 1) + scale_color_gradient2(low = "black", high = "red", mid = 'blue', na.value = "grey") +
    theme_classic() + labs(color = "") + 
    theme(axis.line=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          plot.title=element_text(hjust = 0.5),
          legend.position="none",
          text = element_text(size=15),
          plot.margin = unit(c(0, 0, 0, 0), "lines"))
  ggsave(filename = filename, plot = g, width = 4, height = 4, dpi = 300)
}


for (i in c(1,2,3)){
  plot_eigen_function(anno, i, l)
}


plot_individual_gene <- function(anno, gene, domain){
  pc_tumor_domain <- data.frame(x = anno$x, y = anno$y, pc1 = NA)
  pc_tumor_domain$pc1[anno$area == domain] <- gene_expression_levels[, gene]
  df_clean <- na.omit(pc_tumor_domain)
  filename = paste0("../figure/prostate/", domain, "_", gene,  ".png")
  
  g <- ggplot(df_clean) + geom_point(mapping = aes(x = x, y = y, color = pc1), size = 1.5) + 
    coord_fixed(ratio = 1) + scale_color_gradient2(low = "black", high = "red", mid = 'blue', na.value = "grey") +
    theme_classic() + labs(color = "") + 
    theme(axis.line=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          plot.title=element_text(hjust = 0.5),
          legend.position="none",
          text = element_text(size=15),
          plot.margin = unit(c(0, 0, 0, 0), "lines"))
  ggsave(filename = filename, plot = g, width = 4, height = 4, dpi = 300)
}



# domain 2 (tumor)
genes <- c("TMEFF2", "GOLM1",  "ADGRF1",  "NPY",
           "TAGLN",  "ACTA2",  "FLNA", "MYL9",
           "AMACR",  "CRISP3",  "APOD",  "GOLM1",
           "NR4A1",  "FOS",  "TRIB1",  "FOSB",
           "ACTA2",  "TAGLN",  "FLNA",  "ACTG2",
           "GOLM1",  "TMEFF2",  "ADGRF1",  "OR51E2"
)

for (gene in genes){
  plot_individual_gene(anno, gene, l)
}
