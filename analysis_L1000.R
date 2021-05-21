#######################################################################
#### Overview - Template
#######################################################################

## Goal: 
##       

## Data structures:
## -

#######################################################################
#### SET UP
#######################################################################

#######################################################################
#### Set relevant input and output directories
#######################################################################

## Inputs
L1000dir <- "../data/L1000"

## Outputs
dataDir <- "./data"
resultsDir <- "../results/AALE_KRAS-RIT1"
figsDir <- "../results/Figures"
qcDir <- "./results/QC"


#######################################################################
#### Functions
#######################################################################

source("../scripts/R/load.R")

#######################################################################
#### Load Data
#######################################################################

AALE <- fread(paste0(L1000dir, "/S3_zscores.csv"))
AALE_rows <- which(names(AALE) == "AALE")

AALE_raw <- fread(paste0(L1000dir, "/S3_Normalized.csv"))

#######################################################################
#### Clean AALE z raw normalized expression data -> tidy table
#######################################################################

allele_ids <- as_tibble(t(AALE[1:2, ..AALE_rows]))
names(allele_ids) <- c("Clone_ID", "Allele")
AALE_raw <- as_tibble(AALE_raw)

idx <- c(3, str_which(names(AALE_raw), "AALE"))
AALE_raw <- AALE_raw[,idx]

Sample_IDs <- names(AALE_raw)[-1]
Clone_IDs <- t(AALE_raw[1,-1])
sample_to_clone_IDs <- data.frame(Sample_ID = Sample_IDs,
                                  Clone_ID = Clone_IDs)

AALE_raw <- AALE_raw[-1,]
names(AALE_raw)[1] <- "gene_symbol"

df_AALE_raw <- gather(AALE_raw, "Sample_ID", "expression", -"gene_symbol")
df_AALE_raw <- df_AALE_raw %>%
    full_join(sample_to_clone_IDs, by = "Sample_ID") %>%
    full_join(allele_ids, by = "Clone_ID")

df_AALE_raw <- as_tibble(df_AALE_raw)
df_AALE_raw$expression <- as.numeric(df_AALE_raw$expression)

#######################################################################
#### Clean AALE z scores data
#######################################################################

gene_ids <- AALE[-c(1,2),3]

AALE_rows <- which(names(AALE) == "AALE")
AALE <- AALE[, ..AALE_rows]
AALE <- AALE[-1,]
colnames(AALE) <- as.character(AALE[1,])
AALE <- AALE[-1,]

## ?
for (i in seq(1, ncol(AALE))) {
    AALE[,..i] <- as.numeric(AALE[,..i])
}

#######################################################################
#### Define controls and genes of interest
#######################################################################

oncogenes <- c("ARAF", "BRAF", "EGFR", "MET", "MDM2", "MYC", "NRAS", "RIT1", "KRAS", "TP53", "AURKA", "CTNNB1", "CCND1", "NFE2L2", "PIK3CA")
controls <- c("BFP", "EGFP", "HCRED", "LUCIFERASE")



gene_plot_order <- c("Control", "KRAS", "RIT1", "NRAS", "EGFR", "MET", "ARAF", "BRAF", "MYC", "TP53", "MDM2", "CTNNB1", "AURKA", "CCND1", "PIK3CA", "NFE2L2")
df_plot$gene <- factor(df_plot$gene, levels = gene_plot_order)



#######################################################################
#### Pearson correlation heatmap of normalized gene expression profiles
#######################################################################

## Replicate collapsed (median gene expression)

df <- df_AALE_raw %>%
    group_by(gene_symbol, Allele) %>%
    summarise(avg_expression = median(expression)) %>%
    ungroup() %>%
    mutate(Gene = str_extract(Allele, "^[A-Z0-9]+")) %>%
    filter(Gene %in% c(controls, oncogenes)) %>%
    dplyr::select(-Gene) %>%
    spread(Allele, avg_expression)

df_corr <- cor(df[,-1], method = "pearson")
rownames(df_corr) <- gsub("_", " ", rownames(df_corr))
colnames(df_corr) <- gsub("_", " ", colnames(df_corr))

## hirearchical clustering
hc <- hclust(dist(t(df_corr)))
pdf(file = paste0(resultsDir, "/L1000_AALE_pearson_hclust.pdf"),
    height = 10, width = 16)
plot(hc)
dev.off()
order <- hc$order

## geom_tile heatmap
df_plot <- as_tibble(df_corr) %>%
    mutate(allele_x = rownames(df_corr)) %>%
    gather("allele_y", "pearson", -"allele_x")

df_plot$allele_x <- factor(df_plot$allele_x, levels = rownames(df_corr)[order])
df_plot$allele_y <- factor(df_plot$allele_y, levels = rev(rownames(df_corr)[order]))

g <- ggplot(data = df_plot) +
    geom_tile(aes(x = allele_x, y = allele_y, fill = pearson)) +
    scale_fill_gradient(low = "yellow", high = "red", limits = c(0.7, 1)) +
    scale_x_discrete(position = "top") +
    theme_bw() +
    theme(axis.text.y = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5),
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.ticks.x = element_blank())
dhc <- as.dendrogram(hc)
df_plotclust <- dendro_data(dhc, type = "rectangle")
g_clust <- ggplot(segment(df_plotclust)) +
    geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
    scale_x_continuous(expand = c(0,0)) +
    theme_dendro()
p <- g_clust / g +
    plot_layout(heights = c(1,4))
ggsave(p, filename = "L1000_AALE_pearson_correlations.pdf",
       height = 14, width = 13,
       device = "pdf",
       path = resultsDir)


#######################################################################
#### PCA of normalized gene expression
#######################################################################

## Separate replicates

df <- df_AALE_raw %>%
    dplyr::select(gene_symbol, Sample_ID, expression) %>%
    spread(Sample_ID, expression)

pca <- prcomp(df[,-1])
summary(pca)
df_pca <- data.frame(pca$rotation)

df_pca$Sample_ID <- rownames(df_pca)
df_pca <- df_pca %>%
    full_join(sample_to_clone_IDs, by = "Sample_ID") %>%
    inner_join(allele_ids, by = "Clone_ID")

df_pca$gene <- str_extract(df_pca$Allele, "^[A-Z0-9]+")
df_pca$is_WT <- str_extract(df_pca$Allele, "[A-Za-z0-9]+$")
df_pca <- filter(df_pca, is_WT != "V5")
df_pca$is_WT <- df_pca$is_WT %in% c("WT", controls)

oncogene_subset <- c("KRAS", "RIT1", "ARAF", "EGFR", "MYC")
df_plot <- df_pca %>%
    filter(gene %in% oncogene_subset)
df_controls <- df_pca %>%
    filter(gene %in% controls) %>%
    mutate(gene = "Control")
df_plot <- bind_rows(df_plot, df_controls)

control_color <- "#888888"
names(control_color) <- c("Control")
oncogene_colors <- rev(brewer.pal(length(oncogene_subset), "Set2"))
names(oncogene_colors) <- oncogene_subset
plot_colors <- c(oncogene_colors, control_color)
plot_colors_darken <- darken(plot_colors, 0.2)
names(plot_colors_darken) <- names(plot_colors)

g <- ggplot(data = df_plot, aes(PC1, PC2)) +
    geom_point(aes(fill = gene, color = gene, shape = is_WT),
               size = 3.5, alpha = 0.7) +
    scale_fill_manual(values = plot_colors) +
    scale_color_manual(values = plot_colors_darken) +     
    scale_shape_manual(values = c("TRUE" = 21, "FALSE" = 23),
                       labels = c("TRUE" = "WT", "FALSE" = "Variant"),
                       name = NULL) +
    guides() +
   theme_classic(base_size = 24)    
ggsave(g, filename = "L1000_AALE_PCA_replicate_expression.pdf",
       height = 7, width = 10,
       device = "pdf",
       path = resultsDir)

## Replicate collapsed (median gene expression)

df <- df_AALE_raw %>%
    group_by(gene_symbol, Allele) %>%
    summarise(avg_expression = median(expression)) %>%
    ungroup() %>%
    mutate(Gene = str_extract(Allele, "^[A-Z0-9]+")) %>%
    filter(Gene %in% c(controls, oncogenes)) %>%
    dplyr::select(-Gene) %>%
    spread(Allele, avg_expression)
    
pca <- prcomp(df[,c(-1,-2)])
summary(pca)
df_pca <- data.frame(pca$rotation)
df_pca$allele <- rownames(df_pca)

df_pca$gene <- str_extract(df_pca$allele, "^[A-Z0-9]+")
df_pca$is_WT <- str_extract(df_pca$allele, "[A-Za-z0-9]+$")
df_pca <- filter(df_pca, is_WT != "V5")
df_pca$is_WT <- df_pca$is_WT %in% c("WT", controls)

df_plot <- df_pca %>%
    filter(gene %in% oncogenes)
df_controls <- df_pca %>%
    filter(gene %in% controls) %>%
    mutate(gene = "Control")
df_plot <- bind_rows(df_plot, df_controls)

control_color <- "#333333"
names(control_color) <- c("Control")
oncogene_colors <- rev(brewer.pal(length(oncogenes), "Set2"))
oncogene_colors <- c(oncogene_colors, rev(brewer.pal(length(oncogenes)-8, "Set1")))
names(oncogene_colors) <- oncogenes
plot_colors <- c(oncogene_colors, control_color)
plot_colors_darken <- darken(plot_colors, 0.2)
names(plot_colors_darken) <- names(plot_colors)

# Plot PC1 and PC2
g <- ggplot(data = df_plot, aes(PC1, PC2)) +
    geom_text_repel(data = filter(df_plot, gene %in% c("KRAS", "RIT1") & PC2>0),
                    aes(label = allele),
                    ylim = c(0.2, NA),
                    angle = 90) +
    geom_text_repel(data = filter(df_plot, gene %in% c("KRAS", "RIT1") & PC2<0),
                    aes(label = allele),
                    ylim = c(NA, -0.15),
                    angle = -90) +
    geom_point(aes(fill = gene, color = gene, shape = is_WT),
               size = 3.5, alpha = 0.7) +
    scale_fill_manual(values = plot_colors,
                      name = NULL) +
    scale_color_manual(values = plot_colors_darken,
                       name = NULL) +
    scale_shape_manual(values = c("TRUE" = 21, "FALSE" = 23),
                       labels = c("TRUE" = "WT", "FALSE" = "Variant"),
                       name = NULL) +
    scale_y_continuous(limits = c(-0.25, 0.3)) +
    theme_classic(base_size = 24)
ggsave(g, filename = "L1000_AALE_PCA_expression.pdf",
       height = 7, width = 10,
       device = "pdf",
       path = resultsDir)


#######################################################################
#### PCA of z-scores
#######################################################################

mat <- sapply(AALE, as.numeric)
dim(mat)
    
pca <- prcomp(mat)
df_pca <- data.frame(pca$rotation)
dim(df_pca)
df_pca$allele <- names(AALE)


df_pca$gene <- str_extract(df_pca$allele, "^[A-Z0-9]+")
df_pca$is_WT <- str_extract(df_pca$allele, "[A-Za-z0-9]+$")
df_pca <- filter(df_pca, is_WT != "V5")
df_pca$is_WT <- df_pca$is_WT %in% c("WT", controls)

df_plot <- df_pca %>%
    filter(gene %in% oncogenes)
df_controls <- df_pca %>%
    filter(gene %in% controls) %>%
    mutate(gene = "Control")
df_plot <- bind_rows(df_plot, df_controls)

control_color <- "#333333"
names(control_color) <- c("Control")
oncogene_colors <- rev(brewer.pal(length(oncogenes), "Set2"))
oncogene_colors <- c(oncogene_colors, rev(brewer.pal(length(oncogenes)-8, "Set1")))
names(oncogene_colors) <- oncogenes
plot_colors <- c(oncogene_colors, control_color)
plot_colors_darken <- darken(plot_colors, 0.2)
names(plot_colors_darken) <- names(plot_colors)

gene_plot_order <- c("Control", "KRAS", "RIT1", "NRAS", "EGFR", "MET", "ARAF", "BRAF", "MYC", "TP53", "MDM2", "CTNNB1", "AURKA", "CCND1", "PIK3CA", "NFE2L2")
df_plot$gene <- factor(df_plot$gene, levels = gene_plot_order)

# Plot PC1 and PC2 (12.67% and 5.67% variance explained)
g <- ggplot(data = df_plot, aes(PC1, PC2)) +
    geom_point(aes(fill = gene, color = gene, shape = is_WT),
               size = 3.5, alpha = 0.7) +
    scale_fill_manual(values = plot_colors,
                      name = NULL) +
    scale_color_manual(values = plot_colors_darken,
                       name = NULL) +
    scale_shape_manual(values = c("TRUE" = 21, "FALSE" = 23),
                       labels = c("TRUE" = "WT", "FALSE" = "Variant"),
                       name = NULL) +
    theme_classic(base_size = 24)
ggsave(g, filename = "L1000_AALE_PCA.pdf",
       height = 7, width = 10,
       device = "pdf",
       path = resultsDir)

#######################################################################
## Generate df_AALE with z-scores
#######################################################################
AALE$Geneid <- gene_ids
df_AALE <- tibble(AALE) %>%
    gather("perturbation", "z", -Geneid)
d <- as_tibble(str_split_fixed(df_AALE$perturbation, "_", 2))
colnames(d) <- c("perturb_gene", "perturb_variant")
df_AALE <- bind_cols(d, df_AALE) %>%
    mutate(z = as.numeric(z))

#######################################################################    
## Plot individual gene expression
#######################################################################

"FOSL1" %in% gene_ids$pr_gene_symbol
    
df_plot <- df_AALE %>%
    filter(Geneid == "FOSL1") %>%
    filter(perturb_gene %in% oncogenes) %>%
    filter(!(perturb_variant == "WT_V5"))
df_controls <- df_AALE %>%
    filter(Geneid == "FOSL1") %>%
    filter(perturb_gene %in% controls) %>%
    mutate(perturb_variant = perturb_gene) %>%
    mutate(perturb_gene = "Control")
df_plot <- bind_rows(df_plot, df_controls)

df_plot$pos <- as.numeric(gsub("p.[A-Z]*([0-9]+)[^0-9].*", "\\1", df_plot$perturb_variant))

df_plot <- df_plot %>%
    group_by(perturb_gene) %>%
    arrange(pos) %>%
    ungroup()

perturbations <- unique(df_plot$perturbation)
perturbation_plot_order <- c(controls,
                             str_subset(perturbations, "KRAS"),
                             str_subset(perturbations, "RIT1"),
                             str_subset(perturbations, "NRAS"),
                             str_subset(perturbations, "EGFR"),
                             str_subset(perturbations, "p.D837A"),
                             str_subset(perturbations, "MET"),
                             str_subset(perturbations, "ARAF"),
                             str_subset(perturbations, "p.D429A"),
                             str_subset(perturbations, "BRAF"),
                             str_subset(perturbations, "p.W450L"),
                             str_subset(perturbations, "p.H574N"),                             
                             str_subset(perturbations, "p.D594H"),
                             str_subset(perturbations, "MYC"),
                             str_subset(perturbations, "TP53"),
                             str_subset(perturbations, "MDM2"),
                             str_subset(perturbations, "CTNNB1"),
                             str_subset(perturbations, "AURKA"),
                             str_subset(perturbations, "CCND1"),
                             str_subset(perturbations, "PIK3CA"),
                             str_subset(perturbations, "NFE2L2"))

perturbation_plot_order <- unique(perturbation_plot_order, fromLast = TRUE)
df_plot$perturbation <- factor(df_plot$perturbation, levels = perturbation_plot_order)
gene_plot_order <- c("Control", "KRAS", "RIT1", "NRAS", "EGFR", "MET", "ARAF", "BRAF", "MYC", "TP53", "MDM2", "CTNNB1", "AURKA", "CCND1", "PIK3CA", "NFE2L2")
df_plot$perturb_gene <- factor(df_plot$perturb_gene, levels = gene_plot_order)

## repeats from PCA section
control_color <- "#333333"
names(control_color) <- c("Control")
oncogene_colors <- rev(brewer.pal(length(oncogenes), "Set2"))
oncogene_colors <- c(oncogene_colors, rev(brewer.pal(length(oncogenes)-8, "Set1")))
names(oncogene_colors) <- oncogenes
plot_colors <- c(oncogene_colors, control_color)
plot_colors_darken <- darken(plot_colors, 0.2)
names(plot_colors_darken) <- names(plot_colors)

g <- ggplot(data = df_plot, aes(x = perturbation, y = z)) +
    geom_col(aes(fill = perturb_gene)) +
    geom_hline(yintercept = 1, color = "#222222", linetype = "dashed") +
    geom_hline(yintercept = -1, color = "#222222", linetype = "dashed") +    
    geom_text(data = df_plot %>% filter(z > 0),
              aes(label = perturb_variant), position = position_dodge(width = 0.9),
              angle = 90, hjust = -0.1, vjust = 0.5) +
    geom_text(data = df_plot %>% filter(z < 0),
              aes(label = perturb_variant), position = position_dodge(width = 0.9),
              angle = 90, hjust = 1.1, vjust = 0.5) +
    scale_fill_manual(values = plot_colors,
                      name = NULL) +
    scale_y_continuous(limits = c(-9, 5), name = "Normalized Z-score") +
    scale_x_discrete(labels = NULL, name = NULL) +
    guides(fill = FALSE) +
    facet_grid( . ~ perturb_gene, scales = "free", space = "free",
               switch = "x") +
    theme_cowplot(font_size = 20) +
    theme(axis.text.x = element_text(angle = 75, hjust = 1, vjust = 1),
          axis.ticks.x = element_blank(),
          strip.background = element_blank())
ggsave(g, filename = "L1000_AALE_oncogenes_FOSL1.pdf",
       height = 8, width = 20,
       device = "pdf",
       path = resultsDir)

"NFE2L2" %in% gene_ids$pr_gene_symbol

df_plot <- df_AALE %>%
    filter(Geneid == "NFE2L2") %>%
    filter(perturb_gene %in% c(oncogenes, controls))
df_plot$perturbation <- factor(df_plot$perturbation, levels = perturbation_plot_order)
g <- ggplot(data = df_plot, aes(x = perturbation, y = z)) +
    geom_col(aes(fill = perturb_gene)) +
    scale_fill_manual(values = pca_colors,
                      name = NULL) +
    theme_bw(base_size = 20) +
    theme(axis.text.x = element_text(angle = 75, hjust = 1, vjust = 1))
ggsave(g, filename = "L1000_AALE_oncogenes_NFE2L2.pdf",
       height = 10, width = 25,
       device = "pdf",
       path = resultsDir)


## Plot EMT genes?
