s#######################################################################
#### Overview - Do Differential Expression Analysis
#######################################################################

## Goal: Perform differential expression analysis with edgeR
##       and generate overview data structures for downstream analysis

## Data structures:
## - df_expr
## - df_DE

#######################################################################
#### SET UP
#######################################################################

#######################################################################
#### Set relevant input and output directories
#######################################################################

## Inputs
dataDirs <- c("./")
## Outputs
dataDir <- "./"
resultsDir <- "../results/AALE_KRAS-RIT1"
figsDir <- "../results/AALE_KRAS-RIT1/Figures"
qcDir <- "../results/AALE_KRAS-RIT1"

#######################################################################
#### Functions
#######################################################################


source("../scripts/R/load.R")

# To get metadata for samples
getSampleInfo <- function(dataDir) {
  return(read.table(file=paste0(dataDir, "/sampleInfo.txt"), sep="\t", header=TRUE))
}
# To get list of genes (with annotations)
getGenes <- function(dataDir) {
  data_df <- read.table(paste0(dataDir, "/featureCounts/hg19/counts_all.txt"), header=TRUE)
  data_df <- data_df[1:(dim(data_df)[1]-92),]
  genes_df <- data_df[,1:6]
  return(genes_df)
}
# To get transcript quantification in count form
# Transcript quants were done with featureCounts from the subread package
getCounts <- function(dataDir) {
  sampleInfo <- getSampleInfo(dataDir)
  data_df <- read.table(paste0(dataDir, "/featureCounts/hg19/counts_all.txt"), header=TRUE)
  # Exclude ERCCs
  data_df <- data_df[1:(dim(data_df)[1]-92),]
  counts_df <- data_df[7:dim(data_df)[2]]
  names(counts_df) <- as.vector(sampleInfo$SampleId)
  return(counts_df)
}
# To get perturbation groups info
getGrouping <- function(dataDir, exclude) {
  sampleInfo <- getSampleInfo(dataDir)
  sampleInfo <- sampleInfo[-which(sampleInfo$SampleId %in% exclude),]
  N <- dim(sampleInfo)[1]
  grouping <- vector(length=N)
  for (i in 1:N) {
    sample_name <- as.character(sampleInfo$Notes[i])
    l <- as.vector(strsplit2(sample_name, "_"))
    if ("V5" %in% l) {
      l <- l[-which(l=="V5")]
    }
    grouping[i] <- paste(l[1:(length(l)-3)], collapse = "_")
  }
  return(grouping)
}


#######################################################################
#### Load Data
#######################################################################

## Simple dataframe listing genes
df_genes <- getGenes(dataDirs[1])
## Exclude samples for which RNA-seq failed
## e.g. low coverage
exclude <- c()
## Data frame of transcript counts
counts <- getCounts(dataDirs[1])
if (length(exclude) > 0) {
    counts <- counts[,-which(names(counts) %in% exclude)]
}
df_counts <- counts
rm(counts)
# Get group annotation for Screen 1 first (naming/parsing is different)
sampleInfo <- getSampleInfo(dataDirs[1])
if (length(exclude) > 0) {
    sampleInfo <- sampleInfo[-which(sampleInfo$SampleId %in% exclude),]
}
N <- dim(sampleInfo)[1]

grouping <- sampleInfo$Notes
print(length(grouping))

#######################################################################
#### DE ANALYSIS
#######################################################################

#######################################################################
#### Set up for edgeR
#######################################################################

## TODO: If/else statement for whether DGEList_all.RData already exists

### group samples
group <- factor(grouping)
DEList <- DGEList(counts=df_counts, genes=df_genes, group=group)

### scaling
DEList <- calcNormFactors(DEList)
design <- model.matrix(~0+grouping)

y <- estimateDisp(DEList, design)

save(y, file=paste0(dataDir, "/DGEList.RData"))

load(paste0(dataDir, "/DGEList.RData"))

write.table(y$counts, file = paste0(dataDir, "/transcript_counts_raw.tsv"),
            col.names = NA, row.names = y$genes$Geneid , sep = "\t", quote = FALSE)

#######################################################################
#### Convert counts/rpkms to TPMs
#######################################################################

df_fpkm <- rpkm(y)
df_tpm <- data.frame(apply(df_fpkm, 2, fpkm_to_tpm))
df_tpm$Gene <- y$genes$Geneid
df_tpm <- df_tpm %>% gather("Sample", "TPM", -Gene)

########################n###############################################
#### Check overexpression
#######################################################################

## Box/dot plots 

# cpm_df <- data.frame(cpm(y))
# cpm_df$Gene <- y$genes$Geneid
# cpm_df <- cpm_df %>% gather("Sample", "CPM", -Gene)

plot_overexpression <- function(gene) {
  gene_df <- df_tpm %>% filter(Gene==gene)
  threshold <- gene_df %>% filter(grepl("Renilla", gene_df$Sample)) %>% summarise(mean(log2(TPM))+2*sd(log2(TPM)))
  plot_df <- df_tpm %>% filter(Gene==gene) %>% mutate(gene=grepl(gene, Sample))
  plot_df$Group <- as.vector(group)
  plot_df$Group[which(grepl(gene, plot_df$Sample))] <- gene
  plot_df$Group[which(!grepl(gene, plot_df$Sample))] <- "OTHER"
  plot_df$Group[which(grepl("Renilla", plot_df$Sample))] <- "CONTROL"
  plot_df$Group <- factor(plot_df$Group, levels = c("CONTROL", "OTHER", gene))
  plot_df <- arrange(plot_df, Group)
  colors <- c("grey", brewer.pal(8, "Set1")[c(2,1)])
  names(colors) <- c("CONTROL", "OTHER", gene)
  g <- ggplot(plot_df, aes(x=Group, y=log2(TPM), fill=Group)) +
      geom_boxplot(width = 0.5, alpha = 0.7) +
      geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.5) +
      geom_hline(yintercept = threshold[[1]], linetype = "dashed") +
      theme_classic(base_size = 18) +
      theme(plot.title = element_text(size = 18, hjust = 0.5)) +
      scale_fill_manual(values = colors) + guides(fill=FALSE) +
      ylab(bquote("log"[2]*"TPM" ~ of ~ .(gene))) +


  xlab("Overexpressed Gene") +
      ggtitle(gene)
}

for (gene in c("KRAS", "RIT1")) {
  plot_overexpression(gene)
  ggsave(paste0("overexpression/", gene, ".pdf"),
         device = "pdf", width = 7, height = 5,
         path = qcDir)
}

#######################################################################
#### Overexpression rainbow box/dot plots
#######################################################################


plot_overexpression_rainbow <- function(gene) {
  gene_df <- df_tpm %>% filter(Gene==gene)
  threshold <- gene_df %>%
      filter(grepl("Renilla", gene_df$Sample)) %>%
      summarise(mean(log2(TPM))+2*sd(log2(TPM)))
  plot_df <- df_tpm %>%
      filter(Gene==gene) %>%
      mutate(gene=grepl(gene, Sample)) %>%
      mutate(Group = factor(grouping, levels = c("Renilla", "KRAS-WT", "KRAS-G12V", "KRAS-Q61H", "RIT1-WT", "RIT1-M90I"))) %>%
      mutate(GenePerturb = strsplit2(as.character(Group), split = "-")[,1]) %>%
      mutate(GenePerturb = factor(GenePerturb, levels = c("Renilla", "KRAS", "RIT1")))
  plot_df <- arrange(plot_df, Group)
  colors <- c("#999999", brewer.pal(3, "Dark2"))
  names(colors) <- levels(plot_df$GenePerturb)
  group_names <- gsub("-", " ", levels(plot_df$Group))
  g <- ggplot(plot_df, aes(x=Group, y=log2(TPM))) +
      geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 1,
                   alpha = 0.8,
                   shape = 19) +
      geom_hline(yintercept = threshold[[1]],
                 linetype = "dashed", color = "#999999") +
      scale_fill_manual(name = NULL,
                        values = colors) +
      scale_color_manual(name = NULL,
                         values = darken(colors, 0.3)) +
      scale_x_discrete(labels = group_names,
                       name = "Overexpressed Gene") +
      scale_y_continuous(name = bquote("log"[2]*"TPM" ~ of ~ .(gene))) +
      guides(fill = FALSE, color = FALSE) +
      theme_classic(base_size = 14) +
      theme(plot.title = element_text(size = 16, hjust = 0.5),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                       margin = margin(t = 4, b = 6)),
            axis.text.y = element_text(margin = margin(r = 4, l = 6)))
}

for (gene in c("KRAS", "RIT1")) {
  plot_overexpression_rainbow(gene)
  ggsave(paste0("expression-rainbow/", gene, "2.pdf"), device = "pdf", width = 6, height = 5, path = qcDir)
}

plot_overexpression_bars <- function(gene) {
    gene_df <- df_tpm %>% filter(Gene==gene)
    threshold <- gene_df %>%
        filter(grepl("Renilla", gene_df$Sample)) %>%
        summarise(mean(log2(TPM)))
#        summarise(mean(log2(TPM))+2*sd(log2(TPM)))
    plot_df <- df_tpm %>%
        filter(Gene==gene) %>%
        mutate(gene=grepl(gene, Sample)) %>%
        mutate(Group = factor(grouping, levels = c("Renilla", "RIT1-WT", "RIT1-M90I", "KRAS-WT", "KRAS-G12V", "KRAS-Q61H"))) %>%
        mutate(GenePerturb = strsplit2(as.character(Group), split = "-")[,1]) %>%
        mutate(GenePerturb = factor(GenePerturb, levels = c("Renilla", "KRAS", "RIT1")))
    plot_df <- arrange(plot_df, Group)
    df_mean <- plot_df %>%
        group_by(Group) %>%
        summarize(meanLogTPM = mean(log2(TPM)))
    df_sd <- plot_df %>%
        group_by(Group) %>%
        summarize(sdLogTPM = sd(log2(TPM)))
    df_plot <- inner_join(df_mean, df_sd) %>%
        mutate(lower = meanLogTPM - sdLogTPM) %>%
        mutate(upper = meanLogTPM + sdLogTPM)      
    print(df_plot)
    colors <- rep("#999999", length(unique(plot_df$Sample)))
                                        #  names(colors) <- levels(plot_df$GenePerturb)
    group_names <- gsub("-", "\n", levels(plot_df$Group))
    group_names <- gsub("Renilla", "Vector", group_names)
    g <- ggplot(df_plot, aes(x=Group, y=meanLogTPM)) +
        geom_col(width = 0.7, fill = darken("#999999", 0.2)) +
        geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.3) +
        ## geom_boxplot(data = plot_df, aes(x = Group, y = log2(TPM)),
        ##              width = 0.3) +
        ## geom_point(data = plot_df, aes(x = Group, y = log2(TPM))) +      
                                        #      geom_col(position = position_dodge2(padding = 0.1), width = 0.8) +
        geom_hline(yintercept = threshold[[1]],
                   linetype = "dashed", color = "black") +
        scale_color_manual(name = NULL,
                           values = darken(colors, 0.1)) +
        scale_x_discrete(labels = group_names,
                         name = NULL) +
        scale_y_continuous(name = bquote("log"[2]*"TPM" ~ of ~ .(gene)),
                           breaks = c(0, 2, 4, 6, 8),
                           expand = c(0,0)) +
        guides(fill = FALSE, color = FALSE) +
        theme_classic(base_size = 26) +
        theme(plot.title = element_text(size = 16, hjust = 0.5),
              axis.ticks.x = element_blank(),
              axis.text.x = element_text(#angle = 75, hjust = 1, vjust = 1,
                  margin = margin(t = 4, b = 10)),
              axis.text.y = element_text(margin = margin(r = 4, l = 8)))
}

for (gene in c("KRAS", "RIT1", "NRAS", "HRAS")) {
    plot_overexpression_bars(gene)
    ggsave(paste0("expression-rainbow/", gene, "-bars.pdf"), device = "pdf",
         width = 7, height = 5, path = qcDir)
}

#######################################################################
#### T tests for overexpression of KRAS/RIT1
#######################################################################

perturbations <- c("RIT1-WT", "RIT1-M90I", "KRAS-WT", "KRAS-G12V", "KRAS-Q61H")

do_t_tests <- function(gene) {
    df_gene <- df_tpm %>%
        filter(Gene == gene) %>%
        mutate(log2TPM = log2(TPM)) %>%
        mutate(Sample = gsub("[.]", "-", Sample))
    df_control <- df_gene %>%
        filter(grepl("Renilla", Sample))
    df_gene <- df_gene %>%
        filter(!grepl("Renilla", Sample))
    print(df_gene)
    for (perturb in perturbations) {
        print(perturb)
        df_test <- df_gene %>%
            filter(grepl(perturb, Sample))
        print(df_control$log2TPM)
        print(df_test$log2TPM)
        t <- t.test(df_control$log2TPM, df_test$log2TPM,
                    alternative = c("two.sided"))
        print(t$p.value)
    }
}

for (gene in c("KRAS", "RIT1")) {
    do_t_tests(gene)
}

## Using the above to label figure plots:
## *** for < 0.001 pval
## * for < 0.05 pval

#######################################################################
#### Generate (a simple) sample summaries table
#######################################################################

sample_summaries <- y$samples
sample_summaries$Sample <- rownames(y$samples)

#######################################################################
#### Filtering
#######################################################################

## Remove all transcripts that do not have a cpm of at least 1 in at least 2 samples
# Remove all transcripts that do not have an average log cpm of at least 0.1 across all samples
logCPMc <- cpm(y, log=TRUE) # not actually 'c' since there's no batch effect to correct for
df_plot <- as_tibble(logCPMc)
df_plot$Gene <- y$genes$Geneid
df_plot <- df_plot %>%
    gather("Sample", "logCPM", -Gene) %>%
    group_by(Gene) %>%
    summarise(avgLogCPM = mean(logCPM))
g <- ggplot(df_plot, aes(x = avgLogCPM)) +
    geom_density() +
    theme_bw() +
    labs(x = "Mean Expression (log2CPM)",
         y = "Proportion of Genes")
ggsave("expression_mean_density_prefilter.pdf",
       plot = g,
       device = "pdf",
       width = 7,
       height = 7,
       path = resultsDir)

## FOR edgeR
yfiltered <- y[aveLogCPM(y)>0.1,]
# keep <- rowSums(cpm(y)>1) >= 2
# yfiltered <- y[keep,]
dim(yfiltered)

## Replot after filtering

df_plot <- as_tibble(cpm(yfiltered, log=TRUE))
df_plot$Gene <- yfiltered$genes$Geneid
df_plot <- df_plot %>%
    gather("Sample", "logCPM", -Gene) %>%
    group_by(Gene) %>%
    summarise(avgLogCPM = mean(logCPM))

g <- ggplot(df_plot, aes(x = avgLogCPM)) +
    geom_density() +
    theme_bw() +
    labs(x = "Mean Expression (log2CPM)",
         y = "Proportion of Genes")
ggsave("expression_mean_density_postfilter.pdf",
       plot = g,
       device = "pdf",
       width = 7,
       height = 7,
       path = resultsDir)

#######################################################################
#### PCA
#######################################################################

pca <- prcomp(cpm(yfiltered, log = TRUE), center = TRUE, scale = TRUE)
summary(pca)

df_pca <- data.frame(pca$rotation)
df_pca$sample <- rownames(df_pca)
df_pca$group <- factor(group, levels = c("Renilla", "RIT1-WT", "RIT1-M90I", "KRAS-WT", "KRAS-G12V", "KRAS-Q61H"))
df_pca <- df_pca %>%
    mutate(is_WT = grepl("WT", group) | !grepl("-", group))

perturbation_names <- gsub("-", " ", levels(df_pca$group))
perturbation_names <- gsub("Renilla", "Vector", perturbation_names)
pca_colors <- c("#999999", "#D55E00", "#E69F00", "#0072B2", "#009E73", "#56B4E9")
# Plot PC1 and PC2
g <- ggplot(data=df_pca, aes(PC1, PC2)) +
    geom_point(aes(color = group, fill = group, shape = is_WT),
               size = 6, alpha = 0.8) +
    scale_fill_manual(values = pca_colors,
                      labels = perturbation_names,
                      name = NULL) +
    scale_color_manual(values = darken(pca_colors, 0.3),
                       labels = perturbation_names,
                       name = NULL) +
    scale_shape_manual(values = c("TRUE" = 21, "FALSE" = 23),
                       labels = c("WT", "Variant"),
                       name = NULL) +
    scale_x_continuous(name = "PC1 (96.2% Variance Explained)",
                       limits = c(0.229, 0.239),
                       labels = NULL) +
    scale_y_continuous(name = "PC2 (2.16% Variance Explained)",
                       limits = c(-0.3, 0.5),
                       labels = NULL) +
    theme_classic(base_size = 18) +
    guides(
        fill = guide_legend(
            nrow = 2,
            override.aes = list(shape = c(21, 21, 23, 21, 23, 23))),
        shape = FALSE
    ) +
    theme(
        legend.position = "top",
        legend.justification = "center",
        legend.text = element_text(size = 16),
        legend.box.spacing = unit(0, "pt"),
        axis.ticks = element_blank(),
        axis.title.x = element_text(margin = margin(t = 6)),
        axis.title.y = element_text(margin = margin(r = 6))
    )        
 ggsave(g, filename = "/pca.pdf",
       width = 6, height = 5,
       device = "pdf",
       path = resultsDir)

zscaled_cpm <- scale(t(cpm(yfiltered, log = TRUE)), center = TRUE, scale = TRUE)
zscaled_cpm <- t(zscaled_cpm)

correlation <- cor(zscaled_cpm, method = "pearson")
samples <- rownames(correlation)
df_plot <- correlation
df_plot <- as_tibble(correlation) %>%
    mutate(Sample1 = samples) %>%
    gather("Sample2", "correlation", 1:18)

g <- ggplot(df_plot, aes(Sample1, Sample2)) +
    geom_tile(aes(fill = correlation)) +
    scale_fill_gradientn(limits = c(-1, 1), colors = rev(brewer.pal(9, "RdBu"))) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave(filename = "/sample_correlation_matrix.pdf",
       height = 7, width = 8,
       device = "pdf",
       path = resultsDir)



#######################################################################
#### Consolidate into df_expr comprehensive data structure
#######################################################################

## df_expr generated: joins batch corrected logTPM and logCPM (per gene per sample)

## Get cpm values from the DGEList object (yfiltered)
df_cpm <- data.frame(cpm(yfiltered, log=TRUE))
df_cpm$Gene <- yfiltered$genes$Geneid
df_cpm <- df_cpm %>%
    gather("Sample", "logCPM", -Gene)
## Separately compute rpkm -> tpm
df_fpkm <- rpkm(yfiltered)
df_tpm <- data.frame(apply(df_fpkm, 2, fpkm_to_tpm))
logTPM <- apply(df_tpm, 2, log2)

df_save <- tbl_df(logTPM)
colnames(df_save) <- gsub("[.]", "-", colnames(df_save))
## The different tables in the DGEList object yfiltered
## (e.g. cpm/rpkm/tpm values, genes)
## should match each other in row order to make the following save work

write.csv(df_save, "./transcript_quants_log2TPM.csv",
            row.names = yfiltered$genes$Geneid,
#            col.names = NA,
            quote = FALSE)

df_tpm <- data.frame(logTPM)
df_tpm$Gene <- yfiltered$genes$Geneid
df_tpm <- df_tpm %>% gather("Sample", "logTPM", -Gene)

# Join tpm and cpm values
df_expr <- left_join(df_tpm, df_cpm, by = c("Sample", "Gene")) %>%
    mutate(Sample = str_replace(Sample, "[.]", "-"))
# Join with sample metadata
df_expr <- left_join(df_expr, sample_summaries, by = "Sample")

# Remove all transcripts that do not have an average log cpm of at least 0.1 across all samples
## df_expr_all <- df_expr
## save(df_expr_all, file = paste0(dataDir, "df_expr_all.RData"))
## df_expr <- df_expr_all %>%
##     filter(Gene %in% yfiltered$genes$Geneid)
save(df_expr, file = paste0(dataDir, "df_expr.RData"))

load(file = paste0(dataDir, "df_expr.RData"))

write.csv(df_expr, file = paste0(resultsDir, "/transcript_quants_logTPM-tidy.csv"))


df_save <- 
write.csv(df_save, file = paste0(resultsDir, "/all_logTPM.csv"),
              row.names = FALSE)
rm(df_save)

#######################################################################
#### Calculate means and standard deviations of expression of each gene in each perturbation group
#######################################################################

df_cpm <- df_expr %>%
    group_by(.dots=c("group","Gene")) %>%
    summarise(median(logCPM))
medians_matrix <- spread(df_cpm, group, "median(logCPM)")
dim(medians_matrix)

write.csv(medians_matrix,
          file = paste0(resultsDir, "/Gene-expression_group-median-log2CPM.csv"),
          quote = FALSE,
          col.names = TRUE,
          row.names = FALSE)

medians_matrix <- read.csv(file = paste0(resultsDir, "/Gene-expression_group-median-log2CPM.csv"),
                           header = TRUE)

# rownames(medians_matrix) <- medians_matrix$Gene
# medians_matrix <- medians_matrix[,-1]

#######################################################################
#### Correlation Plots
#######################################################################

## Replicate level
cpm_df <- logCPMc
z_df <- t(cpm_df)
# zscore_df <- 
# (apply(z_df, 2, function (x) {
#     return((x - mean(x))/sd(x))
# }))[1:5, 1:10]
z_df <- scale(z_df, center = TRUE, scale = TRUE)
z_df <- t(z_df)

## With collapsed replicates (median expr levels)
z_medians <- t(medians_matrix)
z_medians<- scale(t(medians_matrix),
                   center = TRUE,
                   scale = TRUE)

z_medians <- t(z_medians)

#######################################################################
#### Spearman Correlation Plot, using ggplot
#######################################################################

correlation <- as.data.frame(cor(z_df, method = "pearson"))
correlation$sample_x <- rownames(correlation)
correlation <- gather(correlation, "sample_y", "correlation", -sample_x)
correlation <- left_join(correlation, dplyr::select(sample_summaries, Sample, group), by = c("sample_x" = "Sample"))
colnames(correlation)[4] <- "group_x"
correlation <- left_join(correlation, dplyr::select(sample_summaries, Sample, group), by = c("sample_y" = "Sample"))
colnames(correlation)[5] <- "group_y"


correlation$group_x_labels <- as.character(correlation$group_x)
correlation$group_x_labels[duplicated(correlation$group_x)] <- ""

correlation$group_y_labels <- as.character(correlation$group_y)
correlation$group_y_labels[duplicated(correlation$group_y)] <- ""

ggplot(correlation, aes(x=sample_x, y=sample_y, fill=correlation)) +
    geom_tile() +
    scale_fill_gradientn(limits = c(-1, 1), colors = rev(brewer.pal(9, "RdBu"))) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    scale_x_discrete(labels = arrange(sample_summaries, Sample)$group[correlation_order], breaks = waiver()) +
    scale_y_discrete(labels = arrange(sample_summaries, Sample)$group[correlation_order], breaks = waiver())

# Clustered spearman correlation plot, using ggplot

correlation <- as.data.frame(cor(z_df, method = "pearson"))
correlation_order <- hclust(dist(cor(z_df, method = "pearson")))$order
sample_levels <- rownames(correlation)[correlation_order]
# correlation <- correlation[correlation_order, correlation_order]
correlation$sample_x <- rownames(correlation)
correlation <- gather(correlation, "sample_y", "correlation", -sample_x)
correlation <- left_join(correlation, dplyr::select(sample_summaries, Sample, group), by = c("sample_x" = "Sample"))
colnames(correlation)[4] <- "group_x"
correlation <- left_join(correlation, dplyr::select(sample_summaries, Sample, group), by = c("sample_y" = "Sample"))
colnames(correlation)[5] <- "group_y"
correlation$sample_x <- factor(correlation$sample_x, levels = sample_levels)
correlation$sample_y <- factor(correlation$sample_y, levels = sample_levels)

g <- ggplot(correlation, aes(x=sample_x, y=sample_y, fill=correlation)) +
    geom_tile() +
    scale_fill_gradientn(limits = c(-1, 1), colors = rev(brewer.pal(9, "RdBu"))) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6)) +
    theme(axis.text.y = element_text( size = 6 )) +
    scale_x_discrete(labels = arrange(sample_summaries, Sample)$group[correlation_order], breaks = waiver()) +
    scale_y_discrete(labels = arrange(sample_summaries, Sample)$group[correlation_order], breaks = waiver())

# ggsave("../results/Screens-all/Correlation-matrix_Z-score.pdf", plot = g, device = "pdf", width = 35, height = 35, units = "in", dpi = 300)

## mat is correlation  matrix
mat <- cor(z_df, method="spearman")
# mat <- mat[apply(mat, 1, function(row) all(row!=-Inf)),]

group <- arrange(sample_summaries, Sample)$group
gene <- arrange(sample_summaries, Sample)$gene
batch <- arrange(sample_summaries, Sample)$batch
## Data frame with column annotations.
matrix_col <- data.frame(group = arrange(sample_summaries, Sample)$group,
                         batch = arrange(sample_summaries, Sample)$batch)
rownames(matrix_col) <- colnames(mat)

## List with colors for each annotation.
matrix_colors <- list(group = colorRampPalette(brewer.pal(9, "Set1"))(length(unique(group))),
                      gene = colorRampPalette(brewer.pal(12, "Set3"))(length(unique(gene))),
                      batch = brewer.pal(5, "Accent"))
# matrix_colors <- list(batch = brewer.pal(5, "Accent"))
names(matrix_colors$group) <- unique(group)
names(matrix_colors$gene) <- unique(gene)
names(matrix_colors$batch) <- unique(batch)

matrix_cluster_cols <- hclust(dist(mat))
matrix_cluster_cols <- sort_hclust(matrix_cluster_cols)

matrix_cluster_rows <- hclust(dist(mat))
matrix_cluster_rows <- sort_hclust(matrix_cluster_rows)
matrix_row <- data.frame(group = arrange(sample_summaries, Sample)$group,
                         gene = arrange(sample_summaries, Sample)$gene,
                         batch = arrange(sample_summaries, Sample)$batch)
rownames(matrix_row) <- rownames(mat)

#######################################################################
#### Spearman Correlation Plot, without clustering
#######################################################################

## First without clustering

# group_x_labels <- as.character(arrange(sample_summaries, Sample)$group)
# group_x_labels[duplicated(group_x_labels)] = ""
# rownames(matrix_row) <- group_x_labels

pheatmap(
  mat               = mat,
  color             = rev(brewer.pal(11, "RdBu")),
  breaks            = seq(-1, 1, length.out = 12),
  border_color      = NA,
  cluster_cols      = FALSE,
  cluster_rows      = FALSE,
  show_colnames     = TRUE,
  show_rownames     = TRUE,
  annotation_col    = matrix_col,
  annotation_colors = matrix_colors,
  annotation_row    = matrix_row,
  annotation_legend = TRUE,
  annotation_names_col = TRUE,
  annotation_names_row = TRUE,
  treeheight_col    = 0,
  treeheight_row    = 0,
  drop_levels       = TRUE,
  display_numbers   = FALSE,
  legend            = TRUE,
  fontsize          = 10,
  main              = "Z-score based correlation"
)

#######################################################################
#### Hierarchically clustered Spearman Correlation Plot, of replicates
#######################################################################

## Then hierarchically clustered (samples on both rows and columns)

pheatmap(
  mat               = mat,
  color             = rev(brewer.pal(11, "RdBu")),
  breaks            = seq(-1, 1, length.out = 12),
  border_color      = NA,
  cluster_cols      = matrix_cluster_cols,
  cluster_rows      = matrix_cluster_rows,
  show_colnames     = TRUE,
  show_rownames     = TRUE,
  annotation_col    = matrix_col,
  annotation_colors = matrix_colors,
  # annotation_row    = matrix_row,
  annotation_names_col = TRUE,
  # annotation_names_row = FALSE,
  # treeheight_col    = 0,
  # treeheight_row    = 0,
  drop_levels       = TRUE,
  display_numbers   = FALSE,
  # legend            = TRUE,
  fontsize          = 12,
  fontsize_row      = 6,
  fontsize_col      = 6,
  width             = 40,
  height            = 35,
  filename          = "./results/Screens-all/Correlation-matrix_Z-score.pdf",
  main              = "Z-score based spearman correlation"
)


#######################################################################
#### Hierarchically clustered Spearman Correlation Plot, of perturbation groups
#######################################################################

mat <- cor(z_medians, method="spearman")
# mat <- mat[apply(mat, 1, function(row) all(row!=-Inf)),]

gene <- arrange(sample_summaries, Sample)$gene

## Data frame with column annotations.
matrix_col <- data.frame(group = arrange(sample_summaries, Sample)$group,
                         batch = arrange(sample_summaries, Sample)$batch)
rownames(matrix_col) <- colnames(mat)

## List with colors for each annotation.
matrix_colors <- list(group = colorRampPalette(brewer.pal(9, "Set1"))(length(unique(group))),
                      gene = colorRampPalette(brewer.pal(12, "Set3"))(length(unique(gene))),
                      batch = brewer.pal(5, "Accent"))
# matrix_colors <- list(batch = brewer.pal(5, "Accent"))
names(matrix_colors$group) <- unique(group)
names(matrix_colors$gene) <- unique(gene)
names(matrix_colors$batch) <- unique(batch)

matrix_cluster_cols <- hclust(dist(mat))
matrix_cluster_cols <- sort_hclust(matrix_cluster_cols)

matrix_cluster_rows <- hclust(dist(mat))
matrix_cluster_rows <- sort_hclust(matrix_cluster_rows)

matrix_row <- data.frame(group = arrange(sample_summaries, Sample)$group,
                         gene = arrange(sample_summaries, Sample)$gene,
                         batch = arrange(sample_summaries, Sample)$batch)
rownames(matrix_row) <- rownames(mat)


pheatmap(
  mat               = mat,
  color             = rev(brewer.pal(11, "RdBu")),
  breaks            = seq(-1, 1, length.out = 12),
  border_color      = NA,
  cluster_cols      = matrix_cluster_cols,
  cluster_rows      = matrix_cluster_rows,
  show_colnames     = TRUE,
  show_rownames     = TRUE,
  # annotation_col    = matrix_col,
  # annotation_colors = matrix_colors,
  # annotation_row    = matrix_row,
  # annotation_names_col = TRUE,
  # annotation_names_row = FALSE,
  # treeheight_col    = 0,
  # treeheight_row    = 0,
  drop_levels       = TRUE,
  display_numbers   = FALSE,
  # legend            = TRUE,
  fontsize          = 16,
  fontsize_row      = 12,
  fontsize_col      = 12,
  width             = 20,
  height            = 20,
  filename          = "./results/Screens-all/Correlation-matrix_Z-score_replicates-collapsed.pdf",
  main              = "Z-score based spearman correlation"
)


#######################################################################
#### Differential expression analysis using edgeR
#######################################################################

# Uses linear model defined above
# And edgeR test glmQLFTest

resultsDir <- "../results/AALE_KRAS-RIT1/diff-expression"
vector_controls <- c("Renilla")
group <- factor(grouping)
conditions <- sort(unique(group))
design <- model.matrix(~0+grouping)
fit <- glmQLFit(yfiltered, design)
# fit <- glmFit(y, design)

#######################################################################
#### Perform edgeR analyses
#######################################################################

## Differential expression analysis comparing each perturbation group to the vector controls which are HCRED, EGFP, and LUCIFERASE.
## Outputs:
## - Volcano plots for each perturbation vs each vector control, saved as PDF
## - DE genes list ordered by p value for each perturbation vs each vector controls, saved as CSV

## number of groups considered in linear model e.g. total of batches and conditions
n_groups <- length(colnames(design))
## get the column indices of control conditions in the design matrix
control_conditions <- sapply(vector_controls, function(x) which(grepl(x, colnames(design))))
rank <- vector()
genes <- vector()

df_DE <- data.frame(Perturbation = factor(),
                    Control = factor(),
                    Geneid = factor(),
                    Chr = factor(),
                    Start = factor(),
                    End = factor(),
                    Strand = factor(),
                    Length = integer(),
                    logFC = double(),
                    logCPM = double(),
                    F = double(),
                    PValue = double(),
                    FDR = double())
df_DE <- as_tibble(df_DE)
lfcCutoff <- 1
tests <- list(c("KRAS-WT", "Renilla"),
              c("KRAS-G12V", "Renilla"),
              c("KRAS-Q61H", "Renilla"),              
              c("RIT1-WT", "Renilla"),
              c("RIT1-M90I", "Renilla"),
              c("KRAS-G12V", "KRAS-WT"),
              c("KRAS-Q61H", "KRAS-WT"),
              c("RIT1-M90I", "RIT1-WT"))
for (test in tests) {
    perturbation = test[1]
    control = test[2]
    print(c(perturbation, control))
    control_idx <- which(grepl(control, colnames(design)))
    perturb_idx <- which(grepl(perturbation, colnames(design)))
    contrast <- rep(0, length(colnames(design)))
    contrast[control_idx] <- -1
    contrast[perturb_idx] <- 1
    print(contrast)
    qlf <- glmQLFTest(fit, contrast =contrast)
    tt <- topTags(qlf, n=dim(fit)[1], sort.by="PValue")
#    print(head(tt$table))
    sigTestResults <- decideTests(qlf, p.value = 0.01, adjust.method = "BH")
    print(head(sigTestResults))
    print(summary(sigTestResults))
    gene <- strsplit2(perturbation, "-")[1]
    gene_regex <- paste0("^", gene, "$")
    df_plot <- as_tibble(tt$table) %>%
        mutate(overexpressed = grepl(gene_regex, Geneid)) %>%
        arrange(overexpressed)
#    names(colors = c("FALSE", "TRUE")
    g <- ggplot(df_plot, aes(x = logFC, y = -log10(FDR))) +
        geom_point(data = df_plot %>% filter(abs(logFC) < lfcCutoff | FDR > 0.05),
                   alpha = 0.5, shape = 16, color = "grey") +
        geom_point(data = df_plot %>% filter(logFC <= -lfcCutoff & FDR <= 0.05),
                   alpha = 0.5, shape = 16, color = "blue") +
        geom_point(data = df_plot %>% filter(logFC >= lfcCutoff & FDR <= 0.05),
                   alpha = 0.5, shape = 16, color = "red") +
        geom_point(data = df_plot %>% filter(overexpressed),
                   alpha = 0.5, shape = 19, color = "purple") +
        geom_vline(xintercept = lfcCutoff, linetype = "dashed", color = "grey") +
        geom_vline(xintercept = -lfcCutoff, linetype = "dashed", color = "grey") +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey") +   
        geom_text_repel(data = filter(df_plot, overexpressed),
                        aes(label = Geneid),
                        nudge_x = 1, nudge_y = 1) +
        xlim(-10, 10) +
        ylim(0, 30) +
        ggtitle(paste(perturbation, "vs", control)) +
        theme_classic()
    ggsave(paste0(perturbation, "_", control, "_volcano.pdf"),
           plot = g,
           device = "pdf", width = 7, height = 6,
           path = resultsDir)
    DE_results <- as_tibble(tt$table)
    up <- DE_results %>%
        filter(logFC > 0)
    down <- DE_results %>%
        filter(logFC < 0)
    write.csv(up, paste0(resultsDir, "/", perturbation, "_", control, "_up.csv"))
    write.csv(down, paste0(resultsDir, "/", perturbation, "_", control, "_down.csv"))
    DE_results$Perturbation <- perturbation
    DE_results$Control <- control
    df_DE <- bind_rows(df_DE, DE_results)
#    print(head(DE_results))
}

save(df_DE, file = paste0(dataDir, "df_DE.RData"))

load(file = paste0(dataDir, "df_DE.RData"))


#### Re-save csv's DE of each perturbation vs Renilla
resultsDir <- "../results/AALE_KRAS-RIT1/diff-expression"
     for (test in tests[1:5]) {
         perturbation = test[1]
         control = test[2]
         df_save <- df_DE %>%
             filter(Perturbation == perturbation) %>%
             filter(Control == control)
         #print(head(df_save))
         write.csv(df_save, paste0(resultsDir, "/", perturbation, "_", control, ".csv"))
         }

#######################################################################
#### Compare KRAS muts to RIT1s
#######################################################################

## KRAS Q61H+G12V vs RIT1 WT+M90I
## Control = KRAS, Perturb = RIT1

KRAS_idx <- c(which(grepl("KRAS-G12V", colnames(design))),
              which(grepl("KRAS-Q61H", colnames(design))))
RIT1_idx <- c(which(grepl("RIT1-WT", colnames(design))),
              which(grepl("RIT1-M90I", colnames(design))))
contrast <- rep(0, length(colnames(design)))
contrast[KRAS_idx] <- -1
contrast[RIT1_idx] <- 1
print(contrast)
qlf <- glmQLFTest(fit, contrast =contrast)
tt <- topTags(qlf, n=dim(fit)[1], sort.by="PValue")
sigTestResults <- decideTests(qlf, p.value = 0.01, adjust.method = "BH")
print(head(sigTestResults))
print(summary(sigTestResults))

gene <- strsplit2(perturbation, "-")[1]
gene_regex <- paste0("^", gene, "$")
df_plot <- as_tibble(tt$table) %>%
    mutate(overexpressed = grepl(gene_regex, Geneid)) %>%
    arrange(overexpressed)

df_DE_KRASvRIT1 <- as_tibble(tt$table)
df_DE_KRASvRIT1$Perturbation <- "RIT1-all"
df_DE_KRASvRIT1$Control <- "KRAS-muts"

save(df_DE_KRASvRIT1, file = paste0(dataDir, "df_DE_KRASvRIT1.RData"))

write.csv(df_DE_KRASvRIT1, paste0(resultsDir, "/RIT1-all_vs_KRAS-muts.csv"))

load(file = paste0(dataDir, "df_DE_KRASvRIT1.RData"))

up <- as_tibble(tt$table) %>% filter(logFC > 0)
down <- as_tibble(tt$table) %>% filter(logFC < 0)
write.csv(down,
          paste0(resultsDir, "/RIT1-all_KRAS-muts_down.csv"))

up_sig <- up %>% filter(abs(logFC) > 1 & FDR < 0.05)
write(as.character(up$Geneid),
          paste0(resultsDir, "/RIT1-all_KRAS-muts_up_sig.txt"))
down_sig <- down %>% filter(abs(logFC) > 1 & FDR < 0.05)
write(as.character(down$Geneid),
          paste0(resultsDir, "/RIT1-all_KRAS-muts_down_sig.txt"))

    
sig <- as_tibble(tt$table) %>% filter(abs(logFC) > 1 & FDR < 0.05)
write(as.character(sig$Geneid),
      paste0(resultsDir, "/RIT1-all_KRAS-muts_sig.txt"))
sig <- sig %>% dplyr::slice(2:201)
write(as.character(sig$Geneid),
          paste0(resultsDir, "/RIT1-all_KRAS-muts_top200sig.txt"))

df_tpm <- df_expr %>%
    dplyr::select(Gene, Sample, logTPM) %>%
    mutate(SampleTPM = paste0("logTPM_", Sample)) %>%
    dplyr::select(Gene, SampleTPM, logTPM) %>%
    spread(SampleTPM, logTPM)
head(df_tpm)

df_logFC <- df_DE %>%
    filter(Control=="Renilla") %>%
    dplyr::select(Perturbation, Control, Geneid, logFC) %>%
    mutate(Test = paste0("logFC_", Perturbation, "_vs_", Control)) %>%
    dplyr::select(Test, Geneid, logFC) %>%
    spread(Test, logFC)

df_save <- left_join( df_DE_KRASvRIT1 , df_logFC ) %>%
    left_join(df_tpm, by = c("Geneid" = "Gene")) %>%
    dplyr::select(-Chr, -Start, -End, -Strand, -Length)

write.csv(df_save, paste0(resultsDir, "/RIT1-all_vs_KRAS-muts_supplemented.csv"))


#######################################################################
#### Write logFCs and FDRs to table for Supplemental Figure
#######################################################################

perturbations <- unique(df_DE$Perturbation)

df <- df_DE %>%
    filter(Control == "Renilla") %>%
    arrange(Geneid)
df_save <- data.frame(Geneid = unique(df$Geneid))
for (perturbation in perturbations) {
    foo <- df %>% filter(Perturbation == perturbation) %>%
        dplyr::select(Geneid, logFC, FDR)
    names(foo) <- c("Geneid", paste0(perturbation, "_Log2FC"), paste0(perturbation, "_FDR"))
    df_save <- left_join(df_save, foo, by = "Geneid")
    print(head(foo))
    print(dim(df_save))
}

## Supp. Table 4
write.csv(df_save, paste0(resultsDir, "/diff-expression/DE-all_vs_Renilla-logFC-FDR.csv"),
          row.names = FALSE)    

df_DE_KRASvM90I <- data.frame(Perturbation = factor(),
                    Control = factor(),
                    Geneid = factor(),
                    Chr = factor(),
                    Start = factor(),
                    End = factor(),
                    Strand = factor(),
                    Length = integer(),
                    logFC = double(),
                    logCPM = double(),
                    F = double(),
                    PValue = double(),
                    FDR = double())
df_DE <- as_tibble(df_DE)
KRAS_muts <- c("KRAS-Q61H", "KRAS-G12V")
M90I_idx <- which(grepl("RIT1-M90I", colnames(design)))
for (KRAS_mut in KRAS_muts) {
    KRAS_idx <- which(grepl(KRAS_mut, colnames(design)))
    contrast <- rep(0, length(colnames(design)))
    contrast[KRAS_idx] <- -1
    contrast[M90I_idx] <- 1
    print(KRAS_mut)
    print(contrast)
    qlf <- glmQLFTest(fit, contrast =contrast)
    tt <- topTags(qlf, n=dim(fit)[1], sort.by="PValue")
    sigTestResults <- decideTests(qlf, p.value = 0.01, adjust.method = "BH")
    print(head(sigTestResults))
    print(summary(sigTestResults))
    gene <- strsplit2(perturbation, "-")[1]
    gene_regex <- paste0("^", gene, "$")
    df_plot <- as_tibble(tt$table) %>%
        mutate(overexpressed = grepl(gene_regex, Geneid)) %>%
        arrange(overexpressed)
    df <- as_tibble(tt$table)
    df$Perturbation <- "RIT1-M90I"
    df$Control <- KRAS_mut
    df_DE_KRASvM90I <- bind_rows(df_DE_KRASvM90I, df)
}

df_G12V <- df_DE_KRASvM90I %>% filter(Control == "KRAS-G12V") %>%
    filter(FDR < 0.05) %>%
    filter(abs(logFC) > 1)
df_Q61H <- df_DE_KRASvM90I %>% filter(Control == "KRAS-Q61H") %>%
    filter(FDR < 0.05) %>%
    filter(abs(logFC) > 1)

df_intersectDE <- inner_join(df_G12V, df_Q61H, by = "Geneid") %>%
    dplyr::select("Geneid", "logFC.x", "FDR.x", "logFC.y", "FDR.y")
dim(df_intersectDE)
