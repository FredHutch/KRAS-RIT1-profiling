#######################################################################
#### Set up
#######################################################################

source("../scripts/R/load.R")

## Inputs
dataDirs <- c("./")
## Outputs
dataDir <- "./"
resultsDir <- "../results/AALE_KRAS-RIT1"
figsDir <- "../results/AALE_KRAS-RIT1/Figures"
qcDir <- "./QC"

#######################################################################
#### Load data
#######################################################################

load(file = paste0(dataDir, "DGEList.RData"))
load(file = paste0(dataDir, "df_expr.RData"))
load(file = paste0(dataDir, "df_DE.RData"))

#######################################################################
#### Generate (a simple) sample summaries table
#######################################################################

sample_summaries <- y$samples
sample_summaries$Sample <- rownames(y$samples)

#######################################################################
#### Summarize read counts
#######################################################################

readcounts <- read.table(paste0(qcDir, "/readcounts/all.txt"))
names(readcounts) <- c("Sample", "Mapped_Reads")

sample_summaries <- left_join(sample_summaries, readcounts, by = "Sample")

## Histogram of library size/read counts
g <- ggplot(sample_summaries) +
    geom_histogram(aes(x = n_reads), binwidth = 10000000) +
    theme_classic()
ggsave(g, filename = paste0(qcDir, "/readcounts_histogram.pdf"),
       width = 7,
       height = 6,
       device = "pdf")

sum(readcounts$n_reads)
median(readcounts$n_reads)


#######################################################################
#### Plot gene body coverages
#######################################################################

#gene_body_coverages <- data.frame(Sample = sample_summaries$Sample)
gene_body_coverages <- data.frame()
for (sample in sample_summaries$Sample) {
    gene_body_coverage <- tbl_df(fread(paste0(qcDir, "/gene-body-coverage/", sample, ".geneBodyCoverage.txt"), header=TRUE))
#    print(gene_body_coverage)
#    print(head(t(gene_body_coverage)))
#    as.numeric(gene_body_coverage)
#    colnames(gene_body_coverage) <- gene_body_coverage[1,]
#    gene_body_coverage <- t(as.numeric(gene_body_coverage))
    gene_body_coverages <- bind_rows(gene_body_coverages, gene_body_coverage)
}

df_plot <- gene_body_coverages
df_plot$Sample <- gene_body_coverages[,1]
df_plot <- df_plot[,-1] %>%
    gather("Percentile", "Reads", -Sample) %>%
    mutate(Percentile = as.numeric(Percentile)) %>%
    mutate(Sample = gsub("KRAS_", "KRAS-", Sample)) %>%
    mutate(Sample = gsub("RIT1_", "RIT1-", Sample)) %>%               
    left_join(sample_summaries, by = "Sample")
df_plot$group <- factor(group, levels = c("Renilla", "RIT1-WT", "RIT1-M90I", "KRAS-WT", "KRAS-G12V", "KRAS-Q61H"))
perturbation_names <- gsub("-", " ", levels(df_plot$group))
perturbation_names <- gsub("Renilla", "Vector", perturbation_names)
plot_colors <- c("#999999", "#D55E00", "#E69F00", "#0072B2", "#009E73", "#56B4E9")
## Make line  plot
g <- ggplot(df_plot) +
    geom_line(aes(x = Percentile, y = Reads,
                  group = Sample, color = group),
              size = 1) +
    scale_color_manual(values = brewer.pal(length(unique(sample_summaries$group)), "Set2"),
                       name = NULL,
                       labels = perturbation_names) +
    scale_y_continuous(name = "Read Coverage",
                       limits = c(0, 8500000)) +
    scale_x_continuous(name = "Gene body percentile (5' -> 3')") +
    guides(
        color = guide_legend(
            nrow = 2)
    ) +
    theme_classic(base_size = 16) +
    theme(
        legend.position = "top",
        legend.justification = "center",
        legend.box.spacing = unit(0, "pt")
    )
ggsave(g, filename = "/gene_body_coverages-all.pdf",
       width = 7,
       height = 6,
       device = "pdf",
       path = qcDir)


#######################################################################
#### Summarize rRNA counts
#######################################################################

rRNAcounts <- tbl_df(fread(paste0(qcDir, "/rRNAcounts/rRNAcounts-all.txt")))
names(rRNAcounts) <- c("count_type", "Sample", "n")
rRNAcounts <- spread(rRNAcounts, count_type, n)
names(rRNAcounts) <- c("Sample", "rRNA_Reads", "Total_Reads")

sample_summaries <- left_join(sample_summaries, rRNAcounts, by = "Sample")

sample_qc <- sample_summaries %>%
    dplyr::select(group, Sample, Total_Reads, Mapped_Reads, rRNA_Reads)

write.csv(sample_qc, paste0(qcDir, "/QC-reads-summary.csv"),
          quote = FALSE, row.names = FALSE)

