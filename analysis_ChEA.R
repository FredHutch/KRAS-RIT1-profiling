#######################################################################
#### Plot ChEA results
#######################################################################

## Goal: Analyze ChEA results and generate analysis plots

## Data structures:
## - df_DE

#######################################################################
#### SET UP
#######################################################################

source("../scripts/R/load.R")

#######################################################################
#### Set relevant input and output directories
#######################################################################

## Inputs
dataDirs <- c("./")
## Outputs
dataDir <- "../results/AALE_KRAS-RIT1/ChEA"
resultsDir <- "../results/AALE_KRAS-RIT1/ChEA"
figsDir <- "../results/AALE_KRAS-RIT1/Figures"

load("df_DE.RData")

perturbations <- unique(df_DE$Perturbation)

#######################################################################
#### Load ChEA3 output files
#######################################################################

libraries <- c("ENCODE_ChIP-seq",
               "ReMap_ChIP-seq",
               "Literature_ChIP-seq",
               "ARCHS4_Coexpression",
               "GTEx_Coexpression",
               "Enrichr_Queries")

for (library in libraries) {
    foo <- fread(paste0(dataDir, "/RIT1-M90I_Renilla_", library, ".tsv"))
    bar <- fread(paste0(dataDir, "/KRAS-G12V_Renilla_", library, ".tsv"))
}

#######################################################################
#### Plot top enriched TFs (Enrichr query)
#######################################################################

df_enrichr <- data.frame(Query_Name = character(),
                         Rank = integer(),
                         Scaled_Rank = double(),
                         Set_name = character(),
                         TF = character(),
                         Intersect = integer(),
                         FET_pvalue = double(),
                         FDR = double(),
                         Odds_Ratio = double(),
                         Library = character())

## Load Enrichr Queries results table
for (perturb in perturbations) {
    foo <- fread(paste0(dataDir, "/", perturb, "_Renilla_Enrichr_Queries.tsv"))
    foo <- as_tibble(foo) %>%
        mutate(Perturbation = perturb)
    names(foo) <- gsub(" ", "_", names(foo))
    names(foo) <- gsub("p-value", "pvalue", names(foo))
    print(head(foo))
    dim(filter(foo, FDR < 0.05))
    df_enrichr <- bind_rows(df_enrichr, foo)
}
df_enrichr <- as_tibble(df_enrichr)


EMT_gene_list_short <- c("HIC1", "TWIST1", "HOXD9", "FOXQ1", "FOSL1")
df_EMT_genes <- data.frame(TF = EMT_gene_list_short, EMT_status = "DB-Confirmed")
EMT_gene_list <- c("PRRX2", "FOXS1", "FOXC2", "BNC1", "FOSB", "HOXC6", "HOXB9", "FOXF1", "SIM2", "FOXL1", "EHF", "EGR2", "ATF3", "ZNF750", "FOXF2", "HOXC8", "JDP2", "HOXA10", "OVOL1", "ELF5", "SOX7", "MSX2", "ASCL2", "RELB", "HEY1")
+df <- data.frame(TF = EMT_gene_list, EMT_status = "Literature")
df_EMT_genes <- bind_rows(df_EMT_genes, df)
df <- data.frame(TF = c("SNAI1", "SNAI2"), EMT_status = "SNAIL")
df_EMT_genes <- bind_rows(df_EMT_genes, df)

## Take top 25 enriched TFs by p-value
for (perturb in perturbations) {
    df_plot <- df_enrichr %>%
        filter(Perturbation == perturb) %>%
        arrange(FET_pvalue) %>%
        dplyr::slice(1:25) %>%
        left_join(df_EMT_genes)
    df_plot$EMT_status[is.na(df_plot$EMT_status)] <- "No"
    xmax <- max(-log10(df_plot$FET_pvalue))
    colors <- c("#de2d26", "#fc9272", "#fee0d2", "#bdbdbd")
    names(colors) <- c("SNAIL", "DB-Confirmed", "Literature", "No")
    ## Plot horizontal bar plot of top 25 TF p-values
    g <- ggplot(df_plot) +
        geom_col(aes(x = reorder(TF, -FET_pvalue),
                     y = -log10(FET_pvalue),
                     fill = EMT_status),
                 width = .75) +
        coord_flip() +
        scale_fill_manual(values = colors,
                          guide = FALSE) +
        scale_x_discrete(name = "Transcription Factor") +    
        scale_y_continuous(name = "-log10(FET p-value)",
                           limits = c(0, xmax+2),
                           expand = c(0,0)) +                           
        theme_classic(base_size = 24) +
        theme(axis.ticks.y = element_blank(),
              axis.text.x = element_text(margin = margin(t = 6, b = 8),
                                         size = 24),
              axis.text.y = element_text(margin = margin(r = 4, l = 6),
                                         size = 20))
    ggsave(plot = g, file = paste0("/", perturb, "_Enrichr-top25-bars.pdf"),
                                   device = "pdf",
                                   height = 8,
                                   width = 10,
                                   path = resultsDir)
}

#######################################################################
#### Compare KRAS and RIT1 Enrichr query results
#######################################################################

foo <- fread(paste0(dataDir, "/RIT1-M90I_Renilla_Enrichr_Queries.tsv"))
bar <- fread(paste0(dataDir, "/KRAS-G12V_Renilla_Enrichr_Queries.tsv"))

df_enrichr <- inner_join(foo, bar,
                         by = c("Query Name", "TF", "Library"),
                         suffix = c("_RIT1-M90I", "_KRAS-G12V"))

df_enrichr %>% filter(TF == "EHF")

df_plot <- df_enrichr
g <- ggplot(df_plot) +
    geom_point(aes(x = -log10(`FET p-value_KRAS-G12V`), y = -log10(`FET p-value_RIT1-M90I`)),
               alpha = 0.5, size = 2, shape = 16) +
    geom_abline(color = "grey") +
    scale_x_continuous(limits = c(0, 47)) +
    scale_y_continuous(limits = c(0, 35)) +    
    theme_classic(base_size = 18)
ggsave(g, file = "/Enrichr_pval_scatterplot.pdf",
       device = "pdf",
       height = 6,
       width = 6,
       path = resultsDir)

foo1 <- mutate(foo, Perturbation = "RIT1-M90I")
bar1 <- mutate(bar, Perturbation = "KRAS-G12V")
df_plot <- bind_rows(foo1, bar1) %>%
    mutate(Perturbation = factor(Perturbation,
                                 levels = c("KRAS-G12V", "RIT1-M90I")))
g <- ggplot(df_plot, aes(x = Perturbation, y = -log10(FET_pvalue), group = TF)) +
    geom_point() +
    stat_summary(geom="line") +
    theme_classic(base_size = 18)
ggsave(g, file = "/Enrichr_difference.pdf",
       device = "pdf",
       height = 7,
       width = 6,
       path = resultsDir)

#######################################################################
#### Compare KRAS and RIT1 overall TF enrichment
#######################################################################

## Score := Mean Integrated Rank
## Library := Rank of TF in each library analysis that could be performed

foo <- as_tibble(fread(paste0(dataDir, "/RIT1-M90I_Renilla_Integrated_meanRank.tsv")))
bar <- as_tibble(fread(paste0(dataDir, "/KRAS-G12V_Renilla_Integrated_meanRank.tsv")))
df_chea <- inner_join(foo, bar,
                      by = c("Query Name", "TF"),
                      suffix = c("_RIT1_M90I", "_KRAS_G12V"))

df_chea %>% filter(TF == "EHF")

df_plot <- df_chea
g <- ggplot(df_plot) +
    geom_point(aes(x = Score_KRAS_G12V, y = Score_RIT1_M90I),
               alpha = 0.5, size = 2, shape = 16) +
    geom_abline(color = "grey") +
    theme_classic(base_size = 18)
ggsave(g, file = "/MeanRank_scatterplot.pdf",
       device = "pdf",
       height = 6,
       width = 6,
       path = resultsDir)

df_chea_top <- df_chea %>%
    filter(Rank_RIT1_M90I < 50 | Rank_KRAS_G12V < 50) %>%
    mutate(score_diff = abs(Score_RIT1_M90I - Score_KRAS_G12V))
#    filter(score_diff > 200)

foo1 <- mutate(foo, Perturbation = "RIT1-M90I")
bar1 <- mutate(bar, Perturbation = "KRAS-G12V")
df_plot <- bind_rows(foo1, bar1) %>%
    mutate(Perturbation = factor(Perturbation,
                                 levels = c("KRAS-G12V", "RIT1-M90I"))) %>%
    filter(TF %in% df_chea_top$TF)
g <- ggplot(df_plot, aes(x = Perturbation, y = Score, group = TF)) +
    geom_point() +
    stat_summary(geom="line", alpha = 0.5) +
    theme_classic(base_size = 18)
ggsave(g, file = "/MeanRank_difference_top.pdf",
       device = "pdf",
       height = 6,
       width = 6,
       path = resultsDir)
