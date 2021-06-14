#######################################################################
#### Analysis of differentially expressed genes
#######################################################################

## Includes:
## - Plotting edgeR results
## - Gene set analysis of DE genes
## - Look at individual genes

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
qcDir <- "../results/AALE_KRAS-RIT1"

#######################################################################
#### Load data
#######################################################################

load(file = paste0(dataDir, "df_expr.RData"))
load(file = paste0(dataDir, "df_DE.RData"))

controls <- c("Renilla", "KRAS-WT")
for (control in controls) {
    df_G12V <- df_DE %>%
        filter(Perturbation == "KRAS-G12V") %>%
        filter(Control == control) %>%
        filter(FDR < 0.05) %>%
        filter(abs(logFC) > 1)
    df_Q61H <- df_DE %>%
        filter(Perturbation == "KRAS-Q61H") %>%
        filter(Control == control) %>%
        filter(FDR < 0.05) %>%
        filter(abs(logFC) > 1)
    foo <- inner_join(df_G12V, df_Q61H, by = c("Geneid", "Control", "Chr", "Start", "End", "Strand", "Length"))
    print(dim(df_G12V))
    print(dim(df_Q61H))
    print(dim(foo))
    listInput <- list(G12V = df_G12V$Geneid, Q61H = df_Q61H$Geneid)
    ## pdf(paste0(resultsDir, "/Venn-vs", control, "-DE.pdf"),
    ##     width = 5,
    ##     height = 5)
    ## venn(listInput)
    ## dev.off()
}

#######################################################################
#### Gene set analysis with goseq
#######################################################################

hallmarks <- read.table("../mSigDB/h.all.v6.2.symbols.gmt", sep = "\t", header = FALSE, fill = TRUE)
hallmarks <- as_tibble(t(hallmarks[, -2]))
colnames(hallmarks) <- hallmarks[1,]
hallmarks <- hallmarks[-1,]
hallmarks <- gather(hallmarks, "term", "gene")

df_GO <- data.frame(Perturbation = character(),
                    Control = character(),
                    Test = character(),
                    category = character(),
#                    direction = character(),
                    over_represented_pvalue = double(),
                    under_represented_pvalue = double(),
                    numDEInCat = integer(),
                    numInCat = integer(),
                    term = character(),
                    ontology = character())
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
        GO.results <- do_goseq_DE(df_DE, perturbation, control, test = "hallmarks")
        GO.results$Perturbation <- perturbation
        GO.results$Control <- control
        GO.results$Test <- "hallmarks"
#        GO.results$direction <- direction
        df_GO <- bind_rows(df_GO, GO.results)
}
df_GO <- df_GO %>%
    mutate(FDR_enriched = p.adjust(over_represented_pvalue, method="BH")) %>%
    mutate(FDR_unenriched = p.adjust(under_represented_pvalue, method="BH"))

df_goseq <- do_goseq_DE(df_DE, "KRAS-G12V", "Renilla", test = "hallmarks")

save(df_GO, file = "./df_GO.RData")

load("./df_GO.RData")

df_plot <- df_GO %>%
    mutate(Perturbation = gsub("-", " ", Perturbation)) %>%
    mutate(category = gsub("HALLMARK_", "", category)) %>%
    mutate(category = gsub("_", " ", category)) %>%
#    filter(direction == "up") %>%
    filter((Control == "Renilla")) %>%
    mutate(significant = FDR_enriched < 0.05)
df_plot <- df_plot %>% mutate(category = factor(category, levels = sort(unique(df_plot$category))))
g <- ggplot(df_plot, aes(x = Perturbation, y = category)) +
    geom_point(aes(size = -log10(FDR_enriched),
                   color = significant)) +
    scale_color_manual(values = c("black", "red"), guide = FALSE) +
    scale_y_discrete(limits = rev(levels(df_plot$category)),
                     name = "mSigDB Hallmark Gene Set") +
    scale_size_continuous(limits = c(0, 5),
#                          oob = squish,
                          name = "-log10(FDR)") +
    coord_equal() +
    theme_classic(base_size = 14) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                     margin = margin(t = 6, b = 5)),
          axis.text.y = element_text(margin = margin(r = 6, l = 4)),
          legend.position = "bottom",
          legend.justification = "center")
ggsave(plot = g, filename = paste0(resultsDir, "/gene-set-analysis/hallmarks_vsRenilla_up_v2.pdf"),
       height = 15,
       width = 8,
       device = "pdf")

df_plot <- df_GO %>%
    mutate(Perturbation = gsub("-", " ", Perturbation)) %>%
    mutate(category = gsub("HALLMARK_", "", category)) %>%
    mutate(category = gsub("_", " ", category)) %>%        
#    filter(direction == "down") %>%
    filter((Control == "Renilla")) %>%
    mutate(significant = FDR_unenriched < 0.05)
df_plot <- df_plot %>% mutate(category = factor(category, levels = rev(sort(unique(df_plot$category)))))
g <- ggplot(df_plot, aes(x = Perturbation, y = category)) +
    geom_count(aes(size = -log10(FDR_unenriched),
                   color = significant)) +
    scale_color_manual(values = c("black", "blue"), guide = FALSE) +
    labs(y = "mSigDB Hallmark Gene Set", size = "-log10(FDR)") +    
    theme_classic()
ggsave(plot = g, filename = paste0(resultsDir, "/gene-set-analysis/hallmarks_vsRenilla_down_v2.pdf"),
       height = 15,
       width = 10,
       device = "pdf")

## SKIPPED
## Only show significantly regulated categories
df_plot <- df_GO %>%
    mutate(Perturbation = gsub("-", " ", Perturbation)) %>%
    mutate(category = gsub("HALLMARK_", "", category)) %>%
    mutate(category = gsub("_", " ", category)) %>%    
    filter(direction == "up") %>%
    filter((Control == "Renilla")) %>%
    mutate(significant = FDR < 0.05) %>%
    group_by(category) %>%
    mutate(significant_category = max(significant)) %>%
    ungroup() %>%
    filter(significant_category == 1) %>%
    mutate(logFDR = log10(FDR))
df_plot <- df_plot %>% mutate(category = factor(category, levels = rev(sort(unique(df_plot$category)))))
g <- ggplot(df_plot, aes(x = Perturbation, y = category)) +
    geom_count(aes(size = -logFDR,
                   color = significant)) +
    scale_color_manual(values = c("black", "red"), guide = FALSE) +
#    scale_size(limits = c(NA, 5), oob=censor) +
    labs(y = "mSigDB Hallmark Gene Set", size = "-log10(FDR)") +
    theme_classic()
ggsave(plot = g, filename = paste0(resultsDir, "/gene-set-analysis/hallmarks_vsRenilla_selected_up.pdf"),
       height = 8,
       width = 9,
       device = "pdf")

df_plot <- df_GO %>%
    mutate(Perturbation = gsub("-", " ", Perturbation)) %>%
    mutate(category = gsub("HALLMARK_", "", category)) %>%
    mutate(category = gsub("_", " ", category)) %>%    
    filter(direction == "down") %>%
    filter((Control == "Renilla")) %>%
    mutate(significant = FDR < 0.05) %>%
    group_by(category) %>%
    mutate(significant_category = max(significant)) %>%
    ungroup() %>%
    filter(significant_category == 1)
df_plot <- df_plot %>% mutate(category = factor(category, levels = rev(sort(unique(df_plot$category)))))
g <- ggplot(df_plot, aes(x = Perturbation, y = category)) +
    geom_count(aes(size = -log10(FDR),
                   color = significant)) +
    scale_color_manual(values = c("black", "blue"), guide = FALSE) +
    labs(y = "mSigDB Hallmark Gene Set") +    
    theme_classic()
ggsave(plot = g, filename = paste0(resultsDir, "/gene-set-analysis/hallmarks_vsRenilla_selected_down.pdf"),
       height = 8,
       width = 9,
       device = "pdf")

## Only KRAS Signaling up/down for figure panel

perturbation_names <- c("RIT1\nWT", "RIT1\nM90I", "KRAS\nWT", "KRAS\nG12V", "KRAS\nQ61H")
df_plot <- df_GO %>%
    mutate(Perturbation = gsub("-", "\n", Perturbation)) %>%
    mutate(category = gsub("HALLMARK_", "", category)) %>%
    mutate(category = gsub("_", " ", category)) %>%    
    filter((Control == "Renilla")) %>%
    mutate(significant = FDR_enriched < 0.05) %>%
    filter(grepl("KRAS", category))
df_plot$Perturbation <- factor(df_plot$Perturbation, levels = perturbation_names)
g <- ggplot(df_plot, aes(x = Perturbation, y = category)) +
    geom_count(aes(size = -log10(FDR_enriched)),
               shape = 19,
               color = "#000000") +
    scale_size_area(limits = c(0, 4), oob = squish,
                    max_size = 10,
                    name = NULL) +
    scale_x_discrete(name = NULL,
                     breaks = c(perturbation_names)) +
    scale_y_discrete(name = NULL,
#                     breaks = c("KRAS SIGNALING DN",
#                                "KRAS SIGNALING UP"),
                     labels = c("KRAS\nSignaling\nDown",
                                "KRAS\nSignaling\nUp")
                     ) +
    coord_fixed() +
    theme_classic(base_size = 20) +
    theme(axis.text.x = element_text(
                                     margin = margin(t = 7)),
          axis.text.y = element_text(margin = margin(r = 7))
          )
ggsave(plot = g, filename = paste0(resultsDir, "/gene-set-analysis/hallmark-KRAS-vsRenilla.pdf"),
       height = 4,
       width = 7,
       useDingbats = FALSE,
       device = "pdf")



gene_set <- fread(file = "../gene-sets/HALLMARK_KRAS_SIGNALING_UP.txt", skip=1)
names(gene_set) <- c("Geneid")

foo <- df_DE %>%
    inner_join(gene_set) %>%
    filter(Perturbation == "KRAS-G12V") %>%
    filter(FDR < 0.05 & abs(logFC) > 1)

bar <- df_DE %>%
    inner_join(gene_set) %>%
    filter(Perturbation == "KRAS-Q61H") %>%
    filter(FDR < 0.05 & abs(logFC) > 1)

df_muts <- inner_join(foo, bar, by = c("Geneid", "Control", "Chr", "Start", "End", "Strand", "Length"))
write.csv(df_muts, paste0(resultsDir, "/gene-set-analysis/genes-up-KRAS-muts-HALLMARK-KRAS-UP.csv"))

bar <- df_DE %>%
    inner_join(gene_set) %>%
    filter(Perturbation == "KRAS-WT") %>%
    filter(FDR < 0.05 & abs(logFC) > 1)

df_wtplus <- inner_join(df, bar, by = c("Geneid", "Control", "Chr", "Start", "End", "Strand", "Length"))
write.csv(df_wtplus, paste0(resultsDir, "/gene-set-analysis/genes-up-KRAS-all-HALLMARK-KRAS-UP.csv"))

#######################################################################
#### Spot check known KRAS-associated genes
#######################################################################

gene_list_kras <- c("KLF4", "ADAM8", "LIF", df_muts$Geneid[1:25])

## From edgeR logFC calculations
for (gene in gene_list_kras) {
    df_plot <- df_DE %>%
        filter(Geneid == gene) %>%
        filter(Control == "Renilla") %>%
        mutate(Perturbation = gsub("-", "\n", Perturbation))
    df_plot$Perturbation <- factor(df_plot$Perturbation, levels = perturbation_names)
    ## df_plot$Perturbation <- as.factor(df_plot$Perturbation)
    ## perturb_names <- gsub("-", "\n", levels(df_plot$Perturbation))
    g <- ggplot(df_plot, aes(Perturbation, logFC)) +
        geom_bar(stat = "identity", width = 0.7) +
        scale_x_discrete(name = NULL) +
        scale_y_continuous(name = bquote("log"[2]*"FC"), expand = c(0,0)) +
        ggtitle(gene) +
        theme_classic(base_size = 26) +
        theme(plot.title = element_text(hjust = 0.5, vjust = 0.5,
                                        margin = margin(b = 20)),,
              axis.ticks.x = element_blank(),
              axis.text.x = element_text(margin = margin(t = 5)),
              axis.text.y = element_text(margin = margin(r = 7)))              
    ggsave(plot = g, filename = paste0("single-gene/", gene, "-logFC.pdf"),
           device = "pdf",
           width = 6.5, height = 6,
           path = resultsDir)
}


plot_expression_bars <- function(gene) {
    perturbations <- c("Renilla", "RIT1-WT", "RIT1-M90I", "KRAS-WT", "KRAS-G12V", "KRAS-Q61H")
    df_plot <- df_tpm %>%
        filter(Gene==gene) %>%
        mutate(gene=grepl(gene, Sample)) %>%
        mutate(Group = factor(grouping, levels = perturbations)) %>%
        mutate(GenePerturb = strsplit2(as.character(Group), split = "-")[,1]) %>%
        mutate(GenePerturb = factor(GenePerturb, levels = c("Renilla", "KRAS", "RIT1")))
    df_plot <- arrange(df_plot, Group)
    df_points <- df_plot
    df_mean <- df_plot %>%
        group_by(Group) %>%
        summarize(meanLogTPM = mean(log2(TPM)))
    df_sd <- df_plot %>%
        group_by(Group) %>%
        summarize(sdLogTPM = sd(log2(TPM)))
    df_plot <- inner_join(df_mean, df_sd) %>%
        mutate(lower = meanLogTPM - sdLogTPM) %>%
        mutate(upper = meanLogTPM + sdLogTPM)      
    print(df_plot)
    print(df_points)
    group_names <- gsub("-", "\n", levels(df_plot$Group))
    group_names <- gsub("Renilla", "Vector", group_names)
    print(max(log2(df_points$TPM)))
    g <- ggplot(df_plot, aes(x=Group, y=meanLogTPM)) +
        geom_col(width = 0.7, fill = darken("#999999", 0.1)) +
        geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.3) +
        geom_point(data = df_points, aes(x = Group, y = log2(TPM)),
                   size = 2, alpha = 0.8) +      
        scale_color_manual(name = NULL,
                           values = darken(colors, 0.1)) +
        scale_x_discrete(labels = group_names,
                          name = NULL) +
        scale_y_continuous(name = bquote("log"[2]*"TPM" ~ of ~ .(gene)),
                           limits = c(0, max(log2(df_points$TPM))+1),
                           breaks = c(0, 2, 4, 6, 8),
                           expand = c(0,0)) +
        guides(fill = FALSE, color = FALSE) +
        theme_classic(base_size = 26) +
        theme(plot.title = element_text(size = 16, hjust = 0.5),
              axis.ticks.x = element_blank(),
              axis.text.x = element_text(#angle = 75, hjust = 1, vjust = 1,
                  margin = margin(t = 4, b = 10)),
              axis.text.y = element_text(margin = margin(r = 4, l = 8)))
    return(g)
}

for (gene in gene_list_kras) {
    g <- plot_expression_bars(gene)
    ggsave(g, filename = paste0("single-gene/", gene, "-logTPM.pdf"),
           device = "pdf",
           width = 7, height = 6,
           path = resultsDir,
           useDingbats = FALSE)
}

## Statistical t-test to directly compare log2TPM values of control vs KRAS perturbed samples
gene_list_kras <- c("ETS1", "HBEGF", "LIF")

for (gene in gene_list_kras) {
    perturbations <- c("Renilla", "RIT1-WT", "RIT1-M90I", "KRAS-WT", "KRAS-G12V", "KRAS-Q61H")    
    df_test <- df_tpm %>%
        filter(Gene == gene) %>%
        mutate(gene = grepl(gene, Sample)) %>%
        mutate(Group = factor(grouping, levels = perturbations)) %>%        
        mutate(GenePerturb = strsplit2(as.character(Group), split = "-")[,1]) %>%
        mutate(GenePerturb = factor(GenePerturb, levels = c("Renilla", "KRAS", "RIT1")))
#    print(df_test)
    for (perturb in perturbations[-1]) {
        df_control <- df_test %>% filter(Group == "Renilla")
        df <- df_test %>% filter(Group == perturb)
        print(c(gene, perturb))
        print(t.test(log2(df_control$TPM), log2(df$TPM))$p.value)
    }
}

## Direct from expression levels
for (gene in gene_list_kras) {
    renilla_expr <- df_expr %>%
        filter(Gene == gene) %>%
        filter(group == "Renilla")
    renilla_expr <- mean(renilla_expr$logTPM)
    df_plot <- df_expr %>%
        filter(Gene == gene) %>%
        mutate(group = gsub("-", "\n", group)) %>%
        mutate(logFC = (logTPM - renilla_expr)/renilla_expr) %>%
        filter(!group=="Renilla") %>%
        group_by(group) %>%
        summarize(log2FC = mean(logFC))
    df_plot$Perturbation <- as.factor(df_plot$Perturbation)
    perturb_names <- gsub("-", "\n", levels(df_plot$Perturbation))
    g <- ggplot(df_plot, aes(reorder(Perturbation, -logFC), logFC)) +
        geom_bar(stat = "identity", width = 0.7) +
        scale_x_discrete(name = NULL) +
        scale_y_continuous(name = bquote("log"[2]*"FC"), expand = c(0,0)) +
        ggtitle(gene) +
        theme_classic(base_size = 24) +
        theme(plot.title = element_text(hjust = 0.5, vjust = 0.5),
              axis.ticks.x = element_blank())
    ggsave(plot = g, filename = paste0("single-gene-logFC/", gene, "-logFC.pdf"),
           device = "pdf",
           width = 6, height = 6,
           path = resultsDir)
}


#######################################################################
#### Gene set analysis with goseq for RIT1 vs KRAS-muts
#######################################################################

load(file = paste0(dataDir, "df_DE_KRASvRIT1.RData"))

df_GO_KRASvRIT1 <- data.frame(Test = character(),
                    category = character(),
                    direction = character(),
                    over_represented_pvalue = double(),
                    under_represented_pvalue = double(),
                    numDEInCat = integer(),
                    numInCat = integer(),
                    term = character(),
                    ontology = character())
## Do the analysis
for (direction in c("up", "down")) {
    ##    for (test in c("GO", "MF", "BP", "hallmarks")) {
    for (test in c("hallmarks")) {    
        d <- do_goseq_DE(df_DE_KRASvRIT1, "RIT1-all", "KRAS-muts", direction = direction, test = test)
        d <- mutate(d, FDR = p.adjust(over_represented_pvalue, method = "BH"))        
        d$direction <- direction
        d$Test <- test
        df_GO_KRASvRIT1 <- bind_rows(df_GO_KRASvRIT1, d)
    }
}

## Save the top results
for (direction in c("up", "down")) {
    for (test in c("GO", "MF", "BP", "hallmarks")) {
        d <- df_GO_KRASvRIT1 %>% filter(direction == direction & Test == test)
        write.csv(d, file = paste0(resultsDir, "/gene-set-analysis/RIT1-all_KRAS-muts_", direction, "_", test, ".csv"),
                  quote = FALSE)
    }
}

## Plot the top results
for (dir in c("up", "down")) {
    for (test in c("MF", "BP")) {
        df_plot <- df_GO_KRASvRIT1 %>% filter(direction == dir & Test == test) %>%
            dplyr::slice(1:10)
        print(df_plot)
        g <- ggplot(df_plot) +
            geom_col(aes(x = reorder(category, -FDR), y = -log10(FDR)),
                     width = 0.7) +
            coord_flip() +
            scale_x_discrete(name = "GO category term") +
            scale_y_continuous(expand = c(0,0)) +
            theme_classic(base_size = 20)
        ggsave(g, file = paste0("/gene-set-analysis/RIT1-all_KRAS-muts_", direction, "_", test, ".pdf"),
               height = 5, width = 9,
               path = resultsDir)
    }
    df_plot <- df_GO_KRASvRIT1 %>% filter(direction == dir & Test == "hallmarks") %>%
        dplyr::slice(1:10)
    print(direction)
    print(df_plot)
    g <- ggplot(df_plot) +
        geom_col(aes(x = reorder(category, -FDR), y = -log10(FDR)),
                 width = 0.7) +
        coord_flip() +
        scale_x_discrete(name = "mSigDB Hallmark") +
        scale_y_continuous(expand = c(0,0)) +
        theme_classic(base_size = 20)
    ggsave(g, file = paste0("/gene-set-analysis/RIT1-all_KRAS-muts_", direction, "_", "hallmarks", ".pdf"),
           height = 5, width = 9,
           path = resultsDir)
}

df_plot <- df_GO_KRASvRIT1 %>%
    filter(Test == "hallmarks")
df_plot <- as_tibble(df_plot) %>%
    mutate(category = gsub("HALLMARK_", "", category)) %>%
    mutate(category = gsub("_", " ", category))
g <- ggplot(df_plot) +
    geom_count(aes(x = direction,
                   y = category,
                   size = -log(FDR))) +
    scale_size_area(limits = c(0, 4), oob = squish,
                    max_size = 10,
                    name = NULL) +
    coord_fixed() +
    theme_classic()
ggsave(g, file = paste0("/gene-set-analysis/RIT1-all_KRAS-muts_hallmarks.pdf"),
       height = 12, width = 6,
       path = resultsDir)

df_plot <- df_GO_KRASvRIT1 %>%
    filter(Test == "hallmarks") %>%
    filter(category %in% c("HALLMARK_INTERFERON_ALPHA_RESPONSE", "HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_INFLAMMATORY_RESPONSE"))

## up = up in RIT1 = down in KRAS, down = down in RIT1 = up in KRAS
df_plot <- as_tibble(df_plot) %>%
    mutate(category = gsub("HALLMARK_", "", category)) %>%
    mutate(category = gsub("_", " ", category))
g <- ggplot(df_plot) +
    geom_count(aes(x = direction, y = category, size = -log(FDR))) +
    coord_fixed() +
    theme_classic()
ggsave(g, file = paste0("/gene-set-analysis/RIT1-all_KRAS-muts_hallmarks_inflammatory.pdf"),
       height = 6, width = 6,
       path = resultsDir)

#### Interferon Gamma Response

ifng_gene_set <- fread(file = "../gene-sets/HALLMARK_INTERFERON_GAMMA_RESPONSE.txt", skip=1)
names(ifng_gene_set) <- "Gene"

df_ifng <- df_DE_KRASvRIT1 %>%
    filter(Geneid %in% ifng_gene_set$Gene) %>%
    filter(FDR < 0.05) %>%
    filter(abs(logFC) > 1)

#### Interferon Alpha Response

ifna_gene_set <- fread(file = "../gene-sets/HALLMARK_INTERFERON_ALPHA_RESPONSE.txt", skip=1)
names(ifna_gene_set) <- "Gene"

df_ifna <- df_DE_KRASvRIT1 %>%
    filter(Geneid %in% ifna_gene_set$Gene) %>%
    filter(FDR < 0.05) %>%
    filter(abs(logFC) > 1)

#### Inflammatory Response

inflammatory_gene_set <- fread(file = "../gene-sets/HALLMARK_INFLAMMATORY_RESPONSE.txt", skip=1)
names(inflammatory_gene_set) <- "Gene"

df_inflam <- df_DE_KRASvRIT1 %>%
    filter(Geneid %in% inflammatory_gene_set$Gene) %>%
    filter(FDR < 0.05) %>%
    filter(abs(logFC) > 1)

#######################################################################
#### Compare with proteomics data
#######################################################################

load(file = "../data/phosphoproteomics/df_proteome.RData")

plot_KRASvRIT1_DE <- function(df, genesetname) {
    df_plot <- df %>%
        mutate(Geneid = factor(Geneid, levels = df$Geneid)) %>%
        filter(Geneid %in% df_proteome$geneSymbol)

    g <- ggplot(df_plot, aes(x = Geneid, y = logFC)) +
        geom_col(width = 0.7) +
        theme_bw(base_size = 24) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
    ggsave(g, file = paste0("/", genesetname, "-RNA-DE.pdf"),
           height = 8, width = 18,
           device = "pdf",
           path = resultsDir)
}

plot_KRASvRIT1_DE(df_ifng, "IFNG-hallmark")
plot_KRASvRIT1_DE(df_ifna, "IFNa-hallmark")
plot_KRASvRIT1_DE(df_inflam, "inflammatory-response-hallmark")

plot_protein_DE <- function(df, genesetname) {
    df_plot <- df_proteome %>%
        filter(geneSymbol %in% df$Geneid) %>%
        mutate(geneSymbol = factor(geneSymbol, levels = df$Geneid)) %>%
        filter(Control == "Renilla") %>%
        mutate(perturbation = paste0(perturb_gene, "_", perturb_variant))

    g <- ggplot(df_plot, aes(x = geneSymbol, y = logFC, fill = perturbation)) +
        geom_col(width = 0.7,
                 position = position_dodge2(padding=0.1)) +
        ## guides(
        ##     fill = guide_legend(
        ##     )) + 
        theme_bw(base_size = 24) +
        theme(
            legend.position = "top",
            legend.justification = "center",
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
    ggsave(g, file = paste0("/", genesetname, "-proteome.pdf"),
           height = 8, width = 18,
           device = "pdf",
           path = resultsDir)
}

plot_protein_DE(df_ifng, "IFNG-hallmark")
plot_protein_DE(df_ifna, "IFNa-hallmark")
plot_protein_DE(df_inflam, "inflammatory-response-hallmark")

plot_expression_bars_multiple_genes <- function(gene_list) {
    perturbations <- c("Renilla", "RIT1-WT", "RIT1-M90I", "KRAS-WT", "KRAS-G12V", "KRAS-Q61H")
    perturbation_colors <- c("#777777", "#D55E00", "#E69F00", "#0072B2", "#009E73", "#56B4E9")        
    df_plot <- df_expr %>%
        filter(Gene %in% gene_list) %>%
        mutate(GenePerturb = strsplit2(as.character(group), split = "-")[,1]) %>%
        mutate(GenePerturb = factor(GenePerturb, levels = c("Renilla", "KRAS", "RIT1")))
    df_plot <- arrange(df_plot, group)
    df_points <- df_plot
    df_mean <- df_plot %>%
        group_by(.dots = c("Gene", "group")) %>%
        summarize(meanLogTPM = mean(logTPM))
    df_sd <- df_plot %>%
        group_by(.dots = c("Gene", "group")) %>%
        summarize(sdLogTPM = sd(logTPM))
    df_plot <- inner_join(df_mean, df_sd) %>%
        mutate(lower = meanLogTPM - sdLogTPM) %>%
        mutate(upper = meanLogTPM + sdLogTPM)      
#    colors <- rep("#999999", length(unique(df_plot$Sample)))
    group_names <- gsub("-", "\n", levels(df_plot$group))
    group_names <- gsub("Renilla", "Vector", group_names)
    print(max(df_points$TPM))
    g <- ggplot(df_plot, aes(x=Gene, y=meanLogTPM, fill = group)) +
        geom_col(width = 0.7,
                 position = position_dodge2(padding=0.1)) +
        geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.3,
                      position = position_dodge2(padding=0.9)) +
        ## geom_boxplot(data = df_points, aes(x = Group, y = log2(TPM)),
        ##             width = 0.3) +
        ## geom_point(data = df_points, aes(x = Group, y = log2(TPM)),
        ##            size = 2, alpha = 0.8) +      
        #      geom_col(position = position_dodge2(padding = 0.1), width = 0.8) +
        scale_fill_manual(name = NULL,
                          values = perturbation_colors,
                          labels = perturbations) +
        scale_x_discrete(labels = group_names,
                         name = NULL) +
        scale_y_continuous(name = bquote("log"[2]*"TPM"),
                           limits = c(0, max(df_points$logTPM)+1),
#                           breaks = c(0, 2, 4, 6, 8),
                           expand = c(0,0)) +
        guides(color = FALSE) +
        theme_classic(base_size = 26) +
        theme(plot.title = element_text(size = 16, hjust = 0.5),
              axis.ticks.x = element_blank(),
              axis.text.x = element_text(#angle = 75, hjust = 1, vjust = 1,
                  margin = margin(t = 4, b = 10)),
              axis.text.y = element_text(margin = margin(r = 4, l = 8)))
    return(g)
}


plot_DE_bars_multiple_genes <- function(gene_list) {
    perturbations <- c("RIT1-WT", "RIT1-M90I", "KRAS-WT", "KRAS-G12V", "KRAS-Q61H")
    perturbation_colors <- c("#D55E00", "#E69F00", "#0072B2", "#009E73", "#56B4E9")    
    df_plot <- df_DE %>%
        filter(Control == "Renilla") %>%
        filter(Geneid %in% gene_list) %>%
        mutate(GenePerturb = strsplit2(as.character(Perturbation), split = "-")[,1]) %>%
        mutate(GenePerturb = factor(GenePerturb, levels = c("KRAS", "RIT1"))) %>%
        mutate(Perturbation = factor(Perturbation, levels = c("RIT1-WT", "RIT1-M90I", "KRAS-WT", "KRAS-G12V", "KRAS-Q61H")))
#    print(df_plot)
#    df_plot <- arrange(df_plot, Group)
#    colors <- rep("#999999", length(unique(df_plot$Sample)))
    df_plot$Geneid <- factor(df_plot$Geneid, levels = levels(gene_list))
    group_names <- gsub("-", " ", perturbations)
    g <- ggplot(df_plot, aes(x=Geneid, y=logFC, fill=Perturbation)) +
        geom_col(width = 0.7,
                 position = position_dodge2(padding=0.1)) +
        #      geom_col(position = position_dodge2(padding = 0.1), width = 0.8) +
        scale_fill_manual(name = NULL,
                           values = perturbation_colors,
                           labels = group_names) +
        scale_x_discrete(name = NULL,
                         breaks = gene_list) +
        scale_y_continuous(name = bquote("log"[2]*"FC"),
                            limits = c(min(df_plot$logFC)-1, max(df_plot$logFC)+1),
                            ) +
        guides(
            fill = guide_legend(
                nrow = 1)
        ) +
        theme_classic(base_size = 28) +
        theme(legend.position = "top",
              legend.justification = "center",
              legend.text = element_text(size = 24),              
              plot.title = element_text(size = 16, hjust = 0.5),
              axis.ticks.x = element_blank(),
              axis.text.x = element_text(#angle = 75, hjust = 1, vjust = 1,
                  margin = margin(t = 4, b = 10)),
              axis.text.y = element_text(margin = margin(r = 4, l = 8)))
    return(g)
}

#######################################################################
#### Inflammatory response genes expression?
#######################################################################

inflammatory_genes <- inflammatory_gene_set$Gene[inflammatory_gene_set$Gene %in% df_proteome$geneSymbol]

g <- plot_DE_bars_multiple_genes(inflammatory_genes[1:5])
ggsave(g, filename = "/diff-expression/inflammatory_genes.pdf",
       height = 8, width = 14,
       path = resultsDir)

#######################################################################
#### MHC I expression?
#######################################################################

HLA_genes <- factor(c("KRAS", "RIT1", "VIM", "HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F", "NLRC5", "IRF1", "IRF2"),
                    levels = c("KRAS", "RIT1", "VIM", "HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F", "NLRC5", "IRF1", "IRF2"))
g <- plot_expression_bars_multiple_genes(HLA_genes)
ggsave(g, filename = "/diff-expression/HLA_genes.pdf",
       height = 8, width = 14,
       path = resultsDir)
g <- plot_DE_bars_multiple_genes(HLA_genes)
ggsave(g, filename = "/diff-expression/HLA_genes-DE.pdf",
       height = 8, width = 14,
       path = resultsDir)

MHC_genes <- factor(c("VIM", "KRAS", "RIT1", "NLRC5", "IRF1", "IRF2"),
                    levels = c("KRAS", "RIT1", "VIM", "IRF1", "IRF2", "NLRC5"))
g <- plot_expression_bars_multiple_genes(MHC_genes)
ggsave(g, filename = "/diff-expression/MHC_genes.pdf",
       height = 8, width = 10,
       path = resultsDir)
g <- plot_DE_bars_multiple_genes(MHC_genes)
ggsave(g, filename = "/diff-expression/MHC_genes-DE.pdf",
       height = 8, width = 10,
       path = resultsDir)

EMT_genes <- factor(c("VIM", "CDH2", "FN1", "KRT19"),
                    levels = c("VIM", "CDH2", "FN1", "KRT19"))
g <- plot_expression_bars_multiple_genes(EMT_genes)
ggsave(g, filename = "/diff-expression/EMT_genes.pdf",
       height = 8, width = 11,
       path = resultsDir)
g <- plot_DE_bars_multiple_genes(EMT_genes)
ggsave(g, filename = "/diff-expression/EMT_genes-DE.pdf",
       height = 8, width = 11,
       path = resultsDir)

