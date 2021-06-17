source("../scripts/R/load.R")
library(Rsamtools)

resultsDir <- "../results/AALE_KRAS-RIT1"

df_mut <- data.frame(Depth = integer(),
                     mut = integer(),
                     Sample = character())
#for (sample in rownames(sample_summaries)) {
for (sample in c("KRAS-G12V_1", "KRAS-G12V_2", "KRAS-G12V_3")) {    
    mpileup <- fread(file = paste0("./STAR/hg19/coverages/", sample, ".mpileup.txt"))
    names(mpileup) <- c("Chromosome", "Pos", "RefBase", "Depth", "ReadBases", "ReadBaseQualities")
    mpileup <- as_tibble(mpileup)
#    print(mpileup[2,4]) # Get total depth
    n_muts <- str_count(mpileup[2, 5], "a")
    df <- data.frame(Sample = sample,
                     Depth = mpileup[2,4],
                     mut = n_muts)
    df_mut <- bind_rows(df_mut, df)
}

#for (mut_oe in c("KRAS-G12V", "KRAS-Q61H", "RIT1-M90I")) {
    mpileup <- fread(file = paste0("./STAR/hg19/coverages/", "KRAS-RIT1-muts", ".txt"))
    names(mpileup) <- c("Chromosome", "Pos", "RefBase", "Depth", "ReadBases", "ReadBaseQualities", "bamfile")
    df_mut <- as_tibble(mpileup) %>%
        mutate(mut = str_count(ReadBases, "a")) %>%
        mutate(ref = Depth - mut) %>%
        gather("allele", "n_reads", c(ref, mut))
#}

for (variant in c("KRAS-G12V", "KRAS-Q61H", "RIT1-M90I")) {
    df_plot <- as_tibble(df_mut) %>%
        filter(grepl(variant, bamfile)) %>%
        mutate(mut = str_count(ReadBases, "[^[:punct:]]")) %>%    
        mutate(ref = Depth - mut) %>%
        # gather("allele", "n_reads", c(ref, mut)) %>%
        mutate(percent_mut = mut/Depth) %>%
        mutate(percent_ref = ref/Depth) %>%
        gather("allele", "percent_reads", c(percent_ref, percent_mut)) %>%
        dplyr::select(-c("mut", "ref", "n_reads")) %>%
        distinct()
    print(variant)
    g <- ggplot(df_plot, aes(x = bamfile, y = percent_reads, fill = allele)) +
        geom_bar(stat = "identity", width = 0.6) +
        scale_fill_manual(values = c("#009E73", "grey"),
                          name = NULL) +
        scale_x_discrete(labels = c("1", "2", "3"),
                         name = "Replicate") +
        scale_y_continuous(name = "% of Reads",
                           labels = c("0", "25", "50", "75", "100"),       
                           expand = c(0,0)) +
        ggtitle(gsub("-", " ", variant)) +
        theme_classic(base_size = 26) +
        theme(legend.position = "none",
              plot.title = element_text(hjust = 0.5,
                                        margin = margin(b = 20)),
              axis.text.x = element_text(margin = margin(t = 4, b = 6)),
              axis.text.y = element_text(margin = margin(r = 7)),
              axis.ticks.x = element_blank())
    ggsave(g, file = paste0("alleles/", variant, ".pdf"),
           device = "pdf",
           width = 5, height = 5,
           path = resultsDir)
}

g <- ggplot(df_M90I, aes(x = bamfile, y = percent_reads, fill = allele)) +
    geom_bar(stat = "identity", width = 0.6) +
    scale_fill_manual(values = c("#009E73", "grey"),
                      labels = c("Mutant", "Wild Type"),
                      name = NULL) +
    scale_x_discrete(labels = c("1", "2", "3"),
                     name = "Replicate") +
    scale_y_continuous(name = "Allele percentage",
                       labels = c("0", "25", "50", "75", "100"),                    
                       expand = c(0,0)) +
    ggtitle("RIT1 M90I") +
    theme_classic(base_size = 24) +
    theme(plot.title = element_text(hjust = 0.5,
                                    margin = margin(b = 20)))
ggsave(g, file = "alleles/alleles_legend.pdf",
       device = "pdf",
       width = 5, height = 6,
       path = resultsDir)


mpileup <- readPileup(file = "./STAR/hg19/coverages/KRAS-G12V.pileup.vcf", variant = "SNP")
