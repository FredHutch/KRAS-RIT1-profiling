#######################################################################
#### Run ChEA3 analysis
#######################################################################

## ChIP-X Enrichment Analysis Version 3
## Prep for and use ChEA3's REST API

#######################################################################
#### Set up
#######################################################################

source("../scripts/R/load.R")
library(httr)
library(jsonlite)

#######################################################################
#### Load DE data/results
#######################################################################

load("df_DE.RData")

perturbations <- unique(df_DE$Perturbation)

#######################################################################
#### Set relevant input and output directories
#######################################################################

## Inputs
dataDirs <- c("./")
## Outputs
dataDir <- "../results/AALE_KRAS-RIT1/ChEA"
resultsDir <- "../results/AALE_KRAS-RIT1/ChEA"
figsDir <- "../results/AALE_KRAS-RIT1/Figures"

#######################################################################
#### Variants vs Renilla for ChEA TF analysis
#######################################################################

## Generate and save gene lists for entering into web portal
## DE of RIT1 M90I vs Renilla
foo <- df_DE %>%
    filter(Perturbation == "RIT1-M90I") %>%
    filter(Control == "Renilla") %>%
    filter(abs(logFC) > 1) %>%
    filter(FDR < 0.05)
write(foo$Geneid, file = paste0(resultsDir, "/diff-expression/RIT1-M90I_Renilla_forChEA.txt"))

bar <- df_DE %>%
    filter(Perturbation == "KRAS-G12V") %>%
    filter(Control == "Renilla") %>%
    filter(abs(logFC) > 1) %>%
    filter(FDR < 0.05)
write(bar$Geneid, file = paste0(resultsDir, "/diff-expression/KRAS-G12V_Renilla_forChEA.txt"))

#### This doesn't work
## DE of KRAS G12V and Q61H vs Renilla
KRAS_idx <- c(which(grepl("KRAS-G12V", colnames(design))),
              which(grepl("KRAS-Q61H", colnames(design))))
control_idx <- c(which(grepl("Renilla", colnames(design))))
contrast <- rep(0, length(colnames(design)))
contrast[KRAS_idx] <- 1
contrast[control_idx] <- -1
print(contrast)
qlf <- glmQLFTest(fit, contrast=contrast)
tt <- topTags(qlf, n=dim(fit)[1], sort.by="PValue")
sigTestResults <- decideTests(qlf, p.value = 0.01, adjust.method = "BH")
print(head(sigTestResults))
print(summary(sigTestResults))

#######################################################################
#### Each vs Renilla for ChEA TF analysis
#######################################################################

df_DE_filtered <- df_DE %>%
    filter(Control == "Renilla") %>%
    filter(abs(logFC) > 1 ) %>%
    filter(FDR < 0.05)

gene_lists <- list()
for (perturb in perturbations) {
    df <- df_DE_filtered %>%
        filter(Perturbation == perturb)
    print(length(df$Geneid))
    print(head(df$Geneid))
    url = "https://amp.pharm.mssm.edu/chea3/api/enrich/"
    encode = "json"
    payload = list(query_name = paste0(perturb, "_query"), gene_set = df$Geneid)
    ## POST to ChEA3 server
    response = POST(url = url,
                    body = payload,
                    encode = encode)
    json = httr::content(response, as = "text")
    ## Convert results to a list of R dataframes
    results <- fromJSON(json)
    print(length(results))
    for (analysis in names(results)) {
        df_save <- results[[analysis]]
#        print(head(df_save))
        # print(gsub("--", "_", analysis))
        filename <- paste0(resultsDir, "/", perturb, "_Renilla_", gsub("--", "_", analysis), ".tsv")
        write.table(df_save, file = filename, sep = "\t", quote = FALSE,
                    row.names = FALSE)
    }
}
