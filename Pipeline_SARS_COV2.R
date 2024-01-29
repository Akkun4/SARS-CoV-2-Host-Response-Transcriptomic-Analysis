library(GOplot)
data = read.delim("analysis.csv")

barplot(
  data$upload_1..fold.Enrichment,
  names.arg = data$GO.biological.process.complete,
  col = "skyblue",
  legend.text = c("Value1", "Value2"),
  xlab = "Pathways",
  ylab = "Enrichment score",
  main = "GO Pathway"
)

ggplot(data, aes(x = data$GO.biological.process.complete, y = data$upload_1..fold.Enrichment, fill = data$GO.biological.process.complete)) +
  geom_bar(stat = "identity", color = "black", width = 0.7) +
  labs(
    title = "GO Pathway",
    x = "Pathway",
    y = "Fold Enrichment"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"  # Remove legend
  )



# Install and load required packages
# install.packages(c("data.table", "ggplot2"))
library(data.table)
library(ggplot2)

# Sample dataset with real gene names
genes <- c(
  "BRCA1", "TP53", "EGFR", "KRAS", "MYC",
  "PTEN", "RB1", "CDKN2A", "PIK3CA", "ERBB2",
  "VHL", "FLT3", "KIT", "MET", "BRAF", "PRKCG", "CXCL3", "CXCL1", "CXCL5", "TNFSF15", "CSF3", "TRAF1", "CXCL16", "SLC40A1",
  "MT1E", "MT1X", "ZCCHC2", "BHLHE41", "DRAM1", "NUAK2", "ROR1-AS1", "RSPO3", "TUBB3",
  "IL32", "OLR1", "SGIP1", "EIF3C", "TMEM158", "NEDD4L", "TNFAIP2", "KYNU", "RGPD5",
  "TRIM47", "FMN1", "H2AC6", "USP12", "CCDC50", "CAMK4", "TRIM36", "ZFP36L2", "PAN3",
  "PLEKHH2", "PFKM", "TMEM123", "IGSF10", "SPEF2", "SPARCL1", "MBNL1", "NADK2", "GPD1L",
  "GJA1", "CCNO", "SLC30A6", "PELO", "RASGRP3", "SAR1B", "CATSPER3", "GPR180", "DYNLT5",
  "DNAI4", "ANKRD22", "FARP1", "IFIT5", "SLC16A12", "PRDM8", "BMP3", "SLC", "HHEX",
  "PTPRK", "GGPS1", "PART1", "LMNTD1", "MARVELD2", "RAD17", "MED21", "STK32B", "ZIC1",
  "ADGRA3", "SREK1IP1", "LGI2", "MR1", "SRP19", "CENPH", "CDYL", "CARHSP1", "BANK1",
  "TXNDC11", "DAB2", "ACMSD", "ANAPC1", "CAST", "SCOC", "CLGN", "CETN3", "SMARCA5",
  "BMP6", "RGPD3", "KLHDC7B", "OMP", "FGF20", "MYOM3", "IL20RA", "Y_RNA", "TAC1", "CPAMD8", "MX2", "RTP4", "VEZTP1",
  "CACNA1I", "THBS4-AS1", "NKAIN4", "AARD", "DBH-AS1", "RSAD2", "PRSS21", "CPA5", "TMEM221",
  "TPD52", "BST2", "CIBAR1-DT", "IYD", "CITED1", "RNA5SP174", "IFIT1", "PLAAT2", "SLC8A2",
  "YBX2", "DNAI7", "STK24-AS1", "RN7SKP26", "IFIT2", "MIR4479", "NRARP", "IFI44L", "CRHBP",
  "ARPC3P4", "CXCL3", "CXCL1", "CXCL5", "RN7SL75P", "MX1", "PTGDR2", "RNA5-8SN2",
  "FBXO39", "INSL6", "CXCL11", "CYP21A2", "MAPK8IP2", "SLC38A4", "SLC", "RN7SL574P", "IFI6",
  "DNAJC19P5", "WFIKKN1", "H3C1", "GABRR2", "BMS1P15", "RPL31P62", "TCF7L1-IT1", "FMR1-IT1",
  "SCT", "DOK7", "CMPK2", "CKM", "BATF2", "LINC02988", "THEGL", "TEKT3", "PSG2", "RPL35P6",
  "RPSAP18", "MPIG6B", "KRT20", "MANEA-DT", "RET", "GJB5", "RNU1-72P", "GRHL1", "OAS2", "CPSF4L",
  "PPP1R14BP2", "CNTNAP2", "RNU6-1048P", "SLX1A", "NPM1P13", "RNA5SP435", "CKS1BP3", "HOXB-AS1",
  "DBH", "DHRS9", "SEMA6A-AS1", "CHUK-DT", "CHRM4", "KRT13", "GCA", "RAMP2", "AARD", "NA",
  "LINC02235", "C1orf226", "KLRD1", "BBOX1", "RXFP4", "PRAL", "SLC4A1", "HAL", "RSF1-IT2",
  "ST6GALNAC2", "KCNT1", "EXOC3L1", "MYCBP2-AS1", "RIMS4", "TMEM225B", "FUT1", "RN7SL3", "MIR4284",
  "BSND", "UBQLNL", "NA", "ANGPTL1", "BIRC7", "ANKRD31", "LCMT1-AS2", "DCST1", "PLIN1", "CDKN1C",
  "CXCL10", "ST6GALNAC1", "FASLG", "RPSAP70", "RAB38", "TEX29", "KRT40", "CPA2", "SLC6A4", "DPYS",
  "AK8", "MIR409", "RPL21P32", "TFR2", "GRM2", "KCNJ5"
)

# Generate synthetic dataset with real gene names
nephropathy_data <- data.table(
  Gene = sample(genes, 100),
  Fold_Enrichment = rnorm(100, mean = 2, sd = 0.5),
  Condition = "Nephropathy"
)
retinopathy_data <- data.table(
  Gene = sample(genes, 100),
  Fold_Enrichment = rnorm(100, mean = 3, sd = 0.7),
  Condition = "Retinopathy"
)
synthetic_data <- rbind(nephropathy_data, retinopathy_data)

# Plot the data
ggplot(synthetic_data, aes(x = Gene, y = Fold_Enrichment, fill = Condition)) +
  geom_bar(stat = "identity", position = "dodge", color = "black", width = 0.7) +
  labs(
    title = "Fold Change for Nephropathy and Retinopathy Genes",
    x = "Gene",
    y = "Fold Enrichment"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Install and load required packages
# install.packages(c("data.table", "ggplot2"))
library(data.table)
library(ggplot2)

# Set a seed for reproducibility
set.seed(123)

# Generate synthetic dataset
genes <- c(paste0("Gene", 1:100), paste0("OtherGene", 1:100))
nephropathy_data <- data.table(
  Gene = sample(genes, 20),
  Fold_Enrichment = rnorm(20, mean = 2, sd = 0.5),
  Condition = "Nephropathy"
)
retinopathy_data <- data.table(
  Gene = sample(genes, 20),
  Fold_Enrichment = rnorm(20, mean = 3, sd = 0.7),
  Condition = "Retinopathy"
)
synthetic_data <- rbind(nephropathy_data, retinopathy_data)

# Plot the data
ggplot(synthetic_data, aes(x = Gene, y = Fold_Enrichment, fill = Condition)) +
  geom_bar(stat = "identity", position = "dodge", color = "black", width = 0.7) +
  labs(
    title = "Fold Change for Nephropathy and Retinopathy Genes",
    x = "Gene",
    y = "Fold Change"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))




# Install and load required packages
# install.packages(c("data.table", "ggplot2"))
library(data.table)
library(ggplot2)

# Define gene lists
nephropathy_genes <- c(
  "PRKCG", "CXCL3", "CXCL1", "CXCL5", "TNFSF15", "CSF3", "TRAF1", "CXCL16", "SLC40A1",
  "MT1E", "MT1X", "ZCCHC2", "BHLHE41", "DRAM1", "NUAK2", "ROR1-AS1", "RSPO3", "TUBB3",
  "IL32", "OLR1", "SGIP1", "EIF3C", "TMEM158", "NEDD4L", "TNFAIP2", "KYNU", "RGPD5",
  "TRIM47", "FMN1", "H2AC6", "USP12", "CCDC50", "CAMK4", "TRIM36", "ZFP36L2", "PAN3",
  "PLEKHH2", "PFKM", "TMEM123", "IGSF10", "SPEF2", "SPARCL1", "MBNL1", "NADK2", "GPD1L",
  "GJA1", "CCNO", "SLC30A6", "PELO", "RASGRP3", "SAR1B", "CATSPER3", "GPR180", "DYNLT5",
  "DNAI4", "ANKRD22", "FARP1", "IFIT5", "SLC16A12", "PRDM8", "BMP3", "SLC", "HHEX",
  "PTPRK", "GGPS1", "PART1", "LMNTD1", "MARVELD2", "RAD17", "MED21", "STK32B", "ZIC1",
  "ADGRA3", "SREK1IP1", "LGI2", "MR1", "SRP19", "CENPH", "CDYL", "CARHSP1", "BANK1",
  "TXNDC11", "DAB2", "ACMSD", "ANAPC1", "CAST", "SCOC", "CLGN", "CETN3", "SMARCA5",
  "BMP6", "RGPD3", "KLHDC7B"
)

retinopathy_genes <- c(
  "OMP", "FGF20", "MYOM3", "IL20RA", "Y_RNA", "TAC1", "CPAMD8", "MX2", "RTP4", "VEZTP1",
  "CACNA1I", "THBS4-AS1", "NKAIN4", "AARD", "DBH-AS1", "RSAD2", "PRSS21", "CPA5", "TMEM221",
  "TPD52", "BST2", "CIBAR1-DT", "IYD", "CITED1", "RNA5SP174", "IFIT1", "PLAAT2", "SLC8A2",
  "YBX2", "DNAI7", "STK24-AS1", "RN7SKP26", "IFIT2", "MIR4479", "NRARP", "IFI44L", "CRHBP",
  "ARPC3P4", "CXCL3", "CXCL1", "CXCL5", "RN7SL75P", "MX1", "PTGDR2", "RNA5-8SN2",
  "FBXO39", "INSL6", "CXCL11", "CYP21A2", "MAPK8IP2", "SLC38A4", "SLC", "RN7SL574P", "IFI6",
  "DNAJC19P5", "WFIKKN1", "H3C1", "GABRR2", "BMS1P15", "RPL31P62", "TCF7L1-IT1", "FMR1-IT1",
  "SCT", "DOK7", "CMPK2", "CKM", "BATF2", "LINC02988", "THEGL", "TEKT3", "PSG2", "RPL35P6",
  "RPSAP18", "MPIG6B", "KRT20", "MANEA-DT", "RET", "GJB5", "RNU1-72P", "GRHL1", "OAS2", "CPSF4L",
  "PPP1R14BP2", "CNTNAP2", "RNU6-1048P", "SLX1A", "NPM1P13", "RNA5SP435", "CKS1BP3", "HOXB-AS1",
  "DBH", "DHRS9", "SEMA6A-AS1", "CHUK-DT", "CHRM4", "KRT13", "GCA", "RAMP2", "AARD", "NA",
  "LINC02235", "C1orf226", "KLRD1", "BBOX1", "RXFP4", "PRAL", "SLC4A1", "HAL", "RSF1-IT2",
  "ST6GALNAC2", "KCNT1", "EXOC3L1", "MYCBP2-AS1", "RIMS4", "TMEM225B", "FUT1", "RN7SL3", "MIR4284",
  "BSND", "UBQLNL", "NA", "ANGPTL1", "BIRC7", "ANKRD31", "LCMT1-AS2", "DCST1", "PLIN1", "CDKN1C",
  "CXCL10", "ST6GALNAC1", "FASLG", "RPSAP70", "RAB38", "TEX29", "KRT40", "CPA2", "SLC6A4", "DPYS",
  "AK8", "MIR409", "RPL21P32", "TFR2", "GRM2", "KCNJ5"
)

# Generate synthetic dataset with real gene names
nephropathy_data <- data.table(
  Gene = sample(nephropathy_genes, 50, replace = TRUE),
  Fold_Enrichment = rnorm(50, mean = 2, sd = 0.5),
  Condition = "Nephropathy"
)

retinopathy_data <- data.table(
  Gene = sample(retinopathy_genes, 50, replace = TRUE),
  Fold_Enrichment = rnorm(50, mean = 3, sd = 0.7),
  Condition = "Retinopathy"
)

synthetic_data <- rbind(nephropathy_data, retinopathy_data)

# Filter common genes
common_genes <- intersect(nephropathy_genes, retinopathy_genes)
common_data <- synthetic_data[Gene %in% common_genes]

# Plot the data for common genes
ggplot(common_data, aes(x = Gene, y = Fold_Enrichment, fill = Condition)) +
  geom_bar(stat = "identity", position = "dodge", color = "black", width = 0.7) +
  labs(
    title = "Fold Enrichment for Common Genes in Nephropathy and Retinopathy",
    x = "Gene",
    y = "Fold Enrichment"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

