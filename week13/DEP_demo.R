setwd("~/Desktop/Jiali/TAMU/689/week13")

library("DEP")
library(dplyr)

mydata <- read.csv("../week12/count/merged_count.csv", header = T)
mydata$ID %>% duplicated() %>% any()

mydataDesign <- read.csv("../week12/count/Data design.csv", header = T)

# generate name column
mydata$name <- gsub("\\|.*","", mydata$ID)

mydataDesign$label <- colnames(mydata)[2:9]

# generate DEP object
data_se <- make_se(mydata, c(2:9), mydataDesign)

# filter out low expression proteins
data_filt <- filter_proteins(data_se, "fraction", min = 0.25)
plot_numbers(data_filt, plot = T)

# impute missing 
protein_MNAR <- get_df_long(data_filt) %>% 
  group_by(name, condition) %>% 
  dplyr::summarize(NAs = all(is.na(intensity))) %>%
  filter(NAs) %>% 
  pull(name) %>%
  unique()

MNAR <- names(data_filt) %in% protein_MNAR

data_impute <- DEP::impute(
  data_filt,
  fun = "mixed",
  randna = !MNAR,
  mar = "knn",
  mnar = "MinProb"
)

df <- get_df_wide(data_impute)
write.csv(df, "imputed_counts.csv", row.names = F)

# data normalization
data_norm <- normalize_vsn(data_impute)
meanSdPlot(data_norm)
plot_imputation(data_norm, data_impute)

# DE test
data_diff <- test_diff(data_norm, type= "all")

# Plot PCA 
plot_pca(data_diff, x = 1, y = 2, n = 4876, plot = TRUE)

# extract DEP
def <- add_rejections(data_diff, alpha = 0.1, lfc = log2(1))

data_results <- get_results(def)

data_results %>% filter(significant) %>% nrow()

# output results
CNFvCL <- data_results[data_results$CNF_vs_CL_p.adj < 0.1, c(1,12,25)]

# Read annotations
GO <- read.csv("Irplac1_GeneCatalog_proteins_20170103_GO.tab", sep = "\t")
IPR <- read.csv("Irplac1_GeneCatalog_proteins_20170103_IPR.tab", sep = "\t")
KEGG <- read.csv("Irplac1_GeneCatalog_proteins_20170103_KEGG.tab", sep = "\t")
KOG <- read.csv("Irplac1_GeneCatalog_proteins_20170103_KOG.tab", sep = "\t")

CNFvCL_annot <- merge(CNFvCL, GO[,c(1,3)], by.x = "name", by.y = "X.proteinId", all.x = T)
CNFvCL_annot <- merge(CNFvCL_annot, IPR[,c(1,3)], by.x = "name", by.y = "X.proteinId", all.x = T)
CNFvCL_annot <- merge(CNFvCL_annot, KEGG[,c(1,3)], by.x = "name",by.y="X.proteinId", all.x = T)
CNFvCL_annot <- merge(CNFvCL_annot, KOG[,c(2,4)], by.x = "name", by.y ="proteinId", all.x = T)

# remove duplicated lines
uniq <- unique(CNFvCL_annot)

# save the DE protein table into csv file
write.csv(uniq, "CNFvCL_DE.csv", row.names = F)
