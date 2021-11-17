setwd("~/Desktop/Jiali/TAMU/689/week12/")
library(dplyr)

# read count files
filename <- list.files("count", pattern = "*txt", full.names = TRUE)
filename

file_list <- filename %>% lapply(read.csv, sep="\t", stringsAsFactors=F)

# combine protein IDs
protein.id <- do.call(rbind, file_list)[,1]
# remove duplicates
protein.id <- protein.id[!duplicated((protein.id))]
# remove contaminant 
protein.id <- protein.id[!(protein.id %in% protein.id[grep("contaminant_", protein.id)])]

protein.id

# Build data table
data <- data.frame(ID = protein.id)
for (i in 1:length(file_list)) {
  sample.id <- gsub("count/|.count.txt", "", filename[i])
  data <- merge(data, file_list[i], by.x="ID", by.y="protein.id", all.x=T)
  colnames(data)[i+1] <- sample.id
}

write.csv(data , "count/merged_count.csv", row.names = F)
# create a design table
colnames(data)
mydesign <- data.frame(label = colnames(data)[-1], condition = "", replicate = "")
write.csv(mydesign, "count/Data design.csv", row.names = F)
