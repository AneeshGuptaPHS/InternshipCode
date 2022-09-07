pacman::p_load("bugphyzz", "curatedMetagenomicData", "EnrichmentBrowser", "SummarizedExperiment", "dplyr", "stringr", "EnhancedVolcano", "dbplyr", "DESeq2", "png", "qusage", "taxize")
source("old_make_signatures.R")
source("deAnaPhyzz.R")

#Opens the data
ye <- curatedMetagenomicData("WirbelJ_2018.relative_abundance", dryrun = F, counts = T, rownames = "long")
ye.se <- ye[[1]]

#Makes the group which is needed to do deAna
grp1 = ifelse(ye.se$study_condition == "control" , 0, 1)
ye.se$GROUP = grp1

#Differential gene expression analysis
ye.se.Ana <- deAna(ye.se, de.method = "DESeq2", filter.by.expr = F)



reg_name <- function(MetaPhlAn){
  temp <- sub(".*g__", "", MetaPhlAn)
  gsub("\\|s__.*", "", temp)
  gsub("_noname", "", temp)
}

.getLast <- function(n)
{
  spl <- unlist(strsplit(n, "\\|"))
  spl[length(spl)]
}

#changes the rownames to the species name
rd <- rowData(ye.se.Ana)
rownames(rd) <- vapply(rownames(rd), reg_name,
                          character(1), USE.NAMES = FALSE)
rownames(rd) <- vapply(rownames(rd), .getLast,
                          character(1), USE.NAMES = FALSE)
rownames(rd) <- substr(rownames(rd), 4, nchar(rownames(rd)))

#filters results to only contain the bacteria of interest
sphingoProducers <- scan("new sps.txt", character(), quote = "")
rd <- as.data.frame(rd)
rd <- filter(rd, genus %in% sphingoProducers | species %in% sphingoProducers, abs(FC) != 30, FC != 0) #gets rid of the bacterial only present on one side
rd[ order((rownames(rd))), ] #clumps genera together and makes the table easier to use overall

#table of significant values
volcanoDataRed <- filter(rd, log2(abs(FC))>1, ADJ.PVAL < .05)
importantColumns <- subset(volcanoDataRed, select = c('species', 'FC', 'ADJ.PVAL'))

#Makes the volcano plot
png(filename = "volcanoWirbelJ.png", res = 300, height = 1800, width = 2500)
EnhancedVolcano(rd,
                lab = rd$genus,
                pCutoff = 0.05,
                title = 'Differential Abundance[DESeq2]',
                x = 'FC',
                y = 'ADJ.PVAL',
                xlim = c(-30, 45),
                ylim = c(0, 25),
                pointSize = 2.5)
dev.off()