pacman::p_load("bugphyzz", "curatedMetagenomicData", "EnrichmentBrowser", "SummarizedExperiment", "dplyr", "stringr", "EnhancedVolcano", "dbplyr", "DESeq2", "png", "qusage", "taxize")
source("old_make_signatures.R")
source("deAnaPhyzz.R")

ye <- curatedMetagenomicData("WirbelJ_2018.relative_abundance", dryrun = F, counts = T, rownames = "long")
ye.se <- ye[[1]]
grp1 = ifelse(ye.se$study_condition == "control" , 0, 1)
ye.se$GROUP = grp1
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


rd <- rowData(ye.se.Ana)
rownames(rd) <- vapply(rownames(rd), reg_name,
                          character(1), USE.NAMES = FALSE)
rownames(rd) <- vapply(rownames(rd), .getLast,
                          character(1), USE.NAMES = FALSE)
rownames(rd) <- substr(rownames(rd), 4, nchar(rownames(rd)))


sphingoProducers <- scan("new sps.txt", character(), quote = "")
rd <- as.data.frame(rd)
rd <- filter(rd, genus %in% sphingoProducers, abs(FC) != 30, FC != 0)
rd[ order((rownames(rd))), ] 

#test <- filter()

#rownames(volcanoData) <- NULL
volcanoDataRed <- filter(rd, log2(abs(FC))>1, ADJ.PVAL < .05)

png(filename = "volcanoWirbelJ.png", res = 300, height = 1800, width = 2500)
EnhancedVolcano(rd,
                lab = rownames(rd),
                pCutoff = 0.05,
                title = 'Differential Abundance[DESeq2]',
                x = 'FC',
                y = 'ADJ.PVAL',
                xlim = c(-30, 45),
                ylim = c(0, 25),
                pointSize = 2.5)
dev.off()