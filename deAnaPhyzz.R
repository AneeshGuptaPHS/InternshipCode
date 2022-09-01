#' 
#' @param study a string, the name of a cmd dataset e.g 
#' @param gmt a gmt 
#' @param condition a string - the study condition e.g. cirrhosis 
#' 
deAnaPhyzz <- function(study, condition){
  study  <- paste(study, "relative_abundance", sep = ".")
  dat.se <- curatedMetagenomicData(study, dryrun = F, counts = T, rownames = "NCBI")[[1]]
  
  grp1 = ifelse(dat.se$study_condition == condition , 1, 0)
  dat.se$GROUP = grp1
  assay(dat.se) <- assay(dat.se) + 1
  
  dat.se.Ana <- deAna(dat.se, de.method = "DESeq2", filter.by.expr = F)
  return(dat.se.Ana)
}

#' 
#' @param sigs a bugphyzz signature
#' 
#' 
makeSpeciesPhyzzGMT <- function(sigs){
  children <- sigs %>% 
    lapply(taxizedb::children, db = "ncbi") %>% 
    lapply(data.table::rbindlist) %>%
    lapply(dplyr::pull,id)
  sigs <- list(children, sigs) %>% 
    unlist(recursive=F)
  sigs <- tapply(unlist(sigs, use.names = FALSE), rep(names(sigs), lengths(sigs)), FUN = c)
  EnrichmentBrowser::writeGMT(sigs, gmt.file = "my.gmt")
  gmt <- qusage::read.gmt("my.gmt")
  return(gmt)
}

#' 
#' @param dat.se.Ana 
#' 
oraGsea <- function(dat.se.Ana, ora_perm = 0, gsea_perm = 1000){
  # Overrepresentation analysis (ORA)
  ora_results <- sbea("ora", dat.se.Ana, gmt, perm = ora_perm)
  ora_results <- gsRanking(ora_results, signif.only = F) %>%
    as.data.frame()
  colnames(ora_results) <- sub("GENE", "MICROBE", colnames(ora_results))
  
  # Gene set enrichment analysis (GSEA)
  gsea_results <- sbea("gsea", dat.se.Ana, gmt, perm = gsea_perm)
  gsea_results <- gsRanking(gsea_results, signif.only = F) %>%
    as.data.frame()
  colnames(gsea_results) <- sub("GENE", "MICROBE", colnames(gsea_results))
  
  results <- list(ora_results, gsea_results)
  names(results) <- c("ora_results", "gsea_results")
  
  return(results)
}
