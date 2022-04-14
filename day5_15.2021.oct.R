library(data.table)
library(tidyverse)
library(bsseq)
library(DSS)
library(dplyr)
library(limma)

setwd("/project/6007512/C3G/projects/pubudu/covid19_new/covid_all")
dss_list <- list.files("DSS_input_wgs_all_corrected", pattern=paste0(".DSS.input.txt$"),recursive=T,full.names=T)



dss_list <- dss_list[!dss_list  %like% "Day0" ]
dss_list <- dss_list[!dss_list  %like% "MED4" ]
  dss_list_df <-  lapply(dss_list, fread)

  sample <- paste0(gsub("\\.DSS.input.txt*","",gsub(".*/","",dss_list)))
  sampledf <- as.data.frame(sample)
  sampledf$order <- 1:length(sample)
  BSobj = makeBSseqData( dss_list_df,
     sample )
  keepLoci = which(rowSums(getCoverage(BSobj)> 2 ) >=4)
  BSobj = BSobj[keepLoci,]

    cols <- fread("Phenotype_all_samples.txt")
  cols$Discharge <- factor(cols$Discharge, levels = c( "Alive", "Deceased"))

  cols <- merge(cols, sampledf, by.x="SampleID", by.y="sample")
  cols$SampleID <- factor(cols$SampleID, levels = sample)

  col2 <- cols
  print("start smoothing")
  BSobj.sm <- BSmooth(
    BSseq = BSobj, 
    BPPARAM = MulticoreParam(workers = 10), 
    verbose = TRUE)
    print("done smoothing")
  col2$Sex_male <- ifelse(col2$Sex=="Male",1,0)
  col2 <- as.data.frame(col2)
	 rownames(col2) <- as.character(col2$SampleID)
  #col2 <- select(col2, c("Discharge","Age_std", "Sex_male","Ordinal_Scores_std"))
  col2 <- select(col2, c("Discharge","Age_std", "Sex_male"))
 
 
 
 
 DMLfit = DMLfit.multiFactor(BSobj, design=col2, formula=~Discharge+Age_std+Sex_male ,smoothing=T)
   print("multifactor done")
 #testing Deceased vs Alive
 #change the direction in the final data frame
 test4 = DMLtest.multiFactor(DMLfit, coef="DischargeDeceased")
 dmrs4 = callDMR(test4, p.threshold=0.05)

dmrs <- dmrs4


save.image(file="mathieu.DSS.Day5and15.smoothing.true.Deceased.vs.alive.DSS.remove.wgs.snps.med8.RData")

  print("save image ")

###################

dmlTest = DMLtest(BSobj, group1=rownames(subset(col2, Discharge=="Deceased")), group2=rownames(subset(col2, Discharge=="Alive")), smoothing=TRUE)

  print("no multifactor start")
dmls = callDML(dmlTest, p.threshold=1)
dmrs.nocov = callDMR(dmlTest, p.threshold=0.05)

  print("done no multifactor")
write.table(dmrs.nocov, file = paste("mathieu.DSS.dmrs.p0.05.day5and15.smoothing.true.Deceased.vs.alive.remove.wgs.snps.nomed4.no.covariates.med8.txt",sep=""), sep = "\t", col.names=T, row.names=F,quote=F)
write.table(dmrs, file = paste("mathieu.DSS.dmrs.p0.05.day5and15.smoothing.true.Deceased.vs.alive.remove.wgs.snps.nomed4.with.covariates.med8.txt",sep=""), sep = "\t", col.names=T, row.names=F,quote=F)

  print("write done")
###Actually these files have both data for with and without covariates
save.image(file="mathieu.DSS.Day5and15.smoothing.true.Deceased.vs.alive.DSS.remove.wgs.snps.med8.RData")


