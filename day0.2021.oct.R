library(data.table)
library(tidyverse)
library(bsseq)
library(DSS)
library(dplyr)
library(limma)

setwd("/project/6007512/C3G/projects/pubudu/covid19_new/covid_all")
dss_list <- list.files("DSS_input_wgs_all_corrected", pattern=paste0(".DSS.input.txt$"),recursive=T,full.names=T)

dss_list <- dss_list[dss_list  %like% "Day0" ]
dss_list <- dss_list[!dss_list  %like% "MED4" ]
  dss_list_df <-  lapply(dss_list, fread)

print("done loading")
  sample <- paste0(gsub("\\.DSS.input.txt*","",gsub(".*/","",dss_list)))
  sampledf <- as.data.frame(sample)
  sampledf$order <- 1:length(sample)
  BSobj = makeBSseqData( dss_list_df,
     sample )
  keepLoci = which(rowSums(getCoverage(BSobj)> 2 ) >=4)
  BSobj = BSobj[keepLoci,]
  cols <- fread("Phenotype_Day0_samples.txt")
  cols$Discharge <- factor(cols$Discharge, levels = c( "Alive", "Deceased"))

  cols <- merge(cols, sampledf, by.x="SampleID", by.y="sample")
  cols$SampleID <- factor(cols$SampleID, levels = sample)

  col2 <- cols

  col2$Sex_male <- ifelse(col2$Sex=="Male",1,0)
  col2 <- as.data.frame(col2)
	 rownames(col2) <- as.character(col2$SampleID)
   col2 <- select(col2, c("Discharge","Age_std", "Sex_male"))


print("start smoothing")
  BSobj.sm <- BSmooth(
    BSseq = BSobj, 
    BPPARAM = MulticoreParam(workers = 10), 
    verbose = TRUE)
	
 DMLfit = DMLfit.multiFactor(BSobj, design=col2, formula=~Discharge+Age_std+Sex_male ,smoothing=T)
 
 print("done multifactor")
 #testing Deceased vs Alive
 #change the direction in the final data frame
 test4 = DMLtest.multiFactor(DMLfit, coef="DischargeDeceased")
 dmrs4 = callDMR(test4, p.threshold=0.05)
 dmrs01 = callDMR(test4, p.threshold=0.01)
 
dmrs <- dmrs4

save.image(file="vini.DSS.Day0.smoothing.true.Deceased.vs.alive.DSS.remove.wgs.snp.newmed8.RData")
dmrs5 = callDMR(test4, p.threshold=0.05, minCG=2)

print("save image")
###################

dmlTest = DMLtest(BSobj, group1=rownames(subset(col2, Discharge=="Deceased")), group2=rownames(subset(col2, Discharge=="Alive")), smoothing=TRUE)

print("done no covariate")
dmls = callDML(dmlTest, p.threshold=1)
dmrs.nocov = callDMR(dmlTest, p.threshold=0.05)

write.table(dmrs.nocov, file = paste("mathieu.DSS.dmrs.p0.05.day0.smoothing.true.Deceased.vs.alive.remove.wgs.snp.nomed4.no.covariates.med8.txt",sep=""), sep = "\t", col.names=T, row.names=F,quote=F)
write.table(dmrs, file = paste("mathieu.DSS.dmrs.p0.05.day0.smoothing.true.Deceased.vs.alive.remove.wgs.snp..nomed4.with.covariates.med8.txt",sep=""), sep = "\t", col.names=T, row.names=F,quote=F)


###Actually these files have both data for with and without covariates
save.image(file="vini.DSS.Day0.smoothing.true.Deceased.vs.alive.DSS.remove.wgs.snp.newmed8.RData")


