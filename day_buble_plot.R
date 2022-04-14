library(GenomicRanges)
library(data.table)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(readr)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(ChIPseeker)
library(gtools)
library(org.Hs.eg.db)
require(biomaRt)
library(ggrepel)
library(eulerr)
library(ggpubr)
library(grid)
library(MAGeCKFlute)
library(clusterProfiler)
library(ReactomePA)
library(splitstackshape)
library(viridis)

file.name="DMRs.DSS.Day0.smoothing.true.Deceased.vs.alive.remove.wgs.snp.no.med4.with.covariates.deltabeta.annotation.med8.txt"

decease.cov <- fread("DMRs.vini.DSS.Day0.smoothing.true.Deceased.vs.alive.remove.wgs.snp.no.med4.with.covariates.deltabeta.med8.txt")

decease.day5_15.cov <- fread("DMRs.vini.DSS.Day5and15.smoothing.true.Deceased.vs.alive.remove.wgs.snp.no.med4.with.covariates.deltabeta.med8.txt")


day0 <- fread(file.name)

#pathway analysis
dff <- day0

#remove 2 dmrs associated with lot of HLA genes and PCDH gene cluster
#dff <- subset(dff, ID != "chr5_140710597_140710750" & ID != "chr6_32637800_32637989")

design="day0.deceased.vs.alive"
dff$anno_group <- dff$annotation


dff$anno_group[ dff$anno_group %like% "exon 1 of" ] <- "First exon"
dff$anno_group[ dff$anno_group %like% "Exon \\(" ] <- "Other exons"

dff$anno_group[ dff$anno_group %like% "intron 1 of" ] <- "First intron"
dff$anno_group[ dff$anno_group %like% "Intron \\(" ] <- "Other introns"

dff$anno_group[ dff$anno_group %like% "Downstream" ] <- "Downstream"

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
GeneHancer <- read_delim("GeneHancer_version_5.0.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

# filter out genes without a gene Symbol
GeneHancer$Genes_filtered <- GeneHancer$Genes
GeneHancer$Genes_filtered <- gsub("ENSG[^,]*,|LOC[^,]*,|GC[0-9][0-9][^,]*,|PIR[0-9][0-9][0-9][0-9][0-9][^,]*,","",GeneHancer$Genes_filtered)
GeneHancer$Genes_filtered[GeneHancer$Genes_filtered==""] <- NA
GeneHancer_GRange <- GRanges(GeneHancer)


# Select only 10kb promoter
DMR_all <- dff
DMR_all$DMR_ID <- paste(DMR_all$seqnames, DMR_all$start, DMR_all$end, sep="_")

# annotate genes in hg38 using GeneHancer dataset
#find all the enhancer genes overlapped with DMRs
DMRs_to_gene <- DMR_all
DMRs_to_gene <- dplyr::select(DMRs_to_gene, c("seqnames", "start", "end", "DMR_ID"))

rownames(DMRs_to_gene) <- 1:dim(DMRs_to_gene)[1]

DMRs_to_gene <- as.data.frame(mergeByOverlaps(query = GRanges(DMRs_to_gene), subject = GeneHancer_GRange))

DMRs_to_gene <- DMRs_to_gene[,c("DMR_ID","Genes","Genes_filtered")] # I kept the filtered list, but it canbe replaced for the non-filtered one if you desire
DMRs_to_gene <- aggregate(Genes_filtered ~ DMR_ID, DMRs_to_gene, function(gene)paste(unique(unlist(strsplit(gene, ","))), collapse = ',')) # merge duplicated DMRsIDs/Genes in a single row


Annotation <- merge(DMR_all, DMRs_to_gene, by = "DMR_ID", all=TRUE)
DMRs_to_gene <- Annotation[,c("seqnames","start","end","SYMBOL","Genes_filtered", "anno_group", "distanceToTSS","diff.Methy", "nCG", "areaStat")]
DMRs_to_gene[is.na(DMRs_to_gene)] <- ""
#Annotation_GRange <-  GRanges(Annotation_GRange)

#background gene list for enrichement
#all the genes in all annotations; i.e all genes and enhancer genes
bg.genelist <- unite(DMRs_to_gene[,4:5], col = "genes", sep=",", remove = TRUE, na.rm = TRUE)

bg.genelist <- as.vector(unique(unlist(strsplit(bg.genelist[,("genes")], ","))))

# merge the anotation For GeneHancer and promoter
promoter_DMRs <- subset(DMRs_to_gene, distanceToTSS > -5000 & distanceToTSS <= 1000 )


promoter_DMRs$genes <- unite(promoter_DMRs[,4:5], col = "genes", sep=",", remove = TRUE, na.rm = TRUE)

#promoter_DMRs$genes <- promoter_DMRs$genes$genes
#Annotation_temp$Targeted[Annotation_temp$Targeted==""]<-NA

# Add genes to the result files
rownames(promoter_DMRs) <- promoter_DMRs$DMR_ID
promoter_DMRs$Genes_filtered <- NULL 

perc_cutoff.promoter <- 90
perc_cutoff <- 90
promoter_DMRs <- promoter_DMRs[order(abs(promoter_DMRs$areaStat)),]

promoter_DMRs <- tail(promoter_DMRs, n=round(dim(promoter_DMRs)[1]*perc_cutoff.promoter/100))
ncg.cutoff=9
promoter_DMRs <- subset(promoter_DMRs, nCG > 3)
#select gene body

#select all the introns and exons
gebebody_DMRs <- subset(DMRs_to_gene, (anno_group %like% 'intron' | anno_group %like% 'exon' | anno_group %like% "UTR" & distanceToTSS > 1000))

gebebody_DMRs$genes <- unite(gebebody_DMRs[,4:5], col = "genes", sep=",", remove = TRUE, na.rm = TRUE)

#gebebody_DMRs$genes <- gebebody_DMRs$genes$genes
#Annotation_temp$Targeted[Annotation_temp$Targeted==""]<-NA

# Add genes to the result files
rownames(gebebody_DMRs) <- gebebody_DMRs$DMR_ID
gebebody_DMRs$Genes_filtered <- NULL 

gebebody_DMRs <- gebebody_DMRs[order(abs(gebebody_DMRs$areaStat)),]
gebebody_DMRs <- tail(gebebody_DMRs, n= round(dim(gebebody_DMRs)[1]*perc_cutoff/100))

gebebody_DMRs <- subset(gebebody_DMRs, nCG > ncg.cutoff)

DMRs_to_gene <- as.data.frame(DMRs_to_gene)

beta = 0
fdr <- 0.05

promoter_DMRs <- as.data.frame(promoter_DMRs)
gebebody_DMRs<- as.data.frame(gebebody_DMRs)



genehancer.hyper.promoter = unique(cSplit(promoter_DMRs, "genes", direction="long"))
genehancer.hyper.promoter.df <- subset(genehancer.hyper.promoter, areaStat > 0)


#gene body scrna
hyper.genebody = unique(cSplit(gebebody_DMRs, "genes", direction="long"))
hyper.genebody.df <- subset(hyper.genebody, areaStat > 0)


data.frame(genes=(as.vector(unique(unlist(strsplit(subset(promoter_DMRs, areaStat > beta  )[,("genes")], ",")))))) -> xx
write.table(xx, file="day0.hyper_promoter_genehancer.genes.txt", sep="\t", col.names=F, row.names=F, quote=F)

data.frame(genes=(as.vector(unique(unlist(strsplit(subset(promoter_DMRs, areaStat > beta  )[,("SYMBOL")], ",")))))) -> yy
write.table(yy, file="day0.hyper_promoter_only.genes.txt", sep="\t", col.names=F, row.names=F, quote=F)
#used https://bioinformatics.psb.ugent.be/cgi-bin/liste/Venn/calculate_venn.htpl to draw venn diagram using Vini's list and this gene list


Genes_GeneHancer <- list(
  "hyper_promoter"  = as.vector(unique(unlist(strsplit(subset(promoter_DMRs, areaStat > beta )[,("SYMBOL")], ",")))),
  "hypo_promoter"  = as.vector(unique(unlist(strsplit(subset(promoter_DMRs, areaStat < -beta )[,("SYMBOL")], ",")))),
  "hyper_promoter_genehancer"  = as.vector(unique(unlist(strsplit(subset(promoter_DMRs, areaStat > beta  )[,("genes")], ",")))),
  "hypo_promoter_genehancer"  = as.vector(unique(unlist(strsplit(subset(promoter_DMRs, areaStat < -beta )[,("genes")], ",")))),
  "hyper_genebody"  = as.vector(unique(unlist(strsplit(subset(gebebody_DMRs, areaStat >  beta )[,("SYMBOL")], ",")))),
  "hypo_genebody"  = as.vector(unique(unlist(strsplit(subset(gebebody_DMRs, areaStat <  -beta )[,("SYMBOL")], ",")))),
  "hyper_genebody_enhancer"  = as.vector(unique(unlist(strsplit(subset(gebebody_DMRs, diff.Methy >  beta )[,("genes")], ",")))),
  "hypo_genebody_enhancer"  = as.vector(unique(unlist(strsplit(subset(gebebody_DMRs, diff.Methy <  -beta )[,("genes")], ","))))
  
  
)

Genes_GeneHancer <- lapply(X = Genes_GeneHancer, FUN = function(x) {x <- TransGeneID(x, fromType = "Symbol", toType = "Entrez", organism = "hsa", ensemblHost = "www.ensembl.org")})


##############################################################################################################
###########                                       KEGG analyses                                    ###########
##############################################################################################################

# Kegg pathway analysis 
KEGG_GeneHancer <- lapply(X = Genes_GeneHancer, FUN = function(x) {x <- enrichKEGG(x, organism='hsa', pAdjustMethod="BH", qvalueCutoff=0.05, pvalueCutoff=0.05, minGSSize = 10)})
KEGG_GeneHancer <- lapply(X = KEGG_GeneHancer, FUN = function(x) {x <- setReadable(x = x, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")})


##############################################################################################################
###########                                  Reactome  analyses                                    ###########
##############################################################################################################

Reactome_GeneHancer <- lapply(X = Genes_GeneHancer, FUN = function(x) {x <- enrichPathway(x, organism='human', pAdjustMethod="BH", qvalueCutoff=0.05, pvalueCutoff=0.05, minGSSize = 10, readable = TRUE)})

##############################################################################################################
###########                                   GO BP analyses                                       ###########
##############################################################################################################

GO_BP_GeneHancer <- lapply(X = Genes_GeneHancer, FUN = function(x) {x <- enrichGO(x, OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod="BH", qvalueCutoff = 0.05, pvalueCutoff = 0.05, minGSSize = 10, readable = TRUE)})

#############################################################################################
#################             Manhattan plot for all tested databases      ##################
#############################################################################################

# GeneHancer
#all peaks
go_list <- list()

for (i in 1:length(Genes_GeneHancer)) {
  
  Manhattan_KEGG <- KEGG_GeneHancer[i][[1]]@result
  Manhattan_Reactome <- Reactome_GeneHancer[i][[1]]@result
  Manhattan_GO <- GO_BP_GeneHancer[i][[1]]@result
  
  # Plot
  Manhattan_KEGG$Approach <- "KEGG"
  Manhattan_Reactome$Approach <- "Reactome"
  Manhattan_GO$Approach <- "GO"
  
  Manhattan_data <- rbind(Manhattan_GO, Manhattan_Reactome, Manhattan_KEGG)
  Manhattan_data <- Manhattan_data[order(Manhattan_data$ID),]
  Manhattan_data <- subset(Manhattan_data, Count >= 5)
  Manhattan_data$BH_combined = p.adjust(Manhattan_data$pvalue, method = "BH")
  Manhattan_data <- Manhattan_data[order(Manhattan_data$BH_combined),]
  
  # Manhattan plot 
  #check whether there are more significant pathways. if there are no significant pathways. there will be no lables
  if(dim(Manhattan_data)[1]>0){
    Selected_terms <- c(if (dim(subset(Manhattan_data, BH_combined < fdr & Approach == "GO"))[1] > 10) subset(Manhattan_data, Approach == "GO")[1:10,"ID"] else subset(Manhattan_data, Approach == "GO")[1:dim(subset(Manhattan_data, BH_combined < fdr & Approach == "GO"))[1],"ID"], 
                        if (dim(subset(Manhattan_data, BH_combined < fdr & Approach == "KEGG"))[1] > 10) subset(Manhattan_data, Approach == "KEGG")[1:5,"ID"] else subset(Manhattan_data, Approach == "KEGG")[1:dim(subset(Manhattan_data, BH_combined < fdr & Approach == "KEGG"))[1],"ID"],
                        if (dim(subset(Manhattan_data, BH_combined < fdr & Approach == "Reactome"))[1] > 10) subset(Manhattan_data, Approach == "Reactome")[1:5,"ID"] else subset(Manhattan_data, Approach == "Reactome")[1:dim(subset(Manhattan_data, BH_combined < fdr & Approach == "Reactome"))[1],"ID"]
    )
    
    
    p <- ggplot(Manhattan_data, aes(x = ID, y = -log10(BH_combined), color = Approach, label = Description)) +
      geom_point(aes(size = Count), alpha = 0.5) +
      ylab("-log10 (FDR)") +
      scale_size(range = c(1,10)) +
      geom_text_repel(data = subset(Manhattan_data, ID %in% Selected_terms), aes(label = Description),
                      segment.color = "grey75", color = "grey25", nudge_y = 2.5, angle = 0, size = 4,
                      show.legend = FALSE, force = 10) +
      geom_hline(yintercept = 1.3, color="black", linetype="dotted", size = 1) +
      ylim(0,25) +
      theme_bw() +
      ggtitle(paste0(names(Genes_GeneHancer)[i]))+
      theme(axis.title = element_text(size = 14, face = "bold"), axis.text.y = element_text(size = 12),
            axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.grid = element_blank())
    print(p)
    print("hello")
    
    Manhattan_data$group <- paste0(names(Genes_GeneHancer)[i])
    go_list[[i]] <- Manhattan_data
  }
  
}





bind_rows(go_list) -> go_df
go_df_methy <- go_df

#uncomment if u want to save the table
#write.table(go_df_methy, file="methylation.enrichment_5k_promoter.txt", sep="\t", col.names=T, row.names=F, quote=F)
#go_df_methy <- subset(go_df_methy, BH_combined < 0.05 | (BH_combined < 0.1 & group=="hyper_promoter_genehancer" | group=="hypo_promoter_genehancer")) 
go_df_methy <- subset(go_df_methy, BH_combined < 0.05) 

go_df_methy_long <- unique(cSplit(go_df_methy, "geneID", direction="long", sep = "/"))

genehancer.hypo.promoter.df <- subset(genehancer.hyper.promoter, areaStat < 0)
hypo_genehancer_interest <- merge(genehancer.hypo.promoter.df, go_df_methy_long, by.x="genes", by.y="geneID")
hyper_genehancer_interest <- merge(genehancer.hyper.promoter.df, go_df_methy_long, by.x="genes", by.y="geneID")
write.table(hypo_genehancer_interest, file="Table_S5-methylation-d0_pathways_after_review.txt", sep="\t", col.names=T, row.names=F, quote=F)

go_df_methy$hla <- str_count(go_df_methy$geneID, "HLA")
go_df_methy$pcdh <- str_count(go_df_methy$geneID, "PCDH")

go_df_methy_ori <- go_df_methy

#go_df_methy$hla_count <- go_df_methy$Count - go_df_methy$hla
#go_df_methy$pcdh_count <- go_df_methy$Count - go_df_methy$pcdh
#go_df_methy <- subset(go_df_methy, pcdh_count > 5 & hla_count > 5 )
go_df_methy <- subset(go_df_methy, (Count-hla > 5) & (Count-pcdh > 5) )



add <- c("GO:0002504", "R-HSA-2132295", "hsa05322", "hsa04658", "GO:0060333", "R-HSA-877300")

go_df_methy_add <- subset(go_df_methy_ori, ID %in% add)
remove <- c("GO:0071230", "hsa00983")

go_df_methy <- subset(go_df_methy, ! ID %in% remove)

go_df_methy <- rbind(go_df_methy, go_df_methy_add)

add_additional_regions <- function(df, group){
  
  
  df$GeneRatio <- "0/1"
  df$BH_combined <- 1
  df$group <- group
  return(df)
}
#uncomment when run aday0
#write.table(go_df_methy, file="day0.methylation.enrichment_5k_promoter_selected.pathways.txt", sep="\t", col.names=T, row.names=F, quote=F)

###map pathways to DMRs
#temp.g.promo.hyper <- subset(go_df_methy, )
#write.table(go_df_methy, file="hyper_genehancer_promoter.txt", sep="\t", col.names=T, row.names=F, quote=F)
#write.table(go_df_methy, file="methylation.enrichment.txt", sep="\t", col.names=T, row.names=F, quote=F)
go_df_methy[order(go_df_methy$BH_combined),] -> go_df_methy

#"R-HSA-2132295", "GO:0002504"- all
#"GO:0060333", "hsa05322", "hsa04658" , "R-HSA-877300" - hypermethy

manually_added_pathways <- c("R-HSA-2132295", "GO:0002504", "GO:0060333", "hsa05322", "hsa04658" , "R-HSA-877300")

temp_manual_df <- subset(go_df_methy, ID %in% manually_added_pathways)

temp_manual_df <- temp_manual_df[c(1,2,4,6,7,9,12,13),]

go_df_methy <- subset(go_df_methy, !ID %in% manually_added_pathways)


#go_df_methy <- go_df_methy[!duplicated(go_df_methy$ID), ]
#add two genebody pathways to represent
add_pathways <- c("GO:0030198", "GO:0000819")
add_pathways <- subset(go_df_methy, ID %in% add_pathways)

go_df_methy_save <- rbind(go_df_methy, temp_manual_df)
write.table(go_df_methy_save, file="Table_S5-methylation-d0_pathways_after_review_2022_01_28.txt", sep="\t", col.names=T, row.names=F, quote=F)
go_df_methy[order(go_df_methy$BH_combined),] %>% head(n=44) -> go_df_methy


go_df_methy <- rbind(go_df_methy, temp_manual_df)


#uset this to save data frame for manuscript
##write.table(go_df_methy_save, file="day0.methylation.enrichment_5k_promoter_selected.pathways.txt", sep="\t", col.names=T, row.names=F, quote=F)
go_df_methy %>% separate(GeneRatio, c("Count","Total"), "/") -> go_df_methy
go_df_methy %>% separate(BgRatio, c("bgCOunt","BGTotal"), "/") -> go_df_methy


add_pathways %>% separate(GeneRatio, c("Count","Total"), "/") -> add_pathways
add_pathways %>% separate(BgRatio, c("bgCOunt","BGTotal"), "/") -> add_pathways

go_df_methy$Percentage <- as.integer(go_df_methy$Count)/ as.integer(go_df_methy$bgCOunt)

add_pathways$Percentage <- as.integer(add_pathways$Count)/ as.integer(add_pathways$bgCOunt)


removable_pathways = c("hsa00982", "GO:0034470", "GO:0060333", "GO:0048704", "GO:0002757", "GO:0048002", "GO:0002478" ,"GO:0019884", "GO:0019882", "R-HSA-5619507", "GO:0050851", "hsa00980", "hsa05164")
go_df_methy <- subset(go_df_methy, !ID %in% removable_pathways)

go_df_methy_temp1 <- subset(go_df_methy, Approach=="GO")
go_df_methy_temp1[order(go_df_methy_temp1$BH_combined),] %>% head(n=8) -> go_df_methy_temp1
go_df_methy_temp2 <- subset(go_df_methy, Approach!="GO")

go_df_methy <- rbind(go_df_methy_temp1, go_df_methy_temp2, add_pathways)

#######create data frme with showing all pathways even without they are not present in a group
reshape2::dcast(unique(go_df_methy), ID  ~ group , value.var="BH_combined", fill=1) -> plot.df.dc1

reshape2::melt(plot.df.dc1,  id.vars = c("ID")) -> plot.df.melt1

colnames(plot.df.melt1) <- c("ID", "group","BH_combined")



plot.df.temp <- merge(unique(dplyr::select(go_df_methy, c("ID", "group", "Percentage"))), plot.df.melt1, by=c("ID", "group") ,all.y=T)

plot.df.temp <- merge(plot.df.temp, unique(dplyr::select(go_df_methy, c("ID", "Approach"))), by="ID")
plot.df.temp <- merge(plot.df.temp, unique(dplyr::select(go_df_methy, c("ID", "Description"))), by="ID")

plot.df.temp$DOC <- ifelse(plot.df.temp$group %like% "hyper", "hypermethylated", ifelse(plot.df.temp$group %like% "hypo", "hypomethylated", "motif"))

# plot.df.temp$group <- factor(plot.df.temp$group, levels=c("hyper_promoter_only","hypo_promoter_only","hyper_promoter_and_genhancer","hypo_promoter_and_genhancer",
#                                                           "hyper_genebody_only", "hypo_genebody_only", "hyper_genebody_and_genhancer",
#                                                           "hypo_genebody_and_genhancer"))
plot.df.temp$group <- factor(plot.df.temp$group, levels=unique(plot.df.temp$group))

plot.df.temp[is.na(plot.df.temp)] <- 0


bubble_plot_methy <- function(approach=NULL){
  ggplot(subset(unique(plot.df.temp), Approach==approach), aes(y = group, x = Description)) + 
    geom_point(aes(fill = -log10(BH_combined), size = Percentage), alpha = 0.75, shape = 21) +
    #scale_fill_gradient2(low = "orange2", mid ="grey85" , high = "skyblue2", limits = c(1.3, 5), midpoint = 3) +
    scale_fill_viridis(option = "A", limits = c(0, 3), name = "-log10\n(FDR)", direction = -1, oob = scales::squish) +
    scale_size_continuous( breaks = c(0,0.05,0.1,0.15,0.2)) +
    ggtitle(approach)+
    facet_grid(rows = vars(DOC), scales = "free", space='free') +
    # scale_fill_gradientn(colours = grDevices::rainbow(n = 10)) +
    
    theme( axis.text.x = element_text(colour = "black", size = 17, face = "bold", angle = -45, vjust = 1.2, hjust = 0), 
           axis.text.y = element_text(colour = "black", face = "bold", size = 16),
           axis.title = element_blank(),
           legend.text = element_text(size = 14, face ="bold", colour ="black"), 
           legend.title = element_text(size = 14, face = "bold"), panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2), 
           legend.position = "right")
}


P1 <- bubble_plot_methy("KEGG")
P2 <- bubble_plot_methy("Reactome")
P3 <- bubble_plot_methy("GO")

genebar.df <- as.data.frame(lengths(Genes_GeneHancer))
#genebar.df <- subset(genebar.df, DOC!="motif")
setDT(genebar.df, keep.rownames = TRUE)[]
colnames(genebar.df) <- c("group", "count")
genebar.df <- subset(genebar.df, group=="hyper_promoter_genehancer" | group=="hypo_promoter_genehancer" | group =="hypo_genebody_enhancer" | group=="hyper_genebody")
genebar.df$group <-  factor(genebar.df$group, level= levels(plot.df.temp$group))



genebar.df$DOC <- ifelse(genebar.df$group %like% "hyper", "hypermethylated", ifelse(genebar.df$group %like% "hypo", "hypomethylated", "motif"))
p4 <- ggplot(genebar.df, aes(x=group, y=count)) +
  geom_bar(stat='identity', width = 0.5) +
  # scale_y_continuous( breaks= seq(0,3000,250), limit=c(0,3000))+
  # theme_bw() +
  facet_grid(rows = vars(DOC), scales = "free", space='free') +
  theme(axis.title = element_text(size = 24, face = "bold"), axis.text.y = element_text(size = 18),
        axis.line= element_blank(), axis.line.x = element_line(colour = "black"),  panel.grid = element_blank(), panel.background = element_blank())+
  coord_flip()

pdf(file=paste0(design,".DSS.DMRs.pathways.at.deltabeta.",beta,".wgs.snp_bublleplot.new.final2.bar.pdf"),height=10,width=15)
ggarrange(P3 + theme(strip.text.y = element_blank()),
          P2 + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank(), strip.text.y = element_blank()),
          P1 + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank()),
          p4+ theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank(), strip.text.y = element_blank()),
          nrow = 1, align = "h", common.legend = TRUE, legend = "right" , widths = c(2,0.8,1,0.2))
#nrow = 1, align = "h", common.legend = TRUE, legend = "right" , widths = c(2.2,1.2,3.7))

# ggarrange(
#   P2 + theme(axis.text.y = element_blank()),
#   P1 + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank(), strip.text.y = element_blank()),
#   p4+ theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank()),
#   nrow = 1, align = "h", common.legend = TRUE, legend = "right" , widths = c(1.2,1.2,0.5))

##4th width value was 0.3

dev.off()


##4th width value was 0.3

dev.off()

