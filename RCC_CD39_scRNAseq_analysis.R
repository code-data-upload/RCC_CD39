 


library(tidyverse); library(dplyr); library(readxl); library(openxlsx); library(ggpubr);
library(gridExtra); library(factoextra); library(circlize); library(seriation); library(dendextend);
library(tibble); library(ggrepel); library(devtools); library(pheatmap); library(xtable);
library(ggsci); library(ggplot2); library(viridis); library(immunarch); library(reshape);
library(psych); library(Seurat); library(patchwork); library(glmGamPoi); library(cowplot);
library(DeconRNASeq); library(clusterProfiler); library(hopach); library(mygene); library(fgsea);
library(R.utils); library(foreach); library(doParallel); library(stringr); library(MCL);
library(scran); library(scDblFinder); library(RColorBrewer); library(harmony); library(msigdbr);
library(ggbreak); library(Nebulosa); library(unikn); library(BiocFileCache); library(Matrix);
library(DDRTree); library(monocle);

rna_counts_matrix<-read.table( "GSE285701_single_cell_rna_data.txt", sep="\t")
adt_counts_matrix<-read.table( "GSE285701_single_cell_adt_data.txt", sep="\t")


### RNA unbiased clustering
#######################################################################################################

sr<-CreateSeuratObject(counts = rna_counts_matrix, names.delim = "_")
DefaultAssay(sr) <- 'RNA'
sr = NormalizeData(sr, normalization.method = "LogNormalize", scale.factor = 10000) ;
sr = FindVariableFeatures(sr, selection.method="vst", nfeatures=2000) ;
sr = ScaleData(sr, features=rownames(sr)) 

sr = RunPCA(sr, features=VariableFeatures(object=sr), npcs=100) ;
sr = RunHarmony(sr, "orig.ident")

ElbowPlot(sr,50)
sr = RunUMAP(sr, reduction="harmony", dims=1:20)
sr = FindNeighbors(sr, reduction="harmony", dims=1:20)

sr = FindClusters(sr, resolution=0.21)
sr$cell_type<-recode(sr$seurat_clusters,"0"="Tex", "1"="Tem", "2"="Temra", "3"="Tcyto", "4"="Tp")

DimPlot(sr, label.size = 7, label = T) + coord_equal()



### classification of ADT CD39+/-PD1+/- cells
#######################################################################################################

sr[["ADT"]] <- CreateAssayObject(counts = adt_counts_matrix)
DefaultAssay(sr) <- 'ADT'

VariableFeatures(sr) <- rownames(sr[["ADT"]])
sr <- NormalizeData(sr, normalization.method = 'CLR', margin = 2) %>%  ScaleData() %>% RunPCA(reduction.name = "apca", approx =FALSE)

sr$anchor<-"anchor"
par(mfrow=c(1,2))
for( tmp_gene in c( "PD-1", "CD39")){
  Da = density(sr@assays$ADT$data[tmp_gene,])
  DeltaY = diff(Da$y)
  Turns = which(DeltaY[-1] * DeltaY[-length(DeltaY)] < 0) +1

    print("########################################");     print(tmp_gene);     print(Da$x[Turns[1:4]])
    plot(Da, xlab="", ylab="", main=tmp_gene)+points(Da$x[Turns[1:4]], Da$y[Turns[1:4]], pch=c(1,16,1,1), col=c("black", "red", "black", "black"))
    text(Da$x[Turns[2]], Da$y[Turns[2]] + 0.03, labels = round(Da$x[Turns[2]], 7), col = "red", cex = 1)
}






tcr_info<-read.table("GSE285701_single_cell_tcr_info.txt")
sr$cell_barcode<-rownames(sr@meta.data)
sr@meta.data<-left_join(sr@meta.data, tcr_info, by = "cell_barcode")

### HLA types and TCR sequences mapped to virus sequences
#######################################################################################################

sr$mhc<-recode(sr$orig.ident,
               "RB183"="HLA-A_11:01;HLA-A_24:02;HLA-B_51:01;HLA-B_54:01;HLA-C_01:02;HLA-C_14:02",
               "RT183"="HLA-A_11:01;HLA-A_24:02;HLA-B_51:01;HLA-B_54:01;HLA-C_01:02;HLA-C_14:02",
               "RB208"= "HLA-A_02:01;HLA-A_24:02;HLA-B_07:02;HLA-B_40:06;HLA-C_01:02;HLA-C_07:67",
               "RT208"="HLA-A_02:01;HLA-A_24:02;HLA-B_07:02;HLA-B_40:06;HLA-C_01:02;HLA-C_07:67",
               "RT199"="HLA-A_02:07;HLA-A_33:03;HLA-B_44:03;HLA-B_46:01;HLA-C_01:03;HLA-C_14:03")
              
hla<-unlist(strsplit(unique(sr$mhc), split = ";"))
tcr_meta<-sr@meta.data[!is.na(sr@meta.data$Sequence),]
tcr_meta$cb<-rownames(tcr_meta)


vdjdb = dbLoad("https://gitlab.com/immunomind/immunarch/raw/dev-0.5.0/private/vdjdb.slim.txt.gz", "vdjdb")
vdjdb<-vdjdb%>%  filter(gene=="TRB", Pathology %in% c("CMV", "EBV", "InfluenzaA"), mhc.class=="MHCI")
vdjdb.TRB$mhc.a<-gsub("[*]","_", vdjdb.TRB$mhc.a)



### DB extraction of patient HLA-types

mhc_matched_db<-c()
for(hh in  hla){
  tmp_mhc_matched_db <- subset(vdjdb.TRB, str_detect(vdjdb.TRB$mhc.a,  hh) )
  if( dim(tmp_mhc_matched_db)[1]!=0 ){
    tmp_mhc_matched_db$from_mhc<-hh
    mhc_matched_db<-rbind(tmp_mhc_matched_db, mhc_matched_db)
  }
}


### cdr3aa matching from our cell cdr3aa to extracted db

final_matched_meta<-c()
for( i in unique(mhc_matched_db$from_mhc)){ 
  one_mhc_matched_db<-subset(mhc_matched_db, from_mhc== i)
  one_mhc_matched_meta<-subset(tcr_meta, str_detect(string=tcr_meta$mhc, i))
  
  for( j in one_mhc_matched_db$cdr3){ 
    tmp_matched_meta<-subset(one_mhc_matched_meta, str_detect(string=one_mhc_matched_meta$CDR3.aa, j))
    
    if( dim(tmp_matched_meta)[1]!=0 ){
      tmp_matched_meta$from_mhc<-i;  tmp_matched_meta$from_cdr3<-j
      tmp_matched_meta$pathology<-paste0(unique(   subset(one_mhc_matched_db, cdr3==j)$Pathology   ),collapse = "_")
      
      if( length(intersect(rownames(final_matched_meta), rownames(tmp_matched_meta)))==0 ){ 
        final_matched_meta<-rbind(tmp_matched_meta, final_matched_meta) }
    }
  }
}





## tumor-specific clonotypes

#######################################################################################################

# identification of tumor-specific cells (scRNAseq) by mapping tumor-reactive clonotypes (bulk seq)

TS_data <- vector()
for (i in 1:length(sr_meta$CDR3.nt)) {
  if ((sum(str_detect(sr_meta$CDR3.nt[i], TS_RT183$data$RT183.clonotypes.TRB$CDR3.nt)) >0) |
      (sum(str_detect(sr_meta$CDR3.nt[i], TS_RT199$data$RT199.clonotypes.TRB$CDR3.nt)) >0) | 
      (sum(str_detect(sr_meta$CDR3.nt[i], TS_RT208$data$RT208.clonotypes.TRB$CDR3.nt)) >0)) {TS_data[i] <- 1
  }  else {TS_data[i] <- 0 }
}

sr_meta$TS <- TS_data
sr$TS<-mapvalues(rownames(sr@meta.data), from = rownames(sr_meta), to=sr_meta$TS)
sr$TS[sr$TS!="0"&sr$TS!="1"]<-NA



# Tumor-specific clonotype shared between blood and tumor

sr_tcr_meta<-subset(sr@meta.data, !is.na(Sequence))
clns_table<-data.frame(table(sr_tcr_meta$Sequence))
clns_table<-clns_table[rev(order(clns_table$Freq)),]

clns_table$TS<-mapvalues(as.character(clns_table$Var1), from = sr_tcr_meta$Sequence, to=sr_tcr_meta$TS)
clns_table_ts1<-subset(clns_table, TS==1); clns_table_ts1

tcr_ts1_rtrb_meta<-sr_tcr_meta%>%filter(Sequence %in% clns_table_ts1$Var1)
tcr_ts1_info<-data.frame(as.data.frame.matrix(table(tcr_ts1_rtrb_meta$Sequence, tcr_ts1_rtrb_meta$sample_type)))
tcr_ts1_info$Var1<-rownames(tcr_ts1_info)
tcr_ts1_info<-merge(x=clns_table_ts1, y=tcr_ts1_info)

# RT vs. RB DEG of tumor-specific clonotype #4, highly represented in both blood and tumor.

rep_clone<-tcr_ts1_info$Var1[tcr_ts1_info$RB>=10& tcr_ts1_info$RT>=10]

sr_tmp_cltp<-subset(sr, Sequence==rep_clone)
Idents(sr_tmp_cltp)<- sr_tmp_cltp$sample_type
rep_clone_rt_vs_rb_deg<-FindMarkers(sr_tmp_cltp, ident.1 = "RT", ident.2 = "RB");
rep_clone_rt_vs_rb_deg$genes<-rownames(rep_clone_rt_vs_rb_deg)







## trajectory analysis
#######################################################################################################



rt_sr<-subset(sr@meta.data, sample_type=="RT")
rb_sr<-subset(sr@meta.data, sample_type=="RB")
intersected_clns<-intersect(unique(rt_sr$Sequence), unique(rb_sr$Sequence))

subset_sr<-sr[,rownames(sr@meta.data%>%filter(Sequence %in% intersected_clns))]; sr
subset_sr<-subset_sr[,rownames(subset_sr@meta.data[!is.na(subset_sr$Sequence),])]
subset_sr<-subset(subset_sr, cell_type!="Tp")

genes<-subset_sr@assays$RNA@counts@Dimnames[[1]]; length(genes)
barc <-subset_sr@assays$RNA@counts@Dimnames[[2]]; length(barc)
mat <- subset_sr@assays$RNA@counts; dim(mat)

fd <- new("AnnotatedDataFrame", data = data.frame("gene_short_name"=genes, row.names=genes))
pd <- new("AnnotatedDataFrame", data = data.frame("barcode"=barc, row.names=barc))
cds <- newCellDataSet(mat, phenoData = pd, featureData = fd);cds

pData(cds)[,c("orig.ident","sample_type", "cell_type", "CD39_PD1", "pathology","TS", "Sequence")]<-
  c(as.factor(subset_sr$orig.ident), as.factor(subset_sr$sample_type), as.factor(subset_sr$cell_type), 
    as.factor(subset_sr$CD39_PD1), as.factor(subset_sr$pathology), as.factor(subset_sr$TS), as.factor(subset_sr$Sequence))


cds <- estimateSizeFactors(cds) 
cds <- estimateDispersions(cds)
cds <- detectGenes(cds, min_expr = 0.1) 
cth <- newCellTypeHierarchy()
cds <- classifyCells(cds, cth, 0.1)

expressed_genes <- row.names(subset(fData(cds), num_cells_expressed > 10)) 
diff_test_res <- differentialGeneTest(cds[expressed_genes,], fullModelFormulaStr = "~TS")
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))
length(ordering_genes) 


for( jj in dir("F:/CRC/ref/removed_geneset/")){
  rm.genes<-read.csv(paste0("F:/notebook_return/CRC/ref/removed_geneset/", jj))
  ordering_genes<-setdiff(ordering_genes, rm.genes[,1])
}

cds <- setOrderingFilter(cds, ordering_genes)
plot_ordering_genes(cds)

cds <- reduceDimension(cds, max_components = 2, method = 'DDRTree')
cds <- orderCells(cds,reverse = F)

plot_cell_trajectory(cds, color_by = "Pseudotime",show_branch_points = F,show_tree = F,cell_size =1 )+coord_equal() 





# top genes highly correlated with pseudotime 

fdata<-fData(cds)
fdata<-fdata[(fdata$use_for_ordering),]

total_cor<-c()
for(gg in fdata$gene_short_name){
  cor_p<-cor(cds$Pseudotime, cds@assayData$exprs[gg,], method = "pearson") 
  tmp_cor<-data.frame(genes=gg, pearson_cor=cor_p)
  total_cor<-rbind(total_cor, tmp_cor)
}

top_cor_p<-subset(total_cor, pearson_cor>= quantile(total_cor$pearson_cor, 0.95))
top_cor_p<-top_cor_p[rev(order(top_cor_p$pearson_cor)),]



