source("~/script/demuxer/src/demux_core_functions.R")


#load C4 sc Data
srt<-readC4('~/data/C4_out/LVW/output/filter_matrix/',0.05,500)
dimplot_publication(srt)
Cells(srt)->valid.cell

#load C4 XY ratio data
demux.info<-read.delim('../C4_out/LVW/output/chrXY.count.final.tsv',header = F)
demux.info<-spread(demux.info,key = V2,value=V3)
colnames(demux.info)[1]<-"barcode"
demux.info<-demux.info[!is.na(demux.info$chrX)&!is.na(demux.info$chrY),]
demux.info$ratioXY<-demux.info$chrX/demux.info$chrY
demux.info$sumXY<-demux.info$chrX+demux.info$chrY
nrow(demux.info)
demux.info[demux.info$barcode %in% valid.cell,]

#observed continuous distribution of ratioXY
density(demux.info[demux.info$ratioXY<=40&demux.info$barcode %in% valid.cell,"ratioXY"]) %>% plot()


#take a look at distribution of ratioXY
rownames(demux.info)<-demux.info$barcode
srt<-AddMetaData(srt,demux.info)
FeaturePlot(srt,"ratioXY",max.cutoff = 30)

#wondering if ratioXY varies across cell types, however, no obvious variation
ggplot(srt@meta.data)+geom_density(aes(x=ratioXY))+facet_wrap(~seurat_clusters)+xlim(c(0,40))

#load independent data from WX
srt_WX<-readC4('../C4_out/WX1/output/filter_matrix/',0.05,1000)
dimplot_publication(srt_WX)
Cells(srt_WX)->valid.cell.WX

demux.info_WX<-read.delim('../C4_out/WX1/output/chrXY.count.final.tsv',header = F)
demux.info_WX<-spread(demux.info_WX,key = V2,value=V3)
colnames(demux.info_WX)[1]<-"barcode"
demux.info_WX<-demux.info_WX[!is.na(demux.info_WX$chrX)&!is.na(demux.info_WX$chrY),]
demux.info_WX$ratioXY<-demux.info_WX$chrX/demux.info_WX$chrY
demux.info_WX$sumXY<-demux.info_WX$chrX+demux.info_WX$chrY

rownames(demux.info_WX)<-demux.info_WX$barcode
srt_WX<-AddMetaData(srt_WX,demux.info_WX)

ggplot(srt_WX@meta.data)+geom_density(aes(x=ratioXY))+facet_wrap(~seurat_clusters)+xlim(c(0,40))


srt<-readC4withXYRatio("../C4_out/LVW/output/",nCount = 200)
srt_WX<-readC4withXYRatio("../C4_out/WX1/output/",nCount = 200)
srt_KD.3.end<-readC4withXYRatio("../C4_out/KD/output/",nCount = 200)
srt_KD.5.end<-readC4withXYRatio("../C4_out/KD1/output/",nCount = 200)
srt_LV2.3.end<-readC4withXYRatio("../C4_out/LV2/output/",nCount = 200)

XYRatioDensityPlot(srt_KD.3.end,"KD_3_end")+XYRatioDensityPlot(srt_KD.5.end,"KD_5_end")+
  XYRatioDensityPlot(srt_WX,"WX_5_end")+XYRatioDensityPlot(srt,"LV+WX_3_end")+
  XYRatioDensityPlot(srt_LV2.3.end,"LV+2_3_end")

library(Rsamtools)
library(GenomicRanges)
library(Biostrings)  # For reverse-complement conversion

# Define BAM file
bamfile <- "data/C4_out/LVW/output/chrM_2380.bam"
bamfile <- "data/C4_out/LVW/output/chrM_5417.bam"

# Define genomic position of interest
target_chrom <- "chrM"
target_pos <- 2380
target_pos <- 5417

which_region <- GRanges(target_chrom, IRanges(target_pos, target_pos))

# Extract reads with sequences, positions, CIGAR, and strand
param <- ScanBamParam(which = which_region, what = c("qname", "seq", "pos", "cigar", "strand"))

bam_data <- scanBam(bamfile, param=param)[[1]]

# Extract fields
read_names <- bam_data$qname
read_sequences <- as.character(bam_data$seq)
read_positions <- bam_data$pos
cigar_strings <- bam_data$cigar
read_strands <- bam_data$strand  # "+" or "-"


# Apply function to each read
bases <- mapply(extract_base_corrected, read_sequences, read_positions, cigar_strings, 
                read_strands, MoreArgs = list(target_genomic_pos = target_pos))

# Create a data frame with the results
results <- data.frame(Read_Name = read_names,Cigar=cigar_strings, Aligned_Start = read_positions, Strand = read_strands, Base_at_Target = bases)

# Print results
results$Base_at_Target %>% table

read2cell_id<-read.delim('data/C4_out/LVW/output/read2cell_id.tsv',header=F,col.names = c("Read_Name","Cell_ID"))
read2cell_id<-read.delim('data/C4_out/LVW/output/read2cell_id.5417.tsv',header=F,col.names = c("Read_Name","Cell_ID"))

results<-dplyr::filter(results,!is.na(Base_at_Target))[c("Read_Name","Base_at_Target")] %>% unique()
results<-left_join(results,read2cell_id)
results<-results %>% group_by(Cell_ID,Base_at_Target) %>% summarise(count=n()) %>% arrange(Cell_ID)


cell_id_reads<-results[results$Base_at_Target %in% c("C","T"),"Cell_ID"] %>% unlist() %>% as.vector()
cell_id_reads<-results[results$Base_at_Target %in% c("G","A"),"Cell_ID"] %>% unlist() %>% as.vector()

length(valid.cell)
intersect(cell_id_reads,valid.cell) %>% length()

results.filtered<-results[results$Base_at_Target %in% c("C","T")&results$Cell_ID %in% intersect(cell_id_reads,valid.cell),]
results.filtered.total_count<-results.filtered %>% group_by(Cell_ID) %>% summarise(total_count=sum(count))
results.filtered<-left_join(results.filtered,results.filtered.total_count)
results.filtered<-dplyr::filter(results.filtered,total_count>=2)
results.filtered$af<-results.filtered$count/results.filtered$total_count
results.filtered<-dplyr::filter(results.filtered,af>0.5)

results.filtered.as.Meta<-results.filtered[c(1,2,5)] %>% as.data.frame() %>% `rownames<-`(results.filtered$Cell_ID)
srt<-AddMetaData(srt,results.filtered.as.Meta)
DimPlot(subset(srt,Base_at_Target %in% c("C","T")),group.by = "Base_at_Target")+ggtitle("base at chrM:2380")
DimPlot(subset(srt,Base_at_Target %in% c("C","T")))


#try to use predetermined SNP to demux LVW
WX_af<-readAFFileWithFiltering("~/data/C4_out/WX1/output/",logR.threshould = NULL)
WX_af$base<-ifelse(WX_af$logR>0,WX_af$alt,WX_af$ref)

LVW_af<-readAFFileWithFiltering("~/data/C4_out/LVW/output/",logR.threshould = NULL)

#LVW_af<-clusterLogRByNumSample(LVW_af,2)
#WX_af<-clusterLogRByNumSample(WX_af,1)

#plotLogRCluster(LVW_af)
#plotLogRCluster(WX_af)

#LVW_af[LVW_af$cluster==2,]
#WX_af[WX_af$cluster==2,]

#load cell level af file
LVW_sc_af<-readAFFileWithFiltering("~/data/C4_out/LVW/output/",logR.threshould = NULL,sc.mode = T)
LVW_sc_af$base<-ifelse(LVW_sc_af$logR>0,LVW_sc_af$alt,LVW_sc_af$ref)
LVW_sc_af<-LVW_sc_af[!is.na(LVW_sc_af$base),]

WX_sc_af<-readAFFileWithFiltering("~/data/C4_out/WX1/output/",logR.threshould = NULL,sc.mode = T)
WX_sc_af$base<-ifelse(WX_sc_af$logR>0,WX_sc_af$alt,WX_sc_af$ref)
WX_sc_af<-WX_sc_af[!is.na(WX_sc_af$base),] 

LVW_af.filter=LVW_af %>% dplyr::filter(.,abs(logR)<0.5) %>% arrange(desc(n.read)) %>% head(20)
WX_af.filter=WX_af %>% arrange(desc(n.read)) %>% head(10)

LVW_sc_af.filter
#test SNP imputation
test.table=constructBackbone(LVW_sc_af[LVW_sc_af$pos %in% LVW_af.filter$pos,],gather.result = T)
plotImputedBackbone(test.table[test.table$pos%in%shared.pos,])+ggtitle("LVW")

test.table.WX<-constructBackbone(WX_sc_af[WX_sc_af$pos %in% LVW_af.filter$pos,],gather.result = T,debug.cell = "CELL1008_N3")
plotImputedBackbone(test.table.WX[test.table.WX$pos%in%shared.pos,])+ggtitle("WX_ref")


test.table.spread<-spread(test.table,key=pos,value=base)
rownames(test.table.spread)<-test.table.spread$cell_id
test.table.spread$cell_id<-NULL
colnames(test.table.spread)<-paste0("position.",colnames(test.table.spread))


test.table.WX.spread<-spread(test.table.WX,key=pos,value=base)
rownames(test.table.WX.spread)<-test.table.WX.spread$cell_id
test.table.WX.spread$cell_id<-NULL
colnames(test.table.WX.spread)<-paste0("position.",colnames(test.table.WX.spread))

sorted_colnames<-colnames(test.table.spread)
sorted_colnames<-sorted_colnames[order(as.numeric(gsub("position.", "", sorted_colnames)))]

test.table.spread<-test.table.spread[sorted_colnames]
test.table.WX.spread<-test.table.WX.spread[sorted_colnames]
rownames(test.table.WX.spread)<-paste0("WX_",rownames(test.table.WX.spread))

combined.SNP.table<-rbind(test.table.spread,test.table.WX.spread[sample(rownames(test.table.WX.spread),1000,F),])

# Compute genetic distance
snp_dist <- dist.gene(combined.SNP.table, method = "percentage")
hc <- hclust(snp_dist, method = "ward.D2")
clusters <- cutree(hc, k = 2)

WX_demuxed<-clusters[clusters==1&!grepl("WX_",names(clusters))] %>% names()
LV_demuxed<-clusters[clusters==2&!grepl("WX_",names(clusters))] %>% names()
WX_spikein<-(names(clusters) %>% grep("WX_",.,value = T) %>% gsub("WX_","",.))

WX_demuxed_table<-LVW_sc_af[LVW_sc_af$cell_id %in% WX_demuxed[1:25],]
WX_demuxed_table$demux<-"WX_demuxed"
LV_demuxed_table<-LVW_sc_af[LVW_sc_af$cell_id %in% LV_demuxed[1:25],]
LV_demuxed_table$demux<-"LV_demuxed"
WX_spikein_table<-WX_sc_af[WX_sc_af$cell_id %in% WX_spikein[1:25],]
WX_spikein_table$demux<-"WX_spikein"

demux.table<-Reduce(rbind,list(WX_demuxed_table,LV_demuxed_table,WX_spikein_table))
demux.table$ident<-"demux"
demux.table$pos<-as.character(demux.table$pos)
ggplot(demux.table[demux.table$pos %in% LVW_af.filter$pos,])+geom_raster(aes(x=pos,y=cell_id,fill=base))+
  geom_raster(aes(x=ident,y=cell_id,fill=demux))+scale_y_discrete(limits=c(WX_spikein[1:25],WX_demuxed[1:25],LV_demuxed[1:25]))+
  theme(axis.text.x=element_text(angle=60,hjust=1))

srt<-subset(srt,cells=c(WX_demuxed,LV_demuxed))
srt@meta.data$demux<-ifelse(rownames(srt@meta.data) %in% WX_demuxed,"WX_demuxed","LV_demuxed")
VlnPlot(srt,features = male_genes,group.by = "demux")

#add module score to cells
srt<-AddModuleScore(srt,features = list(male_genes), name = "Male_Score")
srt<-AddModuleScore(srt,features = list(female_genes), name = "Female_Score")
VlnPlot(srt,features = "Male_Score1",group.by = "orig.ident")
VlnPlot(srt,features = "Female_Score1",group.by = "orig.ident")

#check 10x data
position_with_individual_polymorphism<-NULL
for(sample in TenX_samples){
  table.this=readAFFileWithFiltering(paste0("~/data/cellranger_out/",sample,"/outs/per_sample_outs/",sample,"/count/"),mask.poly = F)
  table.this=table.this %>% dplyr::filter(abs(logR)<abs(log2(0.1)))
  if(is.null(position_with_individual_polymorphism)){
    position_with_individual_polymorphism=table.this
  }else{
    position_with_individual_polymorphism=rbind(position_with_individual_polymorphism,table.this)
  }
}
position_with_individual_polymorphism<-position_with_individual_polymorphism$pos %>% table
position_with_individual_polymorphism<-position_with_individual_polymorphism[position_with_individual_polymorphism>=2] %>% as.data.frame() %>% `colnames<-`(.,c("pos","freq")) 
ggplot(position_with_individual_polymorphism)+geom_bar(aes(x=pos,y=freq),stat="identity")+theme(axis.text.x = element_text(angle=60,hjust=1))+
  ggtitle("positions with polymorphism in 12 individual")

LJQ_af<-readAFFileWithFiltering("~/data/cellranger_out/TFSH190500F_LJQ_3/outs/per_sample_outs/TFSH190500F_LJQ_3/count/")

LJQ_sc_af<-readAFFileWithFiltering("~/data/cellranger_out/TFSH190500F_LJQ_3/outs/per_sample_outs/TFSH190500F_LJQ_3/count/",sc.mode = T)
LJQ_sc_af$cell_id<-gsub("CB:","",LJQ_sc_af$cell_id)
LJQ_sc_af<-dplyr::filter(LJQ_sc_af,pos %in% LJQ_af$pos)

LJQ_af<-callConsensusSeq(LJQ_af)
LJQ_sc_af<-callConsensusSeq(LJQ_sc_af)
LJQ_sc_af.backbone<-constructBackbone(LJQ_sc_af,gather.result = T)
LJQ_sc_af.backbone.false<-LJQ_sc_af.backbone[LJQ_sc_af.backbone$base!=LJQ_sc_af.backbone$consensus.base,]
LJQ_sc_af.backbone.false$pos<-LJQ_sc_af.backbone.false$pos %>% as.numeric()
LJQ_sc_af.backbone.false<-left_join(LJQ_sc_af.backbone.false,LJQ_sc_af,by=c("cell_id","pos"))
LJQ_sc_af.backbone.false[is.na(LJQ_sc_af.backbone.false$base.y),]

standard<-LJQ_af[c("pos","base")]
colnames(standard)[2]<-"consensus.base"
standard$pos<-as.factor(standard$pos)
LJQ_sc_af.backbone<-left_join(LJQ_sc_af.backbone,standard)

LJQ_srt<-readRDS('data/analysis/TFSH190500F_LJQ_3.levelTop.rds')
valid_cells<-Cells(LJQ_srt)
LJQ_sc_af<-dplyr::filter(LJQ_sc_af,cell_id %in% valid_cells)
seurat_cluster_LJQ<-LJQ_srt$seurat_clusters %>% as.data.frame()
colnames(seurat_cluster_LJQ)<-"seurat_cluster"
seurat_cluster_LJQ<-rownames_to_column(seurat_cluster_LJQ,"cell_id")
LJQ_sc_af<-left_join(LJQ_sc_af,seurat_cluster_LJQ)
LJQ_sc_af.by.cluster<-LJQ_sc_af %>% group_by(pos,seurat_cluster) %>% summarise(af=sum(af),baf=sum(baf)) %>% mutate(logR=log2(baf/af))
plotPositionWithBase(LJQ_sc_af.by.cluster,by = "seurat_cluster")
head(seurat_cluster_LJQ)

#mix samples LJQ_ZF
LJQ_ZF_af<-readAFFileWithFiltering("~/data/demux/")
LJQ_ZF_af<-callConsensusSeq(LJQ_ZF_af)
LJQ_sc_af<-readAFFileWithFiltering("~/data/cellranger_out/TFSH190500F_LJQ_3/outs/per_sample_outs/TFSH190500F_LJQ_3/count/",sc.mode = T)
LJQ_sample<-LJQ_sc_af$cell_id %>% unique %>% sample(.,1000)
LJQ_sc_af<-LJQ_sc_af[LJQ_sc_af$cell_id %in% LJQ_sample,]
LJQ_sc_af<-LJQ_sc_af[LJQ_sc_af$pos %in% LJQ_ZF_af$pos,]
LJQ_sc_af$cell_id<-gsub("CB:","LJQ_",LJQ_sc_af$cell_id)


ZF_sc_af<-readAFFileWithFiltering("~/data/cellranger_out/TF_ZF_3/outs/per_sample_outs/TF_ZF_3/count/",sc.mode = T)
ZF_sample<-ZF_sc_af$cell_id %>% unique %>% sample(.,1000)
ZF_sc_af<-ZF_sc_af[ZF_sc_af$cell_id %in% ZF_sample,]
ZF_sc_af<-ZF_sc_af[ZF_sc_af$pos %in% LJQ_ZF_af$pos,]
ZF_sc_af$cell_id<-gsub("CB:","ZF_",ZF_sc_af$cell_id)

merged_LJQ_ZF_sc_af<-rbind(LJQ_sc_af,ZF_sc_af)
merged_LJQ_ZF_sc_af<-callConsensusSeq(merged_LJQ_ZF_sc_af)
plotImputedBackbone(merged_LJQ_ZF_sc_af)
merged_LJQ_ZF_sc_af<-constructBackbone(merged_LJQ_ZF_sc_af,gather.result = T)
plotImputedBackbone(merged_LJQ_ZF_sc_af)

merged_LJQ_ZF_sc_af<-clusterByConsensusSeq(merged_LJQ_ZF_sc_af,n.cluster = 2)
merged_LJQ_ZF_sc_af$true_cluster<-ifelse((grepl("^LJQ_",merged_LJQ_ZF_sc_af$cell_id)&merged_LJQ_ZF_sc_af$cluster==1)|
                                           (grepl("^ZF_",merged_LJQ_ZF_sc_af$cell_id)&merged_LJQ_ZF_sc_af$cluster==2),
                                         "yes","no")
merged_LJQ_ZF_sc_af$true_cluster %>% table()

#try to demux LVW
LVW_sc_af<-readAFFileWithFiltering("~/data/C4_out/LVW/output/",sc.mode = T,min.n.read = 10)
LVW_sc_af<-callConsensusSeq(LVW_sc_af)
getInformativePositionsByCell(LVW_sc_af,min.cell.count = 3)

plotImputedBackbone(LVW_sc_af)
LVW_sc_af.bb<-constructBackbone(LVW_sc_af,gather.result = T)
plotImputedBackbone(LVW_sc_af.bb)

LVW_sc_af.bb<-clusterByConsensusSeq(LVW_sc_af.bb,n.cluster = 2)
LVW_sc_af.bb$cluster %>% table()

#try to demux LV2
LV2_af<-readAFFileWithFiltering("~/data/C4_out/LV2/output/")
LV2_sc_af<-readAFFileWithFiltering("~/data/C4_out/LV2/output/",sc.mode = T,min.n.read = 2)

LV2_af<-callConsensusSeq(LV2_af)
LV2_sc_af<-callConsensusSeq(LV2_sc_af)
LV2_sc_af<-LV2_sc_af[LV2_sc_af$pos %in% LV2_af$pos,]
plotImputedBackbone(LV2_sc_af)
LV2_sc_af.bb<-constructBackbone(LV2_sc_af,gather.result = T)
plotImputedBackbone(LV2_sc_af.bb)

LV2_sc_af.bb<-clusterByConsensusSeq(LV2_sc_af.bb,n.cluster = 2)
LV2_sc_af.bb$cluster %>% table()

#try to demux 4 10x samples
phased_mt<-readAFFileWithFiltering('~/data/demux/LJQ_ZF_ZYX_QGG/LJQ_ZF_ZYX_QGG.mito_reads.downsample.allele.freq.tsv')
sc_merged0<-mixCells10x(TenX_samples[c(1,5,9,12)],positions=phased_mt$pos)
sc_merged<-sc_merged0
phased_mt<-callConsensusSeq(phased_mt)
sc_merged<-callConsensusSeq(sc_merged)
plotImputedBackbone(sc_merged)
sc_merged<-constructBackbone(sc_merged,gather.result = T)
plotImputedBackbone(sc_merged)

#clustering by consensus seq
sc_merged<-clusterByConsensusSeq(sc_merged,n.cluster = 4)

#add sample info
sc_merged$sample<-getField(sc_merged$cell_id,"_",1)

#look for accuracy
(sc_merged[c(1,4,5)] %>% unique())[c(2,3)] %>% table()

#count informative positions
countInformativePositions(sc_merged) %>% length() # 70 informative positions

sc_merged0<-callConsensusSeq(sc_merged0)
plotImputedBackbone(sc_merged0)

#try to demux QGG and ZF
phased_mt<-readAFFileWithFiltering('~/data/demux/QGG_ZF/QGG_ZF.allele.freq.tsv')
sc_merged0<-mixCells10x(TenX_samples[c(5,9)],positions=phased_mt$pos)
sc_merged<-sc_merged0

phased_mt<-callConsensusSeq(phased_mt)
sc_merged<-callConsensusSeq(sc_merged)

plotImputedBackbone(sc_merged)
sc_merged.bb<-constructBackbone(sc_merged,gather.result = T)
plotImputedBackbone(sc_merged.bb)
sc_merged.bb<-clusterByConsensusSeq(sc_merged.bb,n.cluster = 2)
sc_merged.bb$sample<-getField(sc_merged.bb$cell_id,"_",1)
(sc_merged.bb[c(1,4,5)] %>% unique())[c(2,3)] %>% table()

sc_merged$sample<-getField(sc_merged$cell_id,"_",1)
countInformativePositions(sc_merged) %>% length()  # 30 informative positions

#demux 3 samples
phased_mt<-readAFFileWithFiltering('~/data/demux/ZYX_XAH_ZXL/ZYX_XAH_ZXL.allele.freq.tsv')
sc_merged0<-mixCells10x(TenX_samples[grepl("ZYX|XAH|ZXL",TenX_samples)],positions=phased_mt$pos)
sc_merged<-sc_merged0
sc_merged<-callConsensusSeq(sc_merged)
sc_merged$sample<-getField(sc_merged$cell_id,"_",1)
countInformativePositions(sc_merged) %>% length() # 47 informative positions
I.pos<-countInformativePositions(sc_merged)
sc_merged0$pos<-as.factor(sc_merged0$pos)
sc_merged0[sc_merged0$pos %in% I.pos,] %>%  
  ggplot(.)+geom_violin(aes(x=as.factor(pos),y=log10(n.read)))+xlab("informative positions")+
  theme(axis.text.x = element_text(angle=60,hjust=1))+ggtitle("read coverage across informative positions")

#demux 5 samples
phased_mt<-readAFFileWithFiltering('~/data/demux/ZYX_XAH_ZXL_LJQ_WLG/ZYX_XAH_ZXL_LJQ_WLG.allele.freq.tsv')
sc_merged0<-mixCells10x(TenX_samples[grepl("ZYX|XAH|ZXL|LJQ|WLG",TenX_samples)],positions=phased_mt$pos)
sc_merged<-sc_merged0 

sc_merged<-callConsensusSeq(sc_merged)
sc_merged$sample<-getField(sc_merged$cell_id,"_",1)
countInformativePositions(sc_merged) %>% length() # 78 informative positions
I.pos<-countInformativePositions(sc_merged)
sc_merged<-dplyr::filter(sc_merged,pos %in% I.pos)


plotAccuracyByLowerBoundN.pos(sc_merged,5)
sc_merged.bb<-demuxMockSamples(sc_merged,5,10)

(sc_merged.bb[c("cell_id","sample","cluster")] %>% unique)[c("sample","cluster")] %>% table() %>% proportions(.,margin=1) %>% rowMaxs()

plotNumCellByNumPos(sc_merged)
#plot 
plotdata<-data.frame(n.sample=c(2,3,4,5),n.informative.position=c(30,47,70,78))
ccplot(plotdata)+theme_minimal()

ggplot(testdata)+geom_boxplot(aes(x=lower.bound,y=accuracy,group=lower.bound),color="black")+theme_minimal()+
  ylab("accuracy")+xlab("lower bound for number of informative positions")+ylim(c(0,1))

#demux multiplexed 10x official
#demux 5 samples
phased_mt<-readAFFileWithFiltering('~/data/cellranger_out/10xMultiplexedHumanPBMC/outs/allele.freq.tsv')
sc_merged0<-readAFFileWithFiltering('~/data/cellranger_out/10xMultiplexedHumanPBMC/outs/allele.freq.cell.tsv',sc.mode = T)
sc_merged<-sc_merged0

sc_merged<-callConsensusSeq(sc_merged)
sc_merged<-dplyr::filter(sc_merged,n.read>10)
plotImputedBackbone(sc_merged)
sc_merged.bb<-constructBackboneMulti(sc_merged)
plotImputedBackbone(sc_merged.bb)
sc_merged.bb<-clusterByConsensusSeq(sc_merged.bb,3)


sample1_bc<-read.csv('10xMultiplexedHumanPBMC/outs/sample_bcs/SC3_v3_NextGem_DI_CellPlex_Human_PBMC_10K_PBMCs_human_1_count_sample_barcodes.csv',header = F)
sample2_bc<-read.csv('10xMultiplexedHumanPBMC/outs/sample_bcs/SC3_v3_NextGem_DI_CellPlex_Human_PBMC_10K_PBMCs_human_2_count_sample_barcodes.csv',header = F)
sc_merged<-sc_merged0
sc_merged$cell_id<-gsub("CB:","",sc_merged$cell_id)
sc_merged$sample<-case_when(sc_merged$cell_id %in% sample1_bc$V2 ~ "sample1",
                            sc_merged$cell_id %in% sample2_bc$V2 ~ "sample2")
sc_merged<-callConsensusSeq(sc_merged)
plotImputedBackbone(sc_merged,use.pos = phased_mt$pos)


countInformativePositions(sc_merged) %>% length() #oops, the 2 people do not have informative positions
I.pos<-getInformativePositionsByCell(sc_merged)
sc_merged<-dplyr::filter(sc_merged,pos %in% I.pos)
sc_merged.bb<-constructBackboneMulti(sc_merged)
plotImputedBackbone(sc_merged.bb)

#try to demux QGG and ZF with cell level SNP
sc_merged0<-mixCells10x(TenX_samples[c(5,9)])
sc_merged<-sc_merged0
sc_merged<-callConsensusSeq(sc_merged)
sc_merged$sample<-getField(sc_merged$cell_id,"_",1)
I.pos<-getInformativePositionsByCell(sc_merged)$pos %>% unique()
InformativePositions(sc_merged)
sc_merged<-dplyr::filter(sc_merged,pos %in% I.pos)

plotImputedBackbone(sc_merged)
sc_merged.bb<-constructBackboneMulti(sc_merged,gather.result = T)
plotImputedBackbone(sc_merged.bb)
sc_merged.bb<-clusterByConsensusSeq(sc_merged.bb,n.cluster = 2)
sc_merged.bb$sample<-getField(sc_merged.bb$cell_id,"_",1)
(sc_merged.bb[c(1,4,5)] %>% unique())[c(2,3)] %>% table()

sc_merged$sample<-getField(sc_merged$cell_id,"_",1)
countInformativePositions(sc_merged) %>% length()  # 30 informative positions

#mix 10 samples to calculate ED
sc_merged0<-mixCells10x(TenX_samples[1:10])
sc_merged<-sc_merged0
saveRDS(sc_merged0,"data/demux/10individualSCmix.rds")
sc_merged<-callConsensusSeq(sc_merged)
sc_merged$sample<-getField(sc_merged$cell_id,"_",1)
#I.pos<-getInformativePositionsByCell(sc_merged,min.cell.count = 40)$pos %>% unique()
I.pos<-countInformativePositions(sc_merged)
sc_merged<-dplyr::filter(sc_merged,pos %in% I.pos)
plotImputedBackbone(sc_merged,n.sample = 100)
sc_merged.bb<-constructBackboneMulti(sc_merged,gather.result = T)
plotImputedBackbone(sc_merged.bb,n.sample = 50)
sc_merged.bb<-clusterByConsensusSeq(sc_merged.bb,n.cluster = 10)
sc_merged.bb$sample<-getField(sc_merged.bb$cell_id,sep = "_",field = 1)
calculateAccuracyDemuxedMock(sc_merged.bb) 

test.out<-callConsensusSeqIndividual(sc_merged.bb)
test.out=adist(test.out)
test.out=test.out[lower.tri(test.out)]
test.out %>% table() %>% barplot(xlab = "Levenshtein distance",ylab = "# of paired individual",main = "distribution of Levenshtein distance\namong individuals")

#try demux celline
sc_merged0<-readAFFileWithFiltering('~/data/cellranger_out/10xJurkatRaji/outs/',sc.mode = T)
sc_merged<-sc_merged0
sc_merged<-callConsensusSeq(sc_merged)
sc_merged$sample<-getField(sc_merged$cell_id,"_",1)
I.pos<-filterInformativePositionsByCell(sc_merged)
sc_merged<-dplyr::filter(sc_merged,pos %in% I.pos)
I.pos2<-getInformatiavePosByPosTransMat(sc_merged)
sc_merged<-dplyr::filter(sc_merged,pos %in% I.pos2)

plotImputedBackbone(sc_merged)
sc_merged.bb<-constructBackboneMulti(sc_merged,gather.result = T)
plotImputedBackbone(sc_merged.bb)
sc_merged.bb<-clusterByConsensusSeq(sc_merged.bb,n.cluster = 2)

sc_merged.bb[c(1,4)] %>% unique() %>% dplyr::filter(.,cluster==1) %>% 
  write_delim(.,"~/data/cellranger_out/10xJurkatRaji/outs/demuxed1.tsv")
sc_merged.bb[c(1,4)] %>% unique() %>% dplyr::filter(.,cluster==2) %>% 
  write_delim(.,"~/data/cellranger_out/10xJurkatRaji/outs/demuxed2.tsv")

Jurkat<-read.csv('~/data/cellranger_out/10xJurkatRaji/outs/SC3_v3_NextGem_DI_CellPlex_Jurkat_Raji_10K_Jurkat_count_sample_barcodes.csv',header = F)
Raji<-read.csv('~/data/cellranger_out/10xJurkatRaji/outs/SC3_v3_NextGem_DI_CellPlex_Jurkat_Raji_10K_Raji_count_sample_barcodes.csv',header = F)


sc_merged.bb$cell_id<-gsub("CB:","",sc_merged.bb$cell_id)
sc_merged.bb$truth<-case_when(sc_merged.bb$cell_id %in% Jurkat$V2 ~ "Jurkat",
                              sc_merged.bb$cell_id %in% Raji$V2 ~ "Raji")
(sc_merged.bb[c("cell_id","cluster","truth")] %>% unique())[c("cluster","truth")] %>% table() %>% proportions(.,margin=1)


#try demux 8 individual
sc_merged0<-readAFFileWithFiltering('~/data/cellranger_out/SRR8281306/outs/',sc.mode = T)
sc_merged<-sc_merged0
sc_merged<-callConsensusSeq(sc_merged)
#sc_merged$sample<-getField(sc_merged$cell_id,"_",1)
sc_merged<-filterInformativePositionsByCell(sc_merged,min.cell = 1000)

plotImputedBackbone(sc_merged)
sc_merged.bb<-constructBackboneMulti(sc_merged,gather.result = T)
plotImputedBackbone(sc_merged.bb)
sc_merged.bb<-clusterByConsensusSeq(sc_merged.bb,n.cluster = 2)

#try demux 4 peopple 10x big array
sc_merged0<-readAFFileWithFiltering('~/analysis/cellranger/X5_250409MIX01/outs/per_sample_outs/X5_250409MIX01/count/',sc.mode = T)
sc_merged<-sc_merged0

sc_merged$maf.ratio<-pmax(sc_merged$af,sc_merged$baf)/sc_merged$n.read
min.maf.ratio<-sc_merged[sc_merged$n.read>20,] %>% group_by(cell_id) %>% summarise(min.maf.ratio=min(maf.ratio))
doublet.like<-min.maf.ratio[min.maf.ratio$min.maf.ratio<0.8,"cell_id"] %>% unlist %>% as.vector()
sc_merged<-sc_merged[!(sc_merged$cell_id %in%doublet.like),]

sc_merged<-filterDoublet(sc_merged)
sc_merged<-callExactBase(sc_merged)
sc_merged<-filterInformativePositionsByCell(sc_merged,min.cell = 1000)
plotImputedBackbone(sc_merged)
sc_merged.bb<-constructBackboneMulti(sc_merged,gather.result = T)
plotImputedBackbone(sc_merged.bb)
sc_merged.bb<-clusterByConsensusSeq(sc_merged.bb,n.cluster = 4)
(sc_merged.bb[c("cell_id","cluster")] %>% unique)$cluster %>% table()
saveRDS(sc_merged.bb,'analysis/cellranger/X5_250409MIX01/outs/per_sample_outs/X5_250409MIX01/count/demux.rds')
sc_merged.bb<-readRDS('analysis/cellranger/X5_250409MIX01/outs/per_sample_outs/X5_250409MIX01/count/demux.rds')