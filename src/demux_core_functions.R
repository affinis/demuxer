# Lin Lyu, 2025. 4. 28
# this file include all objects for demultiplexing with single cell data in R
# for functions defined in this source code file, to clarify relationship between them, we define following terms:

#1# dependancy: function(s) that this function would call when running, 
### if function A called function B, function B is function A's dependancy

#2# caller: function(s) that would call this function,
### if function A called function B, function A is function B's caller

#3# downstream: function(s) that would utilize the output of this function,
### if function A generates data for function B without been called by function B, function B is the downstream of function A,
### attention that a caller will utilize data from their dependancy, but such case we do not consider caller as dependancy's downstream.

#4# upstream: function(s) that generate the input of this function, 
### if function A generates data for function B without been called by function B, function A is the upstream of function B,
### again, notice that a dependancy will generate data for their caller, but such case we do not consider dependancy as caller's upstream.

#5# NSF: a flag indicating that there is no specific function being this function's caller/dependancy/downstream/upstream IN THIS FILE,
### meaning this function is a top level or bottom level function IN THIS FILE

library(stringr)
library(tidyverse)
library(parallel)
library(ape)

# workflow: readAFFileWithFiltering -> filterDoublet -> callExactBase -> filterInformativePositionsByCell -> 
#           constructBackboneMulti -> clusterByConsensusSeq   

MT_coding_region=c(3307:4262,
                  4470:5511,
                  5904:7445,
                  7586:8269,
                  8366:8572,
                  8527:9207,
                  9207:9990,
                  10059:10404,
                  10470:10766,
                  10760:12137,
                  12337:14148,
                  14149:14673,
                  14747:15887)

position_with_polymorphism<-c(297,300,1670,2442,2617,3572,3577,3583,3605,3611,3727,4264,7428,7526,8219,
                              8243,8303,8425,8463,8466,8469,9980,9983,10413,12105,12113,13761,13762,15547,
                              15666,16183)

position_with_polymorphism_strict<-c(297,300,560,1397,1533,1534,1670,2285,2442,2617,3319,3572,3577,3578,3583,3605,3611,3727,4222,4254,4264,
                                     6023,6668,7428,7526,8219,8269,
                                     8243,8303,8425,8426,8459,8463,8466,8469,8470,9979,9980,9983,10058,10413,11973,12105,12113,12720,13761,13762,
                                     13768,14047,14053,14746,15377,15547,
                                     15666,16183,16247)

# Function: getField
# vec: vector of string
# sep: separator
# field: wanted field
##
# caller: readAFFileWithFiltering
# depandency: NSF
# upstream: NSF
# downstream: NSF
getField<-function(vec,sep,field){
  if(length(field)>1){
    strsplit(vec,sep,fixed=T) %>% lapply(.,`[`,field) %>% lapply(.,paste0,collapse=sep) %>% unlist()
  }else{
    strsplit(vec,sep,fixed=T) %>% lapply(.,`[`,field) %>% unlist()
  }
}

# Function: clusterByConsensusSeq
# imputed.backbone: SNP table processed wtih 'constructBackbone' or 'constructBackboneMulti'
# n.cluster: number of samples in the multiplexed sample
# speed: for a typical demultiplexed sample with 60000-70000 cells, ~2h is needed
##
# caller: NSF
# dependency: NSF
# upstream: constructBackbone, constructBackboneMulti
# downstream: NSF
clusterByConsensusSeq<-function(imputed.backbone,n.cluster){
  imputed.backbone.bkp=imputed.backbone
  imputed.backbone=imputed.backbone[c("cell_id","pos","base")]%>% spread(.,key=pos,value=base)
  rownames(imputed.backbone)=imputed.backbone$cell_id
  imputed.backbone$cell_id=NULL
  snp_dist=dist.gene(imputed.backbone, method = "percentage")
  hc=hclust(snp_dist, method = "ward.D2")
  clusters=cutree(hc, k = n.cluster)
  clusters=clusters %>% as.data.frame()
  colnames(clusters)[1]="cluster"
  clusters=rownames_to_column(clusters,"cell_id")
  imputed.backbone.bkp=left_join(imputed.backbone.bkp,clusters,by="cell_id")
  print(head(imputed.backbone.bkp))
  return(imputed.backbone.bkp)
}

# Function: readAFFileWithFiltering
# af.file.path: string, folder containing file named "allele.freq.tsv" or 'allele.freq.cell.tsv'
# pos.list: vector<integer>, indicating positions wanted in af.file, optional
# logR.threshould: float, logR, i.e., log2(af/baf), set the cutoff for it to filter positions wanted, optional
# sc.mode: logical, if input is an af.file from single cell level or not
# min.n.read: integer, minimum number of read at a position, positions with n.read less than this value will be neglected
# coding.only: logical, indicating if only coding regions on mtDNA are considered
# mask.poly: logical, indicating if positions with individual level polymorphism should be neglected, because these positions will cause ambiguity in demultiplexing. the positions is in experimental vectors called 'position_with_polymorphism' or 'position_with_polymorphism_strict' defined in this file.
# strict.mask: logical, if TRUE, use 'position_with_polymorphism_strict', else use 'position_with_polymorphism'
##
# caller: NSF
# dependency: getField
# downstream: filterDoublet
# upstream: NSF

readAFFileWithFiltering<-function(af.file.path,pos.list=NULL,logR.threshould=NULL,sc.mode=F,min.n.read=5,coding.only=F,mask.poly=T,strict.mask=T){
  if(sc.mode){
    af.file="allele.freq.cell.tsv"
  }else{
    af.file="allele.freq.tsv"
  }

  af.raw=NULL
  if(str_ends(af.file,"/")){
    af.raw=read.delim(paste0(af.file.path,af.file),header=F)
  }else if(str_ends(af.file.path,".tsv")){
    af.raw=read.delim(af.file.path,header=F)
  }else{
    af.raw=read.delim(paste0(af.file.path,"/",af.file),header=F)
  }

  if(sc.mode){
    colnames(af.raw)=c("pos","ref","alt","af","baf","cell_id")
  }else{
    colnames(af.raw)=c("pos","ref","alt","af","baf")
  }

  print(head(af.raw))
  af.mod=NULL
  if(!is.null(pos.list)){
    af.mod=dplyr::filter(af.raw,pos %in% pos.list)
    af.mod=dplyr::filter(af.mod,nchar(af.mod$ref)==1&nchar(af.mod$alt)==1)
  }else{
    af.mod=dplyr::filter(af.raw,nchar(af.raw$ref)==1&nchar(af.raw$alt)==1)
  }
  #remove(af.raw)
  #af.mod$af=getField(af.mod$freqs,",",1) %>% as.numeric()
  #gc()
  #af.mod$baf=getField(af.mod$freqs,",",2) %>% as.numeric()
  #gc()
  af.mod[is.na(af.mod$baf),"baf"]=0
  af.mod$logR=log2(af.mod$baf/af.mod$af)
  if(!sc.mode){
    af.mod=dplyr::filter(af.mod,af!=0)
  }
  if(!is.null(logR.threshould)){
    af.mod=dplyr::filter(af.mod,abs(logR)>=logR.threshould)
  }
  af.mod$n.read=af.mod$af+af.mod$baf
  af.mod=dplyr::filter(af.mod,n.read>=min.n.read)
  if(coding.only){
    af.mod=dplyr::filter(af.mod,pos %in% MT_coding_region)
  }
  if(mask.poly){
    if(strict.mask){
      af.mod=dplyr::filter(af.mod,!(pos %in% position_with_polymorphism_strict))
    }else{
      af.mod=dplyr::filter(af.mod,!(pos %in% position_with_polymorphism))
    }

  }
  print(head(af.mod))
  return(af.mod)
}

# Function: readAFStandard
# read bulkRNASeq af file
##
# caller: NSF
# dependency: NSF
# downstream: 
# upstream: 
readAFStandard<-function(af.file.path,pos.list=NULL){
  af.res=NULL
  for(sample in list.files(af.file.path)){
    af.this=read.delim(file.path(af.file.path,sample,"allele.freq.tsv"),header=F)
    colnames(af.this)=c("pos","ref","alt")
    
    af.mod=NULL
    if(!is.null(pos.list)){
      af.mod=dplyr::filter(af.this,pos %in% pos.list)
    }
    af.mod=dplyr::filter(af.mod,nchar(af.mod$ref)==1&nchar(af.mod$alt)==1)
    remove(af.this)
    af.mod$base=ifelse(af.mod$alt=='.',af.mod$ref,af.mod$alt)
    af.mod$ref=NULL
    af.mod$alt=NULL
    colnames(af.mod)[2]=getField(sample,"_",1)
    if(is.null(af.res)){
      af.res=af.mod
    }else{
      af.res=left_join(af.res,af.mod)
    }
  }
  return(af.res)
}

# Function: callConsensusSeq
# determine the sequence of each demuxed sample
# !! notice !!, a previous function also called 'callConsensusSeq' but replaced by 'callExactBase',
# pay attention to previous code that used 'callConsensusSeq'
##
# caller: NSF
# dependency: NSF
# downstream: NSF
# upstream: clusterByConsensusSeq
callConsensusSeq<-function(clustered.backbone){
  clustered.backbone=clustered.backbone %>% group_by(pos,base,cluster) %>% summarise(count=n())
  clustered.backbone=clustered.backbone %>% group_by(pos,cluster) %>% mutate(major.base=base[which.max(count)])
  clustered.backbone=dplyr::filter(clustered.backbone,base==major.base) %>% unique()
  #return(imputed.backbone)
  consensus.seq=clustered.backbone[c("cluster","pos","base")] %>% spread(key=pos,value=base) %>% as.data.frame()
  consensus.seq$cluster=paste0("individual_",consensus.seq$cluster)
  rownames(consensus.seq)=consensus.seq$cluster
  consensus.seq$cluster=NULL
  consensus.seq=t(consensus.seq) %>% as.data.frame()
  consensus.seq=rownames_to_column(consensus.seq,var = "pos")
  consensus.seq$pos=as.integer(consensus.seq$pos)
  consensus.seq=arrange(consensus.seq,pos)
  #consensus.seq=apply(consensus.seq, 1, paste, collapse = "")
  return(consensus.seq)
}


# Function: clusterLogRByNumSample
# !!deprecated function!!
# logR.table: af table with column logR
# n.sample: cluster of logR values
##
# caller:
# dependency: 
# downstream: 
# upstream: 
clusterLogRByNumSample<-function(logR.table,n.sample){
  n.peak=NULL
  if(n.sample<=3){
    n.peak=n.sample+1
  }else{
    n.peak=4
  }
  logR.table$cluster=kmeans(logR.table$logR,centers=n.peak)$cluster %>% as.factor()
  print(head(logR.table))
  return(logR.table)
}

# Function: plotLogRCluster
# !!deprecated function!!
# logR.table: af table with column logR, should be processed with 'clusterLogRByNumSample' to add column 'cluster'
##
# caller:
# dependency: 
# downstream: 
# upstream: 
plotLogRCluster<-function(logR.table){
  pointplot=ggplot(logR.table)+geom_point(aes(x=pos,y=logR,color=cluster))+
    geom_hline(yintercept = 0,linetype=2)+
    theme_minimal()

  densityplot=ggplot(logR.table)+geom_density(aes(x=logR))+
    coord_flip()+theme(axis.text.y = element_blank(),axis.ticks.y = element_blank())+xlab("")
  ggarrange(pointplot,densityplot,ncol=2,widths=c(3,1),common.legend = T)
}

# Function: callExactBase
# af: af table, if af table is called at single cell level with multiplexed data, be careful to remove doublet
##
# caller: NSF
# dependency: NSF
# downstream: filterInformativePositionsByCell, 
# upstream: filterDoublet
callExactBase<-function(af){
  af$base=ifelse(af$af>af$baf,af$ref,af$alt)
  print(head(af))
  return(af)
}

# callConsensusSeqIndividual
# return: a string indicating individual mt identity
# imputed.backbone: af file processed with constructBackbone, only used for mock data that have column 'sample'
##
# caller: NSF
# dependency: NSF
# downstream: NSF
# upstream: constructBackboneMulti
callConsensusSeqIndividual<-function(imputed.backbone){
 imputed.backbone=imputed.backbone %>% group_by(pos,base,sample) %>% summarise(count=n())
 imputed.backbone=imputed.backbone %>% group_by(pos,sample) %>% mutate(major.base=base[which.max(count)])
 imputed.backbone=dplyr::filter(imputed.backbone,base==major.base) %>% unique()
 #return(imputed.backbone)
 consensus.seq=imputed.backbone[c("sample","pos","base")] %>% spread(key=pos,value=base)
 rownames(consensus.seq)=consensus.seq$sample
 consensus.seq$sample=NULL
 consensus.seq=apply(consensus.seq, 1, paste, collapse = "")
 return(consensus.seq)
}

# Function: constructBackbone
# reference free way to construct original SNP combination with mixed sample
# af.table: should at least contain column "cell_id", "pos" and "base"
# gather.result: logical value, if TRUE, return a table with colnames: "cell_id", "pos" and "base", else it will return a matrix with rownames refer to "cell_id" and colnames refer to "pos"
# dependancy: get_position_transition_matrix, impute_snp_by_position
# speed: ~4min/(1000*100) (cell*pos) for 1 thread
##
# caller:
# dependency: 
# downstream: 
# upstream: 
constructBackbone<-function(af.table,gather.result=T,debug.cell=NULL){
  SNP.table=af.table[c("cell_id","pos","base")] %>% spread(.,key=pos,value=base)
  #cell_id.bkp=SNP.table$cell_id
  rownames(SNP.table)=SNP.table$cell_id
  SNP.table$cell_id=NULL
  # Convert NAs to factors (if not already)
  SNP.table[] <- lapply(SNP.table, factor, levels = c("A", "C", "G", "T","N"))
  print(SNP.table[1:3,1:3])
  message("generating position_transition_matrix...")
  position_transition_matrix=get_position_transition_matrix(SNP.table)
  message("done")
  imputed_snp_matrix=impute_snp_by_position(SNP.table, position_transition_matrix,debug.cell = debug.cell)
  imputed_snp_matrix$cell_id=rownames(imputed_snp_matrix)
  if(gather.result){
    imputed_snp_matrix=gather(imputed_snp_matrix,key="pos",value="base",-cell_id)
    print(head(imputed_snp_matrix))
  }else{
    print(imputed_snp_matrix[1:5,1:5])
  }
  return(imputed_snp_matrix)
}

# Function: constructBackboneMulti
# reference free way to construct original SNP combination with mixed sample using multithread
# return_pos_trans_mat: logical, whether return position transition matrix
# other arguments refer to 'constructBackbone'
# dependancy: get_position_transition_matrix, impute_snp_by_position
# speed: ~4 min/(1000*100) (cell*pos) for 1 thread
##
# caller: NSF
# dependency: get_position_transition_matrix, impute_snp_by_position
# downstream: plotImputedBackbone
# upstream: filterInformativePositionsByCell
constructBackboneMulti<-function(af.table,gather.result=T,debug.cell=NULL,return_pos_trans_mat=F){
  SNP.table=af.table[c("cell_id","pos","base")] %>% spread(.,key=pos,value=base)
  #cell_id.bkp=SNP.table$cell_id
  rownames(SNP.table)=SNP.table$cell_id
  SNP.table$cell_id=NULL
  # Convert NAs to factors (if not already)
  SNP.table[] <- lapply(SNP.table, factor, levels = c("A", "C", "G", "T","N"))
  print(SNP.table[1:3,1:3])
  message("generating position_transition_matrix...")
  position_transition_matrix=get_position_transition_matrix(SNP.table)
  message("done")
  if(return_pos_trans_mat){
    return(return_pos_trans_mat)
  }
  # Set up a cluster
  cl=makeCluster(detectCores() - 1)
  # Export the function to the cluster
  clusterExport(cl, "impute_snp_by_position")
  #split SNP table into list
  row_list=split(SNP.table, seq(nrow(SNP.table)))
  # Use parLapply for parallel processing
  results=parLapply(cl, row_list, impute_snp_by_position,position_transition_matrices=position_transition_matrix)
  # Stop the cluster
  stopCluster(cl)
  imputed_snp_matrix=Reduce(rbind,results)
  imputed_snp_matrix$cell_id=rownames(imputed_snp_matrix)
  if(gather.result){
    imputed_snp_matrix=gather(imputed_snp_matrix,key="pos",value="base",-cell_id)
    print(head(imputed_snp_matrix))
  }else{
    print(imputed_snp_matrix[1:5,1:5])
  }
  return(imputed_snp_matrix)
}

# get_position_transition_matrix
# Function to create position-specific transition matrix
##
# caller: constructBackbone, constructBackboneMulti
# dependency: NSF
# downstream: impute_snp_by_position
# upstream: NSF
get_position_transition_matrix <- function(snp_matrix,transpos=F) {
  if(transpos){
    snp_matrix=t(snp_matrix)
  }
  n_pos <- ncol(snp_matrix)
  transition_matrices <- list()
  
  # Create a transition matrix for each pair of positions, bidirectional
  for (i in 1:n_pos) {
    for (j in 1:n_pos) {
      if(i==j){
        next
      }
      observed_pairs <- table(snp_matrix[, i], snp_matrix[, j])
      pos_i=colnames(snp_matrix)[i]
      pos_j=colnames(snp_matrix)[j]
      # Normalize to probabilities
      transition_matrix <- observed_pairs / rowSums(observed_pairs)
      
      #if no data available, use a random base
      transition_matrix[is.na(transition_matrix)]<-0
      
      # Ensure rownames/colnames match the factor levels of the SNP matrix
      rownames(transition_matrix) <- colnames(transition_matrix) <- levels(snp_matrix[, 1])
      
      # Store both directions (i->j and j->i)
      transition_matrices[[paste0(pos_i, "_to_", pos_j)]] <- transition_matrix
      #transition_matrices[[paste0(pos_j, "_to_", pos_i)]] <- t(transition_matrix)  # transpose for reverse direction
    }
  }
  
  return(transition_matrices)
}

# impute_snp_by_position
# by default cols are positions, rows are cells, if not, use `transpos=T`
##
# caller: constructBackboneMulti
# dependency: NSF
# downstream: NSF
# upstream: get_position_transition_matrix
impute_snp_by_position <- function(snp_matrix, position_transition_matrices,transpos=F,debug.cell=NULL) {
  if(transpos){
    snp_matrix=t(snp_matrix)
  }
  n_pos <- ncol(snp_matrix)
  for (i in 1:nrow(snp_matrix)) {
    if(i%%1000==0){
      message(paste0("processed ",i," cells"))
    }
    if(!is.null(debug.cell)&&rownames(snp_matrix[i, ])==debug.cell){
      message(paste0("Cell: ",debug.cell," started imputing"))
    }
    for (j in 1:n_pos) {
      if (is.na(snp_matrix[i, j])) {
        #loop all available positions to calculate probability score for each base at position j
        scores_this_pos=rep(0,4)
        names(scores_this_pos)=c("A","T","C","G")
        if(!is.null(debug.cell)&&rownames(snp_matrix[i, ])==debug.cell){
          message(paste0("Position: ",colnames(snp_matrix)[j]))
        }
        for (k in 1:n_pos) {
          if (j != k && !is.na(snp_matrix[i, k])) {
            # Create the key for position-pair transition matrix (both directions)
            pos_j=colnames(snp_matrix)[j]
            pos_k=colnames(snp_matrix)[k]
            key.this <- paste0(pos_j, "_to_", pos_k)
            transition_matrix <- position_transition_matrices[[key.this]]
            prev_base <- snp_matrix[i, k]
            probs <- transition_matrix[, prev_base]
            scores_this_pos=scores_this_pos+probs[names(scores_this_pos)]

          }#if (j != k && !is.na(snp_matrix[i, k]))
        }#for (k in 1:n_pos)
        if(!is.null(debug.cell)&&rownames(snp_matrix[i, ])==debug.cell){
          print(scores_this_pos)
        }
        if(scores_this_pos[which.max(scores_this_pos)] %in% scores_this_pos[duplicated(scores_this_pos)]){
          snp_matrix[i, j]="N"
        }else{
          snp_matrix[i, j] <- names(scores_this_pos[which.max(scores_this_pos)])
        }
      }
    }
  }
  return(snp_matrix)
}


# Function: filterInformativePositionsByCell
# filter positions to keep positions that have polymorphism between cells
# mixed.sc.snp.backbone: af table with base at each position determined
# min.cell: min number of cells containing a different base that consider as a polymorphism event
##
# caller: NSF
# dependency: NSF
# downstream: plotImputedBackbone, constructBackboneMulti
# upstream: callExactBase
filterInformativePositionsByCell<-function(mixed.sc.snp.backbone,min.cell=100){
  summarised_table=mixed.sc.snp.backbone %>% group_by(pos,base) %>% summarise(n.max=n())
  divergent.pos=summarised_table$pos[summarised_table$pos %>% duplicated()]
  summarised_table=summarised_table[summarised_table$pos %in% divergent.pos,] %>% 
    group_by(pos) %>% mutate(min.count=n.max[which.min(n.max)])
  #return(summarised_table)
  
  summarised_table=dplyr::filter(summarised_table,min.count>=min.cell)
  I.pos=unique(summarised_table$pos)
  res=dplyr::filter(mixed.sc.snp.backbone,pos %in% I.pos)
  return(res)
}

# Function: plotImputedBackbone
##
# caller: NSF
# dependency: NSF
# downstream: NSF
# upstream: filterInformativePositionsByCell, constructBackboneMulti
plotImputedBackbone<-function(imputed.table,n.sample=25,use.pos=NULL,seed=2025){
  set.seed(seed)
  cell_ids=unique(imputed.table$cell_id)
  if(is.null(use.pos)){
    sorted_positions=unique(imputed.table$pos)
    sorted_positions=sorted_positions[order(as.numeric(sorted_positions))] %>% as.factor()
  }else{
    sorted_positions=use.pos
    sorted_positions=sorted_positions[order(as.numeric(sorted_positions))] %>% as.factor()
  }
  imputed.table$pos=as.factor(imputed.table$pos)
  ggplot(imputed.table[imputed.table$cell_id %in% sample(cell_ids,n.sample),])+
    geom_tile(aes(x=pos,y=cell_id,fill=base),color="black")+scale_x_discrete(limits=sorted_positions)+
    theme(axis.text.x=element_text(angle=60,hjust=1))
}

# Function: filterDoublet
##
# caller: NSF
# dependency: NSF
# downstream: callExactBase
# upstream: readAFFileWithFiltering
filterDoublet<-function(af){
  af$maf.ratio=pmax(af$af,af$baf)/af$n.read
  min.maf.ratio=af[af$n.read>20,] %>% group_by(cell_id) %>% summarise(min.maf.ratio=min(maf.ratio))
  doublet.like=min.maf.ratio[min.maf.ratio$min.maf.ratio<0.8,"cell_id"] %>% unlist %>% as.vector()
  af=af[!(af$cell_id %in%doublet.like),]
  return(af)
}

# Function: calculateDistBetweenSampleAndStandard
# sample.seqs: col1 should be 'pos', then cols following col1 are samples
# std.seqs: col1 should be 'pos', then cols following col1 are samples
##
# caller
# dependency
# upstream: readAFStandard, callConsensusSeq
# downstream
calculateDistBetweenSampleAndStandard<-function(sample.seqs,std.seqs){
  test.seqs=sample.seqs[2:ncol(sample.seqs)] %>% apply(.,2,paste, collapse = "")
  std.seqs=std.seqs[2:ncol(std.seqs)] %>% apply(.,2,paste, collapse = "")
  return(adist(test.seqs,std.seqs))
}
