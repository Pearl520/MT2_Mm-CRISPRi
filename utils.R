#Routine R functions used in RNA-seq/ChIP analyses

inputTEcountfiles <- function(sampleName, simpleName, countDataPath){
  
  suppressMessages(library(dplyr))
  
  count <- data.frame()
  for (i in 1:length(sampleName)){
    inFile <- paste0(countDataPath, sampleName[i], ".cntTable")
    
    count.tmp <- read.table(inFile, sep = "\t", header = T, stringsAsFactors = FALSE)
    colnames(count.tmp) <- c("id", "count")
    #count.tmp <- count.tmp[,c("Gene.ID", "Gene.Name", "FPKM")]
    #for duplicated IDs, just sum up.
    #rpkm.tmp <- aggregate(FPKM~Gene.ID+Gene.Name, data = rpkm.tmp, FUN=sum)
    
    if(nrow(count) == 0){
      count <- count.tmp
      
    } else {
      count <- left_join(count, count.tmp, by = "id")
    }
  }
  colnames(count) <- c("id", as.character(simpleName))
  return(count)  
}

classify_distance <- function(tmp, ignore_direction = F){
  bin <- rep(0, length(tmp))
  
  if(ignore_direction){
    for (i in 1:length(bin)){
      if(abs(tmp[i]) >=0 & abs(tmp[i]) < 10000){ bin[i] <- "0-10kb" }
      else if (abs(tmp[i]) >= 10000 & abs(tmp[i]) < 30000){ bin[i] <- "10-30kb" }
      #else if (tmp[i] >= 20000 & tmp[i] < 30000){ bin[i] <- "20-30kb" }
      else if (abs(tmp[i]) >= 30000 & abs(tmp[i]) < 50000){ bin[i] <- "30-50kb" }
      #else if (tmp[i] >= 40000 & tmp[i] < 50000){ bin[i] <- "40-50kb" }
      else { bin[i] <- ">50kb" } 
    }
  } else{
    for (i in 1:length(bin)){
      if(tmp[i] >=0 & tmp[i] < 10000){ bin[i] <- "0-10kb" }
      else if (tmp[i] >= 10000 & tmp[i] < 30000){ bin[i] <- "10-30kb" }
      #else if (tmp[i] >= 20000 & tmp[i] < 30000){ bin[i] <- "20-30kb" }
      else if (tmp[i] >= 30000 & tmp[i] < 50000){ bin[i] <- "30-50kb" }
      #else if (tmp[i] >= 40000 & tmp[i] < 50000){ bin[i] <- "40-50kb" }
      else if (tmp[i] >= 50000) { bin[i] <- ">50kb" }
      else if (tmp[i] >= -10000 & tmp[i] < 0) { bin[i] <- "-0-10kb" }
      else if (tmp[i] >= -30000 & tmp[i] < -10000) { bin[i] <- "-10-30kb" }
      #else if (tmp[i] >= -30000 & tmp[i] < -20000) { bin[i] <- "-20-30kb" }
      else if (tmp[i] >= -50000 & tmp[i] < -30000) { bin[i] <- "-30-50kb" }
      #else if (tmp[i] >= -50000 & tmp[i] < -40000) { bin[i] <- "-40-50kb"}
      else { bin[i] <- "< -50kb"}
    }
  }
  return(bin)
}

inputStringTieRPKMfiles <- function(sampleName, simpleName, RPKMDataPath){
  
  suppressMessages(library(dplyr))
  
  rpkm <- data.frame()
  for (i in 1:length(sampleName)){
    inFile <- paste0(RPKMDataPath, sampleName[i], ".gene_abund.txt")
    
    rpkm.tmp <- read.table(inFile, sep = "\t", header = T, stringsAsFactors = FALSE)
    
    rpkm.tmp <- rpkm.tmp[,c("Gene.ID", "Gene.Name", "FPKM")]
    #for duplicated IDs, just sum up.
    rpkm.tmp <- aggregate(FPKM~Gene.ID+Gene.Name, data = rpkm.tmp, FUN=sum)
    
    if(nrow(rpkm) == 0){
      rpkm <- rpkm.tmp[,c("Gene.ID", "Gene.Name", "FPKM")] 

    } else {
      rpkm <- left_join(rpkm, rpkm.tmp[,c("Gene.ID", "FPKM")], by = c("Gene.ID"))
    }
  }
  colnames(rpkm) <- c("id", "name",
                      as.character(simpleName))
  return(rpkm)
}

countsToDEseq2FDR <- function(counts, CGroup = 0, TGroup = 0, 
                              min_read = 10){
  
  suppressMessages(library(DESeq2))
  
  counts <- counts[apply(counts, 1, function(x){sum(x)}) > min_read, ]
  groups <- factor(c(rep("CGroup",CGroup),rep("TGroup",TGroup)))
  sampleInfo <- data.frame(groups, row.names = colnames(counts))
  
  dds <- DESeqDataSetFromMatrix(countData = counts, colData = sampleInfo, design = ~ groups)
  dds$groups = relevel(dds$groups,ref="CGroup")
  dds <- DESeq(dds)
  res <- results(dds,independentFiltering=F)
  
  #normalized counts obtained based on DESeq2 "plotCounts" source code
  normCounts <- counts(dds, normalized = TRUE)
  results <- as.data.frame(cbind(normCounts, res))
  results$id <- row.names(results)  
  return (results)
}

classifyDEG <- function(df, ctr.rpkm, trt.rpkm, FDR.col, log2FC.col,  
                        log2FC = 1, RPKM = 1, FDR = 0.05){
  
  df$ctr.mean <- rowMeans(df[ctr.rpkm])
  df$trt.mean <- rowMeans(df[trt.rpkm])
  group <- rep(0, nrow(df)) 
  
  for (i in 1:nrow(df)){
    if(df$ctr.mean[i] > RPKM | df$trt.mean[i] > RPKM){
      if(df[log2FC.col][i,] > log2FC & df[FDR.col][i,] < FDR){
        group[i] <- "up-regulated"
      }
      else if (df[log2FC.col][i,] < -log2FC & df[FDR.col][i,] < FDR){
        group[i] <- "down-regulated"
      }
      else{
        group[i] <- "similar_level"
      }
    } 
    else{
      group[i] <- "low_expression_level"
    }
  }
  return(group)
}

classifyDEG_FConly <- function(df, ctr.rpkm, trt.rpkm, 
                        log2FC = 1, RPKM = 1){
  
  df$ctr.mean <- rowMeans(df[ctr.rpkm])
  df$trt.mean <- rowMeans(df[trt.rpkm])
  group <- rep(0, nrow(df)) 
  df$log2FC.col <- log2(df$trt.mean+0.01) - log2(df$ctr.mean + 0.01)
    
  for (i in 1:nrow(df)){
    if(df$ctr.mean[i] > RPKM | df$trt.mean[i] > RPKM){
      if(df$log2FC.col[i] > log2FC){
        group[i] <- "up-regulated"
      }
      else if (df$log2FC.col[i] < -log2FC){
        group[i] <- "down-regulated"
      }
      else{
        group[i] <- "similar_level"
      }
    } 
    else{
      group[i] <- "low_expression_level"
    }
  }
  return(group)
}

ggScatterplot <- function(df, x, y, group, gene,
                          my.color=c("blue", "grey50", "grey50", "red3"),
                          label.up, label.down, xlab, ylab, title,
                          genes4Label = NULL,
                          FC.line = 2){
  
  suppressMessages(library(ggplot2))
  suppressMessages(library(cowplot))
  suppressMessages(library(ggrepel))
  suppressMessages(library(gridExtra))
  suppressMessages(library(grid))
  
  g <- ggplot(df, aes_string(x=x, y=y, color=group,label=gene)) + 
    geom_point() + scale_color_manual(values=my.color) + 
    xlab(xlab) + ylab(ylab) + ggtitle(title) + theme_cowplot(13) + 
    theme(legend.position = "none") + xlim(0, 25) + ylim(0, 25) + 
    geom_abline(intercept = log2(FC.line), slope = 1, linetype = 2) + 
    geom_abline(intercept = -log2(FC.line), slope = 1, linetype = 2) + 
    geom_abline(intercept = 0, slope = 1, linetype = 1) +
    annotate("text", x = -Inf, y = Inf, label = label.up, 
             col = "red3", size = 3.8, hjust = -0.2, vjust = 1.5) + 
    annotate("text", x= Inf, y =-Inf, label = label.down, 
             col = "blue", size = 3.8, hjust = 1.2, vjust = -1) + 
    geom_text_repel(
      data = subset(df, df[, gene] %in% genes4Label),
      size = 3.8, segment.size = 0.3, segment.color = "black",
      direction = "x", nudge_y = 4.5, nudge_x =-3.5,
      point.padding = 0.25, box.padding = 0.25) +
    geom_point(
      data = subset(df, df[, gene] %in% genes4Label), 
      col = "black", size = rel(1.5))
  return(g)  
}

ggScatterplotSimple <- function(df, x, y, gene,
                                genes4Label = NULL
                                ){
  
  suppressMessages(library(ggplot2))
  suppressMessages(library(cowplot))
  suppressMessages(library(ggrepel))
  suppressMessages(library(gridExtra))
  suppressMessages(library(grid))
  
  g <- ggplot(df, aes_string(x=x, y=y, label=gene)) + 
    geom_point(color = "gray") + 
    xlab(x) + ylab(y) + 
    ggtitle(paste0(y, "_vs_", x)) + theme_cowplot(13) + 
    theme(legend.position = "none") + xlim(0, 20) + ylim(0, 20) + 
    geom_abline(intercept = 0, slope = 1, linetype = 2) + 
    annotate(geom="text", x = -Inf, y = Inf, 
             label=paste0("R=", round(cor(df[, x],
                                          df[, y]), 2)),
             color="blue", size = 3.8, hjust = -0.2, vjust = 1.5) + 
    geom_text_repel(
      data = subset(df, df[, gene] %in% genes4Label),
      size = 3.8, segment.size = 0.3, segment.color = "black",
      direction = "x", nudge_y = 4.5, nudge_x =-3.5,
      point.padding = 0.25, box.padding = 0.25) +
    geom_point(
      data = subset(df, df[, gene] %in% genes4Label), 
      col = "black", size = rel(1.5))
  return(g)  
}

plotGoBars <- function(go_res,top=10, col="red", title="top 10 enriched GO"){
  tmp <- as.data.frame(go_res)
  ord <- order(tmp$p.adjust, decreasing = F)
  tmp <- tmp[ord,]
  tmp <- head(tmp, top)
  tmp$Description <- factor(tmp$Description, levels = tmp$Description)
  g <- ggplot(tmp, aes(x=Description, y=-log10(p.adjust) )) + 
    geom_bar(stat = "identity", fill = col) + coord_flip() + 
    theme_cowplot(13) + ggtitle(title) +  ylab("") +
    scale_x_discrete(limits = rev(levels(tmp$Description))) 
  theme(plot.title = element_text(color="black",face="bold",hjust = 0.5)) 
  g
}

plotAnnoStats <- function(peaks_anno, ttl, simplify=F, category_order=NULL) {
  require(dplyr)
  anno_df <- data.frame()
  
  for(i in 1:length(peaks_anno)){
    tmp <- peaks_anno[[i]]@annoStat
    tmp$Cell <- names(peaks_anno)[i]
    tmp$Feature <- as.character(tmp$Feature)
    
    if(simplify){
      pos <- grep("Exon",tmp$Feature)
      tmp$Feature[pos] <- "Exon"
      
      pos <- grep("Intron",tmp$Feature)
      tmp$Feature[pos] <- "Intron"
      
      pos <- grep("Downstream",tmp$Feature)
      tmp$Feature[pos] <- "Distal Intergenic"
      
      pos <- grep("Promoter",tmp$Feature)
      tmp$Feature[pos] <- "Promoter"
      
      tmp$Feature <- gsub("5' UTR","Exon",tmp$Feature)
      tmp$Feature <- gsub("3' UTR","Exon",tmp$Feature)
      tmp$Feature <- gsub("Distal Intergenic","Intergenic",tmp$Feature)
      tmp$Feature <- factor(tmp$Feature, levels = rev(c("Promoter","Exon","Intron","Intergenic")))
    }
    
    anno_df <- rbind(anno_df, tmp)
  }
  
  anno_df = anno_df %>% group_by(Feature,Cell) %>% summarize(Frequency = sum(Frequency))
  
  anno_df$Feature <- factor(anno_df$Feature)
  
  if(is.null(category_order)){
    anno_df$Cell <- factor(anno_df$Cell)
  }else{
    anno_df$Cell <- factor(anno_df$Cell, levels = category_order)
  }
  
  
  g <- ggplot(anno_df, aes(x=Cell, y=Frequency, fill=Feature)) +
    coord_flip() +
    geom_bar(stat="identity") + ggthemes::theme_base() +
    #scale_fill_brewer(palette = "Set1") +
    scale_fill_npg(breaks=c("Promoter","Exon","Intron","Intergenic")) +
    ggtitle(ttl)+
    xlab("Cell") + ylab("Percentage (%)") + theme(plot.title = element_text(hjust = 0.5,face = "plain",size = 12),
                                                  axis.text = element_text(color="black")) #+ boldText
  g
}

simplifyAnnotation <- function(peaks_anno) {
  #input is a vector of annotation generated by annotatePeak fro ChIPseeker
  #output is a vector of simplified
  tmp <- peaks_anno
  tmp <- gsub("Distal Intergenic","Intergenic", tmp)
  
  pos <- grep("5' UTR",tmp)
  tmp[pos] <- "Intragenic"
  
  pos <- grep("3' UTR",tmp)
  tmp[pos] <- "Intragenic"
  
  pos <- grep("Exon",tmp)
  tmp[pos] <- "Intragenic"
  
  pos <- grep("Intron",tmp)
  tmp[pos] <- "Intragenic"
  
  pos <- grep("Downstream",tmp)
  tmp[pos] <- "Intergenic"
  
  pos <- grep("Promoter",tmp)
  tmp[pos] <- "Promoter"
  
  return (tmp)
}