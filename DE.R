#library(BiocManager)
#file.edit('~/.Rprofile')
#options(BioC_mirror=" https://mirrors.ustc.edu.cn/")
#BiocManager::install("DESeq2")
#BiocManager::install("pheatmap")
#BiocManager::install("apeglm")
#BiocManager::install("limma")
pvalue
####################################################文件读取#####################################################
input = 'J:\\data\\bioinformatic\\chinafir\\new_gtf\\stringtie_add_e_parameter\\gene_count_matrix.csv'
save_dir = "J:\\data\\bioinformatic\\chinafir\\new_gtf\\DE"

group_list=c("26CW","26OW", "2CW", "2OW", "74CW", "74OW", "CK")
num=c(2,2,2,2,2,2,2)
genes <- as.matrix(read.csv(input,row.names="gene_id", header = TRUE))
#############################################分析函数#################################################################
Batch_Deseq_differnece<-function(exprSet, group, num, save_dir="Alldiffenece", save_dir2="NEW_MA",
                                 go_file_path_bp, go_file_path_2_bp, go_file_path_mf, go_file_path_2_mf,
                                 go_file_path_cc, go_file_path_2_cc, 
                                 kegg_file_path, kegg_file_path_2, data_type = 'all'){
  ####加载R包
  library(DESeq2)   
  library(pheatmap)  # 用于作热图的包
  library(ggplot2)
  library(ggrepel)
  library(apeglm)
  library(ggVolcano)
  ########读取文件以及设置创建保存文件夹
  ##create a folder 
  save_pdf<-paste0(save_dir,"/")
  dir.create(save_pdf)
  #ma
  save_dir_ma=paste0(save_pdf,"/MA/")
  dir.create(save_dir_ma)
  #火山图保存
  save_dir_vol=paste0(save_pdf,"/volcano/")
  dir.create(save_dir_vol)
  save_dir_vol_csv=paste0(save_dir_vol,"/csv/")
  dir.create(save_dir_vol_csv)
  save_dir_vol_pdf=paste0(save_dir_vol,"/pdf/")
  dir.create(save_dir_vol_pdf)
  ##########分组创建#########################
  ## creat a group
  group_list= factor(rep(group,num))
  colData=data.frame(row.names = colnames(exprSet),
                     group_list)
  dds_cluster <- DESeqDataSetFromMatrix(countData = exprSet, colData = colData, design = ~group_list)
  dds_cluster_1 <- DESeq(dds_cluster, fitType = 'mean', minReplicatesForReplace = 7, parallel = FALSE) 
  res_cluster <- results(dds_cluster_1)
  res_cluster_1 <- data.frame(res_cluster, stringsAsFactors = FALSE, check.names = FALSE)
  # 依次按照padj值log2FoldChange值进行排序
  res_cluster_1 <- res_cluster_1[order(res_cluster_1$padj, res_cluster_1$log2FoldChange, decreasing = c(FALSE, TRUE)), ]
  # 获取padj（p值经过多重校验校正后的值）小于0.05，表达倍数取以2为对数后大于1或者小于-1的差异表达基因。
  res_cluster_up<- res_cluster_1[which(res_cluster_1$log2FoldChange >= 1 & res_cluster_1$padj < 0.05),]      # 表达量显著上升的基因
  res_cluster_down<- res_cluster_1[which(res_cluster_1$log2FoldChange <= -1 & res_cluster_1$padj < 0.05),]    # 表达量显著下降的基因
  res_cluster_total <- rbind(res_cluster_up,res_cluster_down)
  #################################################样品聚类（检测）######################################################
  library(factoextra)
  save_cluster_pdf<-paste0(save_pdf,"/cluster/")
  dir.create(save_cluster_pdf)
  
  vsd <- vst(dds_cluster, blind = FALSE)
  sampleDists <- dist(t(assay(vsd)))
  sample_cluster <- hcut(sampleDists, k = 2, stand = TRUE)
  # Visualize
  cluster_plot <- fviz_dend(sample_cluster,
                            # 加边框
                            rect = TRUE,
                            # 边框颜色
                            rect_border="cluster",
                            # 边框线条类型
                            rect_lty=2,
                            # 边框线条粗细
                            lwd=1.2,
                            # 边框填充
                            rect_fill = T,
                            # 字体大小
                            cex = 0.4,
                            # 字体颜色
                            color_labels_by_k=T,
                            # 平行放置
                            horiz=T)
  
  plot_name = paste0(save_cluster_pdf,"/","sample_cluster.pdf")
  ggsave(plot_name, cluster_plot, width = 8, height = 6)
  
  #PCA##################################################################
  save_pca_pdf<-paste0(save_pdf,"/pca/")
  dir.create(save_pca_pdf)
  
  vsd <- vst(dds_cluster)
  
  pca_1 <- plotPCA(vsd, intgroup=c("group_list"))
  plot_name <- paste0(save_pca_pdf,"/","pca_type1.pdf")
  ggsave(plot_name, pca_1, width = 6, height = 5)
  
  #PCA-2################################################################
  library("FactoMineR")
  library("factoextra")
  dat=exprSet
  dat=t(dat)
  dat=as.data.frame(dat)
  dat=cbind(dat,group_list)
  ncol(dat)
  dat[,ncol(dat)]
  dat.pca <- PCA(dat[,-ncol(dat)], graph = FALSE)
  plot(dat.pca,choix="ind")
  # Graph of individuals
  pca_2 <- fviz_pca_ind(dat.pca,
                        geom.ind = "point", # show points only (nbut not "text")
                        col.ind = dat$group_list, # color by groups
                        palette = c("#00AFBB", "#E7B800", '#CC00FF', '#FF0099', '#ECA8A9', '#D3E2B7', '#E6E0B0'),
                        addEllipses = TRUE, # Concentration ellipses # 是否圈起来
                        legend.title = "Groups"
  )
  
  plot_name = paste0(save_pca_pdf,"/","pca_type2.pdf")
  ggsave(plot_name, pca_2, width = 6, height = 5)
  
  #总体热图###############################################################
  save_heatmap_pdf<-paste0(save_pdf,"/heatmap/")
  dir.create(save_heatmap_pdf)
  df <- exprSet[intersect(rownames(exprSet),rownames(res_cluster_total)),]    
  # 在原表达矩阵中找到差异表达基因
  df2<- as.matrix(df)
  pdf(paste0(save_heatmap_pdf,"/","total_heatmap.pdf")) 
  heatmap_all <- pheatmap(df2,
                          show_rownames = F,
                          show_colnames = T,
                          cluster_cols = F,
                          cluster_rows=T,
                          height=10,  
                          scale = "row",
                          frontsize = 0.5,
                          angle_col=45, 
                          color =colorRampPalette(c("#8854d0", "#ffffff","#fa8231"))(100),
                          clustering_method = 'single',
  )
  print(heatmap_all)
  dev.off()
  ############################样本一对一比较分析######################################################################
  ## use the Deseq2 to have Diffence analyse
  for (i in 1:length(group)){
    
    name=unique(group)[i]
    print(name)
    colData$group<-relevel(colData$group,ref=name)
    dds=DESeq2::DESeqDataSetFromMatrix(countData = exprSet,
                                       colData = colData,
                                       design = ~group) 
    dds <- dds[ rowSums(DESeq2::counts(dds)) > 10, ]
    dds <- DESeq2::DESeq(dds)
    for (j in 2:length(DESeq2::resultsNames(dds))){
      
      resname=DESeq2::resultsNames(dds)[j]
      
      res=DESeq2::results(dds, name=resname)
      s <- resname
      print(resname)
      s1 <- strsplit(s, "_")
      print(s1)
      
      #if(s1[[1]][2] == 'H.CK'){
      #  if(s1[[1]][4] == 'Q.MPMN'){
      #    next
      #  }
      #}
      #else if (s1[[1]][2] == 'H.MPMN') {
      #  if(s1[[1]][4] != 'Q.MPMN'){
      #    next
      #  }
      #}
      #else if (s1[[1]][2] == 'Q.CK') {
      #  if(s1[[1]][4] != 'Q.MPMN'){
      #    next
      #  }
      #}
      #else{
      #  next
      #}
      resname <- paste0(s1[[1]][2], 'vs', s1[[1]][4])
      ######################################################数据保存#########################################################
      res_s_1 <- data.frame(res, stringsAsFactors = FALSE, check.names = FALSE)
      # 依次按照padj值log2FoldChange值进行排序
      res_s_1 <- res_s_1[order(res_s_1$padj, res_s_1$log2FoldChange, decreasing = c(FALSE, TRUE)), ]
      # 获取padj（p值经过多重校验校正后的值）小于0.05，表达倍数取以2为对数后大于1或者小于-1的差异表达基因。
      res1_s_up<- res_s_1[which(res_s_1$log2FoldChange >= 1 & res_s_1$padj < 0.05),]      # 表达量显著上升的基因
      res1_s_down<- res_s_1[which(res_s_1$log2FoldChange <= -1 & res_s_1$padj < 0.05),]    # 表达量显著下降的基因
      res1_s_total <- rbind(res1_s_up,res1_s_down)
      
      save_dir_csv<-paste0(save_dir_vol_csv,"/", s1[[1]][2], 'vs', s1[[1]][4], '//')
      dir.create(save_dir_csv)
      
      write.csv(res1_s_up, paste0(save_dir_csv,resname,"_up.csv"))
      write.csv(res1_s_down, paste0(save_dir_csv,resname,"_down.csv"))
      write.csv(res1_s_total, paste0(save_dir_csv,resname,"_DE.csv"))
      #####################################################火山图1#######################################################
      res_s_1$color <- ifelse(res_s_1$padj<0.05 & abs(res_s_1$log2FoldChange)>= 1,ifelse(res_s_1$log2FoldChange > 1,'Up','Down'),'Stable')
      color <- c(Up = "red",Stable = "gray",Down = "blue")
      
      p <- ggplot(
        # 指定数据、映射、颜色
        res_s_1, aes(log2FoldChange, -log10(padj), col = color)) +  
        geom_point() +
        theme_bw() +
        ggtitle(resname) +
        scale_color_manual(values = color) +
        # 辅助线
        labs(x="log2 (fold change)",y="-log10 (padj)") +
        geom_hline(yintercept = -log10(0.05), lty=4,col="grey",lwd=0.6) +
        geom_vline(xintercept = c(-1, 1), lty=4,col="grey",lwd=0.6) +
        scale_x_continuous(limits = c(-10, 10)) + 
        # 图例
        theme(legend.position = "none",
              panel.grid=element_blank(),
              axis.title = element_text(size = 16),
              axis.text = element_text(size = 14)+
                # 注释
                geom_text_repel(
                  data = subset(res_s_1, padj < 1e-100 & abs(res_s_1$log2FoldChange) >= 10),
                  aes(label = rownames(res_s_1)),
                  size = 5,
                  box.padding = unit(0.35, "lines"),
                  point.padding = unit(0.3, "lines"))
        )
      save_pdf<-paste0(save_dir_vol_pdf,"/", s1[[1]][2], 'vs', s1[[1]][4], '//')
      dir.create(save_pdf)
      plot_name = paste0(save_pdf,"/",resname,"_volcano_type1.pdf") 
      ggsave(plot_name, p, width = 6, height = 5)
      ########################################火山图2###########################################################
      plot_name_2 = paste0(save_pdf,"/",resname,"_volcano_type2")
      res_s_1 <- add_regulate(res_s_1, log2FC_name = "log2FoldChange",
                              fdr_name = "padj",log2FC = 1, fdr = 0.05)
      colnames(res_s_1)[1] <- "geneName"
      
      res_s_2 <- res_s_1[which(res_s_1$log2FoldChange < 10), ]
      
      res_s_2
      
      # plot
      ggvolcano(res_s_2, x = "log2FoldChange", y = "padj", 
                add_label = FALSE, legend_title = resname,
                  x_lab = 'log2FoldChange', y_lab = '-log10 (padj)',
                  filename= plot_name_2, log2FC_cut=c(-1, 1))
      
      print('volcano1 finish')
      ######################################火山图3#############################################################
      plot_name_3 = paste0(save_pdf,"/",resname,"_volcano_type3")
      gradual_volcano(res_s_2, x = "log2FoldChange", y = "padj",
                      filename=plot_name_3,log2FC_cut=c(1, 1))
      
      
      
      #####################################LFC&&&MA################################################################
      res_lfc <- lfcShrink(dds, coef=j, res=res, type="apeglm")
      res_lfc
      #res=res_lfc
      
      summary(res_lfc)
      summary(res)
      
      #dir.create(save_dir)
      write.csv(res,paste0(save_dir_csv,resname,".csv"))
      
      #dir.create(save_dir2)
      
      save_dir_MA=paste0(save_dir_ma,"/",resname)
      dir.create(save_dir_MA)
      write.csv(res,paste0(save_dir_MA,"/",resname,"_res.csv"))
      write.csv(res_lfc,paste0(save_dir_MA,"/",resname,"_reslfc.csv"))
      pdf(paste0(save_dir_MA,"/",resname,"_MA.pdf")) 
      plotMA(res, ylim=c(-5,5),main=paste0(resname," MA"), alpha=0.05)
      
      dev.off()
      pdf(paste0(save_dir_MA,"/",resname,"_MAlfc.pdf")) 
      xlim <- c(1,1e7); ylim<-c(-5,5)
      plotMA( res_lfc, xlim=xlim, ylim=ylim, main=paste0(resname," apeglm"), alpha=0.05)
      
      dev.off()
      print('ma finish')
      ##########################################热图###############################################################
      each_group <- subset(colData, group_list == s1[[1]][2] | group_list == s1[[1]][4])
      exprSet2 <- t(exprSet)
      df_pre <- exprSet2[intersect(rownames(exprSet2),rownames(each_group)),]  
      df_pre2 <- t(df_pre)
      df <- df_pre2[intersect(rownames(df_pre2),rownames(res1_s_total)),]  
      df = df[apply(df, 1, function(x) sd(x)!=0),] 
      df2<- as.matrix(df)
      pdf(paste0(save_heatmap_pdf,"/",resname,"_heatmap.pdf"))
      pheatmap(df2,
               show_rownames = F,
               show_colnames = T,
               cluster_cols = F,
               cluster_rows=T,
               height=6,  
               scale = "row",
               frontsize = 0.5,
               angle_col=45, 
               color =colorRampPalette(c("#8854d0", "#ffffff","#fa8231"))(100),
               clustering_method = 'single',
      )
      dev.off()
      print('heatmap finish')
    }
  }
}

    
################################################函数自动化过程################################################
    
Batch_Deseq_differnece(genes,group=group_list,num,save_dir = save_dir,
                       go_file_path_bp = go_file_path_bp, go_file_path_2_bp = go_file_path_2_bp,
                       go_file_path_cc = go_file_path_cc, go_file_path_2_cc = go_file_path_2_cc,
                       go_file_path_mf = go_file_path_mf, go_file_path_2_mf = go_file_path_2_mf,
                       kegg_file_path = kegg_file_path, kegg_file_path_2 = kegg_file_path_2)
