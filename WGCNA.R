 
wgcna <- function(
    file_path,
    save_path
  ){
 # BiocManager::install("WGCNA")
  library(WGCNA)
  library(reshape2)
  library(stringr)
  library(DESeq2)
  #####################创建保存文件夹##################
  #保存位置创建
  dir.create(save_path)
  #####################格式转换########################
  data_raw = as.matrix(read.csv(file_path,header = T,row.names="gene_id"))
  add_one <- function(x){
    return (x+1)
  }
  data_raw_1 = apply(data_raw, c(1:2), add_one)
  data_refomat = apply(data_raw_1, c(1:2), log2)
  # data_refomat = varianceStabilizingTransformation(data_raw)
  #################################################
  exprMat <- data_refomat
  
  # 官方推荐 "signed" 或 "signed hybrid"
  # 为与原文档一致，故未修改
  type = "unsigned"
  
  # 相关性计算
  # 官方推荐 biweight mid-correlation & bicor
  # corType: pearson or bicor
  # 为与原文档一致，故未修改
  corType = "pearson"
  
  # corFnc = ifelse(corType=="pearson", cor, bicor)
  # 对二元变量，如样本性状信息计算相关性时，
  # 或基因表达严重依赖于疾病状态时，需设置下面参数
  maxPOutliers = ifelse(corType=="pearson",1,0.05)
  
  # 关联样品性状的二元变量时，设置
  robustY = ifelse(corType=="pearson",T,F)
  
  ##导入数据##
  # dataExpr <- read.table(exprMat, sep='\t', row.names=1, header=T,quote="", comment="", check.names=F)
  dataExpr <- data_refomat
  
  dim(dataExpr)
  
  ## [1] 3600  134
  
  head(dataExpr)[,1:8]
  #############################数据筛选###################
  ## 筛选中位绝对偏差前75%的基因，至少MAD大于0.01
  ## 筛选后会降低运算量，也会失去部分信息
  ## 也可不做筛选，使MAD大于0即可
  m.mad <- apply(dataExpr,1,mad)
  dataExprVar <- dataExpr[which(m.mad >
                                  max(quantile(m.mad, probs=seq(0, 1, 0.25))[2],0.01)),]
  
  ## 转换为样品在行，基因在列的矩阵
  dataExpr <- as.data.frame(t(dataExprVar))
  
  ## 检测缺失值
  gsg = goodSamplesGenes(dataExpr, verbose = 3)
  gsg$allOK
  
  ##  Flagging genes and samples with too many missing values...
  ##   ..step 1
  
  if (!gsg$allOK){
    # Optionally, print the gene and sample names that were removed:
    if (sum(!gsg$goodGenes)>0)
      printFlush(paste("Removing genes:",
                       paste(names(dataExpr)[!gsg$goodGenes], collapse = ",")));
    if (sum(!gsg$goodSamples)>0)
      printFlush(paste("Removing samples:",
                       paste(rownames(dataExpr)[!gsg$goodSamples], collapse = ",")));
    # Remove the offending genes and samples from the data:
    dataExpr = dataExpr[gsg$goodSamples, gsg$goodGenes]
  }
  
  nGenes = ncol(dataExpr)
  nSamples = nrow(dataExpr)
  
  dim(dataExpr)
  
  ## [1]  134 2697
  
  head(dataExpr)[,1:8]
  #########################软阈值筛选#################
  ## 查看是否有离群样品
  sampleTree = hclust(dist(dataExpr), method = "average")
  
  pdf(paste0(save_path,"/sample_cluster.pdf"))
  plot1 <- plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")
  print(plot1)
  dev.off()
  ##############################
  powers = c(c(1:10), seq(from = 12, to=30, by=2))
  sft = pickSoftThreshold(dataExpr, powerVector=powers,
                          networkType=type, verbose=5)
  
  ## pickSoftThreshold: will use block size 2697.
  ##  pickSoftThreshold: calculating connectivity for given powers...
  ##    ..working on genes 1 through 2697 of 2697
  ##    Power SFT.R.sq  slope truncated.R.sq mean.k. median.k. max.k.
  ## 1      1   0.1370  0.825          0.412 587.000  5.95e+02  922.0
  ## 2      2   0.0416 -0.332          0.630 206.000  2.02e+02  443.0
  ## 3      3   0.2280 -0.747          0.920  91.500  8.43e+01  247.0
  ## 4      4   0.3910 -1.120          0.908  47.400  4.02e+01  154.0
  ## 5      5   0.7320 -1.230          0.958  27.400  2.14e+01  102.0
  ## 6      6   0.8810 -1.490          0.916  17.200  1.22e+01   83.7
  ## 7      7   0.8940 -1.640          0.869  11.600  7.29e+00   75.4
  ## 8      8   0.8620 -1.660          0.827   8.250  4.56e+00   69.2
  ## 9      9   0.8200 -1.600          0.810   6.160  2.97e+00   64.2
  ## 10    10   0.8390 -1.560          0.855   4.780  2.01e+00   60.1
  ## 11    12   0.8020 -1.410          0.866   3.160  9.61e-01   53.2
  ## 12    14   0.8470 -1.340          0.909   2.280  4.84e-01   47.7
  ## 13    16   0.8850 -1.250          0.932   1.750  2.64e-01   43.1
  ## 14    18   0.8830 -1.210          0.922   1.400  1.46e-01   39.1
  ## 15    20   0.9110 -1.180          0.926   1.150  8.35e-02   35.6
  ## 16    22   0.9160 -1.140          0.927   0.968  5.02e-02   32.6
  ## 17    24   0.9520 -1.120          0.961   0.828  2.89e-02   29.9
  ## 18    26   0.9520 -1.120          0.944   0.716  1.77e-02   27.5
  ## 19    28   0.9380 -1.120          0.922   0.626  1.08e-02   25.4
  ## 20    30   0.9620 -1.110          0.951   0.551  6.49e-03   23.5
  
  par(mfrow = c(1,2))
  cex1 = 0.9
  # 横轴是Soft threshold (power)，纵轴是无标度网络的评估参数，数值越高，
  # 网络越符合无标度特征 (non-scale)
  pdf(paste0(save_path,"/soft threshold power.pdf"))
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)",
       ylab="Scale Free Topology Model Fit,signed R^2",type="n",
       main = paste("Scale independence"))
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red")
  # 筛选标准。R-square=0.85
  abline(h=0.85,col="red")
  
  # Soft threshold与平均连通性
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
       main = paste("Mean connectivity"))
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers,
       cex=cex1, col="red")
  dev.off()
  
  #################
  power = sft$powerEstimate
  power
  
  ## [1] 6
  ##################################网络构建################
  ##一步法网络构建：One-step network construction and module detection##
  # power: 上一步计算的软阈值
  # maxBlockSize: 计算机能处理的最大模块的基因数量 (默认5000)；
  #  4G内存电脑可处理8000-10000个，16G内存电脑可以处理2万个，32G内存电脑可
  #  以处理3万个
  #  计算资源允许的情况下最好放在一个block里面。
  # corType: pearson or bicor
  # numericLabels: 返回数字而不是颜色作为模块的名字，后面可以再转换为颜色
  # saveTOMs：最耗费时间的计算，存储起来，供后续使用
  # mergeCutHeight: 合并模块的阈值，越大模块越少
  cor <- WGCNA::cor
  net = blockwiseModules(dataExpr, power = power, maxBlockSize = nGenes,
                         TOMType = type, minModuleSize = 50,
                         reassignThreshold = 0, 
                         numericLabels = TRUE, pamRespectsDendro = FALSE,
                         saveTOMs=TRUE, corType = corType,
                         maxPOutliers=maxPOutliers, loadTOMs=TRUE,
                         saveTOMFileBase = paste0('J:\\data\\annotation_final_edition\\expression\\R\\WGCNA\\', 'blockwiseModules', ".tom"),
                         verbose = 3)
  cor<-stats::cor
  
  ##  Calculating module eigengenes block-wise from all genes
  ##    Flagging genes and samples with too many missing values...
  ##     ..step 1
  ##  ..Working on block 1 .
  ##     TOM calculation: adjacency..
  ##     ..will use 47 parallel threads.
  ##      Fraction of slow calculations: 0.000000
  ##     ..connectivity..
  ##     ..matrix multiplication (system BLAS)..
  ##     ..normalization..
  ##     ..done.
  ##    ..saving TOM for block 1 into file WGCNA/LiverFemaleClean.txt.tom-block.1.RData
  ##  ....clustering..
  ##  ....detecting modules..
  ##  ....calculating module eigengenes..
  ##  ....checking kME in modules..
  ##      ..removing 3 genes from module 1 because their KME is too low.
  ##      ..removing 5 genes from module 12 because their KME is too low.
  ##      ..removing 1 genes from module 14 because their KME is too low.
  ##  ..merging modules that are too close..
  ##      mergeCloseModules: Merging modules whose distance is less than 0.25
  ##        Calculating new MEs...
  
  # 根据模块中基因数目的多少，降序排列，依次编号为 `1-最大模块数`。
  # **0 (grey)**表示**未**分入任何模块的基因。
  class(net)
  names(net)
  table(net$colors)
  save(net,file='J:\\data\\annotation_final_edition\\expression\\R\\WGCNA\\net.Rdata')
  ##
  ##   0   1   2   3   4   5   6   7   8   9  10  11  12  13
  ## 135 472 356 333 307 303 177 158 102  94  69  66  63  62
  
  #############################根据TOM值进行聚类#################
  ## 灰色的为**未分类**到模块的基因。
  # Convert labels to colors for plotting
  moduleLabels = net$colors
  moduleColors = labels2colors(moduleLabels)
  # Recalculate MEs with color labels
  MEs0 = moduleEigengenes(dataExpr, moduleColors)$eigengenes
  # 计算根据模块特征向量基因计算模块相异度：
  MEDiss = 1 - cor(MEs0);
  # Cluster module eigengenes
  METree = hclust(as.dist(MEDiss), method = "average");
  pdf(paste0(save_path, '/ME tree.pdf'))
  plot(METree,
       main = "Clustering of module eigengenes",
       xlab = "",
       sub = "")
  # 在聚类图中画出剪切线
  MEDissThres = 0.25
  abline(h = MEDissThres, col = "red")
  dev.off()
  ###############merge后###################
  #
  merge_modules = mergeCloseModules(dataExpr, moduleColors, cutHeight = MEDissThres, 
                                    verbose = 3)
  
  mergedColors = labels2colors(merge_modules$colors)
  merge_MEs0 = moduleEigengenes(dataExpr, mergedColors)$eigengenes
  # 计算根据模块特征向量基因计算模块相异度：
  merge_MEDiss = 1 - cor(merge_MEs0);
  # Cluster module eigengenes
  merge_METree = hclust(as.dist(merge_MEDiss), method = "average");
  
  pdf(paste0(save_path, '/merge_ME_tree.pdf'))
  plot(merge_METree,
       main = "Clustering of module eigengenes",
       xlab = "",
       sub = "")
  # 在聚类图中画出剪切线
  MEDissThres = 0.25
  abline(h = MEDissThres, col = "red")
  dev.off()
  #
  pdf(paste0(save_path, '/module_cluster.pdf'))
  merge_modules = mergeCloseModules(dataExpr, moduleColors, cutHeight = MEDissThres, 
                                    verbose = 3)
  dev.off()
  
  
  mergedMEs = merge_modules$newMEs;
  # Plot the dendrogram and the module colors underneath
  ##############################################
  ############################模块相关性#################
  # module eigengene, 可以绘制线图，作为每个模块的基因表达趋势的展示
  MEs = merge_modules$MEs
  MEs0 = moduleEigengenes(dataExpr, mergedColors)$eigengenes
  # 计算根据模块特征向量基因计算模块相异度：
  MEDiss = 1-cor(MEs0);
  # Cluster module eigengenes
  METree = hclust(as.dist(MEDiss), method = "average")
  # Plot the result
  pdf(paste0(save_path, '/EigengeneNetworks.pdf'))
  plotEigengeneNetworks(MEs0, 
                        "Eigengene adjacency heatmap", 
                        marHeatmap = c(3,4,2,2), 
                        plotDendrograms = FALSE, 
                        xLabelsAngle = 90)
  dev.off()
  ################################TOMplot######################
  # 如果采用分步计算，或设置的blocksize>=总基因数，直接load计算好的TOM结果
  # 否则需要再计算一遍，比较耗费时间
  # TOM = TOMsimilarityFromExpr(dataExpr, power=power, corType=corType, networkType=type)
  load(net$TOMFiles[1], verbose=T)
  
  ## Loading objects:
  ##   TOM
  
  TOM <- as.matrix(TOM)
  
  dissTOM = 1-TOMsimilarityFromExpr(dataExpr, power = power)
  
  nSelect = 400
  set.seed(10)
  select = sample(nGenes, size = nSelect)
  selectTOM = dissTOM[select, select]
  selectTree = hclust(as.dist(selectTOM), method = "average")
  selectColors = moduleColors[select]
  sizeGrWindow(9,9)
  # Transform dissTOM with a power to make moderately strong
  # connections more visible in the heatmap
  plotTOM = dissTOM^7
  # Set diagonal to NA for a nicer plot
  diag(plotTOM) = NA
  # Call the plot function
  plotDiss = selectTOM^7
  diag(plotDiss) = NA
  pdf(paste0(save_path, '/TOM_plot.pdf'))
  TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot",col=gplots::colorpanel(250,'red',"orange",'lemonchiffon'))
  dev.off()
  # 这一部分特别耗时，行列同时做层级聚类
  # TOMplot(plotTOM, net$dendrograms, moduleColors,main = "Network heatmap plot, all genes")
  #######################################导出成cytoscape格式############
  probes = colnames(dataExpr)
  dimnames(TOM) <- list(probes, probes)
  
  # Export the network into edge and node list files Cytoscape can read
  # threshold 默认为0.5, 可以根据自己的需要调整，也可以都导出后在
  # cytoscape中再调整
  cyt = exportNetworkToCytoscape(TOM,
                                 edgeFile = paste0(save_path, "//edges.txt", sep=""),
                                 nodeFile = paste0(save_path, "//nodes.txt", sep=""),
                                 weighted = TRUE, threshold = 0,
                                 nodeNames = probes, nodeAttr = moduleColors)
}
wgcna(
    file_path='J:\\data\\annotation_final_edition\\expression\\stringtie\\gene_tpm_matrix.csv',
    save_path='J:\\data\\annotation_final_edition\\expression\\R\\WGCNA_test'
)
