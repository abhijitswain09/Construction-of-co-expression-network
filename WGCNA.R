# Load only packages used by this analysis
library(WGCNA)
library(ppcor)
args <- commandArgs(trailingOnly = TRUE)
input_file <- if (length(args) >= 1) args[1] else "count_gene_all_rpkm2.csv"
sample_columns <- if (length(args) >= 2) as.integer(strsplit(args[2], ",")[[1]]) else c(8,14,15,5,6,7)
data0 <- read.delim(input_file, header = TRUE, check.names = FALSE)
if (anyNA(sample_columns) || any(sample_columns < 1) || any(sample_columns > ncol(data0)))
  stop("sample_columns must be comma-separated valid 1-based matrix column numbers")
if (length(sample_columns) < 4) stop("WGCNA requires at least four selected samples")
if (!"Gene" %in% names(data0)) stop("Input matrix must contain a Gene column")
data1 <- as.data.frame(t(data0[, sample_columns, drop = FALSE]))
names(data1) <- data0$Gene
rownames(data1) <- names(data0)[sample_columns]
w <- t(data1)
df_new <- w[apply(w != 0, 1, all), , drop = FALSE]
data1 <- t(df_new)
datExpr <- log2(data1 + 1)
names(datExpr) <- colnames(datExpr)
sampleTree <- hclust(dist(datExpr), method = "average")
pdf(file = "N1_1-n-sampleClustering.pdf", width = 12, height = 9)
par(cex = 1.6); par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")
dev.off()
powers <- c(1:20, seq(from=22, to=30, by=2))
sft <- pickSoftThreshold(datExpr, powerVector=powers, verbose=5)
pdf(file="N1_2-n-sft.pdf", width=9, height=5)
par(mfrow=c(1,2)); cex1 <- 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab="Soft Threshold (power)", ylab="Scale Free Topology Model Fit,signed R^2", type="n", main="Scale independence")
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers, cex=cex1, col="red")
abline(h=0.90, col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab="Soft Threshold (power)", ylab="Mean Connectivity", type="n", main="Mean connectivity")
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1, col="red")
dev.off()
datExpr_parco <- pcor(datExpr)
x <- datExpr_parco$estimate
rownames(x) <- colnames(datExpr); colnames(x) <- colnames(datExpr)
y <- (x + t(x))/2
save(datExpr, datExpr_parco, x, y, file="N1_OUT.RData")
adjacency <- adjacency.fromSimilarity(y, type="unsigned", power=12)
TOM <- TOMsimilarity(adjacency); dissTOM <- 1-TOM
save(adjacency, TOM, dissTOM, file="N1_OUT1.RData")
geneTree <- hclust(as.dist(dissTOM), method="average")
pdf(file="Gene_clustering_on_TOM-based_dissimilarity.pdf", width=12, height=9)
plot(geneTree, xlab="", sub="", main="Gene clustering on TOM-based dissimilarity", labels=FALSE, hang=0.04)
dev.off()
dynamicMods <- cutreeDynamic(dendro=geneTree, distM=dissTOM, deepSplit=1, pamRespectsDendro=FALSE, minClusterSize=20)
dynamicColors <- labels2colors(dynamicMods)
pdf(file="Dynamic_tree_cut_N1.pdf", width=8, height=6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", dendroLabels=FALSE, hang=0.03, addGuide=TRUE, guideHang=0.05, main="Gene dendrogram and module colors")
dev.off()
MEList <- moduleEigengenes(datExpr, colors=dynamicColors)
MEs <- MEList$eigengenes
MEDiss <- 1-cor(MEs); METree <- hclust(as.dist(MEDiss), method="average")
pdf(file="clustering_Eigen_gene_N1.pdf", width=7, height=6)
plot(METree, main="Clustering of module eigengenes", xlab="", sub="")
dev.off()
merged <- mergeCloseModules(datExpr, dynamicColors, cutHeight=0.05, verbose=3)
moduleColors <- merged$colors; mergedMEs <- merged$newMEs
pdf(file="merged_dynamic_tree_cut_N1.pdf", width=9, height=6)
plotDendroAndColors(geneTree, cbind(dynamicColors, moduleColors), c("Dynamic Tree Cut","Merged dynamic"), dendroLabels=FALSE, hang=0.03, addGuide=TRUE, guideHang=0.05)
dev.off()
colorOrder <- c("grey", standardColors(150))
moduleLabels <- match(moduleColors, colorOrder)-1
names(moduleLabels) <- colnames(datExpr)
MEs <- mergedMEs
save(MEs, moduleLabels, moduleColors, geneTree, file="N1_OUT2.RData")
modules <- "antiquewhite1"
probes <- names(datExpr)
inModule <- is.finite(match(moduleColors, modules))
modProbes <- probes[inModule]; modTOM <- TOM[inModule, inModule]
cyt <- exportNetworkToCytoscape(modTOM, edgeFile=paste0("CytoscapeInput-edges-",paste(modules,collapse="-"),".txt"), nodeFile=paste0("CytoscapeInput-nodes-",paste(modules,collapse="-"),".txt"), weighted=TRUE, threshold=0.3, nodeNames=modProbes, nodeAttr=moduleColors[inModule])
hub <- chooseOneHubInEachModule(datExpr, moduleColors, numGenes=32, omitColors="grey", power=12, type="unsigned")
intramodularconnectivity <- intramodularConnectivity(adjacency, moduleColors, scaleByMax=FALSE)
netcon <- networkConcepts(adjacency, trait=NULL, networkType="unsigned")
save(cyt, intramodularconnectivity, netcon, file="N1_OUT3.RData")
kME <- signedKME(datExpr, MEs, exprWeights=NULL, MEWeights=NULL, outputColumnName="kME", corFnc="cor", corOptions="use = 'p'")
