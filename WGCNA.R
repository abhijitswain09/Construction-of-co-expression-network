#load library
library(WGCNA)
library(DESeq2) # to perform log2 for normalizing the scale
library(GWENA) #for filtering the raw expression data
library(ppcor) #for doing partial correlation on data
#import the rpkm value of filtered genes
data0 = read.csv("count_gene_all_rpkm2.csv",sep = '\t' ,header=T)
#pick the samples according to the condition
data1=as.data.frame(t(data0[, c(8,14,15,5,6,7)]))
names(data1) = data0$Gene
rownames(data1) = names(data0)[c(8,14,15,5,6,7)]

#data1 = data0[,c(8,14,15,5,6,7)]
#filter out that genes which contains 0 expression data
w = t(data1)
df_new <- w[apply(w!=0, 1, all),]
data1 = t(df_new)
#log2(x+1) value of normalized data
datExpr = log2(data1+1)
#asign names to datExpr for further process
names(datExpr)=colnames(datExpr)
sampleTree = hclust(dist(datExpr), method = "average")
# plot sample tree
pdf(file = "N1_1-n-sampleClustering.pdf", width = 12, height = 9);
par(cex = 1.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="",
     cex.lab = 1.5,cex.axis = 1.5, cex.main = 2)
dev.off()
#Choose soft threshold parameter
powers = c(c(1:20), seq(from = 22, to=30, by=2))
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5) 
# Scale-free topology fit index as a function of the soft-thresholding power
pdf(file = "N1_2-n-sft.pdf", width = 9, height = 5);
par(mfrow = c(1,2));
cex1 = 0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
	xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
	main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red") 
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
	xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
	main = paste("Mean connectivity")) 
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()
#to calculate partial correlation
datExpr_parco = pcor(datExpr)
#extract only the partial correlation data
x= datExpr_parco$estimate
#assign the genes name
rownames(x) = colnames(datExpr)
colnames(x) = colnames(datExpr)
#x is not symmetric thats why we have to change it to symmetric data 
y = (x+t(x))/2
save(datExpr,datExpr_parco, x,y, file="N1_OUT.RData")
#create adjacency matrix from the correlated data
#   if checkSimilarity(y, min = -1, max = 1) its throw an error of  ### Error in checkAdjMat(similarity, min, max) : 
  ####some entries are not between-1and1
  ##then do this   @@@  y[y<0.000001] = 0 @@@@
adjacency = adjacency.fromSimilarity(y,
type = "unsigned",
power = 12)
#Turn adjacency into topological overlap
TOM = TOMsimilarity(adjacency)
dissTOM = 1-TOM
save(adjacency, TOM, dissTOM, file="N1_OUT1.RData")
#Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average")
# Plot the resulting clustering tree (dendrogram)
pdf(file = "Gene_clustering_on_TOM-based_dissimilarity.pdf", width = 12, height = 9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
labels = FALSE, hang = 0.04)
dev.off()
# Module identification using dynamic tree cut
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
deepSplit = 1, pamRespectsDendro = FALSE,
minClusterSize = 20)
#Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
#Plot the dendrogram and colors 
pdf(file = "Dynamic_tree_cut_N1.pdf", width = 8, height = 6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
dendroLabels = FALSE, hang = 0.03,
addGuide = TRUE, guideHang = 0.05,
main = "Gene dendrogram and module colors")
dev.off()
#Calculate eigengenes
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes
#Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs)
#Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average")
#plot cluster of module eigengenes
pdf(file = "clustering_Eigen_gene_N1.pdf", width = 7, height = 6)
plot(METree, main = "Clustering of module eigengenes",
xlab = "", sub = "")
dev.off()
#Plot the cut line into the dendrogram
 abline(h=0.05, col = "red")
#merge the modules
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = 0.05, verbose = 3)
#The merged module colors
moduleColors = merge$colors
#Eigengenes of the new merged modules
mergedMEs = merge$newMEs
#The merged module colors
mergedColors = merge$colors
#plot the merged gene tree
pdf(file = "merged_dynamic_tree_cut_N1.pdf", width = 9, height = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
c("Dynamic Tree Cut", "Merged dynamic"),
dendroLabels = FALSE, hang = 0.03,
addGuide = TRUE, guideHang = 0.05)
dev.off()
#to find how many genes are there in a single module
table(moduleColors)
#Construct numerical labels corresponding to the colors to label
colorOrder = c("grey", standardColors(150)) #if colorOrder having NA values then increase standardColors higher to solve the troubleshooting
moduleLabels = match(moduleColors, colorOrder)-1
#assign the gene names to moduleLabels
names(moduleLabels) = colnames(datExpr)
MEs = mergedMEs
save(MEs, moduleLabels, moduleColors, geneTree, file = "N1_OUT2.RData")
#to find which  genes are in the modules
names(moduleLabels)[moduleColors=="yellow4"]
#export the data modulize to visualize the graph
#Select modules
modules ="antiquewhite1"##if u want to select more than on module then use c("darkseagreen4","yellow4")
probes = names(datExpr)
#Select module probes
inModule = is.finite(match(moduleColors, modules))
modProbes = probes[inModule]
#Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule]
#export the details of nodes and edges to a txt file
cyt = exportNetworkToCytoscape(modTOM,
edgeFile = paste("CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
nodeFile = paste("CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
weighted = TRUE,
threshold = 0.3,
nodeNames = modProbes,
nodeAttr = moduleColors[inModule])
#choose hub gene
hub = chooseOneHubInEachModule(
datExpr,
moduleColors,
numGenes = 32,
omitColors = "grey",
power = 12,
type = "unsigned")
#to find intramodular connectivity
intramodularconnectivity = intramodularConnectivity(adjacency, moduleColors, scaleByMax = FALSE)
#to find network concept list(which contains summary,GS,conformity etc)
netcon= networkConcepts(adjacency,trait = NULL, networkType = "unsigned")
save(cyt,intramodularconnectivity,netcon,file = "N1_OUT3.RData")
#module eigengene significance
kME = signedKME(
datExpr,
MEs,
exprWeights = NULL,
MEWeights = NULL,
outputColumnName = "kME",
corFnc = "cor",
corOptions = "use = 'p'")


