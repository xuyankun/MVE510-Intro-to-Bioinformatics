x = read.table('counts_matrix.txt')
metadata = read.table('metadata.txt',sep='\t',header=TRUE)

rownames(x) # gene
colnames(x) # 80 samples
metadata$patient.id # same as the last command
colnames(x) == metadata$patient.id # all same

# ex1
count.nonzero = function(x){
  x.length = nrow(x)
  counts = vector(length=x.length,mode="double")
  for (i in 1:x.length) {
    num = sum(x[i,]>0)
    counts[i]= num
  }
  return(counts)
}

counts = count.nonzero(x)

x.filtered = x[counts>=20,]  # 38558 left
nrow(x.filtered)/nrow(x) > 0.25 # Ture

# ex2
CPM = function(x){
  x = x+1
  x.col = ncol(x)
  x.row = nrow(x)
  cpm = matrix(0,x.row,x.col)
  for (i in 1:x.col) {
    sum.col = sum(x[,i])
    
    for (j in 1:x.row) {
      cpm[j,i] = (x[j,i]/sum.col) * 1000000 
    }
  }
  return(cpm)
}

x.filtered.cpm = CPM(x.filtered)
x.filtered.logcpm = log(x.filtered.cpm)

# ex3

x.filtered.log = log(1 + x.filtered)

boxplot(x.filtered.logcpm,col = "bisque",ylab="log CPM" ,add = FALSE)
title("Boxplot of all genes in each sample with normalization")

boxplot(x.filtered.log, col = "bisque",ylab="log x.filtered" ,add = FALSE)
title("Boxplot of all genes in each sample without normalization")

plot(x.filtered.logcpm[,1], x.filtered.logcpm[,41], col="blue",
     main="Scatter plot for sample 1 and 41", 
     xlab = "Sample 1", ylab = "Sample 41")
par(new=TRUE)
plot(1:10, 1:10,"l",axes = FALSE,xlab = "", ylab = "")

# ex4

gene1.cpm = t(x.filtered.cpm[1,])
gene1.logcpm = t(x.filtered.logcpm[1,])

boxplot(gene1.cpm[1:40],gene1.cpm[41:80],names = c("not IBD","CD"))
title(" Boxplot of gene1 with CPM data")
boxplot(gene1.logcpm[1:40], gene1.logcpm[41:80],names = c("not IBD","CD"))
title(" Boxplot of gene1 with logCPM data")

diagnosis = as.factor(metadata[,"diagnosis"])
diagnosis = relevel(diagnosis, ref = "Not IBD")
sex = as.factor(metadata[,"Sex"])
age = metadata[,"age.at.diagnosis"]

fit1 = lm(x.filtered.logcpm[1,]~diagnosis)
fit2 = lm(x.filtered.logcpm[1,]~age+sex+diagnosis)

# ex5

gene.length = nrow(x.filtered.logcpm)
fit1.pval = vector(length = gene.length,mode = "double")
fit1.coeff = vector(length = gene.length,mode = "double")
fit2.pval = vector(length = gene.length,mode = "double")
fit2.coeff = vector(length = gene.length,mode = "double")

for (i in 1:gene.length) {
  fit1.i = lm(x.filtered.logcpm[i,]~diagnosis)
  fit2.i = lm(x.filtered.logcpm[i,]~age+sex+diagnosis)
  fit1.pval[i]= summary(fit1.i)$coefficients["diagnosisCD","Pr(>|t|)"]
  fit1.coeff[i]= summary(fit1.i)$coefficients["diagnosisCD","Estimate"]
  fit2.pval[i]= summary(fit2.i)$coefficients["diagnosisCD","Pr(>|t|)"]
  fit2.coeff[i]= summary(fit2.i)$coefficients["diagnosisCD","Estimate"]
}

fit1.padjust = data.frame(cbind(p.adjust(fit1.pval, method = "fdr"), fit1.coeff))
fit2.padjust = data.frame(cbind(p.adjust(fit2.pval, method = "fdr"), fit2.coeff))

rownames(fit1.padjust) = c(rownames(x.filtered))
rownames(fit2.padjust) = c(rownames(x.filtered))

rownames(x.filtered.logcpm) = c(rownames(x.filtered)) 
rownames(x.filtered.cpm) = c(rownames(x.filtered)) 

colnames(fit1.padjust) = c("p-adjust","coefficients")
colnames(fit2.padjust) = c("p-adjust","coefficients")

cutoff = 0.05

fit1.diff.coeff= fit1.padjust[fit1.padjust[,1] < cutoff,2]
fit2.diff.coeff= fit2.padjust[fit2.padjust[,1] < cutoff,2]

num.fit1.diff = length(fit1.diff.coeff)
num.fit2.diff = length(fit2.diff.coeff)

fit1.coeff.up.num = sum(fit1.diff.coeff>0)
fit1.coeff.down.num = num.fit1.diff - fit1.coeff.up.num
fit2.coeff.up.num = sum(fit2.diff.coeff>0)
fit2.coeff.down.num = num.fit2.diff - sum(fit2.diff.coeff>0)

fit1.most.gene = rownames(fit1.padjust)[fit1.padjust[,1] == min(fit1.padjust[,1])]
fit2.most.gene = rownames(fit2.padjust)[fit2.padjust[,1] == min(fit2.padjust[,1])]
# most significant gene: "ENSG00000185499"
coeff.most.gene.fit1 = fit1.padjust["ENSG00000185499",2]
coeff.most.gene.fit2 = fit2.padjust["ENSG00000185499",2]
# almost same

# effect size (log fold-change)
log.fold.change = log(sum(x.filtered.logcpm["ENSG00000185499",41:80]) /
  sum(x.filtered.logcpm["ENSG00000185499",1:40]))
# 0.73
log.fold.change = log(sum(x.filtered.cpm["ENSG00000185499",41:80])/
  sum(x.filtered.cpm["ENSG00000185499",1:40]))
# 2.49

fit1.padjust.order=order(fit1.padjust[,1], decreasing=FALSE)
fit1.padjust.sorted = fit1.padjust[fit1.padjust.order,]
fit2.padjust.order=order(fit2.padjust[,1], decreasing=FALSE)
fit2.padjust.sorted = fit2.padjust[fit2.padjust.order,]

rownames(fit1.padjust.sorted)[1:6]
fit1.padjust.sorted[1:6,2] > 0
# gene
# "ENSG00000185499" "ENSG00000203747" "ENSG00000183010" "ENSG00000150337" "ENSG00000162747" "ENSG00000182240"
# coefficients (up or down - regulate)
# TRUE TRUE TRUE TRUE TRUE TRUE
rownames(fit2.padjust.sorted)[1:6]
fit2.padjust.sorted[1:6,2] > 0
# "ENSG00000185499" "ENSG00000203747" "ENSG00000183010" "ENSG00000162747" "ENSG00000140274" "ENSG00000150337"(old)
# "ENSG00000185499" "ENSG00000203747" "ENSG00000183010" "ENSG00000150337" "ENSG00000162747" "ENSG00000182240"
# TRUE TRUE TRUE TRUE TRUE TRUE


fit2.age.pval = vector(length = gene.length,mode = "double")
fit2.sex.pval = vector(length = gene.length,mode = "double")

for (i in 1:gene.length) {
  fit2.i = lm(x.filtered.logcpm[i,]~age+sex+diagnosis)
  fit2.age.pval[i]= summary(fit2.i)$coefficients["age","Pr(>|t|)"]
  fit2.sex.pval[i] = summary(fit2.i)$coefficients["sexMale","Pr(>|t|)"]
}

fit2.age.padjust = p.adjust(fit2.age.pval)
# all pvalue adjust is 1. none is related to age

fit2.sex.padjust = p.adjust(fit2.sex.pval)
rownames(x.filtered.logcpm)[fit2.sex.padjust==min(fit2.sex.padjust)]
# "ENSG00000067048" most significant for gender


# ex6

xsig = x.filtered.logcpm[rownames(fit1.padjust.sorted)[1:100],]
xsig = as.matrix(xsig)

mypalette = brewer.pal(11,"RdYlBu")
morecols = colorRampPalette(mypalette)
mycols=rev(morecols(255))
column.cols=c("purple","orange")[metadata$diagnosis]
pdf("top100sigGenesHeatmap.pdf",height=9,width=6)
heatmap.2(xsig,trace='none',col=mycols,main='The 100 most significant genes',ColSideColors=column.cols)
dev.off()
# orange is not IBD

# ex7

pca=prcomp(t(x.filtered.logcpm))

PC1 = pca$x[,1]
PC2 = pca$x[,2]
PC3 = pca$x[,3]

plot(PC1,PC2,main = "Scatter plot of PC1 and PC2", xlab = "PC1",
     ylab = "PC2", col=column.cols)

plot(PC1,PC3,main = "Scatter plot of PC1 and PC3", xlab = "PC1",
     ylab = "PC3", col=column.cols)

plot(PC2,PC3,main = "Scatter plot of PC2 and PC3", xlab = "PC2",
     ylab = "PC3", col=column.cols)

# based on gender

sex.cols=c("red","blue")[metadata$Sex]
# red for female, blue for male

plot(PC2,PC3,main = "Scatter plot of PC2 and PC3", xlab = "PC2",
     ylab = "PC3", col=sex.cols)

plot(PC1,PC3,main = "Scatter plot of PC1 and PC3", xlab = "PC1",
     ylab = "PC3", col=sex.cols)

plot(PC1,PC2,main = "Scatter plot of PC1 and PC2", xlab = "PC1",
     ylab = "PC2", col=sex.cols)
