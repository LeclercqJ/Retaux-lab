#THE DETA HAVE BEEN PREPROCESSED IN ANOTHER SCRIPT

#Goal : Run a DESeq2 analysis on the read counts of fixed variant/genes
#and identify Cis vs Trans regulation at play (if any)

#load required packages
library(DESeq2)
library(dplyr)
library(ggplot2)
library(stringr)

#####################################################
######### DESeq2 computation of A and B #############
#Load the matrix of read count for REF and ALT of each replicate 
#(to use for DESeq2 analysis) from FreeBayes proccessed file.
gene_fixed_variant <- read.delim("fixed_variant_gene_DP10.tab", sep = "\t", header = T)

#Keep only REF and ALT count successively for each condition 
#(for parents, REF = SF and ALT = CF total counts (DP))
Fixed_variant_DESeq <- gene_fixed_variant[,c("SF_rep1.DP", "SF_rep2.DP", "SF_rep3.DP", 
                                             "CF_rep1.DP", "CF_rep2.DP", "CF_rep3.DP",
                                             "HybSF_rep1.RO", "HybSF_rep2.RO", "HybSF_rep3.RO", 
                                             "HybSF_rep1.AO", "HybSF_rep2.AO", "HybSF_rep3.AO",
                                             "HybCF_rep1.RO", "HybCF_rep2.RO", "HybCF_rep3.RO",
                                             "HybCF_rep1.AO", "HybCF_rep2.AO", "HybCF_rep3.AO" )]
rownames(Fixed_variant_DESeq) <- gene_fixed_variant$ANN.GENE #rename row with gene names

#Create the sample description file of all samples
cd <- c("Parent", "HybSF", "HybCF")
c <- length(cd) #number of condition : Parent and both hybrids
s <- ncol(Fixed_variant_DESeq) /2 #number of sample (/2 because alt and ref count for each)
r <- s/c #number of replicates
samps <- data.frame(allele=factor(rep(c("ref","alt"), each = r),levels=c("ref","alt")), #Precise ref and alt allele for each replicates
                    condition=factor(rep(cd,each= 2*s/c))) #Specify the condition
rownames(samps) <- colnames(Fixed_variant_DESeq)
print(samps) #check that the file is correct!
rm(c, cd, r, s) #remove these variables to make the environment clearer

#create a function to run the DEseq2 analysis with 3 parameters:
#parent and/or hybrids name and the design for DESeq2
DEanalysis <- function(parent = NULL, hybrid = NULL, design){ 
  #extract from the whole sample info file the subset of sample of interest
  cdt <- ifelse(is.null(parent) == F & is.null(hybrid) ==F, 
                paste(parent, hybrid, sep = "|"), 
                paste0(parent,hybrid))
  sample <- grep(paste0("^(", cdt  ,").*"),rownames(samps))
  
  #Extract the matrix count column matching the sample names of interest 
  Count.Matrix <- Fixed_variant_DESeq[,sample]
  #print information of the selected matrix to check if correct
  print("Count matrix info:")
  print(str(Count.Matrix))
  
  #Extract the sample info matching the sample names of interest
  Samp.design <- samps[sample,] 
  #print information of the selected info file to check if correct
  print("Sample design:")
  print(Samp.design) 
  
  #Create the DESeq Data Set (dds) according to chosen samples and design
  dds <- DESeqDataSetFromMatrix(countData = Count.Matrix , 
                                colData = Samp.design, 
                                design = design)
  
  #Set as reference both condition and allele count (Works for Hybrids vs Parent and/or ALT vs REF only)
  dds$condition <- relevel(dds$condition, ref = "Parent")
  dds$allele <- relevel(dds$allele, ref = "ref")
  
  #run DESeq2 analysis
  dds <- DESeq(dds, fitType = "parametric") 
  
  print("Results:")
  print(resultsNames(dds))
  #output the result outside the function
  if(is.null(parent) == F & is.null(hybrid) == T) {
    #if compare only parent alt vs ref ratio, directly output the result
    return(results(dds, name="allele_alt_vs_ref"))}
  
  else if(is.null(parent) == T & is.null(hybrid) == F) {
    #if compare only alt vs ref ratio in hybrid only, directly output the result
    return(results(dds, name="allele_alt_vs_ref"))}
  
  else if(is.null(parent) == F & is.null(hybrid) == F) {
    #if compare the ratio in both condition, output both results in a list
    return(list(parent_allelic_ratio = results(dds, name="conditionParent.allelealt"), 
                hybrid_allelic_ratio = results(dds, name=paste0("condition", hybrid,".allelealt"))))} 
}

#Run analysis
parent <- DEanalysis(parent = "SF|CF", hybrid = NULL, design= ~allele)
hybSF <- DEanalysis(parent = NULL, hybrid = "HybSF", design= ~allele)
hybCF <- DEanalysis(parent = NULL, hybrid = "HybCF", design= ~allele)


#Save DESeq2 results (fold change and stats)
write.table(parent, "results/DESeq_res_parent.csv", sep = ";", dec = ",", col.names = T, row.names = T)
write.table(hybSF, "results/DESeq_res_hybSF.csv", sep = ";", dec = ",", col.names = T, row.names = T)
write.table(hybCF, "results/DESeq_res_hybCF.csv", sep = ";", dec = ",", col.names = T, row.names = T)

#####################################################

#tHE Following code of Cis/Trans Analysis IS adapted from Bao et al, 2019, Nature Communication: 
#https://doi.org/10.1038/s41467-019-13386-w
#Chunk of code from the github repository (from supplemental): 
#https://github.com/Wendellab/CisTransRegulation/commit/3f8f36adc454251cbd57084f38eeb7f9f2ddf476#
#"cistran.deseq2.r" file; lines ~549-636
#I mainly changed sample names to match ours (i.e. Maxxa and TX2094 to SF and CF), changed plot colors 
#and added a threshold of (adjusted) p-value to the classicCisTrans function to make it easier to change

#####################################################
############Prepare Cis/Trans Analysis ##############
### comparing A-B= 0 is tricky, both are log2FoldChange and its standard error lfcse
# maybe I can compare with t test
#### T test from means and standard errors
# m1, m2: the sample means
# s1, s2: the sample standard deviations
# n1, n2: the same sizes
# se1, se2: the sample standard errors
# se1 <- s1/sqrt(n)
# m0: the null value for the difference in means to be tested for. Default is 0.
# equal.variance: whether or not to assume equal variance. Default is FALSE.
t.test2 <- function(m1,m2,se1,se2,n1,n2,m0=0,equal.variance=FALSE)
{
  if( equal.variance==FALSE )
  {
    # se <- sqrt( (s1^2/n1) + (s2^2/n2) )
    se <- sqrt( (se1^2) + (se2^2) )
    # welch-satterthwaite df
    # df <- ( (s1^2/n1 + s2^2/n2)^2 )/( (s1^2/n1)^2/(n1-1) + (s2^2/n2)^2/(n2-1) )
    df <- ( (se1^2 + se2^2)^2 )/( (se1^2)^2/(n1-1) + (se2^2)^2/(n2-1) )
  } else
  {
    # pooled standard deviation, scaled by the sample sizes
    # se <- sqrt( (1/n1 + 1/n2) * ((n1-1)*s1^2 + (n2-1)*s2^2)/(n1+n2-2) )
    se <- sqrt( (1/n1 + 1/n2) * ((n1-1)*s1^2 + (n2-1)*s2^2)/(n1+n2-2) )
    df <- n1+n2-2
  }
  t <- (m1-m2-m0)/se
  dat <- c(m1-m2, se, t, 2*pt(-abs(t),df))
  names(dat) <- c("Difference of means", "Std Error", "t", "p-value")
  return(dat)
}
x1 = rnorm(3)
x2 = rnorm(3)
# you'll find this output agrees with that of t.test when you input x1,x2
t.test(x1,x2)
t.test2( mean(x1),  mean(x2), sd(x1)/sqrt(3), sd(x2)/sqrt(3), 3,3)
rm(x1,x2)

### Make a categorization table for each F1 condition:
criteria<- as.data.frame(rbind(c("A!=0;B!=0;A=B", "1.Cis only"),
                               c("A!=0;B=0;A!=B", "2.Trans only"),
                               c("A!=0;B!=0;A!=B", "Cis+Trans"),
                               c("A=0;B!=0;A!=B", "5.Compensatory"),
                               c("A=0;B=0;A=B", "6.Conserved") ))
names(criteria) <- c("class","category")
print(criteria)

#Create the function to compute the type of regulation from DESeq2 results
classCisTrans<-function(A.res, B.res, A.n, B.n, log2fc.threshold=0, padj.threshold=0.05)
{
  # A = log2(CF/SF), cis + trans
  A <- A.res[,c("log2FoldChange", "lfcSE", "padj")]
  names(A) <- c("A", "A.SE", "A.padj")
  # B = log2(F1_Hyb_AO/F1_Hyb_RO), cis
  B <- B.res[,c("log2FoldChange", "lfcSE", "padj")]
  names(B) <- c("B", "B.SE", "B.padj")
  
  A<-A[rownames(B),]
  table <- cbind(A,B)
  #A-B, trans
  table$AminusB <- table$A - table$B
  table$AminusB.pvalue <- apply(table,1,function(x) t.test2(m1=x[1],m2=x[4],
                                                            se1=x[2], se2=x[5], 
                                                            n1=A.n, n2=B.n)["p-value"])
  
  table$cisNtrans <- ifelse(abs(table$A) >= log2fc.threshold 
                            & table$A.padj < padj.threshold 
                            & !is.na(table$A.padj), "A!=0", "A=0")
  table$cis <- ifelse(abs(table$B) >= log2fc.threshold 
                      & table$B.padj < padj.threshold 
                      & !is.na(table$B.padj), "B!=0", "B=0")
  table$trans <- ifelse(abs(table$AminusB) >= log2fc.threshold 
                        & table$AminusB.pvalue < padj.threshold, "A!=B", "A=B")
  
  table$class <- paste(table$cisNtrans,table$cis,table$trans,sep=";")
  
  table$category <- as.character(criteria$category[match(table$class,criteria$class)])
  table$category[is.na(table$category)] <- "7.Ambiguous"
  table$category[ table$category=="Cis+Trans" & table$B*table$AminusB >0 ] <- "3.Cis+Trans: enhancing"
  table$category[ table$category=="Cis+Trans" & table$B*table$AminusB <0 ] <- "4.Cis+Trans: compensating"
  
  return(as.data.frame(table))
}

#Function to plot the categories from the different comparison
#Modified a bit from Bao et al because some categories are missing: 
#need to keep the same color code even with fewer categories
plotCisTrans<-function(table, plotTitle="") 
{
  #Choose the list of color:
  colors <- c("red","blue","purple","cyan1","green2","grey25","grey")
  #Create an index from the categories number:
  colors.index <- as.numeric(gsub("\\..*","",levels(as.factor(table$category))))
  #Keep only colors that match the index position:
  colors.plot <- colors[colors.index]
  
  #Make the plot
  p <- ggplot( table, aes(x=A, y=B, color=category)) + geom_point(alpha=0.8) +  
    xlab("Cis + Trans") + ylab("Cis") + ggtitle(plotTitle) + 
    scale_color_manual(values=colors.plot) + geom_vline(xintercept=0) + 
    geom_hline(yintercept=0) +  geom_abline(intercept = 0, slope = 1) + 
    theme_bw()
  
  print(p)
  
  tiff(file = paste0("results/",plotTitle,".tiff"), 
       units="in", width=6, height=6, res=400)
  print(p)
  dev.off()
}

#####################################################
###New from me (end of R script from Bao et al.)#####
#####################################################
############ Cis/Trans SNP vs DEG ###################

#The CisTrans calculcation is Trans = A - B 
#With A the fold change in parents and B the allelic ratio in hybrids
#B comes from number of reads at a single position (SNP) and so is calculate A
#But doing this can introduce discrepancy when comparing with differential expression done using the full length RNA.
#In order to be very confident in the attribution of the category, I also did the CisTrans Analysis using the RNA full length to calculate A (used for classic DEG analysis)
#Then I selected only genes with categories consistent between the two methods of calculation of A (SNP and DEG)


#Load file from DEG analysis at 10hpf (table v6, subtab of 10hpf only)
genes.10hpf <- read.delim("gene_10hpf.csv", header = T, sep = ";", dec = ",") 
genes.10hpf$Previous_NCBI_Name <- gsub("sept-0", "sept", genes.10hpf$Previous_NCBI_Name)

#Load result of DESEQ2 (CF vs SF)
DESeq.DEG <- read.delim("10hpf_DESeq2_results_CFvsSF.tabular", header = F, sep = "\t", 
                        col.names = c("Gene.name", "baseMean", "log2FoldChange", 
                                      "lfcSE", "stat", "pvalue", "padj"))
head(DESeq.DEG)


#### Method for calculation of A : SNP ####
#Threshold of FC and padj set to match the one used in the DEG analysis (used later for A)
cat.HybSF.SNP <- classCisTrans(A.res = parent, B.res = hybSF, A.n = 3, B.n=3, 
                               log2fc.threshold=0.58, padj.threshold=0.01)
cat.HybCF.SNP <- classCisTrans(A.res = parent, B.res = hybCF, A.n = 3, B.n=3, 
                               log2fc.threshold=0.58, padj.threshold=0.01)
table(cat.HybSF.SNP$category)
table(cat.HybCF.SNP$category)

cat.HybSF.SNP <- cbind(rownames(cat.HybSF.SNP),cat.HybSF.SNP)
colnames(cat.HybSF.SNP)[1] <- "Name"
cat.HybCF.SNP <- cbind(rownames(cat.HybCF.SNP),cat.HybCF.SNP)
colnames(cat.HybCF.SNP)[1] <- "Name"

#save result of CisTrans analysis
write.table(cat.HybSF.SNP, "results/catCisTrans_A.snp_HybSF_FC1.5_FDR1percent.csv", 
            sep = ";", dec = ",", col.names = T, row.names = F)
write.table(cat.HybCF.SNP, "results/catCisTrans_A.snp_HybCF_FC1.5_FDR1percent.csv", 
            sep = ";", dec = ",", col.names = T, row.names = F)

#Plot results
plotCisTrans(cat.HybSF.SNP, "HybSF 10hpf SNP")
plotCisTrans(cat.HybCF.SNP, "HybCF 10hpf SNP")


#### Method for calculation of A : DEG ####
#Merge Deseq2 from DEG with gene with fixed variant to get DEG with fixed variant.
DESeq.DEG.A <- merge(gene_fixed_variant, DESeq.DEG, by.x = "ANN.GENE", by.y = "Gene.name", all.x = T)
DESeq.DEG.A[which(is.na(DESeq.DEG.A$baseMean)), c(1,2,39:44)] #Search if non attributed variant -> 0
rownames(DESeq.DEG.A) <- DESeq.DEG.A$ANN.GENE
#Get regulatory categories using A from DEG and B from individual hybrid ratio
#Use the same threshold as DEG (lfc >0.58 & padj<0.01)
#HybSF
cat.HybSF.DEG <- classCisTrans(A.res = DESeq.DEG.A, B.res = hybSF, A.n = 3, B.n=3, 
                               log2fc.threshold=0.58, padj.threshold=0.01)

cat.HybCF.DEG <- classCisTrans(A.res = DESeq.DEG.A, B.res = hybCF, A.n = 3, B.n=3, 
                               log2fc.threshold=0.58,padj.threshold=0.01)


table(cat.HybSF.DEG$category)
table(cat.HybCF.DEG$category)

cat.HybCF.DEG <- cbind(rownames(cat.HybCF.DEG),cat.HybCF.DEG)
colnames(cat.HybCF.DEG)[1] <- "Name"

cat.HybSF.DEG <- cbind(rownames(cat.HybSF.DEG),cat.HybSF.DEG)
colnames(cat.HybSF.DEG)[1] <- "Name"

#save tables for DEG A
write.table(cat.HybSF.DEG, "results/catCisTrans_A.deg_HybSF_FC1.5_FDR1percent.csv", 
            sep = ";", dec = ",", col.names = T, row.names = F)
write.table(cat.HybCF.DEG, "results/catCisTrans_A.deg_HybCF_FC1.5_FDR1percent.csv", 
            sep = ";", dec = ",", col.names = T, row.names = F)

#Plot results
plotCisTrans(cat.HybSF.DEG, "HybSF 10hpf DEG")
plotCisTrans(cat.HybCF.DEG, "HybCF 10hpf DEG")


#### Some Comparisons ####
#Merge info with CisTrans results

Merge <- merge(genes.10hpf[,c(1:8,24,25)], cat.HybSF.SNP[,c(1,14)], 
               by.x = "Previous_NCBI_Name", by.y = "Name", all.y = T)
colnames(Merge)[1] <- "Name"

#Merge results of both hybrids
to_merge <- list(Merge, cat.HybCF.SNP[,c(1,14)], cat.HybSF.DEG[,c(1,14)], cat.HybCF.DEG[,c(1,14)])
Cat.merge <- Reduce(function(...) left_join(..., by= "Name"), to_merge)
colnames(Cat.merge)[11:14] <- c("Cat.HybSF.snp", "Cat.HybCF.snp", 
                                "Cat.HybSF.deg", "Cat.HybCF.deg")

table(Cat.merge$Cat.HybSF.snp)
table(Cat.merge$Cat.HybCF.snp)

table(Cat.merge$Cat.HybSF.deg)
table(Cat.merge$Cat.HybCF.deg)

write.table(Cat.merge, "results/catCisTrans_all-merge_FC1.5_FDR1percent.csv", 
            sep = ";", dec = ",", col.names = T, row.names = F)

Cat.merge <- read.delim("results/catCisTrans_all-merge_FC1.5_FDR1percent.csv", 
                        sep = ";", dec = ",", header = T)

#DEG among genes with fixed variants
deg <- Cat.merge %>% filter(DEG_10hpf == "UP" | DEG_10hpf == "DOWN")
#900 DEG (/!\ not necessarily same annotation)
  
#### Common regulation ####
common.snp.cis <- which(Cat.merge$Cat.HybSF.snp == Cat.merge$Cat.HybCF.snp
                        & Cat.merge$Cat.HybSF.snp == "1.Cis only")
#160 genes with common cis-regulation in both hybrids (A from SNP)
#64.3% (160/249) of HybSF cis-reg genes  & 39.7% (160/403) of HybCF cis-reg genes 

common.deg.cis <- which(Cat.merge$Cat.HybSF.deg == Cat.merge$Cat.HybCF.deg
                        & Cat.merge$Cat.HybSF.deg == "1.Cis only")
#118 genes with common cis-regulation in both hybrids (A from DEG)
#63.1% (118/187) of HybSF cis-reg genes  & 44.7% (118/264) of HybCF cis-reg genes 


#Get common cis-only genes in both hybrids for SNP or DEG
#/!\ Make sure categories have the same thresholds 
#between use of A from SNP or DEG -> padj<0.01 & |lfc|->0.58
true.cis <- Cat.merge %>% filter(Cat.HybSF.snp == "1.Cis only",
                                       Cat.HybCF.snp == "1.Cis only",
                                       Cat.HybSF.deg == "1.Cis only",
                                       Cat.HybCF.deg == "1.Cis only")
#108 "true" cis-only-regulated genes
write.table(true.cis, "cis-only_genes.csv", sep = ";", row.names = F, col.names = T)

#Get genes with common category in both hybrids for SNP or DEG
true.cat <- Cat.merge %>% filter(Cat.HybSF.snp == Cat.HybCF.snp,
                                 Cat.HybSF.snp == Cat.HybSF.deg,
                                 Cat.HybSF.snp == Cat.HybCF.deg)
table(true.cat$Cat.HybSF.snp)
write.table(true.cat, "true_categories_genes.csv", sep = ";", row.names = F, col.names = T)
#4211 genes with same annotation between hybrids and consistent with DEG analysis
#72.9% of consistent annotation btw hyb  and btw SNP/DEG

#Venn diagramm for cis-only and trans-only genes between hybrids
#Create a list with onl the gene belonging to one category (cis-only and trans-only respectively)
#Cis-only list
test.cis.venn <- list()
test.cis.venn["HybSF.deg"] <- Cat.merge %>% 
  filter(Cat.HybSF.deg == "1.Cis only") %>% select(Name)
test.cis.venn["HybSF.snp"] <- Cat.merge %>% 
  filter(Cat.HybSF.snp == "1.Cis only") %>% select(Name)
test.cis.venn["HybCF.deg"] <- Cat.merge %>% 
  filter(Cat.HybCF.deg == "1.Cis only") %>% select(Name)
test.cis.venn["HybCF.snp"] <- Cat.merge %>% 
  filter(Cat.HybCF.snp == "1.Cis only") %>% select(Name)
#Trans-only list
test.trans.venn <- list()
test.trans.venn["HybSF.deg"] <- Cat.merge %>% 
  filter(Cat.HybSF.deg == "2.Trans only") %>% select(Name)
test.trans.venn["HybSF.snp"] <- Cat.merge %>% 
  filter(Cat.HybSF.snp == "2.Trans only") %>% select(Name)
test.trans.venn["HybCF.deg"] <- Cat.merge %>% 
  filter(Cat.HybCF.deg == "2.Trans only") %>% select(Name)
test.trans.venn["HybCF.snp"] <- Cat.merge %>% 
  filter(Cat.HybCF.snp == "2.Trans only") %>% select(Name)

#Plot the Venn diagram and save it
#For saving in tiff
tiff("Venn_cis_allmethod.tiff", units="in", width=6, height=6, res=300)
print(cis)
dev.off()

library(ggvenn)

cis <- ggvenn(
  test.cis.venn, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4
) + ggtitle("Cis-only genes")
pdf("Venn_cis_allmethod.pdf")
print(cis)
dev.off()

cis.2 <- ggvenn(
  test.cis.venn, columns = c("HybSF.deg", "HybCF.deg"),
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4
) + ggtitle("Cis-only genes DEG method")
pdf("Venn_cis_DEGmethod.pdf")
print(cis.2)
dev.off()

cis.3 <- ggvenn(
  test.cis.venn, columns = c("HybSF.snp", "HybCF.snp"),
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4
) + ggtitle("Cis-only genes SNP method")
pdf("Venn_cis_SNPmethod.pdf")
print(cis.3)
dev.off()

trans <-ggvenn(
  test.trans.venn, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4
) + ggtitle("Trans-only genes DEG method")
pdf("Venn_trans_Allmethod.pdf")
print(trans)
dev.off()

trans.2 <- ggvenn(
  test.trans.venn, columns = c("HybSF.deg", "HybCF.deg"),
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4
) + ggtitle("Trans-only genes DEG method")
pdf("Venn_trans_DEGmethod.pdf")
print(trans.2)
dev.off()

trans.3 <- ggvenn(
  test.trans.venn, columns = c("HybSF.snp", "HybCF.snp"),
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4
) + ggtitle("Trans-only genes SNP method")
pdf("Venn_trans_SNPmethod.pdf")
print(trans.3)
dev.off()


