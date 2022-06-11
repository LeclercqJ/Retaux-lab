#Previous step (on Galaxy): RNAseq on cavefish (CF), surface fish (SF) and reciprocal hybrids
#mapping (Hisat2 on Amex2.0 NCBI ref in GFF); Variant calling (Freebayes), 
#SnpSift Filter (QUAL>30), retrieve gene with SnpEff eff (using own custom database with SnpEff build using GTF annotation),
#extract column of interest (SnpSift Extract Field), 
#Remove uncorrect gene id (Text manipulation, should not be an issue when using GTF for SnpEff database)

#Goal:
# 1/ Pre-processing (here)
# 2/ Identify type of regulation (cis, trans...) -> Cf other file
#Cf Landry et al, 2005; McManus et al, 2010; Bell et al, 2013; Bao et al; 2019 
# 3/ Retrieve if gene is DE or not to compare regulation -> other file
# 4/ Create a BED file for the interval restricted to fixed variants only (here)
#(used for cis/trans regulation in /2) 
#and rerun SnpEff annotation to get the impact of only the fixed variants (Galaxy.eu)


#1) Pre-processing:
#load packages required below
library(dplyr)
library(ggplot2)
#################################################################
#load file (SnpEff gene annotation 0bp using the gtf annotation file)
variant <- read.delim("Galaxy_Annotated_variants_GTF_FilterQ30.tabular",header = T, 
                      stringsAsFactors = F, sep = "\t", comment.char = "#")
str(variant)
colnames(variant) <- gsub("GEN\\.|_10hpf", "",colnames(variant))
colnames(variant) <- gsub("\\.+", ".",colnames(variant))
colnames(variant)
str(variant)

variant2 <- variant #671658
#1.1/Remove duplicated gene names/id in a same cell of the tab!
variant2$ANN.GENE <- sapply(variant2$ANN.GENE, function(x) paste(
  unique(unlist(strsplit(x,", "))), collapse = ","))
#1.2/Remove duplicated effect in a same cell of the tab!
variant2$ANN.EFFECT <- sapply(variant2$ANN.EFFECT, function(x) paste(
  unique(unlist(strsplit(x,", "))), collapse = ","))
head(variant2)
str(variant2)

#2/ remove variants in intergenic region (RNAseq: expect reads to be on genes, 
#plus if not, cannot be sure wich gene assign the read to)
variant_gene <- variant2[-which(variant2$ANN.EFFECT == "intergenic_region"),]
#n=503649

#3.1/ Find variant with multiple alternative (ma) and remove them
#if multiple alternative, they are separated with a "," in the ATL column
ma <-grep(",", variant_gene$ALT) #n=7679
variant_no_ma <- variant_gene[-ma,] #503649 - 7679 = 495970 (1.52%)

#3.2/ Remove variant with multiple annotated gene (mg)  
#either overlaping UTR or intron spanning a gene, 
#(i.e last exon of a gene is located after the next gene)
#if multiple annotated gene, they are separated with a "," in the ANN.GENE column
variant_no_mg <- variant_no_ma[-grep(",", variant_no_ma$ANN.GENE),] 

#3.3/ Look for variant without any associated gene and remove them if any
# can be the case when using GFF file with SnpEff (not an issue with GTF)
which(variant_no_mg$ANN.GENE == "")
which(is.na(variant_no_mg$ANN.GENE))
#"-" is used for overlaping/multiple genes, for trna and some genes like nkx2-1
#check if correct (if not remove incorrect annotation, can be the case with GFF)
grep("-", variant_no_mg$ANN.GENE)

#4/#put the AO value as numerical value instead of character 
#(character type due to several alt -> removed in 3.1)
variant_processed <- variant_no_mg
cols.num <- c("SF_rep1.AO", "SF_rep2.AO", "SF_rep3.AO", 
              "CF_rep1.AO", "CF_rep2.AO", "CF_rep3.AO",
              "HybSF_rep1.AO", "HybSF_rep2.AO", "HybSF_rep3.AO",
              "HybCF_rep1.AO", "HybCF_rep2.AO", "HybCF_rep3.AO")
variant_processed[cols.num] <- sapply(variant_processed[cols.num],as.numeric) 
str(variant_processed)

#5/ Keep only variants having reads in all 3 replicates of each condition (SF, CF & both hybrids)
#i.e no NA in any colon
variant_processed2 <- variant_processed[rowSums(is.na(variant_processed)) == 0,]
#rowSums computes column sums across rows that match the condition tested
#Then only the row with 0 NA in any replicates are selected 
#(i.e each variant are present in each replicates). 


#6/ Get fixed variant: i.e. homozygous REF (resp ALT) in all three SF replicates (resp CF rep) #29980
fixed_variant <- variant_processed2[rowSums(variant_processed2[, c("SF_rep1.GT", "SF_rep2.GT", "SF_rep3.GT")] == "0/0") == 3 
                                    & rowSums(variant_processed2[, c("CF_rep1.GT", "CF_rep2.GT", "CF_rep3.GT")] == "1/1") == 3,] 

#7/ merge variant per gene (sum of all fixed variant reads of a gene)
merged_gene <- fixed_variant %>% group_by(ANN.GENE) %>% #group by gene name #7999
  summarise(variant.count = n(), #get the number of variant for each gene
            SF_rep1.DP = sum(SF_rep1.DP, na.rm = T), SF_rep1.RO = sum(SF_rep1.RO, na.rm = T), SF_rep1.AO = sum(SF_rep1.AO, na.rm = T), #merge the read count values of the group (= sum of read / gene for each replicates)
            SF_rep2.DP = sum(SF_rep2.DP, na.rm = T), SF_rep2.RO = sum(SF_rep2.RO, na.rm = T), SF_rep2.AO = sum(SF_rep2.AO, na.rm = T), 
            SF_rep3.DP = sum(SF_rep3.DP, na.rm = T), SF_rep3.RO = sum(SF_rep3.RO, na.rm = T), SF_rep3.AO = sum(SF_rep3.AO, na.rm = T),
            
            CF_rep1.DP = sum(CF_rep1.DP, na.rm = T), CF_rep1.RO = sum(CF_rep1.RO, na.rm = T), CF_rep1.AO = sum(CF_rep1.AO, na.rm = T),
            CF_rep2.DP = sum(CF_rep2.DP, na.rm = T), CF_rep2.RO = sum(CF_rep2.RO, na.rm = T), CF_rep2.AO = sum(CF_rep2.AO, na.rm = T),
            CF_rep3.DP = sum(CF_rep3.DP, na.rm = T), CF_rep3.RO = sum(CF_rep3.RO, na.rm = T), CF_rep3.AO = sum(CF_rep3.AO, na.rm = T),
            
            HybSF_rep1.DP = sum(HybSF_rep1.DP, na.rm = T), HybSF_rep1.RO = sum(HybSF_rep1.RO, na.rm = T), HybSF_rep1.AO = sum(HybSF_rep1.AO, na.rm = T),
            HybSF_rep2.DP = sum(HybSF_rep2.DP, na.rm = T), HybSF_rep2.RO = sum(HybSF_rep2.RO, na.rm = T), HybSF_rep2.AO = sum(HybSF_rep2.AO, na.rm = T),
            HybSF_rep3.DP = sum(HybSF_rep3.DP, na.rm = T), HybSF_rep3.RO = sum(HybSF_rep3.RO, na.rm = T), HybSF_rep3.AO = sum(HybSF_rep3.AO, na.rm = T),
            
            HybCF_rep1.DP = sum(HybCF_rep1.DP, na.rm = T), HybCF_rep1.RO = sum(HybCF_rep1.RO, na.rm = T), HybCF_rep1.AO = sum(HybCF_rep1.AO, na.rm = T),
            HybCF_rep2.DP = sum(HybCF_rep2.DP, na.rm = T), HybCF_rep2.RO = sum(HybCF_rep2.RO, na.rm = T), HybCF_rep2.AO = sum(HybCF_rep2.AO, na.rm = T),
            HybCF_rep3.DP = sum(HybCF_rep3.DP, na.rm = T), HybCF_rep3.RO = sum(HybCF_rep3.RO, na.rm = T), HybCF_rep3.AO = sum(HybCF_rep3.AO, na.rm = T))

#9/ Filter gene by DP >10 in all replicates #5573
merged_gene_DP10 <- merged_gene[which(merged_gene$SF_rep1.DP >=10 & merged_gene$SF_rep2.DP >=10 & merged_gene$SF_rep3.DP >=10
                                      & merged_gene$CF_rep1.DP >=10 & merged_gene$CF_rep2.DP >=10 & merged_gene$CF_rep3.DP >=10
                                      & merged_gene$HybSF_rep1.DP >=10 & merged_gene$HybSF_rep2.DP >=10 & merged_gene$HybSF_rep3.DP >=10
                                      & merged_gene$HybCF_rep1.DP >=10 & merged_gene$HybCF_rep2.DP >=10 & merged_gene$HybCF_rep3.DP >=10),]



#Get some statistics about the number of variants per genes
table(merged_gene_DP10$variant.count)
mean(merged_gene_DP10$variant.count) #4.735331
median(merged_gene_DP10$variant.count) #4
quantile(merged_gene_DP10$variant.count, c(0.25, 0.5, 0.75, 0.9, 0.95, 0.99)) #25% 50% 75% 90% 95% 99% 
max(merged_gene_DP10$variant.count) #49                                       # 2   4   6   9  13  20.28 

p <- ggplot(merged_gene_DP10, aes(x=variant.count)) + 
  geom_histogram(binwidth=1, color="black", fill="white") + 
  geom_vline(aes(xintercept=median(variant.count)),
             color="red", linetype="dashed", size=1) +
  xlab("") + ggtitle("Number of variants per genes") +
  coord_cartesian(ylim = c(0,1050)) + theme_minimal(base_size = 18) +
  annotate(geom="text", x=7, y=820, label="median", color="red")
print(p) 
pdf(file = paste0("Variant_count_plots.pdf")) #save the plot as pdf
print(p)
dev.off()

#Save file (used for Cis/Trans regulation analysis)
write.table(merged_gene_DP10, "DESeq2_analysis_IFB-R-cluster/fixed_variant_gene_DP10.tab", 
            sep = "\t", col.names = T, row.names = F )

#################################################################
#Goal 2 & 3, script in another file (CisTrans analysis)
#################################################################
#Goal 4:
#Create a bed file restricted to fixed variants used for cis/trans analysis
#to use for SnpEff analysis in Galaxy (get impact of variant on proteins)

#A/ get vcf file pre annotation
vcf.preSnpEff <- read.delim("Galaxy82-[SnpSift_Filter_Q30_on_Freebayes_NO_markDup_(data_60)_].vcf",
                            header = F, stringsAsFactors = F, sep = "\t", comment.char = "#")
colnames(vcf.preSnpEff) <- c("CHROM",	"POS",	"ID",	"REF",	"ALT",	"QUAL",	"Filter",	"Info",	"Format",
                             "SF_10hpf_rep3",	"HybCF_10hpf_rep1",	"CF_10hpf_rep2",	"CF_10hpf_rep1",	
                             "SF_10hpf_rep2",	"SF_10hpf_rep1",	"CF_10hpf_rep3",	"HybCF_10hpf_rep3",	
                             "HybSF_10hpf_rep1",	"HybSF_10hpf_rep2",	"HybCF_10hpf_rep2",	
                             "HybSF_10hpf_rep3")


#B/ get fixed variants used (with sum DP >=10) #26390
fixed_variant2 <- fixed_variant[which(fixed_variant$ANN.GENE %in% merged_gene_DP10$ANN.GENE),]
#fixed_variant2 < fixed_variant because genes whith any DP<10 are not kept
fixed_variant2$ANN.GENE <- as.factor(fixed_variant2$ANN.GENE)
summary(fixed_variant2)
levels(fixed_variant2$ANN.GENE)
colnames(fixed_variant2)


#Add the genomic position of each variant (from the VCF)
vcf.preSnpEff.fixed <- semi_join(vcf.preSnpEff[,c(1,2,4:6)], fixed_variant2[,1:5])
vcf.preSnpEff.fixed2 <- merge(vcf.preSnpEff[,c(1,2,4:6)],fixed_variant2[,1:6])
#reorder
vcf.preSnpEff.fixed2b <- vcf.preSnpEff.fixed2[
  order(vcf.preSnpEff.fixed2$CHROM, vcf.preSnpEff.fixed2$POS)
  ,]

#Create a bed file (genomic coordinate) for each variant 
fix.variant.bed <- vcf.preSnpEff.fixed %>% 
  transmute(CHROM = CHROM, POS_START = POS-1, POS_END = POS)
#save and use in Glalaxy
write.table(fix.variant.bed, "Variant_fix-only_interval.bed", 
            col.names = F, row.names = F , sep = "\t")

##################################################################



