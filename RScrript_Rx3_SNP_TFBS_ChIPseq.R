#Rx3 is a gene with change at the cis-regulatory level between cave and surface fish
#We identified 8 SNPs close to the Rx3 gene, that may be responsible for the change in gene expression, 
#for instance by creating or deleting a transcription factor binding site (TFBS)
#Goal: scan DNA sequences for putative TFBS and compare between surface and cave position with fixed SNP (at rx3 level) 
#Use JASPAR database for TFBS and TFBSTools for the analysis

#help from
#https://bioconductor.org/packages/release/bioc/vignettes/TFBSTools/inst/doc/TFBSTools.html
#http://bioconductor.org/packages/release/bioc/manuals/TFBSTools/man/TFBSTools.pdf

#if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
#BiocManager::install("TFBSTools")
#BiocManager::install("JASPAR2020")

library(TFBSTools) #v1.30
library(JASPAR2020) #version 2020 for R v4.1
library(Biostrings) #v2.62.0 #if load FASTA from file
library(dplyr)

opts <- list()
opts[["tax_group"]] <- "vertebrates"
opts[["type"]] <- "ChIP-seq" #ChIP-seq data are more accurate than SELEX according to Hu et al., 2010
PFMatrixList <- getMatrixSet(JASPAR2020, opts)
PFMatrixList

opts2 <- list()
opts2[["tax_group"]] <- "vertebrates"
testPFM <- getMatrixSet(JASPAR2020, opts2)

#transform PFM into PWM (as)
pwmList <- toPWM(PFMatrixList, pseudocounts=0.8) #pseudocounts = 0.8 as recommended in the vignette

#sequences to test ( +/- 15 bp)
subject <- c("GGAATAAGCATAAAAAATAAGCAAGGAATAA", "GGAAGTTGACAAATTCTTCCAATTAAGTTCA", 
             "CTGTGGTTAGAGAATCAGGAATTATGATATA", "TTCAGAGACTGTGACGTTGCACTGTGACAGA",
             "ATGAGCCGAACTGTCCCTGAACACAGAAACT", "CTTGTACAATGCAATAAAAAATGCATATCC", 
             "GAATATAAATATTTTGAATATGATTTCTTTA", "GGACGGTGTCTCTCTCTCTTCCTCTTCCTCTCGCTC",
             "GGAATAAGCATAAAAGATAAGCAAGGAATAA", "GGAAGTTGACAAATTATTCCAATTAAGTTCA", 
             "CTGTGGTTAGAGAATGAGGAATCATGATATA", "TTCAGAGACTGTGACTTTGCACTGTGACAGA",
             "ATGAGCCGAACTGTCACTGAACACAGAAACT", "CTTGTACAATGCAATAAAAAAATGCATATCC", 
             "GAATATAAATATTTTTAATATGATTTCTTTA", "GGACGGTGTCTCTCTCTCTTCCTCTCGCTC") 
seqname <- c("SF_rx3_SNP1", "SF_rx3_SNP2", "SF_rx3_SNP3", "SF_rx3_SNP4", 
             "SF_rx3_SNP5", "SF_rx3_SNP6", "SF_rx3_SNP7", "SF_rx3_SNP8",
             "CF_rx3_SNP1", "CF_rx3_SNP2", "CF_rx3_SNP3", "CF_rx3_SNP4", 
             "CF_rx3_SNP5", "CF_rx3_SNP6", "CF_rx3_SNP7", "CF_rx3_SNP8")

#1/Scan sequences for predicted TFBS
#initialize vectors to store output results
sitesetList <- vector(mode = "list", length = length(subject))
gff3_results <- vector(mode = "list", length = length(subject))
rel_score <- vector(mode = "list", length = length(subject))
pval <- vector(mode = "list", length = length(subject))

for (i in 1:length(subject)) {
  #Scan TFBS
  sitesetList[[i]] <- searchSeq(pwmList, subject[i], seqname = seqname[i],
                           min.score="90%", strand="*")
  cat(paste(seqname[i], "Scanning - Done\n"))
  
  #Extract result in GFF
  gff3_results[[i]] <- writeGFF3(sitesetList[[i]])
  cat(paste(seqname[i], "GFF3 - Done\n"))
  
  #Get relative score
  rel_score[[i]] <- relScore(sitesetList[[i]])
  cat(paste(seqname[i], "Scores - Done\n"))
  
  #Get p-value
  pval[[i]] <-pvalues(sitesetList[[i]], type="TFMPvalue")
  cat(paste(seqname[i], "p-values - Done\n"))
  
  Sys.sleep(1)
}

#gff3_res80 <- gff3_results
#pval80 <- pval
#rel_score80 <- rel_score

#gff3_res85 <- gff3_results
#pval85 <- pval
#rel_score85 <- rel_score

#gff3_res90 <- gff3_results
#pval90 <- pval
#rel_score90 <- rel_score
#gff3_res90_tf <- gff3_results[lenths(gff3_res90) != 0]
#pval90_tf <- pval[lenths(gff3_res90) != 0]
#rel_score90_tf <- rel_score[lenths(gff3_res90) != 0]

#variable for the script below
gff3_res <- gff3_res80
pval <- pval80
rel_score <- rel_score80 

#for 2 sequences (score 80%, ChIP-seq data) -> duplicated entry (same TF with different class attribute), remove it
#gff3_res[[1]] <- gff3_res[[1]][-2,]
#gff3_res[[11]] <- gff3_res[[11]][-17,]

#2/ Compare SF and CF sequences for morph-specific predicted TFBS
#initialize vectors to store output results
gff3_results_SF <- vector(mode = "list", length = length(gff3_res)/2)
gff3_results_CF <- vector(mode = "list", length = length(gff3_res)/2)
SF_only_TF <- vector(mode = "list", length = length(gff3_res)/2)
CF_only_TF <- vector(mode = "list", length = length(gff3_res)/2)  
SF_only_TFinfo <- vector(mode = "list", length = length(gff3_res)/2)
CF_only_TFinfo <- vector(mode = "list", length = length(gff3_res)/2)
SF_sep <- vector(mode = "list", length = length(gff3_res)/2)
CF_sep <- vector(mode = "list", length = length(gff3_res)/2)

for (j in 1:(length(gff3_res)/2)) {
  #add scores and p-values to the gff3
  gff3_results_SF[[j]] <- cbind(TFBS.id = rownames(gff3_res[[j]]),
                                gff3_res[[j]],
                                rel.score = unlist(rel_score[[j]]),
                                p.value = unlist(pval[[j]]))
  gff3_results_CF[[j]] <- cbind(TFBS.id = rownames(gff3_res[[j+length(gff3_res)/2]]),
                                gff3_res[[j+length(gff3_res)/2]],
                                rel.score = unlist(rel_score[[j+length(gff3_res)/2]]),
                                p.value = unlist(pval[[j+length(gff3_res)/2]]))
  
  #extract TF names from the gff3
  SF_attributes <- unlist(gff3_results_SF[[j]]["attributes"])
  SF_TF <- unique(gsub("TF=", "", 
                       as.vector(stri_match(SF_attributes,
                                            regex =  "TF=[^\\;]+"))))
  
  CF_attributes <- unlist(gff3_results_CF[[j]]["attributes"])
  CF_TF <- unique(gsub("TF=", "",  
                       as.vector(stri_match(CF_attributes, 
                                            regex =  "TF=[^\\;]+"))))
  
  #Compare TF names between CF and SF (get only those that are different)
  SF_only_TF[[j]] <- SF_TF[-which(SF_TF %in% CF_TF)]
  CF_only_TF[[j]] <- CF_TF[-which(CF_TF %in% SF_TF)]
  
  #display in the console morph-specific tfbs
  cat(paste0("SF-only TFBS (", unique(gff3_results_SF[[j]]["seqname"]), ") :\n"))
  print(SF_only_TF[[j]])
  cat(paste0("CF-only TFBS (", unique(gff3_results_CF[[j]]["seqname"]), ") :\n"))
  print(CF_only_TF[[j]])
  
  #extract morph-specific values in each gff3
  f1.SF <- function(tf){
    #get the line index for each morph-specific TFBS
    id <- gsub("\\(", "\\\\(", tf)
    id <- gsub("\\)", "\\\\)", id)
    index <- grep(paste0(id,";"), unlist(gff3_results_SF[[j]]["attributes"]))
    
    TFBS_info <- gff3_results_SF[[j]][index,]
    return(TFBS_info)
  }
  SF_only_TFinfo[[j]] <- as.data.frame(t(sapply(SF_only_TF[[j]], f1.SF)))
  f1.CF <- function(tf){
    #get the line index for each morph-specific TFBS
    id <- gsub("\\(", "\\\\(", tf)
    id <- gsub("\\)", "\\\\)", id)
    index <- grep(paste0(id,";"), unlist(gff3_results_CF[[j]]["attributes"]))
    
    TFBS_info <- gff3_results_CF[[j]][index,]
    return(TFBS_info)
  }
  CF_only_TFinfo[[j]] <- as.data.frame(t(sapply(CF_only_TF[[j]], f1.CF)))
  
  #Separate multiple occurrence of one TFBS into individual rows
  if(ncol(SF_only_TFinfo[[j]]) !=0) {
    SF_sep[[j]] <- as.data.frame(matrix(unlist(SF_only_TFinfo[[j]]), ncol = 12))
    }
  if(ncol(CF_only_TFinfo[[j]]) !=0) {
    CF_sep[[j]] <- as.data.frame(matrix(unlist(CF_only_TFinfo[[j]]), ncol = 12))
    }

}
#Merge result by SNP in one single table (per morph)
SF_singledf <- bind_rows(SF_sep[which(lengths(SF_sep) != 0)])
colnames(SF_singledf) <- names(SF_only_TFinfo[[4]])
SF_singledf$TFBS.name <- gsub("TF=", "", as.vector(stri_match(SF_singledf$attributes,
                                                              regex =  "TF=[^\\;]+")))
CF_singledf <- bind_rows(CF_sep[which(lengths(CF_sep) != 0)])
colnames(CF_singledf) <- names(CF_only_TFinfo[[1]])
CF_singledf$TFBS.name <- gsub("TF=", "", as.vector(stri_match(CF_singledf$attributes,
                                                              regex =  "TF=[^\\;]+")))


#for score >= 0.8
SF_only_TF #no predicted specific TFBS for SNP6
unique(SF_singledf$TFBS.name)
CF_only_TF #no predicted specific TFBS for SNP8
unique(CF_singledf$TFBS.name)

unique(SF_singledf$TFBS.name) #23 unique in SF
unique(CF_singledf$TFBS.name) #33 unique in CF
intersect(SF_singledf$TFBS.name, CF_singledf$TFBS.name) #5 shared (but on different sequence/location)
setdiff(SF_singledf$TFBS.name, CF_singledf$TFBS.name) #5 shared
unique(unlist(c(SF_singledf$TFBS.name, CF_singledf$TFBS.name))) #51 unique TFBS 
#33+23 = 56 (5 shared CF/SF). Total 55 unique predicted TFBS in both CF and SF
#so most of them are morph-specific !!!

#Number is the threshold of min score used
#SF80 <- SF_singledf
#CF80 <- CF_singledf



SF85 <- SF80[which(SF80$rel.score > 0.85),]
CF85 <- CF80[which(CF80$rel.score > 0.85),]

SF90 <- gff3_results_SF[which(gff3_results_SF$rel.score > 0.9),] 
#threshold of score>=90%, only 3 TFBS found in SF-only (none for CF-only) 
#-> not very interesting (3 TF involved in heamolymphatic differenciation 
#+ position on SNP8, deletion in CT repeat so can be functional in CF too (score > 0.8))



#FIMO tools (available at https://meme-suite.org/meme/tools/fimo)
#Cite: Charles E. Grant, Timothy L. Bailey, and William Stafford Noble, "FIMO: Scanning for occurrences of a given motif", Bioinformatics, 27(7):1017-1018, 2011.
#Get the PFM of each TFBS
#https://rdrr.io/bioc/TFBSTools/src/R/JASPAR.R
makeFlatFileDir(JASPAR2020)

#load the PFM in a single list
setwd("/shared/ifbstor1/projects/fb_am/variants/DESeq2_variant_analysis/cis_only_analysis/Rx3_TFBS/FlatFileDir")
temp = list.files(pattern="*.pfm")
myfiles = lapply(temp, read.delim, header = FALSE)
names(myfiles) <- gsub("\\.pfm", "",list.files(pattern = ".pfm"))

#select only PFM used above with TFBSTools
myPFM <- myfiles[names(PFMatrixList)]
names(myPFM) <- paste(names(PFMatrixList), name(PFMatrixList))

#save all PFM (+ name/id) in a single txt output
sink("myPFMlist_ChIPseq.txt")
print(myPFM)
sink()

#do the following command to modify the file for FIMO (Unix)
#sed -r -e 's/\$/\>/'  -e 's/`//g' -e '/V[0-9]/d' -e 's/^[1234] //' <myPFMlist_ChIPseq.txt  >myPFMlist_ChIPseq_FIMO.txt

#analyse sequences using FIMO (link above)

#load results
fimo <- read.table("fimo_ChIPseq_pval10-4.tsv", header = T)

fimo.res <- fimo %>% group_by(motif_id, motif_alt_id, gsub("[^0-9]", "", fimo$sequence_name)) %>%
  summarise(n = n()) %>% filter(n == 1) %>% select(motif_id, motif_id)

fimo.unique <- fimo %>% filter(motif_id %in% fimo.res$motif_id)

similar <- fimo.unique[fimo.unique$motif_id %in% c(SF80$TFBS.id, CF80$TFBS.id),]
#only 3 similarly predicted unique TFBS but when look closer, might not be correct because miss an important nucleotide
#PKNOX1, GTT C AGTGACAGTTC (-) wrong C -> should be G
#NFE2L1, ACTGTGACT TT GCACT (+) wrong TT -> should be CA
#Nr5a2, GTGTTCA G GGACAGT (-) wrong G -> should be A


