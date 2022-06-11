#Goal: Identify which missense variants cause a change of aminoacid class

#load packages
library(tidyr)
library(dplyr)
library(stringr)
library(ggplot2)

#load files
#amino acid (aa) table
aa <- read.delim("Amino-acid_class.csv", header = T, sep = ";",
                 col.names = c("Name","Abbreviation","Letter","Class","Polarity","Charge"))

#list of all fixed missense variants (cf pre-processing file)
missense_variants <- read.delim("Galaxy-SnpSift_Extract_Fields_on_MISSENSE_variants_mRNA_ids.tabular", 
                                header = T, sep = "\t")
colnames(missense_variants) <- gsub(pattern = "\\.+", replacement = ".", 
                                     x = colnames(missense_variants))

#remove undesired string
missense_variants$ANN.HGVS_P <- gsub(pattern = "(p\\.)|(, )", replacement = "", 
                                      x = missense_variants$ANN.HGVS_P)
head(missense_variants)

#Seprate reference aa (surface fish, SF) and alternate aa (cavefish, CF) in two columns
missense_variants2 <- missense_variants %>%
  separate(ANN.HGVS_P, into = c("aa_SF", "aa_CF"), 
           sep = "[0-9]+", remove = F)
missense_variants2$aa_change <- gsub(pattern = "[0-9]+", replacement = ">",
                                      x = missense_variants$ANN.HGVS_P)
missense_variants2$aa_change_position <- gsub("[A-Za-z]+", "", missense_variants2$ANN.HGVS_P)

####
#Handle the complex cases (two aa changed side-by-side, e.g AlaVal>IlePro)
#Separate them in individual line per change (e.g Ala>Ile & Val>Pro)
#1/ Get only the line with multiple aa change
mmp <- missense_variants2[!missense_variants2$aa_SF %in% aa$Abbreviation,]
SF_aa_match <- str_match_all(mmp[,"aa_SF"], "[A-Z](.)(.)")
CF_aa_match <- str_match_all(mmp[,"aa_CF"], "[A-Z](.)(.)")

#Get the number of aa in each case (mostly 2 but 3 aa once)
n <- sapply(SF_aa_match, nrow)

#Duplicate the table per row by number of aa changed
idx.dup <- rep(1:nrow(mmp), n)
mmp.dup <- mmp[idx.dup,]

#Get the new number for the second (and third) aa
new_pos <- function(position, n){
  new <- seq(from = position,
             to = position + n-1)
  return(new)
}

#Change colon accordingly
mmp.dup$aa_SF <- grep("[A-Z]", unlist(SF_aa_match), value = T)
mmp.dup$aa_CF <- grep("[A-Z]", unlist(CF_aa_match), value = T)
mmp.dup$aa_change <- paste(mmp.dup$aa_SF, mmp.dup$aa_CF, sep = ">")
mmp.dup$aa_change_position <- unlist(mapply(new_pos, 
                                             as.numeric(mmp$aa_change_position), 
                                             n = n))

#Create the new table from the previous ones
missense_variants2.1 <- missense_variants2[missense_variants2$aa_SF %in% aa$Abbreviation,] 
missense_variants2.2 <- rbind(missense_variants2.1, mmp.dup)
#reorder by genomic coordinate
missense_variants2.2 <- missense_variants2.2[with(missense_variants2.2, order(CHROM, POS)),]

#remove the line with no change (il17ra E335E)
missense_variants2.3 <- missense_variants2.2[-which(missense_variants2.2$aa_SF == missense_variants2.2$aa_CF),]
#and change names with multiple gene names/ transcript id (n=1: LOC111196219/XM_022685429.1)
error_gene <- grep(" ", missense_variants2.3$ANN.GENE)
error_transcript <- grep(" ", missense_variants2.3$ANN.FEATUREID)
missense_variants2.3[error_gene, "ANN.GENE"] <- gsub(",.*", "", missense_variants2.3[error_gene, "ANN.GENE"])
missense_variants2.3[error_transcript, "ANN.FEATUREID"] <- gsub(",.*", "", missense_variants2.3[error_transcript, "ANN.FEATUREID"])

####

missense_variants3 <- missense_variants2.3
SF.colon <- c("Class.SF", "Polarity.SF","Charge.SF")
CF.colon <- c("Class.CF", "Polarity.CF","Charge.CF")
Attribute <- c("Class", "Polarity","Charge")
missense_variants3$Change <- NA
#Get the aa info (class, polarity, charge) for both SF & CF variant
for (i in 1:nrow(missense_variants3)) {
  
  if(missense_variants3[i,"aa_SF"] %in% aa$Abbreviation == F) {next}
  
  else{
    
    missense_variants3[i,SF.colon] <- aa[which(aa$Abbreviation 
                                            == missense_variants3[i,"aa_SF"]),
                                         c(4:6)]
    
    missense_variants3[i,CF.colon] <- aa[which(aa$Abbreviation 
                                                  == missense_variants3[i,"aa_CF"]),
                                            c(4:6)]
    for (j in 1:length(SF.colon)) {
      
      ifelse(missense_variants3[i,SF.colon[j]] != missense_variants3[i,CF.colon[j]],
             ifelse(is.na(missense_variants3[i,"Change"]),
                    missense_variants3[i,"Change"] <- Attribute[j],
                    missense_variants3[i,"Change"] <- paste(missense_variants3[i,"Change"], Attribute[j], sep = ", ")
                    ),
             NA)
    }
    
    ifelse(is.na(missense_variants3[i,"Change"]),
                 missense_variants3[i,"Change"] <- "No Change",
                 NA)
  }
}

#some stats
missense_variants4 <- missense_variants3
missense_variants4$Change <- as.factor(missense_variants4$Change)
missense_variants4$aa_change <- as.factor(missense_variants4$aa_change)
missense_variants4$aa_SF <- as.factor(missense_variants4$aa_SF)
missense_variants4$aa_CF <- as.factor(missense_variants4$aa_CF)

str(missense_variants4)  
summary(missense_variants4)
table(missense_variants4$aa_change)

barplot(table(missense_variants4$aa_change))
barplot

#Results/Conclusion:
#5252 missense variants (-> 5287 when separate MNP/complex case): 
#Change (all) = 4006 (75.76%)
#No change = 1281 (24.24%)

#plot nb of variants per genes (for missens variants only)
variant_perGene <- missense_variants4 %>% group_by(ANN.GENE) %>% #group by gene name #2747
  dplyr::summarise(variant.count = n())

p <- ggplot(variant_perGene, aes(x=variant.count)) + 
  geom_histogram(binwidth=1, color="black", fill="white") + 
  geom_vline(aes(xintercept=median(variant.count)),
             color="red", linetype="dashed", size=1) +
  xlab("") + ggtitle("Number of missense variants per genes") +
  coord_cartesian(ylim = c(0,1600)) + theme_minimal(base_size = 18) +
  annotate(geom="text", x=2.3, y=1050, label="median", color="red")
pdf(file = paste0("Missense_variant_count_plots.pdf")) #save the plot as pdf
print(p)
dev.off()

#Some statistics about the number of variants per genes
table(variant_perGene$variant.count)
mean(variant_perGene$variant.count) #1.924645
median(variant_perGene$variant.count) #1
quantile(variant_perGene$variant.count, c(0.25, 0.5, 0.75, 0.9, 0.95, 0.99)) #25% 50% 75% 90% 95% 99% 
max(variant_perGene$variant.count) #20                                       # 1   1   2   4  5  8 

#Add change using abbreviation i.e P106L instead of Pro106Leu (arbitrary example)
missense_variants4$aa_SF_abbrev <- aa$Letter[match(missense_variants4$aa_SF, aa$Abbreviation)]
missense_variants4$aa_CF_abbrev <- aa$Letter[match(missense_variants4$aa_CF, aa$Abbreviation)]
missense_variants4$aa_change_abbrev <- paste0(missense_variants4$aa_SF_abbrev, 
                                              missense_variants4$aa_change_position, 
                                              missense_variants4$aa_CF_abbrev)

write.table(missense_variants4, "Fixed_missense_variants_Astyanax_10hpf_fulltab.tab", 
            sep = "\t", dec = ".", row.names = F, col.names = T)


########################################################################
############## Get the protein sequence for each variants ##############
########################################################################

#/!\ Look only at missense variants here 
#(possibly other variants type in some of those genes, i.e synonymous, nonsens...)

#Use rentrez package to extract fasta sequences of those protein 
#Use the RNA id from SnpEff -> 1st isoform in NCBI

#Load the required package to access NCBI database
library(rentrez)

#Get the RNA ids for each gene (in SnpEff - from GFF3 Astyanax_mexicanus2.0 annotation file - NCBI)
GeneName <- as.character(unique(missense_variants4$ANN.FEATUREID))
searchID <- function(gene){
#From the NCBI/entrez RNA id, get the Refseq RNA id
id_search <- entrez_search(db="nucleotide", 
                           term= paste("Astyanax mexicanus", gene))
return(id_search$ids)
}
gene_ids <- sapply(GeneName, searchID) #can be long if too many ids
gene_ids

#Create the function to convert RNA Refseq id into petpride id and get the fasta sequence
counter <- 0 #Use a counter to see which gene have wrong ids
getFasta <- function(id){
  linked_seq_ids <- entrez_link(dbfrom="nucleotide", 
                                id = id, 
                                db="protein")
  
  fasta_recs <- entrez_fetch(db="protein", 
                             id=linked_seq_ids$links$nuccore_protein, 
                             rettype="fasta")
  counter <<- counter+1 #to trace back error but keep order
  return(fasta_recs)
}

test <- sapply(gene_ids[59], getFasta)

#try small chunk of ids. Adjust when errors
f1.59 <- sapply(gene_ids[1:59], getFasta)
#f59,179,223,761, 1000, 1010,1075,1855,2367,2485,2610,2724 
#pb no correspondence RNA id - protein id (but exist both in NCBI, manual check)
f61.178 <- sapply(gene_ids[61:178], getFasta) 
f180.222 <- sapply(gene_ids[180:222], getFasta)
f224.300 <- sapply(gene_ids[224:300], getFasta)
f301.500 <- sapply(gene_ids[301:500], getFasta)
f501.760 <- sapply(gene_ids[501:760], getFasta)
f762.999 <- sapply(gene_ids[762:999], getFasta)
f1001.1009 <- sapply(gene_ids[1001:1009], getFasta)
f1011.1074 <- sapply(gene_ids[1011:1074], getFasta)
f1076.1300 <- sapply(gene_ids[1076:1300], getFasta)
f1301.1500 <- sapply(gene_ids[1301:1500], getFasta)
f1501.1854 <- sapply(gene_ids[1501:1854], getFasta)
f1856.2300 <- sapply(gene_ids[1856:2300], getFasta)
f2301.2366 <- sapply(gene_ids[2301:2366], getFasta)
f2368.2484 <- sapply(gene_ids[2368:2484], getFasta)
f2386.2609 <- sapply(gene_ids[2486:2609], getFasta)
f2611.2723 <- sapply(gene_ids[2611:2723], getFasta)
f2725.2747 <- sapply(gene_ids[2725:2747], getFasta)

#to get correct protein id #59,179...: nucc_id +1
#f60,179,223,761, 1000, 1010,1075,1855,2367,2485,2610,2724 
getFasta2 <- function(id){
  prot_id_manual <- as.character(as.numeric(id) + 1)
  
  fasta_recs <- entrez_fetch(db="protein", 
                             id=prot_id_manual, 
                             rettype="fasta")
  
  return(fasta_recs)
}
f60 <- sapply(gene_ids[60], getFasta2)
f179 <- sapply(gene_ids[179], getFasta2)
f223 <- sapply(gene_ids[223], getFasta2)
f761 <- sapply(gene_ids[761], getFasta2)
f1000 <- sapply(gene_ids[1000], getFasta2)
f1010 <- sapply(gene_ids[1010], getFasta2)
f1075 <- sapply(gene_ids[1075], getFasta2)
f1855 <- sapply(gene_ids[1855], getFasta2)
f2367 <- sapply(gene_ids[2367], getFasta2)
f2485 <- sapply(gene_ids[2485], getFasta2)
f2610 <- sapply(gene_ids[2610], getFasta2)
f2724 <- sapply(gene_ids[2724], getFasta2)

#Check that the sum of all retrieve data (i.e. fasta) is correct (n=2747) 
length(c(f1.59, f60, f61.178, f179, f180.222, f223, f224.300, f301.500, f501.760, f761,
         f762.999, f1000, f1001.1009, f1010, f1011.1074, f1075, f1076.1300, f1301.1500, 
         f1501.1854, f1855, f1856.2300, f2301.2366, f2367, f2368.2484, f2485, f2486.2609,
         f2610, f2611.2723, f2724, f2725.2747))

#Merge all fasta records in one file and save it
fasta.all <- c(f1.59, f60, f61.178, f179, f180.222, f223, f224.300, f301.500, f501.760, f761,
               f762.999, f1000, f1001.1009, f1010, f1011.1074, f1075, f1076.1300, f1301.1500, 
               f1501.1854, f1855, f1856.2300, f2301.2366, f2367, f2368.2484, f2485, f2486.2609,
               f2610, f2611.2723, f2724, f2725.2747)

write(fasta.all, file="all_change_peptides.fasta")

########################################################################
############# Change header to meet Mutpred specification ##############
########################################################################
#load fasta file for all genes with at least a missense variant
library(seqinr)
fasta <- read.fasta(file="all_change_peptides.fasta", 
                    seqtype = "AA", as.string = TRUE, set.attributes = FALSE)
names(fasta)
#get duplicated ids
which(duplicated(names(fasta))) #none -> ok

#Get genes not starting with a Methionine
grep("^[^M]", fasta) #(n=9) -> keep them 

####
#Table  with the correspondence between transcript id, gene id and mutation
aa_change <- missense_variants4[, c("ANN.FEATUREID", "ANN.GENE", "aa_change_abbrev")]

#Merge all missense variants in one entry per gene
library(data.table)
change_by_gene <- setDT(aa_change)[,.(aa_change_abbrev=paste(aa_change_abbrev,collapse=" ")),
                         by= c("ANN.FEATUREID","ANN.GENE")]
change_by_gene <- unique(change_by_gene, by="ANN.FEATUREID")
head(change_by_gene)
change_by_gene$ANN.GENE <-as.character(change_by_gene$ANN.GENE)

#Create a new header for the fasta (to use for Mutpred2: id and change separated by space)
mergeID <- change_by_gene %>% 
  unite("MutPred_ID", ANN.GENE, aa_change_abbrev, remove = F, sep = " ")
mergeID

#Change the header of the fasta for mutpred2 specification
mutpred_fasta <- fasta
names(mutpred_fasta) <- mergeID$MutPred_ID
head(mutpred_fasta)
tail(mutpred_fasta)

#save
write.fasta(sequences = as.list(mutpred_fasta), 
            names = names(mutpred_fasta),
            file.out = "Fixed_missense_variants_Astyanax_10hpf_peptides.fasta")

####
#test for all gene if header is correct according to sequence
#prepare the table
#Create a table with the each variant separated by row + corresponding peptide sequence
x <- data.frame(matrix(unlist(mutpred_fasta), nrow=length(mutpred_fasta), byrow=TRUE),stringsAsFactors=FALSE)
x$Mutpred.header <- names(mutpred_fasta)
colnames(x)[1] <- "sequence"
x$gene_name <- gsub("\\ .*", "", x$Mutpred.header)

#Get the peptide sequence per variant + the one letter abbreviation & position of the change
y <- merge(aa_change, x, by.x = "ANN.GENE", by.y = "gene_name", all.x = T)
z <- merge(missense_variants4[,c("ANN.GENE", "ANN.FEATUREID", "aa_change_abbrev", 
                                 "aa_SF_abbrev", "aa_CF_abbrev", "aa_change_position")], 
           y, by = c("ANN.GENE", "ANN.FEATUREID", "aa_change_abbrev"), all.x = T)

#Function to analyse for each variant if the amino acid is found at the correct position
mut_check <- function(aa, pos, seq){
  matching_pos <- gregexpr(aa, seq)[[1]] #get the position of each aa equal to the SF ref
  ifelse(pos %in% matching_pos, #check that the variant is found at the correct position
         Warning_seq <- "no problem",
         Warning_seq <- "WARNING: Incorrect position") #if not, indicate it
  return(Warning_seq)
} 

Warning_seq <- mapply(mut_check, aa = z$aa_SF_abbrev,
                      pos = z$aa_change_position,
                      seq = z$sequence)

pb <- z[which(Warning_seq == "WARNING: Incorrect position"),]
id_error <- unique(pb$ANN.GENE)
#150 variants with wrong position on 75 genes

#Extract the sequence with wrong header  from the whole set
id_error_position <- unlist(sapply(paste0("^",id_error, " "), grep, names(mutpred_fasta)))
id_error_mutpred <- mutpred_fasta[id_error_position]

#remove the X in the sequence (if any)
Xless <- gsub("X", "",id_error_mutpred)
which(gregexpr("X", Xless) >0)
names(Xless) <- names(id_error_mutpred)

#Remake the table with the sequence and one mutation per line (for each gene)
Xless.x <- data.frame(matrix(Xless, nrow=length(Xless), byrow=TRUE),stringsAsFactors=FALSE)
Xless.x$Mutpred.header <- names(id_error_mutpred)
Xless.x$gene_name <- id_error
colnames(Xless.x)[1] <- "sequence"
Xless.y <- merge(aa_change, Xless.x, by.x = "ANN.GENE", by.y = "gene_name", all.y = T) #get also mutation at correct position
Xless.z <- merge(missense_variants4[,c("ANN.GENE", "ANN.FEATUREID", "aa_change_abbrev", 
                                       "aa_SF_abbrev", "aa_CF_abbrev", "aa_change_position")], 
                 Xless.y, by = c("ANN.GENE", "ANN.FEATUREID", "aa_change_abbrev"), all.y = T)

#retest for header-to-sequence position correspondence
Xless.Warning_seq <- mapply(mut_check, aa = Xless.z$aa_SF_abbrev,
                      pos = Xless.z$aa_change_position,
                      seq = Xless.z$sequence)
Xless.pb <- Xless.z[which(Xless.Warning_seq == "WARNING: Incorrect position"),]
unique(Xless.pb$ANN.GENE) #13 genes, 51 mutations with wrong mutation position


####
#correct the sequences detected wrong by Maxime
id_corrected <- which(Xless.x$gene_name %in% unique(Xless.pb$ANN.GENE))
#2 18 19 25 38 47 48 49 57 58 67 68 75
#for the rest, adjust the header manually
names(Xless)[id_corrected] <- c("arhgap19 D276N",
                                "LOC103022313 T92M",
                                "LOC103023024", #remove?
                                "LOC103024947 I19495T G20873S R22323C V22557I A23612V T23622I A23637D V24035M M24146L N25076D P26419S H26438R A26831S K26932N P27457S T27459M I27674V A27687T E27812K R28099C",
                                "LOC103036435 N444S T537S N559I R664Q S716T",
                                "LOC103042045 H621Q",
                                "LOC103042876", #REMOVE,
                                "LOC103045212 A602T D695N L1121Q H1368P I1529L V1596A T1621S I2135V",
                                "LOC111191290 S2549L",
                                "nav2 Q334H A772T S1556A S1820N",
                                "sncaip F612S",                                                  
                                "terf2ip S111A H138N E247V A442T",
                                "znf318 E1516Q V1641A L1676P K1709T")


#final check
x2 <- data.frame(matrix(Xless, nrow=length(Xless), byrow=TRUE),stringsAsFactors=FALSE)
x2$Mutpred.header <- names(Xless)
x2$gene_name <- id_error
x2$Mutation <- gsub("^[A-Za-z0-9]+ ", "", x2$Mutpred.header)
colnames(x2)[1] <- "sequence"

y2 <- x2[-which(x2$gene_name == "LOC103023024" | x2$gene_name == "LOC103042876"),]

z2 <- separate_rows(y2, Mutation)
z2$aa <- gsub("[0-9]+[A-Z]", "", z2$Mutation)
z2$pos <- gsub("[A-Z]", "", z2$Mutation)

#retest for header-to-sequence position correspondence
Xless.Warning_seq2 <- mapply(mut_check, aa = z2$aa,
                            pos = z2$pos,
                            seq = z2$sequence)
Xless.pb2 <- z2[which(Xless.Warning_seq2 == "WARNING: Incorrect position"),]
#no warnings = remove X or changing the header was ok to solve the problem

#remove two genes I could not manually correct
final.corrected <- Xless[-which(names(Xless) == "LOC103023024" | names(Xless) == "LOC103042876")]

#save
write.fasta(sequences = as.list(final.corrected), 
            names = names(final.corrected),
            file.out = "Corrected_variants_Astyanax_10hpf_peptides.fasta")
write.table(z2[,3:6], "Corrected_variants_Astyanax_10hpf_mutation.tab", sep = "\t", col.names = T, row.names = F)

#Correct the initial file with all sequence
all_correct <- x
all_correct[match(z2$gene_name, all_correct$gene_name), ] <- z2[,1:3] #replace wrong header/sequence by the correct one
all_correct <- all_correct[-which(all_correct$gene_name == "LOC103023024" | all_correct$gene_name == "LOC103042876"),] #remove the two genes I could not correct

#save
write.table(all_correct, "Corrected_missense_variants_Astyanax_10hpf_sequence.tab", sep = "\t", col.names = T, row.names = F)


#correct the fasta file & save
write.fasta(sequences = as.list(all_correct$sequence), 
            names = all_correct$Mutpred.header,
            file.out = "Corrected_missense_variants_Astyanax_10hpf_peptides.fasta")

#Then processed as previously done in Policarpo et al., 2021 "Contrasting Gene Decay in Subterranean Vertebrates: Insights from Cavefishes and Fossorial Mammals" DOI: 10.1093/molbev/msaa249

