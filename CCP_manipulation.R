####################################################################################################
##Using the c.DNA position in the Cosmic data to extract the information regarding the alleles
CCP_Std_Mutation<-CCP_Std[which(CCP_Std$Mutation.GRCh37.genome.position!=""),]
CCP_Std_Mutation$Chromosome <- gsub("23", "X",CCP_Std_Mutation$Chromosome)
CCP_Std_Mutation$Chromosome <- gsub("24", "Y",CCP_Std_Mutation$Chromosome)
CCP_Std_Mutation$temp<-paste(CCP_Std_Mutation$Start_Position,CCP_Std_Mutation$End_Position,sep="-")
CCP_Std_Mutation$Mutation.GRCh37.genome.position<-paste(CCP_Std_Mutation$Chromosome,CCP_Std_Mutation$temp,sep=":")
CCP_Std_Mutation$temp<- NULL
row.names(CCP_Std_Mutation)<- NULL


# Divinding the dataframe into the individual mutation segments and doing the proceeding for the further analysis
CCP_Std_Sub<- CCP_Std_Mutation[(grep("(Substitution - coding silent|Substitution - Missense|Substitution - Nonsense)", CCP_Std_Mutation$Mutation.Description)),]
CCP_Std_Sub2<- CCP_Std_Mutation[(grep("(Substitution - coding silent|Substitution - Missense|Substitution - Nonsense)", CCP_Std_Mutation$Mutation.Description, invert = TRUE)),]
rownames(CCP_Std_Sub)<-NULL
rownames(CCP_Std_Sub2)<-NULL

CCP_Std_Sub$Tumor_Seq_Allele <- sapply(strsplit(CCP_Std_Sub$cDNA_Change, ">"), "[", 2)
CCP_Std_Sub$Reference_Allele <- sapply(strsplit(CCP_Std_Sub$cDNA_Change, ">"), "[", 1)
CCP_Std_Sub$Reference_Allele <- sapply(strsplit(CCP_Std_Sub$Reference_Allele, "c."), "[", 2)
CCP_Std_Sub$Reference_Allele <- gsub("[a-z0-9]","", CCP_Std_Sub$Reference_Allele)
CCP_Std_Sub$Reference_Allele <- gsub("[^[:alnum:]]", "",CCP_Std_Sub$Reference_Allele)
CCP_Std_Sub$Reference2Allele <- as.character(with(CCP_Std_Sub, paste(CCP_Std_Sub$Reference_Allele, CCP_Std_Sub$Tumor_Seq_Allele, sep=">")))

Int0gen_Cos_Sub<-CCP_Std_Sub
Int0gen_Cos_Sub$Updated_Accession<- Int0gen_Cos_Sub$Accession.Number
Int0gen_Cos_Sub$Updated_Accession <- gsub("_v68|_v65|_v61|_v70|_v62", "",Int0gen_Cos_Sub$Updated_Accession)
Test_Of_The_result<-Int0gen_Cos_Sub[,c("Chromosome","Start_Position","End_Position","Transcript_Strand","Reference2Allele","ID_sample","ID_tumour","Corrected_Cellline")]
write.table(Test_Of_The_result,"totat.tsv", sep="\t", quote = F, row.names = F)

Consequences<- read.table("consequences.tsv", sep="\t", header=TRUE, stringsAsFactors=F)
CCP_Std_Sub$Reference2Allele <- gsub(">","/",CCP_Std_Sub$Reference2Allele)
Int0gen_Cos_Sub<-CCP_Std_Sub
Int0gen_Cos_Sub$Updated_Accession<- Int0gen_Cos_Sub$Accession.Number
Int0gen_Cos_Sub$Updated_Accession <- gsub("_v68|_v65|_v61|_v70|_v62", "",Int0gen_Cos_Sub$Updated_Accession)
Int0gen_Cos_Sub$Updated_Accession<- sapply(strsplit(Int0gen_Cos_Sub$Updated_Accession, "[.]"), "[", 1)
a<-merge(Int0gen_Cos_Sub,Consequences, by.x= c("Start_Position", "Updated_Accession","Reference2Allele","Transcript_Strand"), 
         by.y= c("START", "TRANSCRIPT_ID","ALLELE","STRAND"), sort= TRUE)

### Repeat this sets of codes twice such as to make update most of NM accessions into ensembl IDs
temp<-unique(as.character(Int0gen_Cos_Sub$Updated_Accession))
temp2<-unique(as.character(a$Updated_Accession))
diff<-setdiff(temp,temp2)
write.table(diff,"Nm2Ent.tsv", sep="\t", quote = F, row.names = F)

lookup2<-read.table("results.txt", sep="\t", header=TRUE, stringsAsFactors=F)
lookup2$X<-NULL
b<-lookup2$Ensembl.Transcript.ID[match(Int0gen_Cos_Sub$Updated_Accession, lookup2$RefSeq.mRNA..e.g..NM_001195597.)]
b[which(is.na(b))] <- Int0gen_Cos_Sub$Updated_Accession[which(is.na(b))]
Int0gen_Cos_Sub$Updated_Accession<-b
a<-merge(Int0gen_Cos_Sub,Consequences, by.x= c("Start_Position", "Updated_Accession","Reference2Allele","Transcript_Strand"), 
         by.y= c("START", "TRANSCRIPT_ID","ALLELE","STRAND"), sort= TRUE)
Annotated_Cosmic_subset1<-a
remove(temp,temp2,diff,b,a,lookup2)

a1<-Annotated_Cosmic_subset1[,c("Data_Set","ID_sample","ID_tumour","CellLines", "Corrected_Cellline","Primary.site","Site.subtype","Primary.histology","Histology.subtype","Tumour.origin","Gene.name","Entrez_Gene_ID","Accession.Number","HGNC.ID","RefSeq_mRNA_Id","RefSeq_prot_ID","SwissProt_acc_Id","SwissProt_entry_Id","Gene.Description","Mutation.ID","Chromosome","Start_Position","End_Position","Transcript_Strand","Mutation.GRCh37.genome.position", "cDNA_Change","Protein_Change","Mutation.Description","Tumor_Seq_Allele","Reference_Allele","Reference2Allele","Updated_Accession")]
require(sqldf)
Unannotated_Cos<- sqldf('SELECT * FROM Int0gen_Cos_Sub EXCEPT SELECT * FROM a1')
Annotated_Cosmic_subset2<-merge(Unannotated_Cos,Consequences, by.x= c("Start_Position","Reference2Allele","Transcript_Strand"), 
                                         by.y= c("START","ALLELE","STRAND"), sort= TRUE)

# The unannotated list of the substituion mutations
temp<-unique(Unannotated_Cos$Start_Position)
temp2<-unique(Annotated_Cosmic_subset2$Start_Position)
diff<-setdiff(temp,temp2)
Unannotated_Substitution<-Int0gen_Cos_Sub[(Int0gen_Cos_Sub$Start_Position %in% diff),]

# The final steps...
Annotated_Cosmic_subset2$Updated_Accession<-NULL
colnames(Annotated_Cosmic_subset2)[34] <- "Unique/Multiple_Updated_TRANSCRIPT_ID"
colnames(Annotated_Cosmic_subset1)[2]<- "Unique/Multiple_Updated_TRANSCRIPT_ID"
remove(a1,diff,temp,temp2)
Annotated_Cosmic_subset2<-Annotated_Cosmic_subset2[,c(1,34,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52)]
colnames(Annotated_Cosmic_subset1) [c(10,11,12,13,14,15,17,18,23,24,27,30)]<-c( "Primary_site","Site_subtype", "Primary_histology","Histology_subtype", "Tumour_origin","Gene_name",
                                                                                "Accession_Number" ,"HGNC_ID","Gene_Description","Mutation_ID","Mutation_GRCh37_genome_position" , "Mutation_Description" )
write.table(Annotated_Cosmic_subset1,"Cosmic_Sub_1.tsv", sep="\t", quote = F, row.names = F)
write.table(Annotated_Cosmic_subset2,"Cosmic_Sub_2.tsv", sep="\t", quote = F, row.names = F)
Annonated_Cosmic_Substitution<-rbind(Annotated_Cosmic_subset1,Annotated_Cosmic_subset2)
write.table(Annonated_Cosmic_Substitution,"Cosmic_Substitution.tsv", sep="\t", quote = F, row.names = F)

remove(Annotated_Cosmic_subset1,Annotated_Cosmic_subset2)

# Inserting the FM-bias values to the original list made of the annotated cosmic list
genes <-read.table("genes.tsv", sep="\t", header=TRUE, stringsAsFactors=F)
Annonated_Cosmic_Substitution$FM_PVALUE<-genes$FM_PVALUE[match(Annonated_Cosmic_Substitution$GENE_ID,genes$GENE_ID)]
Annonated_Cosmic_Substitution$FM_QVALUE<-genes$FM_QVALUE[match(Annonated_Cosmic_Substitution$GENE_ID,genes$GENE_ID)]
Annonated_Cosmic_Substitution$CLUST_ZSCORE<-genes$CLUST_ZSCORE[match(Annonated_Cosmic_Substitution$GENE_ID,genes$GENE_ID)]
Annonated_Cosmic_Substitution$CLUST_PVALUE<-genes$CLUST_PVALUE[match(Annonated_Cosmic_Substitution$GENE_ID,genes$GENE_ID)]
Annonated_Cosmic_Substitution$CLUST_QVALUE<-genes$CLUST_QVALUE[match(Annonated_Cosmic_Substitution$GENE_ID,genes$GENE_ID)]

CCP_Annotated_Substition<-Annonated_Cosmic_Substitution[,c(5,8,9,6,7,10,11,12,13,14,37,18,36,19,17,2,39,38,24,25,1,26,4,3,28,29,35,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57)]
##################################################################################################################
#################################################################################################################
# Dividing the dataframe for the splice site variants and other that are unknown
Not_Substitution_CCP<-CCP_Std_Sub2[(CCP_Std_Sub2$Mutation.Description=="Unknown"),]
rownames(Not_Substitution_CCP)<-NULL
CCP_Std_Splice_site<-Not_Substitution_CCP[(grep("(A>T|A>C|A>G|C>T|C>G|C>A|G>A|G>C|G>T|T>A|T>C|T>G)",
                                                Not_Substitution_CCP$cDNA_Change)),]
rownames(CCP_Std_Splice_site)<-NULL
Unannotated_Unknown<-Not_Substitution_CCP[(grep("(A>T|A>C|A>G|C>T|C>G|C>A|G>A|G>C|G>T|T>A|T>C|T>G)", 
                                                Not_Substitution_CCP$cDNA_Change,invert = TRUE)),]
rownames(Unannotated_Unknown)<-NULL

CCP_Std_Splice_site$Tumor_Seq_Allele <- sapply(strsplit(CCP_Std_Splice_site$cDNA_Change, ">"), "[", 2)
CCP_Std_Splice_site$Reference_Allele <- sapply(strsplit(CCP_Std_Splice_site$cDNA_Change, ">"), "[", 1)
CCP_Std_Splice_site$Reference_Allele  <- gsub("[a-z0-9]","",CCP_Std_Splice_site$Reference_Allele)
CCP_Std_Splice_site$Reference_Allele  <- gsub("[^[:alnum:]]","",CCP_Std_Splice_site$Reference_Allele)
CCP_Std_Splice_site$Reference2Allele <- as.character(with(CCP_Std_Splice_site, paste(CCP_Std_Splice_site$Reference_Allele, CCP_Std_Splice_site$Tumor_Seq_Allele, sep=">")))


Int0gen_Cos_Sub<-CCP_Std_Splice_site
Test_Of_The_result<-Int0gen_Cos_Sub[,c("Chromosome","Start_Position","End_Position","Transcript_Strand","Reference2Allele","ID_sample","ID_tumour","Corrected_Cellline")]
write.table(Test_Of_The_result,"splicesite and others unknown.tsv", sep="\t", quote = F, row.names = F)

Consequences<- read.table("consequences.tsv", sep="\t", header=TRUE, stringsAsFactors=F)
Int0gen_Cos_Sub$Reference2Allele <- gsub(">","/",Int0gen_Cos_Sub$Reference2Allele)
Int0gen_Cos_Sub$Updated_Accession<- Int0gen_Cos_Sub$Accession.Number
Int0gen_Cos_Sub$Updated_Accession <- gsub("_v68|_v65|_v61|_v70|_v62", "",Int0gen_Cos_Sub$Updated_Accession)
Int0gen_Cos_Sub$Updated_Accession<- sapply(strsplit(Int0gen_Cos_Sub$Updated_Accession, "[.]"), "[", 1)
a<-merge(Int0gen_Cos_Sub,Consequences, by.x= c("Start_Position", "Updated_Accession","Reference2Allele","Transcript_Strand"), 
         by.y= c("START", "TRANSCRIPT_ID","ALLELE","STRAND"), sort= TRUE)

### Repeat this sets of codes twice such as to make update most of NM accessions into ensembl IDs
temp1<-unique(Int0gen_Cos_Sub$Updated_Accession)
temp2<-unique(a$Updated_Accession)
diff<-setdiff(temp1,temp2)
write.table(diff, "Nm2Ent_2.tsv", sep="\t", quote = F, row.names = F)
lookup2<-read.table("results_2.txt", sep="\t", header=TRUE, stringsAsFactors=F)

b<-lookup2$Ensembl.Transcript.ID[match(Int0gen_Cos_Sub$Updated_Accession, lookup2$RefSeq.mRNA..e.g..NM_001195597.)]
b[which(is.na(b))] <- Int0gen_Cos_Sub$Updated_Accession[which(is.na(b))]
Int0gen_Cos_Sub$Updated_Accession<-b
a<-merge(Int0gen_Cos_Sub,Consequences, by.x= c("Start_Position", "Updated_Accession","Reference2Allele","Transcript_Strand"), 
         by.y= c("START", "TRANSCRIPT_ID","ALLELE","STRAND"), sort= TRUE)
###################################################################################################################
Annotated_Cosmic_Splice_SUbset1<-a
a1<-Annotated_Cosmic_Splice_SUbset1[,c("Data_Set","ID_sample","ID_tumour","CellLines", "Corrected_Cellline","Primary.site","Site.subtype","Primary.histology","Histology.subtype","Tumour.origin","Gene.name","Entrez_Gene_ID","Accession.Number","HGNC.ID","RefSeq_mRNA_Id","RefSeq_prot_ID","SwissProt_acc_Id","SwissProt_entry_Id","Gene.Description","Mutation.ID","Chromosome","Start_Position","End_Position","Transcript_Strand","Mutation.GRCh37.genome.position", "cDNA_Change","Protein_Change","Mutation.Description","Tumor_Seq_Allele","Reference_Allele","Reference2Allele","Updated_Accession")]
require(sqldf)
Unannotated_cos<- sqldf('SELECT * FROM Int0gen_Cos_Sub EXCEPT SELECT * FROM a1')
Annotated_Cosmic_Splice_subset2<-merge(Unannotated_cos,Consequences, by.x= c("Start_Position","Reference2Allele","Transcript_Strand"), 
                                   by.y= c("START","ALLELE","STRAND"), sort= TRUE)


# The unannotated list of the unknown-splice site mutations
temp<-unique(Unannotated_cos$Start_Position)
temp2<-unique(Annotated_Cosmic_Splice_subset2$Start_Position)
diff<-setdiff(temp,temp2)
Unannotated_Unknown_splicesite<-Int0gen_Cos_Sub[(Int0gen_Cos_Sub$Start_Position %in% diff),]

# The final steps...

Annotated_Cosmic_Splice_subset2$Updated_Accession<-NULL
colnames(Annotated_Cosmic_Splice_subset2)[34] <- "Unique/Multiple_Updated_TRANSCRIPT_ID"
colnames(Annotated_Cosmic_Splice_SUbset1)[2]<- "Unique/Multiple_Updated_TRANSCRIPT_ID"
remove(a1,diff,temp,temp2)
Annotated_Cosmic_Splice_subset2<-Annotated_Cosmic_Splice_subset2[,c(1,34,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52)]
colnames(Annotated_Cosmic_Splice_SUbset1) [c(10,11,12,13,14,15,17,18,23,24,27,30)]<-c( "Primary_site","Site_subtype", "Primary_histology","Histology_subtype", "Tumour_origin","Gene_name",
                                                                                       "Accession_Number" ,"HGNC_ID","Gene_Description","Mutation_ID","Mutation_GRCh37_genome_position" , "Mutation_Description" )
write.table(Annotated_Cosmic_Splice_SUbset1,"Cosmic_Splice_site_1.tsv", sep="\t", quote = F, row.names = F)
write.table(Annotated_Cosmic_Splice_subset2,"Cosmic_Splice_site_2.tsv", sep="\t", quote = F, row.names = F)
Annonated_Cosmic_Splice_Site<-rbind(Annotated_Cosmic_Splice_SUbset1,Annotated_Cosmic_Splice_subset2)
write.table(Annonated_Cosmic_Splice_Site,"Cosmic_Splice_Site.tsv", sep="\t", quote = F, row.names = F)

remove(Annotated_Cosmic_Splice_SUbset1,Annotated_Cosmic_Splice_subset2)

# Inserting the FM-bias values to the original list made of the annotated cosmic list
genes <-read.table("genes.tsv", sep="\t", header=TRUE, stringsAsFactors=F)
Annonated_Cosmic_Splice_Site$FM_PVALUE<-genes$FM_PVALUE[match(Annonated_Cosmic_Splice_Site$GENE_ID,genes$GENE_ID)]
Annonated_Cosmic_Splice_Site$FM_QVALUE<-genes$FM_QVALUE[match(Annonated_Cosmic_Splice_Site$GENE_ID,genes$GENE_ID)]
Annonated_Cosmic_Splice_Site$CLUST_ZSCORE<-genes$CLUST_ZSCORE[match(Annonated_Cosmic_Splice_Site$GENE_ID,genes$GENE_ID)]
Annonated_Cosmic_Splice_Site$CLUST_PVALUE<-genes$CLUST_PVALUE[match(Annonated_Cosmic_Splice_Site$GENE_ID,genes$GENE_ID)]
Annonated_Cosmic_Splice_Site$CLUST_QVALUE<-genes$CLUST_QVALUE[match(Annonated_Cosmic_Splice_Site$GENE_ID,genes$GENE_ID)]


CCP_Annotated_Splice_Site<-Annonated_Cosmic_Splice_Site[,c(5,8,9,6,7,10,11,12,13,14,37,18,36,19,17,2,39,38,24,25,1,26,4,3,28,29,35,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57)]
#####################################################################################################
#####################################################################################################
# Subsetting the orginal mutation datasets
# Insertion variants from the original COSMIC_STD_MUTATION dataframe
temp<-CCP_Std_Sub2[(CCP_Std_Sub2$Mutation.Description=="Insertion - Frameshift"),]
temp2<-CCP_Std_Sub2[(CCP_Std_Sub2$Mutation.Description=="Insertion - In frame"),]
temp3<-CCP_Std_Sub2[(CCP_Std_Sub2$Mutation.Description=="Complex - insertion inframe"),]

CCP_std_INS<-rbind(temp,temp2,temp3)
row.names(CCP_std_INS)<-NULL
str(CCP_std_INS)
remove(temp,temp2,temp3)

CCP_std_INS$Tumor_Seq_Allele<- sapply(strsplit(CCP_std_INS$cDNA_Change, "c."), "[", 2)
CCP_std_INS$Tumor_Seq_Allele<-gsub("[a-z0-9]","", CCP_std_INS$Tumor_Seq_Allele)
CCP_std_INS$Tumor_Seq_Allele<-gsub("[^[:alnum:]]","",CCP_std_INS$Tumor_Seq_Allele)
CCP_std_INS$Reference_Allele <- "-"
CCP_std_INS$Reference2Allele <- as.character(with(CCP_std_INS, paste(CCP_std_INS$Reference_Allele, CCP_std_INS$Tumor_Seq_Allele, sep=">")))

Int0gen_Cos_Sub<-CCP_std_INS
Test_Of_The_result<-Int0gen_Cos_Sub[,c("Chromosome","Start_Position","End_Position","Transcript_Strand","Reference2Allele","ID_sample","ID_tumour","Corrected_Cellline")]
write.table(Test_Of_The_result,"insertion.tsv", sep="\t", quote = F, row.names = F)


Consequences<- read.table("consequences.tsv", sep="\t", header=TRUE, stringsAsFactors=F)
CCP_std_INS$Reference2Allele <- gsub(">","/",CCP_std_INS$Reference2Allele)
Int0gen_Cos_Sub<-CCP_std_INS
Int0gen_Cos_Sub$Updated_Accession<- Int0gen_Cos_Sub$Accession.Number
Int0gen_Cos_Sub$Updated_Accession <- gsub("_v68|_v65|_v61|_v70|_v62", "",Int0gen_Cos_Sub$Updated_Accession)
Int0gen_Cos_Sub$Updated_Accession<- sapply(strsplit(Int0gen_Cos_Sub$Updated_Accession, "[.]"), "[", 1)
a<-merge(Int0gen_Cos_Sub,Consequences, by.x= c("End_Position", "Updated_Accession","Reference2Allele","Transcript_Strand"), 
         by.y= c("START", "TRANSCRIPT_ID","ALLELE","STRAND"), sort= TRUE)

## Repeat these steps twice such that almost NM acessions are converted to ENSEMBL IDs
temp1<-unique(as.character(Int0gen_Cos_Sub$Updated_Accession))
temp2<-unique(as.character(a$Updated_Accession))
diff<-setdiff(temp1,temp2)
write.table(diff,"Nm2Ent_3.tsv", sep="\t", quote = F, row.names = F)
lookup2<-read.table("results_3.txt", sep="\t", header=TRUE, stringsAsFactors=F)
b<-lookup2$Ensembl.Transcript.ID[match(Int0gen_Cos_Sub$Updated_Accession, lookup2$RefSeq.mRNA..e.g..NM_001195597.)]
b[which(is.na(b))] <- Int0gen_Cos_Sub$Updated_Accession[which(is.na(b))]
Int0gen_Cos_Sub$Updated_Accession<-b
a<-merge(Int0gen_Cos_Sub,Consequences, by.x= c("End_Position", "Updated_Accession","Reference2Allele","Transcript_Strand"), 
         by.y= c("START", "TRANSCRIPT_ID","ALLELE","STRAND"), sort= TRUE)

Annotated_Cosmic_Insertion_subset1<-a
remove(temp1,temp2,diff,b,a,lookup2)

a1<-Annotated_Cosmic_Insertion_subset1[,c("Data_Set","ID_sample","ID_tumour","CellLines", "Corrected_Cellline","Primary.site","Site.subtype","Primary.histology","Histology.subtype","Tumour.origin","Gene.name","Entrez_Gene_ID","Accession.Number","HGNC.ID","RefSeq_mRNA_Id","RefSeq_prot_ID","SwissProt_acc_Id","SwissProt_entry_Id","Gene.Description","Mutation.ID","Chromosome","Start_Position","End_Position","Transcript_Strand","Mutation.GRCh37.genome.position", "cDNA_Change","Protein_Change","Mutation.Description","Tumor_Seq_Allele","Reference_Allele","Reference2Allele","Updated_Accession")]
require(sqldf)
Unannotated_cos<- sqldf('SELECT * FROM Int0gen_Cos_Sub EXCEPT SELECT * FROM a1')
Annotated_Cosmic_Insertion_subset2<-merge(Unannotated_cos,Consequences, by.x= c("End_Position","Reference2Allele","Transcript_Strand"), 
                                          by.y= c("START","ALLELE","STRAND"), sort= TRUE)

# The unannotated list of the insertion genomic mutations
temp<-unique(Unannotated_cos$End_Position)
temp2<-unique(Annotated_Cosmic_Insertion_subset2$End_Position)
diff<-setdiff(temp,temp2)
Unannotated_Insertion<-Unannotated_cos[(Unannotated_cos$End_Position %in% diff),]

# The final steps...
Annotated_Cosmic_Insertion_subset2$Updated_Accession<-NULL
colnames(Annotated_Cosmic_Insertion_subset2)[34] <- "Unique/Multiple_Updated_TRANSCRIPT_ID"
colnames(Annotated_Cosmic_Insertion_subset1)[2]<- "Unique/Multiple_Updated_TRANSCRIPT_ID"
Annotated_Cosmic_Insertion_subset2<-Annotated_Cosmic_Insertion_subset2[,c(1,34,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52)]
colnames(Annotated_Cosmic_Insertion_subset1) [c(10,11,12,13,14,15,17,18,23,24,27,30)]<-c( "Primary_site","Site_subtype", "Primary_histology","Histology_subtype", "Tumour_origin","Gene_name",
                                                                                          "Accession_Number" ,"HGNC_ID","Gene_Description","Mutation_ID","Mutation_GRCh37_genome_position" , "Mutation_Description" )
write.table(Annotated_Cosmic_Insertion_subset1,"Cosmic_Insertion_1.tsv", sep="\t", quote = F, row.names = F)
write.table(Annotated_Cosmic_Insertion_subset2,"Cosmic_Insertion_2.tsv", sep="\t", quote = F, row.names = F)
Annonated_Cosmic_Insertion<-rbind(Annotated_Cosmic_Insertion_subset1,Annotated_Cosmic_Insertion_subset2)
write.table(Annonated_Cosmic_Insertion,"Cosmic_Insertion.tsv", sep="\t", quote = F, row.names = F)

remove(Annotated_Cosmic_Insertion_subset1,Annotated_Cosmic_Insertion_subset2)

genes <-read.table("genes.tsv", sep="\t", header=TRUE, stringsAsFactors=F)
Annonated_Cosmic_Insertion$FM_PVALUE<-genes$FM_PVALUE[match(Annonated_Cosmic_Insertion$GENE_ID,genes$GENE_ID)]
Annonated_Cosmic_Insertion$FM_QVALUE<-genes$FM_QVALUE[match(Annonated_Cosmic_Insertion$GENE_ID,genes$GENE_ID)]
Annonated_Cosmic_Insertion$CLUST_ZSCORE<-genes$CLUST_ZSCORE[match(Annonated_Cosmic_Insertion$GENE_ID,genes$GENE_ID)]
Annonated_Cosmic_Insertion$CLUST_PVALUE<-genes$CLUST_PVALUE[match(Annonated_Cosmic_Insertion$GENE_ID,genes$GENE_ID)]
Annonated_Cosmic_Insertion$CLUST_QVALUE<-genes$CLUST_QVALUE[match(Annonated_Cosmic_Insertion$GENE_ID,genes$GENE_ID)]


CCP_Annotated_Insertion<-Annonated_Cosmic_Insertion[,c(5,8,9,6,7,10,11,12,13,14,37,18,36,19,17,2,39,38,24,25,26,1,4,3,28,29,35,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57)]
###########################################################################################################
############################################################################################################

# Deletion dataframes from the original CCP_Std_Sub2 dataset
temp<-CCP_Std_Sub2[(CCP_Std_Sub2$Mutation.Description=="Deletion - Frameshift"),]
temp2<-CCP_Std_Sub2[(CCP_Std_Sub2$Mutation.Description=="Deletion - In frame"),]
temp3<-CCP_Std_Sub2[(CCP_Std_Sub2$Mutation.Description=="Complex - deletion inframe"),]
CCP_Std_Del<-rbind(temp,temp2,temp3)
row.names(CCP_Std_Del)<-NULL

# Dividing the columns of the deletion into 2 segments of the one with the alphabetic character
# after del and other with the numeric character after it
CCP_Std_Del$Reference_Allele<- sapply(strsplit(CCP_Std_Del$cDNA_Change, "c."), "[", 2)
CCP_Std_Del$Reference_Allele <- gsub("del", "|", CCP_Std_Del$Reference_Allele)
CCP_Std_Del$Reference_Allele<- sapply(strsplit(CCP_Std_Del$Reference_Allele,"\\|"), "[", 2)
CCP_Std_Del$Reference_Allele<- toupper(CCP_Std_Del$Reference_Allele)
CCP_Std_Del$Tumor_Seq_Allele<-"-"


CCP_Std_Del$Reference2Allele <- as.character(with(CCP_Std_Del, paste(CCP_Std_Del$Reference_Allele, CCP_Std_Del$Tumor_Seq_Allele, sep=">")))
Int0gen_Cos_Sub<-CCP_Std_Del
Test_Of_The_result<-Int0gen_Cos_Sub[,c("Chromosome","Start_Position","End_Position","Transcript_Strand","Reference2Allele","ID_sample","ID_tumour","Corrected_Cellline")]
write.table(Test_Of_The_result,"deletion.tsv", sep="\t", quote = F, row.names = F)

Consequences<- read.table("consequences.tsv", sep="\t", header=TRUE, stringsAsFactors=F)
CCP_Std_Del$Reference2Allele <- gsub(">","/",CCP_Std_Del$Reference2Allele)
Int0gen_Cos_Sub<-CCP_Std_Del
Int0gen_Cos_Sub$Updated_Accession<- Int0gen_Cos_Sub$Accession.Number
Int0gen_Cos_Sub$Updated_Accession <- gsub("_v68|_v65|_v61|_v70|_v62", "",Int0gen_Cos_Sub$Updated_Accession)
Int0gen_Cos_Sub$Updated_Accession<- sapply(strsplit(Int0gen_Cos_Sub$Updated_Accession, "[.]"), "[", 1)
a<-merge(Int0gen_Cos_Sub,Consequences, by.x= c("Start_Position", "Updated_Accession","Reference2Allele","Transcript_Strand"), 
         by.y= c("START", "TRANSCRIPT_ID","ALLELE","STRAND"), sort= TRUE)

## Repeat the steps twice to get more refined merge of the two datasets
temp1<-unique(as.character(Int0gen_Cos_Sub$Updated_Accession))
temp2<-unique(as.character(a$Updated_Accession))
diff<-setdiff(temp1,temp2)
write.table(diff,"Nm2Ent_4.tsv", sep="\t", quote = F, row.names = F)
lookup2<-read.table("results_4.txt", sep="\t", header=TRUE, stringsAsFactors=F)
b<-lookup2$Ensembl.Transcript.ID[match(Int0gen_Cos_Sub$Updated_Accession, lookup2$RefSeq.mRNA..e.g..NM_001195597.)]
b[which(is.na(b))] <- Int0gen_Cos_Sub$Updated_Accession[which(is.na(b))]
Int0gen_Cos_Sub$Updated_Accession<-b
a<-merge(Int0gen_Cos_Sub,Consequences, by.x= c("Start_Position", "Updated_Accession","Reference2Allele","Transcript_Strand"), 
         by.y= c("START", "TRANSCRIPT_ID","ALLELE","STRAND"), sort= TRUE)

Annotated_Cosmic_Deletion_Subset1<-a
remove(temp1,temp2,diff,b,a,lookup2)

a1<-Annotated_Cosmic_Deletion_Subset1[,c("Data_Set","ID_sample","ID_tumour","CellLines", "Corrected_Cellline","Primary.site","Site.subtype","Primary.histology","Histology.subtype","Tumour.origin","Gene.name","Entrez_Gene_ID","Accession.Number","HGNC.ID","RefSeq_mRNA_Id","RefSeq_prot_ID","SwissProt_acc_Id","SwissProt_entry_Id","Gene.Description","Mutation.ID","Chromosome","Start_Position","End_Position","Transcript_Strand","Mutation.GRCh37.genome.position", "cDNA_Change","Protein_Change","Mutation.Description","Reference_Allele","Tumor_Seq_Allele","Reference2Allele","Updated_Accession")]
require(sqldf)
Unannotated_cos<- sqldf('SELECT * FROM Int0gen_Cos_Sub EXCEPT SELECT * FROM a1')
Annotated_Cosmic_Deletion_Subset2<-merge(Unannotated_cos,Consequences, by.x= c("Start_Position","Reference2Allele","Transcript_Strand"), 
                                         by.y= c("START","ALLELE","STRAND"), sort= TRUE)


# The unannotated list of the deletion genomic mutations
temp<-unique(Unannotated_cos$Start_Position)
temp2<-unique(Annotated_Cosmic_Deletion_Subset2$Start_Position)
diff<-setdiff(temp,temp2)
Unannotated_Deletion<-Unannotated_cos[(Unannotated_cos$Start_Position %in% diff),]
remove(temp,temp2,diff,Unannotated_cos, a1)

# The final steps...
Annotated_Cosmic_Deletion_Subset2$Updated_Accession<-NULL
colnames(Annotated_Cosmic_Deletion_Subset2)[34] <- "Unique/Multiple_Updated_TRANSCRIPT_ID"
colnames(Annotated_Cosmic_Deletion_Subset1)[2]<- "Unique/Multiple_Updated_TRANSCRIPT_ID"
Annotated_Cosmic_Deletion_Subset2<-Annotated_Cosmic_Deletion_Subset2[,c(1,34,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52)]
colnames(Annotated_Cosmic_Deletion_Subset1) [c(10,11,12,13,14,15,17,18,23,24,27,30)]<-c( "Primary_site","Site_subtype", "Primary_histology","Histology_subtype", "Tumour_origin","Gene_name",
                                                                                          "Accession_Number" ,"HGNC_ID","Gene_Description","Mutation_ID","Mutation_GRCh37_genome_position" , "Mutation_Description" )
write.table(Annotated_Cosmic_Deletion_Subset1,"Cosmic_Deletion_1.tsv", sep="\t", quote = F, row.names = F)
write.table(Annotated_Cosmic_Deletion_Subset2,"Cosmic_Deletion_2.tsv", sep="\t", quote = F, row.names = F)
Annonated_Cosmic_Deletion<-rbind(Annotated_Cosmic_Deletion_Subset1,Annotated_Cosmic_Deletion_Subset2)
write.table(Annonated_Cosmic_Deletion,"Cosmic_Deletion.tsv", sep="\t", quote = F, row.names = F)

remove(Annotated_Cosmic_Deletion_Subset1,Annotated_Cosmic_Deletion_Subset2)

genes <-read.table("genes.tsv", sep="\t", header=TRUE, stringsAsFactors=F)
Annonated_Cosmic_Deletion$FM_PVALUE<-genes$FM_PVALUE[match(Annonated_Cosmic_Deletion$GENE_ID,genes$GENE_ID)]
Annonated_Cosmic_Deletion$FM_QVALUE<-genes$FM_QVALUE[match(Annonated_Cosmic_Deletion$GENE_ID,genes$GENE_ID)]
Annonated_Cosmic_Deletion$CLUST_ZSCORE<-genes$CLUST_ZSCORE[match(Annonated_Cosmic_Deletion$GENE_ID,genes$GENE_ID)]
Annonated_Cosmic_Deletion$CLUST_PVALUE<-genes$CLUST_PVALUE[match(Annonated_Cosmic_Deletion$GENE_ID,genes$GENE_ID)]
Annonated_Cosmic_Deletion$CLUST_QVALUE<-genes$CLUST_QVALUE[match(Annonated_Cosmic_Deletion$GENE_ID,genes$GENE_ID)]


CCP_Annotated_Deletion<-Annonated_Cosmic_Deletion[,c(5,8,9,6,7,10,11,12,13,14,37,18,36,19,17,2,39,38,24,25,1,26,4,3,28,29,35,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57)]
########################################################################################################
#######################################################################################################
# Other mutation.descriptions from the CCP_Std_Mutation dataset
temp2<-CCP_Std_Sub2[(CCP_Std_Sub2$Mutation.Description=="Nonstop extension"),]
temp3<-CCP_Std_Sub2[(CCP_Std_Sub2$Mutation.Description=="No detectable mRNA/protein"),]
CCP_Std_other_variants_1<-rbind(temp2,temp3)
row.names(CCP_Std_other_variants_1)<-NULL


CCP_Std_other_variants_1$Tumor_Seq_Allele <- sapply(strsplit(CCP_Std_other_variants_1$cDNA_Change, ">"), "[", 2)
CCP_Std_other_variants_1$Reference_Allele <- sapply(strsplit(CCP_Std_other_variants_1$cDNA_Change, ">"), "[", 1)
CCP_Std_other_variants_1$Reference_Allele <- sapply(strsplit(CCP_Std_other_variants_1$Reference_Allele, "c."), "[", 2)
CCP_Std_other_variants_1$Reference_Allele <- gsub("[a-z0-9]","", CCP_Std_other_variants_1$Reference_Allele)
CCP_Std_other_variants_1$Reference_Allele <- gsub("[^[:alnum:]]", "",CCP_Std_other_variants_1$Reference_Allele)
CCP_Std_other_variants_1$Reference2Allele <- as.character(with(CCP_Std_other_variants_1,
                                                        paste(CCP_Std_other_variants_1$Reference_Allele,
                                                        CCP_Std_other_variants_1$Tumor_Seq_Allele, sep=">")))


Int0gen_Cos_Sub<-CCP_Std_other_variants_1
Test_Of_The_result<-Int0gen_Cos_Sub[,c("Chromosome","Start_Position","End_Position","Transcript_Strand","Reference2Allele","ID_sample","ID_tumour","Corrected_Cellline")]
write.table(Test_Of_The_result,"Other_variants_1.tsv", sep="\t", quote = F, row.names = F)


Consequences<- read.table("consequences.tsv", sep="\t", header=TRUE, stringsAsFactors=F)
CCP_Std_other_variants_1$Reference2Allele <- gsub(">","/",CCP_Std_other_variants_1$Reference2Allele)
Int0gen_Cos_Sub<-CCP_Std_other_variants_1
Int0gen_Cos_Sub$Updated_Accession<- Int0gen_Cos_Sub$Accession.Number
Int0gen_Cos_Sub$Updated_Accession <- gsub("_v68|_v65|_v61|_v70|_v62", "",Int0gen_Cos_Sub$Updated_Accession)
Int0gen_Cos_Sub$Updated_Accession<- sapply(strsplit(Int0gen_Cos_Sub$Updated_Accession, "[.]"), "[", 1)
a<-merge(Int0gen_Cos_Sub,Consequences, by.x= c("Start_Position", "Updated_Accession","Reference2Allele","Transcript_Strand"), 
         by.y= c("START", "TRANSCRIPT_ID","ALLELE","STRAND"), sort= TRUE)

### Repeat this sets of codes twice such as to make update most of NM accessions into ensembl IDs
temp<-unique(as.character(Int0gen_Cos_Sub$Updated_Accession))
temp2<-unique(as.character(a$Updated_Accession))
diff<-setdiff(temp,temp2)
write.table(diff,"Nm2Ent_5.tsv", sep="\t", quote = F, row.names = F)
lookup2<-read.table("results_5.txt", sep="\t", header=TRUE, stringsAsFactors=F)
b<-lookup2$Ensembl.Transcript.ID[match(Int0gen_Cos_Sub$Updated_Accession, lookup2$RefSeq.mRNA..e.g..NM_001195597.)]
b[which(is.na(b))] <- Int0gen_Cos_Sub$Updated_Accession[which(is.na(b))]
Int0gen_Cos_Sub$Updated_Accession<-b
a<-merge(Int0gen_Cos_Sub,Consequences, by.x= c("Start_Position", "Updated_Accession","Reference2Allele","Transcript_Strand"), 
         by.y= c("START", "TRANSCRIPT_ID","ALLELE","STRAND"), sort= TRUE)

Annotated_Cosmic_other_variant1<-a
remove(temp,temp2,diff,b,a,lookup2)

a1<-Annotated_Cosmic_other_variant1[,c("Data_Set","ID_sample","ID_tumour","CellLines", "Corrected_Cellline","Primary.site","Site.subtype","Primary.histology","Histology.subtype","Tumour.origin","Gene.name","Entrez_Gene_ID","Accession.Number","HGNC.ID","RefSeq_mRNA_Id","RefSeq_prot_ID","SwissProt_acc_Id","SwissProt_entry_Id","Gene.Description","Mutation.ID","Chromosome","Start_Position","End_Position","Transcript_Strand","Mutation.GRCh37.genome.position", "cDNA_Change","Protein_Change","Mutation.Description","Tumor_Seq_Allele","Reference_Allele","Reference2Allele","Updated_Accession")]
require(sqldf)
Unannotated_Cos<- sqldf('SELECT * FROM Int0gen_Cos_Sub EXCEPT SELECT * FROM a1')
Annotated_Cosmic_other_variant2<-merge(Unannotated_Cos,Consequences, by.x= c("Start_Position","Reference2Allele","Transcript_Strand"), 
                                by.y= c("START","ALLELE","STRAND"), sort= TRUE)

# The final steps...
Annotated_Cosmic_other_variant2$Updated_Accession<-NULL
colnames(Annotated_Cosmic_other_variant2)[34] <- "Unique/Multiple_Updated_TRANSCRIPT_ID"
colnames(Annotated_Cosmic_other_variant1)[2]<- "Unique/Multiple_Updated_TRANSCRIPT_ID"
remove(a1,diff,temp,temp2)
Annotated_Cosmic_other_variant2<-Annotated_Cosmic_other_variant2[,c(1,34,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52)]
colnames(Annotated_Cosmic_other_variant1) [c(10,11,12,13,14,15,17,18,23,24,27,30)]<-c( "Primary_site","Site_subtype", "Primary_histology","Histology_subtype", "Tumour_origin","Gene_name",
                                                                                "Accession_Number" ,"HGNC_ID","Gene_Description","Mutation_ID","Mutation_GRCh37_genome_position" , "Mutation_Description" )
write.table(Annotated_Cosmic_other_variant1,"Cosmic_other_variant_1.tsv", sep="\t", quote = F, row.names = F)
write.table(Annotated_Cosmic_other_variant2,"Cosmic_other_variant_2.tsv", sep="\t", quote = F, row.names = F)
Annonated_Cosmic_other_variant<-rbind(Annotated_Cosmic_other_variant1,Annotated_Cosmic_other_variant2)
write.table(Annonated_Cosmic_other_variant,"Cosmic_Other_variant.tsv", sep="\t", quote = F, row.names = F)

remove(Annotated_Cosmic_other_variant1,Annotated_Cosmic_other_variant2)

genes <-read.table("genes.tsv", sep="\t", header=TRUE, stringsAsFactors=F)
Annonated_Cosmic_other_variant$FM_PVALUE<-genes$FM_PVALUE[match(Annonated_Cosmic_other_variant$GENE_ID,genes$GENE_ID)]
Annonated_Cosmic_other_variant$FM_QVALUE<-genes$FM_QVALUE[match(Annonated_Cosmic_other_variant$GENE_ID,genes$GENE_ID)]
Annonated_Cosmic_other_variant$CLUST_ZSCORE<-genes$CLUST_ZSCORE[match(Annonated_Cosmic_other_variant$GENE_ID,genes$GENE_ID)]
Annonated_Cosmic_other_variant$CLUST_PVALUE<-genes$CLUST_PVALUE[match(Annonated_Cosmic_other_variant$GENE_ID,genes$GENE_ID)]
Annonated_Cosmic_other_variant$CLUST_QVALUE<-genes$CLUST_QVALUE[match(Annonated_Cosmic_other_variant$GENE_ID,genes$GENE_ID)]



CCP_Annotated_other_variant<-Annonated_Cosmic_other_variant[,c(5,8,9,6,7,10,11,12,13,14,37,18,36,19,17,2,39,38,24,25,1,26,4,3,28,29,35,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57)]
write.table(Annotated_Cosmic_other_variant2,"Cosmic_other_variant_2.tsv", sep="\t", quote = F, row.names = F)

#####################################################################################################################
######################################################################################################################
temp1<-CCP_Std_Sub2[(CCP_Std_Sub2$Mutation.Description=="Complex - frameshift"),]
temp2<-CCP_Std_Sub2[(CCP_Std_Sub2$Mutation.Description=="Complex - compound substitution"),]
Unannotated_CCP_Std_complex<-rbind(temp1,temp2)
row.names(Unannotated_CCP_Std_complex)<-NULL
remove(temp1,temp2)
########################################################################################################################
#######################################################################################################################
CCP_Annotated<-rbind(CCP_Annotated_Substition,CCP_Annotated_Splice_Site,CCP_Annotated_other_variant,CCP_Annotated_Insertion,CCP_Annotated_Deletion)
write.table(CCP_Annotated,"CCP_Total.tsv", sep="\t", quote = F, row.names = F)
