write.table(NM2ENST,"22.tsv", sep="\t", quote = F, row.names = F)
lookup<- read.table("results.txt", sep="\t", header=TRUE, stringsAsFactors=F)

Consequences<- read.table("consequences.tsv", sep="\t", header=TRUE, stringsAsFactors=F)
CCLE_Std$Updated.ENSEMBL.ID<-as.character(CCLE_Std$RefSeq_mRNA_Id)
ex<-lookup$ENSEMBL_acc[match(CCLE_Std$Updated.ENSEMBL.ID, lookup$Ref_mRNA)]
CCLE_Std$Updated.ENSEMBL.ID<-ex
a <-CCLE_Std$Updated.ENSEMBL.ID[is.na(ex)]
test<-CCLE_Std[which(CCLE_Std$Updated.ENSEMBL.ID=="ENST00000377038"),]
test <-CCLE_Std$Updated.ENSEMBL.ID[is.na(CCLE_Std$Updated.ENSEMBL.ID)]
is.na(ex)

CCLE_Std$Reference_Allele<-CCLE$Reference_Allele
identical(CCLE$Tumor_Seq_Allele1,CCLE$Tumor_Seq_Allele2)
CCLE_Std$Tumor_Seq_Allele<-CCLE$Tumor_Seq_Allele1
CCLE_Std$Reference2Allele <- as.character(with(CCLE_Std, paste(CCLE_Std$Reference_Allele, CCLE_Std$Tumor_Seq_Allele, sep="/")))
a<-merge(CCLE_Std,Consequences, by.x= c("Start_Position", "Updated.ENSEMBL.ID","Reference2Allele"), 
         by.y= c("START", "TRANSCRIPT_ID","ALLELE"), sort= TRUE)

a1<-a[c("Data_Set","ID_sample","ID_tumour","CellLines", "Corrected_Cellline","Primary.site","Site.subtype","Primary.histology","Histology.Subtype","Tumour.origin","Gene.name","Entrez_Gene_ID","Accession.Number","HGNC.ID","RefSeq_mRNA_Id","RefSeq_prot_ID","SwissProt_acc_Id","SwissProt_entry_Id","Gene.Description","Mutation.ID","Chromosome","Start_Position","End_Position","Transcript_Strand","Mutation.GRCh37.genome.position", "cDNA_Change","Protein_Change","Mutation.Description","Updated.ENSEMBL.ID","Reference_Allele","Tumor_Seq_Allele","Reference2Allele")]
require(sqldf)
Unannotated <- sqldf('SELECT * FROM CCLE_Std EXCEPT SELECT * FROM a1')

b1<-merge(Unannotated,Consequences, by.x= c("End_Position", "Updated_ENSEMBL_ID","Reference2Allele"), 
         by.y= c("START", "TRANSCRIPT_ID","ALLELE"), sort= TRUE)
a<-a[,c(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,1,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53)]
b1<-b1[,c(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,1,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53)]

colnames(a)[c(1,8,9,10,11,12,13,15,16,21,22,27,30)]<-c("Updated_ENSEMBL_ID","Primary_site","Site_subtype","Primary_histology","Histology_Subtype","Tumour_origin","Gene_name", "Accession_Number","HGNC_ID","Gene_Description","Mutation_ID","Mutation_GRCh37_genome_position","Mutation_Description")
colnames(b1)[c(1,8,9,10,11,12,13,15,16,21,27,30)]<-c("Updated_ENSEMBL_ID","Primary_site","Site_subtype","Primary_histology","Histology_Subtype","Tumour_origin","Gene_name", "Accession_Number","HGNC_ID","Gene_Description","Mutation_GRCh37_genome_position","Mutation_Description")
Annotated<-rbind(a,b1)
colnames(Annotated) [c(1,8,9,10,11,12,13,15,16,21,22,27,30)]<-c("Updated.ENSEMBL.ID","Primary.site","Site.subtype","Primary.histology","Histology.Subtype","Tumour.origin","Gene.name", "Accession.Number","HGNC.ID", "Gene.Description","Mutation.ID", "Mutation.GRCh37.genome.position","Mutation.Description")
a1<-Annotated[c("Data_Set","ID_sample","ID_tumour","CellLines", "Corrected_Cellline","Primary.site","Site.subtype","Primary.histology","Histology.Subtype","Tumour.origin","Gene.name","Entrez_Gene_ID","Accession.Number","HGNC.ID","RefSeq_mRNA_Id","RefSeq_prot_ID","SwissProt_acc_Id","SwissProt_entry_Id","Gene.Description","Mutation.ID","Chromosome","Start_Position","End_Position","Transcript_Strand","Mutation.GRCh37.genome.position", "cDNA_Change","Protein_Change","Mutation.Description","Updated.ENSEMBL.ID","Reference_Allele","Tumor_Seq_Allele","Reference2Allele")]
Unannotated <- sqldf('SELECT * FROM CCLE_Std EXCEPT SELECT * FROM a')


Annotated_2<-merge(Unannotated,Consequences, by.x= c("End_Position", "Transcript_Strand","Reference2Allele"), 
                   by.y= c("START", "STRAND","ALLELE"), sort= TRUE)

temp<-unique(Annotated_2$End_Position)
temp2<-unique(Unannotated$End_Position)
diff<-setdiff(temp2,temp)
Unannotated2<-Unannotated[(Unannotated$End_Position %in% diff),]

Annotated_3<-merge(Unannotated2,Consequences, by.x= c("Start_Position","Reference2Allele"), 
                   by.y= c("START","ALLELE"), sort= TRUE)


temp<-unique(Annotated_3$Start_Position)
temp2<-unique(Unannotated2$Start_Position)
diff<-setdiff(temp2,temp)
Unannotated_FINAL_CCLE<-Unannotated[(Unannotated$Start_Position %in% diff),]



#Standarization of the CCLE annotations
CCLE_Annotated_Standirized<-Annotated[,c(3,6,7,4,5,8,9,10,11,12,38,16,37,17,15,1,40,39,22,23,24,25,35,2,28,29,36,43,44,45,46,47,48,49,50,51,52,53)]
write.table(CCLE_Annotated_Standirized,"CCLE_1.tsv", sep="\t", quote = F, row.names = F)

CCLE_Annotated_Standarized_2<-Annotated_2[,c(4,7,8,5,6,9,10,11,12,13,38,17,37,18,16,35,40,39,23,24,25,1,2,3,27,28,36,43,44,45,46,47,48,49,50,51,52,53)]
colnames(CCLE_Annotated_Standarized_2)[c(6,7,8,9,10,12,15,16,19,23)]<-c("Primary.site","Site.subtype","Primary.histology",
                                                               "Histology.Subtype","Tumour.origin","HGNC.ID","Accession.Number",
                                                               "Updated.ENSEMBL.ID","Mutation.ID","STRAND")
write.table(CCLE_Annotated_Standarized_2,"CCLE_2.tsv", sep="\t", quote = F, row.names = F)

CCLE_Annotated_Standirized_3<-Annotated_3[,c(3,6,7,4,5,8,9,10,11,12,39,16,38,17,15,36,41,40,22,23,1,24,35,2,27,28,37,44,45,46,47,48,49,50,51,52,53,54)]
colnames(CCLE_Annotated_Standirized_3)[c(6,7,8,9,10,12,15,16,19,23)]<-c("Primary.site","Site.subtype","Primary.histology",
                                                                        "Histology.Subtype","Tumour.origin","HGNC.ID","Accession.Number",
                                                                        "Updated.ENSEMBL.ID","Mutation.ID","STRAND")
write.table(CCLE_Annotated_Standirized_3,"CCLE_3.tsv", sep="\t", quote = F, row.names = F)

CCLE_Annotated<-rbind(CCLE_Annotated_Standirized,CCLE_Annotated_Standarized_2,CCLE_Annotated_Standirized_3)
write.table(CCLE_Annotated,"CCLE_total.tsv", sep="\t", quote = F, row.names = F)


genes <-read.table("genes.tsv", sep="\t", header=TRUE, stringsAsFactors=F)
CCLE_Annotated$FM_PVALUE<-genes$FM_PVALUE[match(CCLE_Annotated$GENE_ID,genes$GENE_ID)]
CCLE_Annotated$FM_QVALUE<-genes$FM_QVALUE[match(CCLE_Annotated$GENE_ID,genes$GENE_ID)]
CCLE_Annotated$CLUST_ZSCORE<-genes$CLUST_ZSCORE[match(CCLE_Annotated$GENE_ID,genes$GENE_ID)]
CCLE_Annotated$CLUST_PVALUE<-genes$CLUST_PVALUE[match(CCLE_Annotated$GENE_ID,genes$GENE_ID)]
CCLE_Annotated$CLUST_QVALUE<-genes$CLUST_QVALUE[match(CCLE_Annotated$GENE_ID,genes$GENE_ID)]


a<-colnames(CCP_Annotated)
colnames(CCLE_Annotated)[c(1:43)]<-a
write.table(CCLE_Annotated,"CCLE_total_ALL.tsv", sep="\t", quote = F, row.names = F)
