# Project: Mutation Annotation of Cosmic Cell-Line Project (CCP), Sanger and Cancer 
# Cell-Line Encyclopedia

# Part1: Reading the files in R and doing the usual check ups
CCP<- read.table("CosmicCLP_CompleteExport_v68_070214.tsv", sep="\t", header=TRUE, stringsAsFactors=F)
CCLE<-read.delim("CCLE.maf",sep="\t", header=TRUE, stringsAsFactors=F)

write.table(CCLE,"CCLE.tsv", sep="\t")

# While reading the files the stringasfactor was taken as false so that it would help in the 
# downstreaming analysis of the data. Also we used the read.delim command as the data for the 
# CCLE data set was .maf files.

colnames(CCP)
colnames(CCLE)
head(CCP)
head(CCLE)

#Part2: Standarizing the Celllines names in both datasets
Cellline_CCP<-as.character(CCP[,"Sample.name"])
Cellline_CCLE<- as.character(CCLE[,"Tumor_Sample_Barcode"])

#Using the library reshape to split the both dataframes
library(reshape2)
temp<- colsplit(Cellline_CCLE,"_", c("Celllines","Celltype"))
CCLE_Cellline<-as.character(temp$Celllines)
CCLE_Cellline<-toupper(CCLE_Cellline)
CCLE_Cellline<-gsub("[^[:alnum:]///' ]", "", CCLE_Cellline)
head(CCLE_Cellline)

CCLE$CellLines<-CCLE_Cellline
# Adding the new transformed cellline column in the original Cosmic Cell-lines dataframe
CCLE$Tumor_Sample_Barcode<-NULL
# Removing the original cellline data
CCLE$Celltype_Sample<-temp[,2]
#Adding the sample information in the original dataset

temp<-colsplit(Cellline_CCP,"-", c("part1","part2","part3"))
CCP_Cellline<-as.character(paste(temp[,1],temp[,2],temp[,3],sep=""))
CCP_Cellline<-toupper(CCP_Cellline)
CCP_Cellline<-gsub("[^[:alnum:]///' ]", "", CCP_Cellline)
CCP_Cellline<-gsub(" ", "",CCP_Cellline)
head(CCP_Cellline)

CCP$CellLines<-CCP_Cellline
# Adding the new transformed cellline column in the original Cosmic Cell-lines dataframe
CCP$Sample.name<-NULL
# Removing the original cellline data

remove(Cellline_CCP, Cellline_CCLE)
remove(temp)

cellline_lookup_table <- read.table("cellline_lookup_table.txt",sep="\t", header=TRUE, stringsAsFactors=F)
a <-cellline_lookup_table$STD_Name[match(CCP_Cellline, cellline_lookup_table$cell.line)]
CCP$Corrected_Cellline<-CCP_Cellline
CCP$Corrected_Cellline[!is.na(a)]<-as.character(a[!is.na(a)])
dummy<-CCP[1:1000,]

b<-cellline_lookup_table$STD_Name[match(CCLE_Cellline, cellline_lookup_table$cell.line)]
CCLE$Corrected_Cellline<-CCLE_Cellline
CCLE$Corrected_Cellline[!is.na(b)]<-as.character(b[!is.na(b)])
dummy<-CCLE[1:1000,]

remove(dummy)

## Standarizing and homogenizing both data sets such that they have the necessay information
## Making the 2 data sets homogenous such that both includes the essential variables 
# as have been asked in the project

CCLE_Std <-CCLE[,c(1,2,5,6,7,9,10,34,35,37,39,40,41,42,43,51,52,53)]
CCP_Std<-CCP[,c(1,2,3,4,5,6,7,8,9,11,12,13,14,18,19,22,23,25,26)]

CCLE_Std$Mutation.Description<-paste(CCLE_Std$Variant_Classification,CCLE_Std$Variant_Type, sep="-")
# Since the mutation description in the Cosmic Cell- line data sets is in 2 columns , both 
# are merged into one Combined column as "Mutation. Description"
CCLE_Std$Variant_Classification<- NULL
CCLE_Std$Variant_Type<- NULL

## In the Sanger data set the chromosome number and genomic position are provided as a combined aggregate
## form. Hence these are seperated and given a seperate column in the data set
CCP_Std$Chromosome <- sapply(strsplit(CCP_Std$Mutation.GRCh37.genome.position, ":"), "[", 1)
Position<- sapply(strsplit(CCP_Std$Mutation.GRCh37.genome.position, ":"), "[", 2)
CCP_Std$Start_Position <- sapply(strsplit(Position, "-"), "[", 1)
CCP_Std$End_Position <- sapply(strsplit(Position, "-"), "[", 2)
CCP_Std$Sample.source<- NULL

CCLE_Std$temp<-paste(CCLE_Std$Start_position,CCLE_Std$End_position,sep="-")
CCLE_Std$Mutation.GRCh37.genome.position<-paste(CCLE_Std$Chromosome,CCLE_Std$temp,sep=":")
CCLE_Std$temp<- NULL

#Assigning multiple columns to the both dataframes to get the homogenous columns in the both data sets
CCLE_Std[c("Accession.Number", "HGNC.ID", "ID_Sampe", "ID_Tumor", "Site.subtype", "Primary.histology", "Histology.Subtype","Mutation.ID","Tumour.origin")]<- NA
CCP_Std[c("Entrez_Gene_ID", "RefSeq_mRNA_Id", "RefSeq_prot_ID", "SwissProt_acc_Id", "SwissProt_entry_Id","Gene.Description")]<- NA
CCLE_Std$Data_Set<- "Broad"
CCP_Std$Data_Set<- "Sanger"

##Reordering the whole column 
CCLE_Std<-CCLE_Std[,c("Data_Set","ID_Sampe","ID_Tumor","CellLines", "Corrected_Cellline", "Celltype_Sample","Site.subtype","Primary.histology","Histology.Subtype","Tumour.origin","Hugo_Symbol","Entrez_Gene_Id","Accession.Number","HGNC.ID","Refseq_mRNA_Id","Refseq_prot_Id", "SwissProt_acc_Id", "SwissProt_entry_Id","Description","Mutation.ID","Chromosome","Start_position","End_position","Transcript_Strand", "Mutation.GRCh37.genome.position","cDNA_Change","Protein_Change", "Mutation.Description")]
CCP_Std<-CCP_Std[,c("Data_Set","ID_sample","ID_tumour","CellLines" ,  "Corrected_Cellline", "Primary.site" ,"Site.subtype","Primary.histology", "Histology.subtype","Tumour.origin","Gene.name", "Entrez_Gene_ID" ,"Accession.Number","HGNC.ID", "RefSeq_mRNA_Id","RefSeq_prot_ID" , "SwissProt_acc_Id",  "SwissProt_entry_Id" , "Gene.Description","Mutation.ID","Chromosome" ,"Start_Position","End_Position", "Mutation.GRCh37.strand", "Mutation.GRCh37.genome.position","Mutation.CDS","Mutation.AA", "Mutation.Description")]
colnames(CCLE_Std)[c(2,3,6,11,12,15,16,19,22,23)]<-c("ID_sample","ID_tumour","Primary.site","Gene.name","Entrez_Gene_ID","RefSeq_mRNA_Id","RefSeq_prot_ID","Gene.Description","Start_Position","End_Position")
colnames(CCP_Std)[c(24,26,27)]<-c("Transcript_Strand","cDNA_Change","Protein_Change")

##Saving the files to the external source for further use 
write.table(CCLE_Std, "Broad data set.txt", sep="\t") 
write.table(CCP_Std,"Sanger data set.txt",  sep="\t")
write.table(CCLE, "Coarse Broad data set.txt", sep="\t") 
write.table(CCP,"Coarse Sanger data set.txt",  sep="\t")


## Part3: Finding the common Celllines between two cancer cell-lines project
length(unique(CCP_Std$Corrected_Cellline)) #1008
length(unique(CCLE_Std$Corrected_Cellline)) # 904
Common_Celllines<-intersect((unique(CCP_Std$Corrected_Cellline)),(unique(CCLE_Std$Corrected_Cellline)))
length(Common_Celllines) #594

## Using of the match function to get the common and unique celllines data from the both data frames
Common_Celllines_CCP<-CCP_Std[(CCP_Std$Corrected_Cellline %in% Common_Celllines),] 
Unique_Cellline_CCP<-CCP_Std[(!(CCP_Std$Corrected_Cellline %in% Common_Celllines)),]
Common_Celllines_CCLE<-CCLE_Std[(CCLE_Std$Corrected_Cellline %in% Common_Celllines),] 
Unique_Cellline_CCLE<-CCLE_Std[(!(CCLE_Std$Corrected_Cellline %in% Common_Celllines)),]

write.table(Unique_Cellline_CCLE, "Unique cell line mutation_Board.txt", sep="\t") 
write.table(Unique_Cellline_CCP,"Unique Cell line mutation_Sanger.txt",  sep="\t")


#Part:4 Finding the common genomic postition mutation in common cell lines
Mutonly_Com_CCP<-Common_Celllines_CCP[which(Common_Celllines_CCP$Mutation.GRCh37.genome.position != ""),]
##Mutonly_Com_CCLE<-Common_Celllines_CCLE[which(Common_Celllines_CCLE$Mutation.GRCh37.genome.position != ""),]
Rest_Com_CCP<-Common_Celllines_CCP[which(Common_Celllines_CCP$Mutation.GRCh37.genome.position == ""),]
write.table(Rest_Com_CCP,"Rest Common CL Sanger data set.txt",  sep="\t")

## IMP POINT at this point we selected the cell-lines that had the mutation ID and  genome position in the Sanger Cosmic databases.
## Rest of the data can be found as the  complement of the above code

## doing the same thing for the unique cell lines found in the sanger dataset
Mutonly_Unique_CCP<-Unique_Cellline_CCP[which(Unique_Cellline_CCP$Mutation.GRCh37.genome.position != ""),]
## Mutonly_Unique_CCLE<-Unique_Cellline_CCLE[which(Unique_Cellline_CCLE$Mutation.GRCh37.genome.position != ""),]
Rest_Unique_CCP<-Unique_Cellline_CCP[which(Unique_Cellline_CCP$Mutation.GRCh37.genome.position == ""),]
write.table(Mutonly_Unique_CCP, "Mutation Unique CL Sanger data set.txt", sep="\t") 
write.table(Rest_Unique_CCP,"Rest Unique  CL Sanger data set.txt",  sep="\t")

##length((intersect(Common_Celllines_CCLE$Mutation.GRCh37.genome.position,Mutonly_Com_CCP$Mutation.GRCh37.genome.position)) #30799
##common_Genomic_Variants<-intersect(Common_Celllines_CCLE$Mutation.GRCh37.genome.position,Mutonly_Com_CCP$Mutation.GRCh37.genome.position)
##Mutation_Common<- Mutonly_Com_CCP[match(Mutonly_Com_CCP$Mutation.GRCh37.genome.position,common_Genomic_Variants)]
##x <- Mutonly_Com_CCP[which(Mutonly_Com_CCP$Mutation.GRCh37.genome.position %in% Common_Genomic_Variants),]
##y <- Common_Celllines_CCLE[which(Common_Celllines_CCLE$Mutation.GRCh37.genome.position %in% Common_Genomic_Variants),]

##head(x$Mutation.GRCh37.genome.position)

# Part 5: Making the 2 data sets homogenous such that both includes the essential variables 
# as have been asked in the project
# Creating a unique identifier for the each row of the dataset by fusing the Corrected Cell-lines and Genomic position 

Mutonly_Com_CCP$Unique_identifier<-as.character(paste(Mutonly_Com_CCP$Corrected_Cellline,Mutonly_Com_CCP$Chromosome, Mutonly_Com_CCP$Start_Position,Mutonly_Com_CCP$End_Position, sep =""))
Common_Celllines_CCLE$Unique_identifier<-as.character(paste(Common_Celllines_CCLE$Corrected_Cellline,Common_Celllines_CCLE$Chromosome, Common_Celllines_CCLE$Start_Position,Common_Celllines_CCLE$End_Position, sep =""))
a<- unique(Mutonly_Com_CCP$Unique_identifier)
b<- unique(Common_Celllines_CCLE$Unique_identifier)

# Breaking down of the large data files into smaller chunks of 100,000 rows in order to remove redundancy in the rows and getting unique rows
c<-a[400001:599573]
d<-Mutonly_Com_CCP[(Mutonly_Com_CCP$Unique_identifier%in% c),]

## Used the merge,r function from the other function files.
##Sys.time()
d.list<-split(d,d$Unique_identifier)
##Sys.time()

source("merge.r")
d.result <- lapply(d.list, function(x) mergeData(x))
require(data.table)
m.1<- rbindlist (d.result)
head(m.1)
remove(d,d.list,d.result)


## binding all the elemets of the broken down dataset
Com_CL_Cos<-rbind(i.1,j.1,k.1,l.1,m.1)
remove(i.1,j.1,k.1,l.1,m.1)
Com_CL_Cos$Gene.name<-sapply(strsplit(Com_CL_Cos$Gene.name, "_"), "[", 1)
write.table(Com_CL_Cos, "Common_CL_Sanger dataset.txt", sep="\t") 
 
# gives the number of rows in the final common cellline-common-genomic position
c<-intersect(a,b) #32443 

#Doing the same thing for CCLE dataset
d.list<-split(Common_Celllines_CCLE,Common_Celllines_CCLE$Unique_identifier)
source("merge2.r")
d.result <- lapply(d.list, function(x) mergeData2(x))
require(data.table)
i.2<- rbindlist (d.result)
a<-i.2[(i.2$Unique_identifier=="TT11128638217128638217"),]
## checking whether the merge has worked or not
write.table(Com_CL_Broad, "Common_CL_Broad dataset.txt", sep="\t")
Com_CL_Broad<-i.2
remove(i.2,a,b)


## Fusing the both datasets into common core data
Total_mutation_set<- rbind( Com_CL_Cos,Com_CL_Broad)

## Done upto this point...

temp<-Total_mutation_set[(Total_mutation_set$Unique_identifier%in% c),]
d<- c[20001:32443]
test<-temp[(temp$Unique_identifier%in% d),]
d.list<-split(test,test$Unique_identifier)
source("merge3.r")
d.result<- lapply(d.list, function(x) mergeData3(x))
require(data.table)
Common_gene_variant_Both.3<-rbindlist(d.result)
remove(test,d,d.list,d.result)
remove(temp)

Common_Celline_genomic_variants<-rbind(Common_gene_variant_Both.1,Common_gene_variant_Both.2,Common_gene_variant_Both.3)
remove(Common_gene_variant_Both.1,Common_gene_variant_Both.2,Common_gene_variant_Both.3)

a<-Com_CL_Broad[(!(Com_CL_Broad$Unique_identifier %in% c)),]
b<-Com_CL_Cos[(!(Com_CL_Cos$Unique_identifier %in% c)),]
Sanger_Both_Broad_dataset<-rbind(b,Common_Celline_genomic_variants,a)
write.csv(Sanger_Both_Broad_dataset,"Sanger Broad set.txt", sep="\t")
