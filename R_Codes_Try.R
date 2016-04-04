##################################################################################
############ Finding Metrices that would define the addiction of books ###########
##################################################################################

##################################################################################
########### Doing some exploratory analysis on data ##############################
##################################################################################
bestseller <- read.csv("~/Desktop/bestsellerdata.csv", stringsAsFactors=FALSE)
dim(bestsellerdata)
tail(bestsellerdata)
length(unique(bestsellerdata$visitor_id))
length(unique(bestsellerdata$visit_id))
length(unique(bestsellerdata$event_name))
length(unique(bestsellerdata$time))
length(unique(bestsellerdata$story_id))
length(unique(bestsellerdata$chapter_number))
unique(bestsellerdata$percentage)

## 8589 --> number of visitors
## 17178 --> number of unique visits by these visitors
## 2 --> number of books read by these visitors

##################################################################################
########### Independent Analysis on each book ###################################
################################################################################

book1<-bestsellerdata[bestsellerdata$story_id =="11208",]
row.names(book1)<- NULL
dim (book1)
Visitor_Book1<- as.data.frame(table((book1$visitor_id)))
length(unique(book1$visitor_id))
length(unique(book1$visit_id))

book2<-bestsellerdata[bestsellerdata$story_id =="33946",]
Visitor_Book2<- as.data.frame(table((book2$visitor_id)))
row.names(book2)<- NULL
length(unique(book2$visitor_id))
length (unique(book2$visit_id))


### 8015 user have book1 and 574 have read book2. It is certain that number of visits is also 
## different in each of the cases. It seems 8015 user reads the books in 15671 
## visits while 1157 book visits were made by 574 people.However, it is not certain 
## how many have completed book as that would be one of the yardstick that would
## determine how fasicnating book is.

## testing with few readers (random) to feel how is reading pattern happening
## Case-I
test_subject1 <- book1[(book1$visitor_id == "6c0e7104-66e5-46cc-a943-cf2d51ae4de0"),]
ordered_test_subject1 <- test_subject1[order(test_subject1$chapter_number, test_subject1$percentage),]
test2_subject1 <- book1[(book1$visitor_id == "833d0d8b-a823-4ac0-ada1-dffee73d48c8"),]
ordered_test2_subject1 <- test2_subject1[order(test2_subject1$chapter_number, test2_subject1$percentage),]

## Case-II
test_subject2 <- book2[book2$visitor_id == "33e42755-9c49-40f4-b057-add3c1bc1fc3",]
ordered_test_subject2 <- test_subject2[order(test_subject2$chapter_number, test_subject2$percentage),]
test2_subject2 <- book2[book2$visitor_id == "0bddeb46-8517-48de-83d8-dd01be764971",]
ordered_test2_subject2 <- test2_subject2[order(test2_subject2$chapter_number, test2_subject2$percentage),]

## We observe the reading patterns of people are quite random, in the earlier case 
## only the reader only read last two chapters  while in the later case the reader
## all of the chapters for book1 while in case 2 the readers read upto chapter 4 
## and 6 respectively so it presnts the reading events are quite random

## So the question would arise is how many individuals are there have actually
##  completed the whole book and does that pattern show some indication of addictiveness 
## of the book.

#################################################################################
##############################################################################
range(unique(book1$chapter_number)) ## 1-112
range(unique(book2$chapter_number)) ## 1-32

## Finding out number of chapeters in each of book. Obiviously, book2 is small book
## However, it is unclear whether that would hamper the reading habits. Lets find
## out how many of the above unique readers were able to reach to the end of the book
## starting from chapter 1.

Completedreaders1<-Visitor_Book1[which(Visitor_Book1$Freq > 448),]
row.names(Completedreaders1)= NULL
Completedreaders2<-Visitor_Book2[which(Visitor_Book2$Freq > 128),]
row.names(Completedreaders2)= NULL

Completed_readers1<-book1[which(book1$visitor_id %in% Completedreaders1$Var1),] 
length(unique(Completed_readers1$visitor_id))
row.names(Completed_readers1) = NULL

Completed_readers2<-book2[which(book2$visitor_id %in% Completedreaders2$Var1),] 
length(unique(Completed_readers2$visitor_id))
row.names(Completed_readers2) = NULL

## Here we can see that only 125 readers are found to read all chapters in book1 
## while 62 readers claim to reach chapter all of chapters in book 2. 
## Conidering the original 8015 and 574 unique reader for book1 and book2, only 
## 1.5% complete book1 while 10.8% completed book2. So is this only due to chapter 
## length or some other factors coming to play as well. Furthermore completing book
## is not only the indicator that explains the addictiveness of book. We will delve
## more into it now.

mergeData <- function(x){
  tmp.Accession <- paste(x$Accession.Number, collapse = " | ")
  tmp.Mutation.CDS <- paste(x$Mutation.CDS, collapse = " | ")
  tmp.Mutation.AA <- paste(x$Mutation.AA, collapse = " | ")
  y <- x[1,]
  y$Accession.Number <-  tmp.Accession
  y$Mutation.CDS <- tmp.Mutation.CDS
  y$Mutation.AA <- tmp.Mutation.AA
  return(y)
}



