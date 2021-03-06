#### For one of the book1
book1_reduce <- book1[, c("visit_id","visitor_id","time","chapter_number","percentage")]
book1.list <- split(book1_reduce, book1_reduce$visitor_id)

mergeData <- function(x){
  tmp.visit_id <- length(x$visit_id)
  tmp.consecutive_session <-  length(unique(x$visit_id))
  #tmp.time <- paste(x$time, collapse = " | ")
  tmp.chapter_number <- length(unique(x$chapter_number))
  tmp.highestchapter_number <- max(unique(x$chapter_number))
  tmp.lowestchapter_number <- min(unique(x$chapter_number))
  #tmp.percentage <- paste(x$percentage, collapse = " | ")
  y <- x[1,]
  y$visit_id <-  tmp.visit_id
  y$consecutive_session <- tmp.consecutive_session
  #y$time <- tmp.time
  y$chapter_number <- tmp.chapter_number
  y$highestchapter_number <- tmp.highestchapter_number
  y$lowestchapter_number <- tmp.lowestchapter_number
  #y$percentage <- tmp.percentage
  return(y)
}

book1.list.result <- lapply(book1.list, function(x) mergeData(x))
book1_result <- do.call(rbind, book1.list.result)
row.names(book1_result)<- NULL

Completed_readers1 <- book1_result[which(book1_result$chapter_number==112 & book1_result$highestchapter_number ==112 & book1_result$lowestchapter_number ==1 ),]
## Just 25 readers were able to complete book

### For book2
book2_reduce <- book2[, c("visit_id","visitor_id","time","chapter_number","percentage")]
book2.list <- split(book2_reduce, book2_reduce$visitor_id)

book2.list.result <- lapply(book2.list, function(x) mergeData(x))
book2_result <- do.call(rbind, book2.list.result)
row.names(book1_result)<- NULL

Completed_readers2 <- book2_result[which(book2_result$chapter_number==32 & book2_result$highestchapter_number ==32 & book2_result$lowestchapter_number ==1 ),]
## Just 25 readers were able to complete book


mean(Completed_readers2$consecutive_session) 
## 3.357

mean(Completed_readers1$consecutive_session)
## 19.24

## Doing 2 sample t-test between two consecutive sessions
t.test(Completed_readers1$consecutive_session, Completed_readers2$consecutive_session)


Welch Two Sample t-test

data:  Completed_readers1$consecutive_session and Completed_readers2$consecutive_session
t = 3.186, df = 24.628, p-value = 0.003891
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
  5.607899 26.157815
sample estimates:
  mean of x mean of y 
19.240000  3.357143 
