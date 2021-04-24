todolist <- list.files(pattern="*.txt")
library(stringr)

substrLeft <- function (x, n=8){
  substr(x, 1, nchar(x)-n)
}
list <- c(1:length(todolist))
mydata<-as.data.frame(matrix(nrow=1, ncol=4)) #creates blank matrix to bind all other data to

for (index in list) {
  filename=as.character(todolist[index])
  input <- read.table(filename,sep="\t")
  name <- substrLeft(filename)
  
  totalreads <- as.character(input[1,1])
  totalreads <- sub(" +.*", "", totalreads)
  
  nummapped <- as.character(input[5,1])
  nummapped <- sub(" +.*", "", nummapped)
  
  percentmapped <- as.character(input[5,1])
  percentmapped <- sub(".*\\(", "", percentmapped)
  percentmapped <- sub("%.*", "", percentmapped)
  
  nameddata <- c(name, totalreads, nummapped, percentmapped)
  mydata <-rbind(mydata, nameddata)
}
colnames(mydata) <- c("Subject", "Number of QC-Passed Reads", "Number of Mapped Reads", "Percentage Mapped")
mydata <- mydata[-1,]
write.csv(mydata,file="flagstat_output.csv")