todolist <- list.files(pattern="*.txt")

substrLeft <- function (x, n=13){
  substr(x, 1, nchar(x)-n)
}
list <- c(1:length(todolist))
mydata<-matrix(nrow=1, ncol=12) #creates blank matrix to bind all other data to

for (index in list) {
  filename=as.character(todolist[index])
  input <- read.table(filename,sep="\t")
  name <- substrLeft(as.character(input[1,3]))
  c1data <- as.character(input[,1])
  namedc1data <- c(name,c1data)
  mydata <-rbind(mydata, namedc1data)
}
colnames(mydata) <- c("Name","Basic Statistics","Per base sequence quality","Per sequence quality scores","Per base sequence content","Per sequence GC content","Per base N content","Sequence Length Distribution","Sequence Duplication Levels","Overrepresented sequences","Adapter Content","Kmer Content")
mydata <- mydata[-1,]
write.csv(mydata,file="fastqc_output.csv")
