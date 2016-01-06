args <- commandArgs(trailingOnly = TRUE)
inputfile1<-args[1];
inputfile2<-args[2];

dat000<-read.csv(inputfile1,sep=" ",row.names=NULL);
dat001<-read.csv(inputfile2,sep=" ",row.names=NULL);
MAX<-10000;

png("./steadystate.png",width=1980,height=1280,res=300);
plot(dat001$Gi[1:MAX],dat001$gamma[1:MAX]+1,xlab="Generations",ylab="Richness",col=1,log="y")
points(dat000$Gi[1:MAX],dat000$gamma[1:MAX]+1,col=2)
legend("bottomright",legend=c("cost = 0","cost = 0.001"),col=c(1,2),pch=1,lwd=1,title="costs");
dev.off();
