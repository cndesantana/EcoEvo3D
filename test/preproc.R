dat<-read.csv("./RichnessPerSite.txt",sep=" ")
average<-array(0,dim=c(8,length(names(dat))))
standarddev<-array(0,dim=c(8,length(names(dat))))
nvar<-length(names(dat))
nsite<-8

cat(names(dat),sep=" ",file="./outputsd.dat")
cat(names(dat),sep=" ",file="./outputav.dat")
cat("\n",file="./outputsd.dat",append=TRUE)
cat("\n",file="./outputav.dat",append=TRUE)
for (i in 1:nsite){
	for (j in 1:nvar){
		average[i,j]<-mean(dat[i+8*(0:9),j])
		standarddev[i,j]<-sd(dat[i+8*(0:9),j])
		cat(paste(average[i,j]," ",sep=""),sep=" ",file="./outputav.dat",append=TRUE);
		cat(paste(standarddev[i,j]," ",sep=""),sep=" ",file="./outputsd.dat",append=TRUE);
	}
	cat("\n",file="./outputav.dat",append=TRUE);
	cat("\n",file="./outputsd.dat",append=TRUE);
}
