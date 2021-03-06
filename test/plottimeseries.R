args <- commandArgs(trailingOnly = TRUE)
inputfile<-args[1];

dat<-read.table(inputfile,sep=" ",header=TRUE)
yrange<-c(0,500)
anaG_<-dat$anaG[1]
cost_<-unlist(strsplit(unlist(strsplit(inputfile,'M'))[1],'_'))[7];    
retG_<-unlist(strsplit(unlist(strsplit(inputfile,'M'))[1],'_'))[5];    
retG_<-as.numeric(retG_)/anaG_; 
ri<-1
posrepl<-which(dat$Real==ri)
nsites<-max(dat$Site)
subdat<-dat[posrepl,]
maxGen<-subdat$G[1]
timeseries<-matrix(1,maxGen+1,nsites)
for(i in 1:maxGen){for(j in 1:nsites){timeseries[i+1,j]<-subdat$alpharich[j+nsites*(i-1)];}}

png(paste("./timeseries_anaG_",anaG_,"_retG_",retG_,"_cost_",cost_,".png",sep=""),width=1980,height=1280,res=300);
for(i in 1:nsites){
if(i == 1){plot(timeseries[,i],xlab="Generations",ylab="Beta-Richness",main=paste("Richness per site\n(cost = ",cost_,", retG = ",retG_,")",sep=""),type="b",col=i,ylim=yrange)}
else{points(timeseries[,i],col=i,type="b")}
}
legend("topleft",legend=c("lake1","lake2","lake3","lake4","lake5","lake6","lake7","lake8"),col=1:8,pch=1);
dev.off();
