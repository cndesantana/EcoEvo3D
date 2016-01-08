dat<-read.table("RichnessPerGen_AnaG_10000.0_cost_0.001_MR_0.0101_VR_0.001.txt",sep=" ",header=TRUE)
ri<-1
posrepl<-which(dat$Real==ri)
nsites<-max(dat$Site)
subdat<-dat[posrepl,]
maxGen<-subdat$G[1]
timeseries<-matrix(1,maxGen+1,nsites)
for(i in 1:maxGen){for(j in 1:nsites){timeseries[i+1,j]<-subdat[j+nsites*(i-1),13];}}

png("./timeseries.png",width=1980,height=1280,res=300);
for(i in 1:nsites){
if(i == 1){plot(timeseries[,i],xlab="Generations",ylab="Beta-Richness",main="Richness per site",type="b",col=i,ylim=c(0,150))}
else{points(timeseries[,i],col=i,type="b")}
}
legend("topleft",legend=c("lake1","lake2","lake3","lake4","lake5","lake6","lake7","lake8"),col=1:8,pch=1);
dev.off();
