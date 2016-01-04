dat<-read.csv("./outputav.dat",sep=" ");

png("Speciation_Lakes.png");
lakesnames<-c("lake1","lake2","lake3","lake4","lake5","lake6","lake7","lake8");
lakessizes<-dat$Ji;
orderbysize<-order(dat$Ji);
lakesnames<-lakesnames[orderbysize];
lakessizes<-lakessizes[orderbysize];
mp<-barplot(rbind(dat$SpecANA[orderbysize],dat$SpecCLA[orderbysize],dat$SpecMR[orderbysize],dat$DispersalRich[orderbysize]),main="Origin of species",xlab="Size of the lakes",ylab="Species",col=c("green","darkblue","yellow","red"),legend=c("ANA","CLA","Regional-Mig","Local-Mig"),args.legend=list(x = "topleft"));
text(mp, par("usr")[3] - 0.025, srt = 45, adj = 1,labels = lakessizes,xpd = TRUE, font = 2)
dev.off();

dat<-read.csv("./RichnessPerSite.txt",sep=" ");

png("Richness_Isolation.png");
plot(dat$dT,dat$alpharich,xlab="Isolation level",ylab=expression(paste(alpha,"-rich",sep="")),pch=19,main="Isolation x Richness",lwd=10,log="x");
dev.off();

png("Richness_Size_Isolation.png");
plot(dat$Ji,dat$alpharich,xlab="Size of the lakes",ylab=expression(paste(alpha,"-rich",sep="")),pch=19,main="Size of lakes x Richness x Isolation",lwd=10,log="x",cex=2*(1-dat$dT));
legend("topleft",legend=c(mindt,maxdt),pch=c(19,19),lwd=c(4*(1-mindt),4*(1-maxdt)),title="Isolation level");
dev.off();

