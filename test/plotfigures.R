library(gdata)
library(ggplot2)
library(plyr)
library(reshape2)

error.bar <- function(x, y, upper, mycolor, lower=upper, length=0.1,...){
  if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
    stop("vectors must be same length")
  arrows(x,y+upper, x, y-lower, col=mycolor, angle=90, code=3, length=length, ...)
}

readAreaOfLakes<-function(verticesdatafile){
    dat<-read.csv(verticesdatafile,sep=" ");
    return(dat[,3]);    
}

readVolumeOfLakes<-function(verticesdatafile){
    dat<-read.csv(verticesdatafile,sep=" ");
    return(dat[,2]);    
}


printfigures<-function(inputfile1,inputfile2,inputfile3,nsites,sufix,verticesdatafile){

    cat("Plotting...\n");
    yrange<-c(0,70);
    xrange<-c(0.05,1.25);
    cost<-unlist(strsplit(unlist(strsplit(sufix,'M'))[1],'_'))[5];    
    dat<-read.csv(inputfile1,sep=" ");
    anaG<-dat$anaG[1];    
    nrep<-max(dat$Real);
    figure1<-paste("Richness_Isolation_",sufix,".png",sep=""); 
    png(figure1,width=1980,height=1280,res=300);
    plot(dat$dT,dat$alpharich,xlab="Isolation level",ylab=expression(paste(alpha,"-rich",sep="")),pch=19,main=paste("Isolation x Richness\ncost = ",cost,", anaG = ",anaG,sep=""),lwd=10,log="x",ylim=yrange,xlim=xrange);
    dev.off();
    
    figure2<-paste("Richness_Size_Isolation_",sufix,".png",sep=""); 
    png(figure2,width=1980,height=1280,res=300);
    absolutelakesareas<-readAreaOfLakes(verticesdatafile);#area of the lakes
    absolutelakesvolumes<-readVolumeOfLakes(verticesdatafile);#volume of the lakes
    lakesareas<-rep(absolutelakesareas[order(absolutelakesvolumes[order(unique(dat$Ji))])],nrep)
    mindt<-min(dat$dT);#minimal isolation
    maxdt<-max(dat$dT);#maximal isolation
    plot(lakesareas,dat$alpharich,xlab="Area of the lakes (km^2)",ylab=expression(paste(alpha,"-rich",sep="")),pch=19,main=paste("Size of lakes x Richness x Isolation\ncost = ",cost,", anaG = ",anaG,sep=""),lwd=10,log="x",cex=2*(1.001-((dat$dT-mindt)/(maxdt-mindt) )),ylim=yrange);
    legend("topleft",legend=c(signif(mindt,3),signif(maxdt,3)),pch=c(19,19),pt.cex=c(3.002,1.002),title="Isolation level");
#    text(lakesareas, par("usr")[3] - 0.025, srt = 45, adj = 1,labels = signif(lakesareas,2),xpd = TRUE, font = 2)
    dev.off();

    dat<-read.table(inputfile2,header=TRUE);#average
    figure3<-paste("Speciation_Lakes_",sufix,".png",sep=""); 
    png(figure3,width=1980,height=1280,res=300);
    lakesnames<-paste("lake",1:nsites,sep="");#we should replace this list by the real names of the lakes.
#    absolutelakesareas<-readAreaOfLakes(verticesdatafile);#sizes of the lakes
#    relativelakesareas<-(dat$Ji/dat$J);
#    sumlakesareas<-sum(absolutelakesareas);
#    lakesareas<-absolutelakesareas*sumlakesareas;
    lakesareas<-readAreaOfLakes(verticesdatafile);#areas of the lakes (one value per lake)
    orderbysize<-order(lakesareas);
    lakesnames<-lakesnames[orderbysize];
    mp<-barplot(rbind(dat$SpecANA[orderbysize],dat$SpecCLA[orderbysize],dat$SpecMR[orderbysize],dat$DispersalRich[orderbysize]),main=paste("Origin of species\ncost = ",cost,", anaG = ",anaG,sep=""),xlab="Size of the lakes (Km^2)",ylab="Species",col=c("green","darkblue","yellow","red"),legend=c("ANA","CLA","Regional-Mig","Local-Mig"),args.legend=list(x = "topleft"),ylim=yrange);
    text(mp, par("usr")[3] - 0.75, srt = 45, adj = 1,labels = signif(lakesareas,2),xpd = TRUE, font = 2)
    dev.off();

#    dat1<-read.table("./outputav_RichnessPerSite_AnaG_20000_cost_0.0_MR_3.2e-5_VR_6.09e-6.txt",header=TRUE)
#    dat2<-read.table("./outputsd_RichnessPerSite_AnaG_20000_cost_0.0_MR_3.2e-5_VR_6.09e-6.txt",header=TRUE)
#    dat3<-read.table("./outputav_RichnessPerSite_AnaG_20000_cost_0.001_MR_6.93e-5_VR_2.04e-6.txt",header=TRUE)
#    dat4<-read.table("./outputsd_RichnessPerSite_AnaG_20000_cost_0.001_MR_6.93e-5_VR_2.04e-6.txt",header=TRUE)
#
#    yy1<-dat1$SpecANA;
#    ee1<-dat2$SpecANA;
#    yy2<-dat3$SpecANA;
#    ee2<-dat4$SpecANA;
#    png("./anagenesis.png",width=1980,height=1280,res=300)
#    plot(lakesareas,yy1,xlab="Area of the lakes (km^2)",ylab=expression(paste(beta,"-rich",sep="")),pch=19,lwd=10,log="x",cex=2*(1.001-((dat$dT-mindt)/(maxdt-mindt) )),ylim=c(0,8),main="Anagenesis Speciation");
#    error.bar(lakesareas,yy1,ee1,"black")
#    par(new=T)
#    plot(lakesareas,yy2,xlab="Area of the lakes (km^2)",ylab=expression(paste(beta,"-rich",sep="")),pch=19,lwd=10,log="x",cex=2*(1.001-((dat$dT-mindt)/(maxdt-mindt) )),ylim=c(0,8),col=2);
#    error.bar(lakesareas,yy2,ee2,"red")
#    legend("topleft",legend=c("cost = 0","cost = 0.001"),col=c(1,2),title="costs",pch=c(19,19),pt.cex=2.002);
#    dev.off();
#
#    yy1<-dat1$DispersalRich;
#    ee1<-dat2$DispersalRich;
#    yy2<-dat3$DispersalRich;
#    ee2<-dat4$DispersalRich;
#    png("./dispersal.png",width=1980,height=1280,res=300)
#    plot(lakesareas,yy1,xlab="Area of the lakes (km^2)",ylab=expression(paste(beta,"-rich",sep="")),pch=19,lwd=10,log="x",cex=2*(1.001-((dat$dT-mindt)/(maxdt-mindt) )),ylim=c(0,70),main="Local Migration");
#    error.bar(lakesareas,yy1,ee1,"black")
#    par(new=T)
#    plot(lakesareas,yy2,xlab="Area of the lakes (km^2)",ylab=expression(paste(beta,"-rich",sep="")),pch=19,lwd=10,log="x",cex=2*(1.001-((dat$dT-mindt)/(maxdt-mindt) )),ylim=c(0,70),col=2);
#    error.bar(lakesareas,yy2,ee2,"red")
#    legend("topleft",legend=c("cost = 0","cost = 0.001"),col=c(1,2),title="costs",pch=c(19,19),pt.cex=2.002);
#    dev.off();
}


preprocfunc<-function(inputfile){
    cat("Pre-processing...\n");
    dat<-read.csv(inputfile,sep=" ")
    nsites<-max(dat$Site);
    nrep<-max(dat$Real);
    averageval<-array(0,dim=c(8,length(names(dat))))
    standarddevval<-array(0,dim=c(8,length(names(dat))))
    nvar<-length(names(dat))
    sufix<-ifelse(is.na(as.character(unlist(strsplit(inputfile,"/")))[2]),inputfile,as.character(unlist(strsplit(inputfile,"/")))[2]);
    averageoutputfile<-paste("./outputav_",sufix,sep="") 
    sdoutputfile<-paste("./outputsd_",sufix,sep="") 

    cat(names(dat),sep=" ",file=averageoutputfile)
    cat(names(dat),sep=" ",file=sdoutputfile)
    cat("\n",file=averageoutputfile,append=TRUE)
    cat("\n",file=sdoutputfile,append=TRUE)
    for (j in 1:nvar){
        for (i in 1:nsites){
                possites<-(0:(nrep-1))*nsites+i
    		averageval[i,j]<-mean(dat[possites,j],na.rm=TRUE)
    		standarddevval[i,j]<-sd(dat[possites,j],na.rm=TRUE)
    	}
    }
    cat("Printing the output files...\n");
    for(i in 1:nsites){
        for(j in 1:nvar){
            cat(paste(averageval[i,j]," ",sep=""),sep=" ",file=averageoutputfile,append=TRUE);
            cat(paste(standarddevval[i,j]," ",sep=""),sep=" ",file=sdoutputfile,append=TRUE);
    	}
    	cat("\n",file=averageoutputfile,append=TRUE);
    	cat("\n",file=sdoutputfile,append=TRUE);
    }

    return(list(file1=averageoutputfile,file2=sdoutputfile,nsites=nsites,sufix=sufix));
}

args <- commandArgs(trailingOnly = TRUE)
inputfile<-args[1];
verticesdatafile<-args[2];

preprocres<-preprocfunc(inputfile);
printfigures(inputfile,preprocres$file1,preprocres$file2,preprocres$nsites,preprocres$sufix,verticesdatafile);
