#library(ggplot2)
#library(plyr)
#library(reshape2)

error.bar <- function(x, y, upper, mycolor, lower=upper, length=0.1,...){
  if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
    stop("vectors must be same length")
  arrows(x,y+upper, x, y-lower, col=mycolor, angle=90, code=3, length=length, ...)
}


printfigures<-function(inputfile1,inputfile2,nsites,sufix,verticesdatafile){

    cat("Plotting...\n");
    dat<-read.table(inputfile1,header=TRUE);#average
    dat2<-read.table(inputfile2,header=TRUE);#sd    
    yrange<-range(0,350)*1.2;
    xrange<-c(0.05,1.25);
    anaG<-dat$anaG[1];    
    cost<-unlist(strsplit(unlist(strsplit(sufix,'M'))[1],'_'))[7];    
    retG<-unlist(strsplit(unlist(strsplit(sufix,'M'))[1],'_'))[5];   
    retG<-as.numeric(retG)/anaG; 
    nrep<-max(dat$Real);
    lakesarea<-dat$lakeArea;
    lakesisolation<-dat$dT;    
    figure1<-paste("Richness_Isolation_",sufix,".png",sep=""); 
    png(figure1,width=1980,height=1280,res=300);
    plot(lakesisolation,dat$alpharich,xlab="Isolation level",ylab=expression(paste(alpha,"-rich",sep="")),pch=19,main=paste("Isolation x Richness\ncost = ",cost,", retG = ",retG,sep=""),lwd=10,log="x",ylim=yrange,xlim=xrange);
    error.bar(x = lakesisolation, y = dat$alpharich, upper = dat2$alpharich, lower = dat2$alpharich, mycolor = "black")
    dev.off();
    
    figure2<-paste("Richness_Size_Isolation_",sufix,".png",sep=""); 
    png(figure2,width=1980,height=1280,res=300);
    mindt<-min(dat$dT);#minimal isolation
    maxdt<-max(dat$dT);#maximal isolation
    plot(lakesarea,dat$alpharich,xlab="Area of the lakes (km^2)",ylab=expression(paste(alpha,"-rich",sep="")),pch=19,main=paste("Size of lakes x Richness x Isolation\ncost = ",cost,", retG = ",retG,sep=""),lwd=10,log="x",cex=(1.51-((dat$dT-mindt)/(maxdt-mindt))),ylim=yrange,xlim=range(lakesarea));
    error.bar(x = lakesarea, y = dat$alpharich, upper = dat2$alpharich, lower = dat2$alpharich, mycolor = "black")
    legend("topleft",legend=c(signif(mindt,3),signif(maxdt,3)),pch=c(19,19),pt.cex=c(2.502,1.502),title="Isolation level");
    dev.off();

    figure3<-paste("Speciation_LakesPerArea_",sufix,".png",sep=""); 
    png(figure3,width=1980,height=1280,res=300);
    lakesnames<-paste("lake",1:nsites,sep="");#we should replace this list by the real names of the lakes.
    lakesarea<-dat$lakeArea;
    orderbyarea<-order(lakesarea);
    mp<-barplot(rbind((dat$SpecANA/dat$alpharich)[orderbyarea],(dat$SpecCLA/dat$alpharich)[orderbyarea],(dat$SpecMR/dat$alpharich)[orderbyarea],(dat$DispersalRich/dat$alpharich)[orderbyarea]),main=paste("Origin of species\ncost = ",cost,", retG = ",retG,sep=""),xlab="Area of the lakes (Km^2)",ylab="Proportion of Species",col=c("green","darkblue","yellow","red"),ylim=c(0,1));
    text(mp, par("usr")[3] - 0.01, srt = 45, adj = 1,labels = signif(lakesarea[orderbyarea],2),xpd = TRUE, font = 2)
    dev.off();

    figure4<-paste("Speciation_LakesPerIsolation_",sufix,".png",sep=""); 
    png(figure4,width=1980,height=1280,res=300);
    lakesnames<-paste("lake",1:nsites,sep="");#we should replace this list by the real names of the lakes.
    lakesisolation<-dat$dT;
    orderbyisolation<-order(lakesisolation);
    mp<-barplot(rbind((dat$SpecANA/dat$alpharich)[orderbyisolation],(dat$SpecCLA/dat$alpharich)[orderbyisolation],(dat$SpecMR/dat$alpharich)[orderbyisolation],(dat$DispersalRich/dat$alpharich)[orderbyisolation]),main=paste("Origin of species\ncost = ",cost,", retG = ",retG,sep=""),xlab="Isolation of the lakes",ylab="Proportion of Species",col=c("green","darkblue","yellow","red"),ylim=c(0,1));
    text(mp, par("usr")[3] - 0.01, srt = 45, adj = 1,labels = signif(lakesisolation[orderbyisolation],2),xpd = TRUE, font = 2)
    dev.off();
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
#verticesdatafile<-"./VerticesData.txt"
#printfigures("./outputav_RichnessPerSite_AnaG_500000.0_cost_0.001_MR_7.99e-5_VR_5.11e-5.txt","./outputsd_RichnessPerSite_AnaG_500000.0_cost_0.001_MR_7.99e-5_VR_5.11e-5.txt","RichnessPerSite_AnaG_500000.0_cost_0.001_MR_7.99e-5_VR_5.11e-5.txt",8,"./VerticesData.txt")
printfigures(preprocres$file1,preprocres$file2,preprocres$nsites,preprocres$sufix,verticesdatafile);
