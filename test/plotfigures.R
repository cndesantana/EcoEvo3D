readSizeOfLakes<-function(verticesdatafile){
    dat<-read.csv(verticesdatafile,sep=" ",header=FALSE);
    return(dat[,2]);    
}

printfigures<-function(inputfile1,inputfile2,nsites,sufix,verticesdatafile){

    dat<-read.csv(inputfile1,sep=" ");
    figure1<-paste("Richness_Isolation_",sufix,".png",sep=""); 
    png(figure1,width=1980,height=1280,res=300);
    plot(dat$dT,dat$alpharich,xlab="Isolation level",ylab=expression(paste(alpha,"-rich",sep="")),pch=19,main="Isolation x Richness",lwd=10,log="x");
    dev.off();
    
    figure2<-paste("Richness_Size_Isolation_",sufix,".png",sep=""); 
    png(figure2,width=1980,height=1280,res=300);
    absolutelakessizes<-readSizeOfLakes(verticesdatafile);#sizes of the lakes
    relativelakessizes<-(dat$Ji/dat$J);
    sumlakessizes<-sum(absolutelakessizes);
    lakessizes<-relativelakessizes*sumlakessizes;
    plot(lakessizes,dat$alpharich,xlab="Size of the lakes (km^2)",ylab=expression(paste(alpha,"-rich",sep="")),pch=19,main="Size of lakes x Richness x Isolation",lwd=10,log="x",cex=2*(1-dat$dT));
    mindt<-min(dat$dT);
    maxdt<-max(dat$dT);
    legend("topleft",legend=c(signif(mindt,3),signif(maxdt,3)),pch=c(19,19),lwd=c(4*(1-mindt),4*(1-maxdt)),title="Isolation level");
#    text(lakessizes, par("usr")[3] - 0.025, srt = 45, adj = 1,labels = signif(lakessizes,2),xpd = TRUE, font = 2)
    dev.off();

    dat<-read.table(inputfile2,header=TRUE);
    figure3<-paste("Speciation_Lakes_",sufix,".png",sep=""); 
    png(figure3,width=1980,height=1280,res=300);
    lakesnames<-paste("lake",1:nsites,sep="");#we should replace this list by the real names of the lakes.
#    absolutelakessizes<-readSizeOfLakes(verticesdatafile);#sizes of the lakes
#    relativelakessizes<-(dat$Ji/dat$J);
#    sumlakessizes<-sum(absolutelakessizes);
#    lakessizes<-absolutelakessizes*sumlakessizes;
    lakessizes<-readSizeOfLakes(verticesdatafile);#sizes of the lakes
    orderbysize<-order(lakessizes);
    lakesnames<-lakesnames[orderbysize];
    lakessizes<-lakessizes[orderbysize];
    mp<-barplot(rbind(dat$SpecANA[orderbysize],dat$SpecCLA[orderbysize],dat$SpecMR[orderbysize],dat$DispersalRich[orderbysize]),main="Origin of species",xlab="Size of the lakes (Km^2)",ylab="Species",col=c("green","darkblue","yellow","red"),legend=c("ANA","CLA","Regional-Mig","Local-Mig"),args.legend=list(x = "topleft"));
    text(mp, par("usr")[3] - 0.025, srt = 45, adj = 1,labels = signif(lakessizes,2),xpd = TRUE, font = 2)
    dev.off();
}


preprocfunc<-function(inputfile){
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

    for(i in 1:nsites){
        for(j in 1:nvar){
            cat(paste(averageval[i,j]," ",sep=""),sep=" ",file=averageoutputfile,append=TRUE);
            cat(paste(standarddevval[i,j]," ",sep=""),sep=" ",file=sdoutputfile,append=TRUE);
    	}
    	cat("\n",file=averageoutputfile,append=TRUE);
    	cat("\n",file=sdoutputfile,append=TRUE);
    }

    return(list(file=averageoutputfile,nsites=nsites,sufix=sufix));
}

args <- commandArgs(trailingOnly = TRUE)
inputfile<-args[1];
verticesdatafile<-args[2];

preprocres<-preprocfunc(inputfile);
printfigures(inputfile,preprocres$file,preprocres$nsites,preprocres$sufix,verticesdatafile);
