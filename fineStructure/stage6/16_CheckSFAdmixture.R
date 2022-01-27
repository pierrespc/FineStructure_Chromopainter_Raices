#!/bin/Rscript


require(stringr)
require(scales)
c25 <- c(
  "Caucasus"="#FF7F00",
  "NorthernAfrica"="#E31A1C",
  "SoutheasternMediterraneanSea"="#FDBF6F",
  "Levant"="maroon",
  
  "Iberia"="dodgerblue2",
  "Italia"="skyblue2",
  "NorthernEurope"="blue1",
  "BritishIslands"="steelblue4",
  "Basque"="blue3",
  "CentralEurope"="darkblue",
  
  
  
  "SouthernAfricanHunterGatherer"="mediumseagreen",
  "SouthernBantu"="green1",
  "RainforestHunterGatherer"="darkturquoise",
  "EasternAfrica"="darkgreen",
  "GuineanGulf"="palegreen2",
  "WesternAfrica"="green4",
  
  "NorthernNativeAmerican"="orchid3",
  "SouthAmericanTropicalForests"="deeppink1",
  "NativeAmericanCivilizations"="darkviolet"
)


putSpace<-function(string){
  for(let in LETTERS){
    string=str_replace_all(string,let,paste("\n",let,sep=""))
  }
  string=strsplit(string,split="")[[1]]
  return(paste(string[-1],collapse = ""))
}  

setwd(paste("~/Documents/OtherCollaborations/RAICES/outputsFineStructure/stage6/",sep=""))

SF<-read.table("ProportionsPerIndividual.tsv",header = T,stringsAsFactors =  F,sep="\t")
Admixture<-read.table("../../Admixture/BestRUNperK/Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP_REDUCED.Filtered.pruned.8.AncestryComponentByIndividual.txt",stringsAsFactors = F,header=T,sep="\t")
Admixture$Ind<-paste(Admixture$Population,Admixture$Ind,sep="___")
if(sum(grepl("___",SF$Target))!= sum(SF$Target %in% Admixture$Ind)){
  stop("pb SF and Admixture Ind columns")
}

SF<-SF[grepl("___",SF$Target),]
Admixture<-Admixture[,! names(Admixture) %in% c("Region","latitude","longitude","set","VCFid","cex","Point","Color")]
Admixture<-Admixture[,c("Ind","Population",
                        "brown",
                        "goldenrod","darkorange",
                        "darkolivegreen","palegreen","seagreen","limegreen",
                        "cadetblue")]




pdf("CorrelationsAdmixture_RawSourceFinder2.pdf")
#######comparing all CP/FS/SF clusters to all K in admixture...
colBlue=rgb(red=rep(0,1000),green=rep(0,1000),blue=rep(1,1000),alpha=seq(0,1,length=1000))
colRed=rgb(red=rep(1,1000),green=rep(0,1000),blue=rep(0,1000),alpha=seq(0,1,length=1000))

#for(pop in c(SF$Target[! grepl("___",SF$Target)],"All admixed")){
for(pop in c("PuertoMadryn","All admixed")[1]){
  out<-data.frame(matrix(NA,ncol(Admixture)-2,ncol(SF)-1))
  names(out)<-names(SF)[-1]
  row.names(out)<-names(Admixture)[-c(1,2)]
  plot(0,0,"n",xlab="Admixture components",ylab="SourceFinder Cluster",main=pop,
       xlim=c(-1.5,length(Admixture)-1),
       ylim=c(0,length(SF)+2),
       axes=F)
  y=0
  for(cluster in names(SF)[-1]){
    y=y+1
    x=0
    rect(xleft = x-1.5,xright=x-0.5,ybottom=y-0.5,ytop=y+0.5,col=c25[cluster])
    for(K in names(Admixture)[-c(1,2)]){
      x=x+1
      if(pop=="All admixed"){
        tmp<-merge(SF[,c("Target",cluster)],Admixture[,c("Ind",K)],by.x="Target",by.y="Ind")
      }else{
        tmp<-merge(SF[grepl(paste(pop,"___",sep=""),SF$Target),c("Target",cluster)],Admixture[Admixture$Population==pop,c("Ind",K)],by.x="Target",by.y="Ind")
      }
      cor=cor.test(tmp[,2],tmp[,3])
      corEst<-cor$estimate
      corP<-cor$p.value
      out[K,cluster]<-paste(corEst," (",scientific(corP,digits=3),")",sep="")
      corCol=ifelse(corEst>0,
                    colBlue[round(corEst*1000,digits=0)],
                    colRed[abs(round(corEst*1000,digits=0))])
      rect(xleft = x-0.5,xright=x+0.5,ybottom=y-0.5,ytop=y+0.5,col=corCol)
      text(x=x,y=y,labels = ifelse(corP<1e-4,"4",
                                   ifelse(corP<1e-3,"3",
                                          ifelse(corP<1e-2,"2",
                                                 ifelse(corP<1e-1,"1",""))))
      )
      if(cluster==names(SF)[2]){
        rect(xleft = x-0.5,xright=x+0.5,ybottom=length(SF)+0.5,ytop=length(SF)+1.5,col=K)
      }
    }
  }
}

dev.off()

###Combining Components to better see the correlations between Admixture and SF
AdmixtureGROUP<-read.table("../../Admixture/BestRUNperK/Genotipos_Raices.Plink.Autosomal.HGDP_1KG_SGDP_REDUCED.Filtered.pruned.8.MeanByGroup.txt",stringsAsFactors = F,header=T,sep="\t")

AdmixtureComb<-Admixture[,c("Ind","Population")]

Europe1<-which(AdmixtureGROUP["French",]==max(AdmixtureGROUP["French",]))
AdmixtureComb$Europe1<-Admixture[,names(AdmixtureGROUP)[Europe1]]

Europe2<-which(AdmixtureGROUP["Finnish",]==max(AdmixtureGROUP["Finnish",]))
AdmixtureComb$Europe2<-Admixture[,names(AdmixtureGROUP)[Europe2]]

MiddleEastNorthAfricaCaucasus<-which(AdmixtureGROUP["Bedouin",]==max(AdmixtureGROUP["Bedouin",]))
AdmixtureComb$MiddleEastNorthAfricaCaucasus<-Admixture[,names(AdmixtureGROUP)[MiddleEastNorthAfricaCaucasus]]

NativeAmerican<-which(AdmixtureGROUP["Karitiana",]==max(AdmixtureGROUP["Karitiana",]))
AdmixtureComb$NativeAmerican<-Admixture[,names(AdmixtureGROUP)[NativeAmerican]]

SubsaharanAfrica<-which(AdmixtureGROUP["Subsaharan African",]>0.05)
AdmixtureComb$SubsaharanAfrica<-apply(Admixture[,names(AdmixtureGROUP)[SubsaharanAfrica]],1,sum)

if(length(unique(c(Europe1,Europe2,MiddleEastNorthAfricaCaucasus,NativeAmerican,SubsaharanAfrica)))!=length(AdmixtureGROUP)){
  stop("Admixture combinations not done well")
}

SFcomb<-data.frame(matrix(NA,nrow(SF),length(AdmixtureComb)-1))
names(SFcomb)<-c("Target",names(AdmixtureComb)[-c(1,2)])
SFcomb$Target<-SF$Target

Europe1<-c("Iberia","Italia","Basque")
SFcomb$Europe1<-apply(SF[,Europe1],1,sum)

Europe2<-c("NorthernEurope","BritishIslands","CentralEurope")
SFcomb$Europe2<-apply(SF[,Europe2],1,sum)

MiddleEastNorthAfricaCaucasus<-c("Caucasus","NorthernAfrica","SoutheasternMediterraneanSea","Levant")
SFcomb$MiddleEastNorthAfricaCaucasus<-apply(SF[,MiddleEastNorthAfricaCaucasus],1,sum)

NativeAmerican<-c("NativeAmericanCivilizations","SouthAmericanTropicalForests","NorthernNativeAmerican")
SFcomb$NativeAmerican<-apply(SF[,NativeAmerican],1,sum)

SubsaharanAfrica<-c("SouthernAfricanHunterGatherer","SouthernBantu","RainforestHunterGatherer","EasternAfrica","GuineanGulf","WesternAfrica")
SFcomb$SubsaharanAfrica<-apply(SF[,SubsaharanAfrica],1,sum)

if(length(unique(c(Europe1,Europe2,MiddleEastNorthAfricaCaucasus,NativeAmerican,SubsaharanAfrica)))!=length(SF)-1){
  stop("SourceFinder combinations not done well")
}


c5<-c("Europe1"="goldenrod",
      "Europe2"="darkorange",
      "MiddleEastNorthAfricaCaucasus"="brown",
      "NativeAmerican"="cadetblue",
      "SubsaharanAfrica"="seagreen")

pdf("CorrelationsAdmixture_CombinatedSourceFinder.pdf")
for(pop in c("PuertoMadryn","All admixed")){
  out<-data.frame(matrix(NA,ncol(AdmixtureComb)-2,ncol(SFcomb)-1))
  names(out)<-names(SFcomb)[-1]
  row.names(out)<-names(AdmixtureComb)[-c(1,2)]
  plot(0,0,"n",xlab="Admixture components (Combined)",ylab="SourceFinder Clusters (Combined)",main=pop,
       xlim=c(-1.5,length(AdmixtureComb)-1),
       ylim=c(0,length(SFcomb)+2),
       axes=F)
  y=0
  for(cluster in names(SFcomb)[-1]){
    y=y+1
    x=0
    rect(xleft = x-1.5,xright=x-0.5,ybottom=y-0.5,ytop=y+0.5,col=c5[cluster])
    for(K in names(AdmixtureComb)[-c(1,2)]){
      x=x+1
      if(pop=="All admixed"){
        tmp<-merge(SFcomb[,c("Target",cluster)],AdmixtureComb[,c("Ind",K)],by.x="Target",by.y="Ind")
      }else{
        tmp<-merge(SFcomb[grepl(paste(pop,"___",sep=""),SFcomb$Target),c("Target",cluster)],AdmixtureComb[AdmixtureComb$Population==pop,c("Ind",K)],by.x="Target",by.y="Ind")
      }
      cor=cor.test(tmp[,2],tmp[,3])
      corEst<-cor$estimate
      corP<-cor$p.value
      out[K,cluster]<-paste(corEst," (",scientific(corP,digits=3),")",sep="")
      corCol=ifelse(corEst>0,
                    colBlue[round(corEst*1000,digits=0)],
                    colRed[abs(round(corEst*1000,digits=0))])
      rect(xleft = x-0.5,xright=x+0.5,ybottom=y-0.5,ytop=y+0.5,col=corCol)
      #text(x=x,y=y,labels = ifelse(corP<1e-4,"4",
      #                             ifelse(corP<1e-3,"3",
      #                                    ifelse(corP<1e-2,"2",
      #                                           ifelse(corP<1e-1,"1",""))))
      #)
      text(x=x,y=y,labels = scientific(corP,digits=3),cex=0.5)
      if(cluster==names(SFcomb)[2]){
        print("ho")
        rect(xleft = x-0.5,xright=x+0.5,ybottom=length(SFcomb)+0.5,ytop=length(SFcomb)+1.5,col=c5[K])
        
        
        
      }
    }
  }
  rect(xleft=rep(-1.5,2001),xright=rep(-0.5,2001),
       ybottom=seq(length(SFcomb)+0.5,length(SFcomb)+1.5,length=2001),
       ytop=seq(length(SFcomb)+0.5+1/2001,length(SFcomb)+1.5+1/2001,length=2001),
       col = c(colRed[c(1000:1)],colBlue),border=NA)
  
  text(x=rep(-0.4,9),y=seq(length(SFcomb)+0.5,length(SFcomb)+1.5,length=9),labels = seq(-1,1,0.25),adj = 0,cex=0.4)
  
}

dev.off()



#######combining  combined CP/FS/SF (leaving out those from populations with unclear components in Admixture) 
SFcomb2<-data.frame(matrix(NA,nrow(SF),length(AdmixtureComb)-1))
names(SFcomb2)<-c("Target",names(AdmixtureComb)[-c(1,2)])
SFcomb2$Target<-SF$Target

Europe1<-c("Iberia","Italia")
SFcomb2$Europe1<-apply(SF[,Europe1],1,sum)

Europe2<-c("NorthernEurope","BritishIslands")
SFcomb2$Europe2<-apply(SF[,Europe2],1,sum)

MiddleEastNorthAfricaCaucasus<-c("NorthernAfrica","SoutheasternMediterraneanSea","Levant")
SFcomb2$MiddleEastNorthAfricaCaucasus<-apply(SF[,MiddleEastNorthAfricaCaucasus],1,sum)

NativeAmerican<-c("NativeAmericanCivilizations","SouthAmericanTropicalForests","NorthernNativeAmerican")
SFcomb2$NativeAmerican<-apply(SF[,NativeAmerican],1,sum)

SubsaharanAfrica<-c("SouthernAfricanHunterGatherer","SouthernBantu","RainforestHunterGatherer","EasternAfrica","GuineanGulf","WesternAfrica")
SFcomb2$SubsaharanAfrica<-apply(SF[,SubsaharanAfrica],1,sum)

pdf("CorrelationsAdmixture_CombinatedSourceFinder2.pdf")
for(pop in c("PuertoMadryn","All admixed")){
  out<-data.frame(matrix(NA,ncol(AdmixtureComb)-2,ncol(SFcomb2)-1))
  names(out)<-names(SFcomb2)[-1]
  row.names(out)<-names(AdmixtureComb)[-c(1,2)]
  plot(0,0,"n",xlab="Admixture components (Combined)",ylab="SourceFinder Clusters (Combined clear clusters)",main=pop,
       xlim=c(-1.5,length(AdmixtureComb)-1),
       ylim=c(0,length(SFcomb2)+2),
       axes=F)
  y=0
  for(cluster in names(SFcomb2)[-1]){
    y=y+1
    x=0
    rect(xleft = x-1.5,xright=x-0.5,ybottom=y-0.5,ytop=y+0.5,col=c5[cluster])
    for(K in names(AdmixtureComb)[-c(1,2)]){
      x=x+1
      if(pop=="All admixed"){
        tmp<-merge(SFcomb2[,c("Target",cluster)],AdmixtureComb[,c("Ind",K)],by.x="Target",by.y="Ind")
      }else{
        tmp<-merge(SFcomb2[grepl(paste(pop,"___",sep=""),SFcomb2$Target),c("Target",cluster)],AdmixtureComb[AdmixtureComb$Population==pop,c("Ind",K)],by.x="Target",by.y="Ind")
      }
      cor=cor.test(tmp[,2],tmp[,3])
      corEst<-cor$estimate
      corP<-cor$p.value
      out[K,cluster]<-paste(corEst," (",scientific(corP,digits=3),")",sep="")
      corCol=ifelse(corEst>0,
                    colBlue[round(corEst*1000,digits=0)],
                    colRed[abs(round(corEst*1000,digits=0))])
      rect(xleft = x-0.5,xright=x+0.5,ybottom=y-0.5,ytop=y+0.5,col=corCol)
      #text(x=x,y=y,labels = ifelse(corP<1e-4,"4",
      #                             ifelse(corP<1e-3,"3",
      #                                    ifelse(corP<1e-2,"2",
      #                                           ifelse(corP<1e-1,"1",""))))
      #)
      text(x=x,y=y,labels = scientific(corP,digits=3),cex=0.5)
      if(cluster==names(SFcomb2)[2]){
        print("ho")
        rect(xleft = x-0.5,xright=x+0.5,ybottom=length(SFcomb2)+0.5,ytop=length(SFcomb2)+1.5,col=c5[K])
        
        
        
      }
    }
  }
  rect(xleft=rep(-1.5,2001),xright=rep(-0.5,2001),
       ybottom=seq(length(SFcomb2)+0.5,length(SFcomb2)+1.5,length=2001),
       ytop=seq(length(SFcomb2)+0.5+1/2001,length(SFcomb2)+1.5+1/2001,length=2001),
       col = c(colRed[c(1000:1)],colBlue),border=NA)
  
  text(x=rep(-0.4,9),y=seq(length(SFcomb2)+0.5,length(SFcomb2)+1.5,length=9),labels = seq(-1,1,0.25),adj = 0,cex=0.4)
  
}

dev.off()