# HELP for karyotype function
#You can choose specific chr Ex:ChromosomsNumbers=c(10,5,4,20) or selecting 
#a subset of chr Ex:ChromosomsNumbers=c(1:5)
#chrX & chrY consider as 23 and 24 respectively
# GenesName is a list of genes name will be texting in the plot. it can be empty
# FocusBoundry accept an zooming region. Ex:FocusBoundry=c(9100000,99300000)
#if it left empty hole chr will be ploted.
#plotDir indicate the plot direction. 0 for horizontal and 1 for vertical
# locationBar determine whether the locations shown in the plot or not.(0 and 1)

#setwd("D:/Tehran Uni/Term 1/Mathematical Foundation Of Bioinformatics/HW/HW1/programming")
# install.packages('gridExtra')
# installed.packages('grid')

karyotype=function(ChromosomsNumbers,GenesName,plotDir,FocusBoundry,locationBar)
{

chrdig=400
  cytoBand <- read.delim("cytoBand.txt", header=TRUE)
  genesData=read.csv("Homo_sapiens.GRCh38.p10_genes.csv", header=TRUE,"\t")
  ChromosomsNumbers=sort(ChromosomsNumbers)
  numberOfChrs=length(ChromosomsNumbers)
  scale=10000
  if(length(FocusBoundry)>0 )
  {
    if(length(FocusBoundry)==2)
    {
      fo=sort(FocusBoundry)
      min1=fo[1]/scale
      max1=fo[2]/scale
    }else 
      break
    
  }else
  {
    min1=min(cytoBand[,2])
    max1=(max(cytoBand[,3]))/scale
  }
  colors=colorRampPalette(c("white", "black"))
  palete=c(colors(9),"#FF0000","#2473B4")
  
  partition=(min1+max1)/(numberOfChrs+4)
  plot(min1:max1,min1:max1,type="n",axes=F,ann=F)
  
  Y=0
  for(i in 1:numberOfChrs)
  {
    Y[i]=partition*(i-.5)
    chr=ChromosomsNumbers[i]
    if(chr==23)
      str="chrX"
    else if(chr==24)
    str="chrY"
    else
      str=paste("chr",as.character(chr),sep="")
    x=which(cytoBand[,1]==str)
    site=(cytoBand[x,2:3])/scale
    att=cytoBand[x,5]
    chrNo=dim(site)
    for(j in 1:chrNo[1])
    {
      colorPart=switch(att[j],
             "gneg"=palete[1],
             "gvar"=palete[2],
             "gpos25"=palete[3],
             "gpos33"=palete[4],
             "gpos50"=palete[5],
             "gpos66"=palete[6],
             "gpos75"=palete[7],
             "gpos"=palete[8],
             "gpos100"=palete[9],
             "acen"=palete[10],
             "stalk"=palete[11])
      VectorNumbersPartition=10
      VectorNumbers=(max1-min1)/VectorNumbersPartition  #20 numbers will be shown in plot for chr
      if(site[j,1]<min1 && site[j,2]<=max1 )
      {
        a=min1
        b=site[j,2]
      }else if(site[j,1]>=min1 && site[j,2]>max1)
      {
        a=site[j,1]
        b=max1
      }else if(site[j,1]>=min1 && site[j,2]<=max1)
      {
        a=site[j,1]
        b=site[j,2]
      }
      else
        next
      xD=c(site[j,1],site[j,2],site[j,2],site[j,1])
      yD=c(Y[i],Y[i],Y[i]+chrdig,Y[i]+chrdig)
      if(plotDir==1)
      {
        temp=xD
        xD=yD
        yD=temp
      }
      polygon(x=xD,y=yD,col=colorPart)
    }
    if(locationBar==1)
    {
    if(plotDir==0)
      Maxbound=xD
    else
      Maxbound=yD
    Y_VectorNumbers=Y[i]+chrdig*(-1/4)
    X_VectorNumbers=min1
    for(ii in 1:VectorNumbersPartition)
    {
      if(plotDir==0)
        text(x=X_VectorNumbers,y=Y_VectorNumbers,floor(X_VectorNumbers/100),adj = 1,cex=.6,font=4)
      else
        text(x=Y_VectorNumbers,y=X_VectorNumbers,floor(X_VectorNumbers/100),adj = 1,cex=.6,font=4)
      X_VectorNumbers=floor(VectorNumbers*ii+min1)
     
      if(X_VectorNumbers>Maxbound)
      {
        break
      }
    }
    }
    if(plotDir==1)
    {
      a=Y[i]+500
      b=-500
    }
    else
    {
      a=-300
      b=Y[i]+200
    }
    text(x=a,y=b,str,adj = 1,cex=.8)
  }
# Find Genes Lcations
  colNum=1
  if(plotDir==1)
    colNum=3
  legend("topright", legend=c("gneg","gvar","gpos25", "gpos33","gpos50","gpos66","gpos75","gpos","gpos100","acen","stalk"),y.intersp=0.5,x.intersp =.5,fill=palete,text.font=4,box.lty=0,cex=0.55,ncol=1)
  GeneNum=length(GenesName)
  tmp=length(Y)
  if(tmp>1)
    bound=(Y[2]-Y[1])/1.5
  if(bound>1263)
    bound=1263
  for(i in 1:GeneNum)
  {
    G=which(genesData[,10]==GenesName[i])
    chr=genesData[G,1:3]
    f=strsplit(chr[1,1],"chr")[[1]]
    if(f[2]=="X")
    {
      f[2]="23"
    }else if(f[2]=="Y")
    {
      f[2]="24"
    }
    a=which(ChromosomsNumbers==as.integer(f[2]))
    if(exists("a") && length(a)>0)
    {
      
      locations=dim(chr)
      for(k in 1:locations[1])
      {
        start=chr[k,2]
        end=chr[k,3]
        start=(start+end/2)/scale
        #s2=Y[a]
        #s1=start[1,1]
        if(plotDir==1)
        {
          if(length(ChromosomsNumbers)>18)
          {
            border=250
            cex1=.37
          }
          else
          {
            border=450
            cex1=.55
          }
          adju=1
          y0=start
          y1=y0
          x0=Y[a]+400
          x1=Y[a]+bound-300
        }
        else
        {
          adju=0
          x0=start
          x1=x0
          y0=Y[a]+400
          y1=Y[a]+bound*1.5-300
          cex1=.55
          
        }
        text(x=x1, y=y1,label=GenesName[i],adj = 0,cex=cex1,col="#0000FF")
        segments(x0,y0,x1,y1,col = "#0000FF")
      }
    }
    
  }
}

#-----------------Main--------------------
ChromosomsNumbers=c(1:24)
GenesName=c("SNORD78","BET1L","MIR8078","KRT18P53")
FocusBoundry=c()
plotDir=0
locationBar=0
karyotype(ChromosomsNumbers,GenesName,plotDir,FocusBoundry,locationBar)
