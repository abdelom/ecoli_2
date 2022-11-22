library(HMM)

extrafigures=0

FragmentPosions<-data.frame(recombinant=c(0),rec=c(0),begin=c(0),end=c(0),nb1=c(0),nbpos1=c(0))
i<-1
pdf ("Projet.pdf") 
k_vec = c(1:6, 8, 8, 11:20)
j_vec = c(9:15, 17, 19:28)
for (v in 1:length(j_vec)) {
  j = j_vec[v]
  k=k_vec[v]
  print(k)
  galkpos=769500

  Donor<-seq(0,5000000,1000)*0
  Recipient<-seq(0,5000000,1000)*0
  Donor<-Donor[-1]
  Recipient<-Recipient[-1]
  TabSortedCompiled<-c()
  
  name2=paste("position/", "K12_sample",j, "relgal", k,"_pos.txt",sep="")
  name1=paste("position/","Rel606_sample",j, "relgal", k,"_pos.txt",sep="")
  print(name2)
  print(name1)
  Tab<-read.table(name1)
  Tab2<-read.table(name2)

  
  #exclude problematiq regions
  Tab2<-Tab2[Tab2[,2]<3947000 | Tab2[,2]>=3949000,]
  Tab2<-Tab2[Tab2[,2]<3481000| Tab2[,2]>=3482000,]
  Tab<-Tab[Tab[,1]<3947000 | Tab[,1]>=3949000,]
  Tab<-Tab[Tab[,1]<3481000| Tab[,1]>=3482000,]
  


  #bin data for figures
  a<-hist(Tab[Tab[,2]>=0,1],breaks=seq(0,5000000,1000),plot=F)#,main=paste(j,k))
  #abline(v=769500,lwd=3,col="yellow")
  b<-hist(Tab2[Tab2[,2]>=0,2],breaks=seq(0,5000000,1000),plot=F)#,main=paste(j,k))
  #abline(v=769500,lwd=3,col="yellow")
  # c<-hist(Tab2[,1],breaks=seq(0,5000000,1000),plot=F)#,main=paste(j,k))
  #abline(v=769500,lwd=3,col="yellow")
  
  
  Donor<-Donor+b$counts
  Recipient<-Recipient+a$counts
  
  
  #combine all reads according to the recipient position and attibute to who they matched, keeping only reads that matched both genomes but preferentially one
  TabSorted1<-as.data.frame(Tab[Tab[,2]>=0 & abs(Tab[,1]-Tab[,2])<500000,1])
  TabSorted1$source<-0
  #0 means reacipient match
  dimnames(TabSorted1)[[2]][1]<-"position"
  
  TabSorted2<-as.data.frame(Tab2[Tab2[,2]>=0 & abs(Tab2[,1]-Tab2[,2])<500000,2])
  TabSorted2$source<-1
  #1 means donor match
  dimnames(TabSorted2)[[2]][1]<-"position"
  TabSorted<-rbind(TabSorted1,TabSorted2)
  if (k==1){
    TabSortedCompiled<-TabSorted
    }
  else   
    {TabSortedCompiled<-rbind(TabSortedCompiled,TabSorted)}

  #combine the reads info of both table 0 is for reads that match the recipient, 1 for the ones that match the donor   
  TabSortedCompiled<-TabSortedCompiled[order(TabSortedCompiled$position),]
  # hmm=initHMM(c("MG","REL"),c(0,1),startProbs=c(0.1,0.9),transProbs=matrix(c(0.9999999,0.0000001,0.0000001,0.9999999),2),emissionProbs=matrix(c(0.2,0.8,0.8,0.2),2))
  # hmm=initHMM(c("MG","REL"),c(0,1),startProbs=c(0.1,0.9),transProbs=matrix(c(0.9999999,0.0000001,0.0000001,0.9999999),2),emissionProbs=matrix(c(0.5,0.95,0.5,0.05),2))
  
  # run an hmm to find the boudaries: 3 states Recipient, Donor , and both (some clones have not resolved)
  print("hmm")
  if (j<49) {       
    hmm=initHMM(c("MG","REL","Mix"),c(0,1),startProbs=c(0.1,0.9,0.1),transProbs=matrix(c(0.9999999,0.0000001,0.0000001,0.0000001,0.9999999,0.0000001,0.0000001,0.0000001,0.9999999),3),emissionProbs=matrix(c(0.1,0.9,0.5,0.9,0.1,0.5),3))
    }else
  {
    #  hmm=initHMM(c("MG","REL"),c(0,1),startProbs=c(0.1,0.9),transProbs=matrix(c(0.9999999,0.0000001,0.0000001,0.9999999),2),emissionProbs=matrix(c(0.1,0.9,0.9,0.1),2))
    hmm=initHMM(c("MG","REL","Mix"),c(0,1),startProbs=c(0.1,0.9,0.1),transProbs=matrix(c(0.9999999,0.000000001,0.000000001,0.000000001,0.999999999,0.000000001,0.000000001,0.000000001,0.999999999),3),emissionProbs=matrix(c(0.1,0.9,0.5,0.9,0.1,0.5),3))
  }
  
  
  #posterior=posterior(hmm,as.factor(TabSorted$source))
  viterbihmm=viterbi(hmm,as.factor(TabSortedCompiled$source))
  
  
  TabSortedCompiled$State<-viterbihmm
  TabSortedCompiled$Donor<-as.numeric(TabSortedCompiled$State!="REL")
  
  #plot figure with fit
  print("figure")
  TabSortedCompiled$cumsum<-cumsum(as.numeric(TabSortedCompiled$source*2-1))
  #plot(TabSortedCompiled$position/1e6,TabSortedCompiled$cumsum,pch='.',type='l',main=j,xlim=c(0.5,1),ylim=c(min(TabSortedCompiled[TabSortedCompiled$position/1e6>0.5 & TabSortedCompiled$position/1e6<1,c("cumsum")]),max(TabSortedCompiled[TabSortedCompiled$position/1e6>0.5 & TabSortedCompiled$position/1e6<1,c("cumsum")])))
  # maximim<-min(TabSortedCompiled$cumsum)
  #lines(TabSortedCompiled$position/1e6,TabSortedCompiled$Donor*maximim,col="red",lwd=1)
  # abline(v=769500/1e6,lwd=3,col="yellow")
  
  #plot figure with fit
  plot(a$mids/1e6,log10(Recipient),pch=20,cex=1,main=paste("Recombinant",j),xlab="genome position (Mb)",ylab="log10(coverage)",ylim=c(0,max(log10(Recipient),log10(Donor))),col="#155289",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  # points(b$mids,log10(Donor),pch=20,cex=1,col="blue")
  points(b$mids/1e6,log10(Donor),pch=20,cex=1,col= "#F4E6AA")
  maximim<-max(log10(Donor),log10(Recipient))
  lines(TabSortedCompiled$position/1e6,TabSortedCompiled$Donor*maximim,col="red",lwd=1)
  abline(v=galkpos/1e6,lwd=1,col="yellow")
  
  if(extrafigures==1)
  {
    #plot zoom  with fit
    plot(a$mids/1e6,log10(Recipient),pch=20,cex=1,main=paste("Recombinant",j),xlab="genome position (Mb)",ylab="log10(coverage)",ylim=c(0,2),col="#155289",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,xlim=c(0.7,1))
    # points(b$mids,log10(Donor),pch=20,cex=1,col="blue")
    points(b$mids/1e6,log10(Donor),pch=20,cex=1,col= "#F4E6AA")
    maximim<-2
    lines(TabSortedCompiled$position/1e6,TabSortedCompiled$Donor*maximim,col="red",lwd=1)
    abline(v=galkpos/1e6,lwd=1,col="yellow")
  }
  
  
  #find the length of segments
  print("find length")
  Tabtmp<-TabSortedCompiled
  index=which(Tabtmp$Donor==1)
  rec0=1
  
  while(length(index)>0)
  {
    pd=Tabtmp$position[index[1]]
    Tabtmp<-Tabtmp[c(index[1]:nrow(Tabtmp)),]
    index=which(Tabtmp$Donor==0)
    if (length(index)>0)
    {pe=Tabtmp$position[index[1]]}
    else {pe=4639675}
    FragmentPosions[i,"recombinant"]=j
    
    
    FragmentPosions[i,"rec"]=rec0
    FragmentPosions[i,"begin"]=pd
    FragmentPosions[i,"end"]=pe
    ones<-TabSortedCompiled[TabSortedCompiled$position>=pd & TabSortedCompiled$position<pe,]
    FragmentPosions[i,"nb1"]=sum(ones$source)
    pos<-unique(ones[ones$source==1,]$position)
    #FragmentPosions[i,"nbpos1"]=length(unique(ones[ones$source==1,]$position))
    FragmentPosions[i,"nbpos1"]=length(pos)
    FragmentPosions[i,"length"]=max(pos)-min(pos)
    print(FragmentPosions[i,])
    #plot zoom  with fit  
    if(extrafigures==1)
    {
      plot(a$mids/1e6,log10(Recipient),pch=20,cex=1,main=paste("Recombinant",j,rec0),xlab="genome position (Mb)",ylab="log10(coverage)",ylim=c(0,max(log10(Recipient),log10(Donor))),col="#155289",cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5,xlim=c(pd/1e6-0.05,pe/1e6+0.05))
      # points(b$mids,log10(Donor),pch=20,cex=1,col="blue")
      points(b$mids/1e6,log10(Donor),pch=20,cex=1,col= "#F4E6AA")
      maximim<-max(log10(Recipient),log10(Donor))
      lines(TabSortedCompiled$position/1e6,TabSortedCompiled$Donor*maximim,col="red",lwd=1)
    }
    
    
    i<-i+1
    rec0<-rec0+1
    
    if (length(index)==0){break}
    Tabtmp<-Tabtmp[c(index[1]:nrow(Tabtmp)),]
    index=which(Tabtmp$Donor==1)
  }
  print(j)
}

dev.off()

#FragmentPosions$length<-FragmentPosions$end-FragmentPosions$begin+1
hist(log10(FragmentPosions$length),breaks=20)
DistREL<-hist(log10(FragmentPosions[FragmentPosions$recombinant<=28,]$length),breaks=seq(0,7,0.25))
#DistHS<-hist(log10(FragmentPosions[FragmentPosions$recombinant>28 & FragmentPosions$recombinant>49,]$length),breaks=seq(0,7,0.25))
#Dist536<-hist(log10(FragmentPosions[FragmentPosions$recombinant>48,]$length),breaks=seq(0,7,0.25))

write.csv(FragmentPosions, "fragment_pos")
rm(list=ls(all.names = TRUE))

