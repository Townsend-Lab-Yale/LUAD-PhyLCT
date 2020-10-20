save(newOdd459,file="nearMaff_newOdd459.RData")
nFrame<-data.frame(matrix(ncol=5,nrow=0))
colnames(nFrame)<-c("alt","ref","pos","sample.id","chr")
seqsOnly<-newOdd459[,5:17]
for(i in 1:ncol(seqsOnly)){
  seqsOnly[,i]<-as.character(seqsOnly[,i])
  uniID<-paste0("459_",colnames(seqsOnly)[i])
  for(j in 1:nrow(seqsOnly)){
    alt<-seqsOnly[j,i]
    ref<-as.character(newOdd459[j,"Normal"])
    pos<-newOdd459[j,"Position"]
    chr<-as.character(newOdd459[j,"Chromosome"])
 
    nFrame[nrow(nFrame)+1,]<-c(alt,ref,pos,uniID,chr)
  }
}
nFrame[,"pos"]<-as.numeric(nFrame[,"pos"])
rowsToDrop<-c()
for(i in 1:nrow(nFrame)){
  if(nFrame[i,"ref"]==nFrame[i,"alt"]){rowsToDrop<-c(rowsToDrop,i)}
}
nFrame<-nFrame[-rowsToDrop,]
  


sigs.input<-mut.to.sigs.input(nFrame,sample.id = "sample.id",chr = "chr",pos="pos",ref="ref",alt="alt")
for(i in rownames(sigs.input)){
  res<-whichSignatures(tumor.ref = sigs.input,signatures.ref = signatures.exome.cosmic.v3.may2019,sample.id = i,contexts.needed = T,tri.counts.method = "exome")
  pdf(file= paste0(i,"_anc_signatures.pdf"))
  plotSignatures(res)
  dev.off()
  pdf(file = paste0(i,"_anc_signatures_pie.pdf"))
  makePie(res)
  dev.off()
  save(res,file=paste0(i,"_anc_deconstructSigs.output.RData"))
}