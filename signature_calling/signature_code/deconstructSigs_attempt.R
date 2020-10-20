require("deconstructSigs")
finalFrame<-data.frame(matrix(ncol=6))
colnames(finalFrame)<-c("sample.id","chr","pos","ref","alt","Ancestral_State")
oldFrame<-read.table("~/../Desktop/Ongoing_Work/Cancer/LUAD_Phy_Informed/MAFF/LookedUp_2_met_trees_as_MAFish.csv",sep=" ",header = T)

#Lets normal first, no truncal
nFrame<-oldFrame
nFrame<-nFrame[,c("Tumor_allele","Reference_Allele","pos","Unique_patient_identifier","chr")]
colnames(nFrame)<-c("alt","ref","pos","sample.id","chr")
sigs.input<-mut.to.sigs.input(nFrame,sample.id = "sample.id",chr = "chr",pos="pos",ref="ref",alt="alt")
res435<-whichSignatures(tumor.ref = sigs.input,signatures.ref = signatures.exome.cosmic.v3.may2019,sample.id = "p435",contexts.needed = T,tri.counts.method = "exome")
res459<-whichSignatures(tumor.ref = sigs.input,signatures.ref = signatures.exome.cosmic.v3.may2019,sample.id = "p459",contexts.needed = T,tri.counts.method = "exome")
pdf(file= "p435_whole_non_truncal_signatures.pdf")
plotSignatures(res435)
dev.off()
pdf(file = "p435_whole_non_truncal_signatures_pie.pdf")
makePie(res435)
dev.off()
pdf(file = "p459_whole_non_truncal_signatures.pdf")
plotSignatures(res459)
dev.off()
pdf(file = "p459_whole_non_truncal_signatures_pie.pdf")
makePie(res459)
dev.off()
save(res435,file = "whole_435_non_truncal_deconstructSigs.output.RData")
save(res459,file="whole_459_non_truncal_deconstructSigs.output.RData")

#lets do ero verus nonero each
nFrame<-oldFrame
nFrame<-nFrame[,c("Tumor_allele","Reference_Allele","pos","Unique_patient_identifier","chr","Ero_or_Not")]
colnames(nFrame)<-c("alt","ref","pos","sample.id","chr","Ero_or_Not")
nFrame[,"sample.id"]<-paste0(nFrame[,"sample.id"],nFrame[,"Ero_or_Not"])
sigs.input<-mut.to.sigs.input(nFrame,sample.id = "sample.id",chr = "chr",pos="pos",ref="ref",alt="alt")
for(i in rownames(sigs.input)){
  res<-whichSignatures(tumor.ref = sigs.input,signatures.ref = signatures.exome.cosmic.v3.may2019,sample.id = i,contexts.needed = T,tri.counts.method = "exome")
  pdf(file= paste0(i,"_non_truncal_signatures.pdf"))
  plotSignatures(res)
  dev.off()
  pdf(file = paste0(i,"_non_truncal_signatures_pie.pdf"))
  makePie(res)
  dev.off()
  save(res,file=paste0(i,"_non_truncal_deconstructSigs.output.RData"))
}

###Do it by Ero and Not only
nFrame<-oldFrame
nFrame<-nFrame[,c("Tumor_allele","Reference_Allele","pos","chr","Ero_or_Not")]
colnames(nFrame)<-c("alt","ref","pos","chr","sample.id")
#nFrame[,"sample.id"]<-paste0(nFrame[,"sample.id"],nFrame[,"Ero_or_Not"])
sigs.input<-mut.to.sigs.input(nFrame,sample.id = "sample.id",chr = "chr",pos="pos",ref="ref",alt="alt")
for(i in rownames(sigs.input)){
  res<-whichSignatures(tumor.ref = sigs.input,signatures.ref = signatures.exome.cosmic.v3.may2019,sample.id = i,contexts.needed = T,tri.counts.method = "exome")
  pdf(file= paste0(i,"_non_truncal_signatures.pdf"))
  plotSignatures(res)
  dev.off()
  pdf(file = paste0(i,"_non_truncal_signatures_pie.pdf"))
  makePie(res)
  dev.off()
  save(res,file=paste0(i,"_non_truncal_deconstructSigs.output.RData"))
}
