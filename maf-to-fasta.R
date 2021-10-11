#Defines function in R to turn all MAF files in a directory into their own multi-aligned fasta file, which can be used in BEAST and other programs

#mafs to SNV fasta
mafDirToFasta<-function(direct){
  #get starting Dir
  startDir<-getwd()
  setwd(direct) #change to that directory
  #get all Maf files
  allMafs<-c(list.files(pattern="*.maf"),list.files(pattern="*.MAF"), list.files(pattern="*.Maf"))
  allMafs<-unique(allMafs) #remove redundant matches
  mergeMaf<-NULL
  #combine relevant MAFs into one data frame
  for(i in 1:length(allMafs)){
    thisMaf<-allMafs[i]
    thisMafName<-strsplit(thisMaf,split="")[[1]]
    thisMafName<-paste0(thisMafName[1:(length(thisMafName)-4)],collapse="")
    hold<-read.csv(file = thisMaf,header = T,sep = "\t",stringsAsFactors = F)
    hold<-hold[which(nchar(hold$Reference_Allele)==1),]
    hold<-hold[which(nchar(hold$Tumor_allele)==1),]
    if(is.null(mergeMaf)){mergeMaf<-hold}else{mergeMaf<-rbind(mergeMaf,hold)}
  }
  allMem<-unique(mergeMaf$Unique_patient_identifier) #get all the unique patient IDs
  #generate the name for the normal somatic sample
  normNameTemp<-strsplit(allMem, split = "")
  normName<-paste0(c(normNameTemp[[1]][normNameTemp[[1]]==normNameTemp[[2]]],"N"),collapse = "")
  print(normName)
  #get all unique sites in the merged maf
  allSites<-unique(mergeMaf[,c("chr","Start_Position")])
  seqTable<-data.frame(matrix(nrow = nrow(allSites),ncol=(length(allMem)+3)))
  colnames(seqTable)<-c("chr","pos",normName,allMem)
  #traverse all possible sites in all samples and figure out if that sample has a variant
  for(i in 1:nrow(allSites)){
    thisSite<-allSites[i,]
    allMatches<-mergeMaf[mergeMaf$Start_Position==thisSite$Start_Position & mergeMaf$chr==thisSite$chr,]
    seqTable[i,"chr"]<-thisSite$chr
    seqTable[i,"pos"]<-thisSite$Start_Position
    allNorms<-c()
    for(j in 1:nrow(allMatches)){
      thisMatch<-allMatches[j,]
      allNorms<-unique(c(allNorms,thisMatch$Reference_Allele))
      if(length((allNorms))>1){print(allNorms);stop("Something is wrong; need to split alleles")}
      seqTable[i, thisMatch$Unique_patient_identifier]<-thisMatch$Tumor_allele
    }
    seqTable[i,normName]<-allNorms
  }
  #account for missing data
  seqTableN<-seqTable
  seqTableN[is.na(seqTableN)]<-"N"
  rawMultiFasta<-c()
  #create the multifasta, without regard to the normal sample
  for(k in 3:(ncol(seqTableN))){
    print(paste0((colnames(seqTable)[k])," has this many missing sites relative to other 2: "))
    print(sum(is.na(seqTable[,k]))/nrow(seqTable))
    thisHeader<-paste0(">",colnames(seqTableN)[k])
    thisSeq<-paste0(seqTableN[,k],collapse = "")
    rawMultiFasta<-c(rawMultiFasta,thisHeader,thisSeq)
  }
  #look to the normal reference to add the normal character to any sample without a variant in the MAF
  lookUpMultiFasta<-c()
  for(k in 3:ncol(seqTableN)){
    toChange<-which(seqTableN[,k]=="N")
    seqTableN[toChange,k]<-seqTableN[toChange,normName]
    thisHeader<-paste0(">",colnames(seqTableN)[k])
    thisSeq<-paste0(seqTableN[,k],collapse = "")
    lookUpMultiFasta<-c(lookUpMultiFasta,thisHeader,thisSeq)
  }
  #export all 3 versions of the fasta 
  writeLines(text = rawMultiFasta,con = paste0(c(normNameTemp[[1]][normNameTemp[[1]]==normNameTemp[[2]]],".raw.fa"),collapse = ""))
  writeLines(text = rawMultiFasta,con = paste0(c(normNameTemp[[1]][normNameTemp[[1]]==normNameTemp[[2]]],".looked-up.fa"),collapse = ""))
  write.table(x=seqTable,file = paste0(c(normNameTemp[[1]][normNameTemp[[1]]==normNameTemp[[2]]],".faTSV"),collapse = ""),sep="\t",row.names = F)
  
  setwd(startDir)
}
#apply to all supplied directories
sapply(X = list.dirs(recursive=F),FUN = mafDirToFasta)
