library("ape")
library("seqinr")

args=commandArgs(trailingOnly = T)
fic=args[1]
outfolder=args[2]

#fic=c("refs_cidB_Tn_F1.fas")
#fic=c("refs_test.fa")

sort_refs=function(fic,outfolder){
  name_fic=strsplit(basename(fic),split="[.]")[[1]][1]
  refs=read.dna(fic,format="fasta")
  length_seq=dim(refs)[2]
  n_refs=dim(refs)[1]
  M=as.matrix(dist.dna(refs,model="raw"))*length_seq
  names_refs=colnames(M)
  kept=rep(T,n_refs) #vector saying if ref is kept at the end
  for (i in 1:n_refs){
    for (j in 1:i){
      if (i!=j && M[i,j]<4 && kept[j]!=F){
        kept[i]=F
      }
    }
  }
  #subset the selected sequences
  seq_tot=as.character(refs)
  seq_subset=seq_tot[kept,]
  write.fasta(seq_subset,names=names_refs[kept],file.out=paste(outfolder,name_fic,"_sort.fa",sep=""))
  if (length(which(kept==T))==1){
    write.fasta(seq_subset,names=names_refs[kept],file=paste(outfolder,name_fic,"_sort.fa",sep=""))
  }
  else{
    seq_DNAbin=as.DNAbin(seq_subset)
    write.FASTA(seq_DNAbin,file=paste(outfolder,name_fic,"_sort.fa",sep=""))
    
  }
}

sort_refs(fic,outfolder)

