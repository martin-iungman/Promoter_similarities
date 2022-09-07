#BiocManager::install("JASPAR2020")
# BiocManager::install("TFBSTools")
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)
hg38<-BSgenome.Hsapiens.UCSC.hg38

pfm_tbp <- getMatrixByName(JASPAR2020, "TBP")

seqLogo(toICM(pfm_tbp))

pwm_tbp<-toPWM(pfm_tbp, pseudocounts=0.8)

#desde un string
subject <- DNAString("GAATTCTCTCTTGTTGTAGTCTCTTGACAAAATG")
subject2 <- DNAString("GAATTCTCTATATATACTTGTTGTAGTCTCTTGACAAAATG")
siteset <- searchSeq(pwm_tbp,subject, min.score="50%", strand="*")
head(writeGFF3(siteset))
writeGFF3(siteset)
relScore(siteset)

#desde granges
#primero hay que obtener el dna referido al granges, especificando la especie
#agregarle nombre a los elementos para identificar despues y seguir el mismo camino de antes
library(GenomicRanges)
gr<-GRanges(c("chr2", "chr2"), IRanges(c(100000,500000),c(100300,500500)))
dna<-getSeq(Hsapiens,gr)
names(dna)<-c("seq1","seq2")
siteset2 <- searchSeq(pwm_tbp,dna, min.score="80%", strand="*")
writeGFF3(siteset2)
relScore(siteset2)
