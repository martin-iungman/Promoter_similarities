#GR de ejemplo
gr<-GRanges(c("chr2", "chr2"), IRanges(c(100000,500000),c(100300,500500)))
dna<-getSeq(Hsapiens,gr)
names(dna)<-c("seq1","seq2")
#determino la freq de los dicnucleotido
nuc = oligonucleotideFrequency(dna, width = 2)
nuc = as.data.frame(nuc)

#calculo porcentajes y redondeo. Finalmente me quedo con 
nuc = round(nuc/100,3)
gc_content<-nuc['GC']
