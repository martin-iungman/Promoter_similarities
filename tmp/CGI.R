library(AnnotationHub)
ah<-AnnotationHub()
dm<-query(ah, c("Homo sapiens", "CpG island","hg19"))
ah_cgi<-dm[['AH5086']]
ah_cgi=keepStandardChromosomes(ah_cgi, pruning.mode = "coarse")
countOverlaps(gr,ah_cgi)
fo<-findOverlapPairs(gr, ah_cgi)
width(gr)
values(fo)=DataFrame(length_cgi=width(second(fo)))

library(annotatr)
cgi<-build_annotations("hg19", "hg19_cpgs")
cgi<-cgi[cgi$type=="hg19_cpg_islands"]

gr<-GRanges(seqnames=c("chr1","chr1","chr1"), IRanges(c(29000,28750,135000),c(30000,28850,400000)))
cgi_reg<-annotate_regions(gr, cgi)
cgi_reg=keepStandardChromosomes(ah_cgi, pruning.mode = "coarse")
