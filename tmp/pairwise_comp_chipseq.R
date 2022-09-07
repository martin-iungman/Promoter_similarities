#Script para comparar muestras de ChIPseq

library(compGenomRData)
#datasets disponibles
data_path = system.file('extdata/chip-seq',package='compGenomRData')
chip_files = list.files(data_path, full.names=TRUE)

# load the chromosome info package
library(GenomeInfoDb)

# fetch the chromosome lengths for the human genome
hg_chrs = getChromInfoFromUCSC('hg38')

# find the length of chromosome 21 (en el que se centra el ejemplo)
hg_chrs = subset(hg_chrs, grepl('chr21$',chrom))

#The tileGenome() function from the GenomicRanges package constructs equally sized windows over the genome of interest. The function takes two arguments:
#Firstly, we convert the chromosome lengths data.frame into a named vector.

# downloaded hg_chrs is a data.frame object,
# we need to convert the data.frame into a named vector
seqlengths = with(hg_chrs, setNames(size, chrom))

#Then we construct the windows.
# load the genomic ranges package
library(GenomicRanges)

# tileGenome function returns a list of GRanges of a given width, 
# spanning the whole chromosome
tilling_window = tileGenome(seqlengths, tilewidth=1000)

# unlist converts the list to one GRanges object
tilling_window = unlist(tilling_window)

#We will use the summarizeOverlaps() function from the GenomicAlignments package to count the number of reads in each genomic window.
#The function will do the counting automatically for all our experiments. 
#The summarizeOverlaps() function returns a SummarizedExperiment object.
#The object contains the counts, genomic ranges which were used for the quantification, and the sample descriptions.

# load GenomicAlignments
library(GenomicAlignments)

# fetch bam files from the data folder
bam_files = list.files(
  path       = data_path, 
  full.names = TRUE, 
  pattern    = 'bam$'
)

# use summarizeOverlaps to count the reads
so = summarizeOverlaps(tilling_window, bam_files)

# extract the counts from the SummarizedExperiment
counts = assays(so)[[1]]

#Normalizacion por CPM
# calculate the cpm from the counts matrix
# the following command works because 
# R calculates everything by columns
cpm = t(t(counts)*(1000000/colSums(counts)))

#We remove all tiles which do not have overlapping reads. 
#Tiles with 0 counts do not provide any additional discriminatory power, rather, they introduce artificial similarity between the samples (i.e. samples with only a handful of bound regions will have a lot of tiles with 0 counts, while they do not have to have any overlapping enriched tiles).

# remove all tiles which do not contain reads
cpm = cpm[rowSums(cpm) > 0,]

# change the formatting of the column names
# remove the .chr21.bam suffix
colnames(cpm) = sub('.chr21.bam','',   colnames(cpm))


# remove the GM12878_hg38 prefix
colnames(cpm) = sub('GM12878_hg38_','',colnames(cpm))

# calculates the pearson correlation coefficient between the samples
correlation_matrix = cor(cpm, method='pearson')


## Visualizacion
# load ComplexHeatmap
library(ComplexHeatmap)

# load the circlize package, and define 
# the color palette which will be used in the heatmap
library(circlize)
heatmap_col = circlize::colorRamp2(
  breaks = c(-1,0,1),
  colors = c('blue','white','red')
)

# plot the heatmap using the Heatmap function
Heatmap(
  matrix = correlation_matrix, 
  col    = heatmap_col
)