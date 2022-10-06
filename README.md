# Promoter_similarities

## Preliminary datasets
<h5><a href="scripts/EPD_promoter_analysis.R">EPD_promoter_analysis.R</a></h5>
<h5><a href="scripts/FANTOM3_4_datasets.R">FANTOM3_4_datasets.R</a></h5>


## Correlation of activity between proximal promoters

<h5>1 - <a href="scripts/FANTOM5_prom_ID.R">FANTOM5_prom_ID.R</a></h5>
Identify human "True TSS clusters" from FANTOM5 dataset (Data/TSS_human.bed) as the onea with a minimum score value and associated with an ENSEMBL element (Data/ensembl_genes_grch37.txt: downloaded from BioMart) <br>
Outputs are in Tables/h37_tss_selected_gr.bed and in Tables/clusters_selected_genes.csv <br>

<h5>2 - <a href="scripts/fantom5_links.R">fantom5_links.R</a></h5>
Get the links where the FANATOM5 CAGE files are, using CAGEr package. <br>
The links are in Tables/FANTOM5_files.txt and the whole sample information table of the in Tables/FANTOM5_files_info.tsv 

<h5>3 - <a href="scripts/FANTOM5_files.sh">FANTOM5_files.sh</a></h5>
Bash script to download the files from the links in Tables/FANTOM5_files.txt in parallel, saving the in folder Data/FANTOM5_files/

<h5>4 - <a href="scripts/init_CAGEfiles.R">init_CAGEfiles.R</a></h5>
Loading the CAGE files (in parallel), merging data from the same cell_type and overlapping with the selected TrueTSS regions (Tables/h37_tss_selected_gr.bed). <br>
Output to Tables/true_tss_activ_tissue.rds, with a dataframe for each cell type/tissue in each element from the list. 

<h5>5 - <a href="scripts/tissue_activity.R">tissue_activity.R</a></h5>
With Tables/true_tss_activ_tissue.rds as input, firstly it restricts to the promoters in Tables/clusters_selected_genes.csv. <br>
Secondly, summarizes data from all samples, to filter the active promoters along multiple samples (Tables/active_promoters.tsv) and those highly tissue specific (Tables/tissue_restricted/prom.tsv). 

<h5>6 -<a href="scripts/distances.R">distances.R</a></h5>
Calculate distance to the nearest active promoter in each gene for each promoter in Tables/active_promoters.tsv. <br>
Output in Tables/active_promoters_df.csv 


<h5>7 -<a href="scripts/closed_prom.R">closed_prom.R</a></h5>
Firstly, identify the IC95 of the width of the selected TSS clusters. <br>
Select the pairs of promoters from each gene that are within this genomic distance. <br>
proximal_promoters.rds contains the output.


<h5>8 -<a href="scripts/closed_prom_corr.R">closed_prom_corr.R</a></h5>
Variance of relative activity of promoters in proximal promoters (Tables/proximal_promoters.rds
) along the different human cell types (Tables/true_tss_activ_tissue.rds). <br>
Control is assessed by permuting the pairs of promoters
