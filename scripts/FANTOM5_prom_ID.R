library(rtracklayer)
library(stringi)
library(tidyverse)

human_tss_gr <- import.bed("Data/TSS_human.bed") # importa un archivo bed y lo graba como GRanges
gene_score_2p <- tibble(
  name = mcols(human_tss_gr)$name,
  gene = stri_replace_last_regex(str = mcols(human_tss_gr)$name, pattern = ",\\d.\\d\\d\\d\\d", replacement = ""),
  score = as.numeric(stri_extract_last_regex(str = mcols(human_tss_gr)$name, pattern = "\\d.\\d\\d\\d\\d"))
)
#mcols(human_tss_gr)<-cbind(mcols(human_tss_gr), gene=gene_score_2p$gene)

gene_score_2p <- filter(gene_score_2p, score >= 0.228)
gene_names <- str_split(gene_score_2p$gene, ",") 
gene_names_separated <- map(gene_names, str_split, "@")
head(gene_names_separated, 25)
all_genes <- unlist(gene_names_separated) %>%
  matrix(ncol = 2, byrow = T)
all_genes <- unique(all_genes[,2]) %>% tibble::enframe()
names(all_genes) <- c("id","name")
head(all_genes, 40)
nrow(all_genes)

ENS_txs <- read_tsv(file="Data/ensembl_genes_grch37.txt",col_types = "cccf")
names(ENS_txs) <- c("gene_id","tx_id","gene_symbol", "chr")
head(ENS_txs)
nrow(ENS_txs)

symbol_idx <- which(all_genes$name %in% ENS_txs$gene_symbol)

tx_idx <- which(all_genes$name %in% ENS_txs$tx_id)

all_genes[tx_idx,"gene_symbol"] <- ENS_txs %>%
  select(tx_id, gene_symbol) %>%
  left_join(all_genes[tx_idx,], ., by = c("name" = "tx_id")) %>%
  select(gene_symbol)

all_genes[symbol_idx,"gene_symbol"] <- all_genes[symbol_idx,"name"]

all_genes <- ENS_txs %>%
  select(gene_symbol, gene_id) %>%
  distinct() %>%
  left_join(all_genes, .)

gene_names_matrix <- str_split(gene_score_2p$gene, ",", simplify = T) 
gene_names_matrix <- gsub(".*@","",gene_names_matrix) %>% as_tibble()
colnames(gene_names_matrix) <- paste0("gene_",1:ncol(gene_names_matrix))
gene_names_matrix <- select(gene_score_2p,name) %>% bind_cols(gene_names_matrix)
gene_names_long <- gene_names_matrix %>%
  rename(cluster_id="name") %>%
  pivot_longer(-cluster_id, values_to = "genes") %>%
  select(-name)
gene_names_long

gene_names_long <- gene_names_long %>%
  filter(!genes=="")

gene_names_long <- gene_names_long %>%
  left_join(all_genes, by = c("genes" = "name"))
head(gene_names_long)

clusters_selected_genes <- gene_names_long %>%
  filter(!is.na(gene_id)) %>%
  select(-genes,-id) %>%
  distinct()
length(unique(clusters_selected_genes$cluster_id))
length(unique(clusters_selected_genes$gene_id))
length(unique(clusters_selected_genes$gene_symbol))
# paso 172502 cluster en 22435 genes (gene_id, en gene_symbol son 20935. hay una inconsistencia)

h37_tss_selected_gr <- subset(human_tss_gr, name %in% clusters_selected_genes$cluster_id)
clusters_selected_genes %>%
  group_by(gene_symbol, gene_id) %>%
  mutate(Prom = 1:n()) -> clusters_selected_genes

#clusters_selected_genes<-clusters_selected_genes%>%mutate(cluster_id=str_remove(cluster_id, ",.+"))

genes_alt_prom <- clusters_selected_genes %>%
  pivot_wider(names_from = Prom, values_from = cluster_id)

h37_tss_selected_gr<-export.bed(h37_tss_selected_gr,"Tables/h37_tss_selected_gr.bed")
clusters_selected_genes<-write_csv(clusters_selected_genes,"Tables/clusters_selected_genes.csv")
