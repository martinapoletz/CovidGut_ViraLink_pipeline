# CovidGut_ViraLink_pipeline


## 2_process_a_priori_networks

### Downloading_omnipath_dorothea.R

Script to download directed and signed OmniPath protein-protein interactions and regulatory interactions from DoRothEA (within OmniPath, confidence levels A,B,C)

Input: outdir        
Output: Omnipath and dorothea networks with correct headers

```r
# Capture  messages and errors to a file.
zz <- file("virallink.out", open="a")
sink(zz, type="message", append = TRUE)
message("\nStarting Omnipath/Dorothea download script: Downloading_omnipath_dorothea.R\n")

# Installing packages
if (!requireNamespace("tidyverse", quietly = TRUE)) 
  install.packages("tidyverse", repos = "https://cloud.r-project.org")
if (!requireNamespace("devtools", quietly = TRUE)) 
  install.packages("devtools", repos = "https://cloud.r-project.org")
if(!'OmnipathR' %in% installed.packages()[,"Package"]){
  require(devtools)
  install_github('saezlab/OmnipathR')
}

# Loading the required packages
library(tidyverse)
library(OmnipathR)

# Define parameters
args <- commandArgs(trailingOnly = TRUE)

# Check length of command line parameters
if (length(args) != 1){
  stop("Wrong number of command line input parameters. Please check.")
}

outdir <- args[1]

# Create output directory if doesn't exist
path <- file.path(outdir, "2_process_a_priori_networks", "unprocessed_networks")
dir.create(path, showWarnings = FALSE, recursive=TRUE)

##### PPIs #####

# Get OmniPath PPI interaction network
ia_omnipath <- import_omnipath_interactions() %>% as_tibble()

# If you need higher coverage on PPI add these with `bind_rows` below
ia_ligrec <- import_ligrecextra_interactions() %>% as_tibble()
ia_pwextra <- import_pathwayextra_interactions() %>% as_tibble()
ia_kinaseextra <- import_kinaseextra_interactions() %>% as_tibble()

interactions <- as_tibble(
  bind_rows(
    ia_omnipath %>% mutate(type = 'ppi'),
    ia_pwextra %>% mutate(type = 'ppi'),
    ia_kinaseextra %>% mutate(type = 'ppi'),
    ia_ligrec %>% mutate(type = 'ppi')))

# For the network propagation interactome we need such an interaction table where we have
# the sign of the interactions and the direction of the interactions.
# First we take all the availble databases from Omnipath then we filter those which are ditrected and only inhibition OR stimulation have per one interaction.
# It is a question whether you should use consensus or is_inhibition is_directed columns.

interactions_filtered <- interactions[interactions$consensus_direction==1,]
interactions_filtered <- interactions_filtered[(interactions_filtered $consensus_inhibition +interactions_filtered$consensus_stimulation)==1,]

# Filter for the interactions which are involved in loops or which edges have multiple edges.

human_directed_interactome <- graph_from_data_frame(interactions_filtered, directed = TRUE, vertices = NULL)
human_directed_interactome <- simplify(human_directed_interactome, remove.multiple = FALSE, remove.loops = TRUE)

# The original igraph function deletes the edge attribute so below you can see a slightly changed one.

a = which_multiple(human_directed_interactome)
edge_indexes = E(human_directed_interactome)
human_directed_interactome  = subgraph.edges(human_directed_interactome, edge_indexes[a==FALSE])

# Keep only the giant component.
# First we have it from Csanadi updated for the vertex ids. 
# https://lists.gnu.org/archive/html/igraph-help/2009-08/msg00064.html

giant.component <- function(graph, ...) {
  cl <- clusters(graph, ...)
  induced_subgraph(graph, which(cl$membership == which.max(cl$csize)))
}

human_directed_interactome_giant_component  = giant.component(human_directed_interactome) 

# Write out the igraph both in ncol format

write_graph(human_directed_interactome_giant_component, paste0(path,"/directed_human_ppi_interactome.ncol"), format="ncol")
human_directed_interactome_filtered = as_data_frame(human_directed_interactome_giant_component, what = "edges")
write.table(human_directed_interactome_filtered, paste0(path, "/omnipath_signed_directed.txt"))

##### Regulatory interactions ####

# Download DOROTHEA for transcriptional regulation interactions
# The confidence level will be those TF-target interactions which are currated (a,b,c)

ia_transcriptional <- import_dorothea_interactions(confidence_level = c('A', 'B', 'C')) %>% as_tibble()

# After importing dataframe we go through again the same steps as for the PPIs
# Filtering only interactions where we know the direction and only inhibitory or excitatory.

ia_transcriptional <- ia_transcriptional[ia_transcriptional$consensus_direction==1,]
ia_transcriptional <- ia_transcriptional[(ia_transcriptional$consensus_stimulation+ia_transcriptional$consensus_inhibition)==1,]

# Make a graph and exclude self loops and multiple edges.
# We do not make a giant compnent for the graph, due to the TF-target interactome is not fully discovered and has multiple valuable components. 

ia_transcriptional2 <- graph_from_data_frame(ia_transcriptional, directed = TRUE, vertices = NULL)
ia_transcriptional2 <- simplify(ia_transcriptional2, remove.multiple = FALSE, remove.loops = TRUE)
a = which_multiple(ia_transcriptional2)
edge_indexes = E(ia_transcriptional2)
ia_transcriptional2  = subgraph.edges(ia_transcriptional2, edge_indexes[a==FALSE])

# Write results into files for further work

write_graph(ia_transcriptional2, paste0(path, "/directed_human_TF_targets.ncol"), format="ncol")
human_tf_data_filtered = as_data_frame(ia_transcriptional2, what = "edges")
write.table(human_tf_data_filtered, paste0(path,"/dorothea_abc_signed_directed.txt"))

## reset message sink and close the file connection
sink(type="message")
close(zz)

```

### filter_network_expressed_genes.R

Script to filter networks for expressed genes - runs for regulatory interactions and for ppis

Input: 
Network file, space delimited, gene names in columns "source_genesymbol" "target_genesymbol" or "to" "from"
Table of genes which are expressed, tab delimited, gene names in column "Gene" (ouput from DESeq2)
ID type of the differentially expressed genes - uniprot or gene symbols

Output: Networks in same format as input network (but tab seperated), filtered to include only interactions where source and target node are in the expressed list.


```r
# Capture  messages and errors to a file.
zz <- file("virallink.out", open="a")
sink(zz, type="message", append = TRUE)
message("\nStarting contextualised network reconstruction: filter_network_expressed_genes.R\n")

# Installing packages
if (!requireNamespace("tidyverse", quietly = TRUE)) 
  install.packages("tidyverse", repos = "https://cloud.r-project.org")

# Load packages
library(tidyverse)

# Define parameters
args <- commandArgs(trailingOnly = TRUE)

# Check length of command line parameters
if (length(args) != 5){
  stop("Wrong number of command line input parameters. Please check.")
}

# Output directory
outdir <- args[5]

# ID type of the  expressed genes - uniprot or gene symbols
id_type <- args[4] # symbol or uniprot - the ids in the expression data

# Create output dir if required
path <- file.path(outdir, "2_process_a_priori_networks")
dir.create(path, showWarnings = FALSE, recursive=TRUE)

if (id_type == "symbol"){
  source_col = "source_genesymbol"
  target_col = "target_genesymbol"
} else if (id_type == "uniprot") {
  source_col = "to"
  target_col = "from"
}

# Gene expression file - tab delimited
expressed <- read.csv(args[1], sep = "\t")
dorothea <- args[2]
omnipath <- args[3]

files <- c(dorothea, omnipath)

for (i in files){
  
  # Network file - space delimited 
  network <- read.csv(file.path(i), sep = " ")
  
  # Filter source and target nodes
  network_f <- network %>% filter((get(source_col) %in% expressed$Gene) & (get(target_col) %in% expressed$Gene))
  
  # Get network name for out filename
  file <- strsplit(i, "/")[[1]][4]
  name <- strsplit(file, "_")[[1]][1]
  
  # Save output
  write.table(network_f, file = file.path(path, paste0(name, "_contextualised_network.txt")), sep = "\t", quote = F, row.names = F)
}

# reset message sink and close the file connection
sink(type="message")
close(zz)

```

### get_regulator_deg_network.R

Script to get regulator - differentially expressed genes from the contextualised regulatory network

Input: 
Contextualised regulatory network file (all nodes expressed, tab delimited, output from filter_network_expressed_genes.R
uniprot ids in columns "to" and "from", gene symbols in columns "source_genesymbol" and "target_genesymbol"
Table of differentially expressed genes (csv, output from deseq2)
ID type of the differentially expressed genes - uniprot or gene symbols

Output: Network in same format as input network (but tab seperated), filtered to include only interactions where target nodes are differentially expressed.

```r
##### Set up #####

# Capture  messages and errors to a file.
zz <- file("virallink.out", open="a")
sink(zz, type="message", append = TRUE)
message("\nStarting reg-deg network script: get_regualtor_deg_network.R\n")

# Install required packages
if (!requireNamespace("tidyverse", quietly = TRUE)) 
  install.packages("tidyverse", repos = "https://cloud.r-project.org")

# Load required packages
library(tidyverse)

# Define parameters
args <- commandArgs(trailingOnly = TRUE)

# Check length of command line parameters
if (length(args) != 4){
  stop("Wrong number of command line input parameters. Please check.")
}

# Output directory
outdir <- args[4]

# Create output dir if required
path <- file.path(outdir, "2_process_a_priori_networks")
dir.create(path, showWarnings = FALSE, recursive=TRUE)

# Contextualised regulatory network
reg_net <- read.csv(args[1], sep = "\t")

# Differentially expressed genes (prefiltered)
diff_genes <- read.csv(args[2])

# ID type of the differentially expressed genes - uniprot or gene symbols
id_type <- args[3] # either "uniprot" or "symbol"
  
##### Preprocess #####

# Get column names of network to match to the differentially expressed genes (based on id type)
if(id_type == "symbol"){
  source_col <- "source_genesymbol"
  target_col <- "target_genesymbol"
} else if(id_type == "uniprot"){
  source_col <- "to"
  target_col <- "from"
} else {
  stop("The differential expression data id type is not correctly specified. Should be \"uniprot\" or \"symbol\"")
}

##### Filter network #####

# Filter netowrk so all target genes are differentially expressed
reg_net_f <- reg_net %>% filter(get(target_col) %in% diff_genes$Gene) %>% unique()

# Get list of regulators with number of targeted DEGs
regs <- reg_net_f %>% select(deg_regs = !!source_col) %>% add_count(deg_regs, name = "num_degs") %>% unique()

##### Save #####

write.table(reg_net_f, file = file.path(path,"contextualised_regulator-deg_network.txt"), sep = "\t", quote = F, row.names = F)
write.table(regs, file = file.path(path, "contextualised_regulators_of_degs.txt"), sep = "\t", quote = F, row.names = F)

# reset message sink and close the file connection
sink(type="message")
close(zz)
```

## 3_network_diffusion

### prepare_tiedie_input.R

Script to prepare the input files for TieDIE

Input: 
Contextualised PPI network output from 'filter_network_expressed_genes.R'
TF-DEGs network file output from 'get_regulator_deg_network.R'
Differentially expressed gene file output from 'diff_expression_deseq2.R'
Viral protein- human binding protein file (provided - based on Gordon et al.)

Output: 
pathway.sif - sif format contextualised PPI (Omnipath) network
upstream.input - text file containing the upstream genes with a weight and direction (human binding partners)
downstream.intput - text file containing the downstream genes with a weight and direction (TFs)

```r
##### Set up #####

# Capture  messages and errors to a file.
zz <- file("virallink.out", open="a")
sink(zz, type="message", append = TRUE)
message("\nStarting tiedie input preparation: prepare_tiedie_input.R\n")

# Install required packages
if (!requireNamespace("tidyverse", quietly = TRUE)) 
install.packages("tidyverse", repos = "https://cloud.r-project.org")

# Load packages
library(tidyverse)

# Define parameters
args <- commandArgs(trailingOnly = TRUE)

# Check length of command line parameters
if (length(args) != 5){
  stop("Wrong number of command line input parameters. Please check.")
}

outdir <- args[5]

# Contextualised PPI network (omnipath) output from 'filter_network_expressed_genes.R'
ppis <- read.csv(args[1], sep = "\t")

# TF-DEG interactions output from 'get_regulator_deg_network.R'
tfs <- read.csv(args[2], sep = "\t")

# Filtered differential expression file output from 'diff_expression_deseq2.R'
degs <- read.csv(args[3])

# Human binding proteins of viral proteins
hbps <- read.csv(args[4], sep = "\t")
  
# Create output dir if required
path <- file.path(outdir, "3_network_diffusion", "input_files")
dir.create(path, showWarnings = FALSE, recursive=TRUE)


##### Process pathways #####

# Convert to sif format
ppis2 <- ppis %>% mutate(direction = ifelse(consensus_stimulation == "1", "stimulates>", "inhibits>")) %>%
  select(c(from,direction, to)) %>%
  unique()

# Save pathways
write.table(ppis2, file = file.path(path,"pathway.sif"), sep = "\t", col.names = F, row.names = F, quote = F)

##### Process upstream input #####

# Number of viral proteins bound to each human protein is the weight
# If sign of interaction given then use it, else assume all inhibitory ("-")
if ("sign" %in% colnames(hbps)){
  # Check values in sign column are "-" or "+"
  hbps_f <- hbps %>% filter((sign == "-") | (sign == "+"))
  if (nrow(hbps_f) != nrow(hbps)){
    print("WARNING: Some of the viral-human binding protein interactions were disgarded as the values in the 'sign' column were not '+' or '-'.")
  }
  if (nrow(hbps_f) == 0){
    stop("ERROR: viral-human binding protein interactions do not have the correct values in 'sign' column. They should be '+' or '-'.")
  }
  hbps2 <- hbps %>% select(human_protein, sign) %>% dplyr::rename(direction=sign) %>% group_by(human_protein,direction) %>% summarise(n = n()) %>%
    select(human_protein, n, direction)
} else {
  hbps2 <- hbps %>% select(human_protein) %>% group_by(human_protein) %>% summarise(n = n()) %>% mutate(direction = "-")
}

# Save upstream data
write.table(hbps2, file = file.path(path,"upstream.input"), sep = "\t", col.names = F, row.names = F, quote = F)

##### Process downstream input #####

# (1/#targetgenes) * sum(lfc(targetgene)*signofint)

# Join the tf-deg network with the deg lfc values
tfs2 <- left_join(tfs, degs, by =c("target_genesymbol"="Gene")) %>% select(c(from, to, consensus_stimulation, log2FoldChange)) %>%
  mutate(lfc_sign = ifelse(consensus_stimulation == "1", log2FoldChange, -log2FoldChange))

# Get the sum of all lfc*sign values - and the number of target genes for each tf
tfs3 <- tfs2 %>% select(from, lfc_sign) %>% group_by(from) %>% summarise(sumof = sum(lfc_sign), n = n())

# Divide sumof by n and determine sign (based on sign of the value)
tfs4 <- tfs3 %>% mutate(final_val = sumof/n) %>% mutate(sign = ifelse((final_val >= 0), "+", "-")) %>%
  select(from, final_val, sign)

# Save downstream data
write.table(tfs4, file = file.path(path,"downstream.input"), sep = "\t", col.names = F, row.names = F, quote = F)

# reset message sink and close the file connection
sink(type="message")
close(zz)

```

### tiedie.py

Script to run diffusion analysis using TieDIE

**Minimum Inputs:**
* Separate source/target input heat files: tab-separated, 3 columns each with <gene> <input heat> <sign (+/-)> (as created in the prepare_tiedie_input.R script)
* A search pathway in *.sif* format (geneA <interaction> geneB)

**Outputs:**
* Creates a directory in the current working directory, and writes all output to that
* Information and warnings are logged to standard error

## 4_create_network

### combined_edge_node_tables.R

Script to combine the covid subnetworks into 1 network table and one node table

Input: 
Receptor-tf network output from tiedie
SARS-CoV-2 - human binding protein network (columns 'viral_protein', 'human_protein')
Contextualised regulatory network (TFs - DEGs) (incl. columns 'DEG',"	TF", "consensus_stimulation") output from 'get_regulator_deg_network.R'
Heats values output from TieDie
Sars proteins gene symbol to uniprot conversion table
Differential expression table (unfiltered)

Output: 
Network file where each line represents an interaction
Node table where each lines represents one node in the network

```r
##### Setup #####

# Capture  messages and errors to a file.
zz <- file("virallink.out", open="a")
sink(zz, type="message", append = TRUE)
message("\nStarting network reconstruction: combined_edge_node_tables.R\n")

# Install required packages
if (!requireNamespace("tidyverse", quietly = TRUE)) 
  install.packages("tidyverse", repos = "https://cloud.r-project.org")
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) 
  install.packages("org.Hs.eg.db", repos = "https://cloud.r-project.org")

library(tidyverse)
library(org.Hs.eg.db)

# Define parameters
args <- commandArgs(trailingOnly = TRUE)

# Check length of command line parameters
if (length(args) != 8){
  stop("Wrong number of command line input parameters. Please check.")
}

outdir <- args[8]

# Receptor-tf network from tiedie
rec_tf <- read.csv(args[1], header=F, col.names=c("Source.node","Relationship","Target.node"), sep = "\t")
# Node heats from tiedie
heats <- read.csv(args[2], sep="=")

# Virus-receptor network
hbps <- read.csv(args[3], sep = "\t")
# SARS-Cov2 gene symbols
sars <- read.csv(args[4], sep = "\t")

# tf-deg network
tf_deg <- read.csv(args[5], sep = "\t")

# lfc table (unfiltered)
lfc <- read.csv(args[6])

# ID type of the lfc table (uniprot or symbol)
id_type <- args[7] # "symbol" or "uniprot"

# Create output dir if required
path <- file.path(outdir, "4_create_network")
dir.create(path, showWarnings = FALSE, recursive=TRUE)

##### Process network edges #####

rec_tf2 <- rec_tf %>% mutate(layer= "bindingprot-tf") 

if("sign" %in% colnames(hbps)){
  hbps2 <- hbps %>% mutate(Relationship = ifelse(sign == "-","inhibits>", ifelse(sign=="+", "stimulates>", "unknown")))  %>%
    select(-c(sign)) %>%
    mutate(layer = "cov-bindingprot") %>%
    dplyr::rename(Source.node = viral_protein, Target.node = human_protein) %>%
    filter(Target.node %in% rec_tf2$Source.node)
} else {
  hbps2 <- hbps %>% mutate(Relationship = "unknown", layer = "cov-bindingprot") %>%
    dplyr::rename(Source.node = viral_protein, Target.node = human_protein) %>%
    filter(Target.node %in% rec_tf2$Source.node)
}

tf_deg2 <- tf_deg %>% mutate(layer = "tf-deg") %>%
  dplyr::select(Target.node = to, Source.node = from, Relationship = consensus_stimulation, layer) %>%
  filter(Source.node %in% rec_tf2$Target.node) %>%
  mutate(Relationship = str_replace(Relationship, "1", "stimulates>")) %>%
  mutate(Relationship = str_replace(Relationship, "0", "inhibits>"))

# Join together
whole_net <- rbind(hbps2, rec_tf2, tf_deg2) %>% unique()

# Save
write.table(whole_net, file = file.path(path, "final_network.txt"), sep = "\t", quote = F, row.names = F)

##### Process node table #####

# Get layers of all nodes
nodes1 <- whole_net %>% filter(layer == "cov-bindingprot") %>% mutate(cov_layer = "cov2") %>% dplyr::select(node = Source.node, cov_layer) %>% unique()
nodes2 <- whole_net %>% filter(layer == "cov-bindingprot") %>% mutate(bindingprot_layer = "bindingprot") %>% dplyr::select(node = Target.node, bindingprot_layer) %>% unique()
nodes3 <- whole_net %>% filter(layer == "bindingprot-tf") %>% mutate(ppi1_layer = "bindingprot and/or protein") %>% dplyr::select(node = Source.node, ppi1_layer) %>% unique()
nodes4 <- whole_net %>% filter(layer == "bindingprot-tf") %>% mutate(ppi2_layer = "protein and/or tf") %>% dplyr::select(node = Target.node, ppi2_layer) %>% unique()
nodes5 <- whole_net %>% filter(layer == "tf-deg") %>% mutate(tf_layer = "tf") %>% dplyr::select(node = Source.node, tf_layer) %>% unique()
nodes6 <- whole_net %>% filter(layer == "tf-deg") %>% mutate(deg_layer = "deg") %>% dplyr::select(node = Target.node, deg_layer) %>% unique()

# Join node layers together
all_nodes <- full_join(nodes1, nodes2) %>% full_join(nodes3) %>% full_join(nodes4) %>% full_join(nodes5) %>% full_join(nodes6)

# Combine the ppi layer col
all_nodes <- all_nodes %>% mutate(ppi_layer = ifelse((ppi1_layer == "bindingprot and/or protein" & ppi2_layer == "protein and/or tf"), "protein", "NA")) %>%
  dplyr::select(-c(ppi1_layer, ppi2_layer))

# Combine into 1 column
all_nodes <- all_nodes %>% unite(all_nodes, cov_layer, bindingprot_layer, ppi_layer, tf_layer, deg_layer, sep = ";", remove=FALSE, na.rm = TRUE)
all_nodes[is.na(all_nodes)] <- "NA"
rm(nodes1,nodes2, nodes3, nodes4, nodes5, nodes6)

# Get human node id conversion from Orgdb
id_mapping <- select(org.Hs.eg.db, keys=all_nodes$node, columns=c("SYMBOL","ENTREZID"), keytype="UNIPROT")
# Get only the first mapped id - as need 1:1 mapping later on
id_mapping2 <- id_mapping[!(duplicated(id_mapping$UNIPROT)),]
# Add id conversions to node table
all_nodes <- left_join(all_nodes, id_mapping2, by = c("node"="UNIPROT"))

# Add sars node id conversions
sars_sym <- sars %>%dplyr::select(c(node = Accession, gene_symbol))
all_nodes <- left_join(all_nodes, sars_sym, by = c("node")) %>% mutate(SYMBOL = replace_na(SYMBOL,"")) %>% mutate(gene_symbol = replace_na(as.character(gene_symbol),""))
all_nodes <-all_nodes %>% mutate(gene_symbol = paste0(SYMBOL, gene_symbol)) %>% dplyr::select(-c(SYMBOL))

# Add node heats from tiedie
heats$node <- row.names(heats)
heats$node <- gsub('\\s+', '', heats$node)
all_nodes <- left_join(all_nodes,heats)

# Add lfc and adj p value for all
if (id_type == "symbol"){
  lfc2 <- lfc %>% dplyr::select(gene_symbol = X, log2FoldChange, padj)
} else if(id_type == "uniprot"){
  lfc2 <- lfc %>% dplyr::select(UNIPROT = X, log2FoldChange, padj)
}
all_nodes <- left_join(all_nodes, lfc2)

# Save node table
write.table(all_nodes, file = file.path(path, "node_table.txt"), sep = "\t", quote = F, row.names = F)

# reset message sink and close the file connection
sink(type="message")
close(zz)
```

## 5_betweenness_and_cluster_analysis

### betweenness_and_clustering.R

Script to calculate betweenness centrality and identify clusters (using MCODE) in directed causal network.

REQUIRES CYTOSCAPE PROGRAM TO BE OPEN TO RUN CLUSTERING - if not open the clustering will be skipped

Input: 
Directed causal network output from 'combined_edge_node_tables.R'
Node table output from 'combined_edge_node_tables.R'

Output: 
Updated node table with new column for betweenness centrality and, if MCODE successfully ran, 
3 new columns with MCODE results.

```r

##### Setup #####

# Capture  messages and errors to a file.
zz <- file("virallink.out", open="a")
sink(zz, type="message", append = TRUE)
message("\nStarting betweenness + clustering: betweenness_and_clustering.R\n")

# Install requried packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!requireNamespace("RCy3", quietly = TRUE)) 
  BiocManager::install("RCy3")
if (!requireNamespace("tidyverse", quietly = TRUE)) 
  install.packages("tidyverse", repos = "https://cloud.r-project.org")
if (!requireNamespace("igraph", quietly = TRUE)) 
  install.packages("igraph", repos = "https://cloud.r-project.org")

# Load packages
library(tidyverse)
library(RCy3)
library(igraph)

# Define parameters
args <- commandArgs(trailingOnly = TRUE)

# Check length of command line parameters
if (length(args) != 3){
  stop("Wrong number of command line input parameters. Please check.")
}

outdir <- args[3]

# Directed causal network output from 'combined_edge_node_tables.R'
net <- read.csv(args[1], sep = "\t")
  
# Node table output from 'combined_edge_node_tables.R'
nodes <- read.csv(args[2], sep = "\t")

# Create output dir if required
path <- file.path(outdir, "5_betweenness_and_cluster_analysis")
dir.create(path, showWarnings = FALSE, recursive=TRUE)

##### Betweenness centrality #####

# Convert to igraph network
i_net <- graph_from_data_frame(net, directed=TRUE)

# Calculate node betweenness values
betweenness_res <- betweenness(i_net, v = V(i_net), directed = TRUE, weights = NULL,
            nobigint = TRUE, normalized = FALSE)

# Get into dataframe
betweenness_res2 <- data.frame(betweenness_centrality=sort(betweenness_res, decreasing=TRUE)) %>% tibble::rownames_to_column()

# Add to node table
nodes_2 <- left_join(nodes, betweenness_res2, by = c("node"="rowname"))

##### Clustering #####

cytoscape_func <- function(i_net){
  # Function to carry out MCODE clustering and export MCODE node annotations
  # Only to be run if already checked that the script can connect to cytoscape >= v7
  # Returns node table with MCODE results - or empty table if couldn't get MCODE running
  
  # Install MCODE all
  install_error <- FALSE
  tryCatch({installApp("MCODE")}, error=function(e) {install_error <<- TRUE})
  
  # Only continue if success installing MCODE, else print error
  if (install_error){
    message('Error installing MCODE app - skipping MCODE clustering')
    clust_res <- ""
  } else {

    # Import network into igraph
    createNetworkFromIgraph(i_net,"whole_network")

    # Carry out clustering
    commandsRun('mcode cluster degreeCutoff=2 fluff=FALSE haircut=TRUE nodeScoreCutoff=0.2 kCore=2 maxDepthFromStart=100')
    
    # Wait for Cytoscape to catch up
    Sys.sleep(180)
    
    ###### PROBLEMS ##########
    
    # Get results as table ## POSSIBLE BUG IN RCY3 as the MCODE_Cluster col is mixed up???
    #clust_res <- getTableColumns(table="node",columns=c("name", "MCODE_Node_Status", "MCODE_Score","MCODE_Cluster"))
    
    # Instead get just the clustered/seed nodes and export node table - downside of this appraoch is we can't keep the cluster 
    # values for the unclustered data 
    selectNodes("unclustered", by.col = "MCODE_Node_Status", network="whole_network")
    deleteSelectedNodes(network = "whole_network")
    
    # Creating subnetwork also doesn't appear to work  - another bug?
    #selectEdgesConnectingSelectedNodes(network="whole_network")
    #createSubnetwork(subnetwork.name = "clustered_nodes", network="whole_network")
    
    # Export subnetwork node table
    clust_res <- getTableColumns(table="node",columns=c("name", "MCODE_Node_Status", "MCODE_Score","MCODE_Cluster"))
    
  }
  
  # Close cytoscape session
  closeSession(FALSE)
  
  return(clust_res)
}

# See if cytoscape is open and >= v7
cyto_error <- FALSE
tryCatch( {msg <- cytoscapePing () } , error = function(e) {cyto_error <<- TRUE})

if (!cyto_error){
  if (cytoscapeVersionInfo ()[2] >= 3.7){
    continue = TRUE
    message('Successfully connected to Cytoscape - carrying out MCODE clustering and creating Cytoscape file')
  } else {
    continue = FALSE
    message('Successfully connected to Cytscape BUT version not >= 3.7 - skipping MCODE clustering and creation of Cytoscape file')
  }
} else {
  continue = FALSE
  message('Could not connect to Cytoscape - skipping MCODE clustering and creation of Cytoscape file')
}

# Run clustering if cytoscape open and new enough - and RCy3 app new enough (error getting MCODE results if 'its not')
if(continue & (packageVersion("RCy3") >= "2.6.0")){
  
  # Call function to run MCODE
  results <- cytoscape_func(i_net)

  # If there are results in the table then MCODE ran
  if(results != ""){
    # Un-list the MCODE_cluster column
    results2 <- results %>% unnest(MCODE_Cluster) %>% group_by(name) %>% mutate(MCODE_Cluster = paste0(MCODE_Cluster, collapse = ";")) 
    
    # Join the mcode cluster results to the node table
    nodes_3 <- left_join(nodes_2, results2, by = c("node"="name"))
    
  } else {
    # If MCODE couldn't install then there won't be any results
    nodes_3 <- nodes_2
  }
} else {
  # If Cytoscape wasn't open (or too old) then there won't be any results
  nodes_3 <- nodes_2
}

##### Save output #####

# Save updated node table
write.table(nodes_3, file=file.path(path, "node_table_betweenness_clusters.txt"), row.names = FALSE, quote = FALSE, sep = "\t")

# reset message sink and close the file connection
sink(type="message")
close(zz)

``` 

### cytoscape_visualisation.R

Script to create Cytoscape file containing each of the networks (whole causal network and network clusters - if calculated). 
REQUIRES CYTOSCAPE PROGRAM TO BE OPEN - will skip otherwise.

Input: 
Directed causal network output from 'combined_edge_node_tables.R'
Node table output from 'betweenness_and_clustering.R' (with or without cluster labels)

Output: Cytoscape file containing the whole causal network and cluster networks (if cluster annotations calculated previously)

```r
##### Setup #####

# Capture  messages and errors to a file.
zz <- file("virallink.out", open="a")
sink(zz, type="message", append = TRUE)
message("\nStarting visualisation: cytoscape_visualisation.R\n")

# Install required packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!requireNamespace("RCy3", quietly = TRUE)) 
  BiocManager::install("RCy3", repos = "https://cloud.r-project.org")
if (!requireNamespace("tidyverse", quietly = TRUE)) 
  install.packages("tidyverse", repos = "https://cloud.r-project.org")

# Load packages
library(tidyverse)
library(RCy3)

# Define parameters
args <- commandArgs(trailingOnly = TRUE)

# Check length of command line parameters
if (length(args) != 3){
  stop("Wrong number of command line input parameters. Please check.")
}

# Output directory
outdir <- args[3]

# Directed causal network output from 'combined_edge_node_tables.R'
net <- read.csv(args[1], sep = "\t")

# Node table output from 'betweenness_and_clustering.R'
nodes <- read.csv(args[2], sep = "\t")

# Create output dir if required
path <- file.path(outdir, "5_betweenness_and_cluster_analysis")
dir.create(path, showWarnings = FALSE, recursive=TRUE)

##### Preprocess #####

net2 <- net %>% dplyr::rename("source"= "Source.node", "target"="Target.node", "interaction"="Relationship")
nodes2 <- nodes %>% dplyr::rename("id" = "node")

# Determine if clustering has been carried out
if("MCODE_Cluster" %in% colnames(nodes2)){
  cluster = TRUE
  clus_filt <- nodes2 %>% dplyr::select(MCODE_Cluster) %>% group_by(MCODE_Cluster) %>% count() %>% filter(MCODE_Cluster != "NA")
} else {
  cluster = FALSE
}

# Determine if betweenness centrality calcualtions have been carried out
if("betweenness_centrality" %in% colnames(nodes2)){
  bet = TRUE
} else {
  bet=FALSE 
}

##### Import into Cytoscape #####

# See if cytoscape is open
cyto_error <- FALSE
tryCatch( {msg <- cytoscapePing () } , error = function(e) {cyto_error <<- TRUE})

if (!cyto_error){
  continue = TRUE
  message('Successfully connected to Cytoscape - continuing with visualisation')
} else {
  continue = FALSE
  message('Could not connect to Cytoscape - skipping visualisation')
}

# Run visualisation if cytoscape open
if(continue){
  # Add network to cytoscape - sometimes not working as in 'https://github.com/cytoscape/RCy3/issues/81'
  #createNetworkFromDataFrames(nodes2,net2, title="causal network")
  # Here is a workaround
  createNetworkFromDataFrames(nodes2, edges=net2[,1:2], title="causal network")
  edges <- net2 %>%
    mutate(key=paste(source, "(interacts with)", target))
  loadTableData(edges[,3:5],data.key.column = 'key', table = 'edge')
  
  # Colour by Betweenness centrality if data exists in network
  if (bet){
    setNodeColorMapping("betweenness_centrality",mapping.type="c", network="causal network")
  }
  
  # Subnetwork for each cluster
  if (cluster){
    for (i in clus_filt$MCODE_Cluster){
      clearSelection(type = "both", network = "causal network")
      selectNodes(i, by.col = "MCODE_Cluster", network="causal network")
      selectEdgesConnectingSelectedNodes(network="causal network")
      createSubnetwork(nodes="selected", subnetwork.name = i,network = "causal network")
      layoutNetwork()
    }
  }
  
  # Save cys file
  saveSession(filename = file.path(path, "causal_network.cys"))
}

# reset message sink and close the file connection
sink(type="message")
close(zz)
``` 

### Results ###




## 6_functional_analysis

### cluster_functional_analysis.R

Functional analysis of human genes from clusters inside network
Overrepresentation analysis - reactome and GO BP - using contextualised specific network as background
NB. GO analyses carried out using uniprot IDs for the ppi nodes, but for the Reactome analysis it was necessary to convert to ENTREZ
For the overenrichment analysis (q val <= 0.05) it uses all nodes of the contextualised specific networks (expressed omnipath) as the background.
For the GO analyses it uses simplify with parameter 0.1

Input: 
Node table (csv file) output from cytoscape containing the cluster annotations in the column "MCODE_cluster"
PPI background network file - expressed omnipath network

Output: 
For each cluster (with >= 15 nodes): table of results, dot plot and map plot for GO and Reactome overrepresentation
For all clusters (with >= 15 nodes) together (compareCluster function) dot plot for GO and Reactome overrepresentation

``` r
##### Set up #####

# Capture  messages and errors to a file.
zz <- file("virallink.out", open="a")
sink(zz, type="message", append = TRUE)
message("\nStarting cluster functional analysis: cluster_functional_analysis.R\n")

# Install requried packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!requireNamespace("ReactomePA", quietly = TRUE)) 
  BiocManager::install("ReactomePA", repos = "https://cloud.r-project.org")
if (!requireNamespace("clusterProfiler", quietly = TRUE)) 
  BiocManager::install("clusterProfiler", repos = "https://cloud.r-project.org")
if (!requireNamespace("tidyverse", quietly = TRUE)) 
  install.packages("tidyverse", repos = "https://cloud.r-project.org")
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) 
  BiocManager::install("org.Hs.eg.db", repos = "https://cloud.r-project.org")

# load packages
library(tidyverse)
library(ReactomePA)
library(clusterProfiler)
library(org.Hs.eg.db)

# Define parameters
args <- commandArgs(trailingOnly = TRUE)

# Check length of command line parameters
if (length(args) != 3){
  stop("Wrong number of command line input parameters. Please check.")
}

outdir <- args[3]

# node table containing genes in column gene_symbol and cluster number in 'MCODE_Cluster' - output from 'betweenness_and_clustering.R'
nodes <- read.csv(args[1], sep ="\t")

# Get background network (contextualised specific omnipath)
backgr <- read.csv(args[2], sep = "\t")

# Create output dir if required
path <- file.path(outdir, "6_functional_analysis", "clusters")
dir.create(path, showWarnings = FALSE, recursive=TRUE)

##### Preprocess files #####

preproc_backgr <- function(backgr){
  # Get all nodes in background network
  
  source <- backgr %>% dplyr::select(node = from) %>% unique()
  target <- backgr %>% dplyr::select(node = to) %>% unique()
  nodes <- rbind(source, target) %>% unique()
  
  return(nodes)
}

##### Enrichment analysis #####

go_overrep <- function(net_heats, back_nodes, name, id,folder){
  # Carries out GO overrepresentation analysis.

  # Carry out normal GO gene overrepresentation analysis
  go1 <- enrichGO(gene = as.character(net_heats), OrgDb = "org.Hs.eg.db", keyType = id, universe =as.character(back_nodes$node), ont = 'BP', qvalueCutoff = 0.05)
  
  if (nrow(go1) > 0) {
    # Remove redundancy of GO terms. Cutoff refers to similarity
    go2 <- clusterProfiler::simplify(go1, cutoff=0.1, by = "qvalue", select_fun=min)
    
    # Get as dataframe
    go2_df <- as.data.frame(go2)
    
    if (nrow(go1) >1){
      # Get enrichment map
      map_ora <- emapplot(go2)
      # Save map
      filen <- file.path(folder, paste0(name, "_go_0.1_overrep_map.pdf"))
      pdf(filen)
      print(map_ora)
      dev.off()
    }
    
    # Get dot plot
    dot_plot <- dotplot(go2, showCategory=10, orderBy="qvalue", font.size = 10)
    # Save dot plot
    filep <- file.path(folder, paste0(name, "_go_0.1_overrep_dot.pdf"))
    pdf(filep)
    print(dot_plot)
    dev.off()
   
  } else{
    go2_df <- ""
  }
  
  return(go2_df)
}


reactome_overrep <- function(net_heats_e, back_nodes_e, name, id,folder){
  # Carries out reactome overrepresentation analysis.
  
  # Carry out normal reactome gene overrepresentation analysis (can only use entrez :( )
  re1 <- enrichPathway(gene = as.character(net_heats_e$ENTREZ), organism = "human", universe =as.character(back_nodes_e$ENTREZ), qvalueCutoff = 0.05)
  
  if (nrow(re1) > 0) {
   
    # Get as dataframe
    re1_df <- as.data.frame(re1)
    
    if (nrow(re1) >1){
      # Get enrichment map
      map_ora <- emapplot(re1)
      # Save map
      filen <- file.path(folder, paste0(name, "_reactome_overrep_map.pdf"))
      pdf(filen)
      print(map_ora)
      dev.off()
    }
    
    # Get dot plot
    dot_plot <- dotplot(re1, showCategory=10, orderBy="qvalue", font.size = 10)
    # Save dot plot
    filep <- file.path(folder, paste0(name, "_reactome_overrep_dot.pdf"))
    pdf(filep)
    print(dot_plot)
    dev.off()
  } else{
    re1_df <- ""  
  }
  
  return(re1_df)
}

##### Run all #####

# Preprocess background nodes
back_nodes <- preproc_backgr(backgr)

# Remove unclustered nodes
nodes2 <- nodes %>% filter(MCODE_Cluster != "")

# Un-list the MCODE_cluster column
nodes3 <- nodes2 %>% unnest(MCODE_Cluster) %>% group_by(node) %>% mutate(MCODE_Cluster = paste0(MCODE_Cluster, collapse = ";")) 

# List clusters with <15 nodes
clusters <- nodes3 %>% ungroup() %>% dplyr::select(MCODE_Cluster) %>% group_by(MCODE_Cluster) %>% summarise(n = n()) %>%
  filter(n>= 15)

# Select clusters of interest
nodes3 <- nodes3 %>% filter(MCODE_Cluster %in% clusters$MCODE_Cluster)

# Split df to get seperate clusters
cluster_dfs <- split(nodes3 , f = nodes3$MCODE_Cluster)

compare_ids <- vector(mode = "list")
n<-1
for (i in cluster_dfs){
  
  if (nrow(i) == 0){
    next
  }
  
  # Name of cluster
  nam <- names(cluster_dfs)[n]
  print(nam)
  # Run overenrichment analysis
  go_res <- go_overrep(i$node, back_nodes, nam, "UNIPROT",path)
  write.table(go_res, file = file.path(path, paste0(nam, "_go_0.1_overrep_results.txt")), quote = FALSE, row.names = FALSE, sep = "\t")
  
  # get Entrez ids
  # Convert to entrez genes
  i_e <- bitr(i$node, fromType='UNIPROT', toType='ENTREZID', OrgDb="org.Hs.eg.db")
  back_nodes_e <- bitr(back_nodes$node, fromType='UNIPROT', toType='ENTREZID', OrgDb="org.Hs.eg.db")
  compare_ids[[as.character(n)]] <- i_e$ENTREZID
  
  # Run reactome overenrichment
  reactome_res <- reactome_overrep(i_e, back_nodes_e, nam,"UNIPROT",path)
  write.table(reactome_res, file = file.path(path, paste0(nam, "_reactome_overrep_results.txt")), quote = FALSE, row.names = FALSE, sep = "\t")
  
  n <- n+1
}

# Run together
clusters_go <- compareCluster(compare_ids, fun="enrichGO", ont = 'BP', OrgDb='org.Hs.eg.db',pvalueCutoff=0.1, universe = back_nodes_e$ENTREZID)
#clusters_go2 <- simplify(clusters_go,cutoff = 0.7,by = "p.adjust")
dotplot_g <- dotplot(clusters_go, font.size = 10, title = paste0("NHBE - clustering of PPI nodes"), showCategory=5)
fileg <- file.path(path,"all_clusters_go_overrep_dot.pdf")
pdf(fileg, width=7)
print(dotplot_g)
dev.off()

clusters_react <- compareCluster(compare_ids, fun="enrichPathway",pvalueCutoff=0.1, universe =back_nodes_e$ENTREZID)
#summary(clusters_react)
dotplot_r <- dotplot(clusters_react, font.size = 10, title = paste0("NHBE - clustering of PPI nodes"), showCategory=5)
filer <- file.path(path,"all_clusters_reactome_overrep_dot.pdf")
pdf(filer, width=6)
print(dotplot_r)
dev.off()

# reset message sink and close the file connection
sink(type="message")
close(zz)
```

### network_functional_analysis.R

Functional overrepresentation analysis of causal networks: upstream signalling proteins (binding proteins-TFs inclusive), and differentially expressed genes (DEGs) seperately

Reactome and Gene ontology BP
NB. 
* GO analyses carried out using uniprot IDs, but for the Reactome analysis it is necessary to convert to ENTREZ IDs
* All nodes of the contextualised specific PPI network (expressed omnipath) are used as the background for the upstream signalling proteins
* All TGs of the contextualised specific TF-TG network (expressed dorothea) are used as the background for the DEGs significantly overrepresented functions have q val <= 0.05
* GO analyses use simplify to remove rudundant function (with parameter 0.1)

Input: 
whole network node file (output from combined_edge_node_tables.R)
background network file for ppis (contextualised specific PPI network output from filter_network_expressed_genes.R)
background network file for degs (contextualised specific TF-TG network output from filter_network_expressed_genes.R)

Output: Data table, dot plot and map plot for each overrepresentation analysis

```r
##### Set up #####

# Capture  messages and errors to a file.
zz <- file("virallink.out", open="a")
sink(zz, type="message", append = TRUE)
message("\nStarting functional analysis: network_functional_analysis.R\n")

# Install requried packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!requireNamespace("ReactomePA", quietly = TRUE)) 
  BiocManager::install("ReactomePA", repos = "https://cloud.r-project.org")
if (!requireNamespace("clusterProfiler", quietly = TRUE)) 
  BiocManager::install("clusterProfiler", repos = "https://cloud.r-project.org")
if (!requireNamespace("tidyverse", quietly = TRUE)) 
  install.packages("tidyverse", repos = "https://cloud.r-project.org")
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) 
  BiocManager::install("org.Hs.eg.db", repos = "https://cloud.r-project.org")

#packages
library(tidyverse)
library(ReactomePA)
library(clusterProfiler)
library(org.Hs.eg.db)

# Define parameters
args <- commandArgs(trailingOnly = TRUE)

# Check length of command line parameters
if (length(args) != 4){
  stop("Wrong number of command line input parameters. Please check.")
}

outdir <- args[4]

# Create output dir if required
path1 <- file.path(outdir, "6_functional_analysis", "ppi_layer")
path2 <- file.path(outdir, "6_functional_analysis", "deg_layer")
dir.create(path1, showWarnings = FALSE, recursive=TRUE)
dir.create(path2, showWarnings = FALSE, recursive=TRUE)

##### Preprocess files #####

preproc_backgr_ppi <- function(backgr){
  # Get all nodes in ppi background network
  
  source <- backgr %>% dplyr::select(node = from) %>% unique()
  target <- backgr %>% dplyr::select(node = to) %>% unique()
  nodes <- rbind(source, target) %>% unique()
  
  return(nodes)
}

##### Enrichment analysis #####

go_overrep <- function(net, back_nodes, name, id,folder){
  # Carries out GO overrepresentation analysis.

  # Carry out normal GO gene overrepresentation analysis
  go1 <- enrichGO(gene = as.character(net$node), OrgDb = "org.Hs.eg.db", keyType = id, universe =as.character(back_nodes$node), ont = 'BP', qvalueCutoff = 1, pvalueCutoff = 1)
  
  if (nrow(go1) > 0) {
    # Remove redundancy of GO terms. Cutoff refers to similarity
    go2 <- clusterProfiler::simplify(go1, cutoff=0.1, by = "qvalue", select_fun=min)
    
    # Get as dataframe
    go2_df <- as.data.frame(go2)
    
    if (nrow(go1) >1){
      # Get enrichment map
      map_ora <- emapplot(go2)
      # Save map
      filen <- file.path(folder,paste0(name,"_map_0.1_GO_overrep.pdf"))
      pdf(filen)
      print(map_ora)
      dev.off()
    }
    
    # Get dot plot
    dot_plot <- dotplot(go2, showCategory=10, orderBy="qvalue", font.size = 10)
    # Save dot plot
    filep <- file.path(folder,paste0(name,"_dot_0.1_GO_overrep.pdf"))
    pdf(filep)
    print(dot_plot)
    dev.off()
   
  } else{
    go2_df <- ""
  }
  
  return(go2_df)
}


reactome_overrep <- function(net, back_nodes, name, id,folder){
  # Carries out reactome overrepresentation analysis.
  
  # Convert to entrez genes
  net_e <- bitr(net$node, fromType=id, toType='ENTREZID', OrgDb="org.Hs.eg.db")
  back_nodes_e <- bitr(back_nodes$node, fromType=id, toType='ENTREZID', OrgDb="org.Hs.eg.db")
  
  # Carry out normal reactome gene overrepresentation analysis (can only use entrez :( )
  re1 <- enrichPathway(gene = as.character(net_e$ENTREZ), organism = "human", universe =as.character(back_nodes_e$ENTREZ), qvalueCutoff = 1, pvalueCutoff = 1)
  
  if (nrow(re1) > 0) {
   
    # Get as dataframe
    re1_df <- as.data.frame(re1)
    
    if (nrow(re1) >1){
      # Get enrichment map
      map_ora <- emapplot(re1)
      # Save map
      filen <- file.path(folder,paste0(name,"_map_reactome_overrep.pdf"))
      pdf(filen)
      print(map_ora)
      dev.off()
    }
    
    # Get dot plot
    dot_plot <- dotplot(re1, showCategory=10, orderBy="qvalue", font.size = 10)
    # Save dot plot
    filep <- file.path(folder,paste0(name, "_dot_reactome_overrep.pdf"))
    pdf(filep, width=6.5)
    print(dot_plot)
    dev.off()
  } else{
    re1_df <- ""  
  }
  
  return(re1_df)
}


##### Run all #####

# Read background files
backgr_ppi <- read.csv(args[2], sep = "\t")
backgr_deg <- read.csv(args[3], sep = "\t")
  
# Read nodes file
nodes <- read.csv(args[1], sep = "\t")

# Preprocess background nodes
back_nodes_ppi <- preproc_backgr_ppi(backgr_ppi)
back_nodes_deg <- backgr_deg %>% dplyr::select(node = to) %>% unique()

# Get all ppi nodes
nodes_ppi <- nodes %>% filter(ppi_layer== "protein"|bindingprot_layer== "bindingprot"|tf_layer== "tf") %>% dplyr::select(node)
# Get all deg nodes
nodes_degs <- nodes %>% filter(deg_layer == "deg") %>% dplyr::select(node)

# Run GO overenrichment analysis
go_res_ppi <- go_overrep(nodes_ppi, back_nodes_ppi, "ppis", "UNIPROT",path1)
write.table(go_res_ppi, file = file.path(path1,"ppis_go_0.1_overrep_results.txt"), quote = FALSE, row.names = FALSE, sep = "\t")
go_res_deg <- go_overrep(nodes_degs, back_nodes_deg, "degs", "UNIPROT",path2)
write.table(go_res_deg, file = file.path(path2, "degs_go_0.1_overrep_results.txt"), quote = FALSE, row.names = FALSE, sep = "\t")
  
# Run reactome overenrichment analysis
reactome_res_ppi <- reactome_overrep(nodes_ppi, back_nodes_ppi, "ppis","UNIPROT",path1)
write.table(reactome_res_ppi, file = file.path(path1, "ppis_reactome_overrep_results.txt"), quote = FALSE, row.names = FALSE, sep = "\t")
reactome_res_deg <- reactome_overrep(nodes_degs, back_nodes_deg, "degs","UNIPROT",path2)
write.table(reactome_res_deg, file = file.path(path2, "degs_reactome_overrep_results.txt"), quote = FALSE, row.names = FALSE, sep = "\t")

# reset message sink and close the file connection
sink(type="message")
close(zz)
```

### reformat_functional_result.R

Script to reformat functional enrichment results (output from network_functions_analysis.R) into table of gene annotations

Input: list of overrepresentation result tables to process
Output: Reformatted table for each input table - where rows are genes rather than functions.

```r
### Set up ###

# Capture  messages and errors to a file.
zz <- file("virallink.out", open="a")
sink(zz, type="message", append = TRUE)
message("\nStarting reformatting functional results: reformat_functional_result.R\n")

# Install requried packages
if (!requireNamespace("tidyverse", quietly = TRUE)) 
  install.packages("tidyverse", repos = "https://cloud.r-project.org")
if (!requireNamespace("reshape2", quietly = TRUE)) 
  install.packages("reshape2", repos = "https://cloud.r-project.org")

# Load packages
library(tidyverse)
library(tidyr)
library(reshape2)

# Define parameters
args <- commandArgs(trailingOnly = TRUE)

# Check length of command line parameters
if (length(args) != 1){
  stop("Wrong number of command line input parameters. Please check.")
}

outdir <- args[1]

# Results of the functional analysis
files <- c("6_functional_analysis/deg_layer/degs_go_0.1_overrep_results.txt", "6_functional_analysis/deg_layer/degs_reactome_overrep_results.txt",
           "6_functional_analysis/ppi_layer/ppis_go_0.1_overrep_results.txt", "6_functional_analysis/ppi_layer/ppis_reactome_overrep_results.txt")

# Create output dir if required
path <- file.path(outdir, "6_functional_analysis", "reformatted")
dir.create(path, showWarnings = FALSE, recursive=TRUE)

### ###

# Iterate files
for (f in files){
    
  # Get file name
    if (grepl("/", f, fixed = TRUE)){
      filen_list <- strsplit(f, "/")
      filen <- sapply(filen_list, tail, 1 )
    }
  
    # Open file
    data <- read.csv(file.path(outdir,f), sep = "\t")
    
   if("geneID" %in% colnames(data)){
     coln <- "geneID"
   } else {
     coln <- "core_enrichment"
   }
    
    # Filter for required cols and wide to long
    data2 <- data %>% dplyr::select(c(Description, coln)) %>% separate_rows(coln)
    data2$match <- 1
    
    # Pivot table
    data3 <- data2 %>% pivot_wider(names_from = Description, values_from = match, values_fn = list(match = length))
    
    # convert na to 0
    data3[is.na(data3)] <- 0
    
    # Save
    write.table(data3, file = file.path(path, paste0("reformatted_",filen)), sep = "\t", quote=F, row.names = F)
}

# reset message sink and close the file connection
sink(type="message")
close(zz)
```


