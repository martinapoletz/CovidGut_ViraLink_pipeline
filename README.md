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

## 3_network_diffusion

### prepare_tiedie_input.R

### tiedie.py

## 4_create_network

### combined_edge_node_tables.R

## 5_betweenness_and_cluster_analysis

### betweenness_and_clustering.R

### cytoscape_visualisation.R

## 6_functional_analysis

### cluster_functional_analysis.R

### network_functional_analysis.R

### reformat_functional_result.R




