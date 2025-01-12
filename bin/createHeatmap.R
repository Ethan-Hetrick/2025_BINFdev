################################################
## LOAD LIBRARIES                             ##
################################################
################################################

library(optparse)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)

################################################
################################################
## PARSE COMMAND-LINE PARAMETERS              ##
################################################
################################################
option_list <- list(
  make_option(c("-i", "--input_file"), type="character", default=NULL, metavar="path", help="Input sample file"),
  make_option(c("-g", "--geneFunctions_file"), type="character", default=NULL, metavar="path", help="Gene Functions file."),
  make_option(c("-a", "--annoData_file"), type="character", default=NULL, metavar="path", help="Annotation Data file."),
  make_option(c("-p", "--outprefix"), type="character", default='projectID', metavar="string", help="Output prefix.")
)


opt_parser <- OptionParser(option_list=option_list)
opt        <- parse_args(opt_parser)

sampleInput=opt$input_file
geneInput=opt$geneFunctions_file
annoInput=opt$annoData_file
outprefix=opt$outprefix

testing="Y"
if (testing == "Y"){
  sampleInput="sampleData.csv"
  geneInput="geneFunctions.csv"
  annoInput="annoData.csv"
  outprefix="test"
}


if (is.null(sampleInput)){
  print_help(opt_parser)
  stop("Please provide an input file.", call.=FALSE)
}

################################################
################################################
## READ IN FILES##
################################################
################################################
sampleData=read.csv(sampleInput,row.names=1)
annoData=read.csv(annoInput,row.names=1)
geneFunctions=read.csv(geneInput,row.names=1)

################################################
## Encode Gene Functions for Annotation      ##
################################################

# Convert gene functions to a factor for annotation purposes
geneFunctions$gene_functions <- as.factor(geneFunctions$gene_functions)

# Prepare annotation data for the heatmap (gene functions)
geneFunctions_anno <- geneFunctions[, "gene_functions", drop = FALSE]
colnames(geneFunctions_anno) <- "Gene Functions"

################################################
################################################
## Set colors##
################################################
################################################
annoColors <- list(
  gene_functions = c("Oxidative_phosphorylation" = "#F46D43",
                     "Cell_cycle" = "#708238",
                     "Immune_regulation" = "#9E0142",
                     "Signal_transduction" = "beige", 
                     "Transcription" = "violet"), 
  Group = c("Disease" = "darkgreen",
            "Control" = "blueviolet"),
  Lymphocyte_count = brewer.pal(5, 'PuBu')
)

################################################
## Create a basic heatmap                     ##
################################################

# Ensure only numeric data from sampleData is used for the heatmap (exclude row names)
numericData <- as.matrix(sampleData)  # Convert data to numeric matrix

# Basic Heatmap using numeric data from sampleData
basic_heatmap <- pheatmap(
  numericData, 
  cluster_rows = TRUE, 
  cluster_cols = TRUE, 
  scale = "row", 
  clustering_distance_rows = "euclidean", 
  clustering_distance_cols = "euclidean", 
  clustering_method = "ward.D", 
  show_rownames = TRUE, 
  show_colnames = TRUE,
  filename = paste0("basic_heatmap_", outprefix, ".pdf")
)

################################################
## Create a complex heatmap                   ##
################################################

# Complex Heatmap using numeric data from sampleData and annotations from geneFunctions_anno
complex_heatmap <- pheatmap(
  numericData, 
  cluster_rows = TRUE, 
  cluster_cols = TRUE, 
  scale = "row", 
  clustering_distance_rows = "euclidean", 
  clustering_distance_cols = "euclidean", 
  clustering_method = "ward.D", 
  annotation_row = geneFunctions_anno,  # Annotate rows (genes) based on gene functions
  annotation_col = annoData,  # Use column annotations from the provided annotation data
  annotation_colors = annoColors,  # Use predefined annotation colors
  show_rownames = FALSE, 
  show_colnames = FALSE, 
  legend_breaks = c(min(numericData), median(numericData), max(numericData)),  # Break scale into three bins
  legend_labels = c("low", "medium", "high"), 
  filename = paste0("complex_heatmap_", outprefix, ".pdf")
)
