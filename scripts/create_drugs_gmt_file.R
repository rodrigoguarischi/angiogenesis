### creates_drugs_gmt_file.R ######################################################################
# Read CSV file from drugbank and transform it on a GMT file used on post enrichment analysis
#

### PREAMBLE ######################################################################################

# Load required libraries and datasets

# Define working directory
setwd( "~/angiogenesis_transcriptome/scripts/" );

# Disable automatic conversion of strings to factors
options(stringsAsFactors = FALSE)

### FUNCTIONS #####################################################################################

### DATA ANALYSIS #################################################################################

drugs.gmt <- list();

drugs.gmt$csv <- read.csv("aux_files/drugbank_12.04.2017/all/pharmacologically_active.csv");

for( line.number in 1:nrow(drugs.gmt$csv) ){

  drugs.gmt$drug.id2targets <- rbind(
    drugs.gmt$drug.id2targets,
    cbind( drugs.gmt$csv[line.number,]$Gene.Name, unlist( strsplit( drugs.gmt$csv[line.number,]$Drug.IDs, "; ") ) )
    );
  
  }

colnames( drugs.gmt$drug.id2targets ) <- c("gene.symbol", "drug.id");

drugs.gmt$drug.id2targets <- aggregate( gene.symbol~drug.id, drugs.gmt$drug.id2targets, paste, collapse = "\t" );

write.table(
  x = drugs.gmt$drug.id2targets[,c(1,1,2)],
  file = "../pathway_analysis/post_analysis/drugs2targets.gmt",
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
  );

### SAVE DATA #####################################################################################

### Step 4 - Save Rdata object for later use
save( list = c("drugs.gmt"), file = "aux_files/drugs_gmt.Rdata" );

# # Load expression data info
# load( file = "aux_files/drugs_gmt.Rdata" );
