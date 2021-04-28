### run_DE_analysis.R ############################################################################
# This script consists of four steps: 
# 1) load DESeqDataSet object using build_ddsFull.R script,
# 2) run differential expression (DE) analysis using DESeq2,
# 3) export data on CSV format,
# 4) save Rdata objects for later use

### PREAMBLE ######################################################################################
# Settings for environment

# Define working directory
setwd( "~/angiogenesis_transcriptome/scripts/" );

# Setting options
# Disable automatic conversion of strings to factors
options(stringsAsFactors = FALSE);
# Increase Java max memory. Solve problems when writing large xlsx files
options(java.parameters = "-Xmx6500m");

# Load required libraries and datasets
library( DESeq2 );
library( biomaRt );
library( XML );
library( RCurl );
library( xlsx ); # http://stackoverflow.com/questions/12872699/error-unable-to-load-installed-packages-just-now

### FUNCTIONS #####################################################################################

# Load function mergeDataframesByRownames
source("mergeDataframesByRownames.R");

### Function - add.info() #######################################################################
# Input variables: A vector with EnsemblIDs, desired field name and a biomart object
#
# Output variables: A vector with numeric values
#
# Description: Based on a list of ensemblIDs, query BioMart and return the desired field
add.info <- function( ensemblIDs, field.name, mart ){
  
  # Retrieve information from BioMart
  results <- getBM(
    attributes = c("ensembl_gene_id", field.name),
    filters = "ensembl_gene_id",
    values = ensemblIDs,
    mart = mart
    );
  
  # Collapse multiple rows by ensemblID
  results <- aggregate( .~ensembl_gene_id, na.action = NULL, data = results, paste, collapse = "; " );
  
  # Assign ensemblID as rownames
  rownames( results ) <- results$ensembl_gene_id; 
  
  # Order results to the same order as input
  results <- results[ensemblIDs, ];
  
  return( results[, field.name] );
  }

### DATA ANALYSIS #################################################################################

### Step 1 - Call build_ddsFull.R script to get a DESeqDataSet object
source( "build_ddsFull.R" );

# Initialize object 'de' that will hold data from differential expression
de <- list();

# Initialize object 'tmp' to hold all temporary objects
tmp <- list();

### Step 2 - Run DESeq2 pipeline and get DE genes for each time
de$geneIDs <- list();

for ( current.time in levels(ddsFull$Time) ){
 
  print( paste0( "Calculating DE genes for sample: R", current.time ) );
  
  # subset ddsFull
  tmp$dds <- ddsFull[ , ddsFull$Treat == 'P' | ddsFull$Time == current.time ];

  # If testing R12.5 a error stating "not full rank error" will arise. 
  # Due to experimental design we can consider R12.5 as R12 and this will fix this error
  if( current.time == 12.5 ){
    colData(tmp$dds)$Time <- factor(
      gsub(
        pattern = "12.5",
        replacement = "12",
        x = colData(tmp$dds)$Time
        )
      )
    };

  # "refactor" the factors, in case that levels have been dropped
  tmp$dds$Treat <- droplevels(tmp$dds$Treat);
  tmp$dds$Time  <- droplevels(tmp$dds$Time);
  
  # Make sure that control is the first level in the treatment factor
  tmp$dds$Treat <- relevel(tmp$dds$Treat, "P");
  tmp$dds$Time <- relevel(tmp$dds$Time, "12");
  
  # Remove rows which have 0 or 1 counts for all samples
  tmp$dds <- tmp$dds[ rowSums( counts(tmp$dds) ) > 1, ];

  # Call Deseq2 pipeline
  tmp$dds <- DESeq(tmp$dds);

  tmp$results <- results(
    object = tmp$dds,
    alpha = 0.05,
    lfcThreshold = 1
    );
  
  # Print summary of results
  summary(tmp$results);
  
  de$results[[ paste0("R", current.time) ]] <- list();
  de$results[[ paste0("R", current.time) ]] <- tmp$results;

  de$geneIDs[[ paste0("R", current.time) ]] <- rownames( subset(tmp$results, padj < 0.05 & abs(log2FoldChange) > 1) );
  };

# Count number of DE genes
# lapply( de$geneIDs, length );
# length( unique( unlist(de$geneIDs) ) );

### Get log2 fold-changes using P12 sample as baseline

# Create factor that will be used in contrast
tmp$dds <- ddsFull;
tmp$dds$Group <- factor( paste0( tmp$dds$Treat, tmp$dds$Time ) );
design( tmp$dds ) <- ~ Group;
tmp$dds <- DESeq( tmp$dds );

# List of samples to contrast
tmp$contrast.samples <- c("P15", "P17", "R12", "R12.5", "R15", "R17");

# Get IDs of genes with differential expression
tmp$de.geneIDs.list <- unique( unlist(de$geneIDs) );

de$P12.fold.change <- sapply( 
  tmp$contrast.samples,
  function(sample.name) results(
    object = tmp$dds, 
    contrast = c("Group", sample.name, "P12") )[ tmp$de.geneIDs.list, "log2FoldChange" ]
    );
rownames(de$P12.fold.change) <- tmp$de.geneIDs.list;

# Add column for P12 (fold-change = 0)
de$P12.fold.change <- cbind( data.frame(P12 = 0), de$P12.fold.change);

# Sort data.frame by var(oic.group) - var(control.group)
# Get variance from fold-change observed data
tmp$control.variance <- apply( de$P12.fold.change[, c("P12", "P15", "P17") ], 1, var );
tmp$oir.variance <- apply( de$P12.fold.change[, c("R12", "R12.5", "R15","R17") ], 1, var );  
de$P12.fold.change <- de$P12.fold.change[ order( abs( abs(tmp$oir.variance) - abs(tmp$control.variance) ), decreasing = TRUE ), ];

### Step 3 - Gather information and export result table

### Retrieve information from BioMart about DE genes (does NOT work when in IQ's wired network)

# listMarts()
# listDatasets( useMart("ensembl") )
mouse.mart <- useMart( biomart="ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl", host="www.ensembl.org" );
# listAttributes(mouse.mart) 

de$results.output <- data.frame( ensembl.ids = rownames(de$P12.fold.change) );

# Retrieve information from BioMart
de$results.output$gene.symbol             <- add.info(ensemblIDs = de$results.output$ensembl.ids, field.name = "mgi_symbol", mart = mouse.mart );
de$results.output$entrez.ids              <- add.info(ensemblIDs = de$results.output$ensembl.ids, field.name = "entrezgene", mart = mouse.mart );
de$results.output$description             <- add.info(ensemblIDs = de$results.output$ensembl.ids, field.name = "description", mart = mouse.mart );
de$results.output$human.homologues.ens    <- add.info(ensemblIDs = de$results.output$ensembl.ids, field.name = "hsapiens_homolog_ensembl_gene", mart = mouse.mart );
de$results.output$human.homologues.symbol <- add.info(ensemblIDs = de$results.output$ensembl.ids, field.name = "hsapiens_homolog_associated_gene_name", mart = mouse.mart );

# Include P12 fold change information
de$results.output <- cbind( de$results.output, round(de$P12.fold.change[de$results.output$ensembl.ids,], digits = 3) );
# de$results.output$P12.fold.change <- apply( round(de$P12.fold.change, digits = 3), 1, paste, collapse = "; ");

### Identify genes evaluated by TLDA

# Read the list of TLDA genes
tlda.genes <- read.csv(
  file = "aux_files/TLDA/Biogroup Results.csv",
  header = TRUE, 
  comment.char = "#"
  )$Target.Name;

# Get unique entries and leave only Gene Symbols
tlda.genes <- sub(
  pattern = "-(Mm|Hs).*$",
  replacement = "",
  x = unique(tlda.genes)
  );

de$results.output$on.TLDA <- ifelse( de$results.output$gene.symbol %in% tlda.genes, "X", "");

### Retrieve information about GOs
de$results.output$go.ids            <- add.info(ensemblIDs = de$results.output$ensembl.ids, field.name = "go_id", mart = mouse.mart );

# Check if any of the GOs are angiogenesis-related (following GO ids)
# blood vessel morphogenesis => GO:0048514
# blood vessel development => GO:0001568
# vasculature development => GO:0001944
# angiogenesis => GO:0001525
de$results.output$go.angio <- sapply( strsplit( de$results.output$go.ids, "; " ), function(go.list) sum( c("GO:0048514", "GO:0001568", "GO:0001944", "GO:0001525") %in% go.list ) != 0 );
de$results.output$go.angio <- ifelse( de$results.output$go.angio, "X", "");

### Pull information from PubMed

# Search on pubmed for publications with 'angiogenesis' and respective gene symbol (and human homologue) on Title or Abstract.
# Maximum expected number of PubmedIDs by gene is 100,000 (should be enough)
pubmed.results <- sapply(
  paste( de$results.output$gene.symbol, de$results.output$human.homologues.symbol, sep = "; " ),
  function(gene.symbols){
    query.symbols <- paste( paste0( strsplit(gene.symbols, "; ")[[1]], "[Title/Abstract]"), collapse = "OR+");
    query.url <- paste0( "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&retmax=100000&term=angiogenesis[Title/Abstract])AND(", query.symbols, ")" );
    xmlToList( xmlParse( getURL( query.url ) ) );
    }
  );

de$results.output$pubmed.cits <- sapply( pubmed.results, function(query.results) query.results$Count );
de$results.output$pubmed.ids <- sapply( pubmed.results, function(query.results) paste( query.results$IdList, collapse = "; " ) );

### Add information about possible drugs (drugbank)

drugbank <- list();

# Load drugbank info
for( file in list.files(path = "aux_files/drugbank_12.04.2017/", pattern = "(all|pharmacologically_active).csv", recursive = TRUE, full.names = TRUE) ){
  
  tmp$drug.activity <- sub( pattern = ".csv", replacement = "", x = basename(file) );
  tmp$drug.group <- sub( pattern = "^(.*/)", replacement = "", x = dirname(file) );
  
  drugbank[[tmp$drug.activity]][[tmp$drug.group]] <- read.csv(
    file = file
   );
  
  }

# Append drug infos to resutls data.frame
for( drug.group in names(drugbank$all) ){

  de$results.output[[paste0(drug.group,".drugs")]] <- sapply( 
    strsplit( paste( de$results.output$gene.symbol, de$results.output$human.homologues.symbol, sep = "; " ), "; "), 
    function(gene.symbol) paste0(
      paste( subset( drugbank$pharmacologically_active[[drug.group]], Gene.Name %in% gene.symbol )$Drug.IDs, collapse = "; " ),
      "  /  ",
      paste( subset( drugbank$all[[drug.group]], Gene.Name %in% gene.symbol )$Drug.IDs, collapse = "; " )
      )
    );
  }

### Export results on csv file
rownames(de$results.output) <- NULL;
write.csv( de$results.output, file = "de_results.csv" );

### Background genes (used on enrichment analysis)

# Subset genes with at least 10 reads in 2 distinct samples to use as background
de$background.genes$counts <- as.data.frame( counts(tmp$dds)[ rowSums( counts(tmp$dds) >= 10 ) >= 2, ] );

# Pull IDs information and human homology from BioMart for all expressed genes (background)
de$background.genes$ids <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name", "hsapiens_homolog_ensembl_gene", "hsapiens_homolog_associated_gene_name"),
  filters = "ensembl_gene_id",
  values = rownames(de$background.genes$counts),
  mart = mouse.mart
  );

# Pull IDs information and human homology from BioMart for DEG
de$background.genes$deg.ids <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name", "hsapiens_homolog_ensembl_gene", "hsapiens_homolog_associated_gene_name"),
  filters = "ensembl_gene_id",
  values = de$results.output$ensembl.ids,
  mart = mouse.mart
  );

# Write background files
write.xlsx(
  x = de$background.genes$counts,
  file = "../pathway_analysis/background_genes.xlsx",
  sheetName = "counts_by_gene",
  row.names = TRUE
  );

write.xlsx(
  x = de$background.genes$ids,
  file = "../pathway_analysis/background_genes.xlsx",
  sheetName = "gene_ids_and_homology",
  append = TRUE,
  row.names = TRUE
  );

write.xlsx(
  x = de$background.genes$deg.ids,
  file = "../pathway_analysis/background_genes.xlsx",
  sheetName = "DEG_gene_ids_and_homology",
  append = TRUE,
  row.names = TRUE
  );


### Step 4 - Save Rdata object for later use
save( list = c("de", "ddsFull", "drugbank", "pubmed.results", "mouse.mart", "tlda.genes"), file = "aux_files/expression_data.Rdata" );

# # Load expression data info
# load( file = "aux_files/expression_data.Rdata" );

