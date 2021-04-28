# Define working directory
setwd(dir = "~/angiogenesis_transcriptome/pathway_analysis/post_analysis/cancer_genes_lists/");

# Disable automatic conversion of strings to factors
options(stringsAsFactors = FALSE);

# Load required libraries and datasets
library( VennDiagram );

listofgenes <- list();

listofgenes$cbio <- names( read.delim(file = "cBio_cancer_genes.gmt", sep = "\t")[-c(1,2)] );
listofgenes$cosmic <- read.delim(file = "list_cosmic.txt", sep = "\n", header = FALSE)[,1];

venn.plot <- venn.diagram(
  
  x = list(
    cbio = listofgenes$cbio,
    cosmic = listofgenes$cosmic
    ),
  filename = NULL,
  
  ext.percent = 0,
  lty = "blank",
  fill = c("dodgerblue","red"),
  alpha = 0.35,
  label.col = c("darkblue", "white", "darkred"),
  cex = 2,
  fontface = "bold",
  
  cat.col = c("darkblue", "darkred"),
  cat.cex = 2,
  cat.pos = c(0,0),
  cat.default.pos = "text",
  cat.dist = 0.05
  );

# Remove log file created
file.remove( list.files( pattern="Venn.*.log$") )

#
# Perform Hypergeometric test to access gene overlap significancy
# Assuming number of genes equals to 20k genes
# 

phyper(
  q = sum(listofgenes$cbio %in% listofgenes$cosmic), # Overlap between 2 sets
  m = length(listofgenes$cbio),                      # Number of angiogenesis genes
  n = 20000 - length(listofgenes$cbio),              # Total number of genes - number of angiogenesis genes
  k = length(listofgenes$cosmic),                    # Number of hypoxia genes
  lower.tail = FALSE                                 # What I want is P[X > x].
  );

#
# p = 0
#

# Export image to tiff
tiff(
  filename = "venn_genes.tiff", 
  type = "cairo", 
  width = 15, 
  height = 15, 
  units = 'cm', 
  res = 500, 
  compression = 'lzw'
);

# # Display plot
# grid.newpage();
grid.draw(venn.plot);

dev.off();

