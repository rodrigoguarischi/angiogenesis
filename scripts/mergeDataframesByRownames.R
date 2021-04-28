### Function - mergeDataframesByRownames() ########################################################
# Input variables: Two data frames
#
# Output variables: One merged data frame, by rownames
#
# Description: Merge two dataframes by rowname, replace its rownames and remove duplicated column
mergeDataframesByRownames <- function( dataframe.A, dataframe.B ){
  
  stopifnot( inherits( dataframe.A, "data.frame" ) );
  stopifnot( inherits( dataframe.B, "data.frame" ) );
  
  merged.dataframe <- merge(
    dataframe.A,
    dataframe.B,
    by = "row.names",
    all = TRUE
    );
  
  rownames(merged.dataframe) <- merged.dataframe$Row.names;
  merged.dataframe <- merged.dataframe[, -which(colnames(merged.dataframe) == "Row.names"), drop=FALSE];
  
  return(merged.dataframe);
  }