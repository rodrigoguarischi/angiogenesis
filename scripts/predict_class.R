library( randomForestSRC );

### Function - predict.class() ####################################################################
# Input variables: A rfsrc object and a vector with time cut-points
#
# Output variables: A factor with predicted classes
#
# Description: Assign a group ("low", "intermediate", "high") to each sample of one rfsrc object 
predict.class <- function( rfsrc.obj, breaks ) {
  
  # Fail if object does not inherits rfSRC class
  stopifnot( inherits( rfsrc.obj, "rfsrc" ) );
  
  # You need to set two, and only two, time breaks or breaks does not make sense
  stopifnot( length(breaks) == 2 & breaks[1] < breaks[2] );

  predicted.classes <- ifelse(
    rfsrc.obj$survival[, sum( rfsrc.obj$time.interest < breaks[2] ) ] > 0.5,
    "low",
    ifelse( 
      rfsrc.obj$survival[, sum( rfsrc.obj$time.interest < breaks[1] ) ] > 0.5,
      "intermediate",
      "high"
      )
    );

  return( factor( predicted.classes ) );
  }