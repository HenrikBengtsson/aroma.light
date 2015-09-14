setMethodS3("rowAverages", "matrix", function(X, average=base::mean, deviance=stats::sd, df=function(x, ...) length(if(na.rm) na.omit(x) else x), na.rm=TRUE, ..., asAttributes=TRUE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # 1. Verify the arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument: 'X'
  if (!is.matrix(X))
    stop(paste("Argument 'X' is not a matrix: ", class(X)[1]));

  # Argument: '...', 'average', and 'deviance'.
  args <- list(...);
  args[["average"]] <- average;
  args[["deviance"]] <- deviance;
  args[["df"]] <- df;

  for (kk in seq(along=args)) {
    key <- names(args)[kk];
    arg <- args[[kk]];
    if (!is.null(arg) && !is.function(arg))
      stop(paste("Argument '", key, "' must be a function: ", mode(arg)));
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # 2. Calculate the average and the deviance
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  stats <- list();
  for (kk in seq(along=args)) {
    key <- names(args)[kk];
    arg <- args[[kk]];
    stats[[key]] <- as.matrix(apply(X, MARGIN=1, FUN=arg, na.rm=na.rm));
  }

  if (asAttributes) {
    attrs <- attributes(X);
    X <- stats[["average"]];
    stats[["average"]] <- NULL;
    mostattributes(X) <- c(stats, attrs);
    X;
  } else {
    stats;
  }
}) # rowAverages()

############################################################################
# HISTORY:
# o 2004-05-17
#   Recreated. Made into its own function. This is need by *many* methods.
############################################################################
