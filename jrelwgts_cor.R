#' Johnson Relative Weights - Correlation matrix as input
#'
#' Calculates relative importance weights directly from a correlation matrix.
#' Ideal for polychoric/tetrachoric matrices or pairwise-deletion matrices.
#'
#' @param cor_mat A square, symmetric correlation matrix or data frame. 
#'   The first column/row MUST be the dependent variable (Y).
#'
#' @return A tibble containing the original correlations, raw relative weights, 
#'   and rescaled importance percentages for each predictor.
#' @export
jrelwgts_cor <- function(cor_mat) {
  # check inputs
  cor_mat <- as.matrix(cor_mat)
  if (!isSymmetric(cor_mat, tol = 1e-8)) {
    stop("Input 'cor_mat' must be a symmetric correlation matrix.")
  }
  if (is.null(colnames(cor_mat))) {
    colnames(cor_mat) <- paste0("X", seq_len(ncol(cor_mat)) - 1)
    colnames(cor_mat)[1] <- "Y"
  }
  
  nvar <- ncol(cor_mat)
  
  # parse submatrices
  # rxx = correlation among predictors; rxy = correlation between Y and predictors
  rxx <- cor_mat[2:nvar, 2:nvar, drop = FALSE]
  rxy <- cor_mat[2:nvar, 1, drop = FALSE]
  
  # spectral decomposition (eigenvalues & eigenvectors)
  eigen_decomp <- eigen(rxx)
  evec <- eigen_decomp$vectors
  ev   <- eigen_decomp$values
  
  # ensure eigenvalues aren't slightly negative due to floating-point precision
  ev[ev < 0] <- 0 
  
  # calc delta accouting for possibility of small dimensions
  delta  <- diag(sqrt(ev), nrow = length(ev))
  lambda <- evec %*% delta %*% t(evec)
  
  # relative weights
  lambdasq <- lambda ^ 2
  beta     <- solve(lambda) %*% rxy
  rsquare  <- sum(beta ^ 2)
  
  rawwgt   <- lambdasq %*% (beta ^ 2)
  import   <- (rawwgt / rsquare) * 100
  
  # prep for output
  # drop() or as.vector() ensures matrix structures don't break the tibble
  output_df <- tibble::tibble(
    item        = colnames(cor_mat)[-1],
    correlation = as.vector(rxy),
    raw_weight  = as.vector(rawwgt),
    importance  = as.vector(import)
  )
  
  message(sprintf("Model Fit (R-squared): %0.4f", rsquare))
  
  return(output_df)
}
