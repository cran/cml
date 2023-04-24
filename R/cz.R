cz <- function(w, d, dz)
{
  cz <- -w*d/dz
  cz[which(cz == -Inf)] <- 0
  cz[which(is.na(cz))] <- 0
  cz <- as.matrix(cz)
  cz_diag <- -rowSums(cz)
  cz + diag(cz_diag)
}
