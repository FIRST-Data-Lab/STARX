comb <- function(...) {
  mapply("rbind", ..., SIMPLIFY = FALSE)
}
