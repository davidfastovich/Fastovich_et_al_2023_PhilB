get_confint <- function(x) {
  coefs <- as.data.frame(summary(x)$beta_table)
  lower <- coefs[,'Estimate'] - 1.96 * coefs[, 'Cond. SE']
  upper <- coefs[,'Estimate'] + 1.96 * coefs[, 'Cond. SE']
  coefs$lower <- lower
  coefs$upper <- upper
  coefs$sig <- coefs[,"Estimate"] < coefs[, "lower"] | coefs[,"Estimate"] > coefs[, "upper"]
  return(coefs)
}
