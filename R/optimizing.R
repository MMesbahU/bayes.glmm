
myGWAS_optimizing <- function(model, geno, data0) {

  N1 = nrow(geno)
  vId = rownames(geno)
  pId = c("p", paste0("beta[", 1:2, "]"), paste0("c[", 1:3, "]"))

  y = matrix(nrow = N1, ncol = 8, dimnames = list(vId, c(pId, "se", "lp")))

  for (i in 1:N1) {
    data1 = within(data0, {geno = geno[i, ]})
    fit1 = optimizing(model, verbose = FALSE, hessian = TRUE, algorithm = "LBFGS", data = data1)

    se1 = sqrt(diag(solve(-fit1$hessian)))["p"]
    lp1 = sum(fit1$par[paste0("lp[", 1:576, "]")])

    y[i, ] = c(fit1$par[pId], se1, lp1)
  }

  return(y)
}
