myGWAS_sampling <- function(model, geno, data0, config0) {
  y = list()

  for (i in 1:nrow(geno)) {
    data1 = within(data0, {g = geno[i, ]})

    fit = sampling(model,
                   data = dat1,
                   chains = config0$chains,
                   iter = config0$iter,
                   warmup = config0$warmup,
                   cores = config0$cores)

    y[[i]] = summary(fit, pars = c("p", "beta"))$summary
  }

  names(y) = rownames(geno)

  return(y)
}
