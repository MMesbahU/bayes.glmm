# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'

stan_models <- function(
  model = "dummy",
  rebuild = F
) {

  if(! model %in% c("glm", "glmm")) stop("Error: \"model\" must be \"glm\" or \"glmm\".")

  mydir = find.package("bayes.glmm")

  if(model == "glm") { # fixed model
    if(rebuild == F) {
      mypath = paste0(mydir, "/", "extdata/glm.rdt")
      load(mypath)
      my_stan_model = myglm
    }
    if(rebuild == T) {
      mypath = paste0(mydir, "/", "extdata/glm.stan")
      my_stan_model = stan_model(mypath)
    }
  }

  if(model == "glmm") { # mixed model
    if(rebuild == F) {
      mypath = paste0(mydir, "/", "extdata/glmm.rdt")
      load(mypath)
      my_stan_model = myglmm
    }
    if(rebuild == T) {
      mypath = paste0(mydir, "/", "extdata/glmm.stan")
      my_stan_model = stan_model(mypath)
    }
  }

  return(my_stan_model)

}

myGWAS_fit <- function(
  model = "dummy", # stan model
  method = "dummy", # inference method
  type = "dummy", # model type
  cov = NULL,  # covariate data
  geno = NULL,  # genotype data
  pheno = NULL,  # phenotype data
  prior_non = 1, # non-informative prior
  prior_inf = 0, # informative prior
  iter = 1000, # mcmc iteration
  warmup = 200, # mcmc warnup
  chains = 3, # mcmc chains
  cores = 3 # mcmc cores
) {

  if(class(model) != "stanmodel")
    stop("Error: a stan model is required.")
  if(! method %in% c("optimizing", "sampling"))
    stop("Error: \"method\" must be \"optimizing\" or \"sampling\".")
  if(! type %in% c("linear", "binary", "categorical"))
    stop("Error: \"type\" must be \"linear\", \"binary\", or \"categorical\".")
  if(any(is.na(cov)))
    stop("Error: \"cov\" is required.")
  if(any(is.na(geno)))
    stop("Error: \"geno\" is required.")
  if(any(is.na(pheno)))
    stop("Error: \"pheno\" is required.")

  N = nrow(cov)  # sample number
  D = ncol(cov)  # covariate number
  vId = rownames(geno)  # variant name
  pId = c("a", paste0("beta[", 1:D, "]"), "p") # covariate

  if(type == "linear") K = 1
  if(type != "linear") K = length(unique(pheno))

  dt0 = list(N = N, K = K, D = D, cov = cov, L = L, pheno = pheno)
  dt0$prior_non = prior_non
  dt0$prior_inf = prior_inf

  fit = list()

  for (i in 1:length(vId)) {

    dt1 = within(dt0, {geno = geno[i, ]})

    if(method == "optimizing") {
      out = optimizing(model, verbose = FALSE, hessian = TRUE, algorithm = "LBFGS", data = dt1)
      hessian = out$hessian
      hessian.idx = apply(out$hessian, 1, function(x) all(x == 0))
      hessian = hessian[! hessian.idx, ! hessian.idx]
      se1 = sqrt(diag(solve(-hessian)))["p"]
      fit[[i]] = c(out$par[pId], se1)
      names(fit[[i]]) = c(pId, "se.p")
    }

    if(method == "sampling") {
      out = sampling(model, data = dt1, chains = chains, iter = iter, warmup = warmup, cores = cores)
      fit[[i]] = summary(out, pars = pId)$summary
    }
  }

  names(fit) = vId
  return(fit)
}
