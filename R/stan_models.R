# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'

stan_models <- function(type) {

  if(! type %in% c("glm", "glmm"))
    stop("Error: \"type\" must be \"glm\" or \"glmm\".")

  if(type == "glm")
    model = stan_model("./lib/noK.stan")

  if(type == "glmm")
    model = stan_model("./lib/glmm.stan")

  return(model)

}
