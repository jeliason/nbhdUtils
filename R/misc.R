#' Check a BRMS model using the DHARMa simulated residuals
#'
#' @param model brms model
#' @param integer integer response? (TRUE/FALSE)
#' @param plot make plot?
#' @param ... further arguments for DHARMa::plotResiduals 
#'
#' @return a DHARMa object
#' @export
#'
#' @examples
check_brms <- function(model,             
                       integer = TRUE,   
                       plot = TRUE,       
                       ...                
) {
  
  mdata <- brms::standata(model)
  if (!"Y" %in% names(mdata))
    stop("Cannot extract the required information from this brms model")
  
  dharma.obj <- DHARMa::createDHARMa(
    simulatedResponse = t(brms::posterior_predict(model, ndraws = 1000, re.form = NULL)),
    observedResponse = mdata$Y, 
    fittedPredictedResponse = apply(
      t(brms::posterior_epred(model, ndraws = 1000, re.form = NULL)),
      1,
      mean),
    integerResponse = integer)
  
  if (isTRUE(plot)) {
    plot(dharma.obj, ...)
  }
  
  invisible(dharma.obj)
}