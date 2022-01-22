' shrink_model.R
Shrink linear model from MedReMix.
Based on https://win-vector.com/2014/05/30/trimming-the-fat-from-glm-models-in-r/.

Usage:
    shrink_model.R -i INPUT -o OUTPUT

Options:
    -i --input INPUT            Path to input Rds file
    -o --output OUTPUT          Output path (RDS file)
' -> doc

if (! interactive()) {
    library(docopt)
    args <- docopt(doc, version='Shrink Model v1')
    print(args)
} else {
    message('Running in interactive mode. Be sure to specify args manually.')
}

medremix_model <- readRDS(args[['input']])

stripGlmLR = function(cm) {
  cm$y = c()
  cm$model = c()
  
  cm$residuals = c()
  cm$fitted.values = c()
  cm$effects = c()
  cm$qr$qr = c()  
  cm$linear.predictors = c()
  cm$weights = c()
  cm$prior.weights = c()
  cm$data = c()

  
  cm$family$variance = c()
  cm$family$dev.resids = c()
  cm$family$aic = c()
  cm$family$validmu = c()
  cm$family$simulate = c()
  attr(cm$terms,".Environment") = c()
  attr(cm$formula,".Environment") = c()
  
  cm
}

medremix_model$final_model <- stripGlmLR(medremix_model$final_model)
medremix_model$iteration_models <- lapply(medremix_model$iteration_models, function(z) {stripGlmLR(z)})
medremix_model$zero_model$theta_fit <- stripGlmLR(medremix_model$zero_model$theta_fit)
medremix_model$zero_model$mu_fit <- stripGlmLR(medremix_model$zero_model$mu_fit)

saveRDS(medremix_model, args[['output']])
