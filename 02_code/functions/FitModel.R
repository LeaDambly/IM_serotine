# modified function to take PC priors
FitModel <- function(..., formula = NULL, CovNames = NULL, mesh, spat.ind = "i",
                     predictions = FALSE, tag.pred = "pred",
                     waic = FALSE, cpo = FALSE, dic = FALSE, nthreads = NULL,
                     verbose = FALSE, compute = FALSE, control.fixed = NULL) {
  
  stck <- inla.stack(...)

  if (is.null(CovNames)) {
    CovNames <- unlist(stck$effects$names)
    CovNames <- CovNames[!CovNames %in% c(spat.ind)]
  } else{
    if (!is.null(formula)){
      warning("CovNames and formula are both not NULL")
    }
  }
  mesh <- inla.spde2.matern(Mesh$mesh)
  
  if (!is.null(spat.ind)) {
    CovNames <- c(CovNames, paste0("f(", spat.ind, ", model=mesh)"))
  }
  if (is.null(control.fixed)) {
    control.fixed <- list(mean = 0)
  }
  if (is.null(formula)) {
    Formula <- formula(paste(c("resp ~ 0 ", CovNames), collapse = "+"))
  } else {
    if (is.null(spat.ind)) {
      Formula <- formula
    } else {
      if (any(grepl(paste0("(", spat.ind, ","), formula,
                    fixed = TRUE
      ))) {
        warning(paste0(spat.ind, " already in formula, so will be ignored"))
        Formula <- formula
      } else {
        Formula <- update(formula, paste0(
          " ~ . + f(", spat.ind,
          ", model=mesh)"
        ))
      }
    }
  }

  mod <- inla(Formula,
    family = c("poisson", "binomial"), # likelihood family
    E = inla.stack.data(stck)$e, # component in the mean for the Poisson likelihoods
    Ntrials = inla.stack.data(stck)$Ntrials, # number of trials for the binomial likelihood
    verbose = verbose, # verbose mode for debugging
    control.compute = list(waic = waic, cpo = cpo, dic = dic),
    # control.fixed = control.fixed,
    control.predictor = list(
      A = inla.stack.A(stck),
      compute = compute, link = NULL
    ),
    control.family = list(
      list(control.link = list(model = "log")),
      list(control.link = list(model = "cloglog"))
    ),
    data = inla.stack.data(stck)
  )


  if (predictions) {
    id <- inla.stack.index(stck, tag.pred)$data
    pred <- data.frame(
      mean = mod$summary.fitted.values$mean[id],
      stddev = mod$summary.fitted.values$sd[id]
    )
    res <- list(model = mod, predictions = pred, form = Formula)
  }
  else {
    res <- list(model = mod, form = Formula)
  }
  res
}
