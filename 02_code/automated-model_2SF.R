# NOTE: there might be a lot of warnings regarding spatial objects
# this is because of R's move from proj4 to proj6
# INLA however does NOT yet support proj6, so no point in changing everything - just ignore the warnings
rgdal::set_rgdal_show_exportToProj4_warnings(FALSE)

# model run----
library(INLA)

load("01_data/data_serotine.RData")
#load("01_data/data_serotine_wgs.RData")

files.sources <- paste("02_code/functions_2SF/", list.files("02_code/functions_2SF/"), sep = "")
sapply(files.sources, source)

x <- 5
me <- list(c(5, 150),
           c(25, 150),
           c(45, 150),
           c(65, 150),
           c(85, 150))

of <- rep(list(c(10, 300)),x)

cf <- list(1,
           5,
           9,
           13,
           17)

# spde priors?
sppr <- T

# my attempt
# Probability that the range of the spatial effect is below 7 is 0.01
# Probability that the variation in the spatial effect is larger than 1 is 0.5
if (sppr == T) {
  alpha <- rep(list(2), x)
  rg <- rep(list(c(7, 0.01)), x)
  sg <- rep(list(c(1, 0.5)), x)
}

stcks <- list()
models <- list()
info <- list()

for (i in 2:length(me)) {
  Meshpars <- list(
    max.edge = me[[i]],
    cutoff = cf[[i]],
    offset = of[[i]]
  )
  
  if (sppr == T) {
    Spdepars <- list(
      alpha = alpha[[i]],
      range = rg[[i]],
      sigma = sg[[i]]
    )
    
  } else {
    Spdepars <- NULL
    
  }
  
  Mesh <- MakeSpatialRegion(
    data = NULL,
    bdry = outline,
    meshpars = Meshpars,
    spdepars = Spdepars,
    proj = proj
  )
  
  stk.ip <- MakeIntegrationStack(
    mesh = Mesh$mesh,
    data = covariates,
    area = Mesh$w,
    tag = "ip",
    InclCoords = F
  )
  
  nxy.scale <- 5
  nxy.size <- c(diff(outline@bbox[1, ]), diff(outline@bbox[2, ]))
  nxy <- round(nxy.size / nxy.scale)
  gc()
  
  stk.pred <- MakeProjectionGrid(
    nxy = nxy,
    mesh = Mesh$mesh,
    boundary = outline,
    data = covariates,
    tag = "pred",
    coordnames = c("X", "Y")
  )
  
  
  # Create data stacks for each of the different data types
  # field = presence-absence
  stk.PA <- MakeBinomStack(
    observs = obs$PA,
    data = covariates,
    mesh = Mesh$mesh,
    presname = "NPres",
    trialname = "Ntrials",
    tag = "PA",
    InclCoords = F
  )
  
  # roost = presence only
  stk.PO <- MakePointsStack(
    presences = obs$PO,
    data = covariates,
    mesh = Mesh$mesh,
    tag = "PO",
    InclCoords = F
  )
  
  gc()
  
  form <- formula(resp ~ 0 + broadleaf + arable +
                    grassland + temperature + Intercept +
                    int.PA + int.PO + 
                    f(shared_field, model = Mesh$spde) + 
                    f(i, copy = 'shared_field', fixed = TRUE) + 
                    f(bias_field, model = Mesh$spde))
                  
  # cor <- parallel::detectCores()
  
  stck <- inla.stack(stk.ip, stk.PA, stk.PO, stk.pred$stk)
  
  mod <- inla(formula = form,
              family = c("poisson", "binomial"),
              E = inla.stack.data(stck)$e, # component in the mean for the Poisson likelihoods
              Ntrials = inla.stack.data(stck)$Ntrials, # number of trials for the binomial likelihood
              verbose = FALSE, # verbose mode for debugging
              control.compute = list(waic = TRUE, cpo = TRUE, dic = FALSE),
              # control.fixed = control.fixed,
              control.predictor = list(
                A = inla.stack.A(stck),
                compute = TRUE, link = NULL
              ),
              control.family = list(
                list(control.link = list(model = "log")),
                list(control.link = list(model = "cloglog"))
              ),
              data = inla.stack.data(stck)
  )
  
  id <- inla.stack.index(stck, "pred")$data
  pred <- data.frame(
    mean = mod$summary.fitted.values$mean[id],
    stddev = mod$summary.fitted.values$sd[id]
  )
  model <- list(model = mod, predictions = pred, form = form)
  
  if(sum(model$cpo$failure > 0, na.rm = T) > 0){
    improved.result = inla.cpo(model)
  } else {
    improved.result = "NULL"
  }
  
  if(is.null(Spdepars)){
    Spdepars <- list(
      alpha = "NULL",
      range = c("NULL","NULL"),
      sigma = c("NULL","NULL"))
  }
  
  stcks[[i]] <- list(Mesh, stk.ip, stk.pred, stk.PA, stk.PO)
  models[[i]] <- list(model, improved.result)
  info[[i]] <- list(Meshpars, Spdepars)
  gc()
  
}

#beepr::beep("ping")

save(
  stcks, models, info, 
  file = "03_output/serotine_2SF.RData"
)


