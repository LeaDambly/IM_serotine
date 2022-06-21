# script to run multiple different mesh models in loop

# ignore proj4 warnings
rgdal::set_rgdal_show_exportToProj4_warnings(FALSE)

# model run----
library(INLA)

load("01_data/data_serotine.RData")

files.sources <- paste("02_code/functions/", list.files("02_code/functions/"), sep = "")
sapply(files.sources, source)

# set mesh parameter values
x <- 5 # how many meshes

# max edge
me <- list(c(5, 150),
           c(25, 150),
           c(45, 150),
           c(65, 150),
           c(85, 150))

# offset
of <- rep(list(c(10, 300)),x)

# cutoff
cf <- list(1,
           5,
           9,
           13,
           17)

# spde priors?
sppr <- T

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

for (i in 1:length(me)) {
  # get mesh parameters
  Meshpars <- list(
    max.edge = me[[i]],
    cutoff = cf[[i]],
    offset = of[[i]]
  )

  # mesh parameters if priors are involved
  if (sppr == T) {
    Spdepars <- list(
      alpha = alpha[[i]],
      range = rg[[i]],
      sigma = sg[[i]]
    )

  } else {
    Spdepars <- NULL
    
  }
  # make mesh
  Mesh <- MakeSpatialRegion(
    data = NULL,
    bdry = outline,
    meshpars = Meshpars,
    spdepars = Spdepars,
    proj = proj
  )
  
  # integration points
  stk.ip <- MakeIntegrationStack(
    mesh = Mesh$mesh,
    data = covariates,
    area = Mesh$w,
    tag = "ip",
    InclCoords = F
  )

  nxy.scale <- 5 # scale of projection grid
  nxy.size <- c(diff(outline@bbox[1, ]), diff(outline@bbox[2, ]))
  nxy <- round(nxy.size / nxy.scale)
  gc()

  # projection grid
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

  # NBN = presence only
  stk.PO <- MakePointsStack(
    presences = obs$PO,
    data = covariates,
    mesh = Mesh$mesh,
    tag = "PO",
    InclCoords = F
  )

  gc()

  # model formula
  form <- formula(resp ~ 0 + broadleaf + arable +
    grassland + temperature + Intercept +
    int.PA + int.PO + 
    f(i, model = Mesh$spde))
  
  # fit model
  model <- FitModel(
    stk.ip, stk.PA, stk.PO, stk.pred$stk,
    formula = form,
    CovNames = NULL,
    mesh = Mesh$mesh,
    compute = T,
    predictions = T, # transforming the linear predictors by the inverse of the link function
    waic = T,
    dic  = T,
    cpo = T,
    verbose = F
  )
  
  # recalc cpo if any failed
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


save(
  stcks, models, info, 
  file = "03_output/serotine1.RData"
)
