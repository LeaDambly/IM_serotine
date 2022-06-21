#' Function to get covariate values nearest the data/integration point
#'
#' @param points Points for which we need covariates. Either a SpatialPoints* object or a matrix with coordinates in columns.
#' @param covs Covariates, with coordinates in first and second columns, or a SpatialPointsDataFrame object.
#' @return A SpatialPointsDataFrame with the original points plus the data from the closest point.
#'
#' If there are ties, this will use the first element.
#'
#' @export
GetNearestCovariate <- function(points, covs) {
  require(RANN)
  require(dplyr)
  if(class(points)=="SpatialPointsDataFrame") {
    points <- points@coords
  }
  
  covnames <- names(covs@data)
  pointsXY <- as.data.frame(points)
  covXY <- cbind(covs@coords, covs@data)
  covXY$ID <- seq(nrow(covXY))
  closest <- nn2(covXY[c("X", "Y")], pointsXY[c("X", "Y")], k = 1)
  # we have closest which contains all pointsXY and their nearest covXY ID
  # need to find covariate values for each ID
  # then add those to each point
  closest <- as.data.frame(closest) %>%
    rename(ID = nn.idx)
  joined <- inner_join(closest, covXY, by = "ID") %>%
    dplyr::select(-c(nn.dists, X, Y, ID))
  allpts <- bind_cols(pointsXY, joined)
  
  points.df <- SpatialPointsDataFrame(coords = allpts[c("X", "Y")],
                                      data = allpts[covnames],
                                      proj4string = CRS(proj4string(covs)))
  return(points.df)
}

