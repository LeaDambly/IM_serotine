#' Function to create stack for predictions
#' @param nxy Number of points in x and y directions.
#' @param mesh INLA mesh.
#' @param data Data frame with columns for coordinates, and others are covariates.
#' @param tag Name for tag for the stack (defaults to "points").
#' @param coordnames Names of coorinates (defaults to X and Y)
#' @param boundary Boundary of region to project onto. Defaults to NULL, when the boundary of the mesh will be used. Either of class SpatialPolygons or two columns with the coorinates of the polygon
#' @param intercept Logical: should an intercept be added? Defaults to TRUE
#'
#' @return An INLA stack onto which new data can be projected
#'
#' @export
#' @import INLA

MakeProjectionGrid <- function(nxy, mesh, boundary, data, tag = "pred", coordnames = c("X", "Y"), intercept = TRUE) {
  if ("resp" %in% coordnames) stop("resp cannot be a coordinate name")
  if ("e" %in% coordnames) stop("e cannot be a coordinate name")
  if (is.null(boundary)) stop("you need a boundary")
  if (class(boundary) != "SpatialPolygons") stop("boundary needs to be of class SpatialPolygons")

  projgrid <- inla.mesh.projector(mesh,
    xlim = boundary@bbox["x", ],
    ylim = boundary@bbox["y", ],
    dims = nxy
  )
  # get the points on the grid within the boundary
  xy.in <- splancs::inout(
    projgrid$lattice$loc,
    boundary@polygons[[1]]@Polygons[[1]]@coords
  )

  coord.prd <- projgrid$lattice$loc[xy.in, ]
  colnames(coord.prd) <- coordnames
  A.prd <- projgrid$proj$A[xy.in, ]

  # Extract covariates for points, add intercept and coordinates
  NearestCovs <- GetNearestCovariate(points = coord.prd, covs = data)
  if (intercept) NearestCovs$Intercept <- NA
  NearestCovs@data[, colnames(NearestCovs@coords)] <- NearestCovs@coords

  # stack the predicted data
  stk <- inla.stack(list(
    resp = cbind(NA, rep(NA, nrow(NearestCovs))),
    e = rep(0, nrow(NearestCovs))
  ),
  A = list(1, A.prd), tag = tag,
  effects = list(NearestCovs@data, list(i = 1:mesh$n))
  )

  pred <- list(stk = stk, xy.in = xy.in, predcoords = coord.prd)
  pred
}
