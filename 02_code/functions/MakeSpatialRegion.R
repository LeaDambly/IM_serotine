MakeSpatialRegion <- function(data = NULL, coords = c("X", "Y"), meshpars,
                              spdepars = NULL,
                              bdry = NULL,
                              proj = CRS("+proj=utm")) {
  require(rgeos)

  # extract the boundary of SpatialPolygons object
  region.bdry <- inla.sp2segment(bdry)

  # create the mesh using the domain of a region
  mesh <- inla.mesh.2d(
    boundary = region.bdry, cutoff = meshpars$cutoff,
    max.edge = meshpars$max.edge, offset = meshpars$offset
  )

  # creates a Matern SPDE model
  if (is.null(spdepars)) {
    spde <- inla.spde2.matern(mesh = mesh, alpha = 2)
  } else {
    spde <- inla.spde2.pcmatern(
      mesh = mesh,
      alpha = spdepars$alpha,
      prior.range = spdepars$range,
      prior.sigma = spdepars$sigma
    )
  }

  # compute weights
  # make dual mesh
  dd <- deldir::deldir(mesh$loc[, 1], mesh$loc[, 2])
  tiles <- deldir::tile.list(dd)

  # intersection between domain and dual mesh
  poly.gpc <- as(bdry@polygons[[1]]@Polygons[[1]]@coords, "gpc.poly")

  # w now contains area of voronoi polygons
  # The weights associated to the mesh points are equal to the area of the surrounding polygon
  # in a Voronoi tessellation that falls inside the study region
  w <- sapply(tiles, function(p) rgeos::area.poly(rgeos::intersect(as(cbind(p$x, p$y), "gpc.poly"), poly.gpc)))

  return(list(mesh = mesh, spde = spde, w = w))
}
