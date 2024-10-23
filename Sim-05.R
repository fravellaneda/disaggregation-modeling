#Simulations

library(spatstat)
library(ggplot2)
library(viridis)
library(mvtnorm)
library(dplyr)
library(sf)
library(terra)
require(gridExtra)
library(INLA)
library(SpatialEpi)

set.seed(1234)

roi <- square(1) |>
  as.polygonal() |>
  as.data.frame() |>
  st_as_sf(coords = c("x","y")) |>
  st_combine() |>
  st_cast("POLYGON") |>
  as_Spatial()
boundary <- as(roi,"sf")

n_areas <- 200
fun_gen <- function(x,y) {n_areas*(3/2)*((x)**2+(y-0)**2)}

book.rMatern <- function(n, coords, sigma=1, range, 
                         kappa = sqrt(8*nu)/range, variance = sigma^2, nu=1) {
  m <- as.matrix(dist(coords))
  m <- exp((1-nu)*log(2) + nu*log(kappa*m)-
             lgamma(nu))*besselK(m*kappa, nu)
  diag(m) <- 1
  return(drop(crossprod(chol(variance*m),
                        matrix(rnorm(nrow(coords)*n), ncol=n))))
}

between <- function(x){
  min <- x[1]
  max <- x[2]
  value <- x[3]
  if((value<max) & (value>min)){
    return(1)
  }
  else{
    return(0)
  }
}

alpha <- c(0.1, 0.05, -0.1) 
m.var <- c(0.25,0.05,0.2, 0.15)
range <- c(0.1,0.3,0.2,0.1)

beta <- c(0.3, 0.1, -0.3) 
e.sd <- c(0.1, 0.2, 0.15)

pre1 <- 0.2
pre2 <- 0.1
pre3 <- 0.05  

resolution1 <- 50
x <- seq(1/(2*resolution1), 1-1/(2*resolution1), length.out = resolution1)
dp_pre <- as.matrix(expand.grid(x, x))

theta.ini <- (c(log(range), log((1/m.var)**2), log(1 / e.sd^2),
                beta))[c(1,5,2,6,3,7,4,8,9,10,11,12,13,14)]
theta.ini.e = theta.ini + rnorm(length(theta.ini), 0, 0.1)

varF <- function(m) {
  inla.tmarginal(function(x) 1 / sqrt(exp(x)), m)
}

for (seed in 1:100){
  print(seed)
  #set.seed(seed)
  
  p_gen <- rpoispp(fun_gen,win=square(1))
  p_gen_sf <- st_as_sf(data.frame(x=p_gen$x,y=p_gen$y),coords = c("x","y"))
  voronoi <- terra::voronoi(x = terra::vect(p_gen_sf),bnd = boundary)
  voronoi$id <- 1:nrow(voronoi)
  
  #set.seed(seed)
  z1 <- book.rMatern(1, dp_pre, range = range[1],
                     sigma = sqrt(m.var[1]))
  z2 <- book.rMatern(1, dp_pre, range = range[2],
                     sigma = sqrt(m.var[2]))
  z3 <- book.rMatern(1, dp_pre, range = range[3],
                     sigma = sqrt(m.var[3]))
  z4 <- book.rMatern(1, dp_pre, range = range[4],
                     sigma = sqrt(m.var[4]))
  
  z1 <- z1 - mean(z1)
  z2 <- z2 - mean(z2)
  z3 <- z3 - mean(z3)
  z4 <- z4 - mean(z4)
  
  linear1 <- alpha[1] + z1 + z2 +rnorm(nrow(dp_pre), 0, e.sd[1])
  linear2 <- alpha[2] + beta[1] * z1 + z3 + rnorm(nrow(dp_pre), 0, e.sd[2])
  linear3 <- alpha[3] + beta[2] * z1 + beta[3] * z3 + z4 + 
    rnorm(nrow(dp_pre), 0, e.sd[3])
  
  theta1 <- exp(linear1)
  theta2 <- exp(linear2)
  theta3 <- exp(linear3)
  
  N <- 100000
  dpdf <- as.data.frame(dp_pre)
  names(dpdf) <- c("x","y")
  
  dpdf$linear1 <- linear1
  dpdf$linear2 <- linear2
  dpdf$linear3 <- linear3
  
  dpdf$theta1 <- theta1
  dpdf$theta2 <- theta2
  dpdf$theta3 <- theta3
  
  dpdf$pop <- ((1/resolution1)**2)*N
  
  dpdf$e1 <- dpdf$pop * pre1
  dpdf$e2 <- dpdf$pop * pre2
  dpdf$e3 <- dpdf$pop * pre3
  
  dpdf$lambda1 <- theta1*dpdf$e1
  dpdf$lambda2 <- theta2*dpdf$e2
  dpdf$lambda3 <- theta3*dpdf$e3
  
  dpdf$sir1 <- dpdf$lambda1/dpdf$e1
  dpdf$sir2 <- dpdf$lambda2/dpdf$e2
  dpdf$sir3 <- dpdf$lambda3/dpdf$e3
  
  ## paso simulación
  #set.seed(seed)
  im.1 <- as.im(list(x=x,y=x,z=matrix(dpdf$lambda1 * resolution1**2,
                                      nrow=resolution1)))
  cases.y1 <- st_as_sf(coords(rpoispp(im.1)),coords = c("x","y"))
  cases.y1$y1 <- 1
  cases.y1$y2 <- 0
  cases.y1$y3 <- 0
  
  im.2 <- as.im(list(x=x,y=x,z=matrix(dpdf$lambda2 * resolution1**2,
                                      nrow=resolution1)))
  cases.y2 <- st_as_sf(coords(rpoispp(im.2)),coords = c("x","y"))
  cases.y2$y1 <- 0
  cases.y2$y2 <- 1
  cases.y2$y3 <- 0
  
  im.y3 <- as.im(list(x=x,y=x,z=matrix(dpdf$lambda3 * resolution1**2,
                                       nrow=resolution1)))
  cases.y3 <- st_as_sf(coords(rpoispp(im.y3)),coords = c("x","y"))
  cases.y3$y1 <- 0
  cases.y3$y2 <- 0
  cases.y3$y3 <- 1
  
  cases <- rbind(cases.y1,cases.y2,cases.y3)
  
  voronoi_sf <- st_as_sf(voronoi)
  voronoi_sf$pop <- round(st_area(voronoi_sf)*N,0)
  voronoi_int <- st_intersection(voronoi_sf,cases)
  values <- as.data.frame(voronoi_int)
  values_agg <- aggregate(values[,c("y1","y2","y3")], by=list(values$id), sum)
  values_agg <- values_agg %>% mutate(across(c("y1","y2","y3"), round, 0))
  names(values_agg)[1] <- "id"
  map <- merge(voronoi_sf,values_agg,on="id",all.x=TRUE)
  map[is.na(map)] <- 0
  map$E1 <- expected(population = map$pop,
                     cases = map$y1, n.strata = 1)
  map$E2 <- expected(population = map$pop,
                     cases = map$y2, n.strata = 1)
  map$E3 <- expected(population = map$pop,
                     cases = map$y3, n.strata = 1)
  #SIR
  map$sir.1 <- map$y1/map$E1
  map$sir.2 <- map$y2/map$E2
  map$sir.3 <- map$y3/map$E3
  
  #Points inside each area
  for (i in 1:nrow(map)){
    sample_temp <- st_sample(map[i,],10,type="hexagonal")
    if (i==1){
      loc.d <- cbind(st_coordinates(sample_temp)[, 1],
                     st_coordinates(sample_temp)[, 2])
    }else{
      loc.d <- rbind(loc.d,
                     cbind(st_coordinates(sample_temp)[, 1],
                           st_coordinates(sample_temp)[, 2]))
    }
  }
  
  resolution2 <- 50
  x <- seq(1/(2*resolution2), 1-1/(2*resolution2), length.out = resolution2)
  dp <- as.matrix(expand.grid(x, x))
  bnd <- inla.nonconvex.hull(dp)
  
  mesh <- inla.mesh.2d(
    loc = loc.d, boundary = bnd,
    max.edge = c(0.05,0.1)
  )
  
  mesh_df <- as.data.frame(mesh$loc[,1:2])
  names(mesh_df) <- c("x","y")
  
  spde <- inla.spde2.matern(mesh = mesh, alpha = 2, constr = TRUE)
  
  # Number of mesh nodes
  nv <- mesh$n
  
  # number of groups
  n.pp <- 3
  
  n.ar <- nrow(map)
  p <- map$pop
  
  y.1 <- matrix(NA, nrow = n.ar, ncol = n.pp)
  y.1[, 1] <- map$y1
  e.1 <- map$E1
  
  y.2 <- matrix(NA, nrow = n.ar, ncol = n.pp)
  y.2[, 2] <- map$y2
  e.2 <- map$E2
  
  y.3 <- matrix(NA, nrow = n.ar, ncol = n.pp)
  y.3[, 3] <- map$y3
  e.3 <- map$E3
  
  mesh_sf <- st_as_sf(data.frame(x=mesh$loc[,1],y=mesh$loc[,2],id=1:nv),
                      coords = c("x","y"))
  
  # Projection matrix
  
  A.m <- matrix(0, n.ar, nv)
  
  t1 <-Sys.time()
  for (i in 1:n.ar){
    inter_sf <- st_intersection(mesh_sf,map[i,c()])
    vor <- terra::voronoi(x = terra::vect(inter_sf),bnd = map[i,])
    vor_sf <- st_as_sf(vor)
    vor_int <- st_intersection(vor_sf,map[i,c()])
    vor_int$area <- st_area(vor_int)
    vor_int$prop <- st_area(vor_int)/st_area(map[i,c()])
    vor_int <- vor_int %>% arrange(id)
    idx <- vor_int$id
    A.m[i,idx] <- vor_int$prop
    mesh_sf[mesh_sf$id %in% idx,"area"] <- vor_int$area
  }
  t2 <-Sys.time()
  print("Time matrix construction ")
  print(t2-t1)
  
  mesh_sf[,"pop"] <- N*mesh_sf$area
  mesh_sf[,"pre1"] <- N*pre1*mesh_sf$area
  mesh_sf[,"pre2"] <- N*pre2*mesh_sf$area
  mesh_sf[,"pre3"] <- N*pre3*mesh_sf$area
  
  # Stack for estimatio
  # Create the stack for y1
  stk.y1 <- inla.stack(
    data = list(y = y.1, e = e.1),
    A = list(A.m,A.m,A.m,A.m),
    effects = list(
      spatial.field.sh = 1:nv,
      spatial.field.y1 = 1:nv,
      iid.y1 = 1:nv,
      Intercept.y1 = rep(1, nv)
      #pre = 1/mesh_sf$pre1
    ),
    tag = "y1")
  
  # Create the stack for y2
  stk.y2 <- inla.stack(
    data = list(y = y.2, e = e.2), 
    A = list(A.m,A.m,A.m,A.m), 
    effects = list(
      sh.copy.y2 = 1:nv,
      spatial.field.y2 = 1:nv,
      iid.y2 = 1:nv,
      Intercept.y2 = rep(1, nv)
      #pre = 1/mesh_sf$pre2
    ),
    tag = "y2")
  
  # Create the stack for y3
  stk.y3 <- inla.stack(
    data = list(y = y.3, e = e.3), 
    A = list(A.m,A.m,A.m,A.m,A.m), 
    effects = list(
      sh.copy.y3 = 1:nv,
      y2.copy.y3 = 1:nv,
      spatial.field.y3 = 1:nv,
      iid.y3 = 1:nv,
      Intercept.y3 = rep(1, nv)
      #pre = 1/mesh_sf$pre3
    ),
    tag = "y3")
  
  A.pr <- inla.spde.make.A(mesh = mesh, loc = dp)
  
  # Stack for predicting y1
  stk.y1.pre <- inla.stack(tag = "y1.pre",
                           data = list(y = matrix(NA, nrow=nrow(dp),ncol=n.pp),
                                       e = rep(0,nrow(dp))),
                           A = list(A.pr, A.pr, A.pr, A.pr),
                           effects = list(
                             Intercept.y1 = rep(1, nv),
                             spatial.field.sh = 1:nv,
                             spatial.field.y1 = 1:nv,
                             iid.y1 = 1:nv
                             #pre = 1/dpdf$e1
                           ))
  
  # Stack for predicting y2
  stk.y2.pre <- inla.stack(tag = "y2.pre",
                           data = list(y = matrix(NA, nrow=nrow(dp),ncol=n.pp),
                                       e = rep(0,nrow(dp))),
                           A = list(A.pr, A.pr, A.pr, A.pr),
                           effects = list(
                             Intercept.y2 = rep(1, nv),           
                             sh.copy.y1 = 1:nv,
                             spatial.field.y2 = 1:nv,
                             iid.y2 = 1:nv
                             #pre = 1/dpdf$e2
                           ))
  
  # Stack for predicting y3
  stk.y3.pre <- inla.stack(tag = "y3.pre",
                           data = list(y = matrix(NA, nrow=nrow(dp),ncol=n.pp),
                                       e = rep(0,nrow(dp))),
                           A = list(A.pr, A.pr, A.pr, A.pr, A.pr),
                           effects = list(
                             Intercept.y3 = rep(1, nv),
                             sh.copy.y1 = 1:nv,
                             y2.copy.y3 = 1:nv,
                             spatial.field.y3 = 1:nv,
                             iid.y3 = 1:nv
                             #pre = 1/dpdf$e3
                           ))
  
  # Stack with shared effect
  stk.shared <- inla.stack(tag = "shared",
                           data = list(y = matrix(NA, nrow=nrow(dp),ncol=n.pp)),
                           A = list(A.pr),
                           effects = list(spatial.field.sh = 1:nv))
  
  stk.y1.spec <- inla.stack(tag = "y1.spec",
                            data = list(y = matrix(NA, nrow=nrow(dp),ncol=n.pp)),
                            A = list(A.pr),
                            effects = list(spatial.field.y1 = 1:nv))
  
  stk.y2.spec <- inla.stack(tag = "y2.spec",
                            data = list(y = matrix(NA, nrow=nrow(dp),ncol=n.pp)),
                            A = list(A.pr),
                            effects = list(spatial.field.y2 = 1:nv))
  
  stk.y3.spec <- inla.stack(tag = "y3.spec",
                            data = list(y = matrix(NA, nrow=nrow(dp),ncol=n.pp)),
                            A = list(A.pr),
                            effects = list(spatial.field.y3 = 1:nv))
  
  stk.y1.iid <- inla.stack(tag = "y1.iid",
                           data = list(y = matrix(NA, nrow=nrow(dp),ncol=n.pp)),
                           A = list(A.pr),
                           effects = list(iid.y1 = 1:nv))
  
  stk.y2.iid <- inla.stack(tag = "y2.iid",
                           data = list(y = matrix(NA, nrow=nrow(dp),ncol=n.pp)),
                           A = list(A.pr),
                           effects = list(iid.y2 = 1:nv))
  
  stk.y3.iid <- inla.stack(tag = "y3.iid",
                           data = list(y = matrix(NA, nrow=nrow(dp),ncol=n.pp)),
                           A = list(A.pr),
                           effects = list(iid.y3 = 1:nv))
  
  join.stack <- inla.stack(
    stk.y1, stk.y2, stk.y3,
    stk.y1.pre, stk.y2.pre, stk.y3.pre,
    stk.shared, stk.y1.spec,stk.y2.spec, stk.y3.spec,
    stk.y1.iid, stk.y2.iid, stk.y3.iid)
  
  hyper.eps <- list(prec = list(prior = 'pc.prec', 
                                param = c(1, 0.01)))
  
  hyper <- list(beta = list(prior = 'normal', param = c(0, 10)))
  
  form <- y ~ 0 + Intercept.y1 + Intercept.y2 + Intercept.y3 + #offset(log(pre))+
    f(spatial.field.sh, model = spde) +
    f(spatial.field.y1, model = spde) +
    f(sh.copy.y2, copy = "spatial.field.sh", fixed = FALSE, hyper = hyper) +
    f(spatial.field.y2, model = spde) +
    f(sh.copy.y3, copy = "spatial.field.sh", fixed = FALSE, hyper = hyper) +
    f(y2.copy.y3, copy = "spatial.field.y2", fixed = FALSE, hyper = hyper) +
    f(spatial.field.y3, model = spde) +
    f(iid.y1, model = "iid", hyper = hyper.eps) +
    f(iid.y2, model = "iid", hyper = hyper.eps) +
    f(iid.y3, model = "iid", hyper = hyper.eps) 
  
  
  print("modeling")
  
  t1 <-Sys.time()
  res <- inla(
    formula=form, verbose = FALSE, 
    data = inla.stack.data(join.stack, spde = spde), 
    family = rep("poisson", 3),
    control.predictor = list(A = inla.stack.A(join.stack), 
                             compute = TRUE, link = 1),
    control.compute = list(dic = TRUE,config=TRUE),
    #control.mode = list(theta = theta.ini.e, restart = TRUE),
    #control.inla = list(int.strategy = 'eb'),
    E = inla.stack.data(join.stack)$e
  )
  t2 <- Sys.time()
  print("Time modeling")
  print(t2-t1)
  
  #res$mode$theta
  
  idx.y1 <- inla.stack.index(join.stack, 'y1.pre')$data
  idx.y2 <- inla.stack.index(join.stack, 'y2.pre')$data
  idx.y3 <- inla.stack.index(join.stack, 'y3.pre')$data
  
  idx.sh <- inla.stack.index(join.stack, 'shared')$data
  idx.y1.spec <- inla.stack.index(join.stack, 'y1.spec')$data
  idx.y2.spec <- inla.stack.index(join.stack, 'y2.spec')$data
  idx.y3.spec <- inla.stack.index(join.stack, 'y3.spec')$data
  
  idx.y1.iid <- inla.stack.index(join.stack, 'y1.iid')$data
  idx.y2.iid <- inla.stack.index(join.stack, 'y2.iid')$data
  idx.y3.iid <- inla.stack.index(join.stack, 'y3.iid')$data
  
  # Prediction
  
  ##dpdf$mean.y3 <- res$summary.fitted.values[idx.y3, "mean"]
  
  dpdf$mean.y1 <- res$summary.fitted.values[idx.y1, "mean"]
  dpdf$mean.y2 <- res$summary.fitted.values[idx.y2, "mean"]
  dpdf$mean.y3 <- res$summary.fitted.values[idx.y3, "mean"]
  
  dpdf$mean.l1 <- res$summary.linear.predictor[idx.y1, "mean"]
  dpdf$mean.l2 <- res$summary.linear.predictor[idx.y2, "mean"]
  dpdf$mean.l3 <- res$summary.linear.predictor[idx.y3, "mean"]
  
  dpdf$l1.lower <- res$summary.linear.predictor[idx.y1, "0.025quant"]
  dpdf$l2.lower <- res$summary.linear.predictor[idx.y2, "0.025quant"]
  dpdf$l3.lower <- res$summary.linear.predictor[idx.y3, "0.025quant"]
  
  dpdf$l1.upper <- res$summary.linear.predictor[idx.y1, "0.975quant"]
  dpdf$l2.upper <- res$summary.linear.predictor[idx.y2, "0.975quant"]
  dpdf$l3.upper <- res$summary.linear.predictor[idx.y3, "0.975quant"]
  
  dpdf$shared <- res$summary.linear.predictor[idx.sh, "mean"]
  dpdf$y1.spec <- res$summary.linear.predictor[idx.y1.spec, "mean"]
  dpdf$y2.spec <- res$summary.linear.predictor[idx.y2.spec, "mean"]
  dpdf$y3.spec <- res$summary.linear.predictor[idx.y3.spec, "mean"]
  
  dpdf$y1.iid <- res$summary.linear.predictor[idx.y1.iid, "mean"]
  dpdf$y2.iid <- res$summary.linear.predictor[idx.y2.iid, "mean"]
  dpdf$y3.iid <- res$summary.linear.predictor[idx.y3.iid, "mean"]
  
  p.sd <- lapply(res$internal.marginals.hyperpar[9:11],varF)
  
  tabcrp1 <- cbind(true = alpha, res$summary.fixed[, c(1:3, 5)])
  # Precision of the errors
  tabcrp2 <- cbind(
    true = c(e = e.sd), 
    t(sapply(p.sd, function(m) 
      unlist(inla.zmarginal(m, silent = TRUE))[c(1:3, 7)])))
  colnames(tabcrp2) <- colnames(tabcrp1)
  # Copy parameters 
  tabcrp3 <- cbind(
    true = beta, res$summary.hyperpar[12:14, c(1:3, 5)])
  tabcrp <- rbind(tabcrp1, tabcrp2, tabcrp3)
  
  range.est <- c()
  for (name in c("sh","y1","y2","y3")){
    spde.est <- inla.spde2.result(inla = res, 
                                  name = paste0("spatial.field.",name),
                                  spde = spde, do.transf = TRUE)
    temp <- inla.zmarginal(spde.est$marginals.range.nominal[[1]], 
                           silent = TRUE)
    range.est <- rbind(range.est,unlist(temp)[c(1:3, 7)])
  }
  tabran1 <- cbind(true = c(ran = range),range.est)
  colnames(tabran1) <- colnames(tabcrp1)
  
  m.var.est <- c()
  for (name in c("sh","y1","y2","y3")){
    name = "y3"
    spde.est <- inla.spde2.result(inla = res, 
                                  name = paste0("spatial.field.",name),
                                  spde = spde, do.transf = TRUE)
    temp <- inla.zmarginal(spde.est$marginals.variance.nominal[[1]], 
                           silent = TRUE)
    m.var.est <- rbind(m.var.est,unlist(temp)[c(1:3, 7)])
  }
  tabran2 <- cbind(true = c(var = m.var),m.var.est)
  colnames(tabran2) <- colnames(tabcrp1)
  
  tabran <- rbind(tabran1,tabran2)
  
  #samples <- inla.hyperpar.sample(1000000,res)
  
  allEffects <- as.data.frame(rbind(tabcrp,tabran[c(1,5,2,6,3,7,4,8),]))
  
  write.csv(dpdf,paste0("./s.results/pred_",seed,".csv"), row.names = FALSE)
  write.csv(allEffects,paste0("./s.results/effects_",seed,".csv"), row.names = TRUE)
  st_write(map,paste0("./s.results/maps/map_",seed,".shp"), append = FALSE)
  
  sink(paste0("./s.results/texto_",seed,".txt"))
  print(summary(res))
  cat("Fixed Effects","\n","\n")
  print(tabcrp)
  cat("\n","Random Effects","\n","\n")
  print(tabran)
  cat("\n","Formula","\n","\n")
  print(form)
  cat("\n","time","\n","\n")
  print(t2-t1)
  sink()
  rm(res)
}
