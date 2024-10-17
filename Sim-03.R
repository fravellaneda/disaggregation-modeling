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
library(raster)

roi <- square(1) |>
  as.polygonal() |>
  as.data.frame() |>
  st_as_sf(coords = c("x","y")) |>
  st_combine() |>
  st_cast("POLYGON") |>
  as_Spatial()
boundary <- as(roi,"sf")

n_areas <- 100
fun_gen <- function(x,y) {n_areas*(3/2)*((x)**2+(y-0)**2)}

set.seed(1234)
p_gen <- rpoispp(fun_gen,win=square(1))
p_gen_sf <- st_as_sf(data.frame(x=p_gen$x,y=p_gen$y),coords = c("x","y"))
voronoi <- terra::voronoi(x = terra::vect(p_gen_sf),bnd = boundary)
voronoi$id <- 1:nrow(voronoi)
plot(voronoi)

resolution1 <- 50
x <- seq(1/(2*resolution1), 1-1/(2*resolution1), length.out = resolution1)
dp_pre <- as.matrix(expand.grid(x, x))
plot(dp_pre, asp = 1)

alpha <- c(5, 4, 3) 
m.var <- c(0.5,0.1,0.4, 0.3)
range <- c(0.1,0.4,0.3,0.2)

beta <- c(0.7, 0.5, -0.5) 
e.sd <- c(0.1, 0.2, 0.15)

book.rMatern <- function(n, coords, sigma=1, range, 
                         kappa = sqrt(8*nu)/range, variance = sigma^2, nu=1) {
  m <- as.matrix(dist(coords))
  m <- exp((1-nu)*log(2) + nu*log(kappa*m)-
             lgamma(nu))*besselK(m*kappa, nu)
  diag(m) <- 1
  return(drop(crossprod(chol(variance*m),
                        matrix(rnorm(nrow(coords)*n), ncol=n))))
}

i <- 9
set.seed(i)
z1 <- book.rMatern(1, dp_pre, range = range[1],
                   sigma = sqrt(m.var[1]))
z2 <- book.rMatern(1, dp_pre, range = range[2],
                   sigma = sqrt(m.var[2]))
z3 <- book.rMatern(1, dp_pre, range = range[3],
                   sigma = sqrt(m.var[3]))
z4 <- book.rMatern(1, dp_pre, range = range[4],
                   sigma = sqrt(m.var[4]))


theta1 <- exp(alpha[1] + z1 + z2 +rnorm(nrow(dp_pre), 0, e.sd[1]))
theta2 <- exp(alpha[2] + beta[1] * z1 + z3 + rnorm(nrow(dp_pre), 0, e.sd[2]))
theta3 <- exp(alpha[3] + beta[2] * z1 + beta[3] * z3 + z4 + 
            rnorm(nrow(dp_pre), 0, e.sd[3]))


dpdf2 <- as.data.frame(dp_pre)
names(dpdf2) <- c("x","y")

#Population 

zpop <- book.rMatern(1, dp_pre, range = 0.01,
                   sigma = sqrt(2))
dpdf2$pop <- round(exp(zpop)*100,0)
population <- rast(raster(as.im(list(x=x,y=x,z=matrix(dpdf2$pop,
                                                      nrow=resolution1)))))

ggplot(boundary) + geom_sf() + coord_sf(datum = NA) +
  geom_tile(data = dpdf2[,c("x","y","pop")],
            aes(x = x, y = y, fill = pop)) +
  labs(x = "", y = "") +
  scale_fill_viridis(paste0("RR ",i)) +
  theme_bw()

#Intensity

dpdf2$theta1 <- theta1
dpdf2$theta2 <- theta2
dpdf2$theta3 <- theta3
dpdf2$lambda1 <- dpdf2$pop*theta1
dpdf2$lambda2 <- dpdf2$pop*theta2
dpdf2$lambda3 <- dpdf2$pop*theta3

plot1 <- ggplot(boundary) + geom_sf() + coord_sf(datum = NA) +
  geom_tile(data = dpdf2[,c("x","y","lambda1")],
            aes(x = x, y = y, fill = lambda1)) +
  labs(x = "", y = "") +
  scale_fill_viridis(paste0("RR ",i)) +
  theme_bw()

plot2 <- ggplot(boundary) + geom_sf() + coord_sf(datum = NA) +
  geom_tile(data = dpdf2[,c("x","y","lambda2")],
            aes(x = x, y = y, fill = lambda2)) +
  labs(x = "", y = "") +
  scale_fill_viridis(paste0("RR ",i)) +
  theme_bw()

plot3 <- ggplot(boundary) + geom_sf() + coord_sf(datum = NA) +
  geom_tile(data = dpdf2[,c("x","y","lambda3")],
            aes(x = x, y = y, fill = lambda3)) +
  labs(x = "", y = "") +
  scale_fill_viridis(paste0("RR ",i)) +
  theme_bw()
pplotm1 <- grid.arrange(plot1,plot2, plot3, nrow=1)

#sir continous
dpdf2$e1 <- expected(population = dpdf2$pop,
                     cases = dpdf2$lambda1, n.strata = 1)
dpdf2$e2 <- expected(population = dpdf2$pop,
                     cases = dpdf2$lambda2, n.strata = 1)
dpdf2$e3 <- expected(population = dpdf2$pop,
                     cases = dpdf2$lambda3, n.strata = 1)
dpdf2$sir1 <- dpdf2$lambda1/dpdf2$e1
dpdf2$sir2 <- dpdf2$lambda2/dpdf2$e2
dpdf2$sir3 <- dpdf2$lambda3/dpdf2$e3

## paso simulación
set.seed(1234)
im.1 <- as.im(list(x=x,y=x,z=matrix(dpdf2$lambda1,nrow=resolution1)))
cases.y1 <- st_as_sf(coords(rpoispp(im.1)),coords = c("x","y"))
cases.y1$y1 <- 1
cases.y1$y2 <- 0
cases.y1$y3 <- 0

im.2 <- as.im(list(x=x,y=x,z=matrix(dpdf2$lambda2,nrow=resolution1)))
cases.y2 <- st_as_sf(coords(rpoispp(im.2)),coords = c("x","y"))
cases.y2$y1 <- 0
cases.y2$y2 <- 1
cases.y2$y3 <- 0

im.y3 <- as.im(list(x=x,y=x,z=matrix(dpdf2$lambda3,nrow=resolution1)))
cases.y3 <- st_as_sf(coords(rpoispp(im.y3)),coords = c("x","y"))
cases.y3$y1 <- 0
cases.y3$y2 <- 0
cases.y3$y3 <- 1

cases <- rbind(cases.y1,cases.y2,cases.y3)

#Voronoi + Cases
voronoi_sf <- st_as_sf(voronoi)
voronoi_sf$pop <- terra::extract(population, voronoi_sf,
                                 sum, na.rm = TRUE)$layer
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
plot(map)

#Points inside each area
loc.d <- cbind(st_coordinates(st_point_on_surface(map))[, 1],
               st_coordinates(st_point_on_surface(map))[, 2])
loc2.d <- cbind(st_coordinates(boundary)[, 1],
                st_coordinates(boundary)[, 2])

resolution2 <- 50
x <- seq(1/(2*resolution2), 1-1/(2*resolution2), length.out = resolution2)
dp <- as.matrix(expand.grid(x, x))
bnd <- inla.nonconvex.hull(dp)

summary(dist(loc.d))

mesh <- inla.mesh.2d(
  loc = loc.d, boundary = bnd,
  max.edge = c(0.05,0.1)
)
plot(mesh)
points(loc.d, col = "red")
lines(loc2.d, lwd=2)

mesh_df <- as.data.frame(mesh$loc[,1:2])
names(mesh_df) <- c("x","y")

spde <- inla.spde2.matern(mesh = mesh, alpha = 2, constr = TRUE)

# Number of mesh nodes
nv <- mesh$n

# number of groups
n.pp <- 3

n.ar <- nrow(map)
y.1 <- matrix(NA, nrow = n.ar, ncol = n.pp)
y.1[, 1] <- map$y1
e.1 <- map$E1

y.2 <- matrix(NA, nrow = n.ar, ncol = n.pp)
y.2[, 2] <- map$y2
e.2 <- map$E2

y.3 <- matrix(NA, nrow = n.ar, ncol = n.pp)
y.3[, 3] <- map$y3
e.3 <- map$E3

# Projection matrix

mesh_sf <- st_as_sf(data.frame(x=mesh$loc[,1],y=mesh$loc[,2],id=1:nv),
                    coords = c("x","y"))

A.m <- matrix(0, n.ar, nv)

t1 <-Sys.time()
for (i in 1:n.ar){
  inter_sf <- st_intersection(mesh_sf,map[i,])
  idx <- inter_sf$id
  n.t <- length(idx)
  v <- terra::voronoi(x = terra::vect(mesh_sf[idx,]),
                      bnd = map[i,])
  v_sf <- st_as_sf(v)
  v_int <- st_intersection(v_sf,map[i,"id"])
  v_int$pop <- terra::extract(population, v_int,mean, na.rm = TRUE)$layer
  v_int$propp <- v_int$pop/sum(v_int$pop)
  
  v_int <- v_int %>% arrange(id)
  
  A.m[i,idx] <- v_int$propp
}
t2 <-Sys.time()
print(t2-t1)


# Stack for estimation
# Create the stack for y1
stk.y1 <- inla.stack(
  data = list(y = y.1, e = e.1),
  A = list(A.m,A.m,A.m,1),
  effects = list(
    spatial.field.sh = 1:nv,
    spatial.field.y1 = 1:nv,
    iid.y1 = 1:nv,
    list(Intercept.y1 = rep(1, n.ar))
  ),
  tag = "y1")

# Create the stack for y2
stk.y2 <- inla.stack(
  data = list(y = y.2, e = e.2), 
  A = list(A.m,A.m,A.m,1), 
  effects = list(
    sh.copy.y2 = 1:nv,
    spatial.field.y2 = 1:nv,
    iid.y2 = 1:nv,
    list(Intercept.y2 = rep(1, n.ar))
  ),
  tag = "y2")

# Create the stack for y3
stk.y3 <- inla.stack(
  data = list(y = y.3, e = e.3), 
  A = list(A.m,A.m,A.m,A.m,1), 
  effects = list(
    sh.copy.y3 = 1:nv,
    y2.copy.y3 = 1:nv,
    spatial.field.y3 = 1:nv,
    iid.y3 = 1:nv,
    list(Intercept.y3 = rep(1, n.ar))
  ),
  tag = "y3")

A.pr <- inla.spde.make.A(mesh = mesh, loc = dp)

# Stack for predicting y1
stk.y1.pre <- inla.stack(tag = "y1.pre",
                         data = list(y = matrix(NA, nrow=nrow(dp),ncol=n.pp),
                                     e = rep(0,nrow(dp))),
                         A = list(1, A.pr, A.pr, A.pr),
                         effects = list(
                           data.frame(Intercept.y1 = rep(1, nrow(dp))
                           ),
                           spatial.field.sh = 1:nv,
                           spatial.field.y1 = 1:nv,
                           iid.y1 = 1:nv
                         ))

# Stack for predicting y2
stk.y2.pre <- inla.stack(tag = "y2.pre",
                         data = list(y = matrix(NA, nrow=nrow(dp),ncol=n.pp),
                                     e = rep(0,nrow(dp))),
                         A = list(1, A.pr, A.pr, A.pr),
                         effects = list(
                           data.frame(Intercept.y2 = rep(1, nrow(dp))
                           ),           
                           sh.copy.y1 = 1:nv,
                           spatial.field.y2 = 1:nv,
                           iid.y2 = 1:nv
                           ))

# Stack for predicting y3
stk.y3.pre <- inla.stack(tag = "y3.pre",
                         data = list(y = matrix(NA, nrow=nrow(dp),ncol=n.pp),
                                     e = rep(0,nrow(dp))),
                         A = list(1, A.pr, A.pr, A.pr, A.pr),
                         effects = list(
                           data.frame(Intercept.y3 = rep(1, nrow(dp))
                           ),
                           sh.copy.y1 = 1:nv,
                           y2.copy.y3 = 1:nv,
                           spatial.field.y3 = 1:nv,
                           iid.y3 = 1:nv
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

theta.ini <- (c(c(log(range), log(sqrt(m.var))), log(1 / e.sd^2),
               beta))[c(1,5,2,6,3,7,4,8,9,10,11,12,13,14)]
theta.ini.e = theta.ini + rnorm(length(theta.ini), 0, 0.1)

hyper <- list(beta = list(prior = 'normal', param = c(0, 10)))

form <- y ~ -1 + Intercept.y1 + Intercept.y2 + Intercept.y3 + 
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
  control.mode = list(theta = theta.ini.e, restart = TRUE),
  control.inla = list(int.strategy = 'eb'),
  E = inla.stack.data(join.stack)$e
)
t2 <- Sys.time()
print(t2-t1)


summary(res)

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

dpdf <- data.frame(dp)
names(dpdf) <- c("x","y")

# Prediction

##dpdf$mean.y3 <- res$summary.fitted.values[idx.y3, "mean"]

dpdf$mean.y1 <- res$summary.fitted.values[idx.y1, "mean"]
dpdf$mean.y2 <- res$summary.fitted.values[idx.y2, "mean"]
dpdf$mean.y3 <- res$summary.fitted.values[idx.y3, "mean"]
dpdf$shared <- res$summary.linear.predictor[idx.sh, "mean"]
dpdf$y1.spec <- res$summary.linear.predictor[idx.y1.spec, "mean"]
dpdf$y2.spec <- res$summary.linear.predictor[idx.y2.spec, "mean"]
dpdf$y3.spec <- res$summary.linear.predictor[idx.y3.spec, "mean"]
dpdf$y1.iid <- res$summary.linear.predictor[idx.y1.iid, "mean"]
dpdf$y2.iid <- res$summary.linear.predictor[idx.y2.iid, "mean"]
dpdf$y3.iid <- res$summary.linear.predictor[idx.y3.iid, "mean"]

plot1 <- ggplot(map) +
  geom_tile(data = dpdf[,c("x","y","y2.spec")],
            aes(x = x, y = y, fill = y2.spec)) +
  geom_sf(fill=NA,color="black") + coord_sf(datum = NA) +
  labs(x = "", y = "") +
  scale_fill_viridis("Pred") +
  theme_bw()
plot1

plot2 <- ggplot(map) + 
  geom_tile(data = dpdf2[,c("x","y","sir1")],
            aes(x = x, y = y, fill = sir1)) +
  geom_sf(fill=NA,color="black") + coord_sf(datum = NA) +
  labs(x = "", y = "") +
  scale_fill_viridis("True") +
  theme_bw()

plot3 <- ggplot(map) + 
  geom_tile(data = dpdf2[,c("x","y","pop")],
            aes(x = x, y = y, fill = pop)) +
  geom_sf(fill=NA,color="black") + coord_sf(datum = NA) +
  labs(x = "", y = "") +
  scale_fill_viridis("True") +
  theme_bw()

grid.arrange(plot1,plot2,plot3, nrow=1)

p.sd <- lapply(res$internal.marginals.hyperpar[9:11],
               function(m) {
                 inla.tmarginal(function(x) 1 / sqrt(exp(x)), m)
               })

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
# The complete table
tabcrp <- rbind(tabcrp1, tabcrp2, tabcrp3)
tabcrp

tabran <- cbind(true=theta.ini[1:8], res$summary.hyperpar[1:8, c(1:3, 5)])
tabran
