#Checks if packages are installed and installs them if they are not already installed
list.of.packages <- c("INLA", "geoR", "gstat")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(INLA)
library(geoR)
library(gstat)

#Contains data on 8 variables and 2035 children living in 65 villages.
#Response variable is (pos) and is a binary indicator of the presence of
#malarial parasites in a blood sample.
#Child level covariates are: age (age, in days), usage of bed nets (netuse),
#, and information about whether the bed nets are treated with insecticide
#(treated).
#Village level covariates are: vegetation index (green) and the inclusion
#or not of the village in the primary health care system (phc)
data(gambia)
View(gambia)

#Create coordinates from gambia dataset (in km)
coords <- as.matrix(gambia[,1:2])/1000

#Create index at village level from coords variable
ind <- paste("x", coords[,1], "y", coords[,2], sep="")
#Detect non-duplicated villages
#Notes first row of new coordinate vector
which.nodupl <- which(!duplicated(ind))
#Creates vector with first row NA of length nrow(gambia)
village.index <- c(NA, length=nrow(gambia))
#Populates V1 column in village.index with 1 for first set of coordinates
#2 for next set of duplicated coordinates...up to 64
village.index[1 : (which.nodupl[length(which.nodupl)]-1)] <- rep(1:64, times=as.numeric(diff(which.nodupl)))
#Populates 65th village for final set of coordinates
village.index[which.nodupl[length(which.nodupl)] : nrow(gambia)] <- 65
#Create new column in gambia dataset populated by village.index
gambia$village.index <- village.index

#The SPDE approach is based on the triangulation of the spatial domain.

#Creates a polygon to be used as the boundary of the mesh, extending the domain
#by a distance 'convex' from the given points, while keeping the outside curvature
#radius to at least 'concave' (concave=convex by default)
bnd <- inla.nonconvex.hull(coords, convex=-0.1)
#Creates a triangle mesh based on initial point locations as defined by the boundary.
#The definition of the mesh is a trade-off between the accuracy of the GMRF
#representation and computational costs, both depending on the number of vertices
#used in the triangulation
gambia.mesh <- inla.mesh.2d(boundary=bnd,
                            offset=c(30,60), max.edge=c(20,40))
#Graphical representation of the mesh
#The internal line defines the nonconvex hull
plot(gambia.mesh, main="")
#Include data locations
points(coords, pch=21, bg=1, col="white", cex=1.8)
#Creates Matérn SPDE object (this is also where the amount of vertices is stored: 'n.spde')
gambia.spde <- inla.spde2.matern(mesh=gambia.mesh, alpha=2)
#Creates the sparse weight matrix A by identifying the data locations in the mesh
#and organizing the corresponding values of the basis functions. The sparse matrix A
#maps the GMRF from the G triangulation vertices to the n observation locations
A.est <- inla.spde.make.A(mesh=gambia.mesh, loc=coords)
#Creates all the required indexes for the SPDE model and which is particularly
#useful when replicates or gorups are involved in the model specification
s.index <- inla.spde.make.index(name="spatial.field",
                                n.spde=gambia.spde$n.spde)
#Manages the SPDE objects.
gambia.stack.est <- inla.stack(data=list(y=gambia$pos),
                               A=list(A.est,1,1,1,1,1,1),
                               effects=list(c(s.index, list(Intercept=1)),
                                            list(age=gambia$age/365),
                                            list(treated=gambia$treated),
                                            list(netuse=gambia$netuse),
                                            list(green=gambia$green),
                                            list(phc=gambia$phc),
                                            list(village.index=gambia$village.index)),
                               tag="est")
#Define the formula including all of the effects.
formula <- y ~ -1 + Intercept + treated + netuse + age + green + phc +
  f(spatial.field, model=gambia.spde) +
  f(village.index, model="iid")

gambia.output <- inla(formula,
                      data=inla.stack.data(gambia.stack.est, spde=gambia.spde),
                      family="binomial", Ntrials=1,
                      control.predictor=list(A=inla.stack.A(gambia.stack.est),
                                             compute = TRUE),
                      control.compute=list(dic=TRUE))

fixed.out <- round(gambia.output$summary.fixed[,1:5],3)
View(fixed.out)
#This shows a progressive increase in prevalence with age, and the protective effects
#of bed nets.
#Neither the inclusion in the primary health care system nor the greenness of the
#surrounding space appeared to affect significantly the prevalence of malaria
#The data supports the presence of a spatial effect since a smaller
#DIC (deviance information criterion) are better supported by the data

gambia.output$dic$dic


# Model WITHOUT the spatial effect
gambia.stack.est.noGF <- inla.stack(data=list(y=gambia$pos),
                                    A=list(1, 1, 1, 1, 1, 1, 1),
                                    effects=
                                      list(list(Intercept=rep(1,nrow(gambia))),
                                           list(age=gambia$age/365),
                                           list(treated=gambia$treated),
                                           list(netuse=gambia$netuse),
                                           list(green=gambia$green),
                                           list(phc=gambia$phc),
                                           list(village.index=gambia$village.index)),
                                    tag="est")

formula <- y ~ -1 + Intercept + treated + netuse + age + green + phc + f(village.index,model="iid")

gambia.output.noGF <- inla(formula,
                           data=inla.stack.data(gambia.stack.est.noGF, spde=gambia.spde),
                           family="binomial",Ntrials=1,
                           control.predictor=list(A=inla.stack.A(gambia.stack.est.noGF), compute=TRUE),
                           control.compute=list(dic=TRUE))
#Outputs DIC for non spatial effects model
gambia.output.noGF$dic$dic6