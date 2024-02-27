
### R script to fit gravity-style catchment model to health facility service populations by admin 2 units for DRC (AIM 1)

getwd ()
setwd("C:\\Users\\alyss\\OneDrive\\Documents\\GitHub\\Dissertation\\Data")
Cstack_info()

library(raster)
library(rgdal)
library(rgeos)
library(sp)
library(concaveman)
library(viridis)
library(sf)
library(doParallel)
library(foreach)
library(gdistance)
library(english)

## set up for parallel computation:
# cores_to_use <- 32
# cl <- makeCluster(cores_to_use)
# registerDoParallel(cl)
# stopCluster(cl)

#registerDoParallel(32) # update?
registerDoParallel(cores=20) #current laptop seems to have 24 cores

service_pops <- read.csv("DRC Air de Sante Population from DHIS2 2020.csv")

### source: DHIS2 / Tulane (priv. comm.)

facility.locations <- readOGR("org_units-2021-03-08-09-27.gpkg")
crs(facility.locations) <- "+proj=longlat +datum=WGS84 +no_defs"
fac.csv <- read.csv("org_units-2021-03-08-09-27.csv")

facility.locations$Ref.Ext.parent.1 <- fac.csv$Ref.Ext.parent.1[match(facility.locations$name,fac.csv$Nom)]


additional.facility.locations <- readOGR("health.gpkg")
crs(additional.facility.locations) <- "+proj=longlat +datum=WGS84 +no_defs"
additional.facility.locations <- additional.facility.locations[as.character(additional.facility.locations$amenity) %in% c("doctors","hospital","clinic","health_post"),]
### source: bluesquare: prior comm.
### source: bluesquare: s3://hexa-tulane-playground/GeoHealthAccess/COD/Input/health.gpkg
### remove near matches = suspected duplicates



min.dist <- numeric(length(facility.locations$uuid))
for (i in 1:length(additional.facility.locations$name)) {
  min.dist[i] <- min(sqrt((additional.facility.locations@coords[i,1]-facility.locations@coords[,1])^2+(additional.facility.locations@coords[i,2]-facility.locations@coords[,2])^2)*111*1000)
}
additional.facility.locations <- additional.facility.locations[-which(min.dist < 500),]

##########IMPORTING MAP FRICTION SURFACE
friction <- raster("MAP_friction_car_or_walk.tif")

### source: MAP (1km resolution)


########### IMPORTING BLUESQUARE FICTION SURFACE
#frictionBSQ <-raster("friction_walk_lzw_BSQ.tif")
friction <-raster ("202001_Global_Motorized_Travel_Time_to_Healthcare_COD.tiff")
population <- raster("cod_ppp_2020_UNadj_constrained_1km_matched.tif")
### source: https://www.worldpop.org/geodata/summary?id=49683
### matched to friction grid as below
# population_cod <- raster("cod_ppp_2020_UNadj_constrained.tif")
# pop.vals <- getValues(population_cod)
# match.up <- cellFromXY(friction,xyFromCell(population_cod,which(!is.na(pop.vals))))
# aggregated.matches <- aggregate(pop.vals[!is.na(pop.vals)],list(match.up),sum)
# aggregated.matches <- aggregated.matches[aggregated.matches[,2]>0,]
# population <- friction*NA
# population[aggregated.matches[,1]] <- aggregated.matches[,2]
# writeRaster(population,file="cod_ppp_2020_UNadj_constrained_1km_matched.tif",overwrite=TRUE)

admin <- readOGR("cod_admbnda_adm0_rgc_itos_20190911.shp")

friction[is.na(friction)] <- median(getValues(mask(friction,admin)),na.rm=TRUE)
  

joint.facility.locations.cells <- cellFromXY(friction,rbind(facility.locations@coords[,1:2],additional.facility.locations@coords))
joint.facility.locations.cells.unique <- unique(joint.facility.locations.cells) 
facility.locations.aggregation.list <- match(cellFromXY(friction,facility.locations@coords[,1:2]),joint.facility.locations.cells.unique)
### at 1km resolution we can have a few facilities sharing the same pixel, meaning that
### they're basically indistinguishable in the catchment model
joint.facility.locations.unique <- xyFromCell(friction,joint.facility.locations.cells.unique)
N.unique.fac <- length(joint.facility.locations.unique[,1])

T <- gdistance::transition(friction, function(x) 1/mean(x), 8)
T.GC <- gdistance::geoCorrection(T)
  
in.country.nonzero.pop.pixels <- which(!is.na(getValues(population)))
rowcols <- rowColFromCell(population,1:length(population))
n.rows <- nrow(population)
n.cols <- ncol(population)

############# build a list containing for each health facility the travel times from nearby pixels within a suitable window############
# code isn't working from this point out

#selection= indexing vector or array that specifies which elements of traveltime you're interested in. The values in selection should correspond to the positions or indices in the traveltime vector/array that you want to extract. 
#For example, if selection contains [2, 5, 7], you're interested in extracting the 2nd, 5th, and 7th elements of traveltime.

window.size <- 50 #max(50,min(150,3*as.integer(n.rows*n.cols/N.unique.fac/2)))
catchment.list <- foreach (i = 1:N.unique.fac) %do% {  
  row.selection <- max(rowcols[joint.facility.locations.cells.unique[i],1]-window.size,1):min(n.rows,rowcols[joint.facility.locations.cells.unique[i],1]+window.size)
  col.selection <- max(rowcols[joint.facility.locations.cells.unique[i],2]-window.size,1):min(n.cols,rowcols[joint.facility.locations.cells.unique[i],2]+window.size)
  selection <- cbind(rep(row.selection,times=length(col.selection)),rep(col.selection,each=length(row.selection))) # If selection contains row or column indices that are outside the raster's extent, it could result in an empty extraction
  friction.alt <- friction*NA
  friction.alt[selection] <- friction[selection] # take only pixels lying in DRC from within a window around that facility
  
  postage.stamp.T <- gdistance::transition(friction.alt, function(x) 1/mean(x), 8)
  postage.stamp.T.GC <- gdistance::geoCorrection(postage.stamp.T)
  
  traveltime <- gdistance::accCost(postage.stamp.T.GC, joint.facility.locations.unique[i,])/60/2 #traveltime= vector or array-like object that contains data values relative to friction surface derived travel times to each HF
  traveltime.selection <- traveltime[selection] # uses the indices in selection to extract specific elements from traveltime. The result is a new vector or array that contains only the elements of traveltime at the positions specified by selection.
  inclusion <- which(traveltime.selection<Inf & population[selection] > 0) #ERROR:object is empty/no data
  selection <- selection[inclusion,] #### ERROR: #no data stored
  traveltime.selection <- traveltime.selection[inclusion] #extracted elements are assigned to a new variable named traveltime.selection #ERROR:object is empty/no data stored- seems as if this issue is introduced in line 122?
  catchment.list.sub <- list()
  catchment.list.sub$traveltime <- traveltime.selection
  if (length(selection)>2) {
    catchment.list.sub$rows <- selection[,1]
    catchment.list.sub$cols <- selection[,2]} else {
    catchment.list.sub$rows <- selection[1]
    catchment.list.sub$cols <- selection[2]
  }

  catchment.list.sub
}

N.bound <- 20 # ordinal of furthest facility at which treatment might be sought
  
### summarise the catchment model in terms of the nearest N.bound facilities for each pixel
Nvalid.pix <- length(in.country.nonzero.pop.pixels)
nearest.hf.id.dist <- matrix(10^4,nrow=Nvalid.pix,ncol=2) 
for (i in 1:N.unique.fac) {
  # we want to extract from the catchment list the id number of the nearest health facility 
  # and its corresponding travel time distance 
  cellids <- match(cellFromRowCol(friction,catchment.list[[i]]$rows,catchment.list[[i]]$cols),in.country.nonzero.pop.pixels)
  replacements <- which(catchment.list[[i]]$traveltime < nearest.hf.id.dist[cellids,2])
  nearest.hf.id.dist[cellids[replacements],1] <- i
  nearest.hf.id.dist[cellids[replacements],2] <- catchment.list[[i]]$traveltime[replacements]
}

nlist <- c("",paste0(stringr::str_replace(as.character(ordinal(2:N.bound)),"-",""),"."))
for (name in nlist[-1]) {
  eval(parse(text=paste0(name,"nearest.hf.id.dist <- matrix(10^4,nrow=Nvalid.pix,ncol=2)")))
  for (i in 1:N.unique.fac) {
    cellids <- match(cellFromRowCol(friction,catchment.list[[i]]$rows,catchment.list[[i]]$cols),in.country.nonzero.pop.pixels)
    eval(parse(text=paste0("replacements <- which(catchment.list[[i]]$traveltime < ",name,"nearest.hf.id.dist[cellids,2] & ",paste0(nlist[1:match(name,nlist)-1],"nearest.hf.id.dist[cellids,1]!=i",collapse = " & "),")")))
    eval(parse(text=paste0(name,"nearest.hf.id.dist[cellids[replacements],1] <- i")))
    eval(parse(text=paste0(name,"nearest.hf.id.dist[cellids[replacements],2] <- catchment.list[[i]]$traveltime[replacements]")))
  }
}


#################### SETTING UP BASE ENVIRONMENT TO RUN PYTHON IN R ##################################
#Sys.setenv(RETICULATE_PYTHON = "/home/ecameron/anaconda3/bin/python3")
#Sys.setenv function is used to set environment variables within an R session. Environment variables are key-value pairs that store information about the environment 
#in which a program is running. They are often used to configure the behavior of programs or to store information that needs to be accessed by different parts 
#of a program. Code above is what Ewan used

Sys.setenv(RETICULATE_PYTHON = "C:\\Users\\alyss\\anaconda3\\python.exe") #filepath where python is located 
library(reticulate) #allows use of python script within R
use_python("C:\\Users\\alyss\\anaconda3\\python.exe", required= TRUE) #makes sure that R is calling correct version of Pything

torch <- import("torch") 


### note: torch syntax in R is name$function_name compared with names.function_name in python
### also: R has indices 1:N, python is 0:N; hence there is some care to be taken with this in translating the code below

fac_matches <- match(facility.locations$Ref.Ext.parent.1,service_pops$organisationunitid)
x_service_pops <- service_pops$Population.de.l.Aire.de.Santé[fac_matches]
x_rep_service_pops <- aggregate(x_service_pops,list(facility.locations.aggregation.list),sum)
valid_pops <- as.integer(!is.na(x_rep_service_pops[,2]))
x_rep_service_pops[valid_pops==0,2] <- 1
rep_service_pops <- torch$tensor(x_rep_service_pops[,2],dtype=torch$float)

### build the row-col indexing for a sparse catchment matrix  
vnlist <- paste0(stringr::str_replace(as.character(ordinal(1:N.bound)),"-",""),"_v")
catchment_matrix_rowcollist <- list()
k <- 1
id <- numeric(Nvalid.pix*10)
dist_ordered <- numeric(Nvalid.pix*10)
for (i in 1:Nvalid.pix) { 
    
  for (qq in 1:N.bound) {
    eval(parse(text=paste0("if (",nlist[qq],"nearest.hf.id.dist[i,2]<10^4) {",vnlist[qq]," <- ",nlist[qq],"nearest.hf.id.dist[i,1]-1} else {",vnlist[qq],"<- N.unique.fac+qq-1}")))
  }
    
  v.group <- eval(parse(text=paste0("v.group <- c(",paste0(vnlist,collapse=","),")")))
  rank.ordered <- sort.list(v.group)
  for (qq in 1:N.bound) {
    eval(parse(text=paste0("catchment_matrix_rowcollist[[k+rank.ordered[",qq,"]-1]] <- c(i-1,v.group[",qq,"])")))
  }

  eval(parse(text=paste0("id[c(",paste0("k+",1:N.bound-1,collapse=","),")] <- rank.ordered")))
  eval(parse(text=paste0("dist_ordered[c(",paste0("k+",1:N.bound-1,collapse=","),")] <- c(",paste0(nlist,"nearest.hf.id.dist[i,2]",collapse="+0.1,"),")[rank.ordered]")))
  k <- k+N.bound
}

catchment_matrix_i <- torch$LongTensor(catchment_matrix_rowcollist)
mass_index <- catchment_matrix_i$t()[1]$detach()$numpy() # returns the labeling from health facility ids to column numbers

comb_dist <- torch$tensor(dist_ordered,dtype=torch$float)
pop <- population[in.country.nonzero.pop.pixels]
pop <- torch$tensor(pop,dtype=torch$float)

valid_pops <- torch$tensor(valid_pops,dtype=torch$float)

log_distance_scale_mean <- torch$tensor(log(2),dtype=torch$float,requires_grad=TRUE) # these are the variables of our pytorch model
log_distance_scale <- torch$tensor(rep(0,N.unique.fac+N.bound),dtype=torch$float,requires_grad=TRUE) # these are the variables of our pytorch model
log_masses <- torch$tensor(rep(0,N.unique.fac+N.bound),dtype=torch$float,requires_grad=TRUE)
log_err_scale <- torch$tensor(-2,dtype=torch$float,requires_grad=TRUE)
log_pop_adj <- torch$tensor(-1,dtype=torch$float,requires_grad=TRUE)
logit_threshold <- torch$tensor(-3,dtype=torch$float,requires_grad=TRUE)
logit_competition_factor <- torch$tensor(0,dtype=torch$float,requires_grad=TRUE)

log_masses_expanded <- log_masses[mass_index]
log_distance_scales_expanded <- log_distance_scale[mass_index]$add(log_distance_scale_mean)
  
z <- 1
rel.improvement <- 100000
old.val <- 10000000
mod.factor <- 1  
  
x_opt <- torch$optim$Adam(list(log_distance_scale_mean,log_pop_adj,logit_competition_factor,log_distance_scale,logit_threshold,log_masses),lr=0.01)

while (z < 1000 & abs(rel.improvement) > 0.1) {
    
  log_masses_expanded <- log_masses[mass_index]
  log_distance_scales_expanded <- log_distance_scale[mass_index]$add(log_distance_scale_mean)
    
  raw_attractiveness <- comb_dist$pow(log_distance_scales_expanded$exp()$mul(-1))
  raw_attractiveness_with_mass <- comb_dist$pow(log_distance_scales_expanded$exp()$mul(-1))$mul(log_masses_expanded$exp())
  norm <- raw_attractiveness_with_mass$reshape(list(Nvalid.pix,as.integer(N.bound)))$sum(as.integer(1))$expand(list(as.integer(N.bound),Nvalid.pix))$t()$reshape(list(as.integer(Nvalid.pix*N.bound),as.integer(1)))$squeeze()
  norm_competition_adjusted <- norm$mul(logit_competition_factor$mul(-1)$exp()$add(1)$pow(-1))$add(raw_attractiveness_with_mass$mul(logit_competition_factor$mul(-1)$exp()$add(1)$pow(-1)$mul(-1)$add(1)))
  normed_attractiveness <- raw_attractiveness$div(norm_competition_adjusted)$mul(comb_dist$pow(log_distance_scales_expanded$exp()$mul(-1)))
    
  threshold <- logit_threshold$mul(-1)$exp()$add(1)$pow(-1)
  normed_boundary_attractiveness <- normed_attractiveness$log()$sub(threshold$log())$mul(10)
  normed_boundary_attractiveness <- normed_boundary_attractiveness$div(normed_boundary_attractiveness$abs()$add(1))$add(1)$div(2)
    
  catchment_matrix_v <- normed_boundary_attractiveness
  catchment_matrix <- torch$sparse$FloatTensor(catchment_matrix_i$t(), catchment_matrix_v, torch$Size(c(as.integer(Nvalid.pix),as.integer(N.unique.fac+N.bound))))$t()
    
  pred_pops <- torch$sparse$mm(catchment_matrix,pop$unsqueeze(as.integer(1)))$squeeze()$mul(log_pop_adj$exp())
    
  nll <- torch$distributions$StudentT(5.0,pred_pops$add(1)$log()[sort(unique(facility.locations.aggregation.list))-1],log_err_scale$exp())$log_prob(rep_service_pops$add(1)$log())$mul(valid_pops)$sum()$mul(-1.0)
  nll <- nll$sub(torch$distributions$Normal(0,0.25)$log_prob(log_masses)$sum())
  nll <- nll$sub(torch$distributions$Normal(0,0.1)$log_prob(log_distance_scale)$sum())
  nll <- nll$sub(torch$distributions$Normal(-2,0.5)$log_prob(log_err_scale))
  nll <- nll$sub(torch$distributions$Normal(-2,0.5)$log_prob(logit_threshold))
  nll <- nll$sub(torch$distributions$Normal(log(2),0.1)$log_prob(log_distance_scale_mean))
  nll <- nll$sub(torch$distributions$Normal(0,0.05)$log_prob(log_pop_adj))    

  par(mai=c(1,1,1,1))
  plot(pred_pops[sort(unique(facility.locations.aggregation.list))-1]$detach()$numpy(),rep_service_pops$detach()$numpy(),log='xy',xlim=c(10^2,10^6),ylim=c(10^2,10^6),pch=19,cex=0.5,col=hsv(0,alpha=0.2))
    
  x_opt$zero_grad()
  nll$backward(retain_graph = TRUE)
  x_opt$step()
    
  cat(z," ",nll$detach()$numpy(),"\n")
  cat(sd(log_masses$detach()$numpy()),"\n")
  cat(sd(log_distance_scale$detach()$numpy()),"\n")
  z <- z + 1
  rel.improvement <- old.val-nll$detach()$numpy()
  old.val <- nll$detach()$numpy()
  
}

log_masses_expanded <- log_masses[mass_index]
log_distance_scales_expanded <- log_distance_scale[mass_index]$add(log_distance_scale_mean)
  
raw_attractiveness <- comb_dist$pow(log_distance_scales_expanded$exp()$mul(-1))
raw_attractiveness_with_mass <- comb_dist$pow(log_distance_scales_expanded$exp()$mul(-1))$mul(log_masses_expanded$exp())
norm <- raw_attractiveness_with_mass$reshape(list(Nvalid.pix,as.integer(N.bound)))$sum(as.integer(1))$expand(list(as.integer(N.bound),Nvalid.pix))$t()$reshape(list(as.integer(Nvalid.pix*N.bound),as.integer(1)))$squeeze()
norm_competition_adjusted <- norm$mul(logit_competition_factor$mul(-1)$exp()$add(1)$pow(-1))$add(raw_attractiveness_with_mass$mul(logit_competition_factor$mul(-1)$exp()$add(1)$pow(-1)$mul(-1)$add(1)))
normed_attractiveness <- raw_attractiveness$div(norm_competition_adjusted)$mul(comb_dist$pow(log_distance_scales_expanded$exp()$mul(-1)))
  
threshold <- logit_threshold$mul(-1)$exp()$add(1)$pow(-1)
normed_boundary_attractiveness <- normed_attractiveness$log()$sub(threshold$log())$mul(10)
normed_boundary_attractiveness <- normed_boundary_attractiveness$div(normed_boundary_attractiveness$abs()$add(1))$add(1)$div(2)
  
catchment_matrix_v <- normed_boundary_attractiveness
catchment_matrix <- torch$sparse$FloatTensor(catchment_matrix_i$t(), catchment_matrix_v, torch$Size(c(as.integer(Nvalid.pix),as.integer(N.unique.fac+N.bound))))$t()
  
pred_pops <- torch$sparse$mm(catchment_matrix,pop$unsqueeze(as.integer(1)))$squeeze()$mul(log_pop_adj$exp())

access_per_pix <- torch$sparse$mm(catchment_matrix$t(),torch$tensor(c(rep(1,N.unique.fac),rep(0,N.bound)),dtype=torch$float)$unsqueeze(as.integer(0))$t())
reference.image <- friction*NA
reference.image[in.country.nonzero.pop.pixels] <- as.numeric(access_per_pix$detach()$numpy())
writeRaster(reference.image,file="service_fitted_access_model.tif")

### Admin level aggregation

admin <- readOGR("rdc_zones-de-sante/RDC_Zones de sante.shp")
Nadmin <- length(admin)

psums <- foreach (i=1:Nadmin) %dopar% {
  xpop <- mask(crop(population,admin[i,]),admin[i,])
  xfac <- mask(crop(reference.image,admin[i,]),admin[i,])
  x <- list()
  x[[1]] <- sum(getValues(xpop)*getValues(xfac),na.rm=TRUE)
  x[[2]] <- sum(getValues(xpop),na.rm=TRUE)
  x
}
pops <- numeric(Nadmin)
for (i in 1:Nadmin) {pops[i] <- psums[[i]][[2]]}
mean_facs <- numeric(Nadmin)
for (i in 1:Nadmin) {mean_facs[i] <- psums[[i]][[1]]/psums[[i]][[2]]}

pops[pops<1000] <- NA
mean_facs[is.na(pops)] <- NA
admin$MeanFac <- mean_facs
admin$pop <- pops

n.facs <- numeric(Nadmin)
for (i in 1:Nadmin) {n.facs[i] <- sum(point.in.polygon(joint.facility.locations.unique[,1],joint.facility.locations.unique[,2],admin@polygons[[i]]@Polygons[[1]]@coords[,1],admin@polygons[[i]]@Polygons[[1]]@coords[,2]))}
fac.dense <- n.facs/pops
fac.dense[is.na(pops)] <- NA
admin$facdense <- fac.dense

writeOGR(admin,".","zones_de_sante_service_coverage",driver="ESRI Shapefile",overwrite_layer = TRUE)

