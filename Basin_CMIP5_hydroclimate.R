#************************************************************************************************
# This code brakes down the "CMIP5_MultiModel_World_Historical" to avoid the frequently experienced crashing
# This code performs the analysis for each basin and model and the second code plots roses and Budyko
#************************************************************************************************
# Packages installation and activation
library(calibrate)
library("ncdf4")        # "ncdf4", lib.loc="/Library/Frameworks/R.framework/Versions/3.3/Resources/library")
require(raster)
require(rgdal)
require(maptools)
# options(error=recover)

# Clearing workspace -------------------------------------------------------------------------
rm(list=ls(all=TRUE)) ## Clear Environment
graphics.off()
Sys.time()->start;
# Server directories -------------------------------------------------------------------------
workdir = "~/CMIP5/Workdir/Models/"
basinFold = "~/CMIP5/Basins/"
ResultsFolder = "~/CMIP5/Results/New_PET_World_results/"

# Setting data --------------------------------------------------------------------------------
crs<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
start_cutyear_p<-1910 #end_cutyear_p<-1989
start_cutyear<-2006 #end_cutyear<-2099

# Largest world basins ------------------------------------------------------------------
setwd(basinFold)
basins_name="GRDC_405_basins_from_mouth"
Continents_borders="Continent_World.shp"
basins_shp<-readShapeSpatial(basins_name, proj4string=CRS(as.character(crs)))
World_shp<-readShapeSpatial(Continents_borders, proj4string=CRS(as.character(crs)))
World_line <- as(World_shp, "SpatialLines")
Area<-basins_shp$AREA_CALC
area_min<-10000
numBasins<-length(basins_shp)
Basin<-basins_shp$DRAINAGE
condition_bas<-(Area>=area_min) #(Basin=='KHATANGA')#

# Models ------------------------------------------------------------------------------
model = c("NorESM1-ME","NorESM1-M","MRI-CGCM3","MIROC5","CNRM-CM5","inmcm4",
          "CMCC-CM","IPSL-CM5A-MR","MPI-ESM-MR")
modnum<-length(model)
variable=c("evs","tas","pr","ps","hus","rsd","mrr","mrs")
varnum<-length(variable)

# Run the PETcalc function that computes the Penmen Montieth equation
# using surface air pressure (P), air temperature (T), near-surface specific humidity (Q),
# and surface downwelling shortwave radiation (Rn), assuming G=0 and wind (U=1)
PET_Penman <- function(Ps,T,q,Rn) 
{
  Lv <- 2.501-0.002361*T
  Ga <- 0.0016286*Ps/Lv
  Es <- 0.6108*exp(17.27*T/(T+237.3))
  E <- Ps*q/(0.622+0.378*q)
  D <- 4098*Es/((T+237.3)^2)
  PET_Pen <- ((D/(D+Ga))*Rn)+((Ga/(D+Ga))*(6.43*1.053*(Es-E))/Lv)
  return(PET_Pen)
}

PET_FAO<-function(Ps,T,q,Rn) 
{
  Ga <- 0.000665*Ps
  Es <- 0.6108*exp(17.27*T/(T+237.3))
  E <- Ps*q/(0.622+0.378*q)
  D <- 4098*Es/((T+237.3)^2)
  PET_FAO <- (0.408*D*Rn+Ga*(900/(T+273))*(Es-E))/(D+1.34*Ga)
  return(PET_FAO)
}
PET_MILLY<-function(Rn) {  PET_MIL <- 0.8*Rn
return(PET_MIL)}

# loop for models ----------------------------------------------------------------------------------------
for(mod in 1:modnum){
  print(model[mod])
  # Declare raster stacks
  P_stack <- stack()
  E_stack <- stack()
  T_stack <- stack()
  Ps_stack <- stack()
  q_stack <- stack()
  Rn_stack <- stack()
  R_stack <- stack()
  Sm_stack <- stack()
  var_list <- list()

  setwd(paste0(workdir, model[mod]))
  list_exp = c("historical", "RCP85")#
  numn_exp<-length(list_exp)
  
  for(exp in 1:numn_exp){
    # if (exp==2 && Cond_cell) {break}
    print(paste0("exp=",list_exp[exp]))
    setwd(paste0(workdir, model[mod],"/", list_exp[exp],"/"))
    
    for (i in 1:varnum) {
      list_ts<-list.files(path = ".", pattern = variable[i]) # List all files for a single variable
      num_ts<-length(list_ts) #number of files for a single variable
      
      # Get the file names
      fname_nc<-list_ts
      print(paste0("Variable=",variable[i])) #Check
      var<-stack(fname_nc)
      
      
      # Defineindices to cut the ruster stack---------------------------------------------------------
      if (num_ts>1) {
        # Extract dates and creates time vector
        date_start <- substr(list_ts[1], nchar(list_ts[1])-15, nchar(list_ts[1])-10)
        date_end <- substr(list_ts[num_ts], nchar(list_ts[num_ts])-8, nchar(list_ts[num_ts])-3)
        startyear <- substr(date_start, 1, 4)
        startmonth <- substr(date_start, 5, 6)
        endyear <- substr(date_end, 1, 4)
        endmonth <- substr(date_end, 5, 6)
      } else {
        date_lim <- substr(fname_nc[[1]], nchar(fname_nc[[1]])-15, nchar(fname_nc[[1]])-3)
        dates_split <-strsplit(date_lim, "-")
        date_mat <- matrix(unlist(dates_split), ncol=2, byrow=TRUE)
        startyear <- substr(date_mat[1], 1, 4)
        startmonth <-substr(date_mat[1], 5, 6)
        endyear <- substr(date_mat[2], 1, 4)
        endmonth <- substr(date_mat[2], 5, 6)
      }
      if (list_exp[exp]=="historical") { #Cut the raster stacks for historical sim, from 1910 to 2005
        startindex_cut_p <- abs(start_cutyear_p-as.numeric(startyear))*12+2-as.numeric(startmonth)
        endindex_cut_p <- startindex_cut_p+1151
        
        if(variable[i]=="pr"){
          P_stack <- stack(var[[startindex_cut_p:endindex_cut_p]]*30*24*60*60) #From kg m-2 s-1 to mm/month
        }
        if(variable[i]=="evs"){
          E_stack <- stack(var[[startindex_cut_p:endindex_cut_p]]*30*24*60*60)
        }
        if(variable[i]=="tas"){
          T_stack <- stack(var[[startindex_cut_p:endindex_cut_p]]-273.15) # from K to C
        }
        if(variable[i]=="ps"){
          Ps_stack <- stack(var[[startindex_cut_p:endindex_cut_p]]/1000)
        }
        if(variable[i]=="hus"){ # Specific humidity
          q_stack <- stack(var[[startindex_cut_p:endindex_cut_p]])
        }
        if(variable[i]=="rsd"){ # Solar radiation
          Rn_stack <- stack(var[[startindex_cut_p:endindex_cut_p]]*0.03520667)
        }
        if(variable[i]=="mrr"){ # Runoff
          R_stack <- stack(var[[startindex_cut_p:endindex_cut_p]]*30*24*60*60)
        }
        if(variable[i]=="mrs"){ # Soil moisture
          Sm_stack <- stack(var[[startindex_cut_p:endindex_cut_p]]) #???????????????
        }
      } else {
        # "RCP85": Cut from 2006 to 2099
        startindex_cut <- abs(start_cutyear-as.numeric(startyear))*12+2-as.numeric(startmonth)
        endindex_cut <- startindex_cut+1127
        
        if(variable[i]=="pr"){
          P_stack <- stack(P_stack, var[[startindex_cut:endindex_cut]]*30*24*60*60)
        }
        if(variable[i]=="evs"){
          E_stack <- stack(E_stack, var[[startindex_cut:endindex_cut]]*30*24*60*60)
        }
        if(variable[i]=="tas"){
          T_stack<- stack(T_stack, var[[startindex_cut:endindex_cut]]-273.15)
        }
        if(variable[i]=="ps"){
          Ps_stack <- stack(Ps_stack, var[[startindex_cut:endindex_cut]]/1000)
        }
        if(variable[i]=="hus"){ # Specific humidity
          q_stack <- stack(q_stack, var[[startindex_cut:endindex_cut]]) #adimensionale
        }
        if(variable[i]=="rsd"){ # Solar radiation
          Rn_stack <- stack(Rn_stack, var[[startindex_cut:endindex_cut]]*0.03520667)
        }
        if(variable[i]=="mrr"){ # Runoff
          R_stack <- stack(R_stack, var[[startindex_cut:endindex_cut]]*30*24*60*60)
        }
        if(variable[i]=="mrs"){ # Soil moisture
          Sm_stack <- stack(Sm_stack, var[[startindex_cut:endindex_cut]]) #???????????????
        }
      }
    } 
    rm(var)
    gc()
  }
  # PET
  PET_Pen_stack <- PET_Penman(Ps=Ps_stack, T=T_stack, q=q_stack, Rn=Rn_stack)*30
  PET_FAO_stack <- PET_FAO(Ps=Ps_stack, T=T_stack, q=q_stack, Rn=Rn_stack)*30
  PET_MIL_stack <- PET_MILLY(Rn=Rn_stack)*30

  # Rotare raster stacks and put in a list ---------------------------------------------------
  P_stackR <- rotate(P_stack)
  E_stackR <- rotate(E_stack)
  T_stackR <- rotate(T_stack)
  R_stackR <- rotate(R_stack)
  Sm_stackR <- rotate(Sm_stack)
  PET_Pen_stackR <- rotate(PET_Pen_stack)
  PET_FAO_stackR <- rotate(PET_FAO_stack)
  PET_MIL_stackR <- rotate(PET_MIL_stack)
  
  # Exclude ocean cells ---------------------------------------------------------------------
  ras_col<-Ps_stack@ncols
  ras_rows<-Ps_stack@nrows
  r <- raster(ncol=ras_col, nrow=ras_rows)
  rp <- rasterize(World_line, r)
  values(P_stackR)[!is.na(values(rp))] <- NA
  values(E_stackR)[!is.na(values(rp))] <- NA
  values(T_stackR)[!is.na(values(rp))] <- NA
  values(R_stackR)[!is.na(values(rp))] <- NA
  values(Sm_stackR)[!is.na(values(rp))] <- NA
  values(PET_Pen_stackR)[!is.na(values(rp))] <- NA
  values(PET_FAO_stackR)[!is.na(values(rp))] <- NA
  values(PET_MIL_stackR)[!is.na(values(rp))] <- NA
  
  # Calculate annual mean and cum values
  var_list[[1]]<-stackApply(P_stackR,indices=sort(rep(as.factor(1910:2099),12),decreasing=FALSE),fun=sum,na.rm=TRUE)
  var_list[[2]]<-stackApply(E_stackR,indices=sort(rep(as.factor(1910:2099),12),decreasing=FALSE),fun=sum,na.rm=TRUE)
  var_list[[3]]<-stackApply(R_stackR,indices=sort(rep(as.factor(1910:2099),12),decreasing=FALSE),fun=sum,na.rm=TRUE)
  var_list[[4]]<-stackApply(Sm_stackR,indices=sort(rep(as.factor(1910:2099),12),decreasing=FALSE),fun=mean,na.rm=TRUE)
  var_list[[5]]<-stackApply(T_stackR,indices=sort(rep(as.factor(1910:2099),12),decreasing=FALSE),fun=mean,na.rm=TRUE)
  var_list[[6]]<-stackApply(PET_Pen_stackR,indices=sort(rep(as.factor(1910:2099),12),decreasing=FALSE),fun=sum,na.rm=TRUE)
  var_list[[7]]<-stackApply(PET_FAO_stackR,indices=sort(rep(as.factor(1910:2099),12),decreasing=FALSE),fun=sum,na.rm=TRUE)
  var_list[[8]]<-stackApply(PET_MIL_stackR,indices=sort(rep(as.factor(1910:2099),12),decreasing=FALSE),fun=sum,na.rm=TRUE)
  
  # loop on basins
  for(bas in 1:numBasins){
    if (condition_bas[bas]){
 
      print(paste0("-----------------, bas=",bas))
      ID_bas<-basins_shp$BASIN_ID[bas]
      selected_basin<-basins_shp[basins_shp$BASIN_ID == ID_bas[1],]
      River<-basins_shp$DRAINAGE[bas]
      # print(paste0("------------------------------------------------",River))
      
      ## Variable declaration ------------------------------------------------------------------
      num_var <- length(var_list)
      VAR_bas<-matrix(nrow = 190, ncol = num_var)
      
      for(v in 1:num_var){
        vF<-extract(var_list[[v]][[1]],selected_basin,weights=TRUE,cellnumbers=TRUE)
        
        ## A helper function that tests whether an object is either NULL _or_
        ## a list of NULLs
        is.NullOb <- function(x) is.null(x) | all(sapply(x, is.null))
        
        ## Recursively step down into list, removing all such objects
        rmNullObs <- function(x) {
          x <- Filter(Negate(is.NullOb), x)
          lapply(x, function(x) if (is.list(x)) rmNullObs(x) else x)
        }
        
        vnew<-rmNullObs(vF)
        notNA_position<-which(!is.na(vnew[[1]][,'value'])) #gives the position (row) in the dataframe vnew of a not NA value
        Cond_cell<-(all(is.na(vnew[[1]][,'value']))) # All values are NA
        if(Cond_cell){
          VAR_bas[,v] <- NA
          break
        } else if(all(!is.na(vnew[[1]]))) { # All values are not NA
          fun <- function(x) apply((var_list[[v]][vnew[[x]][,'cell']] * vnew[[1]][,'weight']) / sum(vnew[[1]][,'weight']), 2, sum, na.rm=TRUE)
          VAR_bas[,v]<- sapply(1:length(vnew), fun)
        } else {
          fun <- function(x) apply((var_list[[v]][vnew[[x]][notNA_position,'cell']] * vnew[[1]][notNA_position,'weight'])/ sum(vnew[[1]][notNA_position,'weight']), 2, sum, na.rm=TRUE)
          VAR_bas[,v]<- sapply(1:length(vnew), fun)
        }
      }
      #Save data --------------------------------------------------------------------------------
      dir.create(file.path(ResultsFolder, model[mod],"/"), showWarnings = FALSE)
      setwd(file.path(ResultsFolder, model[mod],"/"))
      
      write.table(VAR_bas, file = paste0(ID_bas,"_",model[mod],".txt"))
      
      # Remove temporary variables
      removeTmpFiles()
      gc()
    }
  }
}

#################################################################################################
# Write result in different files and directories
################################################################################################
# setwd(basinFold)
# basins_name="GRDC_405_basins_from_mouth"
# basins_shp<-readShapeSpatial(basins_name, proj4string=CRS(as.character(crs)))
# numBasins<-length(basins_shp)
# Basin<-basins_shp$DRAINAGE
print(Sys.time()-start)

setwd(paste0(ResultsFolder, model[1],"/"))
list_bas<-list.files() # List all files for a single variable
num_bas<-length(list_bas) #number of files for a single variable
bas_name <- list_bas
# modnum <- 10
for (b in 1:num_bas) {
  if (condition_bas[b]){
    River<-basins_shp$DRAINAGE[b]
    ID_bas<-basins_shp$BASIN_ID[b]
    
    P_files <- matrix(nrow = 190, ncol = modnum)
    E_files <- matrix(nrow = 190, ncol = modnum)
    R_files <- matrix(nrow = 190, ncol = modnum)
    Sm_files <- matrix(nrow = 190, ncol = modnum)
    T_files <- matrix(nrow = 190, ncol = modnum)
    PET_Pen_files <- matrix(nrow = 190, ncol = modnum)
    PET_FAO_files <- matrix(nrow = 190, ncol = modnum)
    PET_MIL_files <- matrix(nrow = 190, ncol = modnum)
    
    for(mod in 1:modnum){
      
      setwd(paste0(ResultsFolder, model[mod],"/"))
      file_name<-list.files(path = ".", pattern = as.character(paste0("^",ID_bas,"_")))
      txt_mod <- read.table(file_name)
      
      P_files[,mod] <- txt_mod[,1]
      E_files[,mod] <- txt_mod[,2]
      R_files[,mod] <- txt_mod[,3]
      Sm_files[,mod] <- txt_mod[,4]
      T_files[,mod] <- txt_mod[,5]
      PET_Pen_files[,mod] <- txt_mod[,6]
      PET_FAO_files[,mod] <- txt_mod[,7]
      PET_MIL_files[,mod] <- txt_mod[,8]
    }
    # WRITE AS TXT
    #Save data --------------------------------------------------------------------------------
    dir.create(file.path("~/CMIP5/Results/R10_NEW/", River,"/"), showWarnings = FALSE)
    setwd(file.path("~/CMIP5/Results/R10_NEW/", River,"/"))
    
    write.table(P_files , file = paste0(River,"_P.txt"))
    write.table(E_files , file = paste0(River,"_E.txt"))
    write.table(R_files , file = paste0(River,"_R.txt"))
    write.table(Sm_files , file = paste0(River,"_Sm.txt"))
    write.table(T_files , file = paste0(River,"_T.txt"))
    write.table(PET_Pen_files , file = paste0(River,"_PET_Pen.txt"))
    write.table(PET_FAO_files , file = paste0(River,"_PET_FAO.txt"))
    write.table(PET_MIL_files , file = paste0(River,"_PET_MIL.txt"))
  }
}

