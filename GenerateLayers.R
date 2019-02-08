library(tidyverse)
library(readxl)
library(raster)

Generated_years <- list.files(path = "DownscaledLayers/Tmin") %>% str_remove("Tmin") %>% str_remove(".rds") %>% as.numeric()

Humans <- read_excel("chile.xlsx", sheet = "Sheet1")


Humans <- Humans %>% dplyr::filter(str_detect(Full.Reference, "Gayo") | str_detect(Full.Reference, "Latorre")| str_detect(Full.Reference, "Santoro"))


Equus <- read_excel("~/Documents/Downscaling/Equus(Dic18).xlsx", sheet = "Sheet1")

Hippidion <- read_excel("~/Documents/Downscaling/Hippidion(Dic18).xlsx")

Years <- sort(unique(c(unique(Equus$age_average),unique(Hippidion$age_average), seq(21000, 8000, -1000))), decreasing = T)

Years <- Years[Years <= 21000]

HumanYears <- round(Humans$C14_date, -1) %>% unique()

Years <- HumanYears[!(HumanYears %in% Generated_years)]

###Profundidades
NewGebco <- read_rds("/home/derek/Documents/Downscaling/NewGebco.rds")

Niveles <- data.frame(YBP = seq(21000, 8000, by = -10), SeaLevel = NA)

Niveles$SeaLevel[1] <- -125

for(i in 2:nrow(Niveles)){
  if (Niveles$YBP[i] > 17000){
    Niveles$SeaLevel[i] = Niveles$SeaLevel[i - 1] + 0.06
  }
  if (Niveles$YBP[i] <= 17000){
    Niveles$SeaLevel[i] = Niveles$SeaLevel[i - 1] + 0.1
  }
  
}

Niveles <- dplyr::filter(Niveles, YBP %in% Years)

##Mascaras de profundidades
Masks <- list()

for(i in 1:nrow(Niveles)){
  Masks[[i]] <- NewGebco
  values(Masks[[i]]) <- ifelse(values(Masks[[i]]) < Niveles$SeaLevel[i], NA,1)
}

Masks <- do.call("stack", Masks)
names(Masks) <- paste("ybp",Years)

#Generar stacks por año

###Capas tmax

TmaxCoeffs <- read_rds("FullCoeficientesTmax2.rds")
TmaxFiles <- list.files(path = "/home/derek/Documents/Downscaling/MaxTemp21_10kyr", pattern = ".asc", recursive = T, full.names = T)
TmaxNow <- getData('worldclim', var='tmax', res=2.5) %>% crop(NewGebco)
TmaxNow <- TmaxNow

####Preds

### Latitud

LatitudTot <- NewGebco
data_matrix <- raster::xyFromCell(NewGebco, 1:ncell(NewGebco))
values(LatitudTot) <- data_matrix[,2]

### Longitud

LongitudTot <- NewGebco
values(LongitudTot) <- data_matrix[,1]

##Creamos stack de predictores

Preds <- stack(NewGebco, LongitudTot, LatitudTot)
names(Preds) <- c("Elev", "Longitud", "Latitud")


STACKS <- list()
for(i in 64:length(Years)){
  TmaxFilesTemp <- (TmaxFiles[str_detect(TmaxFiles, pattern = paste0("grid_data_",Years[i]))])[c(5, 4, 8, 1, 9,7, 6, 2, 12, 11, 10, 3)]
  StackPaleoViewTmax <- list()
  for (j in 1:length(TmaxFilesTemp)){
    StackPaleoViewTmax[[j]] <- raster(TmaxFilesTemp[j]) %>%  crop(NewGebco) %>% resample(NewGebco)
    Added <- (TmaxCoeffs[[j]][[1]] + (NewGebco+Niveles$SeaLevel[i])*TmaxCoeffs[[j]][[2]]  + LongitudTot*TmaxCoeffs[[j]][[4]]+ LatitudTot*TmaxCoeffs[[j]][[3]])
    values(Added) <- ifelse(is.na(values(Added)), values(TmaxNow[[j]]), values(Added))
    Added <- Added/10
    StackPaleoViewTmax[[j]] <- StackPaleoViewTmax[[j]] + Added
  }
  STACKS[[i]] <- do.call("stack", StackPaleoViewTmax)
  STACKS[[i]] <- STACKS[[i]]*Masks[[i]]
  names(STACKS[[i]]) <- paste0("Tmax", 1:12)
  saveRDS(STACKS[[i]], paste0("/home/derek/Documents/Downscaling/DownscaledLayers/Tmax/","Tmax", Years[i], ".rds"))
  message(paste(i, "from", length(Years), "ready"))
  gc()
}

############################################################
############################################################
##################Tmin########################################3
################################################################

###Capas Tmin

TminCoeffs <- read_rds("FullCoeficientesTmin.rds")
TminFiles <- list.files(path = "/home/derek/Documents/Downscaling/MinTemp21_10kyr", pattern = ".asc", recursive = T, full.names = T)
TminNow <- getData('worldclim', var='tmin', res=2.5) %>% crop(NewGebco)

####Preds

### Latitud

LatitudTot <- NewGebco
data_matrix <- raster::xyFromCell(NewGebco, 1:ncell(NewGebco))
values(LatitudTot) <- data_matrix[,2]

### Longitud

LongitudTot <- NewGebco
values(LongitudTot) <- data_matrix[,1]

##Creamos stack de predictores

Preds <- stack(NewGebco, LongitudTot, LatitudTot)
names(Preds) <- c("Elev", "Longitud", "Latitud")


STACKS <- list()
for(i in 1:length(Years)){
  TminFilesTemp <- (TminFiles[str_detect(TminFiles, pattern = paste0("grid_data_",Years[i]))])[c(5, 4, 8, 1, 9,7, 6, 2, 12, 11, 10, 3)]
  StackPaleoViewTmin <- list()
  for (j in 1:length(TminFilesTemp)){
    StackPaleoViewTmin[[j]] <- raster(TminFilesTemp[j]) %>%  crop(NewGebco) %>% resample(NewGebco)
    Added <- (TminCoeffs[[j]][[1]] + (NewGebco+Niveles$SeaLevel[i])*TminCoeffs[[j]][[2]]  + LongitudTot*TminCoeffs[[j]][[4]]+ LatitudTot*TminCoeffs[[j]][[3]])
    values(Added) <- ifelse(is.na(values(Added)), values(TminNow[[j]]), values(Added))
    Added <- Added/10
    StackPaleoViewTmin[[j]] <- StackPaleoViewTmin[[j]] + Added
  }
  STACKS[[i]] <- do.call("stack", StackPaleoViewTmin)
  STACKS[[i]] <- STACKS[[i]]*Masks[[i]]
  names(STACKS[[i]]) <- paste0("Tmin", 1:12)
  saveRDS(STACKS[[i]], paste0("/home/derek/Documents/Downscaling/DownscaledLayers/Tmin/","Tmin", Years[i], ".rds"))
  message(paste(i, "from", length(Years), "ready"))
  gc()
}

########################################################
#######################PRec############################
######################################################3
#############################################3

PrecCoeffs <- read_rds("FullCoeficientesPrec.rds")
PrecFiles <- list.files(path = "/home/derek/Documents/Downscaling/MeanPP21_10kyr", pattern = ".asc", recursive = T, full.names = T)
PrecNow <- getData('worldclim', var='prec', res=2.5) %>% crop(NewGebco)

####Preds

### Latitud

LatitudTot <- NewGebco
data_matrix <- raster::xyFromCell(NewGebco, 1:ncell(NewGebco))
values(LatitudTot) <- data_matrix[,2]

### Longitud

LongitudTot <- NewGebco
values(LongitudTot) <- data_matrix[,1]

##Creamos stack de predictores

Preds <- stack(NewGebco, LongitudTot, LatitudTot)
names(Preds) <- c("Elev", "Longitud", "Latitud")


STACKS <- list()
for(i in 1:length(Years)){
  PrecFilesTemp <- (PrecFiles[str_detect(PrecFiles, pattern = paste0("grid_data_",Years[i]))])[c(5, 4, 8, 1, 9,7, 6, 2, 12, 11, 10, 3)]
  StackPaleoViewPrec <- list()
  for (j in 1:length(PrecFilesTemp)){
    StackPaleoViewPrec[[j]] <- raster(PrecFilesTemp[j]) %>%  crop(NewGebco) %>% resample(NewGebco)
    Added <- exp(PrecCoeffs[[j]][[1]] + (NewGebco+Niveles$SeaLevel[i])*PrecCoeffs[[j]][[2]]  + LongitudTot*PrecCoeffs[[j]][[4]]+ LatitudTot*PrecCoeffs[[j]][[3]])
    values(Added) <- ifelse(is.na(values(Added)), values(PrecNow[[j]]), values(Added))
    Added <- Added
    StackPaleoViewPrec[[j]] <- StackPaleoViewPrec[[j]] + Added
    values(StackPaleoViewPrec[[j]]) <- ifelse(values(StackPaleoViewPrec[[j]]) < 0, 0,values(StackPaleoViewPrec[[j]])) 
  }
  STACKS[[i]] <- do.call("stack", StackPaleoViewPrec)
  STACKS[[i]] <- STACKS[[i]]*Masks[[i]]
  names(STACKS[[i]]) <- paste0("Prec", 1:12)
  saveRDS(STACKS[[i]], paste0("/home/derek/Documents/Downscaling/DownscaledLayers/Prec/","Prec", Years[i], ".rds"))
  message(paste(i, "from", length(Years), "ready"))
  gc()
}


###################Variables biolcimáticas

library(dismo)

TminFilesRds <- list.files("/home/derek/Documents/Downscaling/DownscaledLayers/Tmin", full.names = T)
TmaxFilesRds <- list.files("/home/derek/Documents/Downscaling/DownscaledLayers/Tmax", full.names = T)
PrecFilesRds <- list.files("/home/derek/Documents/Downscaling/DownscaledLayers/Prec", full.names = T)


Bios <- list()

for(i in 5:length(Years)){
  print(paste("Starting", i, "of", length(Years), Sys.time()))
  Tmin <- read_rds(TminFilesRds[str_detect(string = TminFilesRds, pattern = paste0("Tmin", Years[i]))])
  Tmax <- read_rds(TmaxFilesRds[str_detect(string = TmaxFilesRds, pattern = paste0("Tmax", Years[i]))])
  Prec <- read_rds(PrecFilesRds[str_detect(string = PrecFilesRds, pattern = paste0("Prec", Years[i]))])
  Bios[[i]] <- dismo::biovars(prec = Prec, tmin = Tmin, tmax = Tmax)
  saveRDS(Bios[[i]], paste0("/home/derek/Documents/Downscaling/DownscaledLayers/Bio/","Bio", Years[i], ".rds"))
  gc()
}
