# Best downscaling algorithm probado con temperatura minima de Enero LMG

##Load necesary packages

library(tidyverse)
library(raster)
library(rgdal)

## Read the paleoview layer to generate the extent

TminPaleoview <- raster("/home/derek/Documents/Downscaling/TempMinEnero/grid_data_21000BP.asc") %>% trim()

##A list is generated with the 12 layers of a variable (tmin, tmax or prec) in this case prec using the getData function from raster

Prec <- getData('worldclim', var='prec', res=2.5) %>% crop(TminPaleoview) %>% unstack()

### I generate a layer of how the coast was in the last Glacial maximum

NewGebco <- read_rds("/home/derek/Documents/Downscaling/NewGebco.rds")
MaxGlacier <- NewGebco
values(MaxGlacier) <- ifelse(values(MaxGlacier) < -125, NA,values(MaxGlacier))

#### Current day mask

Mask <- Prec[[1]]
values(Mask) <- ifelse(is.na(values(Mask)), NA, 1)

### Create layer of current elevation

Elev <- NewGebco
values(Elev) <- ifelse(values(Elev) < 0, NA,values(Elev))

### Generate a layer of latitude

LatitudTot <- NewGebco
data_matrix <- raster::xyFromCell(NewGebco, 1:ncell(NewGebco))
values(LatitudTot) <- data_matrix[,2]
Latitud <- LatitudTot*Mask

### Longitud

LongitudTot <- NewGebco
values(LongitudTot) <- data_matrix[,1]
Longitud <- LongitudTot*Mask
#generamos los ceoficientes para los 12 meses

##Creamos stack de predictores
Coeficientes <- list()
gc()
for(j in 1:length(Prec)){
  gc()
  Preds <- stack(Elev, Longitud, Latitud, (Prec[[j]] + 1))
  names(Preds) <- c("Elev", "Longitud", "Latitud", "Prec")

  ### Analisis de ventana mobil
  
  w <- c(25,25)

  # Luego generamos las capas para los estimadores

  intercept <- Preds[[1]]
  intercept[] <- NA
  elevationEst <- intercept
  latitudeEst <- intercept
  longitudeEst <- intercept


  for( rl in 1:nrow(Preds)) { 
    v <- getValuesFocal(Preds[[1:4]], row=rl, nrows=1, ngb = w, array = FALSE)
    int <- rep(NA,nrow(v[[1]]))
    x1 <- rep(NA,nrow(v[[1]]))
    x2 <- rep(NA,nrow(v[[1]]))
    x3 <- rep(NA,nrow(v[[1]]))
    x4 <- rep(NA,nrow(v[[1]]))
    for(i in 1:nrow(v[[1]]) ) {
      xy <- na.omit( data.frame(x1=v[[1]][i,],x2=v[[2]][i,], x3=v[[3]][i,],y=v[[4]][i,]) )
      # This tells as that there are more than 170 non-empty cells and less than 625
      if( nrow(xy) > 170 & nrow(xy) <= 624) {
        coefs <- coefficients(glm(as.numeric(xy$y) ~ as.numeric(xy$x1) + as.numeric(xy$x2) + as.numeric(xy$x3), family = poisson))
        int[i] <- coefs[1]
        x1[i] <- coefs[2]
        x2[i] <- coefs[3]
        x3[i] <- coefs[4]
      } else {
        int[i] <- NA
        x1[i] <- NA
        x2[i] <- NA
        x3[i] <- NA 
      }
    }
  
    intercept[rl,] <- int
    elevationEst[rl,] <- x1
    longitudeEst[rl,] <- x2
    latitudeEst[rl,] <- x3
    if(rl %% 100 == 0){
    message(paste(rl, "of", nrow(Preds), "ready"))
    }
  }

  Coeffs <- stack(intercept, elevationEst, latitudeEst, longitudeEst, exp(intercept + Preds$Elev*elevationEst  + Preds$Longitud*longitudeEst+ Preds$Latitud*latitudeEst), Preds$Prec)
  names(Coeffs) <- c("intercept", "elevationEst", "longitudeEst", "latitudeEst", "fitted", "Observed")
  Coeficientes[[j]] <- Coeffs
  saveRDS(Coeficientes, "CoeficientesPrec.rds")
  print(paste("Month", j , "of", length(Prec), "Ready!!", Sys.time()))
  gc()
}

#############Extend the coefficients to the last glacial maximum coastline

library(tidyverse)
library(raster)
library(rgdal)

FullCoeficientes <- list()

w=focalWeight(Mask, 0.5, "circle")

max.it = 1000

for(j in 1:length(Coeficientes)){
  gaps = which(is.na(Coeficientes[[j]][[1]] & is.na(Mask))[])
  intercept.ext = Coeficientes[[j]][[1]]
  elevation.ext = Coeficientes[[j]][[2]]
  latitude.ext = Coeficientes[[j]][[3]]
  longitude.ext = Coeficientes[[j]][[4]]
  for (i in 1:max.it) {
    intercept.ext[gaps] = focal(intercept.ext, w=w, na.rm=TRUE)[gaps]
    elevation.ext[gaps] = focal(elevation.ext, w=w, na.rm=TRUE)[gaps]
    latitude.ext[gaps] = focal(latitude.ext, w=w, na.rm=TRUE)[gaps]
    longitude.ext[gaps] = focal(longitude.ext, w=w, na.rm=TRUE)[gaps]
    if(i %% 10 == 0){
      print(paste(i, "of", max.it))
    }
  }
  intercept.ext = mask(intercept.ext, MaxGlacier)
  elevation.ext = mask(elevation.ext, MaxGlacier)
  latitude.ext = mask(latitude.ext, MaxGlacier)
  longitude.ext = mask(longitude.ext, MaxGlacier)

  ### Prediccion
  
  Coeffs2 <- stack(intercept.ext, elevation.ext, latitude.ext, longitude.ext)
  names(Coeffs2) <- c("intercept", "elevationEst", "longitudeEst", "latitudeEst")
  FullCoeficientes[[j]] <- Coeffs2
  saveRDS(FullCoeficientes, "FullCoeficientesTmax.rds")
  print(paste("Month", j , "of", length(Tmax), "Ready!!", Sys.time()))
  gc()
}
###Relleno
#Coeffs2 <- stack(intercept.ext, elevation.ext, latitude.ext, longitude.ext, (intercept.ext + (MaxGlacier+120)*elevation.ext  + LongitudTot*longitude.ext+ LatitudTot*latitude.ext), TMinWorldclim)
#names(Coeffs2) <- c("intercept", "elevationEst", "longitudeEst", "latitudeEst", "fitted", "Observed")
#TempPred <- Coeffs2$fitted
#values(TempPred) <- ifelse(is.na(values(TempPred)), values(TMinWorldclim), values(TempPred))
