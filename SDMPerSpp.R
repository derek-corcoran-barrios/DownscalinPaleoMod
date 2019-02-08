library(dismo)
library(raster)
library(tidyverse)
library(readxl)

## Hippidion saldiasi, Hippidion principale 

###Leo la base de datos, quito los registros con edades mayoers a 21.000 años y dejo las variables de interes para modelar
Presencias <- read_csv("Equidae.csv") %>% dplyr::filter(age_average <= 21000) %>% dplyr::select(taxon_name, age_average, lng, lat) 

## Creo un dataframe con las especies que tienen mas de 6 presencias
ForFilter <- Presencias %>% group_by(taxon_name) %>%  summarise(n = n()) %>% filter(n > 6 & taxon_name != "Hippidion")

#Filtro las presencias según el data frame for filter y lo divido en una lista de Data.frames según especie

Presencias <- Presencias %>% filter(taxon_name %in% ForFilter$taxon_name) %>% split(.$taxon_name)

SppNames <- names(Presencias)

#Tambien leo la carpeta donde estan las variables bioclimaticas para usarla en el loop

Bios <- list.files("/home/derek/Documents/Downscaling/DownscaledLayers/Bio", full.names = T)

####Primero divido la base de datos por especie segun su edad, este codigo crea una lista, donde cada elemento es un dataframe con los datos de alguna edad en particualar

DFPreds <- Presencias %>%  map(~split(.x, .$age_average)) 

###Ahora con un loop generamos puntos de background (Coherente con Maxent) y de presencia
Backgrounds <- list()

for(x in 1:length(DFPreds)){
  Conditions <- list()
  for(i in 1:length(DFPreds[[x]])){
    #Primero leemos la capa del año con que trabajaremos para usarla de "Mascara" del tiempo i
    Bio <- read_rds(Bios[str_detect(Bios, paste0("Bio", DFPreds[[x]][[i]]$age_average))])
    ###leo las prescencias para trasformarlo en coordenadas
    Pres <- DFPreds[[x]][[i]]
    coordinates(Pres) <- ~lng+lat
    
    #Extraigo las condiciones de las presencias de este tiempo
    
    Condition <- as.data.frame(raster::extract(Bio,Pres))
    Condition <- cbind(DFPreds[[x]][[i]], Condition)
    Condition$Pres <- 1
    #Genero el background con la mascara (Para que no eliga NA de fondo), n = numero de muestras y que
    #no se solape con p (las prescencias)
    set.seed(2018)
    Bg <- as.data.frame(randomPoints(mask = Bio, n = 1000))
    #Extraigo las condiciones del background
    Fondo <- raster::extract(Bio,Bg)
    #Ahora este Background lo transformo en un formato compatible con DFPredHip para poder unirlo a esa base de datos
    Background <- data.frame(taxon_name = "background", age_average = unique(DFPreds[[x]][[i]]$age_average), lng = Bg$x, lat = Bg$y)
    Background <- cbind(Background, Fondo)
    Background$Pres <- 0
    #unimos todo
    
    Conditions[[i]] <- bind_rows(Condition, Background)
    
    ##Mensaje de duracion
    gc()
    message(paste(i, "of", length(DFPreds[[x]]), "ready"))
  }
  
  Conditions <- bind_rows(Conditions) %>% arrange(desc(taxon_name))
  Conditions <- Conditions %>% distinct()
  Backgrounds[[x]] <- Conditions
  print(paste(x, "of", length(DFPreds), "ready"))
}
beepr::beep(8)
saveRDS(Backgrounds, "Backgrounds.rds")

############Maxent usando trinary maps y CV 


library(maxnet)
library(dismo)
library(raster)
library(tidyverse)
library(rasterVis)
library(pROC)
library(trinaryMaps)
library(cmsdm)

Backgrounds <- readRDS("Backgrounds.rds")
Conditions <- list()
for(x in 1:length(Backgrounds)){
  ### set 5 folds for cv in prescences with more than 10 precenses
  if(nrow(filter(Backgrounds[[x]], Pres == 1)) > 10){
    ConditionsPres <- Backgrounds[[x]] %>% dplyr::filter(Pres == 1)
    set.seed(2019)
    ConditionsPres$Fold <- sample(rep(c(1:5), each = ceiling(nrow(ConditionsPres)/5)), nrow(ConditionsPres))
  
  ### set 3 folds for cv in absences
    ConditionsBack <-  Backgrounds[[x]] %>% dplyr::filter(Pres == 0)
    set.seed(2019)
    ConditionsBack$Fold <- sample(rep(c(1:5), each = ceiling(nrow(ConditionsBack)/5)), nrow(ConditionsBack))
  }else{
    ConditionsPres <- Backgrounds[[x]] %>% dplyr::filter(Pres == 1)
    set.seed(2019)
    ConditionsPres$Fold <- sample(rep(c(1:3), each = ceiling(nrow(ConditionsPres)/3)), nrow(ConditionsPres))
    
    ### set 3 folds for cv in absences
    ConditionsBack <-  Backgrounds[[x]] %>% dplyr::filter(Pres == 0)
    set.seed(2019)
    ConditionsBack$Fold <- sample(rep(c(1:3), each = ceiling(nrow(ConditionsBack)/3)), nrow(ConditionsBack))
  }
  #### Join them back together
  
  Conditions[[x]] <- bind_rows(ConditionsPres, ConditionsBack)
}


saveRDS(Conditions, "Conditions.rds")
gc()

##Ahora a hacer los modelos

Evaluacion <- list()
Evals2 <- list()
Modelos <- list()
for(x in 1:length(Conditions)){
Evals <- list()
Models <- list()
  for(i in 1:length(unique(Conditions[[x]]$Fold))){
    Train <- Conditions[[x]] %>% dplyr::filter(Fold != i)
    Test <- Conditions[[x]] %>% dplyr::filter(Fold == i)
    Modelo <- maxnet(Train$Pres, data = Train[,5:23])
    Thresholds <- trinaryMaps::trinaryROCRoots(data.frame(response = Conditions[[x]]$Pres, predictor = as.numeric(predict(Modelo, Conditions[[x]][,5:23], type ="cloglog"))))
    Evals[[i]] <- data.frame(Hi_CI = Thresholds$plotThings$threshHi
                           , Low_CI = Thresholds$plotThings$threshLo, TSS_Thres =Thresholds$plotThings$threshYouden, AUC = as.numeric(Thresholds$plotThings$a.pauc$auc)
  )
  
    Models[[i]] <- Modelo
    message(paste(i, "of", length(unique(Conditions[[x]]$Fold)) ,"ready"))
  }
  Evals <- bind_rows(Evals)
  Evaluacion[[x]] <- Evals
  Evals2[[x]] <- Evals %>% summarise_all(funs(weighted.mean(., w = AUC))) %>% mutate(Spp = SppNames[x])
  Modelos[[x]] <- Models
  print(paste(x, "of", length(Conditions), "ready"))
}

Evals2 <- bind_rows(Evals2)
saveRDS(Evals2, "Evals2.rds")
saveRDS(Modelos, "Modelos.rds")
saveRDS(Evaluacion, "Evaluacion.rds")
beepr::beep(8)


#### Ahora proyectamos especie por especie a cada año

Predicciones <- list()

for(x in 1:length(Modelos)){
Predicted <- list()
Edades <- seq(21000, 8000, by = -1000)
  for (i in 1:length(Edades)){
    Area <- read_rds(Bios[str_detect(Bios, paste0("Bio", Edades[i]))])
    Predict <- list()
    for(j in 1:length(Modelos[[x]])){
      Predict[[j]] <- predict(Area, Modelos[[x]][[j]], type = "cloglog")
    }
    Predicted[[i]] <- weighted.mean(stack(unlist(Predict)), w =  Evaluacion[[x]]$AUC)
    message(paste(i, "of", length(Edades), "ready"))
  }

  Predicted <- do.call("stack", Predicted)
  names(Predicted) <- paste(SppNames[x], Edades)
  Predicciones[[x]] <- Predicted
  print(paste("Species",x, "of", length(SppNames), "ready!!"))
}

beepr::beep(8)

#########Ahora con 0 para ausencias, 1 para presencias seguras y 2 para presencias menos seguras segun trinary masp

#Transformamos los Stack en 2, 1 y 0
PredictedTrinary <- list()
PredictedHi <- list()
PredictedLo <- list()
PredictedTSS <- list()
for(x in 1:length(Predicciones)){
  PredictedTrinary[[x]] <- Predicciones[[x]]
  PredictedLo[[x]] <- Predicciones[[x]]
  PredictedHi[[x]] <- Predicciones[[x]]
  PredictedTSS[[x]] <- Predicciones[[x]]

  values(PredictedTrinary[[x]]) <- ifelse(values(PredictedTrinary[[x]]) < Evals2$Hi_CI[[x]], 0, 
                                   ifelse(values(PredictedTrinary[[x]]) < Evals2$Low_CI[[x]] & values(PredictedTrinary[[x]]) > Evals2$Hi_CI[[x]],1, 2))
  values(PredictedHi[[x]]) <- ifelse(values(PredictedHi[[x]]) < Evals2$Hi_CI[[x]], 0, 1)
  values(PredictedLo[[x]]) <- ifelse(values(PredictedLo[[x]]) < Evals2$Low_CI[[x]], 0, 1)
  values(PredictedTSS[[x]]) <- ifelse(values(PredictedTSS[[x]]) < Evals2$TSS_Thres[[x]], 0, 1)

  names(PredictedTrinary[[x]]) <- paste(SppNames[x], Edades)
  names(PredictedHi[[x]]) <- paste(SppNames[x], Edades)
  names(PredictedLo[[x]]) <- paste(SppNames[x], Edades)
  names(PredictedTSS[[x]]) <- paste(SppNames[x], Edades)
}

beepr::beep(8)


#Calculamos el Area de cada pixel

AreaKm <-raster::area(Predicciones[[1]][[1]])

##Generamos un Dataframe con el area por año para cada especie

Distribuciones <- list()

for(x in 1:length(Predicciones)){
  Distribucion <- list()
  
  for(i in 1:length(Edades)){
    Distribucion[[i]] <- data.frame(Edad = Edades[i], AreaLow = sum(values(AreaKm) * values(PredictedLo[[x]][[i]]), na.rm = T), AreaHi = sum(values(AreaKm) * values(PredictedHi[[x]][[i]]), na.rm = T), Area = sum(values(AreaKm) * values(PredictedTSS[[x]][[i]]), na.rm = T), Spp = SppNames[[x]])
    message(paste(i, "of", length(Edades), "ready"))
  }
  
  Distribucion <- bind_rows(Distribucion)
  Distribuciones[[x]] <- Distribucion
}

Distribuciones <- bind_rows(Distribuciones)
saveRDS(Distribuciones, "Distribuciones.rds")

ggplot(Distribuciones, aes(x = Edad, y = Area)) + geom_line(aes(color = Spp)) + theme_classic() + scale_x_reverse(labels = scales::comma) + theme(legend.position = "bottom") + scale_y_continuous(labels = scales::comma)

library(RColorBrewer)
cols <- brewer.pal(3, "RdBu")
levelplot(PredictedTSS[[4]], par.settings = BuRdTheme, colorkey= FALSE)
levelplot(PredictedTrinary[[4]], par.settings = BuRdTheme, colorkey= FALSE)

Diversity <- Reduce('+', PredictedTSS)
names(Diversity) <- paste("Diversity", Edades)
levelplot(Diversity, par.settings = BuRdTheme)

