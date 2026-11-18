
#Main Script for data analysis for Predator Dominance Analysis  

#PACKAGES----
library("corrplot")
library("dplyr")
library('lme4')
library("AICcmodavg")
library('arm')
library("nlme")
library("DHARMa")
library("glmmTMB")

#Load function:

invisible(lapply(paste("Functions/", dir("Functions"), sep=""), function(x) source(x)))


#DATA----

#All location relevant information has been removed from data to protect privacy

getwd()

point <- read.csv("01. Data/Point_Data.csv")
head(point);dim(point)

field <- read.csv("01. Data/Survey_Data.csv")
head(field);dim(field)

morpho <- read.csv("01. Data/Morpho_Data.csv")
head(morpho);dim(morpho)

obs <- read.csv("01. Data/Observation_Data.csv")
head(obs);dim(obs)


#DATA EXPLORATION----


summary(obs) #removing an unneeded column
obs <- obs[,-which(names(obs) =="Notes")]

length(unique(obs$Obs_ID))

unique(obs$Stage)

length(unique(obs$Morphospecies))
length(unique(morpho$Morphospecies))

length(unique(obs$Site))
length(unique(point$Site)) #this looks right since there where a few sites that I didn't find any insects


summary(morpho) #removing some unneeded columns
morpho <- morpho[,-which(names(morpho) =="NOTES")]
morpho <- morpho[,-which(names(morpho) =="In_Database")]
morpho <- morpho[,-which(names(morpho) =="Duplicates")]


unique(morpho$Trophic)

length(morpho$Trophic[morpho$Trophic == 'Parasite'])
length(morpho$Trophic[morpho$Trophic == 'Predator'])

unique(morpho$Order)

summary(point)

length(unique(point$Survey_Field))
length(field$ID)

hist(point$Plant_Height)
hist(point$Ground_Cover)

summary(field) #removing some unneeded columns
field <- field[,-which(names(field) =="Notes")]
field <- field[,-which(names(field) =="NDVIsum_Field")] #decided this wasn't a good variable since the field size differed
field <- field[,-c(7:16)] #removing climate data which isn't needed for analysis 
field <- field[,-which(names(field) =="X500m_Prop_Water")]
field <- field[,-which(names(field) =="X1km_Prop_Water")]

#exploring some of the variables

hist(field$Field_Area_m2)
hist(field$Crop_Age_Days)
hist(field$Average_Height)
hist(field$Average_GC)

table(field$X500m_Prop_Build) #high proportion of zeros makes it unsuitable for modelling 
field <- field[,-which(names(field) =="X500m_Prop_Build")]

table(field$X500m_Prop_Nat_Graze) #high proportion of zeros makes it unsuitable for modelling
field <- field[,-which(names(field) =="X500m_Prop_Nat_Graze")]

table(field$X500m_Prop_Water)
table(field$X500m_Prop_Crops)

table(field$X1km_Prop_Build)#Not enough variation in data - unsuitable for modelling
field <- field[,-which(names(field) =="X1km_Prop_Build")]

table(field$X1km_Prop_Nat)#high proportion of zeros makes it unsuitable for modelling 
field <- field[,-which(names(field) =="X1km_Prop_Nat")]

table(field$X1km_Rip_Prop)
table(field$X500m_Rip_Prop)
table(field$X1km_Prop_Crops)
table(field$X1km_Prop_Nat_Graze)

hist(point$Ground_Cover)
hist(point$Plant_Height)

#CORRELATION----

cordata2 <- data.frame(ID = point$Survey_Field,Height = point$Plant_Height,GC = point$Ground_Cover,Position = point$Spatial_Position)

cordata2 <- merge(cordata2,field, by = "ID")

cordata2 <- cordata2 %>% dplyr::select(Height, GC, Position,Day_Sampled,Crop_Age_Days,Field_Area_m2,X500m_Prop_Crops,X500m_Rip_Prop,X1km_Prop_Crops,X1km_Rip_Prop,X1km_Prop_Nat_Graze,NDVImean_Field,NDVIsum_500m,NDVImean_500m,NDVIsum_1km,NDVImean_1km)

cordata2$Position[cordata2$Position == 'Outer'] <- 1
cordata2$Position[cordata2$Position == 'Inner'] <- 2
cordata2$Position <- as.numeric(cordata2$Position)

str(cordata2)
cor2 <- cor(cordata2,method = "spearman")
colnames(cor2) <- c("Height", "GC", "Position","Day","Crop Age","Field Area","Crops 500m","Riparian 500m","Crops 1000m","Riparian 1000m","Graze 1000m","NDVI Field Mean","NDVI 500m Sum","NDVI 500m Mean","NDVI 1000m Sum","NDVI 1000m Mean")
rownames(cor2) <- c("Height", "GC", "Position","Day","Crop Age","Field Area","Crops 500m","Riparian 500m","Crops 1000m","Riparian 1000m","Graze 1000m","NDVI Field Mean","NDVI 500m Sum","NDVI 500m Mean","NDVI 1000m Sum","NDVI 1000m Mean")

dev.new(height=8,width=10,dpi=80,pointsize=14,noRStudioGD = T)
corrplot::corrplot(cor2,method="color",  
                   type="upper",addCoef.col = 'black',number.cex = 0.6)
head(cordata2)

#Uncorrelated variables listed below based on correlation plot

#Design variables: Crop age, day sampled and spatial position
#Point Variables: Ground cover, plant height
#Field variables: Field size, field NDVI (mean)
#Landscape variable: Riparian 1km, NDVI (sum) 1km


#DARA PREP----

##prep variables----
variables <- data.frame(Site = point$Site,Height = point$Plant_Height,GC = point$Ground_Cover,Position = point$Spatial_Position,ID = point$Survey_Field,X_Cor = point$X_Cor, Y_Cor = point$Y_Cor)

variables <- merge(variables,field, by = "ID")

head(variables);dim(variables)

variables <- variables %>% dplyr::select(ID,Site,Field,Survey,X_Cor, Y_Cor,Height, GC, Position,Day_Sampled,Crop_Age_Days,Field_Area_m2,X1km_Prop_Crops,X1km_Rip_Prop,NDVImean_Field,NDVIsum_1km)


##prep count data----
head(obs);dim(obs)
head(morpho);dim(morpho)

invert <- data.frame(ID = obs$Site,Morphospecies = obs$Morphospecies)

invert <- merge(invert,morpho,by = "Morphospecies")
head(invert);dim(invert)

#Removing wasps, lepidoptera and Serpintine Leaf Miner to control for chemicals applied during some of the sampling period

length(invert$Trophic[invert$Trophic == 'Parasite'])
invert <- subset(invert, Trophic != "Parasite")

table(invert$Order)
length(invert$Order[invert$Order == 'Lepidoptera'])
invert <- subset(invert, Order != "Lepidoptera")

length(invert$Morphospecies[invert$Morphospecies == 'Serpintine_Leaf_Miner'])
invert <- subset(invert, Morphospecies != "Serpintine_Leaf_Miner")


invert <- invert %>% dplyr::select(ID, Trophic)
head(invert);dim(invert)

#remove tropic levels which aren't predator or herbivore
unique(invert$Trophic)

length(invert$Trophic[invert$Trophic == 'Unknown'])
invert <- subset(invert, Trophic != "Unknown")

length(invert$Trophic[invert$Trophic == 'Ominvore'])
invert <- subset(invert, Trophic != "Ominvore")

length(invert$Trophic[invert$Trophic == 'Pollinator'])
invert <- subset(invert, Trophic != "Pollinator")

length(invert$Trophic[invert$Trophic == 'Fungivore'])
invert <- subset(invert, Trophic != "Fungivore")

length(invert$Trophic[invert$Trophic == 'Hematophagous'])
invert <- subset(invert, Trophic != "Hematophagous")

head(invert);dim(invert)

###Counting herbivores and predators for each site----

herb_count <- invert %>%
  filter(Trophic == "Herbivore") %>%
  group_by(ID) %>%
  summarize(count = n()) 
head(herb_count);dim(herb_count)

herb_count <- data.frame(herb_count)
names(herb_count)[names(herb_count) == "ID"] <- "Site"
str(herb_count)

pred_count <- invert %>%
  filter(Trophic == "Predator") %>%
  group_by(ID) %>%
  summarize(count = n()) 
head(pred_count);dim(pred_count)
pred_count <- data.frame(pred_count)
names(pred_count)[names(pred_count) == "ID"] <- "Site"
str(pred_count)

###Add counts to variables database----
mod_data <- merge(variables,herb_count,by = 'Site',all.x = T)
head(mod_data);dim(mod_data)
mod_data$count[is.na(mod_data$count)] <- 0
names(mod_data)[names(mod_data) == "count"] <- "Herb_Count"
str(mod_data)

mod_data <- merge(mod_data,pred_count,by = 'Site',all.x = T)
head(mod_data);dim(mod_data)
mod_data$count[is.na(mod_data$count)] <- 0
names(mod_data)[names(mod_data) == "count"] <- "Pred_Count"
str(mod_data)

###Predator dominance----

mod_data$Prop_Pred <- ifelse(mod_data$Pred_Count == 0, 0,mod_data$Pred_Count/(mod_data$Herb_Count+mod_data$Pred_Count))

head(mod_data);dim(mod_data)
summary(mod_data$Prop_Pred)
range(mod_data$Prop_Pred)

head(mod_data);dim(mod_data)


#MODELS----

##Step 1 - Modelling design variables----

head(mod_data);dim(mod_data)
mod_data$Sample_Day_Scale <- scale(mod_data$Day_Sampled)
mod_data$Crop_Age_Scale <- scale(mod_data$Crop_Age_Days)

mod_data$Position <- as.factor(mod_data$Position)



Null_Dom <- glmer(Prop_Pred~ 1 + (1|Field/Survey),data = mod_data,family = binomial)

Position <- glmer(Prop_Pred~ Position + (1|Field/Survey),data = mod_data,family = binomial) 
Age <- glmer(Prop_Pred~ Crop_Age_Days + (1|Field/Survey),data = mod_data,family = binomial)
Day <- glmer(Prop_Pred~ Day_Sampled + (1|Field/Survey),data = mod_data,family = binomial)

Position_Age1 <- glmer(Prop_Pred~ Position + Crop_Age_Days + (1|Field/Survey),data = mod_data,family = binomial) 
Position_Age2 <- glmer(Prop_Pred~ Position * Crop_Age_Days + (1|Field/Survey),data = mod_data,family = binomial) 

Position_Day1 <- glmer(Prop_Pred~ Position + Day_Sampled + (1|Field/Survey),data = mod_data,family = binomial) 
Position_Day2 <- glmer(Prop_Pred~ Position * Day_Sampled + (1|Field/Survey),data = mod_data,family = binomial)

Age_Day1 <- glmer(Prop_Pred~ Crop_Age_Days + Day_Sampled + (1|Field/Survey),data = mod_data,family = binomial) 
Age_Day2 <- glmer(Prop_Pred~ Crop_Age_Scale * Sample_Day_Scale + (1|Field/Survey),data = mod_data,family = binomial) #needed re scaled variables

Position_Day_Age1 <- glmer(Prop_Pred~ Position + Day_Sampled + Crop_Age_Days + (1|Field/Survey),data = mod_data,family = binomial) 
Position_Day_Age2 <- glmer(Prop_Pred~ Position * Day_Sampled + Crop_Age_Days + (1|Field/Survey),data = mod_data,family = binomial) 
Position_Day_Age3 <- glmer(Prop_Pred~ Position * Crop_Age_Scale +  Sample_Day_Scale+ (1|Field/Survey),data = mod_data,family = binomial,glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 10000)))
Position_Day_Age4 <- glmer(Prop_Pred~ Position + Crop_Age_Days  *  Day_Sampled + (1|Field/Survey),data = mod_data,family = binomial, glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 10000))) #Didn't converge

aictab(list(null = Null_Dom,position = Position, age = Age, day = Day, position_Age = Position_Age1, positionxage = Position_Age2,position_day = Position_Day1, positionxday = Position_Day2, position_day_age = Position_Day_Age1, positionxday_age =Position_Day_Age2, Positionxage_day = Position_Day_Age3))


##Step 2 - Modelling environmental variables----
#includes fixed effects from top model of the first step

head(mod_data);dim(mod_data)
mod_data$Field_Area_Scaled <- scale(mod_data$Field_Area_m2)

summary(Position_Day_Age3)

#point-level variables
GC <- glmer(Prop_Pred~ Position * Crop_Age_Scale +  Sample_Day_Scale + GC + (1|Field/Survey),data = mod_data,family = binomial) 
Height <- glmer(Prop_Pred~ Position * Crop_Age_Scale +  Sample_Day_Scale + Height + (1|Field/Survey),data = mod_data,family = binomial)


#site-level variables
Field_Area <- glmer(Prop_Pred~ Position * Crop_Age_Scale +  Sample_Day_Scale + Field_Area_Scaled + (1|Field/Survey),data = mod_data,family = binomial) 
NDVI_Field  <- glmer(Prop_Pred~ Position * Crop_Age_Scale +  Sample_Day_Scale + NDVImean_Field + (1|Field/Survey),data = mod_data,family = binomial,glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 10000)))

#landscape-level variables

Prop_Rip  <- glmer(Prop_Pred~ Position * Crop_Age_Scale +  Sample_Day_Scale + X1km_Rip_Prop + (1|Field/Survey),data = mod_data,family = binomial, glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 10000)))
NDVI_1000  <- glmer(Prop_Pred~ Position * Crop_Age_Days   +  Day_Sampled + NDVIsum_1km + (1|Field/Survey),data = mod_data,family = binomial)

aictab(list(null = Null_Dom, positionxage_day = Position_Day_Age3, GC = GC, height = Height, Field_Area = Field_Area, NDVI_Field = NDVI_Field, Proportion_Rip = Prop_Rip, NDVI_Landscape = NDVI_1000))

summary(Position_Day_Age3)

#Extra checks----

#NDVI check

head(mod_data);dim(mod_data)

NDVI_Check_Data <- mod_data[!(mod_data$Survey %in% c(5,7)), ]

NDVI_Field_2 <- glmer(Prop_Pred ~ Position + Crop_Age_Scale +  Sample_Day_Scale + NDVImean_Field + (1|Field/Survey),data = NDVI_Check_Data,family = binomial)
summary(NDVI_Field_2)
summary(NDVI_Field)



#Sensitivity with missing inverts----

##prep data for sensitivity analysis----

Sensitivity <- data.frame(ID = obs$Site,Morphospecies = obs$Morphospecies)

Sensitivity <- merge(Sensitivity,morpho,by = "Morphospecies")
head(Sensitivity);dim(Sensitivity)


Sensitivity <- Sensitivity %>% dplyr::select(ID, Trophic)
head(Sensitivity);dim(Sensitivity)

#remove tropic levels which aren't predator or herbivore
unique(Sensitivity$Trophic)

length(Sensitivity$Trophic[Sensitivity$Trophic == 'Unknown'])
Sensitivity <- subset(Sensitivity, Trophic != "Unknown")

length(Sensitivity$Trophic[Sensitivity$Trophic == 'Ominvore'])
Sensitivity <- subset(Sensitivity, Trophic != "Ominvore")

length(Sensitivity$Trophic[Sensitivity$Trophic == 'Pollinator'])
Sensitivity <- subset(Sensitivity, Trophic != "Pollinator")

length(Sensitivity$Trophic[Sensitivity$Trophic == 'Fungivore'])
Sensitivity <- subset(Sensitivity, Trophic != "Fungivore")

length(Sensitivity$Trophic[Sensitivity$Trophic == 'Hematophagous'])
Sensitivity <- subset(Sensitivity, Trophic != "Hematophagous")

head(Sensitivity);dim(Sensitivity)

#Counting herbivores and predators for each site

herb_count_sense <- Sensitivity %>%
  filter(Trophic == "Herbivore") %>%
  group_by(ID) %>%
  summarize(count = n()) 
head(herb_count_sense);dim(herb_count_sense)

herb_count_sense <- data.frame(herb_count_sense)
names(herb_count_sense)[names(herb_count_sense) == "ID"] <- "Site"
str(herb_count_sense)

Sensitivity[Sensitivity$Trophic == 'Parasite',]$Trophic <- "Predator"
pred_count_sense <- Sensitivity %>%
  filter(Trophic == "Predator") %>%
  group_by(ID) %>%
  summarize(count = n()) 
head(pred_count_sense);dim(pred_count_sense)
pred_count_sense <- data.frame(pred_count_sense)
names(pred_count_sense)[names(pred_count_sense) == "ID"] <- "Site"
str(pred_count_sense)


Sensitivity_1 <- merge(variables,herb_count_sense,by = 'Site',all.x = T)
head(Sensitivity_1);dim(Sensitivity_1)
Sensitivity_1$count[is.na(Sensitivity_1$count)] <- 0
names(Sensitivity_1)[names(Sensitivity_1) == "count"] <- "Herb_Count"
str(Sensitivity_1)

Sensitivity_1 <- merge(Sensitivity_1,pred_count_sense,by = 'Site',all.x = T)
head(Sensitivity_1);dim(Sensitivity_1)
Sensitivity_1$count[is.na(Sensitivity_1$count)] <- 0
names(Sensitivity_1)[names(Sensitivity_1) == "count"] <- "Pred_Count"
str(Sensitivity_1)

#Predator dominance

Sensitivity_1$Prop_Pred <- ifelse(Sensitivity_1$Pred_Count == 0, 0,Sensitivity_1$Pred_Count/(Sensitivity_1$Herb_Count+Sensitivity_1$Pred_Count))

head(Sensitivity_1);dim(Sensitivity_1)
summary(Sensitivity_1$Prop_Pred)
range(Sensitivity_1$Prop_Pred)

head(Sensitivity_1);dim(Sensitivity_1)

##models----
Sensitivity_1$Sample_Day_Scale <- scale(Sensitivity_1$Day_Sampled)
Sensitivity_1$Crop_Age_Scale <- scale(Sensitivity_1$Crop_Age_Days)
Sensitivity_1$Position <- as.factor(Sensitivity_1$Position)
Sensitivity_1$Field_Area_Scaled <-scale(Sensitivity_1$Field_Area_m2)
Sensitivity_1$NDVI1km_Scale <-scale(Sensitivity_1$NDVIsum_1km)

NDVIsum_1km

Null_Sense  <- glmer(Prop_Pred~ 1 + (1|Field/Survey),data = Sensitivity_1,family = binomial)

positionxage_day_Sense <- glmer(Prop_Pred~ Position * Crop_Age_Scale +  Sample_Day_Scale+ (1|Field/Survey),data = Sensitivity_1,family = binomial,glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 10000)))

#point-level variables
GC_Sense <- glmer(Prop_Pred~ Position * Crop_Age_Scale +  Sample_Day_Scale + GC + (1|Field/Survey),data = Sensitivity_1,family = binomial, glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 10000))) 
Height_Sense <- glmer(Prop_Pred~ Position * Crop_Age_Scale +  Sample_Day_Scale + Height + (1|Field/Survey),data = Sensitivity_1,family = binomial,glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 10000)))


#site-level variables
Field_Area_Sense <- glmer(Prop_Pred~ Position * Crop_Age_Scale +  Sample_Day_Scale + Field_Area_Scaled + (1|Field/Survey),data = Sensitivity_1,family = binomial, glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 10000))) 
NDVI_Field_Sense  <- glmer(Prop_Pred~ Position * Crop_Age_Scale +  Sample_Day_Scale + NDVImean_Field + (1|Field/Survey),data = Sensitivity_1,family = binomial,glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 10000)))

#landscape-level variables

Prop_Rip_Sense  <- glmer(Prop_Pred~ Position * Crop_Age_Scale +  Sample_Day_Scale + X1km_Rip_Prop + (1|Field/Survey),data = Sensitivity_1,family = binomial, glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 10000)))
NDVI_1000_Sense  <- glmer(Prop_Pred~ Position * Crop_Age_Scale   +  Sample_Day_Scale + NDVI1km_Scale + (1|Field/Survey),data = Sensitivity_1,family = binomial) # didn't converge

aictab(list(Null = Null_Sense, positionxage_day = positionxage_day_Sense, GC = GC_Sense, height = Height_Sense, Field_Area = Field_Area_Sense, NDVI_Field = NDVI_Field_Sense, Proportion_Rip = Prop_Rip_Sense, NDVI_Landscape = NDVI_1000_Sense))



#SPATIAL AUTOCORRELATION and VIF CHECK----


#VIF CHECK----
library("car")

VIF_Res_1 <- vif(Position_Day_Age3) 

##Top ranked model----

#loop needs a dev.new window open to plot while it runs
dev.new(height=10,width=10,dpi=80,pointsize=14,noRStudioGD = T)


field_numbers <- unique(mod_data$ID)


model_residuals <- simulateResiduals(Position_Day_Age3)
spatial_result <- data.frame(
  field = rep(NA, length(field_numbers)),
  statistic = rep(NA, length(field_numbers)),
  p_value = rep(NA, length(field_numbers)),
  method = rep(NA_character_, length(field_numbers)),
  stringsAsFactors = FALSE)

s <- 1

for (f in field_numbers) {
  
  cat("Field", f, "\n") #What field is loop up to
  
  #Extracting specific residuals for individual fields
  field_indices <- which(mod_data$ID == f)
  field_residuals <- model_residuals
  field_residuals$scaledResiduals <- 
    model_residuals$scaledResiduals[field_indices]
  field_residuals$fittedPredictedResponse <- 
    model_residuals$fittedPredictedResponse[field_indices]
  
  # Test spatial autocorrelation using your grid coordinates
  spatial_test <- testSpatialAutocorrelation(field_residuals, x = mod_data$X_Cor[mod_data$ID == f], y = mod_data$Y_Cor[mod_data$ID == f])
  
  
  spatial_result$field [s] <- f
  spatial_result$statistic [s] <- spatial_test$statistic[1] 
  spatial_result$p_value [s] <- spatial_test$p.value
  spatial_result$method [s] <- spatial_test$method
  
  s <- s + 1
  
}

head(spatial_result);dim(spatial_result)
length(unique(mod_data$ID))

#My own notes/explanation of different bits of results
spatial_test$statistic # The test statistic - note the important one here is observed - want this to be close to zero 
spatial_test$p.value # P-value for the test - indicates if there is significant spatial autocorrelation - statistic above indicates the magnitude and the direction (positive or negative)
spatial_test$method #Method used to get statistic


TopModSpat <- spatial_result
#No Spatial Autocorrelation found in any of the fields/surveys


#PRED AND HERB COUNT MODELS----

#Herbivores
Herb_Mod <- glmmTMB(Herb_Count ~ Position * Crop_Age_Scale +  Sample_Day_Scale + (1|Field/Survey), family = nbinom2, data = mod_data)

Herb_Mod_diagnostic <- simulateResiduals(Herb_Mod, n = 1000)

dev.new(height=10,width=11,dpi=80,pointsize=14,noRStudioGD = T)
testDispersion(Herb_Mod_diagnostic)     #overdispersion

dev.new(height=10,width=11,dpi=80,pointsize=14,noRStudioGD = T)
testZeroInflation(Herb_Mod_diagnostic)  #excess zeros

dev.new(height=10,width=11,dpi=80,pointsize=14,noRStudioGD = T)
testUniformity(Herb_Mod_diagnostic)     #overall distribution fit

#Predators
Pred_Mod <- glmmTMB(Pred_Count ~ Position * Crop_Age_Scale +  Sample_Day_Scale + (1|Field/Survey), family = nbinom2, data = mod_data)

Pred_Mod_diagnostic <- simulateResiduals(Pred_Mod, n = 1000)

dev.new(height=10,width=11,dpi=80,pointsize=14,noRStudioGD = T)
testDispersion(Pred_Mod_diagnostic)     #overdispersion

dev.new(height=10,width=11,dpi=80,pointsize=14,noRStudioGD = T)
testZeroInflation(Pred_Mod_diagnostic)  #excess zeros

dev.new(height=10,width=11,dpi=80,pointsize=14,noRStudioGD = T)
testUniformity(Pred_Mod_diagnostic)    



#PREDICTIONS AND FIGURES----

##Predictions top ranked model----

summary(Position_Day_Age3)

Sample_Day_Pred <- seq(min(mod_data$Sample_Day_Scale),max(mod_data$Sample_Day_Scale),length.out=20) 
Crop_Age_Pred <- seq(min(mod_data$Crop_Age_Scale),max(mod_data$Crop_Age_Scale),length.out=20)
Position_Pred <- c("Outer","Inner")
Position_Pred <- as.factor(Position_Pred)

Design_Pred1 <- expand.grid(Sample_Day_Scale = Sample_Day_Pred, Crop_Age_Scale = Crop_Age_Pred, Position = Position_Pred)

Design_Pred2 <- predictSE(mod = Position_Day_Age3,newdata=Design_Pred1,se.fit = T, type = "link")

Design_Pred3<-data.frame(Design_Pred1,fit.link=Design_Pred2$fit,se.link=Design_Pred2$se.fit)
head(Design_Pred3);dim(Design_Pred3)

Design_Pred3$lci.link<-Design_Pred3$fit.link-(1.96*Design_Pred3$se.link)
Design_Pred3$uci.link<-Design_Pred3$fit.link+(1.96*Design_Pred3$se.link)

Design_Pred3$fit<-invlogit(Design_Pred3$fit.link)
Design_Pred3$se<-invlogit(Design_Pred3$se.link)
Design_Pred3$lci<-invlogit(Design_Pred3$lci.link)
Design_Pred3$uci<-invlogit(Design_Pred3$uci.link)

#checking predctions from conditional vs marginal

library(marginaleffects)

pred_B <- predictions(Position_Day_Age3,
  newdata = datagrid(Sample_Day_Scale = Sample_Day_Pred,
                     Crop_Age_Scale = Crop_Age_Pred[10],
                     Position = "Inner"),
  re.form = NA, type = "response")
head(pred_B)

pred_C <- avg_predictions(Position_Day_Age3, variables = list(Sample_Day_Scale = Sample_Day_Pred), re.form = NA,type = "response")
head(pred_C)

A_pred <- subset(Design_Pred3,Crop_Age_Scale == Crop_Age_Pred[10] & Position == "Inner")

comp <- data.frame(NDVImean_Field = NDVI_Field_Pred, A_predictSE = A_pred$fit,B_datagrid = pred_B$estimate, C_observed = pred_C$estimate)
comp$diff_AB <- comp$A_predictSE - comp$B_datagrid
comp$diff_BC <- comp$B_datagrid - comp$C_observed

max(abs(comp$diff_AB))   #is ~0 which is what we wanted
max(abs(comp$diff_BC))   #how much the marginal definition matters



B_df <- data.frame(Sample_Day_Scale  = pred_B$Sample_Day_Scale ,
                   est = pred_B$estimate, lo = pred_B$conf.low, hi = pred_B$conf.high,
                   version = "B: conditional (fixed covariates)")
C_df <- data.frame(Sample_Day_Scale  = pred_C$Sample_Day_Scale ,
                   est = pred_C$estimate, lo = pred_C$conf.low, hi = pred_C$conf.high,
                   version = "C: marginal (observed data)")\

library('ggplot2')
ggplot(rbind(B_df, C_df), aes(Sample_Day_Scale , est, colour = version, fill = version)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) +
  geom_line()



#Prediction Herb count model----

summary(Herb_Mod)
class(Herb_Mod)

Herb_Pred2 <- predict(object= Herb_Mod,newdata=Design_Pred1,se.fit = T, type = "link",re.form = NA)

Herb_Pred3<-data.frame(Design_Pred1,fit.link=Herb_Pred2$fit,se.link=Herb_Pred2$se.fit)
head(Herb_Pred3);dim(Herb_Pred3)

Herb_Pred3$lci.link<-Herb_Pred3$fit.link-(1.96*Herb_Pred3$se.link)
Herb_Pred3$uci.link<-Herb_Pred3$fit.link+(1.96*Herb_Pred3$se.link)

Herb_Pred3$fit<-exp(Herb_Pred3$fit.link)
Herb_Pred3$se<-exp(Herb_Pred3$se.link)
Herb_Pred3$lci<-exp(Herb_Pred3$lci.link)
Herb_Pred3$uci<-exp(Herb_Pred3$uci.link)


#Prediction Pred count model----

summary(Pred_Mod)

Pred_Pred2 <- predict(object = Pred_Mod,newdata=Design_Pred1,se.fit = T, type = "link",re.form = NA)

Pred_Pred3<-data.frame(Design_Pred1,fit.link=Pred_Pred2$fit,se.link=Pred_Pred2$se.fit)
head(Pred_Pred3);dim(Pred_Pred3)

Pred_Pred3$lci.link<-Pred_Pred3$fit.link-(1.96*Pred_Pred3$se.link)
Pred_Pred3$uci.link<-Pred_Pred3$fit.link+(1.96*Pred_Pred3$se.link)

Pred_Pred3$fit<-exp(Pred_Pred3$fit.link)
Pred_Pred3$se<-exp(Pred_Pred3$se.link)
Pred_Pred3$lci<-exp(Pred_Pred3$lci.link)
Pred_Pred3$uci<-exp(Pred_Pred3$uci.link)


#Figures----

##Top Model----

#to condense code for plotting
xx <- Design_Pred3$Sample_Day_Scale == Sample_Day_Pred[10]
x <- Design_Pred3$Crop_Age_Scale == Crop_Age_Pred[10] & Design_Pred3$Position == "Outer"


OO <- Pred_Pred3$Sample_Day_Scale == Sample_Day_Pred[10] & Pred_Pred3$Position == "Inner"
O_O <- Pred_Pred3$Sample_Day_Scale == Sample_Day_Pred[10] & Pred_Pred3$Position == "Outer"
O <- Pred_Pred3$Crop_Age_Scale == Crop_Age_Pred[10] & Pred_Pred3$Position == "Outer"

HH <- Herb_Pred3$Sample_Day_Scale == Sample_Day_Pred[10] & Herb_Pred3$Position == "Inner"
H_H <- Herb_Pred3$Sample_Day_Scale == Sample_Day_Pred[10] & Herb_Pred3$Position == "Outer"
H <- Herb_Pred3$Crop_Age_Scale == Crop_Age_Pred[10] & Herb_Pred3$Position == "Outer"


dev.new(height=15,width=10,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(5,4,2,2),mfrow=c(3,2),mgp=c(2.5,1,0),xpd = T)

plot(x = mod_data$Crop_Age_Scale,y = mod_data$Prop_Pred,xlab = "Crop Age (Days)",ylab = 'Predator Dominance', ylim = c(0,1),type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,, xaxt = 'n')

polygon(x = c(Design_Pred3$Crop_Age_Scale[xx & Design_Pred3$Position == "Inner"],rev(Design_Pred3$Crop_Age_Scale[xx& Design_Pred3$Position == "Inner"])), y = c(Design_Pred3$lci[xx& Design_Pred3$Position == "Inner"],rev(Design_Pred3$uci[xx& Design_Pred3$Position == "Inner"])),col = rgb(0.7, 0.7, 0.7, 0.7),border = NA)
lines(x=Design_Pred3$Crop_Age_Scale[xx & Design_Pred3$Position == "Inner"],y = Design_Pred3$fit[xx & Design_Pred3$Position == "Inner"],lwd = 2,lty = 1, col = 'grey30')

polygon(x = c(Design_Pred3$Crop_Age_Scale[xx & Design_Pred3$Position == "Outer"],rev(Design_Pred3$Crop_Age_Scale[xx & Design_Pred3$Position == "Outer"])), y = c(Design_Pred3$lci[xx& Design_Pred3$Position == "Outer"],rev(Design_Pred3$uci[xx& Design_Pred3$Position == "Outer"])),col = rgb(0.7, 0.7, 0.7, 0.7),border = NA)
lines(x=Design_Pred3$Crop_Age_Scale[xx & Design_Pred3$Position == "Outer"],y = Design_Pred3$fit[xx & Design_Pred3$Position == "Outer"],lwd = 2,lty = 2, col = 'grey30')
mtext(side=3,line=0,at = -2,'a)',cex=0.8)

axis(side=1, at=seq(from=min(Design_Pred3$Crop_Age_Scale),to=max(Design_Pred3$Crop_Age_Scale),length.out=5),labels=round(seq(from=min(mod_data$Crop_Age_Days),to=max(mod_data$Crop_Age_Days),length.out=5),-1))
mtext("100", side =1, line =1, at =1.3, cex = 0.7)
legend(x = mean(par("usr")[1:2]) + 0.1, 
       y = par("usr")[4] + diff(par("usr")[3:4]) * 0.23, 
       legend = c("Interior", "Edge"), lty = c(1,2), 
       col = "grey30", xpd = NA, bty = "n",
       horiz = TRUE, xjust = 0.4)


plot(x =mod_data$Sample_Day_Scale,y=mod_data$Prop_Pred,xlab = "Day Surveryed",ylab = 'Predator Dominance', ylim = c(0,1),type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2, xaxt = 'n')
mtext(side=1,line=3.3,at = -1.5,'Autumn/Winter',cex=0.6)
mtext(side=1,line=3.3,at = 1.5,'Spring',cex=0.6)
arrows(-0.4,-0.65,1,-0.65, length =0.1)
mtext(side=3,line=0,at = -1.8,'b)',cex=0.8)

polygon(x = c(Design_Pred3$Sample_Day_Scale[x],rev(Design_Pred3$Sample_Day_Scale[x])), y = c(Design_Pred3$lci[x],rev(Design_Pred3$uci[x])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=Design_Pred3$Sample_Day_Scale[x],y = Design_Pred3$fit[x],lwd = 2,col = 'grey30')

axis(side=1, at=seq(from=min(Design_Pred3$Sample_Day_Scale),to=max(Design_Pred3$Sample_Day_Scale),length.out=5),labels=round(seq(from=min(mod_data$Day_Sampled),to=max(mod_data$Day_Sampled),length.out=5),-1))

mtext("120", side =1, line =1, at =0.8, cex = 0.7)


###Predator plots----

plot(x = mod_data$Crop_Age_Scale[mod_data$Pred_Count < 6],y = mod_data$Pred_Count[mod_data$Pred_Count < 6],xlab = "Crop Age (Days)",ylab = 'Predator Abundance', ylim = c(0,6),type = 'p', pch = 16,cex =0.1,col = 'black', las = 1, lwd = 2,, xaxt = 'n')
polygon(x = c(Pred_Pred3$Crop_Age_Scale[OO],rev(Pred_Pred3$Crop_Age_Scale[OO])), y = c(Pred_Pred3$lci[OO],rev(Pred_Pred3$uci[OO])),col = rgb(0.7, 0.7, 0.7, 0.7),border = NA)
lines(x=Pred_Pred3$Crop_Age_Scale[OO],y = Pred_Pred3$fit[OO],lwd = 2,lty = 1, col = 'grey30')
axis(side=1, at=seq(from=min(Pred_Pred3$Crop_Age_Scale),to=max(Pred_Pred3$Crop_Age_Scale),length.out=5),labels=round(seq(from=min(mod_data$Crop_Age_Days),to=max(mod_data$Crop_Age_Days),length.out=5),-1),cex.axis = 0.85)
mtext("c)", side =3, line =0.5, at =-2, cex = 0.8)

polygon(x = c(Pred_Pred3$Crop_Age_Scale[O_O],rev(Pred_Pred3$Crop_Age_Scale[O_O])), y = c(Pred_Pred3$lci[O_O],rev(Pred_Pred3$uci[O_O])),col = rgb(0.7, 0.7, 0.7, 0.7),border = NA)
lines(x=Pred_Pred3$Crop_Age_Scale[O_O],y = Pred_Pred3$fit[O_O],lwd = 2,lty = 2, col = 'grey30')
legend(x = mean(par("usr")[1:2]) + 0.1, 
       y = par("usr")[4] + diff(par("usr")[3:4]) * 0.23, 
       legend = c("Interior", "Edge"), lty = c(1,2), 
       col = "grey30", xpd = NA, bty = "n",
       horiz = TRUE, xjust = 0.4)


plot(x =mod_data$Sample_Day_Scale[mod_data$Pred_Count < 10],y=mod_data$Pred_Count[mod_data$Pred_Count < 10],xlab = "Day Surveyed",ylab = 'Predator Abundance', ylim = c(0,10),type = 'p', pch = 16,cex =0.1,col = 'black', las = 1, lwd = 2, xaxt = 'n')
polygon(x = c(Pred_Pred3$Sample_Day_Scale[O],rev(Pred_Pred3$Sample_Day_Scale[O])), y = c(Pred_Pred3$lci[O],rev(Pred_Pred3$uci[O])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=Pred_Pred3$Sample_Day_Scale[O],y = Pred_Pred3$fit[O],lwd = 2,col = 'grey30')
axis(side=1, at=seq(from=min(Pred_Pred3$Sample_Day_Scale),to=max(Pred_Pred3$Sample_Day_Scale),length.out=5),labels=round(seq(from=min(mod_data$Day_Sampled),to=max(mod_data$Day_Sampled),length.out=5),-1),cex.axis = 0.85)
mtext("160", side =1, line =1, at =1.5, cex = 0.6)
mtext("d)", side =3, line =0.5, at =-1.8, cex = 0.8)
mtext(side=1,line=3.2,at = -1.5,'Autumn/Winter',cex=0.6)
mtext(side=1,line=3.2,at = 1.5,'Spring',cex=0.6)
arrows(-0.40,-6.5,1,-6.5, length =0.1)



###Herbivore plots----
plot(x = mod_data$Crop_Age_Scale[mod_data$Herb_Count < 8],y = mod_data$Herb_Count[mod_data$Herb_Count < 8],xlab = "Crop Age (Days)",ylab = 'Herbivore Abundance', ylim = c(0,8),type = 'p', pch = 16,cex =0.1,col = 'black', las = 1, lwd = 2,, xaxt = 'n')
polygon(x = c(Herb_Pred3$Crop_Age_Scale[HH],rev(Herb_Pred3$Crop_Age_Scale[HH])), y = c(Herb_Pred3$lci[HH],rev(Herb_Pred3$uci[HH])),col = rgb(0.7, 0.7, 0.7, 0.7),border = NA)
lines(x=Herb_Pred3$Crop_Age_Scale[HH],y = Herb_Pred3$fit[HH],lwd = 2,lty = 1, col = 'grey30')
axis(side=1, at=seq(from=min(Herb_Pred3$Crop_Age_Scale),to=max(Herb_Pred3$Crop_Age_Scale),length.out=5),labels=round(seq(from=min(mod_data$Crop_Age_Days),to=max(mod_data$Crop_Age_Days),length.out=5),-1),cex.axis = 0.85)
mtext("f)", side =3, line =0.5, at =-2, cex = 0.8)

polygon(x = c(Herb_Pred3$Crop_Age_Scale[H_H],rev(Herb_Pred3$Crop_Age_Scale[H_H])), y = c(Herb_Pred3$lci[H_H],rev(Pred_Pred3$uci[H_H])),col = rgb(0.7, 0.7, 0.7, 0.7),border = NA)
lines(x=Herb_Pred3$Crop_Age_Scale[H_H],y = Herb_Pred3$fit[H_H],lwd = 2,lty = 2, col = 'grey30')
legend(x = mean(par("usr")[1:2]) + 0.1, 
       y = par("usr")[4] + diff(par("usr")[3:4]) * 0.23, 
       legend = c("Interior", "Edge"), lty = c(1,2), 
       col = "grey30", xpd = NA, bty = "n",
       horiz = TRUE, xjust = 0.4)


plot(x =mod_data$Sample_Day_Scale[mod_data$Herb_Count < 30],y=mod_data$Herb_Count[mod_data$Herb_Count < 30],xlab = "Day Surveyed",ylab = 'Herbivore Abundance', ylim = c(0,30),type = 'p', pch = 16,cex =0.1,col = 'black', las = 1, lwd = 2, xaxt = 'n')
polygon(x = c(Herb_Pred3$Sample_Day_Scale[H],rev(Herb_Pred3$Sample_Day_Scale[H])), y = c(Herb_Pred3$lci[H],rev(Herb_Pred3$uci[H])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=Herb_Pred3$Sample_Day_Scale[H],y = Herb_Pred3$fit[H],lwd = 2,col = 'grey30')
axis(side=1, at=seq(from=min(Herb_Pred3$Sample_Day_Scale),to=max(Herb_Pred3$Sample_Day_Scale),length.out=5),labels=round(seq(from=min(mod_data$Day_Sampled),to=max(mod_data$Day_Sampled),length.out=5),-1),cex.axis = 0.85)
mtext("160", side =1, line =1, at =1.5, cex = 0.6)
mtext("g)", side =3, line =0.5, at =-1.8, cex = 0.8)
mtext(side=1,line=3.2,at = -1.5,'Autumn/Winter',cex=0.6)
mtext(side=1,line=3.2,at = 1.5,'Spring',cex=0.6)
arrows(-0.40,-19.5,1,-19.5, length =0.1)




##Model Diagnostics----
dev.new(height=10,width=15,dpi=80,pointsize=14,noRStudioGD = T)
par(mar = c(4, 4, 3, 2), mfrow = c(2, 3), mgp = c(2.5, 1, 0), cex.lab = 0.8)

testDispersion(Herb_Mod_diagnostic)
mtext("a)", side = 3, line = 0, adj = -0.2, cex = 0.8)

testZeroInflation(Herb_Mod_diagnostic)     
mtext("b)", side = 3, line = 0, adj = -0.25, cex = 0.8)

plotQQunif(Herb_Mod_diagnostic)
mtext("c)", side = 3, line = 0, adj = -0.25, cex = 0.8)

testDispersion(Pred_Mod_diagnostic)  
mtext("d)", side = 3, line = 0, adj = -0.2, cex = 0.8)

testZeroInflation(Pred_Mod_diagnostic)    
mtext("e)", side = 3, line = 0, adj = -0.25, cex = 0.8)

plotQQunif(Pred_Mod_diagnostic)
mtext("f)", side = 3, line = 0, adj = -0.25, cex = 0.8)

##Forest plots----
summary(NDVI_Field)
summary(NDVI_1000)
summary(Height)
summary(Field_Area)
summary(GC)


term_map <- c(
  "(Intercept)" = "Intercept", "NDVImean_Field" = "Field NDVI",
  "Crop_Age_Scale" = "Age", "Sample_Day_Scale" = "Day", 
  "PositionOuter"= "Position (Edge)",
  "PositionOuter:Crop_Age_Scale" ="Position (Edge):Age", 
  "NDVIsum_1km" = "NDVI Sum", 
  "Field_Area_Scaled" = "Field Size", "GC" = "Ground Cover",
  "Crop_Age_Days" = "Age", "Day_Sampled" = "Day",
  "PositionOuter:Crop_Age_Days" = "Position (Edge):Age")


Mods_list <- list(NDVI_Field = NDVI_Field, NDVI_1000 = NDVI_1000, Height=Height, Field_Area= Field_Area, GC=GC)

Mods_list[1]

summary(Mods_list[1]$NDVI_Field)$coefficients
str(Mods_list[1])
names(Mods_list)

Whichmod <- names(Mods_list)
Forest <- vector(mode = "list", length = length(Whichmod))
q <- 1

for (x in Whichmod) {
  
  
  CF <- summary(Mods_list[[x]])$coefficients
  CF <- data.frame(term = row.names(CF),CF)
  row.names(CF) <- 1:nrow(CF)
  colnames(CF) <- c("term","est",'std_err',"z","p")
  
  CF$lci <- CF$est-(1.96*CF$std_err)
  CF$uci <- CF$est+(1.96*CF$std_err)
  
  Forest[[q]] <- CF
  
  q <- q+1
}

Forest[[2]]

Forest <- lapply(Forest, function(df) {
  df$term <- ifelse(df$term %in% names(term_map), term_map[df$term], df$term)
  df
})


dev.new(height=8,width=15,dpi=60,pointsize=14,noRStudioGD = T)
par(mar=c(6,10,1,1),mfrow=c(2,3),mgp=c(2.5,1,0),xpd = F)

plot(x=rev(Forest[[1]]$est),y=1:length(Forest[[1]]$term),ylab="",xlab="Model Estimate",pch=16,yaxt="n",xlim = c(min(Forest[[1]]$lci),max(Forest[[1]]$uci)),cex.axis=0.9,cex.lab=1.2)
axis(2, at=1:length(Forest[[1]]$term),labels = rev(Forest[[1]]$term),las=1,cex.axis=1.1)
arrows(x0 = 0,y0 = 0.5,x1 = 0,y1 = 10,length = 0)
arrows(x0 = rev(Forest[[1]]$lci),y0 = 1:length(Forest[[1]]$term),x1 = rev(Forest[[1]]$uci),y1 = 1:length(Forest[[1]]$term),length = 0)
mtext("a)",side = 3,line = -0.5,at=-16,cex = 0.8)

plot(x=rev(Forest[[2]]$est),y=1:length(Forest[[2]]$term),ylab="",xlab="Model Estimate",pch=16,yaxt="n",xlim = c(min(Forest[[2]]$lci),max(Forest[[2]]$uci)),cex.axis=0.9,cex.lab=1.2)
axis(2, at=1:length(Forest[[2]]$term),labels = rev(Forest[[2]]$term),las=1,cex.axis=1.1)
arrows(x0 = 0,y0 = 0.5,x1 = 0,y1 = 10,length = 0)
arrows(x0 = rev(Forest[[2]]$lci),y0 = 1:length(Forest[[2]]$term),x1 = rev(Forest[[2]]$uci),y1 = 1:length(Forest[[2]]$term),length = 0)
mtext("b)",side = 3,line = -0.5,at=-16,cex = 0.8)

plot(x=rev(Forest[[3]]$est),y=1:length(Forest[[3]]$term),ylab="",xlab="Model Estimate",pch=16,yaxt="n",xlim = c(min(Forest[[3]]$lci),max(Forest[[3]]$uci)),cex.axis=0.9,cex.lab=1.2)
axis(2, at=1:length(Forest[[3]]$term),labels = rev(Forest[[3]]$term),las=1,cex.axis=1.1)
arrows(x0 = 0,y0 = 0.5,x1 = 0,y1 = 10,length = 0)
arrows(x0 = rev(Forest[[3]]$lci),y0 = 1:length(Forest[[3]]$term),x1 = rev(Forest[[3]]$uci),y1 = 1:length(Forest[[3]]$term),length = 0)
mtext("c)",side = 3,line = -0.5,at=-4.5,cex = 0.8)

plot(x=rev(Forest[[4]]$est),y=1:length(Forest[[4]]$term),ylab="",xlab="Model Estimate",pch=16,yaxt="n",xlim = c(min(Forest[[4]]$lci),max(Forest[[4]]$uci)),cex.axis=0.8,cex.lab=1.2)
axis(2, at=1:length(Forest[[4]]$term),labels = rev(Forest[[4]]$term),las=1,cex.axis=1.1)
arrows(x0 = 0,y0 = 0.5,x1 = 0,y1 = 10,length = 0)
arrows(x0 = rev(Forest[[4]]$lci),y0 = 1:length(Forest[[4]]$term),x1 = rev(Forest[[4]]$uci),y1 = 1:length(Forest[[4]]$term),length = 0)
mtext("d)",side = 3,line = -0.5,at=-3.2,cex = 0.8)

plot(x=rev(Forest[[5]]$est),y=1:length(Forest[[5]]$term),ylab="",xlab="Model Estimate",pch=16,yaxt="n",xlim = c(min(Forest[[5]]$lci),max(Forest[[5]]$uci)),cex.axis=0.9,cex.lab=1.2)
axis(2, at=1:length(Forest[[2]]$term),labels = rev(Forest[[5]]$term),las=1,cex.axis=1.1)
arrows(x0 = 0,y0 = 0.5,x1 = 0,y1 = 10,length = 0)
arrows(x0 = rev(Forest[[5]]$lci),y0 = 1:length(Forest[[5]]$term),x1 = rev(Forest[[5]]$uci),y1 = 1:length(Forest[[5]]$term),length = 0)
mtext("e)",side = 3,line = -0.5,at=-5.3,cex = 0.8)


#END----

