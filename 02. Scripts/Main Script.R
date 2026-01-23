
#Main Script for data analysis for Predator Dominance Analysis  

#PACKAGES----
library("corrplot")
library("dplyr")
library('lme4')
library("AICcmodavg")
library('arm')
library("nlme")
library("DHARMa")

#Load function:

invisible(lapply(paste("Functions/", dir("Functions"), sep=""), function(x) source(x)))


#DATA----

#All location relevent information has been removed from data to protect privacy

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

table(field$X1km_Prop_Water)
table(field$X1km_Prop_Crops)
table(field$X1km_Prop_Nat_Graze)

hist(point$Ground_Cover)
hist(point$Plant_Height)

#CORRELATION----

cordata2 <- data.frame(ID = point$Survey_Field,Height = point$Plant_Height,GC = point$Ground_Cover,Position = point$Spatial_Position)

cordata2 <- merge(cordata2,field, by = "ID")

cordata2 <- cordata2 %>% dplyr::select(Height, GC, Position,Day_Sampled,Crop_Age_Days,Field_Area_m2,X500m_Prop_Crops,X500m_Prop_Water,X1km_Prop_Crops,X1km_Prop_Water,X1km_Prop_Nat_Graze,NDVImean_Field,NDVIsum_500m,NDVImean_500m,NDVIsum_1km,NDVImean_1km)

cordata2$Position[cordata2$Position == 'Outer'] <- 1
cordata2$Position[cordata2$Position == 'Inner'] <- 2
cordata2$Position <- as.numeric(cordata2$Position)

str(cordata2)
cor2 <- cor(cordata2,method = "spearman")
colnames(cor2) <- c("Height", "GC", "Position","Day","Crop Age","Field Area","Crops 500m","Water 500m","Crops 1000m","Water 1000m","Graze 1000m","NDVI Field Mean","NDVI 500m Sum","NDVI 500m Mean","NDVI 1000m Sum","NDVI 1000m Mean")
rownames(cor2) <- c("Height", "GC", "Position","Day","Crop Age","Field Area","Crops 500m","Water 500m","Crops 1000m","Water 1000m","Graze 1000m","NDVI Field Mean","NDVI 500m Sum","NDVI 500m Mean","NDVI 1000m Sum","NDVI 1000m Mean")

dev.new(height=8,width=10,dpi=80,pointsize=14,noRStudioGD = T)
corrplot::corrplot(cor2,method="color",  
                   type="upper",addCoef.col = 'black',number.cex = 0.6)
head(cordata2)

#Uncorrelated variables listed below based on correlation plot

#Design variables: Crop age, day sampled and spatial position
#Point Variables: Ground cover, plant height
#Field variables: Field size, field NDVI (mean)
#Landscape variable: water 1km, NDVI (sum) 1km


#DARA PREP----

##prep variables----
variables <- data.frame(Site = point$Site,Height = point$Plant_Height,GC = point$Ground_Cover,Position = point$Spatial_Position,ID = point$Survey_Field)

variables <- merge(variables,field, by = "ID")

head(variables);dim(variables)

variables <- variables %>% dplyr::select(ID,Site,Field,Height, GC, Position,Day_Sampled,Crop_Age_Days,Field_Area_m2,X1km_Prop_Crops,X1km_Prop_Water,NDVImean_Field,NDVIsum_500m,NDVIsum_1km)


##prep count data----
head(obs);dim(obs)
head(morpho);dim(morpho)

invert <- data.frame(ID = obs$Site,Morphospecies = obs$Morphospecies)

invert <- merge(invert,morpho,by = "Morphospecies")
head(invert);dim(invert)

#Removing wasps, lepidoptera and Serpintine Leaf Miner to control for chemicals applied during some of the sampling period

length(invert$Trophic[invert$Trophic == 'Parasite'])
invert <- invert[invert$Trophic != 'Parasite',]

table(invert$Order)
length(invert$Trophic[invert$Order == 'Lepidoptera'])
invert <- invert[invert$Trophic != 'Lepidoptera',]

length(invert$Morphospecies[invert$Morphospecies == 'Serpintine_Leaf_Miner'])
invert <- invert[invert$Morphospecies != 'Serpintine_Leaf_Miner',]


invert <- invert %>% dplyr::select(ID, Trophic)
head(invert);dim(invert)

#remove tropic levels which aren't predator or herbivore
unique(invert$Trophic)

length(invert$Trophic[invert$Trophic == 'Unknown'])
invert <- invert[invert$Trophic != 'Unknown',]

length(invert$Trophic[invert$Trophic == 'Ominvore'])
invert <- invert[invert$Trophic != 'Ominvore',]

length(invert$Trophic[invert$Trophic == 'Pollinator'])
invert <- invert[invert$Trophic != 'Pollinator',]

length(invert$Trophic[invert$Trophic == 'Fungivore'])
invert <- invert[invert$Trophic != 'Fungivore',]

length(invert$Trophic[invert$Trophic == 'Hematophagous'])
invert <- invert[invert$Trophic != 'Hematophagous',]

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


#MODELLING----

##Step 1 - Modelling design variables----

head(mod_data);dim(mod_data)
mod_data$Sample_Day_Scale <- scale(mod_data$Day_Sampled)
mod_data$Crop_Age_Scale <- scale(mod_data$Crop_Age_Days)

mod_data$Position <- as.factor(mod_data$Position)



Null_Dom <- glmer(Prop_Pred~ 1 + (1|Field),data = mod_data,family = binomial)

Position <- glmer(Prop_Pred~ Position + (1|Field),data = mod_data,family = binomial) 
Age <- glmer(Prop_Pred~ Crop_Age_Days + (1|Field),data = mod_data,family = binomial)
Day <- glmer(Prop_Pred~ Day_Sampled + (1|Field),data = mod_data,family = binomial)

Position_Age1 <- glmer(Prop_Pred~ Position + Crop_Age_Days + (1|Field),data = mod_data,family = binomial) 
Position_Age2 <- glmer(Prop_Pred~ Position * Crop_Age_Days + (1|Field),data = mod_data,family = binomial) 

Position_Day1 <- glmer(Prop_Pred~ Position + Day_Sampled + (1|Field),data = mod_data,family = binomial) 
Position_Day2 <- glmer(Prop_Pred~ Position * Day_Sampled + (1|Field),data = mod_data,family = binomial)

Age_Day1 <- glmer(Prop_Pred~ Crop_Age_Days + Day_Sampled + (1|Field),data = mod_data,family = binomial) 
Age_Day2 <- glmer(Prop_Pred~ Crop_Age_Scale * Sample_Day_Scale + (1|Field),data = mod_data,family = binomial) #needed re scaled variables

Position_Day_Age1 <- glmer(Prop_Pred~ Position + Day_Sampled + Crop_Age_Days + (1|Field),data = mod_data,family = binomial) 
Position_Day_Age2 <- glmer(Prop_Pred~ Position * Day_Sampled + Crop_Age_Days + (1|Field),data = mod_data,family = binomial) 
Position_Day_Age3 <- glmer(Prop_Pred~ Position * Crop_Age_Scale +  Sample_Day_Scale+ (1|Field),data = mod_data,family = binomial) #needed re scaled variables
Position_Day_Age4 <- glmer(Prop_Pred~ Position + Crop_Age_Scale  *  Sample_Day_Scale + (1|Field),data = mod_data,family = binomial, glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 10000))) #increased the number of iterations for convergence
Position_Day_Age5 <- glmer(Prop_Pred~ Position * Crop_Age_Scale *  Sample_Day_Scale + (1|Field),data = mod_data,family = binomial, glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 10000))) #increased the number of iterations for convergence

aictab(list(null = Null_Dom,position = Position, age = Age, day = Day, position_Age = Position_Age1, positionxage = Position_Age2,position_day = Position_Day1, positionxday = Position_Day2, position_day_age = Position_Day_Age1, positionxday_age =Position_Day_Age2, positionxage_day = Position_Day_Age3, position_agexday = Position_Day_Age4,positionxagexday = Position_Day_Age5))


##Step 2 - Modelling environmental variables----
#includes fixed effects from top model of the first step

head(mod_data);dim(mod_data)
mod_data$Field_Area_Scaled <- scale(mod_data$Field_Area_m2)

summary(Position_Day_Age3)

#point-level variables
GC <- glmer(Prop_Pred~ Position * Crop_Age_Scale +  Sample_Day_Scale + GC + (1|Field),data = mod_data,family = binomial, glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 10000))) #increased the number of iterations for convergence
Height <- glmer(Prop_Pred~ Position * Crop_Age_Scale +  Sample_Day_Scale + Height + (1|Field),data = mod_data,family = binomial)


#site-level variables
Field_Area <- glmer(Prop_Pred~ Position * Crop_Age_Scale +  Sample_Day_Scale + Field_Area_Scaled + (1|Field),data = mod_data,family = binomial, glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 10000))) #increased the number of iterations for convergence
NDVI_Field  <- glmer(Prop_Pred~ Position * Crop_Age_Scale +  Sample_Day_Scale + NDVImean_Field + (1|Field),data = mod_data,family = binomial)

#landscape-level variables

Prop_Water  <- glmer(Prop_Pred~ Position * Crop_Age_Scale +  Sample_Day_Scale + X1km_Prop_Water + (1|Field),data = mod_data,family = binomial, glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 10000))) #increased the number of iterations for convergence
NDVI_1000  <- glmer(Prop_Pred~ Position * Crop_Age_Days   +  Day_Sampled + NDVIsum_1km + (1|Field),data = mod_data,family = binomial)

aictab(list(null = Null_Dom, positionxage_day = Position_Day_Age3, GC = GC, height = Height, Field_Area = Field_Area, NDVI_Field = NDVI_Field, Proportion_Water = Prop_Water,NDVI_Landscape_1000 = NDVI_1000))

#SPATIAL AUTOCORRELATION CHECK----

##Top ranked model----

#loop needs a dev.new window open to plot while it runs
dev.new(height=10,width=10,dpi=80,pointsize=14,noRStudioGD = T)


field_numbers <- unique(mod_data$ID)


model_residuals <- simulateResiduals(NDVI_Field)
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


##Second ranked model----

#Remember to open a dev.new window so loop will run
dev.new(height=10,width=10,dpi=80,pointsize=14,noRStudioGD = T)


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


SecondModSpat <- spatial_result
#No Spatial Autocorrelation found in any of the fields/surveys


#PREDICTIONS AND FIGURES----

##Predictions top ranked model (Figure 3)----

summary(NDVI_Field)

NDVI_Field_Pred <- seq(min(mod_data$NDVImean_Field ),max(mod_data$NDVImean_Field ),length.out=20)


Enviro_Pred1 <- expand.grid(Sample_Day_Scale = Sample_Day_Pred, Crop_Age_Scale = Crop_Age_Pred, Position = Position_Pred, NDVImean_Field = NDVI_Field_Pred)

Enviro_Pred2 <- predictSE(mod = NDVI_Field,newdata=Enviro_Pred1,se.fit = T, type = "link")

Enviro_Pred3<-data.frame(Enviro_Pred1,fit.link=Enviro_Pred2$fit,se.link=Enviro_Pred2$se.fit)
head(Enviro_Pred3);dim(Enviro_Pred3)

Enviro_Pred3$lci.link<-Enviro_Pred3$fit.link-(1.96*Enviro_Pred3$se.link)
Enviro_Pred3$uci.link<-Enviro_Pred3$fit.link+(1.96*Enviro_Pred3$se.link)

Enviro_Pred3$fit<-invlogit(Enviro_Pred3$fit.link)
Enviro_Pred3$se<-invlogit(Enviro_Pred3$se.link)
Enviro_Pred3$lci<-invlogit(Enviro_Pred3$lci.link)
Enviro_Pred3$uci<-invlogit(Enviro_Pred3$uci.link)

###Figure 3 ----

yy <- Enviro_Pred3$Sample_Day_Scale == Sample_Day_Pred[10]
y <- Enviro_Pred3$Crop_Age_Scale == Crop_Age_Pred[10] & Enviro_Pred3$Position == "Inner" & Enviro_Pred3$NDVImean_Field == NDVI_Field_Pred[10]
yyy<- Enviro_Pred3$Crop_Age_Scale == Crop_Age_Pred[10] & Enviro_Pred3$Position == "Inner" & Enviro_Pred3$Sample_Day_Scale == Sample_Day_Pred[10]

dev.new(height=10,width=10,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,2,2),mfrow=c(2,2),mgp=c(2.5,1,0),xpd = T)

plot(x = mod_data$Crop_Age_Scale,y = mod_data$Prop_Pred,xlab = "Crop Age (Days)",ylab = 'Predator Dominance', ylim = c(0,1),type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,, xaxt = 'n')

polygon(x = c(Enviro_Pred3$Crop_Age_Scale[yy & Enviro_Pred3$Position == "Inner" & Enviro_Pred3$NDVImean_Field == NDVI_Field_Pred[10]],rev(Enviro_Pred3$Crop_Age_Scale[yy& Enviro_Pred3$Position == "Inner" & Enviro_Pred3$NDVImean_Field == NDVI_Field_Pred[10]])), y = c(Enviro_Pred3$lci[yy& Enviro_Pred3$Position == "Inner"& Enviro_Pred3$NDVImean_Field == NDVI_Field_Pred[10]],rev(Enviro_Pred3$uci[yy& Enviro_Pred3$Position == "Inner"& Enviro_Pred3$NDVImean_Field == NDVI_Field_Pred[10]])),col = rgb(0.7, 0.7, 0.7, 0.7),border = NA)
lines(x=Enviro_Pred3$Crop_Age_Scale[yy & Enviro_Pred3$Position == "Inner" & Enviro_Pred3$NDVImean_Field == NDVI_Field_Pred[10]],y = Enviro_Pred3$fit[yy & Enviro_Pred3$Position == "Inner"& Enviro_Pred3$NDVImean_Field == NDVI_Field_Pred[10]],lwd = 2,lty = 1, col = 'grey30')
axis(side=1, at=seq(from=min(Enviro_Pred3$Crop_Age_Scale),to=max(Enviro_Pred3$Crop_Age_Scale),length.out=6),labels=round(seq(from=min(mod_data$Crop_Age_Days),to=max(mod_data$Crop_Age_Days),length.out=6),-1),cex.axis = 1)
mtext("100", side =1, line =1, at =0.7, cex = 0.9)
mtext("a)", side =3, line =0.5, at =-2.75, cex = 0.9)

polygon(x = c(Enviro_Pred3$Crop_Age_Scale[yy & Enviro_Pred3$Position == "Outer" & Enviro_Pred3$NDVImean_Field == NDVI_Field_Pred[10]],rev(Enviro_Pred3$Crop_Age_Scale[yy& Enviro_Pred3$Position == "Outer" & Enviro_Pred3$NDVImean_Field == NDVI_Field_Pred[10]])), y = c(Enviro_Pred3$lci[yy& Enviro_Pred3$Position == "Outer"& Enviro_Pred3$NDVImean_Field == NDVI_Field_Pred[10]],rev(Enviro_Pred3$uci[yy& Enviro_Pred3$Position == "Outer"& Enviro_Pred3$NDVImean_Field == NDVI_Field_Pred[10]])),col = rgb(0.7, 0.7, 0.7, 0.7),border = NA)
lines(x=Enviro_Pred3$Crop_Age_Scale[yy & Enviro_Pred3$Position == "Outer" & Enviro_Pred3$NDVImean_Field == NDVI_Field_Pred[10]],y = Enviro_Pred3$fit[yy & Enviro_Pred3$Position == "Outer"& Enviro_Pred3$NDVImean_Field == NDVI_Field_Pred[10]],lwd = 2,lty = 2, col = 'grey30')
legend('topleft',legend = c('Interior', "Edge"), lty = c(1,2), col = 'grey30',pt.cex = 1)


plot(x =mod_data$Sample_Day_Scale,y=mod_data$Prop_Pred,xlab = "Day Surveyed",ylab = 'Predator Dominance', ylim = c(0,1),type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2, xaxt = 'n')

polygon(x = c(Enviro_Pred3$Sample_Day_Scale[y],rev(Enviro_Pred3$Sample_Day_Scale[y])), y = c(Enviro_Pred3$lci[y],rev(Enviro_Pred3$uci[y])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=Enviro_Pred3$Sample_Day_Scale[y],y = Enviro_Pred3$fit[y],lwd = 2,col = 'grey30')
axis(side=1, at=seq(from=min(Enviro_Pred3$Sample_Day_Scale),to=max(Enviro_Pred3$Sample_Day_Scale),length.out=6),labels=round(seq(from=min(mod_data$Day_Sampled),to=max(mod_data$Day_Sampled),length.out=6),-1),cex.axis = 1)
mtext("120", side =1, line =1, at =0.9, cex = 0.9)
mtext("b)", side =3, line =0.5, at =-2.6, cex = 0.9)
mtext(side=1,line=3,at = -1.5,'Autumn/Winter',cex=0.7)
mtext(side=1,line=3,at = 1.5,'Spring',cex=0.7)
arrows(-0.70,-0.43,1.10,-0.43, length =0.1)



plot(x =mod_data$NDVImean_Field,y=mod_data$Prop_Pred,xlab = "Mean NDVI of Field",ylab = 'Predator Dominance', ylim = c(0,1),type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2, xaxt = 'n')

polygon(x = c(Enviro_Pred3$NDVImean_Field[yyy],rev(Enviro_Pred3$NDVImean_Field[yyy])), y = c(Enviro_Pred3$lci[yyy],rev(Enviro_Pred3$uci[yyy])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=Enviro_Pred3$NDVImean_Field[yyy],y = Enviro_Pred3$fit[yyy],lwd = 2,col = 'grey30')
axis(side=1, at=seq(from=min(Enviro_Pred3$NDVImean_Field),to=max(Enviro_Pred3$NDVImean_Field),length.out=6),labels=round(seq(from=min(Enviro_Pred3$NDVImean_Field),to=max(Enviro_Pred3$NDVImean_Field),length.out=6),2),cex.axis = 0.95)
mtext("0.22", side =1, line =1, at =0.22, cex = 0.8)
mtext("0.35", side =1, line =1, at =0.35, cex = 0.8)
mtext("0.49", side =1, line =1, at =0.49, cex = 0.8)
mtext("c)", side =3, line =0.5, at =0.04, cex = 0.9)


##Predictions second ranked model----

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

##Figure S2 ----
Crop_Age_Pred

xx <- Design_Pred3$Sample_Day_Scale == Sample_Day_Pred[10]
x <- Design_Pred3$Crop_Age_Scale == Crop_Age_Pred[10] & Design_Pred3$Position == "Inner"

head(Design_Pred3,3); dim(Design_Pred3)
head(mod_data,3); dim(mod_data)


dev.new(height=5,width=10,dpi=80,pointsize=14,noRStudioGD = T)
par(mar=c(4,4,2,2),mfrow=c(1,2),mgp=c(2.5,1,0),xpd = T)

plot(x = mod_data$Crop_Age_Scale,y = mod_data$Prop_Pred,xlab = "Crop Age (Days)",ylab = 'Predator Dominance', ylim = c(0,1),type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2,, xaxt = 'n')

polygon(x = c(Design_Pred3$Crop_Age_Scale[xx & Design_Pred3$Position == "Inner"],rev(Design_Pred3$Crop_Age_Scale[xx& Design_Pred3$Position == "Inner"])), y = c(Design_Pred3$lci[xx& Design_Pred3$Position == "Inner"],rev(Design_Pred3$uci[xx& Design_Pred3$Position == "Inner"])),col = rgb(0.7, 0.7, 0.7, 0.7),border = NA)
lines(x=Design_Pred3$Crop_Age_Scale[xx & Design_Pred3$Position == "Inner"],y = Design_Pred3$fit[xx & Design_Pred3$Position == "Inner"],lwd = 2,lty = 1, col = 'grey30')

polygon(x = c(Design_Pred3$Crop_Age_Scale[xx & Design_Pred3$Position == "Outer"],rev(Design_Pred3$Crop_Age_Scale[xx & Design_Pred3$Position == "Outer"])), y = c(Design_Pred3$lci[xx& Design_Pred3$Position == "Outer"],rev(Design_Pred3$uci[xx& Design_Pred3$Position == "Outer"])),col = rgb(0.7, 0.7, 0.7, 0.7),border = NA)
lines(x=Design_Pred3$Crop_Age_Scale[xx & Design_Pred3$Position == "Outer"],y = Design_Pred3$fit[xx & Design_Pred3$Position == "Outer"],lwd = 2,lty = 2, col = 'grey30')
mtext(side=3,line=0,at = -2.5,'a)',cex=1)

axis(side=1, at=seq(from=min(Design_Pred3$Crop_Age_Scale),to=max(Design_Pred3$Crop_Age_Scale),length.out=6),labels=round(seq(from=min(mod_data$Crop_Age_Days),to=max(mod_data$Crop_Age_Days),length.out=6),-1))

legend('topleft',legend = c('Interior', "Edge"), lty = c(1,2), col = 'grey30',pt.cex = 1)


plot(x =mod_data$Sample_Day_Scale,y=mod_data$Prop_Pred,xlab = "Day Surveryed",ylab = 'Predator Dominance', ylim = c(0,1),type = 'p', pch = 16,cex =0.2,col = 'black', las = 1, lwd = 2, xaxt = 'n')
mtext(side=1,line=3,at = -1.5,'Autumn/Winter',cex=0.8)
mtext(side=1,line=3,at = 1.5,'Spring',cex=0.8)
arrows(-1,-0.3,1.25,-0.3, length =0.1)
mtext(side=3,line=0,at = -2.3,'b)',cex=1)

polygon(x = c(Design_Pred3$Sample_Day_Scale[x],rev(Design_Pred3$Sample_Day_Scale[x])), y = c(Design_Pred3$lci[x],rev(Design_Pred3$uci[x])),col = rgb(0.5, 0.5, 0.5, 0.5),border = NA)
lines(x=Design_Pred3$Sample_Day_Scale[x],y = Design_Pred3$fit[x],lwd = 2,col = 'grey30')

axis(side=1, at=seq(from=min(Design_Pred3$Sample_Day_Scale),to=max(Design_Pred3$Sample_Day_Scale),length.out=6),labels=round(seq(from=min(mod_data$Day_Sampled),to=max(mod_data$Day_Sampled),length.out=6),-1))

#END----

