library(raster)
library(randomForest)
library(rgdal)
library(sp)
library(sf)
library(dismo)
library(reshape2)
library(ggplot2)
library(hyperSpec)
library(data.table)
library(caret)
library(dplyr)
library(tidyr)
library(gridExtra)
library(cowplot)
library(RColorBrewer)
library(leaflet)
library(plyr)
library(pROC)
library(plyr)
library(viridis)
#rm(list=ls()) sometimes I need it for dcast, or in general 
#Set working directory where .csv of all ground truth points (classes) and 
# rasters of UAS imagery 
#setwd("Z:/UAV2019/Alpena/Composites11cm")

#CSV with each row as single ground truth point and each column representing
#each band from UAS imagery. each value represented the extracted band value of
#the pixel that intersects with the points. 
data_dsm <- read.csv("AlpenaGT_11cm_DSM.csv")
head(data_dsm)
#these point values were extracted from UAS orthomosaics (11cm spatial resolution) that we collected 
#flying over a wetland in Alpena, MI during summer 2019 
#train is numeric code of vegetation classes (important for models)
#x and y coordinates of each point in UTM 

#NOTE: I sampled these values in ArcMap, and exported a csv 
#we could also do this in R, but it takes longer 

#Subset dataframe and exclude irrelevant columns for this analysis 
aug_dsm <- subset(data_dsm,
                   select=-c(X.1,Green7_12,Red7_12,Red.Edge7_12,NIR7_12,
                             NDVI7_12,DSM7_12,StdVegHeight7_12))
aug_dsm <- drop_na(aug_dsm)


#check how many points we have per class
nrow(aug_dsm[which(aug_dsm$train==1),]) #64 floating veg. points (class 1)
nrow(aug_dsm[which(aug_dsm$train==2),]) #103 duckweed points (class 2) 
nrow(aug_dsm[which(aug_dsm$train==3),]) #66 submergent veg. points (class 3)
nrow(aug_dsm[which(aug_dsm$train==4),]) #82 emergent veg. points (class 4)


####----spectral profile----#### 
#august 
aug_prof_features <- reshape2::melt(aug_dsm,id=c("train"),
                         measure=c("Green8_5","Red8_5","Red.Edge8_5","NIR8_5","NDVI8_5","DSM8_5","StdVegHeight8_5"))
names(aug_prof_features)[3] <- "Reflectance"

#average 
aug_avg<-reshape2::dcast(aug_prof_features,train+variable~.,fun=mean, value.var="Reflectance",na.rm=TRUE)
names(aug_avg)[3]<-"mean_ref"
#se
aug_se<-reshape2::dcast(aug_prof_features, train+variable~., fun=sd, value.var="Reflectance",na.rm=TRUE)
names(aug_se)[3]<-"sd"
table(aug_dsm$train)#order for below command 
#first reps is train sample n, second is reps are for band #'s (4 classes with varying n, now 7 bands in this aug data)
#numerical order, (efb,yv,sv,emergent)
Class_length<-c(rep(64,7),rep(103,7),rep(66,7),rep(82,7))
as.data.frame(Class_length)
aug_se$class_n<-Class_length
#jenky way, but works  
aug_se$se <- aug_se$sd/sqrt(aug_se$class_n)

aug_spectra <- merge(aug_avg, aug_se, by=c("train","variable"), na.rm=TRUE)
aug_spectra$lowerse <- aug_spectra$mean_ref-aug_spectra$se
aug_spectra$upperse <- aug_spectra$mean_ref+aug_spectra$se

aug_spectra$train <- as.factor(aug_spectra$train)
# aug_all <- aug_spectra %>% mutate(train=as.character(train),
#                                       train= if_else(train == 1,"EFB",train),
#                                       train=as.factor(train))  %>% mutate(train=as.character(train),
#                                                                           train= if_else(train == 2,"Yellow Veg.",train),
#                                                                           train=as.factor(train)) %>% mutate(train=as.character(train),
#                                                                                                              train= if_else(train == 3,"Sub.Veg",train),
#                                                                                                              train=as.factor(train))  %>% mutate(train=as.character(train),
#                                                                                                                                                  train= if_else(train == 4,"Emergent",train),
#                                                                                                                                                  train=as.factor(train))
aug_spectra$train <- if_else(aug_spectra$train == 1, "EFB + Lily spp.", as.character(aug_spectra$train))
aug_spectra$train <- if_else(aug_spectra$train == 2, "Duckweed spp.", as.character(aug_spectra$train))
aug_spectra$train <- if_else(aug_spectra$train == 3, "Sub.Veg.", as.character(aug_spectra$train))
aug_spectra$train <- if_else(aug_spectra$train == 4, "Emergent Veg.", as.character(aug_spectra$train))

august_wlegend <- ggplot(aug_spectra,aes(x=variable,y=mean_ref,group=train,color=train)) +
  geom_line(size=1.25)+
  geom_errorbar(aes(ymin=lowerse,ymax=upperse,width=0.05))+
  ylab("Reflectance Value")+
  xlab("Band Name")+
  ggtitle("Veg. Classes Spectral/Textural Profiles in Alpena, MI, August 2019")+
  theme_bw()+
  theme(axis.title.x = element_text(size = 14, face = "bold")) +
  # scale_color_discrete(breaks=c("Frogbit","Lily","SV","Typha + EFB/Litter"),
  #                      labels=c("Frogbit","Lily","SV","Typha + EFB/Litter")) +
  scale_x_discrete(labels=c("Green","Red","Red Edge","NIR","NDVI","StdVegHeight"))

#Actual spectral plot 
aug_spectra_only <- subset(aug_spectra, variable !="DSM8_5") 
aug_spectra_only <- subset(aug_spectra_only, variable !="StdVegHeight8_5") 

aug_spectra_only$train <- factor(aug_spectra_only$train,
                                 levels = c("EFB + Lily spp.","Duckweed spp.","Sub.Veg.","Emergent Veg."))

aug_spectra_plot <- ggplot(aug_spectra_only,aes(x=variable,y=mean_ref,group=train,color=train)) +
  scale_color_viridis(discrete = TRUE,direction=-1) + 
  geom_line(size=1.25)+
  geom_errorbar(aes(ymin=lowerse,ymax=upperse,width=0.05))+
  ylab("% Reflectance")+
  #ggtitle("Class Spectral Profiles")+
  theme_bw()+
  theme(legend.position='none',
        axis.title.x = element_blank(),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  # scale_color_discrete(breaks=c("Frogbit","Lily","SV","Typha + EFB/Litter"),
  #                      labels=c("Frogbit","Lily","SV","Typha + EFB/Litter")) +
  scale_x_discrete(labels=c("Green","Red","Red Edge","NIR","NDVI"))

#rug only 
aug_rugosity <- subset(aug_spectra, variable == "StdVegHeight8_5")
aug_rugosity$train <- factor(aug_rugosity$train,
                                 levels = c("Duckweed spp.","EFB + Lily spp.","Sub.Veg.","Emergent Veg."))

aug_surf_rug_plot <- ggplot(aug_rugosity,aes(x=variable,y=mean_ref,fill=train)) +
  scale_fill_viridis(discrete = "TRUE",direction=-1,
                     labels = c("Duckweed spp.","EFB + Lily spp.","Emergent Veg.","Sub. Veg.")) + 
  geom_bar(position="dodge",stat="identity") +
  geom_errorbar(aes(ymin=lowerse,ymax=upperse,width=0.2), 
                position=position_dodge(0.9)) +
  ylab("SD of Surface Height (m)")+
  #xlab("Band Name")+
  #ggtitle("DSM Texture")+
  theme_bw()+
  theme(legend.position='none',
        axis.title.x=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  # scale_color_discrete(breaks=c("Frogbit","Lily","SV","Typha + EFB/Litter"),
  #                      labels=c("Frogbit","Lily","SV","Typha + EFB/Litter")) +
  scale_x_discrete(labels=c("Rugosity"))

#surface height only 
aug_surface <- subset(aug_spectra, variable == "DSM8_5")
aug_surface$train <- factor(aug_surface$train,
                             levels = c("Duckweed spp.","EFB + Lily spp.","Sub.Veg.","Emergent Veg."))

aug_surf_plot <- ggplot(aug_surface,aes(x=variable,y=mean_ref,fill=train)) +
  scale_fill_viridis(discrete = "TRUE",direction=-1) + 
  geom_bar(position="dodge",stat="identity") +
  geom_errorbar(aes(ymin=lowerse,ymax=upperse,width=0.2), 
                position=position_dodge(0.9)) +
  ylab("Surface Height (m)") +
  #ggtitle("Surface Heights") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(fill="Veg. Class") +
  coord_cartesian(ylim=c(178.5,179.5)) +
  #scale_fill_discrete(name = "Veg. Class", labels = c("EFB + Lily spp.","Emergent Veg.","Duckweed spp.","Sub. Veg.")) +
  scale_x_discrete(labels=c("DSM"))


aug_grid <- grid.arrange(aug_spectra_plot,aug_surf_rug_plot,aug_surf_plot,ncol=3)

pdf(file = "Z:/UAV2019/Alpena/Composites11cm/ClassProfiles_UnlumpedDSM_11cm.pdf",
    width = 5,
    height = 4)

#####----july-----#####
july_prof_features<- reshape2::melt(alpenagt_all,id=c("train"),
                          measure=c("Green7_12","Red7_12","Red.Edge7_12","NIR7_12","NDVI7_12","DSM7_12")) 
names(july_prof_features)[3] <- "Reflectance" 
#average 
july_avg<-reshape2::dcast(july_prof_features,train+variable~.,fun=mean, value.var="Reflectance",na.rm=TRUE)
names(july_avg)[3]<-"mean_ref"
#se
july_se<-reshape2::dcast(july_prof_features, train+variable~., fun=sd, value.var="Reflectance",na.rm=TRUE)
names(july_se)[3]<-"sd"
table(alpenagt_all$train)
#first reps is class sample n, second is reps are for band #'s (4 classes with varying n, 6 bands for each (figure out rugosity issue for this dataset))
Class_length<-c(rep(61,6),rep(98,6),rep(62,6),rep(83,6))
as.data.frame(Class_length)
july_se$class_n<-Class_length
#jenky way, but works  
july_se$se <- july_se$sd/sqrt(july_se$class_n)

july_spectra<-merge(july_avg, july_se, by=c("train","variable"), na.rm=TRUE)
july_spectra$lowerse<-july_spectra$mean_ref-july_spectra$se
july_spectra$upperse<-july_spectra$mean_ref+july_spectra$se

july_spectra$train <- as.factor(july_spectra$train)

july_spectra$train <- as.factor(july_spectra$train)
july_spectra_all <- july_spectra %>% mutate(train=as.character(train),
                                        train= if_else(train == 1,"EFB",train),
                                        train=as.factor(train))  %>% mutate(train=as.character(train),
                                                                            train= if_else(train == 2,"Yellow Veg.",train),
                                                                            train=as.factor(train)) %>% mutate(train=as.character(train),
                                                                                                               train= if_else(train == 3,"Sub.Veg",train),
                                                                                                               train=as.factor(train))  %>% mutate(train=as.character(train),
                                                                                                                                                   train= if_else(train == 4,"Emergent",train),
                                                                                                                                                   train=as.factor(train))
#specral plot 
jul_spectra_only <- subset(july_spectra_all, variable !="DSM7_12")
july_spectra <- ggplot(jul_spectra_only,aes(x=variable,y=mean_ref,group=train,color=train))+
  geom_line(size=1.25)+
  geom_errorbar(aes(ymin=lowerse,ymax=upperse,width=0.05))+
  ylab("% Reflectance")+
  xlab("Band Name")+
  ggtitle("Spectral Profiles, 7/2019")+
  theme_bw()+
  theme(legend.position='none')+
  # scale_color_discrete(breaks=c("Frogbit","Lily","SV","Typha + EFB/Litter"),
  #                      labels=c("Frogbit","Lily","SV","Typha + EFB/Litter")) +
  theme(axis.title.x = element_text(size = 14, face = "bold")) +
  scale_x_discrete(labels=c("Green","Red","Red Edge","NIR","NDVI"))


#dsm, rugosity is funky for this dataset... 
july_surface <- subset(july_spectra_all,variable == "DSM7_12")
july_surf <- ggplot(july_surface,aes(x=variable,y=mean_ref,group=train,fill=train))+
  geom_bar(position="dodge",stat="identity") +
  geom_errorbar(aes(ymin=lowerse,ymax=upperse,width=0.2), 
                position=position_dodge(0.9)) +
  ylab("Surface Height (m)")+
  xlab("Band Name")+
  ggtitle("Surface Heights, 7/2019")+
  theme_bw()+
  labs(fill="Veg. Class")+
  coord_cartesian(ylim=c(179,181.5)) +
  # scale_color_discrete(breaks=c("Frogbit","Lily","SV","Typha + EFB/Litter"),
  #                      labels=c("Frogbit","Lily","SV","Typha + EFB/Litter")) +
  theme(axis.title.x = element_text(size = 14, face = "bold")) +
  scale_x_discrete(labels=c("DSM"))

legend <- get_legend(july_surf)
#as_ggplot(legend)
# #plot side by side 
july_grid <- grid.arrange(july_spectra, july_surf, ncol=2)

#altogether now 
grid_all <- grid.arrange(aug_grid, july_grid,nrow=2)
#pretty good separability with yellow veg 

####-----unlumped chm profiles-----#####
aug_chm <- read.csv("AlpenaGT_11cm_CHM.csv")

aug_chm_prof_features<- reshape2::melt(aug_chm,id=c("train"),
                                   measure=c("Green8_5","Red8_5","Red.Edge8_5","NIR8_5","NDVI8_5","CHM8_5","CHM_Rugosity8_5"))
names(aug_chm_prof_features)[3] <- "Reflectance"

#average 
aug_avg<-reshape2::dcast(aug_chm_prof_features,train+variable~.,fun=mean, value.var="Reflectance",na.rm=TRUE)
names(aug_avg)[3]<-"mean_ref"
#se
aug_se<-reshape2::dcast(aug_chm_prof_features, train+variable~., fun=sd, value.var="Reflectance",na.rm=TRUE)
names(aug_se)[3]<-"sd"
table(aug_chm$train)#order for below command 
#first reps is train sample n, second is reps are for band #'s (4 classes with varying n, now 7 bands in this aug data)
#numerical order, (efb,yv,sv,emergent)
Class_length<-c(rep(60,7),rep(93,7),rep(51,7),rep(81,7))
as.data.frame(Class_length)
aug_se$class_n<-Class_length
#jenky way, but works  
aug_se$se <- aug_se$sd/sqrt(aug_se$class_n)

aug_spectra <- merge(aug_avg, aug_se, by=c("train","variable"), na.rm=TRUE)
aug_spectra$lowerse <- aug_spectra$mean_ref-aug_spectra$se
aug_spectra$upperse <- aug_spectra$mean_ref+aug_spectra$se

aug_spectra$train <- as.factor(aug_spectra$train)
# aug_all <- aug_spectra %>% mutate(train=as.character(train),
#                                   train= if_else(train == 1,"EFB",train),
#                                   train=as.factor(train))  %>% mutate(train=as.character(train),
#                                                                       train= if_else(train == 2,"Yellow Veg.",train),
#                                                                       train=as.factor(train)) %>% mutate(train=as.character(train),
#                                                                                                          train= if_else(train == 3,"Sub.Veg",train),
#                                                                                                          train=as.factor(train))  %>% mutate(train=as.character(train),
#                                                                                                                                              train= if_else(train == 4,"Emergent",train),
#                                                                                                                                              train=as.factor(train))
# 


aug_spectra$train <- if_else(aug_spectra$train == 1, "EFB + Lily spp.", as.character(aug_spectra$train))
aug_spectra$train <- if_else(aug_spectra$train == 2, "Duckweed spp.", as.character(aug_spectra$train))
aug_spectra$train <- if_else(aug_spectra$train == 3, "Sub.Veg.", as.character(aug_spectra$train))
aug_spectra$train <- if_else(aug_spectra$train == 4, "Emergent Veg.", as.character(aug_spectra$train))


august_wlegend <- ggplot(aug_spectra,aes(x=variable,y=mean_ref,group=train,color=train)) +
  geom_line(size=1.25)+
  geom_errorbar(aes(ymin=lowerse,ymax=upperse,width=0.05))+
  ylab("Reflectance Value")+
  xlab("Band Name")+
  ggtitle("Veg. Classes Spectral/Textural Profiles in Alpena, MI, August 2019")+
  theme_bw()+
  theme(axis.title.x = element_text(size = 14, face = "bold")) +
  # scale_color_discrete(breaks=c("Frogbit","Lily","SV","Typha + EFB/Litter"),
  #                      labels=c("Frogbit","Lily","SV","Typha + EFB/Litter")) +
  scale_x_discrete(labels=c("Green","Red","Red Edge","NIR","NDVI","StdVegHeight"))

#Actual spectral plot 
aug_spectra_only <- subset(aug_spectra, variable !="CHM8_5") 
aug_spectra_only <- subset(aug_spectra_only, variable !="CHM_Rugosity8_5") 

# aug_spectra <- ggplot(aug_spectra_only,aes(x=variable,y=mean_ref,group=train,color=train)) +
#   geom_line(size=1.25)+
#   geom_errorbar(aes(ymin=lowerse,ymax=upperse,width=0.05))+
#   ylab("% Reflectance")+
#   ggtitle("Spectral Profiles, 8/2019")+
#   theme_bw()+
#   theme(legend.position='none',
#         axis.title.x = element_blank()) +
#   # scale_color_discrete(breaks=c("Frogbit","Lily","SV","Typha + EFB/Litter"),
#   #                      labels=c("Frogbit","Lily","SV","Typha + EFB/Litter")) +
#   scale_x_discrete(labels=c("Green","Red","Red Edge","NIR","NDVI"))

#rug only 
aug_rugosity <- subset(aug_spectra, variable == "CHM_Rugosity8_5")
aug_rugosity$train <- factor(aug_rugosity$train,
                            levels = c("Duckweed spp.","EFB + Lily spp.","Sub.Veg.","Emergent Veg."))

aug_chm_rug_plot <- ggplot(aug_rugosity,aes(x=variable,y=mean_ref,fill=train)) +
  scale_fill_viridis(discrete = "TRUE",direction=-1) + 
  geom_bar(position="dodge",stat="identity") +
  geom_errorbar(aes(ymin=lowerse,ymax=upperse,width=0.2), 
                position=position_dodge(0.9)) +
  ylab("SD of Canopy Height (m)")+
  #xlab("Band Name")+
  #ggtitle("CHM Texture")+
  theme_bw()+
  theme(legend.position='none',
        axis.title.x=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  # scale_color_discrete(breaks=c("Frogbit","Lily","SV","Typha + EFB/Litter"),
  #                      labels=c("Frogbit","Lily","SV","Typha + EFB/Litter")) +
  scale_x_discrete(labels=c("Rugosity"))

#surface height only 
aug_canopy <- subset(aug_spectra, variable == "CHM8_5")
aug_canopy$train <- factor(aug_canopy$train,
                            levels = c("Duckweed spp.","EFB + Lily spp.","Sub.Veg.","Emergent Veg."))

aug_canopy_plot <- ggplot(aug_canopy,aes(x=variable,y=mean_ref,fill=train)) +
  scale_fill_viridis(discrete = "TRUE", direction=-1) + 
  geom_bar(position="dodge",stat="identity") +
  geom_errorbar(aes(ymin=lowerse,ymax=upperse,width=0.2), 
                position=position_dodge(0.9)) +
  ylab("Canopy Height (m)") +
  #ggtitle("Canopy Heights") +
  theme_bw() +
  theme(axis.title.x = element_blank()) +
  labs(fill="Veg. Class") +
  coord_cartesian(ylim=c(0,0.16)) +
  theme(legend.position='none',
        axis.title.x=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  #scale_fill_discrete(name = "Veg. Class", labels = c("Floating Veg.","Yellow Veg.","Sub.Veg.","Emergent Veg.")) +
  scale_x_discrete(labels=c("CHM"))


aug_grid <- plot_grid(aug_spectra_plot,aug_chm_rug_plot,aug_surf_rug_plot,aug_canopy_plot,aug_surf_plot, 
                         ncol=3,nrow=2,labels=c("a)","b)","c)","d)","e)"))
aug_grid
ggsave("FigureS1.jpg",width=9,height=6,units=c("in"),dpi=300)


#####---- using all predictor vars across dates (will keep efb/lil and yv separate for now)-----######
all_vars <- subset(alpenagt_all,
                   select=-c(X,Y))
all_vars <- drop_na(all_vars)

model_allpred <- train(as.factor(train)~., # Class is a function of the variables we decided to include
                       data = all_vars, # Use the train data frame as the training data
                       method = 'rf',# Use the 'random forest' algorithm
                       trControl = trainControl(method = 'repeatedcv', # Use cross-validation
                                                number = 5,
                                                repeats = 50,
                                                classProbs = T,
                                                savePredictions = T)) 
model_allpred

confusionMatrix.train(model_allpred)
plot(varImp(model_allpred,scale=F))

#alternative Error Matrix 
sub_all <- subset(model_allpred$pred,model_allpred$pred$mtry==model_allpred$bestTune$mtry)
caret::confusionMatrix(table(sub_all$pred,sub_all$obs))

#create own df of confm values
Floating <- c(2008,822,71,349)
YellowVeg <- c(420,4342,288,300)
SubVeg <- c(135,472,893,250)
Emergent <- c(324,27,85,3714)

#Ref <- c("Floating","SubVeg","Emergent")
conf_df <- data.frame(Floating,YellowVeg,SubVeg,Emergent)
con_matx <- as.matrix(conf_df)
con_diag <- diag(con_matx)
n <- sum(con_matx)

avg_OA <- sum(con_diag) / n
avg_OA
# 75.6%

# observed (true) cases per class
rowsums <- apply(con_matx, 1, sum)

# predicted cases per class
colsums <- apply(con_matx, 2, sum)

# Producer accuracy
avg_PA <- con_diag / colsums
# User accuracy
avg_UA <- con_diag / rowsums
outAcc_all <- data.frame(producerAccuracy = avg_PA, userAccuracy = avg_UA)
outAcc_all




####----aug only unlumped -----######
all_vars_aug <- subset(aug_dsm,
                       select=-c(X.1,Green7_12,Red7_12,Red.Edge7_12,NIR7_12,NDVI7_12,DSM7_12, StdVegHeight7_12))
all_vars_aug <- drop_na(all_vars_aug)

tunegrid <- expand.grid(.mtry = (1:7)) 

model_aug <- train(make.names(train)~., # Class is a function of the variables we decided to include
                   data = all_vars_aug, # Use the train data frame as the training data
                   method = 'rf',
                   tuneGrid = tunegrid, # Use the 'random forest' algorithm
                   trControl = trainControl(method = 'repeatedcv', # Use cross-validation
                                            number = 5,
                                            repeats=50,
                                            classProbs = T,
                                            savePredictions=T)) 
model_aug

confusionMatrix.train(model_aug)
plot(varImp(model_aug,scale=F))

#alternative Error Matrix 
sub_aug <- subset(model_aug$pred,model_aug$pred$mtry==model_aug$bestTune$mtry)
caret::confusionMatrix(table(sub_aug$pred,sub_aug$obs))

#for precision, recall, and f1 score 
caret::confusionMatrix(table(sub_aug$pred,sub_aug$obs), mode="prec_recall")


#create own df of confm values
Floating <- c(4317,276,407)
SubVeg <- c(947,2090,263)
Emergent <- c(366,45,3689)


#Ref <- c("Floating","SubVeg","Emergent")
conf_df <- data.frame(Floating,SubVeg,Emergent)
con_matx <- as.matrix(conf_df)
con_diag <- diag(con_matx)
n <- sum(con_matx)

avg_OA <- sum(con_diag) / n
avg_OA
#81.3

# observed (true) cases per class
rowsums <- apply(con_matx, 1, sum)

# predicted cases per class
colsums <- apply(con_matx, 2, sum)

# Producer accuracy
avg_PA <- con_diag / colsums
# User accuracy
avg_UA <- con_diag / rowsums
outAcc_aug <- data.frame(producerAccuracy = avg_PA, userAccuracy = avg_UA)
outAcc_aug

rocs_aug <- llply(unique(model_aug$pred$obs), function(cls) {
  roc(response = model_aug$pred$obs==cls, 
      predictor = model_aug$pred[,as.character(cls)])
})

plot(rocs_aug[[1]], col="green", legacy.axes=TRUE, xlim=c(1,0), ylim=c(0,1),
     xlab="False Positive Rate", ylab="True Positive Rate")
l_ply(rocs_aug[2], plot, col=c("blue"), add=T)
l_ply(rocs_aug[3], plot, col=c("salmon"), add=T)
axis(1,seq(1,0,0.2))

legend("bottomright",title="Class ROC", legend=c("Floating","Submergent","Emergent"),
       col=c("green","blue","salmon"), lwd=1)

aucs <- laply(rocs_aug, function(x) as.numeric(x$auc))
summary(aucs)
mean(aucs) #all are one-versus-others, so take mean to get overall model performance 

text(x=1.000000e-01,y=0.5,"Mean AUC = 0.88")

eers <- laply(rocs_aug, function(x) x$sensitivities[which.min(abs(x$sensitivities-x$specificities))])
summary(eers)

conf <- ldply(unique(model_aug$pred$obs), function(cls) data.frame(tp=sum(model_aug$pred$pred[model_aug$pred$obs==cls]==cls), 
                                                                   fn=sum(model_aug$pred$pred[model_aug$pred$obs==cls]!=cls), 
                                                                   tn=sum(model_aug$pred$pred[model_aug$pred$obs!=cls]!=cls),
                                                                   fp=sum(model_aug$pred$pred[model_aug$pred$obs!=cls]==cls)))
conf$tpr <- with(conf, tp/(tp+fn))
conf$tnr <- with(conf, tn/(tn+fp))
summary(conf[,5:6])

#####---- all pred vars floating/yv lumped----###### 
aug_dsm <- read.csv("AlpenaGT_11cm_DSM.csv")
aug_dsm$train <- if_else(aug_dsm$train == 2, 1, 
                         as.double(aug_dsm$train))

#take out a good random portion of floating, way too many
(floatingindex <- which(aug_dsm$train == '1'))
(deleteindex <- sample(floatingindex, length(floatingindex) - 100)) 
aug_dsm <- aug_dsm[-deleteindex,]

aug_dsm$train <- if_else(aug_dsm$train == 1, "Floating Veg.", as.character(aug_dsm$train))
aug_dsm$train <- if_else(aug_dsm$train == 3, "Sub.Veg.", as.character(aug_dsm$train))
aug_dsm$train <- if_else(aug_dsm$train == 4, "Emergent Veg.", as.character(aug_dsm$train))

####----new plots with lumped classes----#### 
aug_prof_features<- reshape2::melt(aug_dsm,id=c("train"),
                                   measure=c("Green8_5","Red8_5","Red.Edge8_5","NIR8_5","NDVI8_5","DSM8_5","StdVegHeight8_5"))
names(aug_prof_features)[3] <- "Reflectance"

#average 
aug_avg<-reshape2::dcast(aug_prof_features,train+variable~.,fun=mean, value.var="Reflectance",na.rm=TRUE)
names(aug_avg)[3]<-"mean_ref"
#se
aug_se<-reshape2::dcast(aug_prof_features, train+variable~., fun=sd, value.var="Reflectance",na.rm=TRUE)
names(aug_se)[3]<-"sd"
table(aug_dsm$train)#order for below command 
#first reps is train sample n, second is reps are for band #'s (4 classes with varying n, now 7 bands in this aug data)
#numerical order, (emergent,floating, sub veg)
Class_length<-c(rep(82,7),rep(100,7),rep(66,7))
as.data.frame(Class_length)
aug_se$class_n<-Class_length
#jenky way, but works  
aug_se$se <- aug_se$sd/sqrt(aug_se$class_n)

aug_spectra <- merge(aug_avg, aug_se, by=c("train","variable"), na.rm=TRUE)
aug_spectra$lowerse <- aug_spectra$mean_ref-aug_spectra$se
aug_spectra$upperse <- aug_spectra$mean_ref+aug_spectra$se

aug_spectra$train <- as.factor(aug_spectra$train)
# aug_all <- aug_spectra %>% mutate(train=as.character(train),
#                                   train= if_else(train == 1,"Floating Veg.",train),
#                                   train=as.factor(train))  %>% mutate(train=as.character(train),
#                                                                       train= if_else(train == 2,"Sub.Veg.",train),
#                                                                       train=as.factor(train)) %>% mutate(train=as.character(train),
#                                                                                                          train= if_else(train == 3,"Emergent Veg.",train),
#                                                                                                          train=as.factor(train)) 

aug_spectra$train <- if_else(aug_spectra$train == 1, "Floating Veg.", as.character(aug_spectra$train))
aug_spectra$train <- if_else(aug_spectra$train == 2, "Sub.Veg.", as.character(aug_spectra$train))
aug_spectra$train <- if_else(aug_spectra$train == 3, "Emergent Veg.", as.character(aug_spectra$train))

aug_all <- aug_spectra

#Actual spectral plot 
aug_spectra_only <- subset(aug_all, variable !="DSM8_5") 
aug_spectra_only <- subset(aug_spectra_only, variable !="StdVegHeight8_5") 

aug_spectra_only$train <- factor(aug_spectra_only$train,
                                 levels = c("Emergent Veg.","Floating Veg.","Sub.Veg."))

aug_spectra_plot <- ggplot(aug_spectra_only,aes(x=variable,y=mean_ref,group=train,color=train)) +
  scale_color_viridis(discrete = TRUE,direction=-1) + 
  geom_line(size=1.25)+
  geom_errorbar(aes(ymin=lowerse,ymax=upperse,width=0.05))+
  ylab("% Reflectance")+
  #ggtitle("Class Spectral Profiles")+
  theme_bw()+
  theme(legend.position='none',
        axis.title.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  # scale_color_discrete(breaks=c("Frogbit","Lily","SV","Typha + EFB/Litter"),
  #                      labels=c("Frogbit","Lily","SV","Typha + EFB/Litter")) +
  scale_x_discrete(labels=c("Green","Red","Red Edge","NIR","NDVI"))


#rug only 
aug_rugosity <- subset(aug_all, variable == "StdVegHeight8_5")
aug_rugosity$train <- factor(aug_rugosity$train,
                                 levels = c("Emergent Veg.","Floating Veg.","Sub.Veg."))


aug_rug_dsm_lumped <- ggplot(aug_rugosity,aes(x=variable,y=mean_ref,fill=train)) +
  scale_fill_viridis(discrete = TRUE,direction=-1) + 
  geom_bar(position="dodge",stat="identity") +
  geom_errorbar(aes(ymin=lowerse,ymax=upperse,width=0.2), 
                position=position_dodge(0.9)) +
  ylab("SD of Surface Height (m)")+
  #ggtitle("DSM Texture, 8/2019")+
  theme_bw()+
  theme(legend.position='none',
        axis.title.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  # scale_color_discrete(breaks=c("Frogbit","Lily","SV","Typha + EFB/Litter"),
  #                      labels=c("Frogbit","Lily","SV","Typha + EFB/Litter")) +
  scale_x_discrete(labels=c("Rugosity"))

#surface height only 
aug_surface <- subset(aug_all, variable == "DSM8_5")
aug_surface$train <- factor(aug_surface$train,
                             levels = c("Emergent Veg.","Floating Veg.","Sub.Veg."))

aug_surf_lumped <- ggplot(aug_surface,aes(x=variable,y=mean_ref,fill=train)) +
  scale_fill_viridis(discrete = TRUE,direction=-1) +
  geom_bar(position="dodge",stat="identity") +
  geom_errorbar(aes(ymin=lowerse,ymax=upperse,width=0.2), 
                position=position_dodge(0.9)) +
  ylab("Surface Height (m)")+
  #ggtitle("Surface Heights, 8/2019")+
  theme_bw()+
  theme(axis.title.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + 
  labs(fill="Veg. Class")+
  coord_cartesian(ylim=c(178.6,179.4)) +
  #scale_fill_discrete(name = "Veg. Class", labels = c("Floating Veg.","Sub.Veg.","Emergent Veg.")) +
  scale_x_discrete(labels=c("DSM"))

aug_grid_lumped <- grid.arrange(aug_spectra_plot,aug_rug_dsm_lumped,
                         aug_surf_lumped,ncol=3)


######-----july-----####### 
july_prof_features <- reshape2::melt(alpenagt_all,id=c("train"),
                                    measure=c("Green7_12","Red7_12","Red.Edge7_12","NIR7_12","NDVI7_12","DSM7_12","StdVegHeight7_12")) 
names(july_prof_features)[3] <- "Reflectance" 
#average 
july_avg<-reshape2::dcast(july_prof_features,train+variable~.,fun=mean, value.var="Reflectance",na.rm=TRUE)
names(july_avg)[3]<-"mean_ref"
#se
july_se<-reshape2::dcast(july_prof_features, train+variable~., fun=sd, value.var="Reflectance",na.rm=TRUE)
names(july_se)[3]<-"sd"
table(alpenagt_all$train)
#first reps is class sample n, second is reps are for band #'s (4 classes with varying n, 7 bands for each)
Class_length<-c(rep(100,7),rep(66,7),rep(82,7))
as.data.frame(Class_length)
july_se$class_n<-Class_length
#jenky way, but works  
july_se$se <- july_se$sd/sqrt(july_se$class_n)

july_spectra<-merge(july_avg, july_se, by=c("train","variable"), na.rm=TRUE)
july_spectra$lowerse<-july_spectra$mean_ref-july_spectra$se
july_spectra$upperse<-july_spectra$mean_ref+july_spectra$se

july_spectra$train <- as.factor(july_spectra$train)

july_spectra$train <- as.factor(july_spectra$train)
july_spectra_all <- july_spectra %>% mutate(train=as.character(train),
                                            train= if_else(train == 1,"Floating Veg.",train),
                                            train=as.factor(train))  %>% mutate(train=as.character(train),
                                                                                train= if_else(train == 2,"Sub.Veg.",train),
                                                                                train=as.factor(train)) %>% mutate(train=as.character(train),
                                                                                                                   train= if_else(train == 3,"Emergent Veg.",train),
                                                                                                                   train=as.factor(train))
#specral plot 
jul_spectra_only <- subset(july_spectra_all, variable !="StdVegHeight7_12")
jul_spectra_only <- subset(jul_spectra_only, variable !="DSM7_12")

july_spectra <- ggplot(jul_spectra_only,aes(x=variable,y=mean_ref,group=train,color=train))+
  geom_line(size=1.25)+
  geom_errorbar(aes(ymin=lowerse,ymax=upperse,width=0.05))+
  ylab("% Reflectance")+
  ggtitle("Spectral Profiles, 7/2019")+
  theme_bw()+
  theme(legend.position='none', 
        axis.title.x = element_blank())+
  # scale_color_discrete(breaks=c("Frogbit","Lily","SV","Typha + EFB/Litter"),
  #                      labels=c("Frogbit","Lily","SV","Typha + EFB/Litter")) +
  scale_x_discrete(labels=c("Green","Red","Red Edge","NIR","NDVI"))

#rug only 
july_rugosity <- subset(july_spectra_all, variable == "StdVegHeight7_12")
july_rug_lumped <- ggplot(july_rugosity,aes(x=variable,y=mean_ref,fill=train)) +
  geom_bar(position="dodge",stat="identity") +
  geom_errorbar(aes(ymin=lowerse,ymax=upperse,width=0.2), 
                position=position_dodge(0.9)) +
  ylab("Rugosity (Std Dev. of Surface Height)")+
  ggtitle("Textural Profiles, 7/2019")+
  theme_bw()+
  theme(legend.position='none',
        axis.title.x = element_blank()) +
  # scale_color_discrete(breaks=c("Frogbit","Lily","SV","Typha + EFB/Litter"),
  #                      labels=c("Frogbit","Lily","SV","Typha + EFB/Litter")) +
  scale_x_discrete(labels=c("Rugosity"))

#dsm 
july_surface <- subset(july_spectra_all,variable == "DSM7_12")
july_surf <- ggplot(july_surface,aes(x=variable,y=mean_ref,group=train,fill=train))+
  geom_bar(position="dodge",stat="identity") +
  geom_errorbar(aes(ymin=lowerse,ymax=upperse,width=0.2), 
                position=position_dodge(0.9)) +
  ylab("Surface Height (m)")+
  ggtitle("Surface Heights, 7/2019")+
  theme_bw()+
  theme(axis.title.x = element_blank()) + 
  labs(fill="Veg. Class")+
  coord_cartesian(ylim=c(178.75,179.4)) +
  # scale_color_discrete(breaks=c("Frogbit","Lily","SV","Typha + EFB/Litter"),
  #                      labels=c("Frogbit","Lily","SV","Typha + EFB/Litter")) +
  scale_x_discrete(labels=c("DSM"))

# #plot side by side 
july_grid <- grid.arrange(july_spectra, july_rug_lumped, july_surf, ncol=3)

#altogether now 
grid_all <- grid.arrange(aug_grid, july_grid,nrow=2)
 

#####-----lumped chm plots ----####
aug_chm <- read.csv("AlpenaGT_11cm_CHM.csv")

aug_chm$train <- if_else(aug_chm$train == 2, 1, 
                         as.double(aug_chm$train))

#take out a good random portion of floating, way too many
(floatingindex <- which(aug_chm$train == '1'))
(deleteindex <- sample(floatingindex, length(floatingindex) - 100)) 
aug_chm <- aug_chm[-deleteindex, ]

aug_chm$train <- if_else(aug_chm$train == 1, "Floating Veg.", as.character(aug_chm$train))
aug_chm$train <- if_else(aug_chm$train == 3, "Sub. Veg.", as.character(aug_chm$train))
aug_chm$train <- if_else(aug_chm$train == 4, "Emergent Veg.", as.character(aug_chm$train))

aug_chm_prof_features<- reshape2::melt(aug_chm,id=c("train"),
                                       measure=c("Green8_5","Red8_5","Red.Edge8_5","NIR8_5","NDVI8_5","CHM8_5","CHM_Rugosity8_5"))
names(aug_chm_prof_features)[3] <- "Reflectance"

#average 
aug_avg<-reshape2::dcast(aug_chm_prof_features,train+variable~.,fun=mean, value.var="Reflectance",na.rm=TRUE)
names(aug_avg)[3]<-"mean_ref"
#se
aug_se<-reshape2::dcast(aug_chm_prof_features, train+variable~., fun=sd, value.var="Reflectance",na.rm=TRUE)
names(aug_se)[3]<-"sd"
table(aug_chm$train)#order for below command 
#first reps is train sample n, second is reps are for band #'s (4 classes with varying n, now 7 bands in this aug data)
#numerical order, (emergent, flloating, sv)
Class_length<-c(rep(81,7),rep(100,7),rep(51,7))
as.data.frame(Class_length)
aug_se$class_n<-Class_length
#jenky way, but works  
aug_se$se <- aug_se$sd/sqrt(aug_se$class_n)

aug_spectra <- merge(aug_avg, aug_se, by=c("train","variable"), na.rm=TRUE)
aug_spectra$lowerse <- aug_spectra$mean_ref-aug_spectra$se
aug_spectra$upperse <- aug_spectra$mean_ref+aug_spectra$se

aug_spectra$train <- as.factor(aug_spectra$train)

#Actual spectral plot 
aug_spectra_only <- subset(aug_spectra, variable !="CHM8_5") 
aug_spectra_only <- subset(aug_spectra_only, variable !="CHM_Rugosity8_5") 

aug_spectra_only$train <- factor(aug_spectra_only$train,
                                 levels = c("Emergent Veg.","Floating Veg.","Sub. Veg."))

# aug_spectra_plot <- ggplot(aug_spectra_only,aes(x=variable,y=mean_ref,group=train,color=train)) +
#   scale_color_viridis(discrete = TRUE,direction=-1) + 
#   geom_line(size=1.25)+
#   geom_errorbar(aes(ymin=lowerse,ymax=upperse,width=0.05))+
#   ylab("% Reflectance")+
#   ggtitle("Spectral Profiles, 8/2019")+
#   theme_bw()+
#   theme(legend.position='none',
#         axis.title.x = element_blank()) +
#   # scale_color_discrete(breaks=c("Frogbit","Lily","SV","Typha + EFB/Litter"),
#   #                      labels=c("Frogbit","Lily","SV","Typha + EFB/Litter")) +
#   scale_x_discrete(labels=c("Green","Red","Red Edge","NIR","NDVI"))

#rug only 
aug_rugosity <- subset(aug_spectra, variable == "CHM_Rugosity8_5")
aug_rugosity$train <- factor(aug_rugosity$train,
                             levels = c("Emergent Veg.","Floating Veg.","Sub. Veg."))

aug_rug_canopy_plot <- ggplot(aug_rugosity,aes(x=variable,y=mean_ref,fill=train)) +
  scale_fill_viridis(discrete = TRUE,direction=-1) + 
  geom_bar(position="dodge",stat="identity") +
  geom_errorbar(aes(ymin=lowerse,ymax=upperse,width=0.2), 
                position=position_dodge(0.9)) +
  ylab("SD of Canopy Height (m)")+
  #xlab("Band Name")+
  #ggtitle("CHM Texture, 8/2019")+
  theme_bw()+
  theme(legend.position='none',
        axis.title.x=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  # scale_color_discrete(breaks=c("Frogbit","Lily","SV","Typha + EFB/Litter"),
  #                      labels=c("Frogbit","Lily","SV","Typha + EFB/Litter")) +
  scale_x_discrete(labels=c("Rugosity"))

#surface height only 
aug_canopy <- subset(aug_spectra, variable == "CHM8_5")
aug_canopy$train <- factor(aug_canopy$train,
                            levels = c("Emergent Veg.","Floating Veg.","Sub. Veg."))

aug_canopy_plot <- ggplot(aug_canopy,aes(x=variable,y=mean_ref,fill=train)) +
  scale_fill_viridis(discrete = TRUE,direction=-1) +
  geom_bar(position="dodge",stat="identity") +
  geom_errorbar(aes(ymin=lowerse,ymax=upperse,width=0.2), 
                position=position_dodge(0.9)) +
  ylab("Canopy Height (m)") +
  #ggtitle("Canopy Heights, 8/2019") +
  theme_bw() +
  theme(legend.position='none',
         axis.title.x=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(fill="Veg. Class") +
  coord_cartesian(ylim=c(0,0.16)) +
  scale_x_discrete(labels=c("CHM"))


aug_grid <- plot_grid(aug_spectra_plot,aug_rug_canopy_plot,aug_rug_dsm_lumped, aug_canopy_plot,aug_surf_lumped,
                      ncol=3, nrow=2, labels=c("a)","b)","c)","d)","e)"))

#####----first rf model----#####
all_vars_lumped <- subset(alpenagt_all,
                          select=-c(X,Y))
all_vars_lumped <- drop_na(all_vars_lumped)


#recursive feature elimination 
ctrl <- rfeControl(functions=rfFuncs,
                   method="repeatedcv",
                   repeats=5,
                   verbose = FALSE)
subsets <- c(1:14)
rfProfile <- rfe(all_vars_lumped[2:15],
                 all_vars_lumped$train,
                 sizes=subsets,
                 rfeControl = ctrl)
predictors(rfProfile)
#consistent with var importance plots, run model with 12 predictors? 
#elimnate both rugosity metrics 
all_vars_subset <- subset(all_vars_lumped,
                          select=c("train","Red8_5","Green8_5","NIR8_5","NDVI7_12","NDVI8_5","Red.Edge8_5","Red7_12","NIR7_12",
                                   "DSM8_5","DSM7_12","Green7_12","Red.Edge7_12"))
model_subset <- train(make.names(train)~., # Class is a function of the variables we decided to include
                          data = all_vars_subset, # Use the train data frame as the training data
                          method = 'rf',# Use the 'random forest' algorithm
                          trControl = trainControl(method = 'repeatedcv', # Use cross-validation
                                                   number = 5,
                                                   repeats = 50,
                                                   classProbs = T, #need to include for AUC/ROC, took out as.factor() on train 
                                                   savePredictions = T)) #get warning with these params
model_subset

# model_all_lumped <- train(as.factor(train)~., # Class is a function of the variables we decided to include
#                           data = all_vars_lumped, # Use the train data frame as the training data
#                           method = 'rf',# Use the 'random forest' algorithm
#                           trControl = trainControl(method = 'repeatedcv', # Use cross-validation
#                                                    number = 5,
#                                                    repeats = 50,
#                                                    classProbs = T, #need to include for AUC/ROC, took out as.factor() on train 
#                                                    savePredictions = T,
#                                                    verboseIter = T)) #get warning with these params
# model_all_lumped

#try this way 
#levels(all_vars_lumped$train) <- c("Floating","Submerged","Emergent")

control <- trainControl(method='repeatedcv', 
                        number=5, 
                        repeats=50, 
                        search='grid',
                        classProbs = T, #need to include for AUC/ROC, took out as.factor() on train 
                        savePredictions = T)

tunegrid <- expand.grid(.mtry = (1:14)) 

rf_gridsearch <- train(make.names(train)~., 
                       data = all_vars_lumped,
                       method = 'rf',
                       metric = 'Accuracy',
                       tuneGrid = tunegrid)
print(rf_gridsearch)
plot(rf_gridsearch)

library(e1071)
bestMtry <- tuneRF(all_vars_lumped[2:15],all_vars_lumped$train, stepFactor = 1.5, improve = 1e-5, ntree = 500)

model_all_lumped <- train(make.names(train)~., # Class is a function of the variables we decided to include
                          data = all_vars_lumped, # Use the train data frame as the training data
                          method = 'rf',# Use the 'random forest' algorithm
                          trControl = trainControl(method = 'repeatedcv', # Use cross-validation
                                                   number = 5,
                                                   repeats = 50,
                                                   classProbs = T, #need to include for AUC/ROC, took out as.factor() on train 
                                                   savePredictions = T)) #get warning with these params
model_all_lumped

confusionMatrix.train(model_all_lumped)
plot(varImp(model_all_lumped,scale=F))

#alternative Error Matrix 
sub_all_lumped <- subset(model_all_lumped$pred,model_all_lumped$pred$mtry==model_all_lumped$bestTune$mtry)
caret::confusionMatrix(table(sub_all_lumped$pred,sub_all_lumped$obs))

#good explanation of performance statistics here 
# https://topepo.github.io/caret/measuring-performance.html#measures-for-predicted-classes

#for precision, recall, and f1 score 
caret::confusionMatrix(table(sub_all_lumped$pred,sub_all_lumped$obs), mode="prec_recall")

#create own df of confm values
Floating <- c(4369,324,307)
SubVeg <- c(838,1957,305)
Emergent <- c(129,67,3954)


#Ref <- c("Floating","SubVeg","Emergent")
conf_df <- data.frame(Floating,SubVeg,Emergent)
con_matx <- as.matrix(conf_df)
con_diag <- diag(con_matx)
n <- sum(con_matx)

avg_OA <- sum(con_diag) / n
avg_OA
# 84.4%. 

# observed (true) cases per class
rowsums <- apply(con_matx, 1, sum)

# predicted cases per class
colsums <- apply(con_matx, 2, sum)

# Producer accuracy
avg_PA <- con_diag / colsums
# User accuracy
avg_UA <- con_diag / rowsums
outAcc_all_lumped <- data.frame(producerAccuracy = avg_PA, userAccuracy = avg_UA)
outAcc_all_lumped

#ROC and AUC 
library(pROC)
library(plyr)

rocs_all <- llply(unique(model_all_lumped$pred$obs), function(cls) {
   roc(response = model_all_lumped$pred$obs==cls, 
       predictor = model_all_lumped$pred[,as.character(cls)])
  })

plot(rocs_all[[1]], col="green", legacy.axes=TRUE, xlim=c(1,0), ylim=c(0,1),
     xlab="False Positive Rate", ylab="True Positive Rate")
l_ply(rocs_all[2], plot, col=c("blue"), add=T)
l_ply(rocs_all[3], plot, col=c("salmon"), add=T)
axis(1,seq(1,0,0.2))

legend("bottomright",title="Class ROC", legend=c("Floating","Submergent","Emergent"),
       col=c("green","blue","salmon"), lwd=1)

aucs <- laply(rocs_all, function(x) as.numeric(x$auc))
summary(aucs)
mean(aucs) #all are one-versus-others, so take mean to get overall model performance 

text(x=1.000000e-01,y=0.5,"Mean AUC = 0.90")

eers <- laply(rocs_all, function(x) x$sensitivities[which.min(abs(x$sensitivities-x$specificities))])
summary(eers)

conf <- ldply(unique(model_all_lumped$pred$obs), function(cls) data.frame(tp=sum(model_all_lumped$pred$pred[model_all_lumped$pred$obs==cls]==cls), 
                                                           fn=sum(model_all_lumped$pred$pred[model_all_lumped$pred$obs==cls]!=cls), 
                                                           tn=sum(model_all_lumped$pred$pred[model_all_lumped$pred$obs!=cls]!=cls),
                                                           fp=sum(model_all_lumped$pred$pred[model_all_lumped$pred$obs!=cls]==cls)))
conf$tpr <- with(conf, tp/(tp+fn))
conf$tnr <- with(conf, tn/(tn+fp))
summary(conf[,5:6])

####----aug only bothdsm -----######
all_vars_aug <- subset(aug_dsm,
                       select=-c(X.1,X,Y,Green7_12,Red7_12,Red.Edge7_12,NIR7_12,NDVI7_12,DSM7_12, StdVegHeight7_12))
all_vars_aug <- drop_na(all_vars_aug)

tunegrid <- expand.grid(.mtry = (1:7)) 

model_aug <- train(make.names(train)~., # Class is a function of the variables we decided to include
                   data = all_vars_aug, # Use the train data frame as the training data
                   method = 'rf',
                   tuneGrid = tunegrid, # Use the 'random forest' algorithm
                   trControl = trainControl(method = 'repeatedcv', # Use cross-validation
                                            number = 5,
                                            repeats=50,
                                            classProbs = T,
                                            savePredictions=T)) 
model_aug

confusionMatrix.train(model_aug)

tiff(filename="FigS9.tif",width=8,height=6,units="in",res=600)
plot(varImp(model_aug,scale=F),xlab="Importance (% Increase MSE)",
     ylab="Bands")
dev.off()

#alternative Error Matrix 
sub_aug <- subset(model_aug$pred,model_aug$pred$mtry==model_aug$bestTune$mtry)
caret::confusionMatrix(table(sub_aug$pred,sub_aug$obs))

#for precision, recall, and f1 score 
caret::confusionMatrix(table(sub_aug$pred,sub_aug$obs), mode="prec_recall")


#create own df of confm values
Floating <- c(4317,276,407)
SubVeg <- c(947,2090,263)
Emergent <- c(366,45,3689)


#Ref <- c("Floating","SubVeg","Emergent")
conf_df <- data.frame(Floating,SubVeg,Emergent)
con_matx <- as.matrix(conf_df)
con_diag <- diag(con_matx)
n <- sum(con_matx)

avg_OA <- sum(con_diag) / n
avg_OA
#81.3

# observed (true) cases per class
rowsums <- apply(con_matx, 1, sum)

# predicted cases per class
colsums <- apply(con_matx, 2, sum)

# Producer accuracy
avg_PA <- con_diag / colsums
# User accuracy
avg_UA <- con_diag / rowsums
outAcc_aug <- data.frame(producerAccuracy = avg_PA, userAccuracy = avg_UA)
outAcc_aug

rocs_aug <- llply(unique(model_aug$pred$obs), function(cls) {
  roc(response = model_aug$pred$obs==cls, 
      predictor = model_aug$pred[,as.character(cls)])
})

plot(rocs_aug[[1]], col="green", legacy.axes=TRUE, xlim=c(1,0), ylim=c(0,1),
     xlab="False Positive Rate", ylab="True Positive Rate")
l_ply(rocs_aug[2], plot, col=c("blue"), add=T)
l_ply(rocs_aug[3], plot, col=c("salmon"), add=T)
axis(1,seq(1,0,0.2))

legend("bottomright",title="Class ROC", legend=c("Floating","Submergent","Emergent"),
       col=c("green","blue","salmon"), lwd=1)

aucs <- laply(rocs_aug, function(x) as.numeric(x$auc))
summary(aucs)
mean(aucs) #all are one-versus-others, so take mean to get overall model performance 

text(x=1.000000e-01,y=0.5,"Mean AUC = 0.88")

eers <- laply(rocs_aug, function(x) x$sensitivities[which.min(abs(x$sensitivities-x$specificities))])
summary(eers)

conf <- ldply(unique(model_aug$pred$obs), function(cls) data.frame(tp=sum(model_aug$pred$pred[model_aug$pred$obs==cls]==cls), 
                                                                          fn=sum(model_aug$pred$pred[model_aug$pred$obs==cls]!=cls), 
                                                                          tn=sum(model_aug$pred$pred[model_aug$pred$obs!=cls]!=cls),
                                                                          fp=sum(model_aug$pred$pred[model_aug$pred$obs!=cls]==cls)))
conf$tpr <- with(conf, tp/(tp+fn))
conf$tnr <- with(conf, tn/(tn+fp))
summary(conf[,5:6])

####----july only bothdsm -----######
all_vars_july <- subset(all_vars_lumped,
                       select=-c(Green8_5,Red8_5,Red.Edge8_5,NIR8_5,NDVI8_5,DSM8_5, StdVegHeight8_5))
all_vars_july <- drop_na(all_vars_july)

model_july <- train(make.names(train)~., # Class is a function of the variables we decided to include
                   data = all_vars_july, # Use the train data frame as the training data
                   method = 'rf',# Use the 'random forest' algorithm
                   trControl = trainControl(method = 'repeatedcv', # Use cross-validation
                                            number = 5,
                                            repeats=50,
                                            classProbs = T,
                                            savePredictions=T)) 
model_july

confusionMatrix.train(model_july)
plot(varImp(model_july,scale=F))

#alternative Error Matrix 
sub_july <- subset(model_july$pred,model_july$pred$mtry==model_july$bestTune$mtry)
caret::confusionMatrix(table(sub_july$pred,sub_july$obs))

#for precision, recall, and f1 score 
caret::confusionMatrix(table(sub_july$pred,sub_july$obs), mode="prec_recall")


#create own df of confm values
Floating <- c(4264,339,397)
SubVeg <- c(917,1891,292)
Emergent <- c(249,138,3763)


#Ref <- c("Floating","SubVeg","Emergent")
conf_df <- data.frame(Floating,SubVeg,Emergent)
con_matx <- as.matrix(conf_df)
con_diag <- diag(con_matx)
n <- sum(con_matx)

avg_OA <- sum(con_diag) / n
avg_OA
#81.0

# observed (true) cases per class
rowsums <- apply(con_matx, 1, sum)

# predicted cases per class
colsums <- apply(con_matx, 2, sum)

# Producer accuracy
avg_PA <- con_diag / colsums
# User accuracy
avg_UA <- con_diag / rowsums
outAcc_july <- data.frame(producerAccuracy = avg_PA, userAccuracy = avg_UA)
outAcc_july

rocs_july <- llply(unique(model_july$pred$obs), function(cls) {
  roc(response = model_july$pred$obs==cls, 
      predictor = model_july$pred[,as.character(cls)])
})

plot(rocs_july[[1]], col="green", legacy.axes=TRUE, xlim=c(1,0), ylim=c(0,1),
     xlab="False Positive Rate", ylab="True Positive Rate")
l_ply(rocs_july[2], plot, col=c("blue"), add=T)
l_ply(rocs_july[3], plot, col=c("salmon"), add=T)
axis(1,seq(1,0,0.2))

legend("bottomright",title="Class ROC", legend=c("Floating","Submergent","Emergent"),
       col=c("green","blue","salmon"), lwd=1)

aucs <- laply(rocs_july, function(x) as.numeric(x$auc))
summary(aucs)
mean(aucs) #all are one-versus-others, so take mean to get overall model performance 

text(x=1.000000e-01,y=0.5,"Mean AUC = 0.88")

eers <- laply(rocs_july, function(x) x$sensitivities[which.min(abs(x$sensitivities-x$specificities))])
summary(eers)

conf <- ldply(unique(model_july$pred$obs), function(cls) data.frame(tp=sum(model_july$pred$pred[model_july$pred$obs==cls]==cls), 
                                                                   fn=sum(model_july$pred$pred[model_july$pred$obs==cls]!=cls), 
                                                                   tn=sum(model_july$pred$pred[model_july$pred$obs!=cls]!=cls),
                                                                   fp=sum(model_july$pred$pred[model_july$pred$obs!=cls]==cls)))
conf$tpr <- with(conf, tp/(tp+fn))
conf$tnr <- with(conf, tn/(tn+fp))
summary(conf[,5:6])

####----multidate spectral only----##### 
all_vars_so <- subset(all_vars_lumped,
                      select=-c(DSM8_5,StdVegHeight8_5,DSM7_12,StdVegHeight7_12))
all_vars_so <- drop_na(all_vars_so)

model_all_so <- train(make.names(train)~., # Class is a function of the variables we decided to include
                      data = all_vars_so, # Use the train data frame as the training data
                      method = 'rf',# Use the 'random forest' algorithm
                      trControl = trainControl(method = 'repeatedcv', # Use cross-validation
                                               number = 5,
                                               repeats=50,
                                               classProbs=T,
                                               savePredictions=T)) 
model_all_so

confusionMatrix.train(model_all_so)
plot(varImp(model_all_so,scale=F))

#alternative Error Matrix 
sub_all_so <- subset(model_all_so$pred,model_all_so$pred$mtry==model_all_so$bestTune$mtry)
caret::confusionMatrix(table(sub_all_so$pred,sub_all_so$obs))

#for precision, recall, and f1 score 
caret::confusionMatrix(table(sub_all_so$pred,sub_all_so$obs), mode="prec_recall")

#create own df of confm values
Floating <- c(4241,356,403)
SubVeg <- c(891,1918,291)
Emergent <- c(167,70,3913)


#Ref <- c("Floating","SubVeg","Emergent")
conf_df <- data.frame(Floating,SubVeg,Emergent)
con_matx <- as.matrix(conf_df)
con_diag <- diag(con_matx)
n <- sum(con_matx)

avg_OA <- sum(con_diag) / n
avg_OA
# 82.2%

# observed (true) cases per class
rowsums <- apply(con_matx, 1, sum)

# predicted cases per class
colsums <- apply(con_matx, 2, sum)

# Producer accuracy
avg_PA <- con_diag / colsums
# User accuracy
avg_UA <- con_diag / rowsums
outAcc_all_so <- data.frame(producerAccuracy = avg_PA, userAccuracy = avg_UA)
outAcc_all_so

rocs_so <- llply(unique(model_all_so$pred$obs), function(cls) {
  roc(response = model_all_so$pred$obs==cls, 
      predictor = model_all_so$pred[,as.character(cls)])
})

plot(rocs_so[[1]], col="green", legacy.axes=TRUE, xlim=c(1,0), ylim=c(0,1),
     xlab="False Positive Rate", ylab="True Positive Rate")
l_ply(rocs_so[2], plot, col=c("blue"), add=T)
l_ply(rocs_so[3], plot, col=c("salmon"), add=T)
axis(1,seq(1,0,0.2))

legend("bottomright",title="Class ROC", legend=c("Floating","Submergent","Emergent"),
       col=c("green","blue","salmon"), lwd=1)

aucs <- laply(rocs_so, function(x) as.numeric(x$auc))
summary(aucs)
mean(aucs) #all are one-versus-others, so take mean to get overall model performance 

text(x=1.000000e-01,y=0.5,"Mean AUC = 0.9")

eers <- laply(rocs_so, function(x) x$sensitivities[which.min(abs(x$sensitivities-x$specificities))])
summary(eers)

conf <- ldply(unique(model_all_so$pred$obs), function(cls) data.frame(tp=sum(model_all_so$pred$pred[model_all_so$pred$obs==cls]==cls), 
                                                                   fn=sum(model_all_so$pred$pred[model_all_so$pred$obs==cls]!=cls), 
                                                                   tn=sum(model_all_so$pred$pred[model_all_so$pred$obs!=cls]!=cls),
                                                                   fp=sum(model_all_so$pred$pred[model_all_so$pred$obs!=cls]==cls)))
conf$tpr <- with(conf, tp/(tp+fn))
conf$tnr <- with(conf, tn/(tn+fp))
summary(conf[,5:6])


####----august spectral only----##### 
aug_so <- subset(aug_dsm,
                      select=-c(X.1,X,Y,DSM8_5,StdVegHeight8_5,Green7_12,Red7_12,Red.Edge7_12,NIR7_12,NDVI7_12,DSM7_12,StdVegHeight7_12))
aug_so <- drop_na(aug_so)

tunegrid <- expand.grid(.mtry = (1:5)) 

model_aug_so <- train(make.names(train)~., # Class is a function of the variables we decided to include
                      data = aug_so, # Use the train data frame as the training data
                      method = 'rf',
                      tuneGrid = tunegrid, # Use the 'random forest' algorithm
                      trControl = trainControl(method = 'repeatedcv', # Use cross-validation
                                               number = 5,
                                               repeats=50,
                                               classProbs=T,
                                               savePredictions=T)) 
model_aug_so

confusionMatrix.train(model_aug_so)

tiff(filename="FigS11.tif",width=8,height=6,units="in",res=600)
plot(varImp(model_aug_so,scale=F),xlab="Importance (% Increase MSE)",
     ylab="Bands")
dev.off()


#alternative Error Matrix 
sub_aug_so <- subset(model_aug_so$pred,model_aug_so$pred$mtry==model_aug_so$bestTune$mtry)
caret::confusionMatrix(table(sub_aug_so$pred,sub_aug_so$obs))

#for precision, recall, and f1 score 
caret::confusionMatrix(table(sub_aug_so$pred,sub_aug_so$obs), mode="prec_recall")

#create own df of confm values
Floating <- c(4245,374,381)
SubVeg <- c(964,2102,234)
Emergent <- c(318,108,3674)


#Ref <- c("Floating","SubVeg","Emergent")
conf_df <- data.frame(Floating,SubVeg,Emergent)
con_matx <- as.matrix(conf_df)
con_diag <- diag(con_matx)
n <- sum(con_matx)

avg_OA <- sum(con_diag) / n
avg_OA
# 80.8%

# observed (true) cases per class
rowsums <- apply(con_matx, 1, sum)

# predicted cases per class
colsums <- apply(con_matx, 2, sum)

# Producer accuracy
avg_PA <- con_diag / colsums
# User accuracy
avg_UA <- con_diag / rowsums
outAcc_aug_so <- data.frame(producerAccuracy = avg_PA, userAccuracy = avg_UA)
outAcc_aug_so

rocs_so <- llply(unique(model_aug_so$pred$obs), function(cls) {
  roc(response = model_aug_so$pred$obs==cls, 
      predictor = model_aug_so$pred[,as.character(cls)])
})

plot(rocs_so[[1]], col="green", legacy.axes=TRUE, xlim=c(1,0), ylim=c(0,1),
     xlab="False Positive Rate", ylab="True Positive Rate")
l_ply(rocs_so[2], plot, col=c("blue"), add=T)
l_ply(rocs_so[3], plot, col=c("salmon"), add=T)
axis(1,seq(1,0,0.2))

legend("bottomright",title="Class ROC", legend=c("Floating","Submergent","Emergent"),
       col=c("green","blue","salmon"), lwd=1)

aucs <- laply(rocs_so, function(x) as.numeric(x$auc))
summary(aucs)
mean(aucs) #all are one-versus-others, so take mean to get overall model performance 

text(x=1.000000e-01,y=0.5,"Mean AUC = 0.9")

eers <- laply(rocs_so, function(x) x$sensitivities[which.min(abs(x$sensitivities-x$specificities))])
summary(eers)

conf <- ldply(unique(model_aug_so$pred$obs), function(cls) data.frame(tp=sum(model_aug_so$pred$pred[model_aug_so$pred$obs==cls]==cls), 
                                                                      fn=sum(model_aug_so$pred$pred[model_aug_so$pred$obs==cls]!=cls), 
                                                                      tn=sum(model_aug_so$pred$pred[model_aug_so$pred$obs!=cls]!=cls),
                                                                      fp=sum(model_aug_so$pred$pred[model_aug_so$pred$obs!=cls]==cls)))
conf$tpr <- with(conf, tp/(tp+fn))
conf$tnr <- with(conf, tn/(tn+fp))
summary(conf[,5:6])



alpenagt_all <- alpenagt_all %>% mutate(train=as.double(train),
                                        train= if_else(train== 2,1,train),
                                        train=as.double(train)) #floating 
alpenagt_all <- alpenagt_all%>% mutate(train=as.double(train),
                                       train= if_else(train== 3,2,train),
                                       train=as.double(train)) #submergent
alpenagt_all <- alpenagt_all%>% mutate(train=as.double(train),
                                       train= if_else(train== 4,3,train),
                                       train=as.double(train)) #emergent 

nrow(alpenagt_all[which(alpenagt_all$train==1),]) #159 efb/floating 
nrow(alpenagt_all[which(alpenagt_all$train==2),]) #62 sub 
nrow(alpenagt_all[which(alpenagt_all$train==3),]) #83 emergent



#take out a good random portion of floating, way too many
(floatingindex <- which(alpenagt_all$train == '1'))
(deleteindex <- sample(floatingindex, length(floatingindex) - 100)) 
alpenagt_all <- alpenagt_all[-deleteindex, ]

#check if works 
nrow(alpenagt_all[which(alpenagt_all$train==1),]) #100 floating

####----new plots with lumped classes----#### 
aug_prof_features<- reshape2::melt(alpenagt_all,id=c("train"),
                                   measure=c("Green8_5","Red8_5","Red.Edge8_5","NIR8_5","NDVI8_5","DSM8_5","StdVegHeight8_5"))
names(aug_prof_features)[3] <- "Reflectance"

#average 
aug_avg<-reshape2::dcast(aug_prof_features,train+variable~.,fun=mean, value.var="Reflectance",na.rm=TRUE)
names(aug_avg)[3]<-"mean_ref"
#se
aug_se<-reshape2::dcast(aug_prof_features, train+variable~., fun=sd, value.var="Reflectance",na.rm=TRUE)
names(aug_se)[3]<-"sd"
table(alpenagt_all$train)#order for below command 
#first reps is train sample n, second is reps are for band #'s (4 classes with varying n, now 7 bands in this aug data)
#numerical order, (efb,yv,sv,emergent)
Class_length<-c(rep(100,7),rep(66,7),rep(82,7))
as.data.frame(Class_length)
aug_se$class_n<-Class_length
#jenky way, but works  
aug_se$se <- aug_se$sd/sqrt(aug_se$class_n)

aug_spectra <- merge(aug_avg, aug_se, by=c("train","variable"), na.rm=TRUE)
aug_spectra$lowerse <- aug_spectra$mean_ref-aug_spectra$se
aug_spectra$upperse <- aug_spectra$mean_ref+aug_spectra$se

aug_spectra$train <- as.factor(aug_spectra$train)
aug_all <- aug_spectra %>% mutate(train=as.character(train),
                                  train= if_else(train == 1,"Floating Veg.",train),
                                  train=as.factor(train))  %>% mutate(train=as.character(train),
                                                                      train= if_else(train == 2,"Sub.Veg.",train),
                                                                      train=as.factor(train)) %>% mutate(train=as.character(train),
                                                                                                         train= if_else(train == 3,"Emergent Veg.",train),
                                                                                                         train=as.factor(train)) 



#Actual spectral plot 
aug_spectra_only <- subset(aug_all, variable !="DSM8_5") 
aug_spectra_only <- subset(aug_spectra_only, variable !="StdVegHeight8_5") 


aug_spectra_lumped <- ggplot(aug_spectra_only,aes(x=variable,y=mean_ref,group=train,color=train)) +
  geom_line(size=1.25)+
  geom_errorbar(aes(ymin=lowerse,ymax=upperse,width=0.05))+
  ylab("% Reflectance")+
  ggtitle("Spectral Profiles, 8/2019")+
  theme_bw()+
  theme(legend.position='none')+
  theme(axis.title.x = element_blank()) +
  # scale_color_discrete(breaks=c("Frogbit","Lily","SV","Typha + EFB/Litter"),
  #                      labels=c("Frogbit","Lily","SV","Typha + EFB/Litter")) +
  scale_x_discrete(labels=c("Green","Red","Red Edge","NIR","NDVI"))

#rug only 
aug_rugosity <- subset(aug_all, variable == "StdVegHeight8_5")
aug_rug_lumped <- ggplot(aug_rugosity,aes(x=variable,y=mean_ref,fill=train)) +
  geom_bar(position="dodge",stat="identity") +
  geom_errorbar(aes(ymin=lowerse,ymax=upperse,width=0.2), 
                position=position_dodge(0.9)) +
  ylab("Rugosity (Std Dev. of Surface Height)")+
  ggtitle("Textural Profiles, 8/2019")+
  theme_bw()+
  theme(legend.position='none',
        axis.title.x = element_blank()) +
  # scale_color_discrete(breaks=c("Frogbit","Lily","SV","Typha + EFB/Litter"),
  #                      labels=c("Frogbit","Lily","SV","Typha + EFB/Litter")) +
  scale_x_discrete(labels=c("Rugosity"))

#surface height only 
aug_surface <- subset(aug_all, variable == "DSM8_5")
aug_surf_lumped <- ggplot(aug_surface,aes(x=variable,y=mean_ref,fill=train)) +
  geom_bar(position="dodge",stat="identity") +
  geom_errorbar(aes(ymin=lowerse,ymax=upperse,width=0.2), 
                position=position_dodge(0.9)) +
  ylab("Surface Height (m)")+
  ggtitle("Surface Heights, 8/2019")+
  theme_bw()+
  theme(axis.title.x = element_blank()) + 
  labs(fill="Veg. Class")+
  coord_cartesian(ylim=c(178.7,179.4)) +
  # scale_color_discrete(breaks=c("Frogbit","Lily","SV","Typha + EFB/Litter"),
  #                      labels=c("Frogbit","Lily","SV","Typha + EFB/Litter")) +
  scale_x_discrete(labels=c("DSM"))

aug_grid <- grid.arrange(aug_spectra_lumped,aug_rug_lumped,
                         aug_surf_lumped,ncol=3)
ggsave("FigureS2.jpg",width=9,height=6,units=c("in"),dpi=300)


#july
july_prof_features <- reshape2::melt(alpenagt_all,id=c("train"),
                                     measure=c("Green7_12","Red7_12","Red.Edge7_12","NIR7_12","NDVI7_12","DSM7_12","StdVegHeight7_12")) 
names(july_prof_features)[3] <- "Reflectance" 
#average 
july_avg<-reshape2::dcast(july_prof_features,train+variable~.,fun=mean, value.var="Reflectance",na.rm=TRUE)
names(july_avg)[3]<-"mean_ref"
#se
july_se<-reshape2::dcast(july_prof_features, train+variable~., fun=sd, value.var="Reflectance",na.rm=TRUE)
names(july_se)[3]<-"sd"
table(alpenagt_all$train)
#first reps is class sample n, second is reps are for band #'s (4 classes with varying n, 7 bands for each)
Class_length<-c(rep(100,7),rep(66,7),rep(82,7))
as.data.frame(Class_length)
july_se$class_n<-Class_length
#jenky way, but works  
july_se$se <- july_se$sd/sqrt(july_se$class_n)

july_spectra<-merge(july_avg, july_se, by=c("train","variable"), na.rm=TRUE)
july_spectra$lowerse<-july_spectra$mean_ref-july_spectra$se
july_spectra$upperse<-july_spectra$mean_ref+july_spectra$se

july_spectra$train <- as.factor(july_spectra$train)

july_spectra$train <- as.factor(july_spectra$train)
july_spectra_all <- july_spectra %>% mutate(train=as.character(train),
                                            train= if_else(train == 1,"Floating Veg.",train),
                                            train=as.factor(train))  %>% mutate(train=as.character(train),
                                                                                train= if_else(train == 2,"Sub.Veg.",train),
                                                                                train=as.factor(train)) %>% mutate(train=as.character(train),
                                                                                                                   train= if_else(train == 3,"Emergent Veg.",train),
                                                                                                                   train=as.factor(train))
#specral plot 
jul_spectra_only <- subset(july_spectra_all, variable !="StdVegHeight7_12")
jul_spectra_only <- subset(jul_spectra_only, variable !="DSM7_12")

july_spectra <- ggplot(jul_spectra_only,aes(x=variable,y=mean_ref,group=train,color=train))+
  geom_line(size=1.25)+
  geom_errorbar(aes(ymin=lowerse,ymax=upperse,width=0.05))+
  ylab("% Reflectance")+
  ggtitle("Spectral Profiles, 7/2019")+
  theme_bw()+
  theme(legend.position='none', 
        axis.title.x = element_blank())+
  # scale_color_discrete(breaks=c("Frogbit","Lily","SV","Typha + EFB/Litter"),
  #                      labels=c("Frogbit","Lily","SV","Typha + EFB/Litter")) +
  scale_x_discrete(labels=c("Green","Red","Red Edge","NIR","NDVI"))

#rug only 
july_rugosity <- subset(july_spectra_all, variable == "StdVegHeight7_12")
july_rug_lumped <- ggplot(july_rugosity,aes(x=variable,y=mean_ref,fill=train)) +
  geom_bar(position="dodge",stat="identity") +
  geom_errorbar(aes(ymin=lowerse,ymax=upperse,width=0.2), 
                position=position_dodge(0.9)) +
  ylab("Rugosity (Std Dev. of Surface Height)")+
  ggtitle("Textural Profiles, 7/2019")+
  theme_bw()+
  theme(legend.position='none',
        axis.title.x = element_blank()) +
  # scale_color_discrete(breaks=c("Frogbit","Lily","SV","Typha + EFB/Litter"),
  #                      labels=c("Frogbit","Lily","SV","Typha + EFB/Litter")) +
  scale_x_discrete(labels=c("Rugosity"))

#dsm 
july_surface <- subset(july_spectra_all,variable == "DSM7_12")
july_surf <- ggplot(july_surface,aes(x=variable,y=mean_ref,group=train,fill=train))+
  geom_bar(position="dodge",stat="identity") +
  geom_errorbar(aes(ymin=lowerse,ymax=upperse,width=0.2), 
                position=position_dodge(0.9)) +
  ylab("Surface Height (m)")+
  ggtitle("Surface Heights, 7/2019")+
  theme_bw()+
  theme(axis.title.x = element_blank()) + 
  labs(fill="Veg. Class")+
  coord_cartesian(ylim=c(178.75,179.4)) +
  # scale_color_discrete(breaks=c("Frogbit","Lily","SV","Typha + EFB/Litter"),
  #                      labels=c("Frogbit","Lily","SV","Typha + EFB/Litter")) +
  scale_x_discrete(labels=c("DSM"))

# #plot side by side 
july_grid <- grid.arrange(july_spectra, july_rug_lumped, july_surf, ncol=3)

#altogether now 
grid_all <- grid.arrange(aug_grid, july_grid,nrow=2)



#####----aug chm lumped----#####

#again mutate() suddenly not working,gotta love R 
# alpenagt_all_chm <- alpenagt_all_chm %>% mutate(train=as.double(train),
#                                         train= if_else(train== 2,1,train),
#                                         train=as.double(train)) #floating 
alpenagt_all_chm$train <- if_else(alpenagt_all_chm$train == 2, 1, 
                                  as.double(alpenagt_all_chm$train))

# alpenagt_all_chm <- alpenagt_all_chm%>% mutate(train=as.double(train),
#                                        train= if_else(train== 3,2,train),
#                                        train=as.double(train)) #submergent
alpenagt_all_chm$train <- if_else(alpenagt_all_chm$train == 3, 2, 
                                  as.double(alpenagt_all_chm$train))

# alpenagt_all_chm <- alpenagt_all_chm%>% mutate(train=as.double(train),
#                                        train= if_else(train== 4,3,train),
#                                        train=as.double(train)) #emergent 
alpenagt_all_chm$train <- if_else(alpenagt_all_chm$train == 4, 3, 
                                  as.double(alpenagt_all_chm$train))

nrow(alpenagt_all_chm[which(alpenagt_all_chm$train==1),]) #153 efb/floating 
nrow(alpenagt_all_chm[which(alpenagt_all_chm$train==2),]) #51 sub 
nrow(alpenagt_all_chm[which(alpenagt_all_chm$train==3),]) #81 emergent

#take out a good random portion of floating, way too many
(floatingindex <- which(alpenagt_all_chm$train == '1'))
(deleteindex <- sample(floatingindex, length(floatingindex) - 100)) 
alpenagt_all_chm <- alpenagt_all_chm[-deleteindex, ]

#check if works 
nrow(alpenagt_all_chm[which(alpenagt_all_chm$train==1),]) #100 floating

all_vars_lumped <- subset(alpenagt_all_chm,
                          select=-c(X,Y))

all_vars_lumped <- drop_na(all_vars_lumped)

tunegrid <- expand.grid(.mtry = (1:7)) 

model_aug_chm <- train(make.names(train)~., # Class is a function of the variables we decided to include
                      data = all_vars_lumped, # Use the train data frame as the training data
                      method = 'rf',
                      tuneGrid=tunegrid,# Use the 'random forest' algorithm
                      trControl = trainControl(method = 'repeatedcv', # Use cross-validation
                                               number = 5,
                                               repeats = 50,
                                               classProbs = T, #need to include for AUC/ROC, took out as.factor() on train 
                                               savePredictions = T)) #get warning with these params
model_aug_chm


tiff(filename="FigS12.tif",width=8,height=6,units="in",res=600)
plot(varImp(model_aug_chm,scale=F),xlab="Importance (% Increase MSE)",
     ylab="Bands")
dev.off()

# model_all_lumped <- train(as.factor(train)~., # Class is a function of the variables we decided to include
#                           data = all_vars_lumped, # Use the train data frame as the training data
#                           method = 'rf',# Use the 'random forest' algorithm
#                           trControl = trainControl(method = 'repeatedcv', # Use cross-validation
#                                                    number = 5,
#                                                    repeats = 50,
#                                                    classProbs = T, #need to include for AUC/ROC, took out as.factor() on train 
#                                                    savePredictions = T,
#                                                    verboseIter = T)) #get warning with these params
# model_all_lumped

#try this way 
#levels(all_vars_lumped$train) <- c("Floating","Submerged","Emergent")

control <- trainControl(method='repeatedcv', 
                        number=5, 
                        repeats=50, 
                        search='grid',
                        classProbs = T, #need to include for AUC/ROC, took out as.factor() on train 
                        savePredictions = T)

tunegrid <- expand.grid(.mtry = (1:7)) 

rf_gridsearch <- train(make.names(train)~., 
                       data = all_vars_lumped,
                       method = 'rf',
                       metric = 'Accuracy',
                       tuneGrid = tunegrid)
print(rf_gridsearch)
plot(rf_gridsearch)
#actually maybe I should do it this way, I get best result with mtry of 1.... 

library(e1071)
bestMtry <- tuneRF(all_vars_lumped[2:8],all_vars_lumped$train, stepFactor = 1.5, improve = 1e-5, ntree = 500)

model_all_lumped <- train(make.names(train)~., # Class is a function of the variables we decided to include
                          data = all_vars_lumped, # Use the train data frame as the training data
                          method = 'rf',
                          tuneGrid = tunegrid, # Use the 'random forest' algorithm
                          trControl = trainControl(method = 'repeatedcv', # Use cross-validation
                                                   number = 5,
                                                   repeats = 50,
                                                   classProbs = T, #need to include for AUC/ROC, took out as.factor() on train 
                                                   savePredictions = T)) #get warning with these params
model_all_lumped

confusionMatrix.train(model_all_lumped)
plot(varImp(model_all_lumped,scale=F))

#alternative Error Matrix 
sub_all_lumped <- subset(model_all_lumped$pred,model_all_lumped$pred$mtry==model_all_lumped$bestTune$mtry)
caret::confusionMatrix(table(sub_all_lumped$pred,sub_all_lumped$obs))

#good explanation of performance statistics here 
# https://topepo.github.io/caret/measuring-performance.html#measures-for-predicted-classes

#for precision, recall, and f1 score 
caret::confusionMatrix(table(sub_all_lumped$pred,sub_all_lumped$obs), mode="prec_recall")

#create own df of confm values
Floating <- c(4091,329,580)
SubVeg <- c(932,1226,392)
Emergent <- c(384,64,3602)


#Ref <- c("Floating","SubVeg","Emergent")
conf_df <- data.frame(Floating,SubVeg,Emergent)
con_matx <- as.matrix(conf_df)
con_diag <- diag(con_matx)
n <- sum(con_matx)

avg_OA <- sum(con_diag) / n
avg_OA
#76.9 %. 

# observed (true) cases per class
rowsums <- apply(con_matx, 1, sum)

# predicted cases per class
colsums <- apply(con_matx, 2, sum)

# Producer accuracy
avg_PA <- con_diag / colsums
# User accuracy
avg_UA <- con_diag / rowsums
outAcc_all_lumped <- data.frame(producerAccuracy = avg_PA, userAccuracy = avg_UA)
outAcc_all_lumped

#ROC and AUC 
library(pROC)
library(plyr)

rocs_all <- llply(unique(model_all_lumped$pred$obs), function(cls) {
  roc(response = model_all_lumped$pred$obs==cls, 
      predictor = model_all_lumped$pred[,as.character(cls)])
})

plot(rocs_all[[1]], col="green", legacy.axes=TRUE, xlim=c(1,0), ylim=c(0,1),
     xlab="False Positive Rate", ylab="True Positive Rate")
l_ply(rocs_all[2], plot, col=c("blue"), add=T)
l_ply(rocs_all[3], plot, col=c("salmon"), add=T)
axis(1,seq(1,0,0.2))

legend("bottomright",title="Class ROC", legend=c("Floating","Submergent","Emergent"),
       col=c("green","blue","salmon"), lwd=1)

aucs <- laply(rocs_all, function(x) as.numeric(x$auc))
summary(aucs)
mean(aucs) #all are one-versus-others, so take mean to get overall model performance 

text(x=1.000000e-01,y=0.5,"Mean AUC = 0.90")

eers <- laply(rocs_all, function(x) x$sensitivities[which.min(abs(x$sensitivities-x$specificities))])
summary(eers)

conf <- ldply(unique(model_all_lumped$pred$obs), function(cls) data.frame(tp=sum(model_all_lumped$pred$pred[model_all_lumped$pred$obs==cls]==cls), 
                                                                          fn=sum(model_all_lumped$pred$pred[model_all_lumped$pred$obs==cls]!=cls), 
                                                                          tn=sum(model_all_lumped$pred$pred[model_all_lumped$pred$obs!=cls]!=cls),
                                                                          fp=sum(model_all_lumped$pred$pred[model_all_lumped$pred$obs!=cls]==cls)))
conf$tpr <- with(conf, tp/(tp+fn))
conf$tnr <- with(conf, tn/(tn+fp))
summary(conf[,5:6])



####----classified raster----####
# rasterPath <- raster("Z:/UAV2019/Alpena/Composites11cm/AlpenaBothDSM_11cmMasked.tif")
# plot(rasterPath)
#open raster and check n bands
r<-stack("Z:/UAV2019/Alpena/Composites11cm/MoreGT11cm/AlpenaBothDSM_11cmMasked.tif")
r
r@layers
object.size(r)
#plot NDVI (even with weird -10000 value form pix4D?)
hist(r[[5]]) #renders even though I cannot see values for each layer, weird masking issue in Arc 

# make sp df
coordinates(alpenagt_all) <- ~X + Y
proj4string(alpenagt_all) <- "+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs"

rbPal <- colorRampPalette(c("orange","blue","green"))
alpenagt_all$col <- rbPal(3)[as.numeric(cut(alpenagt_all$train,breaks=3))]

pal <- brewer.pal(n=6,name="PiYG")
cuts <- c(-0.25,0,0.25,0.5,0.75,1)

x1 <- c()

plot(r[[5]], breaks=cuts,
     col=pal,
     xlab="UTM Easting (m)",
     ylab="UTM Northing (m)",
     legend.args=list(text="NDVI"))
points(alpenagt_all,
       pch=19,
       cex=0.75,
       col=alpenagt_all$col)
legend(x=306900, y=4994440,c("Emergent Veg.","Floating Veg.","Submergent Veg."),pch=19,cex=0.75,col=alpenagt_all$col)

r_brick <- brick(r)
hist(r_brick[[5]])
#for nice looking figure
aug_ndvi <- r_brick[[5]]
#must mask out werid neg values in the thousands  (just white areas in between wetland veg?)
pal <- brewer.pal(n=6,name="PiYG")
cuts <- c(-1,0,0.25,0.5,0.75,1)
plot(r_brick[[5]],breaks=cuts,col=pal)

ndvi_df <- rasterToPoints(r[[5]])
ndvi_df <- data.frame(ndvi_df)

ggplot() +
  geom_raster(data = ndvi_df,
              aes(x = x, y = y, alpha = AlpenaBothDSM_11cmMasked.5)) + 
  coord_equal()

####-----rf classification from package----##### 
#rename raster bands to match those in classification model 
names(r) <- c("Green8_5","Red8_5","Red.Edge8_5","NIR8_5","NDVI8_5","DSM8_5","StdVegHeight8_5",
              "Green7_12","Red7_12","Red.Edge7_12","NIR7_12","NDVI7_12","DSM7_12","StdVegHeight7_12")

#k fold split with extracted values and classes 
j <- kfold(alpenagt_all, k = 5, by=alpenagt_all$train)
table(j)

x <- list() 
for (k in 1:5) {
  train <- alpenagt_all[j!= k,]
  test <- alpenagt_all[j == k,]
  rf <- randomForest(x=train[,4:17], y=factor(train[,"train"]),ntree=501,proximity = TRUE, importance=TRUE)
  pclass <- predict(rf, test, type='class')
  # create a data.frame using the reference and prediction
  x[[k]] <- cbind(test$Class, as.integer(pclass))}
#nested for loop somehow for repeated cv?  

#obtain variable importance of rf object
importance(rf)

#create classified raster based on rf predictions (from randomForest object in for loop)
# classified_rast <- predict(r,rf,filename="ClassPred.img",type="response",
#                            na.rm=TRUE,overwrite=TRUE,progress="window")

####-----classified raster from caret model-----#####
#call in raster of model predictions from previous session 
classified_rast <- raster("Z:/UAV2019/Alpena/Composites11cm/ClassPred.img")


#simple plot
# pal <- colorRampPalette(c("darkolivegreen1","darkolivegreen","darkgreen"))
# cuts <- c(0,1,2,3)
# plot(classified_rast,breaks=cuts,col=pal(3))

##ggplot for relabeling 
r_points <- rasterToPoints(classified_rast)
r_df <- data.frame(r_points)

#convert numbers to class names 
r_df$ClassPred <- if_else(r_df$ClassPred == 1, "Floating Veg.", as.character(r_df$ClassPred))
r_df$ClassPred <- if_else(r_df$ClassPred == 2, "Submergent Veg.", as.character(r_df$ClassPred))
r_df$ClassPred <- if_else(r_df$ClassPred == 3, "Emergent Veg.", as.character(r_df$ClassPred))

r_df$ClassPred <- factor(r_df$ClassPred,
                             levels = c("Emergent Veg.","Floating Veg.","Submergent Veg."))



# r_df <- r_df %>% mutate(ClassPred=as.character(ClassPred),
#                         ClassPred= if_else(ClassPred== 1,'Floating Veg.',ClassPred), 
#                         ClassPred=as.factor(ClassPred))
# r_df <- r_df %>% mutate(ClassPred=as.character(ClassPred),
#                         ClassPred= if_else(ClassPred== 2,'Submergent Veg.',ClassPred), 
#                         ClassPred=as.factor(ClassPred))
# r_df <- r_df %>% mutate(ClassPred=as.character(ClassPred),
#                         ClassPred= if_else(ClassPred== 3,'Emergent Veg.',ClassPred), 
#                         ClassPred=as.factor(ClassPred))

names(r_df$x) <- "long"
names(r_df$y) <- "lat"

library(stringr)
library(ggsn)

class.fig <- ggplot(data=r_df) +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  geom_tile(aes(x=x,y=y,fill=ClassPred)) +
  scale_fill_viridis(discrete = TRUE,direction=-1) + 
  #scale_color_manual(values=c("Floating Veg"="#00B81F","Submergent Veg."="#00B0F6","Emergent Veg."="#F8766D")) + 
  #ggtitle("Vegetation Classification of UAS Imagery")+
  ylab("UTM Northing (m)") +
  xlab("UTM Easting (m)") +
  labs(fill= "Predicted Vegetation Classes") +
  theme(legend.key.height= unit(5,"mm")) + 
  north(scale=0.12, location= "bottomright", symbol = 12, x.min=306775, x.max=306950,
        y.min=4994365, y.max=4994500) +
  scalebar(location = "bottomright", dist=50, dist_unit = "m", 
           x.min=306850, x.max=306900,
           y.min=4994375, y.max=4994525, transform = FALSE) +
  ggplot(data=as.data.frame(rosie_pts)) +
  geom_point(geom_point(aes(x=X,y=Y,color=veg_class))) +
  scale_fill_viridis(discrete = TRUE,direction=-1)
  

grid.arrange(class.fig)

# transform = TRUE, model = "WGS84"
setwd("Z:/UAV2019/Alpena")

rosie_pts <- read.csv("RosiePts_wClassPredictions_V2.csv")
rosie_pts <- drop_na(rosie_pts)

coordinates(rosie_pts) <- ~X + Y
proj4string(rosie_pts) <- "+proj=longlat +datum=WGS84"
#proj4string(alpenagt_all) <- "+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs"
rosie_pts <- spTransform(rosie_pts, CRS("+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs"))

rosie_pts$veg_class <- if_else(rosie_pts$veg_class == "F", "Floating Vegetation", as.character(rosie_pts$veg_class))
rosie_pts$veg_class <- if_else(rosie_pts$veg_class == "MV", "Mixed Vegetation", as.character(rosie_pts$veg_class))
rosie_pts$veg_class <- if_else(rosie_pts$veg_class == "OW", "Open Water", as.character(rosie_pts$veg_class))
rosie_pts$veg_class <- if_else(rosie_pts$veg_class == "S", "Submergent Vegetation", as.character(rosie_pts$veg_class))
rosie_pts$veg_class <- if_else(rosie_pts$veg_class == "T", "Typha x glauca", as.character(rosie_pts$veg_class))

rosie_pts <- as.data.frame(rosie_pts)
r_df <- as.data.frame(r_df)


multispec_11cm <- stack("Z:/UAV2019/Alpena/Composites11cm/AlpenaStackDSM_11cm_Clip.tif")
multispec_11cm <- multispec_11cm[[1:4]]
plotRGB(multispec_11cm,4,1,2,stretch="lin")
#plot(rose_fig,add=T)

ms_df <- as.data.frame(multispec_11cm,xy=TRUE)

e <- extent(multispec_11cm)

library(RStoolbox)
png("Figure3.png",width=900,height=600,res=600)
rose_fig <- ggplot() +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  ggRGB(multispec_11cm,r=4,g=1,b=2,stretch="lin",ggLayer=TRUE) +
  geom_tile(aes(x=x,y=y,fill=ClassPred),data=r_df) +
  scale_fill_viridis(discrete = TRUE,direction=-1,
                     name="RF-Predicted Veg. Classes") + 
  geom_point(aes(x=X,y=Y,color=veg_class),data=rosie_pts,
             color="black", size=2.5, show.legend = FALSE )+ #hacky way to get black borders on pts 
   scale_size( guide = "none" )+
  geom_point(aes(x=X,y=Y,color=veg_class),data=rosie_pts,
             size=2) + #stroke=2 
  scale_color_viridis(discrete=TRUE, option="magma",
                      name="Designated Veg. Zone") +
  ylab("UTM Northing (m)") +
  xlab("UTM Easting (m)") +
  #labs(color= "Field Determined Vegetation Classes", x=NULL,y=NULL) +
  theme(legend.key.height= unit(5,"mm"),
        legend.position = "right",
        legend.background = element_rect(fill = "lightgrey")) + 
  north(scale=0.15, location= "bottomright", symbol = 12, x.min=306775, x.max=306950,
         y.min=4994365, y.max=4994500) +
  scalebar(location = "bottomright", dist=50, dist_unit = "m", st.size=2,
           x.min=306850, x.max=306900,
            y.min=4994375, y.max=4994525, transform = FALSE) + 
  xlim(c(306625.3,e[2])) + ylim(e[3:4]) +
  coord_quickmap()
rose_fig
ggsave("Figure3.jpg",width=8,height=6,units=c("in"),dpi=600)


plot(class.fig)
plot(rose_fig,add=TRUE)

####----alt sf way for final fig-----######## 
#back to spdf 
coordinates(r_df) <- ~x + y
#proj4string(rosie_pts) <- "+proj=longlat +datum=WGS84"
proj4string(r_df) <- "+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs"
#rosie_pts <- spTransform(rosie_pts, CRS("+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs"))

st_geometry(r_df)
ggplot() + 
  geom_sf(data=st_geometry(r_df),aes(fill=ClassPred),geometry=geometry) +
  #scale_fill_viridis(discrete = TRUE,direction=-1) +
  coord_sf()
  
