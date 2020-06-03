#Prediction Kebara 2 from N=33 H. sapiens females: MODEL C
require(Morpho)
require(geomorph)
require(Rvcg)
require(ape)
require(rgl)
library(doParallel)
source("pls_function.R")

#Load landmark data----
Trunk_data <- read.morphologika("D:/R progress/Trunk_analyses_R/Kebara2_prediction/R1/Kebara_2_prediction/N60_coordinates.txt")

Pelves_Kebara2 <- read.morphologika("D:/R progress/Trunk_analyses_R/Kebara2_prediction/Kebara2_pelves_reconstructions.txt")
Kebara2_thorax_reconstructions <- read.morphologika("D:/R progress/Trunk_analyses_R/Kebara2_prediction/R1/Kebara_2_prediction/Kebara2_thorax_reconstructions.txt")

Trunk_females <- bindArr(Trunk_data[,,1:4],Trunk_data[,,11],Trunk_data[,,14:15],Trunk_data[,,18:24],Trunk_data[,,37:54],Trunk_data[,,58],along=3)

#Resliding trunk data ----
#Define landmarks and semilandmarks----
fix <- c(1:134)
Curve1 <- c(1,135:147,46)
Curve2 <- c(23,148:160,68)
Curve3 <- c(3,161:173,47)
Curve4 <- c(25,174:186,69)
Curve5 <- c(5,187:199,49)
Curve6 <- c(27,200:212,71)
Curve7 <- c(7,213:225,51)
Curve8 <- c(29,226:238,73)
Curve9 <- c(9,239:251,53)
Curve10 <- c(31,252:264,75)
Curve11 <- c(11,265:277,55)
Curve12 <- c(33,278:290,77)
Curve13 <- c(13,291:303,57)
Curve14 <- c(35,304:316,79)
Curve15 <- c(15,317:329,59)
Curve16 <- c(37,330:342,81)
Curve17 <- c(17,343:355,61)
Curve18 <- c(39,356:368,83)
Curve19 <- c(19,369:381,63)
Curve20 <- c(41,382:394,85)
Curve21 <- c(21,395:407,65)
Curve22 <- c(43,408:420,87)
Curve23 <- c(2,421:433,45)
Curve24 <- c(24,434:446,67)
Curve25 <- c(4,447:459,48)
Curve26 <- c(26,460:472,70)
Curve27 <- c(6,473:485,50)
Curve28 <- c(28,486:498,72)
Curve29 <- c(8,499:511,52)
Curve30 <- c(30,512:524,74)
Curve31 <- c(10,525:537,54)
Curve32 <- c(32,538:550,76)
Curve33 <- c(12,551:563,56)
Curve34 <- c(34,564:576,78)
Curve35 <- c(14,577:589,58)
Curve36 <- c(36,590:602,80)
Curve37 <- c(16,603:615,60)
Curve38 <- c(38,616:628,82)
Curve39 <- c(18,629:641,62)
Curve40 <- c(40,642:654,84)
Curve41 <- c(20,655:667,64)
Curve42 <- c(42,668:680,86)
Curve43 <- c(22,681:693,66)
Curve44 <- c(44,694:706,88)
Curve45 <- c(89,707:722,128)
Curve46 <- c(106,723:737,108)
Curve47 <- c(107,738:752,109)
Curve48 <- c(110,753:760,112)
Curve49 <- c(111,761:768,113)
Curve50 <- c(114,769:776,118)
Curve51 <- c(115,777:784,119)
Curve52 <- c(108,785:790,124)
Curve53 <- c(109,791:796,125)
Curve54 <- c(116,797:802,122)
Curve55 <- c(117,803:808,123)
Curve56 <- c(114,809:826,120)
Curve57 <- c(115,827:844,121)
Curve58 <- c(126,845:849,118)
Curve59 <- c(127,850:854,119)
Curve60 <- c(118,855:862,133)
Curve61 <- c(119,863:870,134)

Curves <- list(Curve1,Curve2,Curve3,Curve4,Curve5,Curve6,Curve7,Curve8,Curve9,Curve10,Curve11,Curve12,Curve13,Curve14,Curve15,Curve16,Curve17,Curve18,Curve19,Curve20,Curve21,Curve22,Curve23,Curve24,Curve25,Curve26,Curve27,Curve28,Curve29,Curve30,Curve31,Curve32,Curve33,Curve34,Curve35,Curve36,Curve37,Curve38,Curve39,Curve40,Curve41,Curve42,Curve43,Curve44,Curve45,Curve46,Curve47,Curve48,Curve49,Curve50,Curve51,Curve52,Curve53,Curve54,Curve55,Curve56,Curve57,Curve58,Curve59,Curve60,Curve61)

Surfaces <- c(1:1030)[-c(fix,Curve1,Curve2,Curve3,Curve4,Curve5,Curve6,Curve7,Curve8,Curve9,Curve10,Curve11,Curve12,Curve13,Curve14,Curve15,Curve16,Curve17,Curve18,Curve19,Curve20,Curve21,Curve22,Curve23,Curve24,Curve25,Curve26,Curve27,Curve28,Curve29,Curve30,Curve31,Curve32,Curve33,Curve34,Curve35,Curve36,Curve37,Curve38,Curve39,Curve40,Curve41,Curve42,Curve43,Curve44,Curve45,Curve46,Curve47,Curve48,Curve49,Curve50,Curve51,Curve52,Curve53,Curve54,Curve55,Curve56,Curve57,Curve58,Curve59,Curve60,Curve61)]

proc_trunk <- procSym(Trunk_females)

relaxarray <- Trunk_females

for (i in 1:33) {
  print(i)
  
  relaxarray[,,i] <- relaxLM(Trunk_females[,,i],proc_trunk$mshape,SMvector=fix, deselect=TRUE, surp=Surfaces,outlines=Curves)
  
}

Todo_trunk <- relaxarray[,,1:33]

# Thorax and pelvis n=27 males H. sapiens
Thorax_temp <- read.morphologika("D:/R progress/Trunk_analyses_R/Kebara2_prediction/Thorax_temp_mlk_stern.txt")
Thorax_indices <- vcgKDtree(Trunk_data[,,26],Thorax_temp[,,1],k=1)$index
Thorax_females <- Todo_trunk[Thorax_indices,,]
Pelvis_temp <- read.morphologika("D:/R progress/Trunk_analyses_R/Kebara2_prediction/Pelvis_temp_mlk.txt")
Pelvis_indices <- vcgKDtree(Trunk_data[,,26],Pelvis_temp[,,1],k=1)$index
Pelvis_females <- Todo_trunk[Pelvis_indices,,]

Thorax_mesh <- vcgImport("D:/R progress/Trunk_analyses_R/Kebara2_prediction/Thorax_temp.ply", updateNormals = TRUE, readcolor = FALSE, clean = TRUE,silent = FALSE)
Pelvis_mesh <- vcgImport("D:/R progress/Trunk_analyses_R/Kebara2_prediction/Pelvis_temp.ply", updateNormals = TRUE, readcolor = FALSE, clean = TRUE,silent = FALSE)

#GPA N=33 females H. sapiens

Thorax_females_GPA <- procSym(Thorax_females)
Pelvis_females_GPA <- procSym(Pelvis_females)

#MODEL C: PLS n=33 males Homo sapiens ----
PLS_ModelC <- pls2B(Thorax_females_GPA$rotated,Pelvis_females_GPA$rotated, rounds=999,cv=T)
PLS_ModelC

col_polyg <- c("black")
col_plot <- c(rep("skyblue",33))
col_leg <- c("skyblue")

PLS_ModelC_comm <- getPLSCommonShape(PLS_ModelC)

plot(PLS_ModelC_comm$XscoresScaled[,1],PLS_ModelC_comm$YscoresScaled[,1], xlab="Block 1 (Thorax)",ylab="Block 2 (Pelvis)",main="TWO BLOCKS PLS ANALYSIS (1st latent variable)",pch=18,col=col_plot,cex=2,xlim=c(-0.15,0.15),ylim=c(-0.20,0.20))
abline(h=0.00,v=0.00,col="black")
legend("bottomright", inset=.005, legend = c("Females"), horiz=FALSE, cex=1.2,bty="n",pt.cex=1.7,col=col_leg,pch=18)

groups33thorax <- read.table("D:/R progress/Trunk_analyses_R/Kebara2_prediction/R1/Kebara_2_prediction/Groups_N33females.txt",header=TRUE)
if (require(car)) {
  for(ii in 1:length(levels(groups33thorax$SPECIES))){
    dataEllipse(PLS_ModelC_comm$XscoresScaled[groups33thorax$SPECIES==levels(groups33thorax$SPECIES)[ii],1],PLS_ModelC_comm$YscoresScaled[groups33thorax$SPECIES==levels(groups33thorax$SPECIES)[ii],1],add=TRUE,center.pch = FALSE,levels=c(.95), plot.points = F,col=col_polyg[ii],fill=F, fill.alpha=0.10,robust=TRUE)}
}


Covar_ModelC <- plsCoVar(PLS_ModelC,i=1,sdx=1,sdy=1)

PLS1_ModelC_neg_x<- tps3d(Thorax_mesh,Thorax_temp[,,1],Covar_ModelC$x[,,1])
shade3d(PLS1_ModelC_neg_x, col="skyblue",back="lines",box=FALSE,axes=FALSE,specular="white")

PLS1_ModelC_pos_x<- tps3d(Thorax_mesh,Thorax_temp[,,1],Covar_ModelC$x[,,2])
shade3d(PLS1_ModelC_pos_x, col="skyblue",back="lines",box=FALSE,axes=FALSE,specular="white")

PLS1_ModelC_neg_y<- tps3d(Pelvis_mesh,Pelvis_temp[,,1],Covar_ModelC$y[,,1])
shade3d(PLS1_ModelC_neg_y,col="skyblue",back="lines",box=FALSE,axes=FALSE,specular="white")

PLS1_ModelC_pos_y<- tps3d(Pelvis_mesh,Pelvis_temp[,,1],Covar_ModelC$y[,,2])
shade3d(PLS1_ModelC_pos_y,col="skyblue",back="lines",box=FALSE,axes=FALSE,specular="white")

#Error de predicción

truth_thorax <- vecx(Thorax_females_GPA$rotated) #Convert the 3D array into a matrix (Thorax)

errorC <- NULL
for (i in 1:dim(PLS_ModelC$predicted.x)[3]) # loop for cross-validation
  errorC <- c(errorC,mean((truth_thorax-PLS_ModelC$predicted.x[,,i])^2))

NoCVerror <- NULL
for (i in 1:dim(PLS_ModelC$predicted.x)[3]) # loop for MSE non-validated
  NoCVerror <- c(NoCVerror,mean((truth_thorax-vecx(predictPLSfromData(PLS_ModelC,y=Pelvis_females_GPA$rotated,ncomp=i)))^2))

# Cross-validation plot----

plot(1:31,errorC,main="Cross-validated prediction error",pch=19,xlab="Latent Variables used for Prediction",ylab="MSE",ylim=c(0,0.000010),col="skyblue")
lines(1:31,errorC,col="skyblue")
points(1:31,NoCVerror,col="black",pch=15)                  
lines(1:31,NoCVerror,col="black")
legend(x="topright",y=NULL,legend=c("Cross-validated", "Nonvalidated"), col = c("skyblue","black"),pch=c(19,15),lty=1)


#VALIDACIÓN DE LA PREDICCIÓN EN LA MUESTRA n=33 females

cv <- 1:dim(Trunk_females)[3]
predC <- Trunk_females[Thorax_indices,,]*0
for (i in cv) {
  plsC <- pls2B(Thorax_females_GPA$rotated[,,-i],Pelvis_females_GPA$rotated[,,-i],rounds = 99)
  predC[,,i] <- predictPLSfromData(plsC,y=Pelvis_females_GPA$rotated[,,i],ncomp=5)
}

predC[,,1:33]

actual_pred_thorax_N33 <- bindArr(Thorax_females[,,1:33],predC[,,1:33],along=3)

r2morphologika(actual_pred_thorax_N33,"actual_pred_thorax_N33.txt")

#Best prediction Model C

Best_pred <- tps3d(Thorax_mesh,Thorax_temp[,,1],predC[,,27])
shade3d(Best_pred,col="azure",back="lines",box=FALSE,axes=FALSE,specular="white")
Orig_rot <- rotonto(Thorax_females_GPA$mshape,Thorax_females[,,27],scale = T)$yrot
Orig <- tps3d(Thorax_mesh,Thorax_temp[,,1],Orig_rot)
shade3d(Orig,col="darkolivegreen1",back="lines",box=FALSE,axes=FALSE,specular="white")

#Worst prediction Model C

Worst_pred <- tps3d(Thorax_mesh,Thorax_temp[,,1],predC[,,4])
shade3d(Worst_pred,col="azure",back="lines",box=FALSE,axes=FALSE,specular="white")
Orig_rot <- rotonto(Thorax_females_GPA$mshape,Thorax_females[,,4],scale = T)$yrot
Orig <- tps3d(Thorax_mesh,Thorax_temp[,,1],Orig_rot)
shade3d(Orig,col="mistyrose2",back="lines",box=FALSE,axes=FALSE,specular="white")

#PREDICTIONS based on MODEL C----

#Prediction Kebara 2 thorax from Rak & Arensburg pelvis: Model C
Kebara2RA_toGPA_ModelC <- rotonto(Pelvis_females_GPA$mshape,Pelves_Kebara2[,,1],scale = T)$yrot
Pred_Kebara2RA_ModelC <- predictPLSfromData(PLS_ModelC,y=Kebara2RA_toGPA_ModelC,ncomp=5) #number of LVs to perform the prediction
deformGrid3d(Thorax_females_GPA$mshape,Pred_Kebara2RA_ModelC)
deformGrid3d(Pred_Kebara2RA_ModelC,Pred_Kebara2RA_ModelC)
Thorax_Kebara2RA_ModelC <- tps3d(Thorax_mesh,Thorax_temp[,,1],Pred_Kebara2RA_ModelC)
shade3d(Thorax_Kebara2RA_ModelC,col="azure",back="lines",box=FALSE,axes=FALSE,specular="white")

#Superimposition of prediction Kebara 2 thorax from Rak & Arensburg pelvis (Model C) with thorax reconstruction of G-O et al. 2018
shade3d(Thorax_Kebara2RA_ModelC,col="azure",back="lines",box=FALSE,axes=FALSE,specular="white")
Thorax_Kebara2GO_mesh_ModelC <- vcgImport("D:/R progress/Trunk_analyses_R/Kebara2_prediction/Kebara2_thorax_GO.ply", updateNormals = TRUE, readcolor = FALSE, clean = TRUE,silent = FALSE)
Thorax_Kebara2GO_toGPA_ModelC <- rotonto(Thorax_females_GPA$mshape,Kebara2_thorax_reconstructions[,,1],scale = T)$yrot
Thorax_Kebara2GO_ModelC <- tps3d(Thorax_Kebara2GO_mesh_ModelC,Kebara2_thorax_reconstructions[,,1],Thorax_Kebara2GO_toGPA_ModelC)
shade3d(Thorax_Kebara2GO_ModelC,col="deepskyblue",back="lines",box=FALSE,axes=FALSE,specular="white")

#Superimposition of prediction Kebara 2 thorax from Rak & Arensburg pelvis (Model C) with thorax reconstruction Sawyer and Maley (2005)
shade3d(Thorax_Kebara2RA_ModelC,col="azure",back="lines",box=FALSE,axes=FALSE,specular="white")
Thorax_Kebara2S_M_mesh_ModelC <- vcgImport("D:/R progress/Trunk_analyses_R/Kebara2_prediction/R1/Kebara_2_prediction/Kebara2_thorax_SM.ply", updateNormals = TRUE, readcolor = FALSE, clean = TRUE,silent = FALSE)
Thorax_Kebara2S_M_toGPA_ModelC <- rotonto(Thorax_females_GPA$mshape,Kebara2_thorax_reconstructions[,,2],scale = T)$yrot
Thorax_Kebara2S_M_ModelC <- tps3d(Thorax_Kebara2S_M_mesh_ModelC,Kebara2_thorax_reconstructions[,,2],Thorax_Kebara2S_M_toGPA_ModelC)
shade3d(Thorax_Kebara2S_M_ModelC,col="palegreen",back="lines",box=FALSE,axes=FALSE,specular="white")

#Prediction Kebara 2 thorax from Sawyer and Maley (2005) pelvis: Model C
Kebara2SM_toGPA_ModelC <- rotonto(Pelvis_females_GPA$mshape,Pelves_Kebara2[,,2],scale = T)$yrot
Pred_Kebara2SM_ModelC <- predictPLSfromData(PLS_ModelC,y=Kebara2SM_toGPA_ModelC,ncomp=5) #number of LVs to perform the prediction
deformGrid3d(Thorax_females_GPA$mshape,Pred_Kebara2SM_ModelC)
deformGrid3d(Pred_Kebara2SM_ModelC,Pred_Kebara2SM_ModelC)
Thorax_Kebara2SM_ModelC <- tps3d(Thorax_mesh,Thorax_temp[,,1],Pred_Kebara2SM_ModelC)
shade3d(Thorax_Kebara2SM_ModelC,col="azure",back="lines",box=FALSE,axes=FALSE,specular="white")

#Superimposition of prediction Kebara 2 thorax from Sawyer and Maley (2005) pelvis (Model C) with thorax reconstruction of G-O et al. 2018
shade3d(Thorax_Kebara2SM_ModelC,col="azure",back="lines",box=FALSE,axes=FALSE,specular="white")
Thorax_Kebara2GO_mesh_ModelC <- vcgImport("D:/R progress/Trunk_analyses_R/Kebara2_prediction/Kebara2_GO.ply", updateNormals = TRUE, readcolor = FALSE, clean = TRUE,silent = FALSE)
Thorax_Kebara2GO_toGPA_ModelC <- rotonto(Thorax_females_GPA$mshape,Kebara2_thorax_reconstructions[,,1],scale = T)$yrot
Thorax_Kebara2GO_ModelC <- tps3d(Thorax_Kebara2GO_mesh_ModelC,Kebara2_thorax_reconstructions[,,1],Thorax_Kebara2GO_toGPA_ModelC)
shade3d(Thorax_Kebara2GO_ModelC,col="deepskyblue",back="lines",box=FALSE,axes=FALSE,specular="white")

#Superimposition of prediction Kebara 2 thorax from Sawyer and Maley (2005) pelvis (Model C) with thorax reconstruction of Sawyer and Maley (2005)
shade3d(Thorax_Kebara2SM_ModelC,col="azure",back="lines",box=FALSE,axes=FALSE,specular="white")
Thorax_Kebara2S_M_mesh_ModelC <- vcgImport("D:/R progress/Trunk_analyses_R/Kebara2_prediction/R1/Kebara_2_prediction/Kebara_SM.ply", updateNormals = TRUE, readcolor = FALSE, clean = TRUE,silent = FALSE)
Thorax_Kebara2S_M_toGPA_ModelC <- rotonto(Thorax_females_GPA$mshape,Kebara2_thorax_reconstructions[,,2],scale = T)$yrot
Thorax_Kebara2S_M_ModelC <- tps3d(Thorax_Kebara2S_M_mesh_ModelC,Kebara2_thorax_reconstructions[,,2],Thorax_Kebara2S_M_toGPA_ModelC)
shade3d(Thorax_Kebara2S_M_ModelC,col="palegreen",back="lines",box=FALSE,axes=FALSE,specular="white")
