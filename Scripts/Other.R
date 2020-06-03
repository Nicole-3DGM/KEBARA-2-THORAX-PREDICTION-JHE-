#Comparar cada modelo con su muestra asociada

freqA_pred <- read.table("D:/R progress/Trunk_analyses_R/Kebara2_prediction/R1/Kebara_2_prediction/FreqA_pred.txt",header = T)
den_A <- plot(density(freqA_pred$pred),col="dodgerblue3",lwd=3,xlim=c(0,0.30),ylim=c(0,22))
freqA_actual <- read.table("D:/R progress/Trunk_analyses_R/Kebara2_prediction/R1/Kebara_2_prediction/FreqA_actual.txt",header = T)
lines(density(freqA_actual$actual), col = "black",lwd=3,lty=3)
legend(x="topright",y=NULL,legend=c("Model A", "n = 27 H. sapiens males"), lty=c(1,3), lwd=3,col = c("dodgerblue3","black"))

freqB_pred <- read.table("D:/R progress/Trunk_analyses_R/Kebara2_prediction/R1/Kebara_2_prediction/FreqB_pred.txt",header = T)
den_B <- plot(density(freqB_pred$pred),col="red",lwd=3,xlim=c(0,0.30),ylim=c(0,22))
freqB_actual <- read.table("D:/R progress/Trunk_analyses_R/Kebara2_prediction/R1/Kebara_2_prediction/FreqB_actual.txt",header = T)
lines(density(freqB_actual$actual), col = "black",lwd=3,lty=3)
legend(x="topright",y=NULL,legend=c("Model B", "n = 60 H. sapiens"), lty=c(1,3), lwd=3,col = c("red","black"))

freqC_pred <- read.table("D:/R progress/Trunk_analyses_R/Kebara2_prediction/R1/Kebara_2_prediction/FreqC_pred.txt",header = T)
den_C <- plot(density(freqC_pred$pred),col="skyblue",lwd=3,xlim=c(0,0.30),ylim=c(0,22))
freqC_actual <- read.table("D:/R progress/Trunk_analyses_R/Kebara2_prediction/R1/Kebara_2_prediction/FreqC_actual.txt",header = T)
lines(density(freqC_actual$actual), col = "black",lwd=3,lty=3)
legend(x="topright",y=NULL,legend=c("Model C", "n = 33 H. sapiens females"), lty=c(1,3), lwd=3,col = c("skyblue","grey"))

#Comparar modelos entre sí

den_A <- plot(density(freqA_pred$pred),col="dodgerblue3",lwd=3,xlim=c(0,0.30),ylim=c(0,22))
lines(density(freqB_pred$pred), col = "red",lwd=3)
lines(density(freqC_pred$pred), col = "skyblue",lwd=3)

legend(x="topright",y=NULL,legend=c("Model A", "Model B","Model C"), lwd=3,col = c("dodgerblue3","red","skyblue"))

#COMPARAR CROSS-VALIDATIONS

plot(1:58,errorB,main="Cross-validated prediction error",pch=19,xlab="Latent Variables used for Prediction",ylab="MSE",ylim=c(0,0.000010),col="red")
lines(1:58,errorB,col="red")
points(1:25,errorA,col="dodgerblue3",pch=19) 
lines(1:25,errorA,col="dodgerblue3")
points(1:31,errorC,col="skyblue",pch=19)                  
lines(1:31,errorC,col="skyblue")
legend(x="topright",y=NULL,legend=c("Model A", "Model B", "Model C"), col = c("dodgerblue3","red","skyblue"),pch=19,lty=1)
