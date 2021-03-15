

#procedure from http://web.missouri.edu/~huangf/data/mvnotes/Documents/pca_in_r_2.html



#import data, teat end shape converted to numeric (F = 0, R = 1, P = 2)
library(readr)
PCA <- read_csv("udder-and-teat-traits-for-PCA.csv")

#Double check variables recognized as numeric
PCA$FTL = as.numeric(PCA$FTL)
PCA$FTW = as.numeric(PCA$FTW)
PCA$FTS = as.numeric(PCA$FTS)
PCA$RTL = as.numeric(PCA$RTL)
PCA$RTW = as.numeric(PCA$RTW)
PCA$RTS = as.numeric(PCA$RTS)
PCA$UD = as.numeric(PCA$UD)
PCA$UH = as.numeric(PCA$UH)
PCA$UW = as.numeric(PCA$UW)
PCA$FTP = as.numeric(PCA$FTP)
PCA$RTP = as.numeric(PCA$RTP)
PCA$UC = as.numeric(PCA$UC)
PCA$FUA = as.numeric(PCA$FUA)


###############################################################################################
###############################################################################################
#PCA on all udder and teat traits together
#
#


#important parameters
R <- cor(PCA[,4:16], use = "complete.obs")     #saves correlation matrix
p <- ncol(PCA[,4:16])    #number of observations

#Bartlett's test of sphericity
bart<-function(dat){
  R<-cor(dat, use = "complete.obs") #saves corr matrix
  p<-ncol(dat) #number of variables
  n<-nrow(dat) #number of observations
  chi2<- -((n-1)-((2*p)+5)/6) * log(det(R)) #this is the formula
  df <-(p*(p-1)/2)
  crit<-qchisq(.95,df) #critical value
  p<-pchisq(chi2,df,lower.tail = F) #pvalue
  cat("Bartlett's test of sphericity: X2(",df,")=",chi2,", p=",round(p,5),sep="" )
}

bart(PCA[,4:16]) #p<.001; PCA is appropriate for these data

e <- eigen(R) #Sovling for the eigenvalues and eigenvectors 
str(e)

L <- e$values #placing eigenvalues in L
Vm <- matrix(0,nrow=p,ncol=p) #creating pxp matrix with 0s
Vm

#Vm is orthogonal matrix since all correlations between variables are 0
diag(Vm) <- L #Putting eigenvalues in the diagonals 
e$vectors # These are the eigenvectors -- are the standardized regression weights 

loadings <- e$vectors %*% sqrt(Vm) #these are the loadings 
#or the correlation of the component variables with the original variables -- 
   #sometimes referred to as the P matrix. and PP` is the original correlation matrix 
#SPSS refers to this as the component matrix 

#original correlation matrix e$vectors %*% Vm %*% t(e$vectors) # V L V` 

L/length(L) #proportion of variance accouted for by each PC 

zdat<-scale(PCA[,4:16])
pca.scores<-zdat%*%-e$vectors #negative changes the directions so a larger PC score corresponds with large size and with the corrected loading weights in loadings2
colnames(pca.scores)<-c('PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9',
                        'PC10','PC11','PC12','PC13')
head(pca.scores)
#NOTE: these scores are scaled such that they have a mean of zero and the variance comes out to eigenvalue for that component
round(colMeans(pca.scores),2)
apply(pca.scores,2,var)
e$values

head(scale(pca.scores)[,1])

pv<-e$values/p #percentage of variance accounted for by each variable

round(cor(pca.scores, use = "complete.obs"),2) #0s Not an error, this is a correlation matrix, and the PC scores are not correlated with each other, completely orthogonal 

cor(PCA[,4:16],pca.scores[,1], use = "complete.obs") # this is correlating the variables with the first pc

#flipping the loadings sign

sign<-vector(mode="numeric",length=p)
sign <- sign(colSums(loadings))
loadings2<-loadings %*% diag(sign)
head(loadings2)

#Merging the pc scores with the measurement data 
Udder_Teat_PCs<-cbind(PCA,pca.scores[,1:4])

write.csv(Udder_Teat_PCs,file = "Udder_Teat_PCs.csv")    #use this later for GWA
write.csv(loadings2,file = "Udder_Teat_loadings2.csv")  #need this for loadings plot 

#Plot PC variance explained
plot(pv, pch = 16, ylab="Percent variance explained", xlab="PCs")
abline(h=0.0769)

#Plot loadings per PC
#just because it's easier, open Udder_Teat_loadings2.csv in excel and change headers to PCs and row labels to corresponding trait names, then save and import

loadingsU_T <- read_csv("udder&teat PCA/loadingsU&T.csv")

library(ggplot2)
ggplot(data = loadingsU_T, aes(x=reorder(Trait, -PC1), y=PC1)) + 
  geom_bar(stat="identity", fill="maroon4") + 
  coord_flip() + theme_bw()

ggplot(data = loadingsU_T, aes(x=reorder(Trait, -PC2), y=PC2)) + 
  geom_bar(stat="identity", fill="orchid4") + 
  coord_flip() + theme_bw()

ggplot(data = loadingsU_T, aes(x=reorder(Trait, -PC3), y=PC3)) + 
  geom_bar(stat="identity", 
  fill="darkorchid4") + 
  coord_flip() + theme_bw()

ggplot(data = loadingsU_T, aes(x=reorder(Trait, -PC4), y=PC4)) + 
  geom_bar(stat="identity", 
  fill="royalblue4") + 
  coord_flip() + theme_bw()

ggplot(data = loadingsU_T, aes(x=reorder(Trait, -PC5), y=PC5)) + 
  geom_bar(stat="identity", 
  fill="steelblue4") + 
  coord_flip() + theme_bw()


###############################################################################################
###############################################################################################
#repeat, but only on udder traits
#
#


#important parameters
R <- cor(PCA[,12:16], use = "complete.obs")     #saves correlation matrix
p <- ncol(PCA[,12:16])    #number of observations

#Bartlett's test of sphericity
bart<-function(dat){
  R<-cor(dat, use = "complete.obs") #saves corr matrix
  p<-ncol(dat) #number of variables
  n<-nrow(dat) #number of observations
  chi2<- -((n-1)-((2*p)+5)/6) * log(det(R)) #this is the formula
  df <-(p*(p-1)/2)
  crit<-qchisq(.95,df) #critical value
  p<-pchisq(chi2,df,lower.tail = F) #pvalue
  cat("Bartlett's test of sphericity: X2(",df,")=",chi2,", p=",round(p,5),sep="" )
}

bart(PCA[,12:16]) #p<.001, yay we can do PCA 


e<-eigen(R) #Sovling for the eigenvalues and eigenvectors 
str(e)

L<-e$values #placing eigenvalues in L
Vm<-matrix(0,nrow=p,ncol=p) #creating pxp matri with 0s
Vm
#Vm is orthogonal matrix since all correlations between variables are 0
diag(Vm)<-L #Putting eigenvalues in the diagonals 
e$vectors # These are the eigenvectors -- are the standardized regression weights 

loadings<-e$vectors %*% sqrt(Vm) #these are the loadings 
#or the corerlation of the component variables with the original variables -- sometimes referred to as the P matris. and PP` is the original correlation matrix 
#SPSS referes to this as the component matrix 

#original correlation matrix e$vectors %*% Vm %*% t(e$vectors) # V L V` 

L/length(L) #proportion of variance accouted for by each PC 

zdat<-scale(PCA[,12:16])
pca.scores<-zdat%*%-e$vectors #negative changes the directions so a larger PC score corresponds with large size and with the corrected loading weights in loadings 2
colnames(pca.scores)<-c('PC1','PC2','PC3','PC4','PC5')
head(pca.scores)
#NOTE: these scores are scaled such that they have a mean of zero and the variance comes out to eigenvalue for that component
round(colMeans(pca.scores),2)
apply(pca.scores,2,var)
e$values

head(scale(pca.scores)[,1])

pv<-e$values/p #percentage of variance accounted for by each variable

round(cor(pca.scores, use = "complete.obs"),2) #0s Not an error, this is a correlation matrix, and the PC scores are not correlated with each other, completely orthogonal 

cor(PCA[,12:16],pca.scores[,1], use = "complete.obs") # this is correlating the variables with the first pc

#flipping the loadings sign

sign<-vector(mode="numeric",length=p)
sign <- sign(colSums(loadings))
loadings2<-loadings %*% diag(sign)
head(loadings2)



#Merging the pc scores with the measurement data 
Udder_PCs<-cbind(PCA,pca.scores[,1:2])

write.csv(Udder_PCs,file = "Udder_PCs.csv")

write.csv(loadings2,file = "loadingsU_only.csv")

#just because it's easier, open Udder_Teat_loadings2.csv in excel and change headers to PCs and row labels to corresponding trait names, then save and import
loadingsU_only <- read_csv("loadingsU_only.csv")


#Plot PC variance explained
plot(pv, pch = 16, ylab="Percent variance explained", xlab="PCs")
abline(h=(1/5))

#Plot loadings per PC

library(ggplot2)
ggplot(data = loadingsU_only, aes(x=reorder(Trait, -PC1), y=PC1)) + geom_bar(stat="identity", 
  fill="darkseagreen4") + 
  coord_flip() + theme_bw()

ggplot(data = loadingsU_only, aes(x=reorder(Trait, -PC2), y=PC2)) + geom_bar(stat="identity", 
  fill="darkslateblue") + 
  coord_flip() + theme_bw()

###############################################################################################
###############################################################################################
#repeat, but only on udder traits
#
#



#important parameters
R <- cor(PCA[,4:11], use = "complete.obs")     #saves correlation matrix
p <- ncol(PCA[,4:11])    #number of observations

#Bartlett's test of sphericity
bart<-function(dat){
  R<-cor(dat, use = "complete.obs") #saves corr matrix
  p<-ncol(dat) #number of variables
  n<-nrow(dat) #number of observations
  chi2<- -((n-1)-((2*p)+5)/6) * log(det(R)) #this is the formula
  df <-(p*(p-1)/2)
  crit<-qchisq(.95,df) #critical value
  p<-pchisq(chi2,df,lower.tail = F) #pvalue
  cat("Bartlett's test of sphericity: X2(",df,")=",chi2,", p=",round(p,5),sep="" )
}

bart(PCA[,4:11]) #p<.001, yay we can do PCA 


e<-eigen(R) #Sovling for the eigenvalues and eigenvectors 
str(e)

L<-e$values #placing eigenvalues in L
Vm<-matrix(0,nrow=p,ncol=p) #creating pxp matri with 0s
Vm
#Vm is orthogonal matrix since all correlations between variables are 0
diag(Vm)<-L #Putting eigenvalues in the diagonals 
e$vectors # These are the eigenvectors -- are the standardized regression weights 

loadings<-e$vectors %*% sqrt(Vm) #these are the loadings 
#or the corerlation of the component variables with the original variables -- sometimes referred to as the P matris. and PP` is the original correlation matrix 
#SPSS referes to this as the component matrix 

#original correlation matrix e$vectors %*% Vm %*% t(e$vectors) # V L V` 

L/length(L) #proportion of variance accouted for by each PC 

zdat<-scale(PCA[,4:11])
pca.scores<-zdat%*%-e$vectors #negative changes the directions so a larger PC score corresponds with large size and with the corrected loading weights in loadings 2
colnames(pca.scores)<-c('PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8')
head(pca.scores)
#NOTE: these scores are scaled such that they have a mean of zero and the variance comes out to eigenvalue for that component
round(colMeans(pca.scores),2)
apply(pca.scores,2,var)
e$values

head(scale(pca.scores)[,1])

pv<-e$values/p #percentage of variance accounted for by each variable

round(cor(pca.scores, use = "complete.obs"),2) #0s Not an error, this is a correlation matrix, and the PC scores are not correlated with each other, completely orthogonal 

cor(PCA[,4:11],pca.scores[,1], use = "complete.obs") # this is correlating the variables with the first pc

#flipping the loadings sign

sign<-vector(mode="numeric",length=p)
sign <- sign(colSums(loadings))
loadings2<-loadings %*% diag(sign)
head(loadings2)


#Merging the pc scores with the measurement data 
Teat_PCs<-cbind(PCA,pca.scores[,1:3])

write.csv(Teat_PCs,file = "Teat_PCs.csv")
write.csv(loadings2,file = "loadingsT_only.csv")

#Plot PC variance explained
plot(pv, pch = 16, ylab="Percent variance explained", xlab="PCs")
abline(h=(1/8))

#Plot loadings per PC
#just because it's easier, open Udder_Teat_loadings2.csv in excel and change headers to PCs and row labels to corresponding trait names, then save and import
loadingsT_only <- read_csv("loadingsT_only.csv")

library(ggplot2)
ggplot(data = loadingsT_only, aes(x=reorder(Trait, -PC1), y=PC1)) + geom_bar(stat="identity", 
  fill="darkolivegreen3") + 
  coord_flip() + theme_bw()

ggplot(data = loadingsT_only, aes(x=reorder(Trait, -PC2), y=PC2)) + geom_bar(stat="identity", 
  fill="skyblue3") + 
  coord_flip() + theme_bw()

ggplot(data = loadingsT_only, aes(x=reorder(Trait, -PC3), y=PC3)) + geom_bar(stat="identity", 
  fill="rosybrown3") + 
  coord_flip() + theme_bw()

ggplot(data = loadingsT_only, aes(x=reorder(Trait, -PC4), y=PC4)) + geom_bar(stat="identity", 
  fill="navajowhite3") + 
  coord_flip() + theme_bw()



###############################################################################################
###############################################################################################
#repeat, but only on risk traits ID'd in Miles et al 2019 Prev Vet Med 
#
#


#import data with udder and teat traits only, teat end shape converted to numeric
PCA <- read.csv("PCA_risk-traits-only.csv")

#designate all variables as quantitative
PCA$RTW = as.numeric(PCA$RTW)
PCA$RTS = as.numeric(PCA$RTS)
PCA$UH = as.numeric(PCA$UH)
PCA$FUA = as.numeric(PCA$FUA)

#important parameters
R <- cor(PCA[,4:7], use = "complete.obs")     #saves correlation matrix
p <- ncol(PCA[,4:7])    #number of observations

#Bartlett's test of sphericity
bart<-function(dat){
  R<-cor(dat, use = "complete.obs") #saves corr matrix
  p<-ncol(dat) #number of variables
  n<-nrow(dat) #number of observations
  chi2<- -((n-1)-((2*p)+5)/6) * log(det(R)) #this is the formula
  df <-(p*(p-1)/2)
  crit<-qchisq(.95,df) #critical value
  p<-pchisq(chi2,df,lower.tail = F) #pvalue
  cat("Bartlett's test of sphericity: X2(",df,")=",chi2,", p=",round(p,5),sep="" )
}

bart(PCA[,4:7]) #p<.001, yay we can do PCA 


e<-eigen(R) #Sovling for the eigenvalues and eigenvectors 
str(e)

L<-e$values #placing eigenvalues in L
Vm<-matrix(0,nrow=p,ncol=p) #creating pxp matri with 0s
Vm
#Vm is orthogonal matrix since all correlations between variables are 0
diag(Vm)<-L #Putting eigenvalues in the diagonals 
e$vectors # These are the eigenvectors -- are the standardized regression weights 

loadings<-e$vectors %*% sqrt(Vm) #these are the loadings 
#or the corerlation of the component variables with the original variables -- sometimes referred to as the P matris. and PP` is the original correlation matrix 
#SPSS referes to this as the component matrix 

#original correlation matrix e$vectors %*% Vm %*% t(e$vectors) # V L V` 

L/length(L) #proportion of variance accouted for by each PC 

zdat<-scale(PCA[,4:7])
pca.scores<-zdat%*%-e$vectors #negative changes the directions so a larger PC score corresponds with large size and with the corrected loading weights in loadings 2
colnames(pca.scores)<-c('PC1','PC2','PC3','PC4')
head(pca.scores)
#NOTE: these scores are scaled such that they have a mean of zero and the variance comes out to eigenvalue for that component
round(colMeans(pca.scores),2)
apply(pca.scores,2,var)
e$values

head(scale(pca.scores)[,1])

pv<-e$values/p #percentage of variance accounted for by each variable

round(cor(pca.scores, use = "complete.obs"),2) #0s Not an error, this is a correlation matrix, and the PC scores are not correlated with each other, completely orthogonal 

cor(PCA[,4:7],pca.scores[,1], use = "complete.obs") # this is correlating the variables with the first pc

#flipping the loadings sign

sign<-vector(mode="numeric",length=p)
sign <- sign(colSums(loadings))
loadings2<-loadings %*% diag(sign)
head(loadings2)
#Merging the pc scores with the measurement data 
Risk.Trait.Only_PCs<-cbind(PCA,pca.scores[,1:4])

write.csv(Risk.Trait.Only_PCs, file = "Risk-Trait-Only_PCs.csv")
write.csv(loadings2,file = "loadingsRisk_only.csv")

#Plot PC variance explained
plot(pv, pch = 16, ylab="Percent variance explained", xlab="PCs")
abline(h=1/4)

#Plot loadings per PC

#just because it's easier, open Udder_Teat_loadings2.csv in excel and change headers to PCs and row labels to corresponding trait names, then save and import
loadingsRisk_only <- read_csv("loadingsRisk_only.csv")

library(ggplot2)
ggplot(data = loadingsRisk_only, aes(x=reorder(Trait, -PC1), y=PC1)) + geom_bar(stat="identity", fill="maroon4") + 
  coord_flip() + theme_bw()

ggplot(data = loadingsRisk_only, aes(x=reorder(Trait, -PC2), y=PC2)) + geom_bar(stat="identity", fill="orchid4") + 
  coord_flip() + theme_bw()

ggplot(data = loadingsRisk_only, aes(x=reorder(Trait, -PC3), y=PC3)) + geom_bar(stat="identity", 
                                                                        fill="darkorchid4") + 
  coord_flip() + theme_bw()

ggplot(data = loadingsRisk_only, aes(x=reorder(Trait, -PC4), y=PC4)) + geom_bar(stat="identity", 
                                                                        fill="royalblue4") + 
  coord_flip() + theme_bw()
















