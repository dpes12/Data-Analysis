dat <- read.csv("Malware.csv",na.strings="")

set.seed(1)
# Creating 350 Observations
selected.rows <- sample(1:nrow(dat),size=350,replace=FALSE)

mydata <- dat[selected.rows, 2:11]

dim(mydata)#Checking the dimensions
library(dummies)
dummy <- fastDummies::dummy_cols(mydata)#Creating Dummy Variables
knitr::kable(dummy)


dummy <- fastDummies::dummy_cols(mydata, remove_first_dummy = TRUE)#Deleting one dummy column
knitr::kable(dummy)

write.csv(dummy,"mydata.csv")#Writing new data into "mydata.csv"
#Installing all the required Packages
install.packages("scatterplot3d")
library(scatterplot3d)

install.packages("scales")
library(scales)

install.packages("ggpubr")
library(ggpubr)

install.packages("FactoMineR")
library(FactoMineR)

install.packages("factoextra")
library(factoextra)

install.packages("remotes")
remotes::install_github("vqv/ggbiplot")
library(ggbiplot)

data(mydata)
data(dummy)
str(mydata)
str(dummy)

pca.dummy <- prcomp(dummy[c(1 ,7 ,9 , 11:16)], #Performing PCA analysis with numeric variables excluding Verified.as.Malware
                   scale=TRUE)#Standardised the data
summary(pca.dummy)
pca.dummy$rotation

plot(pca.dummy, type="l", main="Scree plot - dummy") #Creating scree plot


fviz_pca_biplot(pca.dummy,
                axes = c(1,2), #Specifying the PCs to be plotted.
                #Parameters for samples
                col.ind=dummy$Verified.as.Malware, #Outline colour of the shape
                fill.ind=dummy$Verified.as.Malware, #fill colour of the shape
                alpha=0.5, #transparency of the fill colour
                pointsize=4, #Size of the shape
                pointshape=21, #Type of Shape
                #Parameter for variables
                col.var="red", #Colour of the variable labels
                label="var", #Show the labels for the variables only
                repel=TRUE, #Avoid label overplotting
                addEllipses=TRUE, #Add ellipses to the plot
                legend.title=list(colour="Verified As Malware",fill="Verified as Malware",alpha="Verified As Malware"))

#cross tabulating inc.executable with verified as malware

prop.table(xtabs ( ~ inc.excutable + Verified.as.Malware , data = dummy), margin = 2)*100

#cross tabulating outside.network with verified as malware

prop.table(xtabs ( ~ outside.network + Verified.as.Malware , data = dummy), margin = 2)*100

##cross tabulating URL count with verified as malware
prop.table(xtabs ( ~ URL.count + Verified.as.Malware , data = dummy), margin = 2)*100


