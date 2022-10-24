### Script for PCA
### Customized by mperdomo 2022 for Vandy LATTE Fall 2022

### you must first download each package using the function install.packages()
library(devtools)
install_github("vqv/ggbiplot")
install.pacakages(tidyverse)
install.pacakages(FactoMineR)
install.pacakages(factoextra)

### upload necessary libraries
library(ggbiplot)
library(tidyverse)
library(FactoMineR)
library(factoextra)

### this function clears the global environment
rm(list = ls())

filename = "PCA"
today.date = "10212022"

### set the working directory
### this is the folder in which your source files can be found
setwd("/Users/YOURCOMPUTERUSERNAME/Desktop/PCA Practice/texts")
TAALES <- read.csv("TAALES.csv")
Groups <- read.csv("groups.csv")

### Convert the first column to from the data frame 
### into a row name so we can run the PCA
TAALES <- TAALES %>% 
  remove_rownames %>% 
  column_to_rownames(var="Filename")

### Select the indices that will be used in the PCA 
### For an exploratory study, you may choose the included indices
### after controlling for normality. Otherwise, you make your indices
### selection based on theoretical grounds or previous work.
TAALES_clean <- select(TAALES, c("KF_Freq_AW",
                           "TL_Freq_AW",
                           "Brown_Freq_AW",
                           "MRC_Familiarity_AW",
                           "MRC_Concreteness_AW",
                           "MRC_Imageability_AW",
                           "MRC_Meaningfulness_AW",
                           "Kuperman_AoA_AW",
                           "Brysbaert_Concreteness_Combined_AW",
                           "SUBTLEXus_Freq_AW",
                           "SUBTLEXus_Range_AW",
                           "BNC_Written_Freq_AW",
                           "BNC_Written_Range_AW",
                           "BNC_Spoken_Freq_AW",
                           "BNC_Spoken_Range_AW",
                           "BNC_Spoken_Bigram_Normed_Freq",
                           "BNC_Written_Bigram_Freq_Normed",
                           "BNC_Spoken_Trigram_Normed_Freq",
                           "BNC_Written_Trigram_Freq_Normed",
                           "BNC_Spoken_Bigram_Proportion",
                           "BNC_Written_Bigram_Proportion",
                           "BNC_Spoken_Trigram_Proportion",
                           "BNC_Written_Trigram_Proportion",
                           "COCA_Academic_Range_AW",
                           "COCA_Academic_Frequency_AW",
                           "COCA_fiction_Range_AW",
                           "COCA_fiction_Frequency_AW",
                           "COCA_magazine_Range_AW",
                           "COCA_magazine_Frequency_AW",
                           "COCA_news_Range_AW",
                           "COCA_news_Frequency_AW",
                           "COCA_spoken_Range_AW",
                           "COCA_spoken_Frequency_AW",
                           "COCA_Academic_Bigram_Frequency",
                           "eat_types",
                           "eat_tokens",
                           "USF",
                           "McD_CD",
                           "Sem_D",
                           "All_AWL_Normed",
                           "All_AFL_Normed",
                           "Spoken_AFL_Normed",
                           "Written_AFL_Normed",
                           "Freq_HAL",
                           "Log_Freq_HAL",
                           "Ortho_N",
                           "Phono_N",
                           "OG_N",
                           "Freq_N",
                           "OLD",
                           "OLDF",
                           "PLD",
                           "PLDF",
                           "BG_Mean",
                           "LD_Mean_RT",
                           "WN_Mean_RT",
                           "BG_Mean_CW",
                           "content_poly",
                           "hyper_noun_S1_P1",
                           "hyper_verb_S1_P1",
                           "hyper_verb_noun_s1_p1"))



# PCA ---------------------------------------------------------------------
### https://www.statmethods.net/advstats/factor.html

### One option is using factoMineR 
### https://rdrr.io/cran/FactoMineR/f/vignettes/FactoMineR.Rmd

#################### prcomp() ########################
######################################################
### http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/118-principal-component-analysis-in-r-prcomp-vs-princomp/

#### Take a peak at PCA
TAALES.pca <- prcomp(TAALES_clean, center = TRUE, scale. = TRUE)
summary(TAALES.pca)

#### Loadings (correlations between variables and factor)
#### These are simply coefficients of principal components (recall that principal components are linear combinations of original variables) that are equal to eigenvectors entries
head(TAALES.pca$rotation)
#### by individual rows
head(TAALES.pca$x)

### Choose the number of components to retain based on the eigenvalue
ncomp <- 10 #change the number here

### Extract the factor rotation
### Rotation transforms the original factors into new ones that are easier to interpret
### https://stats.stackexchange.com/questions/59213/how-to-compute-varimax-rotated-principal-components-in-r
rawLoadings     <- TAALES.pca$rotation[,1:ncomp] %*% diag(TAALES.pca$sdev, ncomp, ncomp)
rotatedLoadings <- varimax(rawLoadings)$loadings
invLoadings     <- t(pracma::pinv(rotatedLoadings))
scores          <- scale(TAALES_clean) %*% invLoadings
print(scores[1:5,])                   # Scores computed via rotated loadings
rotatedLoadings
pcs <- as.data.frame(scores)

### Save the loadings file with varimax rotations
### Varimax rotations maximize the sum of the variance of the squared loadings
### This rotation is orthogonal, where we assume the variables are not correlated

write.csv(unclass(loadings(varimax(rawLoadings))), file = "pca.loadings.csv")

## Eigenvalues
## The eigenvalues measure the amount of variation retained by each principal component. Eigenvalues are large for the first PCs and small for the subsequent PCs.
## An eigenvalue > 1 indicates that PCs account for more variance than accounted by one of the original variables in standardized data. This is commonly used as a cutoff point for which PCs are retained. This holds true only when the data are standardized.
eig.val <- get_eigenvalue(TAALES.pca)
eig.val

#Another way of deciding how many components to retain is to use Kaiser’s criterion: that we should only retain principal components for which the variance is above 1 (when principal component analysis was applied to standardised data). We can check this by finding the variance of each of the principal components. (look at eigenva)
  
## Screeplot
fviz_screeplot(TAALES.pca, addlabels = TRUE, ncp=12)

## A simple method to extract the results, for variables, from a PCA output is to use the function get_pca_var(). This function provides a list of matrices containing all the results for the active variables (coordinates, correlation between variables and axes, squared cosine and contributions)
var <- get_pca_var(TAALES.pca)
var

# Coordinates
# Correction by square root of eigenvalue is done to standardize the variance of PC scores to 1 and therefore to allow for correlation interpretation of loadings. These standardized loadings are sometimes called loadings as well. See for example PCA function from FactoMineR package. It never uses a word loadings, it uses word coordinates for standardized loadings.
head(var$coord)
#write.csv(var$coord, "/Volumes/My Passport/Digital Text Exposure/BOX/coordinates.csv")
website_pcs <- TAALES.pca$ind$coord
#write.csv(website_pcs, "/Volumes/My Passport/Digital Text Exposure/BOX/BOX NLP 2021/websites.csv")

# Cos2: quality on the factore map
head(var$cos2)

# Contributions to the principal components: contains the contributions (in percentage) of the variables to the principal components. The contribution of a variable (var) to a given principal component is (in percentage) : (var.cos2 * 100) / (total cos2 of the component)
head(var$contrib)
#write.csv(var$contrib, "/Volumes/My Passport/Digital Text Exposure/BOX/contributions.csv")


### Look at correlations
#devtools::install_github("vsimko/corrplot")
library("corrplot")
corrplot(var$contrib, is.corr=FALSE, tl.cex = 0.4)

# Contributions of variables to PC1
fviz_contrib(TAALES.pca, choice = "var", axes = 1, top = 30)

##The most important (or, contributing) variables can be highlighted on the correlation plot as follow:
fviz_pca_var(TAALES.pca, col.var = "contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07")
)

# Change the transparency by contrib values
fviz_pca_var(TAALES.pca, alpha.var = "contrib")

#### Plot PCA
ggbiplot(TAALES.pca)
ggbiplot(TAALES.pca, labels=rownames(TAALES))

# top 5 contributing individuals and variable
fviz_pca_biplot(TAALES.pca, select.ind = list(contrib = 20), 
                select.var = list(contrib = 5),
                ggtheme = theme_minimal())

#### PC1 v. PC2
#### add/remove var.axes=FALSE for vectors
ggbiplot(TAALES.pca, obs.scale = 1, var.scale = 1, choices=c(1,2),
              groups = Groups$website, ellipse = TRUE, circle = TRUE, var.axes=FALSE) + 
  scale_color_discrete(name = '')

#### PC3 v. PC4
ggbiplot(TAALES.pca, obs.scale = 1, var.scale = 1, choices=c(3,4),
         groups = Groups$website, ellipse = TRUE, circle = TRUE, var.axes=FALSE) + 
  scale_color_discrete(name = '')

#### PC5 v. PC6
ggbiplot(TAALES.pca, obs.scale = 1, var.scale = 1, choices=c(5,6),
         groups = Groups$website, ellipse = TRUE, circle = TRUE, var.axes=FALSE) + 
  scale_color_discrete(name = '')

#### PC7 v. PC8
ggbiplot(TAALES.pca, obs.scale = 1, var.scale = 1, choices=c(7,8),
         groups = Groups$website, ellipse = TRUE, circle = TRUE, var.axes=FALSE) + 
  scale_color_discrete(name = '')

#### PC9 v. PC10
ggbiplot(TAALES.pca, obs.scale = 1, var.scale = 1, choices=c(9,10),
         groups = Groups$website, ellipse = TRUE, circle = TRUE, var.axes=FALSE) + 
  scale_color_discrete(name = '')



### More info
### https://aaronschlegel.me/factor-analysis-principal-component-method-r-part-two.html

