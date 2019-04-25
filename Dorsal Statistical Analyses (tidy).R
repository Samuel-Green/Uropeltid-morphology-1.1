### Hi

# This script assumes that data will already have been generated from other scripts and ready to analyse

# Set wprking directory first of course

setwd("F:/Elliot Samuel Shayle/Box Sync/Elliot Shayle/From Filipa's Dropbox/Elliot GM project/R_files_GM")

setwd("C:/Users/ellio/OneDrive - University College London/From Filipa's Dropbox/Elliot GM project/R_files_GM") ### On my Surface w/ OneDrive)

setwd("C:/Users/ellis2/OneDrive - University College London/From Filipa's Dropbox/Elliot GM project/R_files_GM") ### NHM computer W/ OneDrive

load("Data_Needed_for_Dorsal_Analysis_Ver1.0.Rdata")

# In line with Page and Cooper (2017), I should try doing the analyses for species means and also for individual specimens

#############################################

#########################
# 0: Preamble           #
#########################

# Creating the Geomorph data frame

dorsal.snake.gpa <- gpagen(mydata.new, curves = mysliders, ProcD = F) # alignment with bending energy
snake.gdf <- geomorph.data.frame(dorsal.snake.gpa)

#########################
# 1: Regression         #
#########################

# Does body shape change with body size?

#################### Using individual level data ####################

Snake.lm <- procD.lm(coords ~ log(Csize), iter = 9999, data = snake.gdf)
summary(Snake.lm) 

### Based on the tiny R^2 value (0.00382), I can say that only 0.382% of variation in morphology is due to the size of the snake
### and that due to the massive P value, it is also a non significant result

#################### Using mean species level data ####################

### Error: I am unsure how to create a mean size, and how to attach to the species means dataset

########################################
# 2: Kmult test of phylogenetic signal #
########################################

### Methods/code borrowed from Page 2017 et al and adapted for my data

library(ape)
library(rgl)
library(geomorph)
library(geiger)
library(maps)
library(phytools)
library(BAMMtools)
library(convevol)

### Remember to have species data ready for analysis, as well as the tree https://github.com/NaturalHistoryMuseum/river-dolphin-convergence/blob/master/analyses/01-crania-sppmean.R
### For further reading on the topic: https://rfunctions.blogspot.com/2014/02/measuring-phylogenetic-signal-in-r.html

phy.data.frame <- geomorph.data.frame(coords = dorsal.mean.PCA.data$coords,
                                      phy = tree_ult,
                                      group = dorsal.mean.PCA.data$genus)

### change 'group' in phy.data.frame for Natalie's MANOVA as well where necessary

check <- name.check(tree_ult, pca, data.names = pca$species)
print(check)

phy1 <- physignal(dorsal.mean.PCA.data$coords, tree_ult)
summary(phy1)

### There is a significant effect from phylogenetic signal

##################################################################
# 3: MANOVA between genus and coordinates W/ phylogenetic signal #
##################################################################

### Again, from Natalie's paper

manova <- procD.pgls(coords ~ group, phy = phy, data = phy.data.frame, iter =  19999) ### Only worked after I removed 'RRPP = TRUE'
summary(manova) ### Despite removing 'RRPP = TRUE' the summary says it still implemented this anyway

### Non significant differences between species means by genus (surprising if I'm interpreting it right)
### It seems that there is little difference in mean shape based upon genus when phylogenetic signal is accounted for
### Or maybe this isn't surprising, becuase they share an evolutionary history, so perhaps they just didn't have much time/need to diverge

##############################################################################
# 4: MANOVA between genus/species and shape (Filipa's script recommendation) #
##############################################################################

library(geomorph)

### I will need a GPA before running the analysis. For example, 'dorsal.snake.gpa'.

### Single-Factor MANOVA ###
#create geomorph data frame
### Using species and country as a grouping variable
snake.gdf <- geomorph.data.frame(dorsal.snake.gpa, gp = as.factor(paste(info$species, info$Country, sep=", ")))
mean.snake.anova.sp.co <- procD.lm(coords ~ gp, iter = 19999, data = snake.gdf)
summary(mean.snake.anova.sp.co)
plot(mean.snake.anova.sp.co) # diagnostic plots

### Using Genus and country as a grouping variable
snake.gdf <- geomorph.data.frame(dorsal.snake.gpa, gp = as.factor(paste(info$Genus, info$Country, sep=", ")))
mean.snake.anova.Ge.co <- procD.lm(coords ~ gp, iter = 19999, data = snake.gdf)
summary(mean.snake.anova.Ge.co)
plot(mean.snake.anova.Ge.co) # diagnostic plots

### Using Genus only as a grouping variable
snake.gdf <- geomorph.data.frame(dorsal.snake.gpa, gp = as.factor(info$Genus))
mean.snake.anova.Ge <- procD.lm(coords ~ gp, iter = 19999, data = snake.gdf)
summary(mean.snake.anova.Ge)
plot(mean.snake.anova.Ge) # diagnostic plots

### Using country only as a grouping variable
snake.gdf <- geomorph.data.frame(dorsal.snake.gpa, gp = as.factor(info$Country))
mean.snake.anova.co <- procD.lm(coords ~ gp, iter = 19999, data = snake.gdf)
summary(mean.snake.anova.co)
plot(mean.snake.anova.co) # diagnostic plots

### Using species only as a grouping variable
snake.gdf <- geomorph.data.frame(dorsal.snake.gpa, gp = as.factor(info$species))
mean.snake.anova.sp <- procD.lm(coords ~ gp, iter = 19999, data = snake.gdf)
summary(mean.snake.anova.sp)
plot(mean.snake.anova.sp) # diagnostic plots

### For test outcomes, refer to my OneNote document

### There was also an option to run a 'factorial MANOVA' with more than 1 factor
### I don't think this is relevant since I don't have that many factors to compare between
### But willing to be proved wrong

###################################################
# 5: 2B-PLS for analysing lateral and dorsal view #
###################################################

### Make sure to run both dorsal and lateral tests first

### Note that this analysis used the individual sample level data

dorsal.lateral.2BPLS <- two.b.pls(dorsal.snake.gpa$coords, lateral.snake.gpa$coords, iter = 999, print.progress = T) 

### Error ^: Unequal number of specimens in each sample
### Find a script that can drop non matching entries...

print(dorsal.lateral.2BPLS)
plot(dorsal.lateral.2BPLS)

### Note that this analysis used mean data

dorsal.lateral.means.2BPLS <- two.b.pls(dorsal.mean.PCA.data$coords, lateral.mean.PCA.data$coords, iter = 999, print.progress = T)
print(dorsal.lateral.means.2BPLS)
plot(dorsal.lateral.means.2BPLS)

### For results and analysis, please refer to my OneNote document