# PCA with species means
# march 2019

##################

# Load packages
library(devtools)
library(geomorph) # geometric morphometric analyses for 2D/3D data
library(dplyr)

##################

# see loaded objects in workspace
ls()
#clean R workspace
rm(list=ls())

# Set working directory

setwd("C:/Users/ellis2/OneDrive - University College London/Combined Photo + TPS library for David/Lateral") ### NHM computer W/ OneDrive

setwd("F:/Elliot Samuel Shayle/Box Sync/Elliot Shayle/Combined Photo + TPS library for David/Lateral") ### Main desktop W/ Box

##################

## Reading data into geomorph

# Read landmark data from tps file (exported from TPSdig)
mydata = readland.tps("LateralExpandedHeadLMDataVer1.1.tps", specID = "ID", readcurves=TRUE )
print(mydata) # prints all data in console

dim(mydata) #retrieves dimensions of an object
# x landmarks, 2 dimensions, y specimens

# Read in text file (for landmarks or other variables)
# Could have info on genera, species, habitat, food, etc.
info <- read.csv("ExpandedLateralClassifierVer1.3.csv", header=TRUE, row.names=1)
print(info)

#need to delete a few curve points that fall on top of fixed landmarks
#check which ones need to be deleted
mysliders=define.sliders(mydata, write.file = FALSE)

## Delete landmarks ##
delete = c(6, 35, 36, 65) # make a vector of the landmark numbers to delete
# Remember, 'delete' is now the name of a vector

#delete selected landmarks
mydata.new <- mydata[-delete,,] 


### Define sliding landmarks,
#https://www.rdocumentation.org/packages/geomorph/versions/3.0.7/topics/define.sliders
?define.sliders

#mysliders <- define.sliders(mydata.new, nsliders = 56, write.file = TRUE) # Manual definition of sliders

mysliders <- rbind(define.sliders(c(2,6:33,3)), define.sliders(c(1,34:61,3))) # Automatic definition of sliders (semi-landmarks)

dim(mysliders)

print(mysliders)

##################

# Running GPA - with fixed and sliding-landmarks

#In geomorph, gpagen() will align specimens with fixed landmarks, or a combination of
#fixed and sliding semilandmarks. For the latter, one must specify which landmarks are 
#semilandmarks (and for semilandmarks on curves, which landmarks they slide in-between).

?gpagen #check default parameters before running

### Procrustes analysis. There are 2 methods available

lateral.snake.gpa <- gpagen(mydata.new, curves = mysliders, ProcD = T) #GPA-alignment using Procrustes distance
## how to run GPA with bending energy, as recommended by DC Adams
lateral.snake.gpa <- gpagen(mydata.new, curves = mysliders, ProcD=FALSE) #Using bending energy

## Descriptive data
summary(lateral.snake.gpa)
str(lateral.snake.gpa)
plotAllSpecimens(lateral.snake.gpa$coords)    #GPA-aligned data

# Calculate means with multiple categories
# which labels/colours/groups would it be useful to represent in the plot?
# here, we're only creating 2 groups - species and genera, but could easily create
# more with country (India/SL), etc

# transform data into 2d array before calculating means, 
y = two.d.array(lateral.snake.gpa$coords)
proc2D <- as.data.frame(y) #transform into dataframe
proc2D$species <- info$species #add group info
proc2D$Genus <- info$Genus
proc2D$Country <- info$Country
str(proc2D)

#After grouping by species and genus, calculate means
specmeans <- as.data.frame(dplyr::summarise_each(dplyr::group_by(proc2D, species, Genus, Country), dplyr::funs(mean)))
str(specmeans)

# PCA can only use numeric data, so let's remove groups info and add in the species names as row names
final <- data.frame(specmeans[,-c(1:3)], row.names = specmeans[,1])
str(final)

# Data matrix needs to be put back in a 3D array (see 'arrayspecs')
specmeans3d <- arrayspecs(final, 61, 2)
str(specmeans3d)


# Put it all together for further analyses
lateral.mean.PCA.data <- list(coords = specmeans3d,
                     species = specmeans$species,
                     genus = specmeans$Genus)

############################

# Graphics and Visualization
### PCA with landmarks and semi-landmarks

### Basic PCA with landmarks and semi-landmarks

# PCA plot (tangent space) 
?plotTangentSpace

lateral.PCA.means <- plotTangentSpace(lateral.mean.PCA.data$coords, axis1=1, axis2 = 2,
                              label = lateral.mean.PCA.data$species,
                              groups = lateral.mean.PCA.data$genus,
                              legend = F,
                              add = T)

print(lateral.PCA.means) # has info on proportion of variance for PC1 and 2

### Stylish new PCA

#https://www.r-bloggers.com/tips-tricks-7-plotting-pca-with-tps-grids/

#group information stored in the object "specmeans"
#let's create a vector to colour the groups (genera)
col.gp <- rainbow(length(levels(specmeans$Genus))) # generates a set of different colour over the rainbow spectrum
# note! - so assign individual colours to each genus, see website above

#This vector needs to have dimension labels
names(col.gp) <- levels(specmeans$Genus)

#Using match() we can generate a vector of length(n) assigning a colour to each specimen
col.gp <- col.gp[match(specmeans$Genus, names(col.gp))]

# your PCA results are stored in "PCA.means"

lateral.PCA.means$pc.summary$importance ##PCA.means is a list containing the pc.summary
#this  can be used to make the new axis labels
xlab = paste("PC1", " (", round(lateral.PCA.means$pc.summary$importance[2,1]*100, 1), "%)", sep="")
ylab = paste("PC2", " (", round(lateral.PCA.means$pc.summary$importance[2,2]*100, 1), "%)", sep="")

#use the advanced plotting functions layout() and par() to build up this multipart figure
mat <- matrix(c(4,5,0,1,1,2,1,1,3), 3)
#divided up the plotting window into a 3 by 3 grid, and the numbers correspond to the order 
#and location of each item being plotted

layout(mat, widths=c(1,1,1), heights=c(1,1,1))# set the size of the rows and columns

# Item 1 to plot, the graph of PC1 vs PC2
par(mar=c(4, 4, 1, 1)) # sets the margins

#PCA.means also contains the pc scores in PCA.means$pc.scores.
#Here we take the first two columns for PC1 and PC2
plot(lateral.PCA.means$pc.scores[,1], lateral.PCA.means$pc.scores[,2], 
     pch=21, cex=2, bg=col.gp, xlab=xlab, ylab=ylab, asp=T, cex.lab=1.2) #cex.lab is to change font size of the axis label
legend("bottomright", legend= unique(specmeans$Genus), pch=19,  col=unique(col.gp))
text(x=lateral.PCA.means$pc.scores[,1], y=lateral.PCA.means$pc.scores[,2], labels = lateral.mean.PCA.data$species, pos = 3)

#Finally we add the TPS grids to demonstrate the shape changes along each axis. PCA contains the shape matrices PCA$pc.shapes.
meanshape <- mshape(specmeans3d)# assign mean shape for use with plotRefToTarget below
# Find the mean species
meanspec <- findMeanSpec(specmeans3d)
ref <- specmeans3d[, , meanspec]

ref <- mshape(lateral.mean.PCA.data$coords)# assign mean shape for use with plotRefToTarget below

# plot thin plate splines
par (mar = c (0, 0, 0, 0)) # sets the margins
PC1min=plotRefToTarget(ref, lateral.PCA.means$pc.shapes$PC1min)
PC1max=plotRefToTarget(ref, lateral.PCA.means$pc.shapes$PC1max)
PC2min=plotRefToTarget(ref, lateral.PCA.means$pc.shapes$PC2min)
PC2max=plotRefToTarget(ref, lateral.PCA.means$pc.shapes$PC2max)

#Done!

# Make a new dataframe
pca <- data.frame(lateral.PCA.means$pc.scores,
                  species = lateral.mean.PCA.data$species,
                  genus = lateral.mean.PCA.data$genus)

# Sort PCA data based on the tip labels of the phylogeny
pca <- pca[match(tree_ult$tip.label, pca$species), ]


### Old, dodgy phylomorphospaces

# Plot the phylomorphospaces 
col.node <- c("chartreuse", "blue4")

phylomorphospace(tree_ult, pca[, 1:2], label= "horizontal", control= col.node, add= F) # My attempt

phylomorphospace(tree_ult, X=pca[,c(1,2)], label="horizontal", lwd=1, xlab="PC 1",ylab="PC 2", axes=T) # Filipa's attempt

project.phylomorphospace(tree_ult, pca[, 1:2], nsteps = 200, sleep = 0.1,
                         direction = "both") # Animated phylomorphospace

### Stylish, new phylomorphospaces

library(phytools)
plotTree(tree_ult,node.numbers=T)

cols<-rep("grey",length(tree_ult$tip.label)+tree_ult$Nnode)
names(cols)<-1:length(cols)
# we want to plot  nodes descended from the following with different colours:
# Rhinophis - node 47
cols[getDescendants(tree_ult,47)]<-"#0092ff"
# Plactypletrurus 44
cols[getDescendants(tree_ult,44)]<-"#ffdb00"
# Plectrurus 45
cols[getDescendants(tree_ult,45)]<-"#48ff00"
# Uropeltis 41
cols[getDescendants(tree_ult,41)]<-"#ff00db"
# Teretrurus 39
cols[getDescendants(tree_ult,39)]<-"#4900ff"
# Platypectrurus 38
cols[getDescendants(tree_ult,38)]<-"#ffdb00"
# Melanophidium 33
cols[getDescendants(tree_ult,33)]<-"#ff0000"

#you can put whatever colours you'd like
#http://sape.inf.usi.ch/quick-reference/ggplot2/colour

#plot
phylomorphospace(tree_ult, pca[,c(1,2)],control=list(col.node=cols),
                 label="horizontal", lwd=1, 
                 xlab="PC 1",ylab="PC 2",
                 xlim=c(-0.2,.2))

#plot ultrametric tree and phylomorphospace side by side
layout(matrix(c(1,2),1,2))
plotTree(tree_ult, ftype = "i", mar=c(3.2,0.1,1.4,0.1))
par(mar=c(4.1,4.1,2.1,0.5))
phylomorphospace(tree_ult, pca[,c(1,2)],control=list(col.node=cols),
                 label="horizontal", lwd=1, 
                 xlab="PC 1",ylab="PC 2",
                 xlim=c(-0.2,.2))
legend("topright", legend= unique(specmeans$Genus), pch=19,  col=unique(col.gp))

# tree_ult comes from Tree_Script_FS.R file

#to change the colour you'll need to look up how to do that 


### Trying to plot PCA & Morphospace

multi.tree <- multi2di(tree_ult,random = T) # to make the tree 'fully bifurcating'

plotGMPhyloMorphoSpace(multi.tree, pca[,1:2], tip.labels = F, node.labels = F, ancStates = F, xaxis = 1, yaxis = 2, zaxis = NULL) # Error in match.names(clabs, names(xi)) : names do not match previous names

###############################################################

### Maybe Natalie Cooper's code would help with making pretty phylomorphospaces

# ----------------------------------
# STEP 4: Plot the phylomorphospaces
# ----------------------------------
# Create a colour vector (in the order of the tips, then nodes)
# River dolphin tips will be blue, others grey and nodes black
plot.colours <- c("gray56", "deepskyblue3")
cols <- c(plot.colours[pca$river], rep("black", whaletree$Nnode))
names(cols)<- 1:(length(whaletree$tip.label) + whaletree$Nnode)

# Plot the phylomorphospaces 
# Exclude taxon names to make plots clearer
phylomorphospace(whaletree, pca[, 1:2], control = list(col.node = cols), label = "off")
phylomorphospace(whaletree, pca[, c(1, 3)], control = list(col.node = cols), label ="off")
phylomorphospace(whaletree, pca[, 2:3], control = list(col.node = cols), label = "off")

##############################################################################################

#-----------------------------
# Plot thin plate splines
# ----------------------------

# Calculate the mean shape using mshape
meanshape <- mshape(specmeans3d)

# Find the mean species
meanspec <- findMeanSpec(specmeans3d)
ref <- specmeans3d[, , meanspec]

# plot thin plate splines
par (mar = c (0, 0, 0, 0))
plotRefToTarget(ref, lateral.PCA.means$pc.shapes$PC1min)
plotRefToTarget(ref, lateral.PCA.means$pc.shapes$PC1max)
plotRefToTarget(ref, lateral.PCA.means$pc.shapes$PC2min)
plotRefToTarget(ref, lateral.PCA.means$pc.shapes$PC2max)
plotRefToTarget(ref, lateral.PCA.means$pc.shapes$PC3min)
plotRefToTarget(ref, lateral.PCA.means$pc.shapes$PC3max)
plotRefToTarget(ref, lateral.PCA.means$pc.shapes$PC4min)
plotRefToTarget(ref, lateral.PCA.means$pc.shapes$PC4max)


##################
## Another way of doing it, if we only needed to add species info
# create geomorph dataframe
dorsal_df = geomorph.data.frame(lateral.snake.gpa, genus = info$Genus, ID = info$Accession_ID, species = info$species)
attributes(dorsal_df)

# calculate means per species
y <- two.d.array(lateral.snake.gpa$coords)
means <- aggregate(y ~ dorsal_df$species, FUN=mean)
means
str(means)

# pca can only use numeric data, so let's remove species info and add in the species names as row names
final <- data.frame(means[,-c(1)], row.names = means[,1])
str(final)

# Data matrix needs to be put back in a 3D array (see 'arrayspecs')
specmeans3d <- arrayspecs(final, 65, 2)
str(specmeans3d)

############################

# Graphics and Visualization

# PCA plot (tangent space) 

lateral.PCA.means <- plotTangentSpace(specmeans3d, axis1=1, axis2 = 2,
                              label = means[,1],
                              groups = proc2D$genus,
                              legend = F)

?plotTangentSpace

print(lateral.PCA.means) # has info on proportion of variance for PC1 and 2


