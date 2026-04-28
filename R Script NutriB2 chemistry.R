library(multcomp)
library(vegan)
library(pairwiseAdonis)
library(ggplot2)
library(ape)
library(thematic)
library(RColorBrewer)
library(vegan3d)

data <- read.table(file.choose(), header=T, sep=";") 
data1 <- read.table(file.choose(), header=T, sep=";")
data1 <- data.frame(data1[,-1])

dataOsmia <- subset(data, Genus=="Andrena")
rownames(dataOsmia) <- dataOsmia[,1]
data1 <- data.frame(dataOsmia[,-1])
data1 <- data.frame(data1[,-1])
data1 <- data.frame(data1[,-1])
data1 <- subset(data1, select=-Sum)

y <- dataOsmia$Species 
table(y)

dataOsmia <- subset(dataOsmia, dataOsmia$Species != "Bombus_humilis" & dataOsmia$Species != "Bombus_lucorum" & dataOsmia$Species != "Bombus_soroeensis")  #to delete species with low n, do not forget to change data1 accordingly

x <- data$CV_FA
boxplot(x ~ data$avg_female_size , adj=1, col ="grey", xlab = "Average female size", ylab = "CV for fatty acids")


ggplot(data = data, aes(x = Genus, y = PLratio)) +
  geom_boxplot(aes(group = Genus), width = 0.5, size = 0.4, position = position_dodge(0.8)) +
  geom_point(aes(shape= Genus, fill= Species), position = position_dodge(0.8)) +
  scale_shape_manual(values = c(21, 20, 22, 23, 24, 25))+
  theme_classic()

model1 <- glm(x ~ proboiscis_length, na.action = na.omit, data = data)

anova (model1, test="F")
library(emmeans)
em <- emmeans(model1, "Genus")
contrast(em, "pairwise", adjust = "Tukey")

### P:L-ratio plot
lines <- data.frame(
  intercept = rep(0, 7),
  slope = c(0.1, 0.25, 0.5, 1, 2, 4, 0.75))


ggplot(data, aes(x=AA, y=FA, color=Speciescode, shape=Genus, fill=Speciescode)) + geom_point(size=5)+
  geom_abline(aes(intercept = intercept, slope = slope),
              linetype = "dashed", data = lines)+
  scale_shape_manual(values = c(21, 20, 22, 23, 24, 25,0,1,2,3,4,5,6))+
  theme_classic() 

# Change the low and high colors
# Sequential color scheme
sp2+scale_color_gradient(low="blue", high="red")
# Diverging color scheme
mid<-mean(data$PLratio)
sp2+scale_color_gradient2(midpoint=mid, low="blue", mid="white",
                          high="red", space ="Lab" )

### NMDS
mds2 <- metaMDS(data1, distance="bray", k=2, trymax=5, autotransform=TRUE, noshare=0.1, expand=TRUE, trace=1)
mds2
mds2$stress 
clr <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499")
pts <- c(15, 16, 17, 18, 3, 4, 8, 7, 9)
cdx <- data$Genuscode
ordiplot(mds2, display="sites", type="none") |> points("sites", pch=pts[cdx], col=clr[cdx], cex=2)
legend(x="bottomright", legend=levels(as.factor(data$Genus)), col=clr, pch=pts, cex=1.3)
ordiellipse(mds2, groups = data$Species, draw = "polygon", lty = 1, col = clr)

mds3D <- metaMDS(data1, distance="bray", k=3, trymax=5, autotransform=TRUE, noshare=0.1, expand=TRUE, trace=1)
mds3D
mds3D$stress 
ordiplot3d(mds3D, col = clr[cdx], ax.col= "black", pch = pts[cdx])
legend(x="topright", legend=levels(as.factor(data$Species)), col=clr, pch=pts, cex=0.5)

distm1 <- vegdist(data1, method="bray")
write.table(distm1,"distm1.csv", sep=";")

data.ado <- adonis2(distm1 ~ data$Species, permutations=10000) 
data.ado

pairwise.adonis<-pairwise.adonis2(distm1 ~ data$Species, data = data)
pairwise.adonis

### Environmental fits
model<-envfit(mds2,data1,perm=100) ## fits environmental vectors or factors(i.e. sugars/AAs) onto an ordination.
## The projections of points onto vectors have maximum correlation with
## corresponding environmental variables, and the factors show the averages
## of factor levels.
model2<-envfit(mds2~Species,data,perm=100) ## fits families into ordination
model2<-envfit(mds2~plantspecies,data,perm=100)
model 
model2
plot(model) 

### Dynamic range boxes
library(tidyverse)
library(dynRB)
library(ggplot2)
library(reshape2)
library(vegan)

AA.test <- read.table(file.choose(), header=T, sep=";")
r <- dynRB_VPa(AA.test, steps = 201)

Overlaps = 
  r$result %>%
  select(V1, V2, port_mean) %>%
  filter(!(V1 == V2))

p.overlap = 
  ggplot(Overlaps,
         aes(x = V1, y = V2)) +
  geom_tile(aes(fill = port_mean)) +
  theme_classic()
p.overlap

V.individual = 
  r$result %>%
  select(V1, vol_V1_mean) %>%
  mutate(V2 = V1) %>% ## to fit into plot of overlaps
  group_by(V1, V2) %>%
  summarize(niche.volume = median(vol_V1_mean))
V.individual

p.overlap +
  geom_point(
    aes(size = niche.volume),
    data = V.individual) +
  xlab("Overlapper") +
  ylab ("Overlappee") +
  scale_fill_continuous(name = "Overlap", low="#DDCC77", high="#CC6677")

##### Checking for collinearity #####
library(corrplot) # correlation plot
library(Hmisc) # correlation plot
library(GGally)

data <- read.table(file.choose(), header=T, sep=";")
data <- data.frame(data[,-1])

cor(data, method = "pearson") 
ggpairs(data)

res <- cor(data, method = "pearson") #compute correlation matrix
#round(res, 2)
res2 <- rcorr(as.matrix(data))
# p values of the correlations
res2
write.table(as.data.frame(res2$P), "p-values correlations fatty acids.csv")
#function to join all info
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}
corremat<-flattenCorrMatrix(res2$r, res2$P)
#matrix with correlation coefficients
corrplot(res, type = "upper", order = "hclust",
         tl.col = "black", tl.cex=0.7)
