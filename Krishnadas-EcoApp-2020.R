#######################################################################################
## WESTERN GHATS (KODAGU) DATA ##
## BETA DIVERSITY PATTERNS BETWEEN FRAGMENT AND INTACT FOREST
## Data from: Seed size predicts community composition and carbon storage potential 
## of tree communities in rainforest fragments in India's Western Ghats. 
## Dryad Digital Repository. doi:/10.5061/dryad.7s7r1
#######################################################################################

## set your working directory
# wd = file.path("")
# setwd(wd)

# load required packages
library(nlme)
library(vegan)
library(ggplot2)
library(ggord)
library(reshape)
library(reshape2)
library(ggsignif)
library(grid)
library(gridExtra)
library(SpadeR)
library(betapart)
library(ecodist)
library(adespatial)
#for varpart plots
library(limma)
library(tidyverse)
library(tidyr)
library(dplyr)
library(ggforce)
library(hierDiversity)

## set ggplot theme to be used for all plots
thm = theme_bw()+theme(plot.background = element_blank(), text = element_text(size=14),# family = "sans serif"), 
                       panel.grid.minor = element_blank(), panel.grid.major = element_blank(),  
                       strip.background= element_blank())

## INPUT DATA
kod.dat = read.csv("Osuri_Sanakran_2016_JAE_plot_data.csv", header=T) #retained filename from Dryad
kod.dat = kod.dat[-which(is.na(kod.dat$species)),] #remove UID species
kod.loc = read.csv("Coorg_sites_UTM.csv", header = T)

soil.dat = read.csv("soildata.csv", header = T)
soil = aggregate(C ~ uniqid, data = soil.dat, FUN = mean, na.rm=F)
soil$N = aggregate(N ~ uniqid, data = soil.dat, FUN = mean)[,2]
kod.loc = merge(kod.loc, soil, by = "plot.code", by.y = "uniqid", all = T)
kod.loc = kod.loc[-c(22:25),] # gets repeated
kod.loc$CN = with(kod.loc, C/N)

kod.dat = merge(kod.dat, kod.loc, by = "plot.code", all = T)
kod.dat = kod.dat[-which(is.na(kod.dat$plot.id)),]
kod.dat$plot.all = with(kod.dat, paste(plot.code, type, sep = "-"))

# remove sites that are spatially disjunct from others
kod.dat = kod.dat[-which(kod.dat$lat > 12.3),]
kod.dat = droplevels(kod.dat[which(kod.dat$site!='Aabilu'),])
kod.dat = droplevels(kod.dat[which(kod.dat$site!='Biruga'),])
str(kod.dat)

## Table S1: species occurrences b/w control and fragment
mat.kod = as.data.frame(cast(kod.dat, species ~ type, value = 'diameter_cm', fun = length))
rownames(mat.kod) = mat.kod$species
View(mat.kod)
mat.kod$Occurrence = with(mat.kod, ifelse(CONTROL==0, "F",
                                          ifelse(FRAGMENT==0, "C", "B")))
write.csv(mat.kod, file = "TableS1-species occurences.csv")

###############################################################################################
## GAMMA DIVERSITY: CHAO-SHEN INDICES OF ENTROPY
###############################################################################################
#install_github('AnneChao/SpadeR')
## SUBSET DATA BY HABITAT
tree.c = droplevels(subset(kod.dat, type=="CONTROL"))
tree.f = droplevels(subset(kod.dat, type=="FRAGMENT"))

## COMMUNITY MATRICES (SITE X SPECIES) FOR FRAGMENT AND CONTROL
mat.f = with(tree.f, table(plot.code, species))
mat.c = with(tree.c, table(plot.code, species))

cs.f = ChaoSpecies(mat.f,"abundance", k=5, conf=0.95)
cs.c = ChaoSpecies(mat.c,"abundance", k=5, conf=0.95)

cs.all = as.data.frame(rbind(cs.f$Species_table[5,], cs.c$Species_table[5,]))
cs.all$type = c("Fragment","Contiguous")

gp.cs = ggplot(data=cs.all, aes(x=type, y=Estimate)) +
  geom_point(size=2) + ggtitle(expression(paste("a. ",gamma,"-diversity (Richness: Chao-Shen)"))) +
  geom_errorbar(data = cs.all, aes(x=type, ymin=`95%Lower`, ymax=`95%Upper`, 
                                   width=0.3, alpha = 0.4)) +
  thm + theme_classic() + theme(legend.position = "none", 
                                axis.title.x = element_blank(),
                                axis.text.x = element_text(size=12))
## SPECIES ACCUMULATION
accu.f = specaccum(mat.f, method = "exact", permutations = 100,
                   conditioned=TRUE, gamma = "chao")
accu.c = specaccum(mat.c, method = "exact", permutations = 100,
                   conditioned=TRUE, gamma = "chao")
# Fig S2
png("FigS2-specaccum.png", height=3, width=3.5, 
    pointsize=6, units="in", res=600)

plot(accu.c, xlab="Sites", ci = 2, ci.type = "bar", lty = "solid", ci.col = 'black', 
     ci.lty = 1, ylab = "Rarefied species richness")
plot(accu.f, add=T, method="rarefaction", ci = 2, ci.type = "bar", 
     lty = "dashed", ci.col = 'black', ci.lty = 1, 
     ylab = "", ylim=c(0,180), xlim=c(0,30))
leg=c("Fragment","Contiguous")
legend(x=1, y=150,legend=leg, lty=c("dashed","solid"), bty="n")
dev.off()

##############################################################################################
## ALPHA DIVERSITY: RAREFIED SPECIES RICHNESS
##############################################################################################
mat.all = cast(kod.dat, plot.all ~ species, value = 'diameter_cm', fun = length)
row.names(mat.all) = mat.all[,1]
mat.all = mat.all[,-c(1,ncol(mat.all))]

plots  =  strsplit(rownames(mat.all), split = "-")
location  =  sapply(plots, "[", 1)
type  =  sapply(plots, "[", 2)

locs = cbind(location, type)
dat.loc = merge(locs, kod.loc, by.x = 'location', by.y = 'plot.code')

rare.full = rarefy(mat.all, sample = 20, se=T)
div.full = data.frame(rare = rare.full[1,], location = location, type = type)
div.full.2 = merge(div.full, dat.loc, by = "location")
div.full.2$type = with(div.full.2, ifelse(type.x=="CONTROL","Contiguous","Fragment"))

gp.rr = ggplot(data = div.full.2, aes(x=type, y=rare)) + geom_boxplot() + 
  thm + theme(axis.text.x = element_blank(), axis.title = element_blank()) +
  geom_signif(comparisons = list(c("Contiguous", "Fragment")), map_signif_level=TRUE) +
  ggtitle(expression(paste("b. ",alpha,"-diversity (Rarefied richness)"))) + 
  theme_classic() + theme(axis.title = element_blank(), axis.text.x = element_text(size=12))

## Figure 1: ALPHA AND GAMMA 
gp.div = grid.arrange(arrangeGrob(gp.cs, gp.rr, ncol=2))
ggsave(gp.div, file="0-Fig1-alpha&gamma.png", height=4, width=7, units="in", dpi=600)

## HIERARCHICAL PARTITIONING OF DIVERSITY
groups  =  as.matrix(cbind(sort(type), rep("total", times=nrow(groups))))

## partitioning
hier.tree  =  hierDiversity(dat=as.matrix(mat.all)[order(type,location),],
                            group=groups)
names(hier.tree)
hier.div = as.data.frame(rbind(t(hier.tree$lev2$total),
                               t(hier.tree$lev1$CONTROL), 
                               t(hier.tree$lev1$FRAGMENT)))
#t(hier.tree$lev2$total)))
hier.div$index = c('a. Alpha','b. Beta','c. Gamma','c. Turnover','d. Homogeneity')
rownames(hier.div) = NULL
hier.div$type = rep(c("Landscape","Contiguous","Fragment"), each = 5)

hier.div2 = subset(hier.div, index!='c. Turnover'& 
                     index!='d. Homogeneity' & type!='Landscape')
# Figure S3
gp.hier = ggplot(data = hier.div2, aes(x=type, y=total)) +
  geom_point(size = 2) + facet_wrap(~index, scales = "free_y", nrow = 1) +
  xlab("Diversity index") + ylab("Index value") + 
  geom_errorbar(data = hier.div2, aes(x=type, ymin=total-1.96*SE, ymax=total+1.96*SE, 
                                      width=0.3, alpha = 0.4)) +
  thm + theme(legend.position = "none", axis.text.x = element_text(size = 8))
gp.hier

#############################################################################################
## BETA DIVERSITY AND BETA DEVIATION (NULL MODELS)
#############################################################################################
## CUSTOM CODE MODIFIED FROM KRAFT ET AL. 2011, TELLO ET AL.; downloaded from: 
# http://science.sciencemag.org/content/333/6050/1755/tab-figures-data
# http://jsebastiantello.weebly.com/r-code.html
#############################################################################################
## USE COMMON SPECIES POOL FOR FRAGMENT & CONTROL
#############################################################################################
mat.kod = as.matrix(cast(kod.dat, type ~ species, value = 'height_m', fun = length))
spp.pool = colSums(mat.kod)[-ncol(mat.kod)]

## BRAY-CURTIS DISSIMILARITY
## FRAGMENT
## Create random matrices 
rand.f <- assemble.from.pool.randA(compo=mat.f, fix.local.abund=T, fix.rSAD=T, 
                                   rand.N=999, save.output=FALSE, pool.type = "user.defined",
                                   spp.pool = spp.pool)
null.f <- rand.f$rand.datasets

# Partition diversity for the empirical and null matrices
bray.null.f <- compo.dists.nullM(compo=mat.f, null.matrices=null.f, 
                                 dist.method="bray", binary.data=FALSE, print.progress=FALSE)

obsvec.bray.f = melt(as.matrix(bray.null.f$empirical.values), varnames = c("row", "col"))
obsvec.bray.f = obsvec.bray.f[obsvec.bray.f$row != obsvec.bray.f$col,]

sesvec.bray.f = melt(bray.null.f$summary.table[,,"SES"], varnames = c("row", "col"))
sesvec.bray.f = sesvec.bray.f[sesvec.bray.f$row != sesvec.bray.f$col,]

## CONTROL
## Create random matrices 
rand.c <- assemble.from.pool.randA(compo=mat.c, fix.local.abund=T, fix.rSAD=T, 
                                   rand.N=999, save.output=FALSE, pool.type = "user.defined",
                                   spp.pool = spp.pool)
null.c <- rand.c$rand.datasets

bray.null.c <- compo.dists.nullM(compo=mat.c, null.matrices=null.c, 
                                 dist.method="bray", binary.data=FALSE, print.progress=T)
obsvec.bray.c = melt(as.matrix(bray.null.c$empirical.values), varnames = c("row", "col"))
obsvec.bray.c = obsvec.bray.c[obsvec.bray.c$row != obsvec.bray.c$col,]

sesvec.bray.c = melt(bray.null.c$summary.table[,,"SES"], varnames = c("row", "col"))
sesvec.bray.c = sesvec.bray.c[sesvec.bray.c$row != sesvec.bray.c$col,]

## combine fragment and control results
bray.dat = data.frame(bray.obs = c(obsvec.bray.f[,3], obsvec.bray.c[,3]), 
                      bray.disp = c(sesvec.bray.f[,3], sesvec.bray.c[,3]), 
                      type = rep(c("Fragment","Contiguous"), 
                                 c(nrow(sesvec.bray.f),nrow(sesvec.bray.c))))
## models
mod.bcraw = lm(bray.obs ~ type, data = bray.dat)
summary(mod.bcraw)

mod.bcdisp = lm(bray.disp ~ type, data = bray.dat)
summary(mod.bcdisp)

## beta observed
gp.all.bray.obs = ggplot(data = bray.dat, aes(x=type, y=bray.obs)) +
  geom_boxplot() + thm + geom_abline(intercept = 0, slope = 0, linetype = "dashed") +
  xlab("") + ggtitle(expression(paste("a. ",beta,"-diversity (Bray-Curtis)"))) + 
  theme_classic() + theme(axis.title = element_blank(), axis.text = element_text(size = 12)) +
  geom_signif(comparisons = list(c("Contiguous", "Fragment")), map_signif_level=TRUE)

## beta deviation
gp.all.bray.ses = ggplot(data = bray.dat, aes(x=type, y=bray.disp)) +
  geom_boxplot() + thm + geom_abline(intercept = 0, slope = 0, linetype = "dashed") +
  xlab("") + ggtitle(expression(paste("b. ",beta,"-deviation (Bray-Curtis)"))) + 
  theme_classic() + theme(axis.title = element_blank(), axis.text = element_text(size = 12)) +
  geom_signif(comparisons = list(c("Contiguous", "Fragment")), map_signif_level=TRUE)

##############################################
## SORENSEN'S INDEX
##############################################
## FRAGMENT
## Create random matrices 
rand.f <- assemble.from.pool.randA(compo=mat.f, fix.local.abund=T, fix.rSAD=T, 
                                   rand.N=999, save.output=FALSE, pool.type = "user.defined",
                                   spp.pool = spp.pool)
null.f <- rand.f$rand.datasets

whit.null.f <- compo.dists.nullW(compo=mat.f, null.matrices=null.f, print.progress=T)

obsvec.whit.f = melt(as.matrix(whit.null.f$empirical.values), varnames = c("row", "col"))
obsvec.whit.f = obsvec.whit.f[obsvec.whit.f$row != obsvec.whit.f$col,]

sesvec.whit.f = melt(whit.null.f$summary.table[,,"SES"], varnames = c("row", "col"))
sesvec.whit.f = sesvec.whit.f[sesvec.whit.f$row != sesvec.whit.f$col,]

## CONTROL
## Create random matrices 
rand.c <- assemble.from.pool.randA(compo=mat.c, fix.local.abund=T, fix.rSAD=T, 
                                   rand.N=999, save.output=FALSE, pool.type = "user.defined",
                                   spp.pool = spp.pool)
null.c <- rand.c$rand.datasets

whit.null.c <- compo.dists.nullW(compo=mat.c, null.matrices=null.c, print.progress=T)

obsvec.whit.c = melt(as.matrix(whit.null.c$empirical.values), varnames = c("row", "col"))
obsvec.whit.c = obsvec.whit.c[obsvec.whit.c$row != obsvec.whit.c$col,]

sesvec.whit.c = melt(whit.null.c$summary.table[,,"SES"], varnames = c("row", "col"))
sesvec.whit.c = sesvec.whit.c[sesvec.whit.c$row != sesvec.whit.c$col,]

## combine fragment and control results
whit.dat = data.frame(whit.obs = c(obsvec.whit.f[,3], obsvec.whit.c[,3]), 
                      whit.disp = c(sesvec.whit.f[,3], sesvec.whit.c[,3]), 
                      type = rep(c("Fragment","Contiguous"), c(nrow(sesvec.whit.f),nrow(sesvec.whit.c))))
## models
mod.wraw = aov(whit.obs ~ type, data = whit.dat)
summary(mod.wraw)
TukeyHSD(mod.wraw)

mod.wdisp = aov(wdisp ~ type, data = whit.all.vec)
summary(mod.wdisp)
TukeyHSD(mod.wdisp)## beta observed
gp.all.whit.obs = ggplot(data = whit.dat, aes(x=type, y=whit.obs)) +
  geom_boxplot() + thm + geom_abline(intercept = 0, slope = 0, linetype = "dashed") +
  xlab("") + ggtitle(expression(paste("c. ",beta,"-diversity (Sorensen)"))) + 
  theme_classic() + theme(axis.title = element_blank(), axis.text = element_text(size = 12)) +
  geom_signif(comparisons = list(c("Contiguous", "Fragment")), map_signif_level=TRUE)

## beta deviation
gp.all.whit.ses = ggplot(data = whit.dat, aes(x=type, y=whit.disp)) +
  geom_boxplot() + thm + geom_abline(intercept = 0, slope = 0, linetype = "dashed") +
  xlab("") + ggtitle(expression(paste("d. ",beta,"-deviation (Sorensen)"))) + 
  theme_classic() + theme(axis.title = element_blank(), axis.text = element_text(size = 12)) +
  geom_signif(comparisons = list(c("Contiguous", "Fragment")), map_signif_level=TRUE)

## Figure 2
gp.bray.whit = grid.arrange(arrangeGrob(gp.all.bray.obs, gp.all.bray.ses, 
                                        gp.all.whit.obs, gp.all.whit.ses, ncol=2),
                            left = textGrob(expression(paste("More similar" %<->% "Less similar")), 
                                            rot = 90, vjust = 1, gp=gpar(fontsize=16)))

###############################################################################################
## Variance partitioning of alpha diversity
###############################################################################################
pcnm.all = pcnm(dist(dat.loc[,c('xcoord','ycoord')]/1000))
pcnm.all.vec = pcnm.all$vectors[, which(pcnm.all$values>0)]
fwd.sel.rr = forward.sel(as.matrix(div.full.2$rare), pcnm.all.vec)
pcnm.sel.rr = pcnm.all.vec[,colnames(pcnm.all.vec) %in% fwd.sel.rr$variables]
vp.all.rr = varpart(as.matrix(div.full.2$rare), as.factor(div.full.2$'type'),  
                     pcnm.sel.rr, div.full.2[,c("MAP","elevation","CN")])

rda.s.rr = rda(as.matrix(div.full.2$rare) ~ pcnm.sel.rr, data = dat.loc)
RsquareAdj(rda.s.rr)

rda.clim.rr = rda(as.matrix(div.full.2$rare) ~ MAP + elevation + CN, data = dat.loc)
RsquareAdj(rda.clim.rr)

rda.frag.rr = rda(as.matrix(div.full.2$rare) ~ type , data = dat.loc)
RsquareAdj(rda.frag.rr)

## plot varpart
df.venn <- data.frame(x = c(0, 1.25, -1.25),
                      y = c(1, -0.8, -0.8))
df.vp.rr = data.frame(rsq = round(vp.all.rr$part$indfract[-8,3], digits = 2),
                   x = c(0, 1.5, -1.4, 0.75, 0, -0.8, 0),
                   y = c(1.5, -1, -1, 0.5, -1.5, 0.5, 0))
df.vp.rr$rsq = ifelse(df.vp.rr$rsq<0, NA, df.vp.rr$rsq)
vars = data.frame(vars = c("F","S","E"),
                  x = c(0, 1.5, -1.5), y = c(2.5, -2, -2))
gp.vp.rr = ggplot(df.venn, aes(x0 = x, y0 = y, r = 2)) +
  geom_circle(alpha = .3, size = 1) + theme_void() + 
  coord_fixed() + ggtitle(expression(paste("a. ",alpha,"-diversity"))) +
  geom_text(aes(label="Residuals=0.44", x = 0, y = -3), size = 4.5) +
  theme(legend.position = 'bottom', plot.title = element_text(hjust = 0.5, size = 16)) +
  annotate("text", x = df.vp.rr$x, y = df.vp.rr$y, label = df.vp.rr$rsq, size = 5) +
  annotate("text", x = vars$x, y = vars$y, label = vars$vars, size = 5) 

#############################################################################################
## COMPOSITIONAL VARIATION
#############################################################################################
## PERMANOVA
ado.tree1 = adonis(mat.all ~ type, datamethod="bray") 
ado.tree2 = adonis(mat.all ~ type, method="jaccard")

## NMDS
nmds.tree = metaMDS(vegdist(mat.all), distance="bray")

## Figure S4
png("FigS4-NMDS-BW.png", width = 6, height = 5, units = 'in', res = 600)

op  =  ordiplot(nmds.tree, type ="none") 
points(op, what="sites", pch=as.numeric(factor(type)))
ordihull(nmds.tree, groups=type, draw="polygon",label=F)
#ordispider(nmds.tree, display="sites", groups=sort(type), spiders = "centroid")
legend(x=0.4, y=0.4, legend = c("Contiguous","Fragment"), pch = c(1,2), bty="n")
text(x=0.6, y=-0.4, labels = "Stress = 0.17")
dev.off()

###########################################################################################
## VARIANCE PARTITIONING for composition
###########################################################################################
pcnm.all = pcnm(dist(dat.loc[,c('xcoord','ycoord')]/1000))
pcnm.all.vec = pcnm.all$vectors[, which(pcnm.all$values>0)]
fwd.sel.vp = forward.sel(mat.all, pcnm.all.vec)
pcnm.sel = pcnm.all.vec[,colnames(pcnm.all.vec) %in% fwd.sel.vp$variables]
vp.all.mat = varpart(mat.all, dat.loc[,c('type')],  pcnm.sel, dat.loc[,c("MAP","elevation","CN")])
plot(vp.all.mat, cutoff=0.005)

rda.clim = rda(mat.all ~ MAP + elevation + CN, data = dat.loc)
RsquareAdj(rda.clim)

rda.frag = rda(mat.all ~ type , data = dat.loc)
RsquareAdj(rda.frag)

rda.sp = rda(mat.all ~ pcnm.sel, data = dat.loc)
RsquareAdj(rda.sp)

rda.tot = rda(mat.all ~ pcnm.sel + type + MAP + elevation + CN, data = dat.loc)
RsquareAdj(rda.tot)
anova(rda.fclim, rda.tot)

# plot
df.venn = data.frame(x = c(0, 1.25, -1.25), y = c(1, -0.8, -0.8))
df.vp = data.frame(rsq = round(vp.all.mat$part$indfract[-8,3], digits = 2),
                   x = c(0, 1.5, -1.4, 0.75, 0, -0.8, 0),
                   y = c(1.5, -1, -1, 0.5, -1.5, 0.5, 0))
df.vp$rsq = ifelse(df.vp$rsq<0, NA, df.vp$rsq)
vars = data.frame(vars = c("F","S","E"),
                  x = c(0, 1.5, -1.5), y = c(2.5, -2, -2))
gp.vp.full = ggplot(df.venn, aes(x0 = x, y0 = y, r = 2)) +
  geom_circle(alpha = .3, size = 1) + theme_void() + 
  coord_fixed() + ggtitle("b. Species composition") +
  geom_text(aes(label="Residuals=0.41", x = 0, y = -3), size = 4.5) +
  theme(legend.position = 'bottom', plot.title = element_text(hjust = 0.5, size = 16)) +
  annotate("text", x = df.vp$x, y = df.vp$y, label = df.vp$rsq, size = 5) +
  annotate("text", x = vars$x, y = vars$y, label = vars$vars, size = 5) 

gp.rda.vp = grid.arrange(arrangeGrob(gp.vp.rr, gp.vp.full, ncol=2))

####################################################################################################
## PARTIAL CORRELOGRAMS: wtI staistic
####################################################################################################
loc.env.all = unique(kod.dat[,c('plot.all',"xcoord","ycoord","elevation","MAP")])

manr.map = mantel(dist(loc.env.all[c("MAP")]) ~ dist(loc.env.all[c("ycoord","xcoord")]/1000))
manr.ele = mantel(dist(loc.env.all[c("elevation")]) ~ dist(loc.env.all[c("ycoord","xcoord")]/1000))

## BETA DIVERSITY: BRAY CURTIS
pman.dist = pmgram(vegdist(mat.all, method = "bray"), 
                   space = dist(loc.env.all[c("ycoord","xcoord")]/1000),
                   partial = dist(loc.env.all[c("MAP","elevation")]))
pmanr.dist = mantel.partial(vegdist(mat.all, method = "bray"), dist(loc.env.all[c("ycoord","xcoord")]/1000),
                            dist(scale(loc.env.all[c("MAP","elevation")])))
pman.env = pmgram(vegdist(mat.all, method = "bray"), dist(scale(loc.env.all[c("MAP","elevation")])),
                  partial = dist(loc.env.all[c("ycoord","xcoord")]/1000))
pmanr.env = mantel.partial(vegdist(mat.all, method = "bray"), dist(scale(loc.env.all[c("MAP","elevation")])),
                           dist(loc.env.all[c("ycoord","xcoord")]/1000))
pman.map = pmgram(vegdist(mat.all, method = "bray"), dist(loc.env.all[c("MAP")]),
                  partial = dist(loc.env.all[c("ycoord","xcoord")]/1000))
pmanr.map = mantel.partial(vegdist(mat.all, method = "bray"), dist(loc.env.all[c("MAP")]),
                           dist(loc.env.all[c("ycoord","xcoord")]/1000))
pman.ele = pmgram(vegdist(mat.all, method = "bray"), dist(loc.env.all[c("elevation")]),
                  partial = dist(loc.env.all[c("ycoord","xcoord")]/1000))
pmanr.ele = mantel.partial(vegdist(mat.all, method = "bray"), dist(loc.env.all[c("elevation")]),
                           dist(loc.env.all[c("ycoord","xcoord")]/1000))

## plot correlogram
mantel.all = as.data.frame(cbind(rbind(pman.dist$mgram, pman.ele$mgram, pman.map$mgram)))
mantel.all$vars = as.factor(c(rep("c. Space", nrow(pman.dist$mgram)), 
                              rep("d. Elevation", nrow(pman.ele$mgram)), rep("e. MAP", nrow(pman.map$mgram))))
mantel.all$sig = with(mantel.all, ifelse(pval<0.05, 16, 1))
vars = levels(mantel.all$vars)

gp.mgram.all = ggplot(data = mantel.all, aes(x=round(lag, digits=0), y=wtI)) +
  geom_point(size = 1.5, shape=mantel.all$sig) + geom_line() + 
  facet_wrap(~ vars, ncol = 3, scales = "free") + 
  xlab("Euclidean distance") + ylab("Correlation (wtI)") + 
  geom_hline(yintercept = 0, linetype = "dashed") + thm +
  theme(strip.text.x = element_text(size = 16))

#Figure 3: combine varpart + pmgram
gp.vp.mcg = grid.arrange(arrangeGrob(gp.rda.vp, gp.mgram.all, nrow=2))
ggsave(gp.vp.mcg, file="Fig3-vp-pmcg.png", height = 6.5, width = 7.5, units="in", dpi=600)

#############################################################################################
## PER HABITAT: DETECTING SIGNATURES OF ASSEMBLY MECHANISMS
###########################################################################################
###################################################################################################
## MANTEL CORRELOGRAMS WITHIN FRAGMENT AND CONTIGUOUS FOREST
###################################################################################################
## environmental variation for fragment and control
loc.env.c = unique(tree.c[,c('plot.all',"xcoord","ycoord","elevation","MAP")])
loc.env.f = unique(tree.f[,c('plot.all',"xcoord","ycoord","elevation","MAP")])

## composition vs. geodist
mantel.bray.c = pmgram(vegdist(mat.c, method = "bray"), 
                       dist(loc.env.c[, c("xcoord","ycoord")]/1000),
                       partial = vegdist(scale(loc.env.c[c("MAP","elevation")]),
                                         method = "gower"))
mantel.bray.f = pmgram(vegdist(mat.f, method = "bray"), 
                       dist(loc.env.f[, c("xcoord","ycoord")]/1000),
                       partial = vegdist(scale(loc.env.f[c("MAP","elevation")]), 
                                         method = 'gower'))
plot(mantel.bray.c)
plot(mantel.bray.f)

## composition vs. env dist
mantel.bray.map.c = pmgram(vegdist(mat.c, method = "bray"), 
                           dist(loc.env.c[c("MAP")]), 
                          partial = dist(loc.env.c[, c("xcoord","ycoord")]/1000))
mantel.bray.ele.c = pmgram(vegdist(mat.c, method = "bray"), dist(loc.env.c[c("elevation")]),
                           partial = dist(loc.env.c[, c("xcoord","ycoord")]/1000))
mantel.bray.map.f = pmgram(vegdist(mat.f, method = "bray"), dist(loc.env.f[c("MAP")]),
                          partial = dist(loc.env.f[, c("xcoord","ycoord")]/1000))
mantel.bray.ele.f = pmgram(vegdist(mat.f, method = "bray"), dist(loc.env.f[c("elevation")]),
                           partial = dist(loc.env.f[, c("xcoord","ycoord")]/1000))
## env vs geodist
mantel.map.c = pmgram(dist(loc.env.c[c("MAP")]), dist(loc.env.c[c("ycoord","xcoord")]/1000),
                      partial = dist(loc.env.c[c("elevation")]))
mantel.map.f = pmgram(dist(loc.env.f[c("MAP")]), dist(loc.env.f[c("ycoord","xcoord")]/1000),
                     partial = dist(loc.env.f[c("elevation")]))
mantel.ele.c = pmgram(dist(loc.env.c[c("elevation")]), dist(loc.env.c[c("ycoord","xcoord")]/1000),
                     partial = dist(loc.env.c[c("MAP")]))
mantel.ele.f = pmgram(dist(loc.env.f[c("elevation")]), dist(loc.env.f[c("ycoord","xcoord")]/1000),
                     partial = dist(loc.env.f[c("MAP")]))

manr.map.c = mantel(dist(loc.env.c[c("MAP")]) ~ dist(loc.env.c[c("ycoord","xcoord")]/1000))
manr.map.f = mantel(dist(loc.env.f[c("MAP")]) ~ dist(loc.env.f[c("ycoord","xcoord")]/1000))

manr.ele.c = mantel(dist(loc.env.c[c("elevation")]) ~ dist(loc.env.c[c("ycoord","xcoord")]/1000))
manr.ele.f = mantel(dist(loc.env.f[c("elevation")]) ~ dist(loc.env.f[c("ycoord","xcoord")]/1000))

## plot correlograms per habitat
## make vartable
# vs. geodist
times = nrow(mantel.bray.c$mgram)+nrow(mantel.bray.f$mgram) 
mantel.hab = as.data.frame(rbind(mantel.bray.c$mgram, mantel.bray.f$mgram, 
                                 mantel.map.c$mgram, mantel.map.f$mgram, 
                                 mantel.ele.c$mgram, mantel.ele.f$mgram))
mantel.hab$vars = as.factor(c(rep("Composition", times), rep("MAP", times), 
                              rep("Elevation", times)))
mantel.hab$type = as.factor(rep(c(rep("Contiguous",9), rep("Fragment",8)), times = 3))
mantel.hab$sig = with(mantel.hab, ifelse(is.na(pval), 1, ifelse(pval<0.05,16,1)))
mantel.hab$wtI = with(mantel.hab, ifelse(is.na(wtI), 0, wtI))

vars = data.frame(expand.grid(levels(mantel.hab$type), levels(mantel.hab$vars)))
colnames(vars) = c("type", "vars")
dat_text = data.frame(x=rep(0.4,6), y =1, lab = c('a.','b.','c.','d.','e.','f.'), vars)

gp.mgram.geo = ggplot(data = mantel.hab, aes(x=round(lag, digits=0), y=wtI)) +
  geom_point(size = 1.5, shape = mantel.hab$sig) + geom_line() + 
  facet_grid(vars ~ type, scales = "free") + 
  xlab("Geographical distance (km)") + ylab("Correlation (wtI)") + 
  geom_hline(yintercept = 0, linetype = "dashed") + thm +
  geom_text(data=dat_text, aes(x=x, y=y, label = lab), show.legend = FALSE, colour="black") 
gp.mgram.geo
## vs. envdist
## plot correlograms
mantel.bray.env = as.data.frame(cbind(rbind(mantel.bray.map.c$mgram, mantel.bray.ele.c$mgram, 
                                            mantel.bray.map.f$mgram, mantel.bray.ele.f$mgram)))
mantel.bray.env$vars = as.factor(c(rep("MAP", nrow(mantel.bray.map.c$mgram)), 
                                   rep("Elevation", nrow(mantel.bray.ele.c$mgram)),
                                   rep("MAP", nrow(mantel.bray.map.f$mgram)),
                                   rep("Elevation", nrow(mantel.bray.ele.f$mgram))))
mantel.bray.env$type = as.factor(c(rep("Contiguous", 2*nrow(mantel.bray.map.c$mgram)),
                                   rep("Fragment", 2*nrow(mantel.bray.map.f$mgram))))
mantel.bray.env$point = with(mantel.bray.env, ifelse(pval<0.05,16,1))

mantel.map = subset(mantel.bray.env, vars=="MAP")
mantel.ele = subset(mantel.bray.env, vars=="Elevation")

dat_text = data.frame(type = levels(mantel.ele$type), lab = c('g.','h.'))
gp.mgram.ele = ggplot(data = mantel.ele, aes(x=round(lag, digits=0), y=wtI)) +
  geom_point(size = 1.5, shape = mantel.ele$point) + geom_line() + 
  facet_grid(~ type, scales = "free") + 
  xlab("Difference in elevation") + ylab("Correlation (wtI)") + 
  geom_hline(yintercept = 0, linetype = "dashed") + thm +
  geom_text(data=dat_text, aes(x=0.5, y=0.5, label = lab), show.legend = FALSE, colour="black") 
gp.mgram.ele

dat_text = data.frame(type = levels(mantel.map$type), lab = c('i.','j.'))
gp.mgram.map = ggplot(data = mantel.map, aes(x=round(lag, digits=0), y=wtI)) +
  geom_point(size = 1.5, shape = mantel.map$point) + geom_line() + 
  facet_grid(~ type, scales = "free") + 
  xlab("Difference in MAP") + ylab("Correlation (wtI)") + 
  geom_hline(yintercept = 0, linetype = "dashed") + thm +
  geom_text(data=dat_text, aes(x=0.5, y=0.5, label = lab), show.legend = FALSE, colour="black") 
gp.mgram.map

vars = data.frame(expand.grid(levels(mantel.bray.env$type), levels(mantel.bray.env$vars)))
colnames(vars) = c("type","vars")
dat_text = data.frame(x=rep(0.5,4), y =0.5, lab = c('g.','h.','i.','j.'), vars)

gp.mgram.env = ggplot(data = mantel.bray.env, aes(x=round(lag, digits=0), y=wtI)) +
  geom_point(size = 1.5, shape = mantel.bray.env$point) + geom_line() + 
  facet_grid(vars ~ type, scales = "free") + 
  xlab("Environmental distance") + ylab("Correlation (wtI)") + 
  geom_hline(yintercept = 0, linetype = "dashed") + thm +
  geom_text(data=dat_text, aes(x=x, y=y, label = lab), show.legend = FALSE, colour="black") 

# Figure 4
gp.mgram.full = grid.arrange(gp.mgram.geo, gp.mgram.ele, gp.mgram.map,
                             heights=unit(c(4.5,2.25,2.25), c("in")))

## BETA DEVIATION: BRAY CURTIS
## composition vs. geodist
mantel.bdev.c = pmgram(as.dist(bray.null.c$summary.table[,,"SES"]), 
                       dist(loc.env.c[, c("xcoord","ycoord")]/1000),
                       partial = vegdist(scale(loc.env.c[c("elevation","MAP")]), method = "gower"))
mantel.bdev.f = pmgram(as.dist(bray.null.f$summary.table[,,"SES"]), 
                       dist(loc.env.f[, c("xcoord","ycoord")]/1000),
                      partial = vegdist(scale(loc.env.f[c("elevation","MAP")]), method = "gower"))
## composition vs. env dist
mantel.bdev.map.c = pmgram(as.dist(bray.null.c$summary.table[,,"SES"]), 
                           dist(loc.env.c[c("MAP")]),# method = "gower"),
                           partial = dist(loc.env.c[, c("xcoord","ycoord")]/1000))
mantel.bdev.ele.c = pmgram(as.dist(bray.null.c$summary.table[,,"SES"]), 
                           dist(loc.env.c[c("elevation")]),
                           partial = dist(loc.env.c[, c("xcoord","ycoord")]/1000))
mantel.bdev.map.f = pmgram(as.dist(bray.null.f$summary.table[,,"SES"]), 
                           dist(loc.env.f[c("MAP")]),
                           partial = dist(loc.env.f[, c("xcoord","ycoord")]/1000))
mantel.bdev.ele.f = pmgram(as.dist(bray.null.f$summary.table[,,"SES"]), 
                           dist(loc.env.f[c("elevation")]),
                           partial = dist(loc.env.f[, c("xcoord","ycoord")]/1000))

## plot correlograms per habitat
## make vartable
mantel.hab = as.data.frame(rbind(mantel.bdev.c$mgram, mantel.bdev.f$mgram))
mantel.hab$type = as.factor(c(rep("Contiguous",9), rep("Fragment",8)))
mantel.hab$sig = with(mantel.hab, ifelse(is.na(pval), 1, ifelse(pval<0.05, 16, 1)))
mantel.hab$wtI = with(mantel.hab, ifelse(is.na(wtI), 0, wtI))

dat_text = data.frame(x=rep(0.5,2), y =1, lab = c('a.','b.'), type=levels(mantel.hab$type))

gp.mgram.bdev = ggplot(data = mantel.hab, aes(x=round(lag, digits=0), y=wtI)) +
  geom_point(size = 1.5, shape = mantel.hab$sig) + geom_line() + 
  facet_grid(~type, scales = "free") + 
  xlab("Geographical distance (km)") + ylab("Correlation (wtI)") + 
  # geom_errorbar(aes(x=round(lag, digits=0), ymin = llim, ymax = ulim), 
  #               width=0.3, alpha = 0.4) +
  geom_hline(yintercept = 0, linetype = "dashed") + thm +
  geom_text(data=dat_text, aes(x=x, y=y, label = lab), show.legend = FALSE, colour="black") 

## vs. envdist
mantel.bdev.env = as.data.frame(cbind(rbind(mantel.bdev.map.c$mgram, mantel.bdev.ele.c$mgram, 
                                            mantel.bdev.map.f$mgram, mantel.bdev.ele.f$mgram)))
mantel.bdev.env$vars = as.factor(c(rep("MAP", nrow(mantel.bdev.map.c$mgram)), 
                                   rep("Elevation", nrow(mantel.bdev.ele.c$mgram)),
                                   rep("MAP", nrow(mantel.bdev.map.f$mgram)),
                                   rep("Elevation", nrow(mantel.bdev.ele.f$mgram))))
mantel.bdev.env$type = as.factor(c(rep("Contiguous", 2*nrow(mantel.bdev.map.c$mgram)),
                                   rep("Fragment", 2*nrow(mantel.bdev.map.f$mgram))))
mantel.bdev.env$point = with(mantel.bdev.env, ifelse(pval<0.05,16,1))

vars = data.frame(expand.grid(levels(mantel.bdev.env$type), levels(mantel.bdev.env$vars)))
colnames(vars) = c("type","vars")
dat_text = data.frame(x=rep(0.5,4), y=0.4, lab = c('c.','d.','e.','f.'), vars)

gp.mgram.bdev.env = ggplot(data = mantel.bdev.env, aes(x=round(lag, digits=0), y=wtI)) +
  geom_point(size = 1.5, shape = mantel.bdev.env$point) + geom_line() + 
  facet_grid(vars ~ type, scales = "free") + #ylim(c(-0.5,0.5)) +
  xlab("Environmental distance") + ylab("Correlation (wtI)") + 
  geom_hline(yintercept = 0, linetype = "dashed") + thm +
  geom_text(data=dat_text, aes(x=x, y=y, label = lab), show.legend = FALSE, colour="black") 
# Figure S5
gp.mgram.bdev.full = grid.arrange(arrangeGrob(gp.mgram.bdev, gp.mgram.bdev.env))

## BETA DIVERSITY: SORENSEN
## composition vs. geodist
mantel.whit.c = pmgram(vegdist(mat.c, method = "bray", binary=T), 
                       dist(loc.env.c[, c("xcoord","ycoord")]/1000),
                      partial = vegdist(scale(loc.env.c[c("elevation","MAP")]), method = "gower"))
mantel.whit.f = pmgram(vegdist(mat.f, method = "bray", binary=T), 
                       dist(loc.env.f[, c("xcoord","ycoord")]/1000),
                      partial = vegdist(scale(loc.env.f[c("elevation","MAP")]), method = "gower"))
## composition vs. env dist
mantel.whit.map.c = pmgram(vegdist(mat.c, method = "bray", binary=T), 
                           dist(loc.env.c[c("MAP")]),
                           partial = dist(loc.env.c[, c("xcoord","ycoord")]/1000))
mantel.whit.ele.c = pmgram(vegdist(mat.c, method = "bray", binary=T), 
                           dist(loc.env.c[c("elevation")]),
                           partial = dist(loc.env.c[, c("xcoord","ycoord")]/1000))
mantel.whit.map.f = pmgram(vegdist(mat.f, method = "bray", binary=T), 
                           dist(loc.env.f[c("MAP")]),
                           partial = dist(loc.env.f[, c("xcoord","ycoord")]/1000))
mantel.whit.ele.f = pmgram(vegdist(mat.f, method = "bray", binary=T), 
                           dist(loc.env.f[c("elevation")]),
                           partial = dist(loc.env.f[, c("xcoord","ycoord")]/1000))
## plot correlograms per habitat
## make vartable
mantel.hab = as.data.frame(rbind(mantel.whit.c$mgram, mantel.whit.f$mgram))
mantel.hab$type = as.factor(c(rep("Contiguous",9), rep("Fragment",8)))
mantel.hab$sig = with(mantel.hab, ifelse(is.na(pval), 1, ifelse(pval<0.05, 16, 1)))
mantel.hab$wtI = with(mantel.hab, ifelse(is.na(wtI), 0, wtI))

dat_text = data.frame(x=rep(0.5,2), y =0.6, lab = c('a.','b.'), type=levels(mantel.hab$type))

gp.mgram.sor = ggplot(data = mantel.hab, aes(x=round(lag, digits=0), y=wtI)) +
  geom_point(size = 1.5, shape = mantel.hab$sig) + geom_line() + 
  facet_grid(~type, scales = "free") + 
  xlab("Geographical distance (km)") + ylab("Correlation (wtI)") + 
  geom_hline(yintercept = 0, linetype = "dashed") + thm +
  geom_text(data=dat_text, aes(x=x, y=y, label = lab), show.legend = FALSE, colour="black") 

## vs. envdist
mantel.whit.env = as.data.frame(cbind(rbind(mantel.whit.map.c$mgram, mantel.whit.ele.c$mgram, 
                                            mantel.whit.map.f$mgram, mantel.whit.ele.f$mgram)))
mantel.whit.env$vars = as.factor(c(rep("MAP", nrow(mantel.whit.map.c$mgram)), 
                                   rep("Elevation", nrow(mantel.whit.ele.c$mgram)),
                                   rep("MAP", nrow(mantel.whit.map.f$mgram)),
                                   rep("Elevation", nrow(mantel.whit.ele.f$mgram))))
mantel.whit.env$type = as.factor(c(rep("Contiguous", 2*nrow(mantel.whit.map.c$mgram)),
                                   rep("Fragment", 2*nrow(mantel.whit.map.f$mgram))))
mantel.whit.env$point = with(mantel.whit.env, ifelse(pval<0.05,16,1))

vars = data.frame(expand.grid(levels(mantel.whit.env$type), levels(mantel.whit.env$vars)))
colnames(vars) = c("type","vars")
dat_text = data.frame(x=rep(0.5,4), y=0.4, lab = c('c.','d.','e.','f.'), vars)

gp.mgram.sor.env = ggplot(data = mantel.whit.env, aes(x=round(lag, digits=0), y=wtI)) +
  geom_point(size = 1.5, shape = mantel.whit.env$point) + geom_line() + 
  facet_grid(vars ~ type, scales = "free") + #ylim(c(-0.5,0.5)) +
  xlab("Environmental distance") + ylab("Correlation coefficient (wtI)") + 
  geom_hline(yintercept = 0, linetype = "dashed") + thm +
  geom_text(data=dat_text, aes(x=x, y=y, label = lab), show.legend = FALSE, colour="black") 

# Figure S6
gp.mgram.sor.full = grid.arrange(arrangeGrob(gp.mgram.sor, gp.mgram.sor.env))

## BETA DEVIATION: SORENSEN
## composition vs. geodist
mantel.wdev.c = pmgram(as.dist(whit.null.c$summary.table[,,"SES"]), 
                       dist(loc.env.c[, c("xcoord","ycoord")]/1000),
                      partial = vegdist(scale(loc.env.c[c("elevation","MAP")]), method = "gower"))
mantel.wdev.f = pmgram(as.dist(whit.null.f$summary.table[,,"SES"]), 
                       dist(loc.env.f[, c("xcoord","ycoord")]/1000),
                      partial = vegdist(scale(loc.env.f[c("elevation","MAP")]), method = 'gower'))
## composition vs. env dist
mantel.wdev.map.c = pmgram(as.dist(whit.null.c$summary.table[,,"SES"]), 
                           dist(loc.env.c[c("MAP")]),
                           partial = dist(loc.env.c[, c("xcoord","ycoord")]/1000))
mantel.wdev.ele.c = pmgram(as.dist(whit.null.c$summary.table[,,"SES"]), dist(loc.env.c[c("elevation")]),
                           partial = dist(loc.env.c[, c("xcoord","ycoord")]/1000))
mantel.wdev.map.f = pmgram(as.dist(whit.null.f$summary.table[,,"SES"]), dist(loc.env.f[c("MAP")]),
                           partial = dist(loc.env.f[, c("xcoord","ycoord")]/1000))
mantel.wdev.ele.f = pmgram(as.dist(whit.null.f$summary.table[,,"SES"]), dist(loc.env.f[c("elevation")]),
                           partial = dist(loc.env.f[, c("xcoord","ycoord")]/1000))
## plot correlograms per habitat
## make vartable
mantel.hab = as.data.frame(rbind(mantel.wdev.c$mgram, mantel.wdev.f$mgram))
mantel.hab$type = as.factor(c(rep("Contiguous",9), rep("Fragment",8)))
mantel.hab$sig = with(mantel.hab, ifelse(is.na(pval), 1, ifelse(pval<0.05, 16, 1)))
mantel.hab$wtI = with(mantel.hab, ifelse(is.na(wtI), 0, wtI))

dat_text = data.frame(x=rep(0.5,2), y =0.6, lab = c('a.','b.'), type=levels(mantel.hab$type))

gp.mgram.wdev = ggplot(data = mantel.hab, aes(x=round(lag, digits=0), y=wtI)) +
  geom_point(size = 1.5, shape = mantel.hab$sig) + geom_line() + 
  facet_grid(~type, scales = "free") + 
  xlab("Geographical distance (km)") + ylab("Correlation (wtI)") + 
  geom_hline(yintercept = 0, linetype = "dashed") + thm +
  geom_text(data=dat_text, aes(x=x, y=y, label = lab), show.legend = FALSE, colour="black") 

## vs. envdist
mantel.wdev.env = as.data.frame(cbind(rbind(mantel.wdev.map.c$mgram, mantel.wdev.ele.c$mgram, 
                                            mantel.wdev.map.f$mgram, mantel.wdev.ele.f$mgram)))
mantel.wdev.env$vars = as.factor(c(rep("MAP", nrow(mantel.wdev.map.c$mgram)), 
                                   rep("Elevation", nrow(mantel.wdev.ele.c$mgram)),
                                   rep("MAP", nrow(mantel.wdev.map.f$mgram)),
                                   rep("Elevation", nrow(mantel.wdev.ele.f$mgram))))
mantel.wdev.env$type = as.factor(c(rep("Contiguous", 2*nrow(mantel.wdev.map.c$mgram)),
                                   rep("Fragment", 2*nrow(mantel.wdev.map.f$mgram))))
mantel.wdev.env$point = with(mantel.wdev.env, ifelse(pval<0.05,16,1))

vars = data.frame(expand.grid(levels(mantel.wdev.env$type), levels(mantel.wdev.env$vars)))
colnames(vars) = c("type","vars")
dat_text = data.frame(x=rep(0.5,4), y=0.4, lab = c('c.','d.','e.','f.'), vars)

gp.mgram.wdev.env = ggplot(data = mantel.wdev.env, aes(x=round(lag, digits=0), y=wtI)) +
  geom_point(size = 1.5, shape = mantel.wdev.env$point) + geom_line() + 
  facet_grid(vars ~ type, scales = "free") + #ylim(c(-0.5,0.5)) +
  xlab("Environmental distance") + ylab("Correlation (wtI)") + 
  geom_hline(yintercept = 0, linetype = "dashed") + thm +
  geom_text(data=dat_text, aes(x=x, y=y, label = lab), show.legend = FALSE, colour="black") 

# Figure S7
gp.mgram.wdev.full = grid.arrange(arrangeGrob(gp.mgram.wdev, gp.mgram.wdev.env))

###########################################################################################
## WITHIN HABITAT MRM
###########################################################################################
####################################
## BRAY-CURTIS
####################################
## FRAGMENT
env.f = droplevels(subset(dat.loc, type=='FRAGMENT'))
mrm.ele.f = MRM(dist(env.f[, c("elevation")]) ~  
                      dist(env.f[, c("xcoord","ycoord")]/1000), nperm = 1000)
mrm.map.f = MRM(dist(env.f[, c("MAP")]) ~  
                      dist(env.f[, c("xcoord","ycoord")]/1000), nperm = 1000)
# beta diversity
mrm.braw.f.geo = MRM(as.dist(bray.null.f$empirical.values) ~ 
                       dist(env.f[, c("xcoord","ycoord")]/1000), nperm = 1000)
mrm.braw.f.env = MRM(as.dist(bray.null.f$empirical.values) ~ dist(env.f[, c("elevation")]) + 
                       dist(env.f[, c("MAP")]), nperm = 1000) #+ dist(env.f[, c("CN")]
mrm.braw.f.all = MRM(as.dist(bray.null.f$empirical.values) ~ dist(env.f[, c("elevation")]) + 
                   dist(env.f[, c("MAP")]) + 
                     dist(env.f[, c("xcoord","ycoord")]/1000), nperm = 1000)

## variance explained by different fractions
## env alone
mrm.braw.f.all$r.squared[1] - mrm.braw.f.geo$r.squared[1]

## space alone
mrm.braw.f.all$r.squared[1] - mrm.braw.f.env$r.squared[1]

## jointly by env & geo
(mrm.braw.f.geo$r.squared[1] + mrm.braw.f.env$r.squared[1]) - mrm.braw.f.all$r.squared[1]

# beta deviation
mrm.bdev.f.geo = MRM(as.dist(bray.null.f$summary.table[,,"SES"]) ~ 
                       dist(env.f[, c("xcoord","ycoord")]/1000), nperm = 100)
mrm.bdev.f.env = MRM(as.dist(bray.null.f$summary.table[,,"SES"]) ~ dist(env.f[, c("elevation")]) + 
                       dist(env.f[, c("MAP")]), nperm = 100)
mrm.bdev.f.all = MRM(as.dist(bray.null.f$summary.table[,,"SES"]) ~ dist(env.f[, c("elevation")]) + 
                   dist(env.f[, c("MAP")]) + #dist(env.f[, c("CN")]) + 
                     dist(env.f[, c("xcoord","ycoord")]/1000), nperm = 100)

## variance explained by different fractions
## env alone
mrm.bdev.f.all$r.squared[1] - mrm.bdev.f.geo$r.squared[1]

## space alone
mrm.bdev.f.all$r.squared[1] - mrm.bdev.f.env$r.squared[1]

## jointly by env & geo
(mrm.bdev.f.geo$r.squared[1] + mrm.bdev.f.env$r.squared[1]) - mrm.bdev.f.all$r.squared[1]

## CONTIGUOUS
env.c = droplevels(subset(dat.loc, type=='CONTROL'))
mrm.ele.c = MRM(dist(env.c[, c("elevation")]) ~  
                      dist(env.c[, c("xcoord","ycoord")]/1000), nperm = 1000)
mrm.map.c = MRM(dist(env.c[, c("MAP")]) ~  
                      dist(env.c[, c("xcoord","ycoord")]/1000), nperm = 1000)
# beta diversity
mrm.braw.c.geo = MRM(as.dist(bray.null.c$empirical.values) ~ 
                       dist(env.c[, c("xcoord","ycoord")]/1000), nperm = 1000)
mrm.braw.c.env = MRM(as.dist(bray.null.c$empirical.values) ~ dist(env.c[, c("elevation")]) + 
                       dist(env.c[, c("MAP")]) + dist(env.c[, c("CN")]), nperm = 1000)
mrm.braw.c.all = MRM(as.dist(bray.null.c$empirical.values) ~ dist(env.c[, c("elevation")]) + 
                       dist(env.c[, c("MAP")]) + dist(env.c[, c("CN")]) + 
                       dist(env.c[, c("xcoord","ycoord")]/1000), nperm = 1000)
## variance explained by different fractions
## env alone
mrm.braw.c.all$r.squared[1] - mrm.braw.c.geo$r.squared[1]

## space alone
mrm.braw.c.all$r.squared[1] - mrm.braw.c.env$r.squared[1]

## jointly by env & geo
(mrm.braw.c.geo$r.squared[1] + mrm.braw.c.env$r.squared[1]) - mrm.braw.c.all$r.squared[1]

# beta deviation
mrm.bdev.c.geo = MRM(as.dist(bray.null.c$summary.table[,,"SES"]) ~ 
                       dist(env.c[, c("xcoord","ycoord")]/1000), nperm = 100)
mrm.bdev.c.env = MRM(as.dist(bray.null.c$summary.table[,,"SES"]) ~ dist(env.c[, c("elevation")]) + 
                       dist(env.c[, c("MAP")]) + dist(env.c[, c("CN")]), nperm = 100)
mrm.bdev.c.all = MRM(as.dist(bray.null.c$summary.table[,,"SES"]) ~ dist(env.c[, c("elevation")]) + 
                       dist(env.c[, c("MAP")]) + dist(env.c[, c("CN")]) + 
                       dist(env.c[, c("xcoord","ycoord")]/1000), nperm = 100)
## variance explained by different fractions
## env alone
mrm.bdev.c.all$r.squared[1] - mrm.bdev.c.geo$r.squared[1]

## space alone
mrm.bdev.c.all$r.squared[1] - mrm.bdev.c.env$r.squared[1]

## jointly by env & geo
(mrm.bdev.c.geo$r.squared[1] + mrm.bdev.c.env$r.squared[1]) - mrm.bdev.c.all$r.squared[1]

##################################
## SORENSEN's
##################################
## FRAGMENT
# beta diversity
mrm.wraw.f.geo = MRM(as.dist(whit.null.f$empirical.values) ~ 
                       dist(env.f[, c("xcoord","ycoord")]/1000), nperm = 1000)
mrm.wraw.f.env = MRM(as.dist(whit.null.f$empirical.values) ~ dist(env.f[, c("elevation")]) + 
                       dist(env.f[, c("MAP")]) + dist(env.f[, c("CN")]), nperm = 1000)
mrm.wraw.f.all = MRM(as.dist(whit.null.f$empirical.values) ~ dist(env.f[, c("elevation")]) + 
                       dist(env.f[, c("MAP")]) + dist(env.f[, c("CN")]) +
                       dist(env.f[, c("xcoord","ycoord")]/1000), nperm = 1000)

## variance explained by different fractions
## env alone
mrm.wraw.f.all$r.squared[1] - mrm.wraw.f.geo$r.squared[1]

## space alone
mrm.wraw.f.all$r.squared[1] - mrm.wraw.f.env$r.squared[1]

## jointly by env & geo
(mrm.wraw.f.geo$r.squared[1] + mrm.wraw.f.env$r.squared[1]) - mrm.wraw.f.all$r.squared[1]

# beta deviation
mrm.wdev.f.geo = MRM(as.dist(whit.null.f$summary.table[,,"SES"]) ~ 
                       dist(env.f[, c("xcoord","ycoord")]/1000), nperm = 100)
mrm.wdev.f.env = MRM(as.dist(whit.null.f$summary.table[,,"SES"]) ~ dist(env.f[, c("elevation")]) + 
                       dist(env.f[, c("MAP")]) + dist(env.f[, c("CN")]), nperm = 100)
mrm.wdev.f.all = MRM(as.dist(whit.null.f$summary.table[,,"SES"]) ~ dist(env.f[, c("elevation")]) + 
                       dist(env.f[, c("MAP")]) + dist(env.f[, c("CN")]) + 
                       dist(env.f[, c("xcoord","ycoord")]/1000), nperm = 100)
## variance explained by different fractions
## env alone
mrm.wdev.f.all$r.squared[1] - mrm.wdev.f.geo$r.squared[1]

## space alone
mrm.wdev.f.all$r.squared[1] - mrm.wdev.f.env$r.squared[1]

## jointly by env & geo
(mrm.wdev.f.geo$r.squared[1] + mrm.wdev.f.env$r.squared[1]) - mrm.wdev.f.all$r.squared[1]

## CONTIGUOUS
# beta diversity
mrm.wraw.c.geo = MRM(as.dist(whit.null.c$empirical.values) ~ 
                       dist(env.c[, c("xcoord","ycoord")]/1000), nperm = 1000)
mrm.wraw.c.env = MRM(as.dist(whit.null.c$empirical.values) ~ dist(env.c[, c("elevation")]) + 
                       dist(env.c[, c("MAP")]) + dist(env.c[, c("CN")]), nperm = 1000)
mrm.wraw.c.all = MRM(as.dist(whit.null.c$empirical.values) ~ dist(env.c[, c("elevation")]) + 
                       dist(env.c[, c("MAP")]) + dist(env.c[, c("CN")]) + 
                       dist(env.c[, c("xcoord","ycoord")]/1000), nperm = 1000)
## variance explained by different fractions
## env alone
mrm.wraw.c.all$r.squared[1] - mrm.wraw.c.geo$r.squared[1]

## space alone
mrm.wraw.c.all$r.squared[1] - mrm.wraw.c.env$r.squared[1]

## jointly by env & geo
(mrm.wraw.c.geo$r.squared[1] + mrm.wraw.c.env$r.squared[1]) - mrm.wraw.c.all$r.squared[1]

# beta deviation
mrm.wdev.c.geo = MRM(as.dist(whit.null.c$summary.table[,,"SES"]) ~ 
                       dist(env.c[, c("xcoord","ycoord")]/1000), nperm = 100)
mrm.wdev.c.env = MRM(as.dist(whit.null.c$summary.table[,,"SES"]) ~ dist(env.c[, c("elevation")]) + 
                       dist(env.c[, c("MAP")]) + dist(env.c[, c("CN")]), nperm = 100)
mrm.wdev.c.all = MRM(as.dist(whit.null.c$summary.table[,,"SES"]) ~ dist(env.c[, c("elevation")]) + 
                       dist(env.c[, c("MAP")]) + dist(env.c[, c("CN")]) + 
                       dist(env.c[, c("xcoord","ycoord")]/1000), nperm = 100)
## variance explained by different fractions
## env alone
mrm.wdev.c.all$r.squared[1] - mrm.wdev.c.geo$r.squared[1]

## space alone
mrm.wdev.c.all$r.squared[1] - mrm.wdev.c.env$r.squared[1]

## jointly by env & geo
(mrm.wdev.c.geo$r.squared[1] + mrm.wdev.c.env$r.squared[1]) - mrm.wdev.c.all$r.squared[1]

######################
## end
######################