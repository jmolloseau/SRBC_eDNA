library(readxl)
SRBC_macros <- read_excel("data/raw/SRBC_macros.xlsx")

library(dplyr)
library(tidyr)
library(tibble)

#converting excel file from SRBC to site X taxa table w/counts
SRBC_by_site <- SRBC_macros %>%
  group_by(Site, Taxa) %>% 
  summarise(total = sum(Count), .groups = "drop") %>% 
  pivot_wider(
    names_from = Taxa,
    values_from = total,
    values_fill = 0
  )
SRBC_site_macro <- data.frame (SRBC_by_site, row.names="Site")
#presence absence data frame
SRBC_PA <- replace (SRBC_site_macro,SRBC_site_macro>0,1 )
#IBI scores and location info
SRBC_IBIs <- read.csv("data/raw/49_sites.csv")
SRBC_site_IBIs <- data.frame (SRBC_IBIs, row.names="Site")

#alpha diversity
SRBC_shan <- diversity(SRBC_site_macro, index= "shannon")
  SRBC_alpha <- as.data.frame(SRBC_shan)
  SRBC_alpha_fin <- rownames_to_column(SRBC_alpha, var='Site')

ggplot(data=SRBC_alpha_fin, mapping = aes(x=Site, y=SRBC_shan)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(title= "Alpha Diversity",
       x= "Site",
       y= "Shannon Diversity") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))

ggplot(data=SRBC_IBIs, mapping = aes(x=Site, y=IBI_Score)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(title= "Indices of Biotic Integrity Scores",
       x= "Site",
       y= "Score") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))


#beta diversity
library(betapart)
SRBC_beta<- beta.multi.abund(SRBC_site_macro, index.family="bray")
SRBC_beta

SRBC_beta_bray<- SRBC_beta$beta.BRAY

#turnover vs. nestedness
SRBC_TN<- beta.multi(SRBC_PA, index.family="sorensen")
SRBC_TN

SRBC_turn<- SRBC_TN$beta.SIM

SRBC_nest<- SRBC_TN$beta.SNE

SRBC_beta_full<- as.data.frame(cbind(SRBC_beta_bray,SRBC_turn,SRBC_nest))


library(ggplot2)
df_long <- pivot_longer(
  SRBC_beta_full,
  cols = everything(),
  names_to = "category",
  values_to = "value"
)

ggplot(data=df_long, mapping = aes(x=category, y=value)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(title= "Comparison of Beta Diversity, Nestedness, and Turnover",
  x= "Beta Diversity Measure",
  y= "Value") +
theme_minimal()

#NODF
SRBC_NODF <- nestednodf(SRBC_site_macro, order = TRUE, weighted = FALSE, wbinary = FALSE)
SRBC_NODF
plot(SRBC_NODF, col = "black", names = TRUE)
#overall not very nested, which isn't what I would expect (thought that human disturbance would lead to nestedness)
#could be because we are looking at genus, which could be more dependent on dispersal than environmental factors
#would I observe more nestedness if I looked at the family level instead?

#clustering
SRBC.norm <- decostand(SRBC_site_macro, "normalize")
#normalizing because of differences in abundances between sites
SRBC.ch <- vegdist(SRBC.norm, "euc")

SRBC.ch.single <- hclust(SRBC.ch, method = "single")
plot(SRBC.ch.single,labels = rownames(SRBC_site_macro),main = "Chord - Single linkage")


SRBC.ch.complete <- hclust(SRBC.ch, method="complete")
plot(SRBC.ch.complete, labels=rownames(SRBC_site_macro),main= "Chord - Complete Linkage")

SRBC.ch.centroid <- hclust(SRBC.ch, method="centroid")
plot(SRBC.ch.centroid, labels=rownames(SRBC_site_macro),main= "Chord - Centroid Linkage")

SRBC.ch.ward <- hclust(SRBC.ch, method="ward.D2")
plot(SRBC.ch.ward, labels=rownames(SRBC_site_macro),main= "Chord- Ward Linkage")

library(dendextend)

class(SRBC.ch.ward)
dend1 <- as.dendrogram(SRBC.ch.ward)
class(dend1)
dend2 <- as.dendrogram(SRBC.ch.complete)
dend12 <- dendlist(dend1, dend2)
tanglegram(
  untangle(dend12),
  sort = TRUE, common_subtrees_color_branches = TRUE, main_left = "Ward method",
  main_right = "Complete linkage",
  margin_inner = 6,
  margin_outer = 5
)

#PCA
SRBC_PCA <- rda(SRBC_site_macro, scale=TRUE)
summary(SRBC_PCA)
screeplot(SRBC_PCA, bstick = TRUE, npcs = length(SRBC_PCA$CA$eig))

source("cleanplot.pca.R")
par(mfrow = c(1, 2))
cleanplot.pca(SRBC_PCA, ax1=1, ax2=2, scaling = 1) #species scores
cleanplot.pca(SRBC_PCA, ax1=1, ax2=2, scaling = 2) #site scores
#kind of a complete mess

#PERMANOVA
community_ordered <- SRBC_by_site[order(SRBC_by_site$Site),]
  community_matrix <- data.frame(community_ordered, row.names="Site")
metadata_matrix <- SRBC_IBIs[order(SRBC_IBIs$Site),]

SRBC_perm1 <- adonis2(community_matrix ~ IBI_Score+Location+IBI_Score:Location, data=metadata_matrix, permutations = 999, method="bray")
SRBC_perm1

#NMDS
SRBC.nmds <- metaMDS(community_matrix, distance = "bray", try=1000)
SRBC.nmds
SRBC.nmds$stress
stressplot(SRBC.nmds)
plot(SRBC.nmds)
#adding in IBI score and cardinal direction data
scores <-data.frame(SRBC.nmds$points)
#create bins with IBI scores
IBI_bins <- cut(metadata_matrix$IBI_Score, breaks=5)
NMDS_scores <-cbind(SRBC_site_IBIs,scores, IBI_bins)
#create hulls with IBI bins
NMDS_vlowIBI <- NMDS_scores[NMDS_scores$IBI_bins=="23.4,38.2",][chull(NMDS_scores[NMDS_scores$IBI_bins=="23.4,38.2", c("MDS1","MDS2")]),]
NMDS_lowIBI <- NMDS_scores[NMDS_scores$IBI_bins=="38.2,52.9",][chull(NMDS_scores[NMDS_scores$IBI_bins=="38.2,52.9", c("MDS1","MDS2")]),]
NMDS_midIBI <- NMDS_scores[NMDS_scores$IBI_bins=="52.9,67.6",][chull(NMDS_scores[NMDS_scores$IBI_bins=="52.9,67.6", c("MDS1","MDS2")]),]
NMDS_highIBI <- NMDS_scores[NMDS_scores$IBI_bins=="67.6,82.4",][chull(NMDS_scores[NMDS_scores$IBI_bins=="67.6,82.4", c("MDS1","MDS2")]),]
NMDS_vhighIBI <- NMDS_scores[NMDS_scores$IBI_bins=="82.4,97.2",][chull(NMDS_scores[NMDS_scores$IBI_bins=="82.4,97.2", c("MDS1","MDS2")]),]

NMDS_IBI_hull <- rbind(NMDS_vlowIBI,NMDS_lowIBI,NMDS_midIBI,NMDS_highIBI,NMDS_vhighIBI)

SRBC_NMDS_plot <- ggplot() + geom_point(NMDS_scores, mapping=aes(MDS1, MDS2, color=IBI_bins, shape=Location)) +
  geom_polygon(NMDS_IBI_hull, mapping=aes(MDS1, MDS2, fill=IBI_bins, group=IBI_bins), alpha=0.05)
  #creates polygons around the group hulls
SRBC_NMDS_plot
