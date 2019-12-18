###############################
### TMP Metagenome Analysis ###
###############################

# Ryan Blaustein <ryan.blaustein@northwestern.edu>
# last edited: 12-18-2019

#####################

### call packages
library(Nonpareil)
library(ggplot2)
library(gplots)
library(gridExtra)
library(gtools)
library(reshape2)
library(scales)
library(vegan)
library(VennDiagram)


#####
##### LOAD DATA #####
#####

### set wd to */Toothbrush_Microbiome_Project

## taxa tables
# genus table
genus_table = read.table('data/taxa_tables/genus/genus_merge_tb_hmp_oreg_gerl_clean.txt',
                         header = T, sep = '\t', row.names = 1, comment.char = '')
# species table
species_table = read.table('data/taxa_tables/species/species_merge_tb_hmp_oreg_gerl_clean.txt',
                           header = T, sep = '\t', row.names = 1, comment.char = '')
# species table: HMP oral subset
species_hmp_subset = read.table('data/taxa_tables/species/hmp_oral_subset_metaphlan_species.txt',
                                header = T, sep = '\t', row.names = 1, comment.char = '')

## map table
map_table = read.table('data/map_tables/mapping_tb_oreg_gerl_hmp_clean.txt',
                       sep='\t',h=T,row.names=1,check=F,comment='')
# remove neg controls; set Env factor
species_table = species_table[,-c(grep("neg", map_table$Env))]
genus_table = genus_table[,-c(grep("neg", map_table$Env))]
map_table = map_table[-c(grep("neg", map_table$Env)),]
map_table$Env = as.character(map_table$Env)
map_table$Env = factor(map_table$Env,
                       levels = c("toothbrush", "building_dust", "Oral", "Skin", "Gut", "Vaginal"))

## sourcetracker
st_predictions_genus = read.table('data/sourcetracker/sink_predictions.txt',
                                  header = T, sep = '\t', row.names = 1, comment.char = '')

## shortbred
sb_u90 = read.table('data/biobakery/shortbred/sb_CARD_u90_merge.txt',
                    sep='\t', h=T, row.names=1, check=F, comment='',quote="\"")

## halla
# similarity scores: taxa
halla_species = read.table('data/biobakery/halla/similarity_scores_taxa.txt',
                           header = T, sep = '\t', row.names = 1, comment.char = '')
# similarity scores: ARGs
halla_ARGs = read.table('data/biobakery/halla/similarity_scores_ARGs.txt',
                        header = T, sep = '\t', row.names = 1, comment.char = '')

## redcap data
survey_table_clean = read.table('data/metadata/redcap_data_clean.txt',
                                sep='\t',h=T,row.names=1,check=F,comment='')


#####
##### METAGENOME COVERAGE #####
#####

### set wd to */Toothbrush_Microbiome_Project/data/nonpareil

## vector for calling sample outputs
np_list = c('A10.npo','A11.npo','A16.npo','A17.npo','A18.npo',
            'A19.npo','A22.npo','A2.npo','A5.npo',
            'A7.npo','A9.npo','B12.npo','B28.npo','B29.npo',
            'B30.npo','B32.npo','B34.npo','B35a.npo','B35b.npo',
            'B37.npo','B39.npo','B42.npo','B44.npo',
            'B49.npo','B51.npo','B6.npo','B8.npo','C53.npo',
            'C54.npo','C57.npo','C60.npo','C62.npo')

## plot curves (FIGURE S1)
Nonpareil.legend(Nonpareil.curve.batch(np_list, 
                                       plot.opts = list(legend = NA, 
                                                        #plot.model = F, 
                                                        plot.diversity = F,
                                                        main = NA)), cex = 0.5)

## generate summary
summary(Nonpareil.curve.batch(np_list, plot = F))

## compute summary stats: mean, med, sd, se
mean(summary(Nonpareil.curve.batch(np_list, plot = F))[,2])
median(summary(Nonpareil.curve.batch(np_list, plot = F))[,2])
sd(summary(Nonpareil.curve.batch(np_list, plot = F))[,2])
sd(summary(Nonpareil.curve.batch(np_list, plot = F))[,2])/sqrt(32)


######################
##### MICROBIOTA #####
######################

### set wd back to */Toothbrush_Microbiome_Project

###
### Taxa across all sample types
###

### Frequency table
freq_table = matrix(nrow = dim(species_table)[1], ncol = length(levels(map_table$Env)))
freq_table = as.data.frame(freq_table)
colnames(freq_table) = levels(map_table$Env)
rownames(freq_table) = rownames(species_table)

for (x in 1:dim(species_table)[1]) {
  freq_table[x,1] = length(which(species_table[x, grep("toothbrush", map_table$Env)] > 0))/length(grep("toothbrush", map_table$Env))
  freq_table[x,2] = length(which(species_table[x, grep("building_dust", map_table$Env)] > 0))/length(grep("building_dust", map_table$Env))
  freq_table[x,3] = length(which(species_table[x, grep("Oral", map_table$Env)] > 0))/length(grep("Oral", map_table$Env))
  freq_table[x,4] = length(which(species_table[x, grep("Skin", map_table$Env)] > 0))/length(grep("Skin", map_table$Env))
  freq_table[x,5] = length(which(species_table[x, grep("Gut", map_table$Env)] > 0))/length(grep("Gut", map_table$Env))
  freq_table[x,6] = length(which(species_table[x, grep("Vaginal", map_table$Env)] > 0))/length(grep("Vaginal", map_table$Env))
}
freq_table_tb = freq_table[order(freq_table$toothbrush, decreasing = T),]

###
### Alpha diversity by sample type
###

# compute shannon and simpson diversity metrics
diversity_vec = matrix(nrow = dim(species_table)[2], ncol = 2)
diversity_vec = as.data.frame(diversity_vec)
for (a in 1:dim(species_table)[2]) {
  diversity_vec[a,1] = diversity(species_table[,a], index = "shannon")
  diversity_vec[a,2] = diversity(species_table[,a], index = "simpson")
}
colnames(diversity_vec) = c("Shannon", "Simpson")

# add environment factor
diversity_vec$Env = map_table$Env
diversity_vec$Env = factor(diversity_vec$Env,
                           levels = c("toothbrush", "Oral", "Skin", "Gut", "Vaginal", "building_dust"))

## boxplots (FIGURE 1B)
ggplot(diversity_vec, aes(x = Env, y = Shannon, fill = Env)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  #geom_violin(alpha = 0.8) +
  geom_jitter(size = 0.2, width = 0.1, alpha = 0.35) +
  scale_fill_manual(values=c("blue", "lightgreen", "yellow", "brown", "violet", "gray"),
                    labels = c("Toothbrush", "HMP-Oral", "HMP-Skin", "HMP-Gut", "HMP-Vaginal", "Building dust")) +
  theme_bw() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 14))

## ANOVA stats
summary(aov(Shannon ~ Env, diversity_vec))
TukeyHSD(aov(Shannon ~ Env, diversity_vec))

###
### Beta diversity between sample types
###

# beta-diversity measure
beta <- vegdist(t(species_table), 'jaccard', binary = T)
beta.mat <- as.matrix(beta)

# projection
pcoa <- cmdscale(beta, k = 4, eig = TRUE)

# cleanup
ord <- as.data.frame(pcoa$points)
names(ord) <- c('pcoa1', 'pcoa2', 'pcoa3', 'pcoa4') ## rename coordinates

# Percent explained variation
eig <- eigenvals(pcoa)
100*head(eig/sum(eig))

# add metadata
ord$Category = map_table$Env

# re-order for plot overlay
ord = rbind(ord[c(grep("Oral|Gut|Skin|Vaginal", ord$Category),
                  grep("dust", ord$Category),
                  grep("toothbrush", ord$Category)),])
ord = as.data.frame(ord)

# set vectors for point shape, size, and stroke
vec_size = rep(2.7, 2547)
vec_size[grep("toothbrush", ord$Category)] = rep(5, 34)

vec_shape = rep(21, 2547)
vec_shape[grep("toothbrush", ord$Category)] = rep(23, 34)

vec_stroke = rep(0.5, 2547)
vec_stroke[grep("toothbrush", ord$Category)] = rep(0.7, 34)

#re-order legend
ord$Category = factor(ord$Category,
                      levels = c("toothbrush", "Oral", "Skin", "Gut", "Vaginal", "building_dust"))

## plot PCoA (FIGURE 1A)
ggplot(data = ord, aes(x = pcoa1, y = pcoa2, fill = Category)) +
  geom_point(alpha = 0.6, 
             stroke = vec_stroke,
             shape = vec_shape,
             size = vec_size) +
  theme_bw() +
  xlab('PCo1 (24.4%)') +
  ylab('PCo2 (10.2%)') +
  scale_fill_manual(values=c("blue", "lightgreen", "yellow", "brown", "violet", "gray"),
                    labels = c("Toothbrush", "HMP-Oral", "HMP-Skin", "HMP-Gut", "HMP-Vaginal", "Building dust")) +
  theme(axis.text = element_text(size=15, color = "black"),
        axis.title = element_text(size=16, color = "black"),
        legend.text = element_text(size=15, color = "black"),
        legend.title = element_blank()) +
  guides(fill = guide_legend(override.aes=list(shape=c(23, 21, 21, 21, 21, 21),
                                               size = 4))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

## PERMANOVA stats

# effect of sample type
beta <- vegdist(t(species_table), 'jaccard', binary = T)
adonis(beta ~ Env, data = map_table, permutations = 999)

# toothbrush vs. oral
beta <- vegdist(t(species_table[,grep("toothbrush|Oral",map_table$Env)]), 'jaccard', binary = T)
adonis(beta ~ Env, data = map_table[grep("toothbrush|Oral",map_table$Env),], permutations = 999)

# toothbrush vs. skin
beta <- vegdist(t(species_table[,grep("toothbrush|Skin",map_table$Env)]), 'jaccard', binary = T)
adonis(beta ~ Env, data = map_table[grep("toothbrush|Skin",map_table$Env),], permutations = 999)

# toothbrush vs. gut
beta <- vegdist(t(species_table[,grep("toothbrush|Gut",map_table$Env)]), 'jaccard', binary = T)
adonis(beta ~ Env, data = map_table[grep("toothbrush|Gut",map_table$Env),], permutations = 999)

# toothbrush vs. vaginal
beta <- vegdist(t(species_table[,grep("toothbrush|Vaginal",map_table$Env)]), 'jaccard', binary = T)
adonis(beta ~ Env, data = map_table[grep("toothbrush|Vaginal",map_table$Env),], permutations = 999)

# toothbrush vs. dust
beta <- vegdist(t(species_table[,grep("toothbrush|dust",map_table$Env)]), 'jaccard', binary = T)
adonis(beta ~ Env, data = map_table[grep("toothbrush|dust",map_table$Env),], permutations = 999)

## Distance between samples

# toothbrush
mean(as.matrix(vegdist(t(species_table[,grep("toothbrush",map_table$Env)]), 
                       'jaccard', binary = T)))
sd(as.matrix(vegdist(t(species_table[,grep("toothbrush",map_table$Env)]), 
                     'jaccard', binary = T)))/sqrt(34)

# oral
mean(as.matrix(vegdist(t(species_table[,grep("Oral",map_table$Env)]), 
                       'jaccard', binary = T)))
sd(as.matrix(vegdist(t(species_table[,grep("Oral",map_table$Env)]), 
                     'jaccard', binary = T)))/sqrt(length(grep("Oral",map_table$Env)))

###
### Toothbrush-associated taxa
###

### Frequency-abundance

# subest table
species_toothbrush = species_table[,grep("TMP", colnames(species_table))]

# add average
species_toothbrush = data.frame(species_toothbrush, avg = apply(species_toothbrush,1,mean))
species_toothbrush_avg_freq = data.frame(species_toothbrush, freq_table)
species_toothbrush_avg_freq = species_toothbrush_avg_freq[-c(which(species_toothbrush_avg_freq$avg == 0)),]
species_toothbrush_avg_freq = species_toothbrush_avg_freq[order(species_toothbrush_avg_freq$toothbrush, decreasing = T),]

## toothbrush microbiota frequency-abundance (FIGURE 1C)
ggplot(species_toothbrush_avg_freq, aes(x = avg, #size = avg, 
                                        y = toothbrush, fill = Oral)) +
  geom_point(alpha = 0.7, shape = 21, stroke = 0.5, size = 4) +
  theme_bw() +
  scale_fill_gradient(low = "white", high = "black") +
  scale_x_log10() +
  geom_hline(yintercept=0.5, lty = 2) +
  #scale_size_continuous(range=c(1.5, 7), breaks=c(0.1, 1, 5, 10, 15, 50)) +
  ylab('Toothbrush frequency') +
  xlab('Relative abundance (%)') +
  labs(fill = "HMP-Oral Freq.") +
  theme(axis.text.y = element_text(size=17, color = "black"), 
        axis.text.x = element_text(size=17, color = "black"),
        axis.title = element_text(size=18, color = "black"),
        legend.title = element_text(size=13, color = "black"),
        legend.text = element_text(size=13, color = "black"),
        legend.position = c(0.24,0.76))

### Heatmap of conserved taxa

# re-name by species
species_toothbrush_avg_freq$taxonomy = rownames(species_toothbrush_avg_freq)
rownames(species_toothbrush_avg_freq) = unlist(strsplit(rownames(species_toothbrush_avg_freq), "|", T))[grep("s__", unlist(strsplit(rownames(species_toothbrush_avg_freq), "|", T)))]
rownames(species_toothbrush_avg_freq) = substr(rownames(species_toothbrush_avg_freq), start = 4, 
                                               stop = nchar(rownames(species_toothbrush_avg_freq)))

# set color palette
mypalette <- colorRampPalette(c("white","black"))(n = 100)
mycolors = c(0:100)

# plot key for heat colors
hist(c(0:100), breaks = 100, col = mypalette, 
     border = mypalette, ylim = c(0,0.5),
     xlab = NULL, ylab = NULL, 
     axes = FALSE,
     main = NA)
axis(side = 1, labels = c("0.00", "0.25", "0.50", "0.75", "1.00"),
     at = c(0,25,50,75,100), las = 1, cex = 4)

# sidebars for phylum
phyla_list = unlist(strsplit(species_toothbrush_avg_freq$taxonomy, "|", T))[grep("p__", unlist(strsplit(species_toothbrush_avg_freq$taxonomy, "|", T)))]
phyla_list = substr(phyla_list, start = 4, stop = nchar(phyla_list))
phyla_list = phyla_list[1:length(which(species_toothbrush_avg_freq$toothbrush >= 0.5))]
phyla_col = data.frame(phyla_list, col = rep("white", length(phyla_list)))
phyla_col$col = as.character(phyla_col$col)
phyla_col$col[grep("Actinobacteria", phyla_col$phyla_list)] = rep("purple", length(phyla_col$col[grep("Actinobacteria", phyla_col$phyla_list)]))
phyla_col$col[grep("Bacteroidetes", phyla_col$phyla_list)] = rep("green", length(phyla_col$col[grep("Bacteroidetes", phyla_col$phyla_list)]))
phyla_col$col[grep("Firmicutes", phyla_col$phyla_list)] = rep("skyblue", length(phyla_col$col[grep("Firmicutes", phyla_col$phyla_list)]))
phyla_col$col[grep("Fusobacteria", phyla_col$phyla_list)] = rep("red", length(phyla_col$col[grep("Fusobacteria", phyla_col$phyla_list)]))
phyla_col$col[grep("Proteobacteria", phyla_col$phyla_list)] = rep("gold", length(phyla_col$col[grep("Proteobacteria", phyla_col$phyla_list)]))

# plot phyla legend
plot(1:50, col = "white")
legend(1,49, legend=c("Actinobacteria", "Bacteroidetes", "Firmicutes", "Fusobacteria", "Proteobacteria"),
       c("purple", "green", "skyblue", "red", "gold"), cex=1)

# "clean" species names
heat_data = as.matrix(species_toothbrush_avg_freq[1:length(which(species_toothbrush_avg_freq$toothbrush >= 0.5)),
                                                  c(36,38,39,40,41,37)])
rownames(heat_data) = c("Rothia mucilaginosa", "Streptococcus sanguinis", "Veillonella (unclassified)", "Streptococcus parasanguinis", "Rothia dentocariosa", 
                        "Veillonella parvula", "Streptococcus mitis", "Klebsiella oxytoca", "Stenotrophomonas maltophilia", "Streptococcus salivarius",
                        "Actinomyces viscosus", "Actinomyces oris", "Rothia aeria", "Enterobacter cloacae", "Streptococcus gordonii",           
                        "Haemophilus parainfluenzae", "Actinomyces odontolyticus", "Granulicatella unclassified", "Streptococcus australis", "Granulicatella adiacens", 
                        "Streptococcus infantis", "Kocuria rhizophila", "Abiotrophia defectiva", "Veillonella atypica", "Lautropia mirabilis",
                        "Klebsiella (unclassified)", "Pseudomonas (unclassified)", "Actinomyces massiliensis", "Oribacterium sinus", "Veillonella dispar",                 
                        "Neisseria (unclassified)", "Corynebacterium matruchotii", "Prevotella melaninogenica", "Gemella haemolysans", "Gemella (unclassified)", 
                        "Stomatobaculum longum", "Leptotrichia (unclassified)")

## plot heatmap (FIGURE 1D)
heatmap.2(heat_data,
          dendrogram = "row",
          Colv = FALSE,
          distfun = function(x) vegdist(x, method = 'jaccard'),
          key=FALSE, 
          symkey=FALSE, 
          density.info='none',
          trace = "none",
          cexCol = 1.6,
          labCol = c("Toothbrush", "HMP-Oral", "HMP-Skin", "HMP-Gut", "HMP-Vaginal", "Building dust"),
          labRow = as.expression(lapply(rownames(heat_data), function(a) bquote(italic(.(a))))),
          srtCol = 45,
          cexRow = 1.2,
          col = mypalette,
          #breaks = mycolors,
          lhei = c(0.5,20),
          lwid = c(0.8,5),
          offsetRow = 0,
          offsetCol = 0,
          RowSideColors = as.character(phyla_col$col),
          margins = c(8, 20))

# fraction of conserved toothbrush species (>50% samples) that are also conserved in oral microbiome
length(which(as.data.frame(heat_data)$Oral > 0.5))/length(as.data.frame(heat_data)$Oral)

###
### Sourcetracker ###
###

# quick stats
apply(100*st_predictions_genus, 2, mean)
apply(100*st_predictions_genus, 2, sd)/sqrt(34)

# proportions
st_predictions_genus_melt = melt(st_predictions_genus)
st_predictions_genus_melt$variable = factor(st_predictions_genus_melt$variable,
                                            levels = c("Gut", "Oral", "Skin", "Vaginal", "Unknown"))

## violin plot(FIGURE S2)
ggplot(st_predictions_genus_melt, 
       aes(x = variable, y = value, fill = variable)) +
  geom_violin(scale="width") +
  #geom_boxplot(outlier.shape = NA, width = 0.1, alpha = 0.6) +
  geom_jitter(width = 0.12, size = 1.5, alpha = 0.5) +
  stat_summary(fun.y = mean, geom = 'point', fill = 'white', shape = 21, size = 3.5) +
  theme_bw() +
  ylim(-0.001,1.05) +
  ylab("Fraction of toothbrush microbiota") +
  xlab("Putative source") +
  scale_fill_manual(values=c("brown", "lightgreen", "yellow", "violet", "orange"),
                    labels = c("HMP-Gut", "HMP-Oral", "HMP-Skin", "HMP-Vaginal", "Unknown")) +
  theme(axis.text.y = element_text(size=11, color = "black"), 
        axis.title = element_text(size=12, color = "black"),
        axis.text.x = element_text(size=12, color = "black", angle = 45, hjust = 1)) +
  theme(panel.border = element_rect(colour = "black")) +
  scale_x_discrete(limits = c("Gut", "Oral", "Skin", "Vaginal", "Unknown"), 
                   labels = c("HMP-Gut", "HMP-Oral", "HMP-Skin", "HMP-Vaginal", "Unknown")) +
  theme(legend.position = 'none')


###########################
##### HMP-ORAL SUBSET #####
###########################

# subset oral metagenomes
HMP_oral = t(species_table[,grep("Oral", map_table$Env)])

# beta-diversity measure
beta <- vegdist(HMP_oral, 'jaccard', binary = T)
beta.mat <- as.matrix(beta)

# projection
pcoa <- cmdscale(beta, k = 4, eig = TRUE)

# cleanup
ord <- as.data.frame(pcoa$points)
names(ord) <- c('pcoa1', 'pcoa2', 'pcoa3', 'pcoa4') ## rename coordinates

# add clusters (34 centroids)
ord$cluster <- kmeans(ord[, c('pcoa1', 'pcoa2')], centers = 34)$cluster
ord$cluster <- as.character(ord$cluster)

# get centroid closest neighbor
centroids = kmeans(ord[, c('pcoa1', 'pcoa2')], centers = 34)$centers
ord_clust = rbind(ord[,c(1,2,5)],
                  data.frame(centroids,
                             cluster = rep("centroid", 34)))
ord_cent_dist = as.matrix(dist(ord_clust[,c(1,2)]))
ord_cent_dist = as.data.frame(ord_cent_dist)
min_d = apply(ord_cent_dist[c(1:1259), c((1293-33):1293)],
              2, min)
vec_centroid = c(1:34)
for (i in 1:34) {
  vec_centroid[i] = which(ord_cent_dist[,1259 + i] == min_d[i])
}
ord$cent_neighbor = rep("point", length(dim(ord)[1]))
ord$cent_neighbor[vec_centroid] = rep("centroid", length(vec_centroid))

# Percent explained variation
eig <- eigenvals(pcoa)
100*head(eig/sum(eig))

# plot PCoA (FIGURE S3)
ggplot(data = ord, aes(x = pcoa1, y = pcoa2, fill = cent_neighbor)) +
  geom_point(size = 3, alpha = 0.7, pch = 21) +
  theme_bw() +
  xlab('PCo1 (19.1%)') +
  ylab('PCo2 (13.4%)') +
  scale_fill_manual(values = c("red", "gray")) +
  theme(axis.text = element_text(size=15, color = "black"),
        axis.title = element_text(size=16, color = "black"),
        legend.position = "none")

# SRA accession list (TABLE S3)
data.frame(number = c(1:34),
           name = rownames(ord[vec_centroid,]))


######################
##### RESISTOMES #####
######################

# trim table to contian only ARGs detected in the toothbrush and oral samples
sb_u90_clean = sb_u90[which(apply(sb_u90[,1:68],1,sum) > 0),]

# add summary stats
sb_u90_clean$tb_avg = apply(sb_u90_clean[,1:34],1,mean)
sb_u90_clean$hmp_avg = apply(sb_u90_clean[,35:68],1,mean)
sb_u90_clean$tb_se = apply(sb_u90_clean[,1:34],1,sd)/sqrt(34)
sb_u90_clean$hmp_se = apply(sb_u90_clean[,35:68],1,sd)/sqrt(34)

# compute ARG frequency in respective sample types
freq_table = matrix(nrow = dim(sb_u90_clean)[1], ncol = 3)
freq_table = as.data.frame(freq_table)
for (i in 1:dim(sb_u90_clean)[1]) {
  freq_table[i,1] = length(which(sb_u90_clean[i,1:68] > 0))/68
  freq_table[i,2] = length(which(sb_u90_clean[i,1:34] > 0))/34
  freq_table[i,3] = length(which(sb_u90_clean[i,35:68] > 0))/34
}
colnames(freq_table) = c("freq_all", "freq_tmp", "freq_hmp")

# add frequencies and subset data by marker length
sb_u90_clean_freq = data.frame(sb_u90_clean, freq_table)
sb_u90_sub = sb_u90_clean_freq
sb_u90_sub = sb_u90_clean_freq[c(which(sb_u90_clean_freq$marker_length >= 30)),]

# write-table to export (edit for TABLE S2)
#write.table(sb_u90_sub, 'data/biobakery/shortbred/sb_u90_30aa.txt', sep="\t", col.names=NA, quote = FALSE)

# re-order based on frequency
sb_u90_sub = sb_u90_sub[order(sb_u90_sub$freq_hmp, decreasing = T),]
sb_u90_sub = sb_u90_sub[order(sb_u90_sub$freq_tmp, decreasing = T),]
sb_u90_sub = sb_u90_sub[order(sb_u90_sub$freq_all, decreasing = T),]

###
### Resistome enrichments
###

## RPKM Counts (t-test, wilcox-test)
RPKM_count = data.frame(t(sb_u90_sub[,grep("count", colnames(sb_u90_sub))]),
                        study = c(rep("Toothbrush", 34), rep("HMP", 34)))
# test for significance
p_vals = matrix(nrow=176, ncol=2)
p_vals = as.data.frame(p_vals)
for(i in 1:176) {
  p_vals[i,1] = t.test(log(RPKM_count[grep("Toothbrush", RPKM_count$study),i]+1),
                       log(RPKM_count[grep("HMP", RPKM_count$study),i]+1))$p.val
  p_vals[i,2] = wilcox.test(log(RPKM_count[grep("Toothbrush", RPKM_count$study),i]+1),
                            log(RPKM_count[grep("HMP", RPKM_count$study),i]+1))$p.val
}
colnames(p_vals) = c("p_ttest", "p_wilcox")
rownames(p_vals) = colnames(RPKM_count)[1:176] 
# result
length(which(p.adjust(p_vals[,1], method = "bonferroni") < 0.05))
rownames(p_vals)[which(p.adjust(p_vals[,1], method = "bonferroni") < 0.05)]

## RPKM Counts (remove 0s -- wilcox-test)
# remove 0's
RPKM_count_no0 = RPKM_count
for (i in 1:176) {
  RPKM_count_no0[,i][which(RPKM_count_no0[,i] == 0)] = rep(NA, length(which(RPKM_count_no0[,i] == 0)))
}
# test for significant difference
med_p_val = c(1:176)
for(i in 1:176) {
  med_p_val[i] = wilcox.test(log(RPKM_count_no0[grep("Toothbrush", RPKM_count$study),i]+1),
                             log(RPKM_count_no0[grep("HMP", RPKM_count$study),i]+1))$p.val
}
names(med_p_val) = colnames(RPKM_count)[1:176] 
# result
length(which(p.adjust(med_p_val, method = "bonferroni") < 0.05))

## RPKM Counts (binomial -- GLM)
# set binomial table
RPKM_count_binom = RPKM_count
for (i in 1:176) {
  RPKM_count_binom[,i][which(RPKM_count[,i] > 0)] = rep(1, length(which(RPKM_count[,i] > 0)))
}
# test for significance
p.glm = c(1:176)
for(i in 1:176) {
  p.glm[i] = unlist(anova(glm(as.numeric(RPKM_count_binom[,i]) ~ RPKM_count_binom$study,
                              family = binomial(link='logit')), test = "Chisq"))[10]
}
names(p.glm) = colnames(RPKM_count)[1:176] 
# result
length(which(p.adjust(p.glm, method = "bonferroni") < 0.05))
p.adjust(p.glm, method = "bonferroni")[which(p.adjust(p.glm, method = "bonferroni") < 0.05)]

# stats for TABLE S2
data.frame(GLM_p = p.glm, Bon_q = p.adjust(p.glm, method = "bonferroni"))
#write.table(data.frame(GLM_p = p.glm, Bon_q = p.adjust(p.glm, method = "bonferroni")), 
#            'data/biobakery/shortbred/enrichment_stats.txt',sep="\t", col.names=NA, quote = FALSE)

###
### Resistome summary stats
###

# total ARGs
dim(RPKM_count_binom)
head(RPKM_count_binom)

# toothbrush ARGs per sample
apply(RPKM_count_binom[c(1:34),-c(177)], 1, sum)
median(apply(RPKM_count_binom[c(1:34),-c(177)], 1, sum))
mean(apply(RPKM_count_binom[c(1:34),-c(177)], 1, sum))
sd(apply(RPKM_count_binom[c(1:34),-c(177)], 1, sum))/sqrt(34)

# toothbrush ARGs detected
length(which(apply(RPKM_count_binom[c(1:34),-c(177)], 2, sum) > 0))

# HMP ARGs per sample
apply(RPKM_count_binom[c(35:68),-c(177)], 1, sum)
median(apply(RPKM_count_binom[c(35:68),-c(177)], 1, sum))
mean(apply(RPKM_count_binom[c(35:68),-c(177)], 1, sum))
sd(apply(RPKM_count_binom[c(35:68),-c(177)], 1, sum))/sqrt(34)

# HMP ARGs detected
length(which(apply(RPKM_count_binom[c(35:68),-c(177)], 2, sum) > 0))

# t-test and wilcox test
t.test(apply(RPKM_count_binom[c(1:34),-c(177)], 1, sum),
       apply(RPKM_count_binom[c(35:68),-c(177)], 1, sum))
wilcox.test(apply(RPKM_count_binom[c(1:34),-c(177)], 1, sum),
            apply(RPKM_count_binom[c(35:68),-c(177)], 1, sum))

###
### Alpha diversity: resistomes and HMP oral subset (taxa)
###

## alpha diversity: HMP oral subset (taxa)

# shannon and simpson diversity metrics
diversity_hmp_sub = matrix(nrow = 34, ncol = 2)
diversity_hmp_sub = as.data.frame(diversity_hmp_sub)
for (a in 1:34) {
  diversity_hmp_sub[a,1] = diversity(species_hmp_subset[,a], index = "shannon")
  diversity_hmp_sub[a,2] = diversity(species_hmp_subset[,a], index = "simpson")
}
colnames(diversity_hmp_sub) = c("Shannon", "Simpson")
rownames(diversity_hmp_sub) = colnames(species_hmp_subset)

# stats
median(diversity_hmp_sub[,1])
mean(diversity_hmp_sub[,1])

# compare to tb diversity
t.test(diversity_hmp_sub[,1], diversity_vec[grep("toothbrush", diversity_vec$Env),1])
wilcox.test(diversity_hmp_sub[,1], diversity_vec[grep("toothbrush", diversity_vec$Env),1])

## alpha diversity: resistomes

# shannon and simpson diversity metrics
diversity_vec_res = matrix(nrow = 68, ncol = 2)
diversity_vec_res = as.data.frame(diversity_vec_res)
for (a in 1:68) {
  diversity_vec_res[a,1] = diversity(RPKM_count[a,c(1:176)], index = "shannon")
  diversity_vec_res[a,2] = diversity(RPKM_count[a,c(1:176)], index = "simpson")
}
colnames(diversity_vec_res) = c("Shannon", "Simpson")
rownames(diversity_vec_res) = rownames(RPKM_count)
diversity_vec_res$Env = RPKM_count[,177]

# stats
t.test(diversity_vec_res$Shannon[1:34], diversity_vec_res$Shannon[35:68])
wilcox.test(diversity_vec_res$Shannon[1:34], diversity_vec_res$Shannon[35:68])
t.test(diversity_vec_res$Simpson[1:34], diversity_vec_res$Simpson[35:68])
wilcox.test(diversity_vec_res$Simpson[1:34], diversity_vec_res$Simpson[35:68])

## alpha diversity: toothbrush taxa vs. toothbrush ARGs

# order ARG 
rownames(RPKM_count)[1:34]
rownames(RPKM_count)[1:34][mixedorder(rownames(RPKM_count)[1:34])]
diversity_vec_res$Shannon[1:34]
diversity_vec_res$Shannon[1:34][mixedorder(rownames(RPKM_count)[1:34])]

# order taxa
colnames(species_table)[1:34]
diversity_vec[1:34,1]

# correlation (spearman)
cor.test(diversity_vec_res$Shannon[1:34][mixedorder(rownames(RPKM_count)[1:34])], diversity_vec[1:34,1], method = "spearman")

###
### Resistome plots
###

### alpha diversity boxplot

# set data frame
ad_taxa_res = data.frame(Shannon = c(diversity_vec[grep("toothbrush", diversity_vec$Env),1],
                                     diversity_hmp_sub[,1], diversity_vec_res$Shannon),
                         Sample_type = c(rep("Toothbrush", 34), rep("HMP-Oral", 34),
                                         rep("Toothbrush", 34), rep("HMP-Oral", 34)),
                         Data_type = c(rep("Species", 68), rep("ARGs", 68)))

# set factors
ad_taxa_res$Data_type = factor(ad_taxa_res$Data_type,
                               levels = c("Species", "ARGs"))

## plot (FIGURE 2B) 
ggplot(ad_taxa_res, aes(x = Data_type, y = Shannon, fill = Sample_type)) +
  geom_boxplot(position = position_dodge(0.8), alpha = 0.6, outlier.shape = NA) +
  geom_jitter(alpha = 0.3, size = 0.4,
              position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.1)) +
  scale_fill_manual(values = c("green", "blue")) +
  theme_bw() +
  ylim(0, 4.25) +
  theme(legend.title = element_blank(),
        #legend.position = "bottom",
        legend.position = c(0.18, 0.12), 
        legend.text = element_text(size = 8),
        axis.title.x = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14))

### frequency and RPKMS

# subset by detection in >5% samples
sb_u90_sub = sb_u90_sub[c(which(sb_u90_sub$freq_all > 0.05)),]

# re-order based on frequency
sb_u90_sub = sb_u90_sub[order(sb_u90_sub$freq_hmp, decreasing = T),]
sb_u90_sub = sb_u90_sub[order(sb_u90_sub$freq_tmp, decreasing = T),]
sb_u90_sub = sb_u90_sub[order(sb_u90_sub$freq_all, decreasing = T),]

# write-table to export (add drug class from CARD database)
#write.table(sb_u90_sub, 'data/biobakery/shortbred/sb_card_u90_30aa_0.05freq.txt', sep="\t", col.names=NA, quote = FALSE)

# load table with gene category manually added
sb_u90_sub = read.table('data/biobakery/shortbred/sb_card_u90_30aa_0.05freq_edit.txt', sep='\t', h=T, row.names=1, check=F, comment='',quote="\"")

# re-order based on frequency
sb_u90_sub = sb_u90_sub[order(sb_u90_sub$freq_hmp, decreasing = F),]
sb_u90_sub = sb_u90_sub[order(sb_u90_sub$freq_tmp, decreasing = F),]
sb_u90_sub = sb_u90_sub[order(sb_u90_sub$freq_all, decreasing = F),]

# add layer for colors by drug class
sb_u90_sub$cols = rep("gray", length(sb_u90_sub$drug_class))
sb_u90_sub$cols[grep("amino", sb_u90_sub$drug_class)] = rep("darkblue", length(grep("amino", sb_u90_sub$drug_class)))
sb_u90_sub$cols[grep("beta", sb_u90_sub$drug_class)] = rep("red", length(grep("beta", sb_u90_sub$drug_class)))
sb_u90_sub$cols[grep("fluor", sb_u90_sub$drug_class)] = rep("darkgreen", length(grep("fluor", sb_u90_sub$drug_class)))
sb_u90_sub$cols[grep("fosf", sb_u90_sub$drug_class)] = rep("purple", length(grep("fosf", sb_u90_sub$drug_class)))
sb_u90_sub$cols[grep("macro", sb_u90_sub$drug_class)] = rep("darkorange", length(grep("macro", sb_u90_sub$drug_class)))
sb_u90_sub$cols[grep("multi", sb_u90_sub$drug_class)] = rep("gold", length(grep("multi", sb_u90_sub$drug_class)))
sb_u90_sub$cols[grep("peptide", sb_u90_sub$drug_class)] = rep("brown", length(grep("peptide", sb_u90_sub$drug_class)))
sb_u90_sub$cols[grep("tet", sb_u90_sub$drug_class)] = rep("pink", length(grep("tet", sb_u90_sub$drug_class)))
sb_u90_sub$cols[grep("tric", sb_u90_sub$drug_class)] = rep("skyblue", length(grep("tric", sb_u90_sub$drug_class)))
sb_u90_sub$cols[grep("sulf", sb_u90_sub$drug_class)] = rep("magenta", length(grep("sulf", sb_u90_sub$drug_class)))

# use only > 10%
sb_u90_sub = sb_u90_sub[c(which(sb_u90_sub$freq_all > 0.1)),]

# prep data
sb_u90_stats = data.frame(name = rep(rownames(sb_u90_sub), 2),
                          gene = rep(sb_u90_sub$gene, 2),
                          drug_class = rep(sb_u90_sub$drug_class, 2),
                          freq = c(sb_u90_sub$freq_tmp, -sb_u90_sub$freq_hmp),
                          study = c(rep("Toothbrush", length(sb_u90_sub$tb_avg)),
                                    rep("HMP-Oral", length(sb_u90_sub$hmp_avg))))
sb_u90_stats$name = factor(sb_u90_stats$name, levels = rownames(sb_u90_sub))

## frequency plot (FIGURE 3)
ggplot(sb_u90_stats,
       aes(x = name, y = freq, fill = study)) +
  geom_bar(stat = "identity", color = "darkgray") +
  geom_point(aes(x = name, y = 0),
             shape = 21, size = 2.2, color = "black",
             fill = rep(sb_u90_sub$cols, 2)) +
  xlab("ARG Protein Family") +
  ylab("Sample Frequency") +
  theme_bw() +
  scale_x_discrete(labels = as.character(sb_u90_sub$gene)) +
  ylim(-1.05,1.05) +
  geom_hline(yintercept = 0.5, lty = 2, col = "black") +
  geom_hline(yintercept = -0.5, lty = 2, col = "black") +
  scale_fill_manual(values = c("green", "blue")) +
  theme(legend.position = c(0.85,0.85),
        legend.text = element_text(size = 14),
        legend.title = element_blank()) +
  theme(axis.title = element_text(size = 14),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 6,
                                   face = "italic"),
        legend.position = "none") + 
  coord_flip() +
  theme(axis.ticks.y = element_line(size = 1.5, lineend = 'round',
                                    color = sb_u90_sub$cols))

### RPKM dot-plot

# prep data
sb_u90_dot = data.frame(t(sb_u90_sub[grep("count", colnames(sb_u90_sub))]),
                        study = c(rep("Toothbrush", 34), rep("HMP-Oral", 34)))
sb_u90_dot = melt(sb_u90_dot, "study")

## plot medians (FIGURE 3)
ggplot(sb_u90_dot,
       aes(x = variable, y = log(value+1), fill = study)) +
  geom_jitter(shape = 21, width = 0, size = 1, alpha = 0.2) +
  #geom_point(aes(x = variable, y = 0), 
  #          fill = "black", size = 0.8, shape = 21) +
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
               geom = "crossbar", width = 0.8, fatten = 2,
               aes(col = study)) +
  #stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
  #             geom = "crossbar", width = 0.7, fatten = 0.01, col = "darkgray") +
  theme_bw() +
  ylim(0.0001, 7.5) +
  ylab("log(RPKM + 1)") +
  scale_x_discrete(labels = as.character(sb_u90_sub$gene)) +
  scale_fill_manual(values = c("green", "blue")) +
  scale_color_manual(values = c("green", "blue")) +
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 14),
        #axis.text.y = element_text(size = 7),
        axis.text.y = element_blank(),
        legend.position = "none") + 
  coord_flip() +
  theme(axis.ticks.y = element_line(size = 1.5, lineend = 'round',
                                    color = sb_u90_sub$cols))

## plot (for legend only)
sb_u90_stats$drug_class = factor(sb_u90_stats$drug_class,
                                 levels = levels(sb_u90_stats$drug_class)[c(1:5,8:11,6,7)])
ggplot(sb_u90_stats,
       aes(x = name, y = freq, col = study)) +
  geom_bar(stat = "identity", color = "darkgray") +
  geom_point(aes(x = name, y = 0, fill = drug_class),
             shape = 21, size = 4, color = "black") +
  xlab("ARG Protein Family") +
  ylab("Sample Frequency") +
  theme_bw() +
  scale_x_discrete(labels = as.character(sb_u90_sub$gene)) +
  ylim(-1.05,1.05) +
  geom_hline(yintercept = 0.5, lty = 2, col = "black") +
  geom_hline(yintercept = -0.5, lty = 2, col = "black") +
  scale_fill_manual(values = c("red", "darkgreen", "purple", "darkorange", "brown", "pink", "skyblue", "gold", "gray")) +
  theme(axis.title = element_text(size = 14),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 6)) + 
  theme(legend.text = element_text(size = 14),
        legend.title = element_blank()) +
  coord_flip() +
  theme(axis.ticks.y = element_line(size = 1.5, lineend = 'round',
                                    color = sb_u90_sub$cols))

### core/conserved microbiota and resistome

## microbial taxa
tb_core = rownames(species_toothbrush_avg_freq[which(species_toothbrush_avg_freq$toothbrush >= .5),])
dust_core = rownames(species_toothbrush_avg_freq[which(species_toothbrush_avg_freq$building_dust >= .5),])
oral_core = rownames(species_toothbrush_avg_freq[which(species_toothbrush_avg_freq$Oral >= .5),])
gut_core = rownames(species_toothbrush_avg_freq[which(species_toothbrush_avg_freq$Gut >= .5),])
skin_core = rownames(species_toothbrush_avg_freq[which(species_toothbrush_avg_freq$Skin >= .5),])
vaginal_core = rownames(species_toothbrush_avg_freq[which(species_toothbrush_avg_freq$Vaginal >= .5),])

# print (FIGURE 2A)
venn.diagram(list(tb_core, oral_core), 
             category.names = c("Toothbrush", "Oral"),
             filename = 'data/core/conserved_taxa.png',
             fill = c("blue", "green"),
             cex = 4,
             #scaled = 0,
             cat.cex = 0,
             cat.pos = 0)

## resistomes
tb_core_res = sb_u90_sub$gene[which(sb_u90_sub$freq_tmp >= 0.5)]
oral_core_res = sb_u90_sub$gene[which(sb_u90_sub$freq_hmp >= 0.5)]

# print (FIGURE 2A)
venn.diagram(list(tb_core_res, oral_core_res), 
             category.names = c("Toothbrush", "Oral"),
             filename = 'data/core/conserved_ARG.png',
             fill = c("blue", "green"),
             cex = 4,
             #scaled = 0,
             cat.cex = 0,
             cat.pos = 0) 


####################
##### METADATA #####
####################

## plot theme
blank_theme <- theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_blank(),
        panel.grid=element_blank(),
        axis.ticks = element_blank(),
        plot.title=element_text(size=14, face="bold"))

## plot function
gg_survey = function (X, Z) {
  ggplot(X, 
         aes(x="", y=value, fill=Var1))+
    geom_bar(width = 0.5, stat = "identity") +
    coord_polar("y", start=0) +
    scale_fill_brewer(palette = "Blues") +
    blank_theme +
    theme(axis.text.x = element_blank(),
          legend.title = element_blank(),
          plot.title = element_text(size = 10, hjust = 0.5),
          legend.text = element_text(size = 6),
          legend.position = c(1, 0.5)) +
    geom_text(aes(y = value/3 + c(0, cumsum(value)[-length(value)]), 
                  label = percent(value/(34-Z))), size=2.5) +
    theme(legend.key.size = unit(3, "mm")) +
    theme(legend.text = element_blank())
}

## Survey questions
colnames(survey_table_clean)

# Gender
RC_gender = melt(tapply(rep(1, 34), survey_table_clean$Gender, sum))
RC_gender = as.data.frame(rbind(RC_gender[2,],
                                RC_gender[1,]))
S1 = gg_survey(RC_gender, 0) +
  ggtitle("Gender identification")

# Cavities
RC_cavities = melt(tapply(rep(1, 34), survey_table_clean$Cavities, sum))
RC_cavities = as.data.frame(rbind(RC_cavities[1,],
                                  RC_cavities[3,],
                                  RC_cavities[2,]))
RC_cavities$Var1 = factor(RC_cavities$Var1,
                          levels = c("0", "1-5", ">5"))
S2 = gg_survey(RC_cavities, 2) +
  ggtitle("Cavities")

# Oral_Infections_P_G
RC_OI = melt(tapply(rep(1, 34), survey_table_clean$Oral_Infections_P_G, sum))
RC_OI$Var1 = factor(RC_OI$Var1,
                    levels = rev(levels(RC_OI$Var1)))
S3 = gg_survey(RC_OI, 0) +
  ggtitle("Any oral infections")

# Current_dental_devices_NG_R
RC_DD = melt(tapply(rep(1, 34), survey_table_clean$Current_dental_devices_NG_R, sum))
RC_DD$Var1 = factor(RC_DD$Var1,
                    levels = rev(levels(RC_DD$Var1)))
length(which(is.na(survey_table_clean$Current_dental_devices_NG_R)))
S4 = gg_survey(RC_DD, 3) +
  ggtitle("Dental device")

# Nonantibiotic_oral_medications
RC_med = melt(tapply(rep(1, 34), survey_table_clean$Nonantibiotic_oral_medications, sum))
RC_med$Var1 = factor(RC_med$Var1,
                     levels = rev(levels(RC_med$Var1)))
S5 = gg_survey(RC_med, 0) +
  ggtitle("Nonantibiotic oral meds")

# Missing_teeth
RC_teeth = melt(tapply(rep(1, 34), survey_table_clean$Missing_teeth, sum))
RC_teeth$Var1 = factor(RC_teeth$Var1,
                       levels = rev(levels(RC_teeth$Var1)))
S6 = gg_survey(RC_teeth, 0) +
  ggtitle("Any missing teeth")

# Adenoids_or_tonsils_removed
RC_AT = melt(tapply(rep(1, 34), survey_table_clean$Adenoids_or_tonsils_removed, sum))
RC_AT$Var1 = factor(RC_AT$Var1,
                    levels = rev(levels(RC_AT$Var1)))
S7 = gg_survey(RC_AT, 0) +
  ggtitle("Adenoids/tonsils removed")

# Alcohol
RC_alc = melt(tapply(rep(1, 34), survey_table_clean$Alcohol_Consumption, sum))
RC_alc = as.data.frame(rbind(RC_alc[2,],
                             RC_alc[1,]))
S8 = gg_survey(RC_alc, 0) +
  ggtitle("Alcohol consumption")

# Coffee_Consumption
RC_coff = melt(tapply(rep(1, 34), survey_table_clean$Coffee_Consumption, sum))
RC_coff = as.data.frame(rbind(RC_coff[2,],
                              RC_coff[1,]))
S9 = gg_survey(RC_coff, 0) +
  ggtitle("Coffee consumption")

# Dietary_Restrictions
RC_diet = melt(tapply(rep(1, 34), survey_table_clean$Dietary_Restrictions, sum))
RC_diet$Var1 = factor(RC_diet$Var1,
                      levels = rev(levels(RC_diet$Var1)))
S10 = gg_survey(RC_diet, 0) +
  ggtitle("Dietary restrictions")

# Toothpaste_Brand
RC_tp = melt(tapply(rep(1, 34), survey_table_clean$Toothpaste_Brand, sum))
RC_tp = as.data.frame(rbind(RC_tp[3,],
                            RC_tp[2,],
                            RC_tp[1,]))
S11 = gg_survey(RC_tp, 0) +
  ggtitle("Toothpaste brand")

# Toothbrush age
RC_tb_age = melt(tapply(rep(1, 34), survey_table_clean$Current_Toothbrush_Duration, sum))
RC_tb_age = as.data.frame(rbind(RC_tb_age[2,],
                                RC_tb_age[4,],
                                RC_tb_age[3,], 
                                RC_tb_age[1,]))
RC_tb_age$Var1 = factor(RC_tb_age$Var1,
                        levels = c("<1 mo.", "1-2 mo.", "2-3 mo.", ">3 mo."))
S12 = gg_survey(RC_tb_age, 0) +
  ggtitle("Toothbrush age")

# Brushing_Teeth_Frequency
RC_tb_freq = melt(tapply(rep(1, 34), survey_table_clean$Brushing_Teeth_Frequency, sum))
RC_tb_freq = as.data.frame(rbind(RC_tb_freq[1,],
                                 RC_tb_freq[3,],
                                 RC_tb_freq[2,]))
RC_tb_freq$Var1 = factor(RC_tb_freq$Var1,
                         levels = c("1x", "2x", ">2x"))
S13 = gg_survey(RC_tb_freq, 0) +
  ggtitle("Toothbrushing per day")

# Toothbrush_Location
RC_tb_loc = melt(tapply(rep(1, 34), survey_table_clean$Toothbrush_Location, sum))
RC_tb_loc$Var1 = factor(RC_tb_loc$Var1,
                        levels = rev(levels(RC_tb_loc$Var1)))
length(which(is.na(survey_table_clean$Toothbrush_Location)))
S14 = gg_survey(RC_tb_loc, 4) +
  ggtitle("Toothbrush storage")

# Nearby toothbrush
RC_tb_near = melt(tapply(rep(1, 34), survey_table_clean$Nearby_Toothbrush, sum))
RC_tb_near$Var1 = factor(RC_tb_near$Var1,
                         levels = rev(levels(RC_tb_near$Var1)))
S15 = gg_survey(RC_tb_near, 0) +
  ggtitle("Nearby toothbrush")

# Floss_Frequency
RC_floss = melt(tapply(rep(1, 34), survey_table_clean$Floss_Frequency, sum))
RC_floss = as.data.frame(rbind(RC_floss[3,],
                               RC_floss[1,],
                               RC_floss[2,]))
RC_floss$Var1 = factor(RC_floss$Var1,
                       levels = c("Daily", "1-3x per wk", "Less freq."))
S16 = gg_survey(RC_floss, 0) +
  ggtitle("Floss use")

# Mouthwash_Frequency
RC_mw = melt(tapply(rep(1, 34), survey_table_clean$Mouthwash_Frequency, sum))
RC_mw = as.data.frame(rbind(RC_mw[3,],
                            RC_mw[1,],
                            RC_mw[2,]))
RC_mw$Var1 = factor(RC_mw$Var1,
                    levels = c("Daily", "1-3x per wk", "Less freq."))
S17 = gg_survey(RC_mw, 0) +
  ggtitle("Mouthwash use")

# Bathroom_Window
RC_window = melt(tapply(rep(1, 34), survey_table_clean$Bathroom_Window, sum))
RC_window$Var1 = factor(RC_window$Var1,
                        levels = rev(levels(RC_window$Var1)))
S18 = gg_survey(RC_window, 1) +
  ggtitle("Bathroom window")

# Window_Open_Frequency
RC_window_freq = melt(tapply(rep(1, 34), survey_table_clean$Window_Open_Frequency, sum))
RC_window_freq = as.data.frame(rbind(RC_window_freq[3,],
                                     RC_window_freq[2,],
                                     RC_window_freq[4,],
                                     RC_window_freq[1,]))
RC_window_freq$Var1 = factor(RC_window_freq$Var1,
                             levels = c("Daily", "Weekly", "Monthly", "Never"))
S19 = gg_survey(RC_window_freq, 1) +
  ggtitle("Window open frequency")

# Pets
RC_pets = melt(tapply(rep(1, 34), survey_table_clean$Pets, sum))
RC_pets$Var1 = factor(RC_pets$Var1,
                      levels = rev(levels(RC_pets$Var1)))
S20 = gg_survey(RC_pets, 0) +
  ggtitle("Domestic pets")

## plot all (edit for FIGURE 4)
grid.arrange(S1, S13, S12, S11, S16, S17, S14, 
             S15, S18, S19, S2, S4, S6, S3,S7, 
             S5, S9, S8, S10, S20, ncol = 4)

###
### metadata correlations with alpha diversity
###

## add diversity metrics to metadata table
survey_table_clean$AD_taxa = diversity_vec[1:34,1]
survey_table_clean$AD_ARGs = diversity_vec_res$Shannon[1:34][mixedorder(rownames(RPKM_count)[1:34])]

## stats

# taxa
colnames(survey_table_clean)[1:21]
div_meta = as.data.frame(matrix(nrow = 21, ncol = 2))
for (i in 1:21) {
  div_meta[i,1] = unlist(summary(aov(survey_table_clean$AD_taxa ~ survey_table_clean[,i])))[9]
  div_meta[i,2] = as.numeric(unlist(kruskal.test(survey_table_clean$AD_taxa ~ survey_table_clean[,i]))[3])
}
colnames(div_meta) = c("aov", "kw")
rownames(div_meta) = colnames(survey_table_clean)[1:21]

# ARGs
colnames(survey_table_clean)[1:21]
div_res_meta = as.data.frame(matrix(nrow = 21, ncol = 2))
for (i in 1:21) {
  div_res_meta[i,1] = unlist(summary(aov(survey_table_clean$AD_ARGs ~ survey_table_clean[,i])))[9]
  div_res_meta[i,2] = as.numeric(unlist(kruskal.test(survey_table_clean$AD_ARGs ~ survey_table_clean[,i]))[3])
}
colnames(div_res_meta) = c("aov", "kw")
rownames(div_res_meta) = colnames(survey_table_clean)[1:21]

## outputs
div_meta
div_res_meta

###
### factors that correlated with microbiome diversity
###

### taxa

# gender
taxa_gen = ggplot(survey_table_clean,
                  aes(x = Gender, y = AD_taxa)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.9, fill = "gray") +
  geom_jitter(width = 0.1, alpha = 0.5, size = 0.8) +
  theme_bw() +
  ylab("Shannon index") +
  ylim(0.3, 3.5) +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        #axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none")

# adenoids and tonsils
taxa_AT = ggplot(survey_table_clean,
                 aes(x = Adenoids_or_tonsils_removed, y = AD_taxa)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.9, fill = "gray") +
  geom_jitter(width = 0.1, alpha = 0.5, size = 0.8) +
  theme_bw() +
  ylab("Shannon index") +
  ylim(0.3, 3.5) +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        #axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none")

# mouthwash frequency
taxa_MF = ggplot(survey_table_clean,
                 aes(x = Mouthwash_Frequency, y = AD_taxa)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.9, fill = "gray") +
  geom_jitter(width = 0.1, alpha = 0.5, size = 0.8) +
  theme_bw() +
  ylab("Shannon index") +
  ylim(0.3, 3.5) +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        #axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none")

## taxa boxplots (FIGURE S4; add p-values)
grid.arrange(taxa_gen, taxa_MF, taxa_AT, nrow = 1)

### ARGs

# gender
ARG_gen = ggplot(survey_table_clean,
                 aes(x = Gender, y = AD_ARGs)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.9, fill = "gray") +
  geom_jitter(width = 0.1, alpha = 0.5, size = 0.8) +
  theme_bw() +
  ylab("Shannon index") +
  ylim(0.3, 4.25) +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        #axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none")

# toothbrushing frequency
ARG_tb_freq = ggplot(survey_table_clean,
                     aes(x = Brushing_Teeth_Frequency, y = AD_ARGs)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.9, fill = "gray") +
  geom_jitter(width = 0.1, alpha = 0.5, size = 0.8) +
  theme_bw() +
  ylab("Shannon index") +
  ylim(0.3, 4.25) +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        #axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none")

# nonantimicrobial medication
survey_table_clean$Nonantibiotic_oral_medications = factor(survey_table_clean$Nonantibiotic_oral_medications,
                                                           levels = c("Yes", "No"))
ARG_meds = ggplot(survey_table_clean,
                  aes(x = Nonantibiotic_oral_medications, y = AD_ARGs)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.9, fill = "gray") +
  geom_jitter(width = 0.1, alpha = 0.5, size = 0.8) +
  theme_bw() +
  ylab("Shannon index") +
  ylim(0.3, 4.25) +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        #axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none")

# toothbrush storage
df_tb_storage = survey_table_clean[-c(which(is.na(survey_table_clean$Toothbrush_Location))),]
ARG_tb_storage = ggplot(df_tb_storage,
                        aes(x = Toothbrush_Location, y = AD_ARGs)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.9, fill = "gray") +
  geom_jitter(width = 0.1, alpha = 0.5, size = 0.8) +
  theme_bw() +
  ylab("Shannon index") +
  ylim(0.3, 4.25) +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        #axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none")

# bathroom window
df_bw = survey_table_clean[-c(which(is.na(survey_table_clean$Bathroom_Window))),]
df_bw$Bathroom_Window = factor(df_bw$Bathroom_Window,
                               levels = c("Yes", "No"))
ARG_bw = ggplot(df_bw,
                aes(x = Bathroom_Window, y = AD_ARGs)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.9, fill = "gray") +
  geom_jitter(width = 0.1, alpha = 0.5, size = 0.8) +
  theme_bw() +
  ylab("Shannon index") +
  ylim(0.3, 4.25) +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        #axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none")

# bathroom window open frequency
df_bw_open = survey_table_clean[-c(which(is.na(survey_table_clean$Window_Open_Frequency))),]
df_bw_open$Window_Open_Frequency = factor(df_bw_open$Window_Open_Frequency,
                                          levels = c("Daily", "Weekly", "Monthly", "Never"))
ARG_bw_open = ggplot(df_bw_open,
                     aes(x = Window_Open_Frequency, y = AD_ARGs)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.9, fill = "gray") +
  geom_jitter(width = 0.1, alpha = 0.5, size = 0.8) +
  theme_bw() +
  ylab("Shannon index") +
  ylim(0.3, 4.25) +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        #axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none")

## ARG boxplots (FIGURE S6; add p-values)
grid.arrange(ARG_gen, ARG_tb_freq, ARG_meds, 
             ARG_tb_storage, ARG_bw, ARG_bw_open, nrow = 2)

###
### oral hygiene stats
###

### taxa

## toothbrushing per day

# all
summary(aov(AD_taxa ~ Brushing_Teeth_Frequency, survey_table_clean))
kruskal.test(AD_taxa ~ Brushing_Teeth_Frequency, survey_table_clean)

# pairwise
wilcox.test(survey_table_clean$AD_taxa[grep("1x", survey_table_clean$Brushing_Teeth_Frequency)],
            survey_table_clean$AD_taxa[grep("^2x$", survey_table_clean$Brushing_Teeth_Frequency)])
wilcox.test(survey_table_clean$AD_taxa[grep("1x", survey_table_clean$Brushing_Teeth_Frequency)],
            survey_table_clean$AD_taxa[grep(">2x", survey_table_clean$Brushing_Teeth_Frequency)])
wilcox.test(survey_table_clean$AD_taxa[grep("^2x$", survey_table_clean$Brushing_Teeth_Frequency)],
            survey_table_clean$AD_taxa[grep(">2x", survey_table_clean$Brushing_Teeth_Frequency)])

## floss frequency

# all
summary(aov(AD_taxa ~ Floss_Frequency, survey_table_clean))
kruskal.test(AD_taxa ~ Floss_Frequency, survey_table_clean)

# pairwise
wilcox.test(survey_table_clean$AD_taxa[grep("Daily", survey_table_clean$Floss_Frequency)],
            survey_table_clean$AD_taxa[grep("Less", survey_table_clean$Floss_Frequency)])
wilcox.test(survey_table_clean$AD_taxa[grep("Daily", survey_table_clean$Floss_Frequency)],
            survey_table_clean$AD_taxa[grep("1-3", survey_table_clean$Floss_Frequency)])
wilcox.test(survey_table_clean$AD_taxa[grep("Less", survey_table_clean$Floss_Frequency)],
            survey_table_clean$AD_taxa[grep("1-3", survey_table_clean$Floss_Frequency)])

## mouthwash frequency

# all
summary(aov(AD_taxa ~ Mouthwash_Frequency, survey_table_clean))
kruskal.test(AD_taxa ~ Mouthwash_Frequency, survey_table_clean)

# pairwise
wilcox.test(survey_table_clean$AD_taxa[grep("Daily", survey_table_clean$Mouthwash_Frequency)],
            survey_table_clean$AD_taxa[grep("Less", survey_table_clean$Mouthwash_Frequency)])
wilcox.test(survey_table_clean$AD_taxa[grep("Daily", survey_table_clean$Mouthwash_Frequency)],
            survey_table_clean$AD_taxa[grep("1-3", survey_table_clean$Mouthwash_Frequency)])
wilcox.test(survey_table_clean$AD_taxa[grep("Less", survey_table_clean$Mouthwash_Frequency)],
            survey_table_clean$AD_taxa[grep("1-3", survey_table_clean$Mouthwash_Frequency)])

### ARGs

## toothbrushing per day

# all
summary(aov(AD_ARGs ~ Brushing_Teeth_Frequency, survey_table_clean))
kruskal.test(AD_ARGs ~ Brushing_Teeth_Frequency, survey_table_clean)

# pairwise
wilcox.test(survey_table_clean$AD_ARGs[grep("1x", survey_table_clean$Brushing_Teeth_Frequency)],
            survey_table_clean$AD_ARGs[grep("^2x$", survey_table_clean$Brushing_Teeth_Frequency)])
wilcox.test(survey_table_clean$AD_ARGs[grep("1x", survey_table_clean$Brushing_Teeth_Frequency)],
            survey_table_clean$AD_ARGs[grep(">2x", survey_table_clean$Brushing_Teeth_Frequency)])
wilcox.test(survey_table_clean$AD_ARGs[grep("^2x$", survey_table_clean$Brushing_Teeth_Frequency)],
            survey_table_clean$AD_ARGs[grep(">2x", survey_table_clean$Brushing_Teeth_Frequency)])

## floss frequency

# all
summary(aov(AD_ARGs ~ Floss_Frequency, survey_table_clean))
kruskal.test(AD_ARGs ~ Floss_Frequency, survey_table_clean)

# pairwise
wilcox.test(survey_table_clean$AD_ARGs[grep("Daily", survey_table_clean$Floss_Frequency)],
            survey_table_clean$AD_ARGs[grep("Less", survey_table_clean$Floss_Frequency)])
wilcox.test(survey_table_clean$AD_ARGs[grep("Daily", survey_table_clean$Floss_Frequency)],
            survey_table_clean$AD_ARGs[grep("1-3", survey_table_clean$Floss_Frequency)])
wilcox.test(survey_table_clean$AD_ARGs[grep("Less", survey_table_clean$Floss_Frequency)],
            survey_table_clean$AD_ARGs[grep("1-3", survey_table_clean$Floss_Frequency)])

## mouthwash frequency

# all
summary(aov(AD_ARGs ~ Mouthwash_Frequency, survey_table_clean))
kruskal.test(AD_ARGs ~ Mouthwash_Frequency, survey_table_clean)

# pairwise
wilcox.test(survey_table_clean$AD_ARGs[grep("Daily", survey_table_clean$Mouthwash_Frequency)],
            survey_table_clean$AD_ARGs[grep("Less", survey_table_clean$Mouthwash_Frequency)])
wilcox.test(survey_table_clean$AD_ARGs[grep("Daily", survey_table_clean$Mouthwash_Frequency)],
            survey_table_clean$AD_ARGs[grep("1-3", survey_table_clean$Mouthwash_Frequency)])
wilcox.test(survey_table_clean$AD_ARGs[grep("Less", survey_table_clean$Mouthwash_Frequency)],
            survey_table_clean$AD_ARGs[grep("1-3", survey_table_clean$Mouthwash_Frequency)])

###
### oral hygiene plots
###

## set factors
survey_table_clean$Brushing_Teeth_Frequency = factor(survey_table_clean$Brushing_Teeth_Frequency,
                                                     levels = c(">2x", "2x", "1x"))
survey_table_clean$Floss_Frequency = factor(survey_table_clean$Floss_Frequency,
                                            levels = c("Daily", "1-3x per wk", "Less freq."))
survey_table_clean$Mouthwash_Frequency = factor(survey_table_clean$Mouthwash_Frequency,
                                                levels = c("Daily", "1-3x per wk", "Less freq."))

## plots

# taxa
tb_tax = ggplot(survey_table_clean,
                aes(x = Brushing_Teeth_Frequency, y = AD_taxa)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.9, fill = "gray") +
  geom_jitter(width = 0., alpha = 0.7, size = 0.6) +
  theme_bw() +
  xlab("Toothbrushing per day") +
  ylab("Shannon (Taxa)") +
  ylim(0.3, 3.5) +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        #axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none")

fl_tax = ggplot(survey_table_clean,
                aes(x = Floss_Frequency, y = AD_taxa)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.9, fill = "gray") +
  geom_jitter(width = 0., alpha = 0.7, size = 0.6) +
  theme_bw() +
  xlab("Toothbrushing per day") +
  ylab("Shannon (Taxa)") +
  ylim(0.3, 3.5) +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        #axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none")

mw_tax = ggplot(survey_table_clean,
                aes(x = Mouthwash_Frequency, y = AD_taxa)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.9, fill = "gray") +
  geom_jitter(width = 0., alpha = 0.7, size = 0.6) +
  theme_bw() +
  xlab("Toothbrushing per day") +
  ylab("Shannon (Taxa)") +
  ylim(0.3, 3.5) +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        #axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none")

# ARGs
tb_arg = ggplot(survey_table_clean,
                aes(x = Brushing_Teeth_Frequency, y = AD_ARGs)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.9, fill = "gray") +
  geom_jitter(width = 0., alpha = 0.7, size = 0.6) +
  theme_bw() +
  xlab("Toothbrushing per day") +
  ylab("Shannon (ARGs)") +
  ylim(0.3, 4.25) +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        #axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none")

fl_arg = ggplot(survey_table_clean,
                aes(x = Floss_Frequency, y = AD_ARGs)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.9, fill = "gray") +
  geom_jitter(width = 0., alpha = 0.7, size = 0.6) +
  theme_bw() +
  xlab("Toothbrushing per day") +
  ylab("Shannon (ARGs)") +
  ylim(0.3, 4.25) +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        #axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none")

mw_arg = ggplot(survey_table_clean,
                aes(x = Mouthwash_Frequency, y = AD_ARGs)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.9, fill = "gray") +
  geom_jitter(width = 0., alpha = 0.7, size = 0.6) +
  theme_bw() +
  xlab("Toothbrushing per day") +
  ylab("Shannon (ARGs)") +
  ylim(0.3, 4.25) +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        #axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none")

# plot (FIGURE 5; add significant p-values)
grid.arrange(tb_tax, fl_tax, mw_tax,
             tb_arg, fl_arg, mw_arg, nrow = 2)

###
### metadata correlations with particular (groups of) taxa or ARGs: halla
###

## set heat colors
mypalette <- colorRampPalette(c("white","violet"))(n = 99)
mycolors = seq(0, 0.5, length.out = 100)

# plot key for heat colors
hist(c(0:100), breaks = 100, col = mypalette, 
     border = mypalette, ylim = c(0,0.5),
     xlab = NULL, 
     ylab = NULL, 
     axes = FALSE,
     main = NA)
axis(side = 1, labels = c("0.0", "0.1", "0.2", "0.3", "0.4", "0.5"),
     at = c(0,20,40,60,80,100), las = 1, cex = 4)

## taxon associations heatmap (FIGURE 6; edit by adding 'association boxes')
heatmap.2(as.matrix(halla_species),
          Rowv = FALSE,
          Colv = FALSE,
          dendrogram="none",
          key=FALSE, 
          symkey=FALSE, 
          density.info='none',
          trace = "none",
          cexCol = 1.3,
          srtCol = 65,
          cexRow = 1.2,
          labRow = as.expression(lapply(rownames(halla_species), function(a) bquote(italic(.(a))))),
          labCol = c("Nonantibiotic oral meds", "Mouthwash frequency", "Dietary restrictions", "Alcohol consumption",
                     "Coffee consumption", "Toothbrush storage", "Bathroom window", "Window open frequency",
                     "Domestic pets", "Missing teeth", "Dental device", "Adenoids or tonsils removed"),
          lhei = c(0.5,20),
          lwid = c(0.8,5),
          col = mypalette, 
          breaks = mycolors,
          offsetRow = 0,
          offsetCol = 0,
          margins = c(15, 15))

## ARG associations heatmap (FIGURE 6; edit by adding 'association boxes')
heatmap.2(as.matrix(halla_ARGs[,1:7]),
          Rowv = FALSE,
          Colv = FALSE,
          dendrogram="none",
          key=FALSE, 
          symkey=FALSE, 
          density.info='none',
          trace = "none",
          cexCol = 1.6,
          srtCol = 65,
          cexRow = 1.4,
          labRow = as.expression(lapply(rownames(halla_ARGs), function(a) bquote(italic(.(a))))),
          labCol = c("Dietary restrictions", "Alcohol consumption", "Coffee consumption", 
                     "Toothbrush age", "Bathroom window", "Domestic pets", "Missing teeth"),
          RowSideColors = as.character(halla_ARGs$color),
          lhei = c(0.5,20),
          lwid = c(0.8,5),
          col = mypalette, 
          breaks = mycolors,
          offsetRow = 1.5,
          offsetCol = 0,
          margins = c(15, 15))

# plot ARG legend
plot(1:50, col = "white")
legend(1,49, legend=c("fluoroquinolone", "fosfomycin", "macrolide", "peptide", 
                      "sulfonamide", "tetracycline", "triclosan", "multidrug", "other"),
       c("darkgreen", "purple", "orange", "brown",
         "darkblue", "pink", "skyblue", "yellow", "gray"), cex=1)
