library(phyloseq)
library(vegan)
library(mixOmics)
library(ggplot2)
library(tictoc)
library(here)


packageVersion("phyloseq")
packageVersion("mixOmics")
packageVersion("ggplot2")
packageVersion("tictoc")
packageVersion("here")

ps <- readRDS(here("mouse_ps_filt.rds"))

##
# Alpha diversity
##
alpha_metrics <- c("Observed", "Shannon", "Simpson", "Fisher")

# boxplots
plot_richness(ps, x = "age_binned", measures = alpha_metrics) + geom_boxplot()

# violin plots
## color violins by binned age
plot_richness(ps, x = "age_binned", measures = alpha_metrics) + 
  geom_violin(aes(fill = factor(age_binned)), trim = FALSE, draw_quantiles = c(0.25, 0.5, 0.75)) + 
  geom_jitter(height = 0, width = 0.1)

## color jittered points by age
plot_richness(ps, x = "age_binned", measures = alpha_metrics) + 
  geom_violin(trim = FALSE, draw_quantiles = c(0.25, 0.5, 0.75)) + 
  geom_jitter(aes(colour=factor(age)), height = 0, width = 0.1)

# perform Kruskal-Wallis test
alphas <- estimate_richness(ps, measures = alpha_metrics)
kw <- t(sapply(alphas, function(x) unlist(kruskal.test(x~sample_data(ps)$age_binned)[c("p.value","statistic")])))
kw

# apply log(1+x) transform
qplot(log10(rowSums(otu_table(ps))),binwidth=0.2) +
  xlab("Logged counts-per-sample")
pslog <- transform_sample_counts(ps, function(x) log(1 + x))

##
# Beta diversity
##

# calculate and plot weighted UniFrac 
wuf <- distance(pslog, method = "wunifrac")
wuf.pcoa <- ordinate(pslog, method = "MDS", distance = wuf)
wuf.pcoa.egvals <- wuf.pcoa$values$Eigenvalues

plot_ordination(pslog, wuf.pcoa, color = "age_binned") + 
  labs(col = "Binned Age") + 
  coord_fixed(sqrt(wuf.pcoa.egvals[2]/wuf.pcoa.egvals[1]))

# calculate and plot Bray-Curtis
bc <- distance(pslog, method = "bray")
bc.pcoa <- ordinate(pslog, method = "MDS", distance = bc)
bc.egvals <- bc.pcoa$values$Eigenvalues

plot_ordination(pslog, bc.pcoa, color = "age_binned") + 
  labs(col = "Binned Age") + 
  coord_fixed(sqrt(bc.egvals[2]/bc.egvals[1]))

# perform PERMANOVA test on Bray-Curtis DM
adonis(bc ~ age_binned, as(sample_data(pslog), "data.frame"))



# random forest
library(caret)
library(randomForest)
library(magrittr)

sample_data(pslog)$age2 <- cut(sample_data(pslog)$age, c(0, 100, 400))
dataMatrix <- data.frame(age = sample_data(pslog)$age2, otu_table(pslog))

# take 8 mice at random to be the training set, and the remaining 4 the test set
trainingMice <- sample(unique(sample_data(pslog)$host_subject_id), size = 8)
inTrain <- which(sample_data(pslog)$host_subject_id %in% trainingMice)
training <- dataMatrix[inTrain,]
testing <- dataMatrix[-inTrain,]

rfFit <- train(age ~ ., data = training, method = "rf",
               preProc = "center", proximity = TRUE)
rfClasses <- predict(rfFit, newdata = testing)
table(rfClasses, testing$age)

# plot results
rf_prox <- cmdscale(1 - rfFit$finalModel$proximity) %>%
  data.frame(sample_data(pslog)[inTrain, ])

ggplot(rf_prox) +
  geom_point(aes(x = X1, y = X2, col = age_binned),
             size = 3, alpha = 0.7) +
  scale_color_manual(values = c("#A66EB8", "#238DB5", "#748B4F")) +
  guides(col = guide_legend(override.aes = list(size = 4))) +
  labs(col = "Binned Age", x = "Axis1", y = "Axis2")

# check taxa importance
as.vector(tax_table(ps)[which.max(importance(rfFit$finalModel)), c("Family", "Genus")])

impOtu <- as.vector(otu_table(pslog)[,which.max(importance(rfFit$finalModel))])
maxImpDF <- data.frame(sample_data(pslog), abund = impOtu)
ggplot(maxImpDF) +   geom_histogram(aes(x = abund), bins = 30) +
  facet_grid(age2 ~ .) +
  labs(x = "Abundance of discriminative bacteria", y = "Number of samples")
