library(ggplot2)
library(dplyr)

options(echo=TRUE)

args <- commandArgs(trailingOnly = TRUE)

gwas <- read.table(args[1], header=F, col.names=c("Populations", "Population", "BayesFactor", "Type"))

gwas$Collection <- args[2]

replicated.distributions <- filter(gwas, Type == "replicated") %.% group_by(Collection, Populations, Population) %.% summarize(MeanBayesFactor=mean(BayesFactor), StdDevBayesFactor=sd(BayesFactor))

gwas <- inner_join(gwas, replicated.distributions) %.% mutate(ZScore = (BayesFactor - MeanBayesFactor) / StdDevBayesFactor)

min.value <- min(gwas$ZScore)
max.value <- max(gwas$ZScore)

min.abs <- min(abs(min.value), abs(max.value))
max.abs <- max(abs(min.value), abs(max.value))

z.scale <- function(x) {
  (x - min.value) / (max.value - min.value)
}

z.colors <- c("#d7191c", "#d7191c", "#fdae61", "#aaaa4a", "#abd9e9", "#2c7bb6", "#2c7bb6")
z.values <- c(-1000000, z.scale(-10), z.scale(-5), z.scale(0), z.scale(5), z.scale(10), 1000000)

dataset.mins <- gwas %.% group_by(Collection) %.% summarize(min = min(BayesFactor), max = max(BayesFactor)) %.% mutate(lower = min + 0.9 * (max - min)) %.% select(Collection, lower)

real.zscores <- filter(gwas, Type=="real")
mean.std.zscore <- real.zscores %.% group_by(Collection, Populations)  %.% summarize(meanZ = mean(ZScore), sdZ = sd(ZScore))
real.zscores <- inner_join(real.zscores, mean.std.zscore)

bayes.factors <- real.zscores %.% group_by(Collection, Populations)  %.%
  summarize(BayesFactor = max(0, mean(2 * log(dnorm(ZScore, meanZ, sdZ) / dnorm(ZScore)))))
bayes.factors$Stars <- symnum(bayes.factors$BayesFactor, corr = FALSE, na = FALSE, 
                           cutpoints = c(0, 2, 6, 10.0, 10000.0), 
                           symbols = c(" ", "*", "**", "***"))
bayes.factors <- inner_join(bayes.factors, dataset.mins)

vlines <- rbind(gwas %.% group_by(Collection, Populations) %.% summarize(XIntercept = 0.5),
       	    	    gwas %.% group_by(Collection, Populations) %.% summarize(XIntercept = mean(Populations) + 0.5))

p <- ggplot(filter(gwas, Type=="replicated"), aes(Population + 1, BayesFactor)) + 
  geom_text(data=bayes.factors, aes(label=Stars, x=(Populations + 1)/ 2, y=lower), size=14, color="gray70") + 
  geom_point(alpha=0.1, color="gray50") + facet_wrap(~ Collection, nrow=2, scales="free") + 
  theme_bw() + 
  geom_vline(data=vlines, aes(xintercept=XIntercept), alpha=0.0) +
  xlab("Population") +
  scale_x_continuous(breaks=1:6)
  #xlim(0, Populations)
  

p + geom_point(data=filter(gwas, Type=="real"), aes(Population + 1 - (0.04 * Populations), BayesFactor, color=ZScore), size=3.5, alpha=1.0) +
  scale_color_gradientn(colours=z.colors, values=z.values)
  
ggsave("gwas.pdf", width=6, height=4)
