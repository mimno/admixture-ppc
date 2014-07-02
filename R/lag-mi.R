library(ggplot2)
library(dplyr)

options(echo=TRUE)

args <- commandArgs(trailingOnly = TRUE)

lag.data <- read.table(args[1], header=F, col.names=c("Populations", "Population", "Lag", "MI", "Type"))

lag.data$Collection <- args[2]

replicated.distributions <- filter(lag.data, Type == "replicated") %.% group_by(Collection, Lag, Population) %.% summarize(MeanMI=mean(MI), StdDevMI=sd(MI))

lag.data <- inner_join(lag.data, replicated.distributions) %.% mutate(ZScore = (MI - MeanMI) / StdDevMI)

min.value <- min(lag.data$ZScore)
max.value <- max(lag.data$ZScore)

min.abs <- min(abs(min.value), abs(max.value))
max.abs <- max(abs(min.value), abs(max.value))

z.scale <- function(x) {
  (x - min.value) / (max.value - min.value)
}

dataset.min <- lag.data %.% summarize(min = min(MI), max = max(MI)) %.% mutate(lower = min + 0.95 * (max - min), Key = 0)

real.zscores <- filter(lag.data, Type=="real")
mean.std.zscore <- real.zscores %.% group_by(Collection, Lag)  %.% summarize(meanZ = mean(ZScore), sdZ = sd(ZScore))
real.zscores <- inner_join(real.zscores, mean.std.zscore)

real.zscores$Key <- 0

bayes.factors <- real.zscores %.% group_by(Key, Lag, Collection)  %.%
  summarize(BayesFactor = max(0, mean(2 * log(dnorm(ZScore, meanZ, sdZ) / dnorm(ZScore)))))
bayes.factors$Stars <- symnum(bayes.factors$BayesFactor, corr = FALSE, na = FALSE, 
                              cutpoints = c(0, 2, 6, 10.0, Inf), 
                              symbols = c(" ", "*", "**", "***"))
print(bayes.factors)
bayes.factors <- inner_join(bayes.factors, dataset.min)

vlines <- data.frame(XIntercept = c(7.5, 12.5, 17.5, 22.5, 27.5))

lag.data <- lag.data %.% mutate(JitteredLag = Lag + 8 * 0.1 * Population - (8 * 0.05 * Populations))

z.colors <- c("#d7191c", "#d7191c", "#fdae61", "#aaaa4a", "#abd9e9", "#2c7bb6", "#2c7bb6")
z.values <- c(-1000000, z.scale(-10), z.scale(-5), z.scale(0), z.scale(5), z.scale(10), 1000000)

p <- ggplot(filter(lag.data, Type=="replicated"), aes(JitteredLag, MI)) + 
  geom_text(data=bayes.factors, aes(label=Stars, x=Lag, y=lower), size=8, color="gray70") + 
  geom_point(alpha=0.1, color="gray50") + facet_wrap(~ Collection, nrow=1) + 
  theme_bw() +
  geom_vline(data=vlines, aes(xintercept=XIntercept), color="gray70") + 
  xlab("Lag") + xlim(3, 32)

p + geom_line(data=filter(lag.data, Type=="real"), aes(JitteredLag, MI, group=Population), color="gray70") +
  geom_point(data=filter(lag.data, Type=="real"), aes(JitteredLag, MI, color=ZScore), size=3.5, alpha=0.7) +
  scale_color_gradientn(colours=z.colors, values=z.values)


ggsave("lag-mi.pdf", width=11, height=5)

