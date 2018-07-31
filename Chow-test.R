####
# Created by Mateu Menendez-Serra
####

setwd("insert the path of your directory")

#install.packages("strucchange")
#install.packages("ggplot2")
#install.packages("gridExtra")
library(strucchange) #https://CRAN.R-project.org/package=strucchange
library(ggplot2)
library(gridExtra)

# an example data is provided in the following file: "clustering\_index\_vs\_env\_gradient.csv" 
data <- read.table(file = "clustering\_index\_vs\_env\_gradient.csv", sep = ";", header = T)
data$Salinity <- log10(data$Salinity)
data <- data[order(data$Salinity), ] # is mandatory to order the data according to the salinity values 
                                    # (increasing or decreasing, indistinctly) 

# Execute the Chow test with the "sctest" function from the package "strucchange" to each point of the salinity 
# gradient, excluding the tails (the function needs a minimum number of 3 points on both sides of the tested value, 
# in this case we have eliminated 5 values of each tail), and extract the p.value and the F statistic value.

for (i in (nrow(data)-5):5){ 
  data[i, 3] <- sctest(Clustering_index ~ Salinity, type = "Chow", point = i, data = data)$p.value   # Add the results to the existing
  data[i, 4] <- sctest(Clustering_index ~ Salinity, type = "Chow", point = i, data = data)$statistic # data object a new column
}

colnames(data)[3:4] <- c("pvalue", "Statistic")

max.statistic <- max(data$Statistic, na.rm = T) # Get the maximum value of the statistic
max.observation <- which.max(data$Statistic) # Get observation with the maximum value of the statistic
pvalue <- data[max.observation, 3] # Get the p.value corresponding to the maximum value of the statistic
threshold.value <- data[max.observation, 1] # Get the salinity value corresponding to the maximum value of the statistic 
thershold.salinity.value <- 10^data[max.observation, 1]

# Plot the variation of the "Response variable" as well as the value of the F statistic
# and the p.value along the salinity gradient. 

g_chow.sal <- ggplot(data = data, aes(x = 10^Salinity, y = Clustering_index))+
  theme_bw()+
  theme(aspect.ratio = 4/4, legend.position = "bottom")+ 
  labs(x="Max. Salinity (% dissolved salts w/v)",  y ="Clustering index", caption = "")+
  geom_point(alpha = 0.2, size = 4, colour = "black")+
  scale_x_log10()+
  geom_smooth(method ="loess", span = 0.4, se = T, colour = "green", size = 1.5)+
  geom_vline(xintercept = thershold.salinity.value, colour = "orange", linetype = "dotted", size = 0.75)
g_chow.sal

g_chow.stat <- ggplot(data = data, aes(x = 10^Salinity, y = Statistic))+
  theme_bw()+
  theme(aspect.ratio = 4/4, legend.position = "bottom")+
  labs(x="Max. Salinity (% dissolved salts w/v)", y = "Structural Change Test F Statistic\n -Chow test-", caption = paste("Threshold salinity value:", thershold.salinity.value, "% (w/V)"))+
  geom_point(alpha = 0.2, size = 4, colour = "black")+
  scale_x_log10()+
  geom_smooth(method = "loess", span = 0.2, se = T, colour = "blue", size = 1.5)+
  geom_vline(xintercept = thershold.salinity.value, colour = "orange", linetype = "dotted", size = 0.75)
g_chow.stat

g_chow.sal.pval <- ggplot(data, aes(x = 10^Salinity, y = log10(pvalue)))+
  theme_bw()+
  theme(aspect.ratio = 4/4, legend.position = "bottom")+
  labs(x="Max. Salinity (% dissolved salts w/v)", y="log transformed p-values", caption = "")+
  ylim(c(-16, -6))+
  geom_point(alpha=0.2, size=4, colour="red")+
  scale_x_log10()+
  geom_vline(xintercept = thershold.salinity.value, colour="orange", linetype = "dotted", size=0.75)
g_chow.sal.pval

# Aggregation of the three plots 
grid.arrange(g_chow.sal, g_chow.stat, g_chow.sal.pval, nrow=1, top = "Method RTCC: Randomized trait community clustering index\nEstimating structural changes in regressions models - Chow test to infer thresholds")
