###########################################################################################
################# Max vs hourly max wind speed and cumulative wind speed ##################
###########################################################################################

library(ggtext)
library(RColorBrewer)
library(ggplot2)
require(gridExtra)

BDobject <- readRDS(file = paste0(dir, "IIDIPUS_Input/BDobjects_v3/BD_TC20200404VUT_7325"))

display.brewer.pal(n = 11, name = "RdYlBu")
cbPalette <- brewer.pal(name = "RdYlBu", n = 11) # Defining the colour blind palette to be used in the plots

# Plot the relative frequencies of gradings
p <- ggplot(data.frame(BDobject@data$grading), aes(x=BDobject@data$grading)) +
  geom_bar(fill = cbPalette[9]) + 
  scale_fill_brewer(palette="Blues") + 
  theme_minimal()

p +
  labs(title = "Frequency of damage gradings",
       x = "Grading", y = "Frequency") +
  theme(plot.title = element_markdown(lineheight = 1.1, size = 18), #face = "bold"),
        plot.margin = unit(c(2, 1, 1, 1), "lines"))

# Create a dataframe of all the hazmeans as columns and all the buildings as rows
hazMeanDF <- data.frame(BDobject@data[6:59])

#####################################################################
##################### Total maximum wind speed ######################
#####################################################################

#
### Get a dataframe with 1. the maximum wind speeds over each row in the hazMean dataframe and 2. the grading for each row
#

maxwindspeedDF <- data.frame(maxwindspeed = rep(NA, nrow(hazMeanDF)),
                             grading = BDobject@data$grading)
for(i in 1:length(maxwindspeedDF$maxwindspeed)){
  maxwindspeedDF$maxwindspeed[i] <- max(hazMeanDF[i,], na.rm = T)
}

#
### Violin plot for each level of grading (first all of the gradings on the same plot)
#

theme_update(plot.title = element_text(hjust = 0.5)) # Centering plot titles in ggplot

p <- ggplot(maxwindspeedDF, aes(x=grading, y=maxwindspeed, fill=grading)) +
  geom_violin(trim=FALSE, scale="area", width=1.2, position = position_dodge(width = 0.2)) + 
  scale_fill_brewer(palette="Blues") + 
  theme_minimal() + 
  stat_summary(fun=mean, colour = cbPalette[2])
p +
  theme(legend.position="none") +
  labs(title = "Distribution and <b style='color:#D73027'>mean</b> of gradings",
       x = "Grading", y = "Max wind speed over all time points") +
  theme(plot.title = element_markdown(lineheight = 1.1, size = 18), #face = "bold"),
        plot.margin = unit(c(2, 1, 1, 1), "lines"))

#
### Now each with its own individual plot compared to the complement
#

# Damaged

Damaged <- ggplot(maxwindspeedDF, aes(x=grading=="Damaged", y=maxwindspeed, fill=grading=="Damaged")) +
  geom_violin(trim=FALSE, scale="area", width=1.2, position = position_dodge(width = 0.2)) + 
  scale_fill_brewer(palette="Blues") + 
  theme_minimal() + 
  stat_summary(fun=mean, colour = cbPalette[2]) +
  labs(title = "Grading==Damaged",
       y = "Max wind speed") +
  theme(plot.title = element_markdown(lineheight = 1.1, size = 14), #face = "bold"),
        plot.margin = unit(c(1, 0.5, 0.5, 0.5), "lines")) +
  theme(legend.position="none",
        axis.title.x=element_blank())

# destroyed

destroyed <- ggplot(maxwindspeedDF, aes(x=grading=="destroyed", y=maxwindspeed, fill=grading=="destroyed")) +
  geom_violin(trim=FALSE, scale="area", width=1.2, position = position_dodge(width = 0.2)) + 
  scale_fill_brewer(palette="Blues") + 
  theme_minimal() + 
  stat_summary(fun=mean, colour = cbPalette[2]) +
  labs(title = "Grading==destroyed",
       y = "Max wind speed") +
  theme(plot.title = element_markdown(lineheight = 1.1, size = 14), #face = "bold"),
        plot.margin = unit(c(1, 0.5, 0.5, 0.5), "lines")) +
  theme(legend.position="none",
        axis.title.x=element_blank())

# moderate

moderate <- ggplot(maxwindspeedDF, aes(x=grading=="moderate", y=maxwindspeed, fill=grading=="moderate")) +
  geom_violin(trim=FALSE, scale="area", width=1.2, position = position_dodge(width = 0.2)) + 
  scale_fill_brewer(palette="Blues") + 
  theme_minimal() + 
  stat_summary(fun=mean, colour = cbPalette[2]) +
  labs(title = "Grading==moderate",
       y = "Max wind speed") +
  theme(plot.title = element_markdown(lineheight = 1.1, size = 14), #face = "bold"),
        plot.margin = unit(c(1, 0.5, 0.5, 0.5), "lines")) +
  theme(legend.position="none",
        axis.title.x=element_blank())

# possible

possible <- ggplot(maxwindspeedDF, aes(x=grading=="possible", y=maxwindspeed, fill=grading=="possible")) +
  theme(legend.position="none") +
  geom_violin(trim=FALSE, scale="area", width=1.2, position = position_dodge(width = 0.2)) + 
  scale_fill_brewer(palette="Blues") + 
  theme_minimal() + 
  stat_summary(fun=mean, colour = cbPalette[2]) +
  labs(title = "Grading==possible",
       y = "Max wind speed") +
  theme(plot.title = element_markdown(lineheight = 1.1, size = 14), #face = "bold"),
        plot.margin = unit(c(1, 0.5, 0.5, 0.5), "lines")) +
  theme(legend.position="none",
        axis.title.x=element_blank())

grid.arrange(Damaged, destroyed, moderate, possible, ncol=2, nrow=2)

#####################################################################
################## Cumulative maximum wind speed ####################
#####################################################################

#
### Get the cumulative max wind speed rather than the maximum of the maximum 
# i.e. sum each row in the hazMean dataframe rather than taking the maximum of each row

cltwindspeedDF <- data.frame(cltwindspeed = rep(NA, nrow(hazMeanDF)),
                             grading = BDobject@data$grading)
for(i in 1:length(cltwindspeedDF$cltwindspeed)){
  cltwindspeedDF$cltwindspeed[i] <- sum(hazMeanDF[i,], na.rm = T)
}

p <- ggplot(cltwindspeedDF, aes(x=grading, y=cltwindspeed, fill=grading)) +
  geom_violin(trim=FALSE, scale="area", width=1.3, position = position_dodge(width = 0.2)) + 
  scale_fill_brewer(palette="Blues") + 
  theme_minimal() + stat_summary(fun=mean, colour = cbPalette[2])

p +
  theme(legend.position="none") +
  labs(title = "Distribution and <b style='color:#D73027'>mean</b> of gradings",
       x = "Grading", y = "Cumulative wind speed over all time points") +
  theme(plot.title = element_markdown(lineheight = 1.1, size = 18), #face = "bold"),
        plot.margin = unit(c(2, 1, 1, 1), "lines"))

#
### Each individually
#

# Damaged

Damaged <- ggplot(cltwindspeedDF, aes(x=grading=="Damaged", y=cltwindspeed, fill=grading=="Damaged")) +
  geom_violin(trim=FALSE, scale="area", width=1.2, position = position_dodge(width = 0.2)) + 
  scale_fill_brewer(palette="Blues") + 
  theme_minimal() + 
  stat_summary(fun=mean, colour = cbPalette[2]) +
  labs(title = "Grading==Damaged",
       y = "Cumulative wind speed") +
  theme(plot.title = element_markdown(lineheight = 1.1, size = 14), #face = "bold"),
        plot.margin = unit(c(1, 0.5, 0.5, 0.5), "lines")) +
  theme(legend.position="none",
        axis.title.x=element_blank())

# destroyed

destroyed <- ggplot(cltwindspeedDF, aes(x=grading=="destroyed", y=cltwindspeed, fill=grading=="destroyed")) +
  geom_violin(trim=FALSE, scale="area", width=1.2, position = position_dodge(width = 0.2)) + 
  scale_fill_brewer(palette="Blues") + 
  theme_minimal() + 
  stat_summary(fun=mean, colour = cbPalette[2]) +
  labs(title = "Grading==destroyed",
       y = "Cumulative wind speed") +
  theme(plot.title = element_markdown(lineheight = 1.1, size = 14), #face = "bold"),
        plot.margin = unit(c(1, 0.5, 0.5, 0.5), "lines")) +
  theme(legend.position="none",
        axis.title.x=element_blank())

# moderate

moderate <- ggplot(cltwindspeedDF, aes(x=grading=="moderate", y=cltwindspeed, fill=grading=="moderate")) +
  geom_violin(trim=FALSE, scale="area", width=1.2, position = position_dodge(width = 0.2)) + 
  scale_fill_brewer(palette="Blues") + 
  theme_minimal() + 
  stat_summary(fun=mean, colour = cbPalette[2]) +
  labs(title = "Grading==moderate",
       y = "Cumulative wind speed") +
  theme(plot.title = element_markdown(lineheight = 1.1, size = 14), #face = "bold"),
        plot.margin = unit(c(1, 0.5, 0.5, 0.5), "lines")) +
  theme(legend.position="none",
        axis.title.x=element_blank())

# possible

possible <- ggplot(cltwindspeedDF, aes(x=grading=="possible", y=cltwindspeed, fill=grading=="possible")) +
  theme(legend.position="none") +
  geom_violin(trim=FALSE, scale="area", width=1.2, position = position_dodge(width = 0.2)) + 
  scale_fill_brewer(palette="Blues") + 
  theme_minimal() + 
  stat_summary(fun=mean, colour = cbPalette[2]) +
  labs(title = "Grading==possible",
       y = "Cumulative wind speed") +
  theme(plot.title = element_markdown(lineheight = 1.1, size = 14), #face = "bold"),
        plot.margin = unit(c(1, 0.5, 0.5, 0.5), "lines")) +
  theme(legend.position="none",
        axis.title.x=element_blank())

grid.arrange(Damaged, destroyed, moderate, possible, ncol=2, nrow=2)
