library(ggplot2)

################################################################################
#############  ***Dictyota***  ##############

dictyota_pie_1 <- data.frame(
  group = c("Female_Biased", "Male_Biased", "Unbiased"),
  value = c(2.5,3.2,94.3))

# Barplot
bp_Dictyota<- ggplot(dictyota_pie_1, aes(x="", y=value, fill=group))+
  geom_bar(width = 1, stat = "identity") +
  ggtitle("Dictyota") +
  theme(plot.title = element_text(face="bold", hjust = 0.5))
bp_Dictyota

#create pie
pie_Dictyota <- bp_Dictyota + coord_polar("y", start=0)
pie_Dictyota

## Add FCs ######
dictyota_pie_2<- data.frame(
  group = c("FB-FC > 4", "FB-FC > 16", "MB-FC > 4", "MB-FC > 16","Unbiased"),
  value = c(1.71,0.86,2.34,1.25,94.3))

Barplot
bp_Dictyota2<- ggplot(dictyota_pie_2, aes(x="", y=value, fill=group))+
  geom_bar(width = 1, stat = "identity") +
  ggtitle("Dictyota") +
  theme(plot.title = element_text(face="bold", hjust = 0.5))
bp_Dictyota2

#create pie
pie_Dictyota2 <- bp_Dictyota2 + coord_polar("y", start=0)
pie_Dictyota2

################################################################################
############## ***Myriotrichia***  ############

Myrio_pie_1 <- data.frame(
  group = c("Female_Biased", "Male_Biased", "Unbiased"),
  value = c(9.0,5.4,86.1))

# Barplot
bp_Myrio<- ggplot(Myrio_pie_1, aes(x="", y=value, fill=group))+
  geom_bar(width = 1, stat = "identity") +
  ggtitle("Myriotrichia") +
  theme(plot.title = element_text(face="bold", hjust = 0.5))
bp_Myrio

#create pie
pie_Myrio <- bp_Myrio + coord_polar("y", start=0)
pie_Myrio

## Add FCs ######
Myrio_pie_2<- data.frame(
  group = c("FB-FC > 4", "FB-FC > 16", "MB-FC > 4", "MB-FC > 16","Unbiased"),
  value = c(5.3,2.7,3.2,1.6,86.1))

#Barplot
bp_Myrio2<- ggplot(Myrio_pie_2, aes(x="", y=value, fill=group))+
  geom_bar(width = 1, stat = "identity") +
  ggtitle("Myriotrichia") +
  theme(plot.title = element_text(face="bold", hjust = 0.5))

bp_Myrio2

#create pie
pie_Myrio2 <- bp_Myrio2 + coord_polar("y", start=0)
pie_Myrio2

################################################################################
############## ***Desmarestia***  ############

Desma_pie_1 <- data.frame(
  group = c("Female_Biased", "Male_Biased", "Unbiased"),
  value = c(8.7,14.1,77.1))

#Barplot
bp_Desma<- ggplot(Desma_pie_1, aes(x="", y=value, fill=group))+
  geom_bar(width = 1, stat = "identity") +
  ggtitle("Desmarestia") +
  theme(plot.title = element_text(face="bold", hjust = 0.5))
bp_Desma

#create pie
pie_Desma <- bp_Desma + coord_polar("y", start=0)
pie_Desma

## Add FCs ######
Desma_pie_2<- data.frame(
  group = c("FB-FC > 4", "FB-FC > 16", "MB-FC > 4", "MB-FC > 16","Unbiased"),
  value = c(2.3,0.6,3.7,1.1,77.1))

#Barplot
bp_Desma2<- ggplot(Desma_pie_2, aes(x="", y=value, fill=group))+
  geom_bar(width = 1, stat = "identity") +
  ggtitle("Desmarestia") +
  theme(plot.title = element_text(face="bold", hjust = 0.5))
bp_Desma2

#create pie
pie_Desma2 <- bp_Desma2 + coord_polar("y", start=0)
pie_Desma2

################################################################################
############## ***Sacchoriza***  ############

Saccho_pie_1 <- data.frame(
  group = c("Female_Biased", "Male_Biased", "Unbiased"),
  value = c(10.5,12.0,78.0))

#Barplot
bp_Saccho<- ggplot(Saccho_pie_1, aes(x="", y=value, fill=group))+
  geom_bar(width = 1, stat = "identity") +
  ggtitle("Sacchoriza") +
  theme(plot.title = element_text(face="bold", hjust = 0.5))
bp_Saccho

#create pie
pie_Saccho <- bp_Saccho + coord_polar("y", start=0)
pie_Saccho

## Add FCs ######
Saccho_pie_2<- data.frame(
  group = c("FB-FC > 4", "FB-FC > 16", "MB-FC > 4", "MB-FC > 16","Unbiased"),
  value = c(4.6,2.3,5.5,3.3,78))

#Barplot
bp_Saccho2<- ggplot(Saccho_pie_2, aes(x="", y=value, fill=group))+
  geom_bar(width = 1, stat = "identity") +
  ggtitle("Sacchoriza") +
  theme(plot.title = element_text(face="bold", hjust = 0.5))
bp_Saccho2

#create pie
pie_Saccho2 <- bp_Saccho2 + coord_polar("y", start=0)
pie_Saccho2
