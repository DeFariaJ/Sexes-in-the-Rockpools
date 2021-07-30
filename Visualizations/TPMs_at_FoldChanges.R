library(tidyr)
filter_df3$FoldChange <- "Unbiased"

filter_df3$FoldChange[filter_df3$padj<.05 & abs(filter_df3$log2FoldChange) >1 & abs(filter_df3$log2FoldChange) <2] <- "2-4"
filter_df3$FoldChange[filter_df3$padj<.05 & abs(filter_df3$log2FoldChange) >=2 & abs(filter_df3$log2FoldChange) <2.59] <- "4-6"
filter_df3$FoldChange[filter_df3$padj<.05 & abs(filter_df3$log2FoldChange) >=2.59 & abs(filter_df3$log2FoldChange) <3] <- "6-8"
filter_df3$FoldChange[filter_df3$padj<.05 & abs(filter_df3$log2FoldChange) >=3 & abs(filter_df3$log2FoldChange) <3.32] <- "8-10"
filter_df3$FoldChange[filter_df3$padj<.05 & abs(filter_df3$log2FoldChange) >=3.32 & abs(filter_df3$log2FoldChange) <4] <- "10-16"
filter_df3$FoldChange[filter_df3$padj<.05 & abs(filter_df3$log2FoldChange) >=4 & abs(filter_df3$log2FoldChange) <4.32] <- "16-20"
filter_df3$FoldChange[filter_df3$padj<.05 & abs(filter_df3$log2FoldChange) >=4.32] <- ">20"

filter_df3_bp <- filter_df3[!grepl("Unbiased", filter_df3$FoldChange),]

filter_df3_bp$FoldChange<- factor(filter_df3_bp$FoldChange,levels = c("2-4", "4-6", "6-8", "8-10","10-16","16-20",">20"))
binary_heatMap_Dictyo$FoldChange <- filter_df3_bp$FoldChange

#Bar plot
ggplot(filter_df3_bp, aes(x=FoldChange)) + 
  geom_bar() +
  theme(axis.title.x = element_text(face="bold", size=16),
        axis.title.y = element_text(face="bold", size=16),
        axis.text=element_text(size=14),
        legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())


#**************************************
xxx <- binary_heatMap_Dictyo %>% 
  group_by(FoldChange) %>%
  summarise(meanF = mean(Females), meanM = mean(Males), sdF = sd(Females), sdM = sd(Males))
#**************************************

ggplot(xxx) +
  geom_point(aes(x = FoldChange, y = meanF),colour="black", size = 2) +
  geom_point(aes(x= FoldChange, y = meanM), colour = "black", size =2) +
  geom_line(aes(x= FoldChange, y = meanF), group= 1, colour = "orange", size = 1.2) +
  geom_line(aes(x= FoldChange, y = meanM), group= 1, colour = "blue", size = 1.2) +
  geom_errorbar(aes(x = FoldChange, ymin= meanF-sdF, ymax=meanF+sdF), width=0.4, size = 0.4)+
  ylab("Log2TPM (mean)")+
  ggtitle("Dictyota") +
  theme(axis.title.x = element_text(face="bold", size=16),
        axis.title.y = element_text(face="bold", size=16),
        axis.text=element_text(size=14),
        plot.title = element_text(face="bold", size=18,hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())





 
