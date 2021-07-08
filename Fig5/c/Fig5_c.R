library(ggplot2)

df <- read.csv("ELISpot_Data.txt", header = TRUE, sep = "\t")

positions <- c("MNRRPILTIIN", "NRRPILTIINT", "RRPILTIINTG", "RPILTIINTGR", "PILTIINTGRL", "ILTIINTGRLQ", "LTIINTGRLQW", "NRRPILTIIN", "RRPILTIINT", "RPILTIINTG", "PILTIINTGR", "ILTIINTGRL", "LTIINTGRLQ", "TIINTGRLQW", "RRPILTIIN", "RPILTIINT", "PILTIINTG", "ILTIINTGR", "LTIINTGRL", "TIINTGRLQ", "IINTGRLQW")

pdf('BarPlot.pdf', height = 5, width = 8, useDingbats = FALSE)

p<-ggplot(data=df, aes(x=Peptide, y=Average)) +
  geom_bar(stat="identity", fill="steelblue") + 
  theme_minimal() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  scale_x_discrete(limits = positions) + 
  scale_y_continuous(name="IFNg spots per 5 x 104 cells", limits=c(0, 100)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  geom_hline(yintercept = 50, colour = "red", linetype = "dashed")
p

dev.off()