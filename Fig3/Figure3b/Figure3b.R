
rm(list=ls())

library(ggplot2)
library(ggpubr)
library(dplyr)
library(ggthemes)
library(ggupset)


IO<-read.delim(file="IO_signatures.txt", header = TRUE, sep = "\t", row.names = 1)
IO$Date <- as.Date(IO$Date, format = "%d-%b-%y")


RESP= read.table(file="Response_duration.txt", header = TRUE, sep = "\t")
RESP$Start <- as.Date(RESP$Start, format = "%Y-%m-%d")
RESP$End <- as.Date(RESP$End, format = "%Y-%m-%d")

# Subset chest wall data
IO_thoracic<-subset(IO, Tissue %in% "chest wall")


# Plot longitudinals
a <-ggplot(data = IO_thoracic, mapping = aes(x = Date, y = APM))+
  scale_x_date(breaks=as.Date(c("2013-11-28","2015-01-28","2017-01-17","2017-07-04", "2018-02-02","2018-06-15")), date_labels = c("d373","d799","d1,519","d1,687","d1,900", "d2,033"))+
  geom_point(aes(y = APM, colour = "APM"))+
  stat_summary(fun = mean, geom = 'line', colour = "blue", aes(y = APM, group = 1))+
  geom_point(aes(y = IFN.Gamma, color = "IFN-gamma signaling"))+
  stat_summary(fun = mean, geom = 'line', colour = "black", aes(y = IFN.Gamma, group = 1))+
  geom_point(aes(y = Prolif, color = "Proliferation"))+
  stat_summary(fun = mean, geom = 'line', colour = "red", aes(y = Prolif, group = 1))+
  scale_linetype_manual(name ="", values = c('dashed')) +
  geom_rect(data = RESP[2:9,], aes(xmin = as.Date(Start,"%Y-%m-%d"), xmax = as.Date(End,"%Y-%m-%d"), fill = Response), 
            ymin = -Inf, ymax = Inf, alpha = 0.2,
            inherit.aes = FALSE, show.legend = FALSE)+
  scale_fill_manual(values = c("green","grey", "red", "orange", "yellow" )) +
  scale_color_manual(values = c(
    'APM' = 'blue',
    'IFN-gamma signaling' = 'black',
    'Proliferation' = 'red'),
    labels = c(
      'APM' = 'APM',
      'IFN-gamma signaling' = expression(IFN~gamma~ signaling),
      'Proliferation' = 'Tumor proliferation'
    )) +
  labs(color = 'Signatures scores', y = "Signatures scores")+
  geom_hline(aes(yintercept = 6.5, linetype= "TIS, primary TNBC (TCGA)"),colour = 'grey47')+
  theme_classic()+
  theme(text = element_text(family="Helvetica"), legend.position="right",axis.title.x=element_blank(),
        axis.text.x=element_blank(),legend.text.align = 0, plot.margin=margin(t = 1, b = 1, unit = "pt"))

b<-ggplot(data = IO_thoracic, mapping = aes(x = Date, y = CD8.T.Cells)) +
  scale_x_date(breaks=as.Date(c("2013-11-28","2015-01-28","2017-01-17","2017-07-04", "2018-02-02","2018-06-15")), date_labels = c("d373","d799","d1,519","d1,687","d1,900", "d2,033"))+
  geom_point(aes(y = CD8.T.Cells, color = "CD8 T Cells"))+
  geom_point(aes(y = Th1.Cells, color = "Th1 Cells"))+
  geom_point(aes(y = NK.Cells, color = "NK Cells"))+
  stat_summary(fun = mean, geom = 'line', colour = "darkred", aes(y = CD8.T.Cells, group = 1))+
  stat_summary(fun = mean, geom = 'line', colour = "darkorange", aes(y = Th1.Cells, group = 1))+
  stat_summary(fun = mean, geom = 'line', colour = "gray47", aes(y = NK.Cells, group = 1))+
  geom_rect(data = RESP[2:9,], aes(xmin = as.Date(Start,"%Y-%m-%d"), xmax = as.Date(End,"%Y-%m-%d"), fill = Response), 
            ymin = -Inf, ymax = Inf, alpha = 0.2,
            inherit.aes = FALSE, show.legend = FALSE)+
  scale_fill_manual(values = c("green","grey", "red", "orange", "yellow" ))+
  scale_color_manual(values = c(
    'CD8 T Cells' = 'darkred',
    'Th1 Cells' = 'darkorange',
    'NK Cells' = 'gray47'))+
  labs(color = 'CTL abundance', y = "CTL abundance scores")+
  theme_classic()+
  theme(text = element_text(family="Helvetica"),axis.title.x=element_blank(),legend.text.align = 0, plot.margin=margin(t = 1, b = 1, unit = "pt"))

# Facet longitudinal plots
figure3b <- ggarrange(a, b,
                    ncol = 1, nrow = 2, align = "hv") 
pdf("longitudinal3b.pdf", useDingbats = FALSE, height = 5, width = 8)
figure3b
dev.off()

