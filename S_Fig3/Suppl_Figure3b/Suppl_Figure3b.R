rm(list=ls())

library(corrplot)
library(Hmisc)
library(ggplot2)
library(pheatmap)


IO<-read.delim(file="IO_signatures.txt", header = TRUE, sep = "\t", row.names = 1)
IO$Date <- as.Date(IO$Date, format = "%d-%b-%y")


RESP= read.table(file="Response_duration.txt", header = TRUE, sep = "\t")
RESP$Start <- as.Date(RESP$Start, format = "%Y-%m-%d")
RESP$End <- as.Date(RESP$End, format = "%Y-%m-%d")

IO_thoracic<-subset(IO,(Tissue %in% "chest wall"))

# Plot longitudinals
a <-ggplot(data = IO_thoracic, aes(x = Date, y = TIS))+ylim(4, 8)+
  scale_x_date(breaks=as.Date(c("2013-11-28","2015-01-28","2017-01-17","2017-07-04", "2018-02-02","2018-06-15")), date_labels = c("d373","d799","d1,519","d1,687","d1,900", "d2,033"))+
  geom_point(aes(y = TIS, colour = "Tumor inflammation signature"))+
  stat_summary(fun = mean, geom = 'line', colour = "black", aes(y = TIS, group = 1))+
  geom_point(aes(x=Date,y = Prolif, color = "Tumor proliferation"))+
  geom_line(aes(x=Date, y=Prolif), color="red", na.rm = TRUE)+
  geom_hline(aes(yintercept = 6.5, linetype= "TIS, primary TNBC (TCGA)"),colour = 'grey47')+
  scale_linetype_manual(name ="", values = c('dashed')) +
  scale_color_manual(values = c(
    'Tumor inflammation signature' = 'black',
    'Tumor proliferation' = 'red')) +
  labs(color = 'Signatures scores')+
  geom_rect(data = RESP[2:9,], aes(xmin = as.Date(Start,"%Y-%m-%d"), xmax = as.Date(End,"%Y-%m-%d"), fill = Response), 
            ymin = -Inf, ymax = Inf, alpha = 0.2,
            inherit.aes = FALSE, show.legend = FALSE)+
  scale_fill_manual(values = c("green","grey", "red", "orange", "yellow" )) +
  theme_classic()+
  theme(text = element_text(family="Helvetica"), legend.position="right",axis.title.x=element_blank(),axis.text.x=element_blank(),axis.title.y=element_blank(), plot.margin=margin(t = 1, b = 1, unit = "pt"))


b<-ggplot(data = IO_thoracic, mapping = aes(x = Date, y = CD8.T.Cells)) +
  scale_x_date(breaks=as.Date(c("2013-11-28","2015-01-28","2017-01-17","2017-07-04", "2018-02-02","2018-06-15")), date_labels = c("d373","d799","d1,519","d1,687","d1,900", "d2,033"))+
  geom_point(aes(y = CD8.T.Cells, color = "CD8 T Cells, T Cells, NK Cells, \n NK CD56dim Cells"))+
  geom_point(aes(y = T.Cells, color = "CD8 T Cells, T Cells, NK Cells, \n NK CD56dim Cells"))+
  geom_point(aes(y = NK.Cells, color = "CD8 T Cells, T Cells, NK Cells, \n NK CD56dim Cells"))+
  geom_point(aes(y = NK.CD56dim.Cells, color = "CD8 T Cells, T Cells, NK Cells, \n NK CD56dim Cells"))+
  stat_summary(fun = mean, geom = 'line', colour = "lightgray", aes(y = CD8.T.Cells, group = 1))+
  stat_summary(fun = mean, geom = 'line', colour = "lightgray", aes(y = T.Cells, group = 1))+
  stat_summary(fun = mean, geom = 'line', colour = "lightgray", aes(y = NK.Cells, group = 1))+
  stat_summary(fun = mean, geom = 'line', colour = "lightgray", aes(y = NK.CD56dim.Cells, group = 1))+
  geom_point(aes(y = Th1.Cells, color = "Th1 Cells"))+
  geom_point(aes(y =CD45, colour = "CD45"))+
  geom_point(aes(y = Cytotoxic.Cells, color = "Cytotoxic Cells"))+
  geom_point(aes(y = Inflammatory.Chemokines, color = "Inflammatory Chemokines"))+
  stat_summary(fun = mean, geom = 'line', colour = "gray47", aes(y = Th1.Cells, group = 1))+
  stat_summary(fun = mean, geom = 'line', colour = "red", aes(y = CD45, group = 1))+
  stat_summary(fun = mean, geom = 'line', colour = "#3366CC", aes(y = Cytotoxic.Cells, group = 1))+
  stat_summary(fun = mean, geom = 'line', colour = "darkgreen", aes(y = Inflammatory.Chemokines, group = 1))+
  geom_rect(data = RESP[2:9,], aes(xmin = as.Date(Start,"%Y-%m-%d"), xmax = as.Date(End,"%Y-%m-%d"), fill = Response), 
            ymin = -Inf, ymax = Inf, alpha = 0.2,
            inherit.aes = FALSE, show.legend = FALSE)+
  scale_fill_manual(values = c("green","grey", "red", "orange", "yellow" ))+
  scale_color_manual(values = c(
    'CD8 T Cells, T Cells, NK Cells, \n NK CD56dim Cells' = 'lightgray',
    'Th1 Cells' = 'gray47',
    'CD45' = 'red',
    'Cytotoxic Cells' = '#3366CC',
    'Inflammatory Chemokines' = 'darkgreen'))+
  labs(color = 'Effector immune cell abundance')+
  theme_classic()+
  theme(text = element_text(family="Helvetica"), axis.title.y=element_blank(), axis.title.x=element_blank(),axis.text.x=element_blank(), plot.margin=margin(t = 1, b = 1, unit = "pt"))

c<-ggplot(data = IO_thoracic, mapping = aes(x = Date, y = APM)) +labs(x = "Date")+
  scale_x_date(breaks=as.Date(c("2013-11-28","2015-01-28","2017-01-17","2017-07-04", "2018-02-02","2018-06-15")), date_labels = c("d373","d799","d1,519","d1,687","d1,900", "d2,033"))+
  geom_point(aes(y = DC, color = "DC and Macrophages"))+
  geom_point(aes(y = Macrophages, color = "DC and Macrophages"))+
  geom_point(aes(y = B.Cells, color = "B Cells"))+
  stat_summary(fun = mean, geom = 'line', colour = "lightgrey", aes(y = DC, group = 1))+
  stat_summary(fun = mean, geom = 'line', colour = "lightgrey", aes(y = Macrophages, group = 1))+
  stat_summary(fun = mean, geom = 'line', colour = "gray47", aes(y = B.Cells, group = 1))+
  geom_point(aes(y = APM, colour = "MHC Class I (APM)"))+
  geom_point(aes(y = Immunoproteasome, color = "Immunoproteasome"))+
  geom_point(aes(y = MHC2, color = "MHC Class II"))+
  stat_summary(fun = mean, geom = 'line', colour = "magenta", aes(y = Immunoproteasome, group = 1))+
  stat_summary(fun = mean, geom = 'line', colour = "darkred", aes(y = APM, group = 1))+
  stat_summary(fun = mean, geom = 'line', colour = "#000066", aes(y = MHC2, group = 1))+
  geom_rect(data = RESP[2:9,], aes(xmin = as.Date(Start,"%Y-%m-%d"), xmax = as.Date(End,"%Y-%m-%d"), fill = Response), 
            ymin = -Inf, ymax = Inf, alpha = 0.2,
            inherit.aes = FALSE, show.legend = FALSE)+
  scale_fill_manual(values = c("green","grey", "red", "orange", "yellow" ))+
  scale_color_manual(values = c(
    'DC and Macrophages' = 'lightgrey',
    'B Cells' = 'gray47',
    'Immunoproteasome' = 'magenta',
    'MHC Class I (APM)' = 'darkred',
    'MHC Class II' = '#000066'))+
  labs(color = 'Antigen processing/presentation')+
  theme_classic()+
  theme(text = element_text(family="Helvetica"),axis.title.y=element_blank(),axis.title.x=element_blank(),axis.text.x=element_blank(),plot.margin=margin(t = 1, b = 1, unit = "pt"))


# Make longitudinal of CD4 and CD8 IHC to add as d

suppl.data<-read.delim("Suppl.Table2.txt")

cd4.cd8<- subset(suppl.data[, c(1:4,18,19,21:26)], Specimen %in% c("FFPE","Frozen_tumor/FFPE"))
cd4.cd8$CD4<-gsub("NE","0",cd4.cd8$CD4)
cd4.cd8$CD8<-gsub("NE","0",cd4.cd8$CD8)
cd4.cd8=type.convert(cd4.cd8)
cd4.cd8$Date<-as.Date(cd4.cd8$Date, format = "%d-%b-%y") 

d<-ggplot(data = cd4.cd8[2:33,], mapping = aes(x = Date, y = CD4)) +labs(x = "Date")+
  scale_x_date(breaks=as.Date(c("2013-11-28","2015-01-28","2017-01-17","2017-07-04", "2018-02-02","2018-06-15")), date_labels = c("d373","d799","d1,519","d1,687","d1,900", "d2,033"))+
  geom_point(aes(y = CD4, colour = "CD4"))+
  geom_point(aes(y = CD8, color = "CD8"))+
  stat_summary(fun = mean, geom = 'line', colour = "blue", aes(y = CD4, group = 1))+
  stat_summary(fun = mean, geom = 'line', colour = "black", aes(y = CD8, group = 1))+
  geom_rect(data = RESP[2:9,], aes(xmin = as.Date(Start,"%Y-%m-%d"), xmax = as.Date(End,"%Y-%m-%d"), fill = Response), 
            ymin = -Inf, ymax = Inf, alpha = 0.2,
            inherit.aes = FALSE, show.legend = FALSE)+
  scale_fill_manual(values = c("green","grey", "red", "orange", "yellow" ))+
  scale_color_manual(values = c(
    'CD4' = 'blue',
    'CD8' = 'black'))+
  labs(color = 'Cell counts (Immunohistochemistry)', y= "Signature scores")+
  theme_classic()+
  theme(text = element_text(family="Helvetica"),axis.title.x=element_blank(),plot.margin=margin(t = 1, b = 1, unit = "pt"))


# Facet longitudinal plots
figure <- ggarrange(a, b, c, d,
                    ncol = 1, nrow = 4, align = "hv", label.y = 1) 
pdf("longitudinal_Suppl3.pdf", useDingbats = FALSE, height = 7, width = 8)

figure
dev.off()
