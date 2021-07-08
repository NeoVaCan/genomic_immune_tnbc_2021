library(dplyr)
library("ggpubr")
library(multcomp)
library(rstatix)

df <- read.csv(file = 'Neoantigens_Prop.txt', sep = "\t")
df["class"] <- as.factor(df$class)

res.aov <- aov(p.late ~ class, data = df)
summary(glht(res.aov, linfct = mcp(class = "Tukey")))

group1 <- c("Inter", "Late", "Late")
group2 <- c("Early", "Early", "Inter")
p <- c(0.8799, 0.0159, 0.0276)
y.position <- c(0.28,0.82,0.40)

stat.test <- data.frame(group1, group2, p, y.position)

plot2 <- ggboxplot(df, x = "class", y = "p.late",
color = "class", palette = c("chartreuse", "blue1", "darkviolet"),
order = c("Early", "Inter", "Late"),
ylab = "p.late", xlab = "MolClocl_Timing")
plot2 <- plot2 + stat_pvalue_manual(
  stat.test,  label = "p", tip.length = 0
  )

#Plot both in one PDF

pdf('BoxPlot.pdf', width = 8, height = 4)
grid.arrange(plot1, plot2, ncol=2)
dev.off()