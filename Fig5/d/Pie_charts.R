#Pie charts 

pdf('Pie_Charts.pdf', width = 5, height = 5)

slices <- c(92, 14, 85, 48)
lbls <- c("clonal [early]", "clonal [NA]", "clonal [late]", "subclonal")
pie(slices, labels = lbls, col=c("gray10", "gray30", "gray50", "yellow"), main="M5A CLS")

slices <- c(128, 9, 41, 43)
lbls <- c("clonal [early]", "clonal [NA]", "clonal [late]", "subclonal")
pie(slices, labels = lbls, col=c("gray10", "gray30", "gray50", "yellow"), main="M5B CLS")

slices <- c(88, 9, 110, 93)
lbls <- c("clonal [early]", "clonal [NA]", "clonal [late]", "subclonal")
pie(slices, labels = lbls, col=c("gray10", "gray30", "gray50", "yellow"), main="M5C CLS")

slices <- c(89, 19, 109, 3)
lbls <- c("clonal [early]", "clonal [NA]", "clonal [late]", "subclonal")
pie(slices, labels = lbls, col=c("gray10", "gray30", "gray50", "yellow"), main="M6A CLS")

slices <- c(97, 24, 140, 6)
lbls <- c("clonal [early]", "clonal [NA]", "clonal [late]", "subclonal")
pie(slices, labels = lbls, col=c("gray10", "gray30", "gray50", "yellow"), main="M6B CLS")

slices <- c(129, 2, 65, 24)
lbls <- c("clonal [early]", "clonal [NA]", "clonal [late]", "subclonal")
pie(slices, labels = lbls, col=c("gray10", "gray30", "gray50", "yellow"), main="M7 CLS")

slices <- c(124, 5, 24, 21)
lbls <- c("clonal [early]", "clonal [NA]", "clonal [late]", "subclonal")
pie(slices, labels = lbls, col=c("gray10", "gray30", "gray50", "yellow"), main="M8A CLS")

slices <- c(91, 71, 26)
lbls <- c("clonal [early]", "clonal [late]", "subclonal")
pie(slices, labels = lbls, col=c("gray10", "gray50", "yellow"), main="M8B CLS")

slices <- c(103, 24, 79, 18)
lbls <- c("clonal [early]", "clonal [NA]", "clonal [late]", "subclonal")
pie(slices, labels = lbls, col=c("gray10", "gray30", "gray50", "yellow"), main="M9 CLS")

slices <- c(106, 16, 139, 14)
lbls <- c("clonal [early]", "clonal [NA]", "clonal [late]", "subclonal")
pie(slices, labels = lbls, col=c("gray10", "gray30", "gray50", "yellow"), main="M10C CLS")

slices <- c(105, 56, 88)
lbls <- c("clonal [early]", "clonal [late]", "subclonal")
pie(slices, labels = lbls, col=c("gray10", "gray50", "yellow"), main="M10D CLS")

slices <- c(106, 9, 52, 67)
lbls <- c("clonal [early]", "clonal [NA]", "clonal [late]", "subclonal")
pie(slices, labels = lbls, col=c("gray10", "gray30", "gray50", "yellow"), main="M10E CLS")

slices <- c(69, 6, 9, 3)
lbls <- c("clonal [early]", "clonal [NA]", "clonal [late]", "subclonal")
pie(slices, labels = lbls, col=c("gray10", "gray30", "gray50", "yellow"), main="M11A CLS")

slices <- c(70, 12, 34)
lbls <- c("clonal [early]", "clonal [NA]", "clonal [late]")
pie(slices, labels = lbls, col=c("gray10", "gray30", "gray50"), main="M11B CLS")

slices <- c(108, 20, 94, 51)
lbls <- c("clonal [early]", "clonal [NA]", "clonal [late]", "subclonal")
pie(slices, labels = lbls, col=c("gray10", "gray30", "gray50", "yellow"), main="M12A CLS")

slices <- c(3, 18)
lbls <- c("clonal [early]", "subclonal")
pie(slices, labels = lbls, col=c("gray10", "yellow"), main="M12B CLS")

slices <- c(98, 2, 21)
lbls <- c("clonal [early]", "clonal [NA]", "clonal [late]")
pie(slices, labels = lbls, col=c("gray10", "gray30", "gray50"), main="M15 CLS")

slices <- c(27, 3, 18)
lbls <- c("clonal [early]", "clonal [NA]", "clonal [late]")
pie(slices, labels = lbls, col=c("gray10", "gray30", "gray50"), main="M16 CLS")

slices <- c(118, 24, 62)
lbls <- c("clonal [early]", "clonal [NA]", "clonal [late]")
pie(slices, labels = lbls, col=c("gray10", "gray30", "gray50"), main="M17 CLS")

slices <- c(98, 64, 39, 7)
lbls <- c("clonal [early]", "clonal [NA]", "clonal [late]", "subclonal")
pie(slices, labels = lbls, col=c("gray10", "gray30", "gray50", "yellow"), main="M18 CLS")

dev.off()