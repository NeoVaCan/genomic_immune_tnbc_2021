library(clonevol)
library(fishplot)

x <- read.table("ClonEvol_DataSet.txt", header=TRUE, sep="\t")

# shorten vaf column names as they will be
vaf.col.names <- grep('.vaf', colnames(x), value=T)
var.col.names <- grep('.var', colnames(x), value=T)
ref.col.names <- grep('.ref', colnames(x), value=T)
sample.names <- gsub('.vaf', '', vaf.col.names)
x[, sample.names] <- x[, vaf.col.names] 
vaf.col.names <- sample.names
# prepare sample grouping
#sample.groups <- c("M1", "M5A", "M5B", "M5C", "M6A", "M6B", "M7", "M8A", "M8B", "M9", "M10C", "M10D", "M10E", "M11A", "M11B", "M12A", "M12B", "M15", "M16", "M17", "M18"); names(sample.groups) <- vaf.col.names
sample.groups <- c("Chest_Wall", "Chest_Wall", "M1", "Chest_Wall", "Lung_Right", "Lung_Right", "Lung_Right", "Diaphragm", "Lung_Left", "Lung_Right", "Pericardium", "Lung_Right", "Diaphragm", "Diaphragm", "Lung_Left", "Lung_Left", "Lung_Left", "Pericardium", "Lymph_Nodes", "Lymph_Nodes", "Chest_Wall"); names(sample.groups) <- vaf.col.names
# setup the order of clusters to display in various plots (later)
x <- x[order(x$cluster),]
#If you want ClonEvol to choose colors for you, simply set it to NULL, like this:
clone.colors <- c("black", "yellow3", "#56B4E9", "#009E73", "green2", "red", "gray50", "magenta" , "orange", "maroon")
#clone.colors <- NULL
#pdf('BoxPlot.pdf', width = 10, height = 30, useDingbats = FALSE, title='') 
pp <- plot.variant.clusters(x,
    cluster.col.name = 'cluster', 
    show.cluster.size = FALSE, 
    cluster.size.text.color = 'blue', 
    vaf.col.names = vaf.col.names, 
    vaf.limits = 100, 
    sample.title.size = 20,
    violin = FALSE,
    box = TRUE,
    jitter = TRUE,
    jitter.shape = 1,
    jitter.color = clone.colors, 
    jitter.size = 1,
    jitter.alpha = 1, 
    jitter.center.method = 'median', 
    jitter.center.size = 1, 
    jitter.center.color = 'darkgray', 
    jitter.center.display.value = 'none', 
    highlight = 'is.driver', 
    highlight.shape = 21,
    highlight.color = 'white', 
    highlight.fill.color = 'white', 
    highlight.note.col.name = 'gene', 
    highlight.note.size = 2, 
    order.by.total.vaf = FALSE)
#dev.off()
plot.pairwise(x, col.names = vaf.col.names, out.prefix = 'Plots/AllMets_variants.pairwise.plot',
colors = clone.colors)
pdf('FlowPlot.pdf', width=20, height=10, useDingbats=FALSE, title='') 
plot.cluster.flow(x, vaf.col.names = vaf.col.names,
                  #sample.names = c("M5A", "M5B", "M5C", "M7", "M8A", "M8B", "M9"),
                  sample.names = c("M1", "M5A", "M5B", "M5C", "M6A", "M6B", "M7", "M8A", "M8B", "M9", "M10C", "M10D", "M10E", "M11A", "M11B", "M12A", "M12B", "M15", "M16", "M17", "M18"), 
                  colors = clone.colors, y.title = "CCF (%)")
dev.off()
y = infer.clonal.models(variants = x, 
    cluster.col.name = 'cluster',
    vaf.col.names = vaf.col.names, 
    sample.groups = sample.groups,
    #sample.groups = NULL,
    cancer.initiation.model='monoclonal', 
    subclonal.test = 'bootstrap', 
    subclonal.test.model = 'non-parametric', 
    num.boots = 1000,
    founding.cluster = 1, 
    cluster.center = 'mean', 
    ignore.clusters = 2,
    clone.colors = clone.colors,
    #JBH adjusted original 0.01
    min.cluster.vaf = 0.01,
    p.value.cutoff = 0.05,
    # min probability that CCF(clone) is non-negative
    #JBH adjusted original 0.05
    sum.p = 0.05,
    # alpha level in confidence interval estimate for CCF(clone) 
    #JBH adjusted original 0.05
    alpha = 0.05,
    random.seed = 666,
    score.model.by = "probability",
    vaf.in.percent = TRUE)
y <- transfer.events.to.consensus.trees(y, 
    x[x$is.driver,],
    cluster.col.name = 'cluster', 
    event.col.name = 'gene')

y <- convert.consensus.tree.clone.to.branch(y, branch.scale = 'sqrt')

plot.clonal.models(y,
    # box plot parameters
    box.plot = FALSE,
    fancy.boxplot = TRUE,
    fancy.variant.boxplot.highlight = 'is.driver', 
    fancy.variant.boxplot.highlight.shape = 21, 
    fancy.variant.boxplot.highlight.fill.color = 'red', 
    fancy.variant.boxplot.highlight.color = 'black', 
    fancy.variant.boxplot.highlight.note.col.name = 'gene', 
    fancy.variant.boxplot.highlight.note.color = 'blue', 
    fancy.variant.boxplot.highlight.note.size = 2,
    fancy.variant.boxplot.jitter.alpha = 1, 
    fancy.variant.boxplot.jitter.center.color = 'grey50', 
    fancy.variant.boxplot.base_size = 12, 
    fancy.variant.boxplot.plot.margin = 1, 
    fancy.variant.boxplot.vaf.suffix = '.VAF',
    # bell plot parameters
    clone.shape = 'bell',
    bell.event = TRUE, 
    bell.event.label.color = 'blue', 
    bell.event.label.angle = 60, 
    clone.time.step.scale = 1, 
    bell.curve.step = 2,
    # node-based consensus tree parameters
    merged.tree.plot = TRUE, 
    tree.node.label.split.character = NULL, 
    tree.node.shape = 'circle', 
    tree.node.size = 30, 
    tree.node.text.size = 0.5, 
    merged.tree.node.size.scale = 1.25, 
    merged.tree.node.text.size.scale = 2.5, 
    merged.tree.cell.frac.ci = FALSE,
    # branch-based consensus tree parameters
    merged.tree.clone.as.branch = TRUE, 
    mtcab.event.sep.char = ',', 
    mtcab.branch.text.size = 1, 
    mtcab.branch.width = 0.75, 
    mtcab.node.size = 3, 
    mtcab.node.label.size = 1, 
    mtcab.node.text.size = 1.5,
    # cellular population parameters
    cell.plot = TRUE,
    num.cells = 100,
    cell.border.size = 0.25, 
    cell.border.color = 'black', 
    clone.grouping = 'horizontal', 
    #meta-parameters 
    scale.monoclonal.cell.frac = TRUE, 
    show.score = FALSE,
    cell.frac.ci = TRUE, 
    disable.cell.frac = FALSE, 
    # output figure parameters 
    out.dir = 'Plots', 
    out.format = 'pdf', 
    overwrite.output = TRUE,
    width = 30,
    height = 20,
    max.num.models.to.plot = 35,
    # vector of width scales for each panel from left to right 
    #panel.widths = c(3,4,2,4,2))
    panel.widths = c(3,2,2,2))

input_names = sample.names

fishplot_input <- generateFishplotInputs(y, rescale = TRUE)

vafs = clonevol:::estimate.clone.vaf(y$variants, vaf.col.names=vaf.col.names)
scales = vafs[vafs$cluster == 1, vaf.col.names]/max(vafs[vafs$cluster == 1, vaf.col.names])
for (i in 1:length(fishplot_input$cell.fractions)){
fishplot_input$cell.fractions[[i]] = (fishplot_input$cell.fractions[[i]] *
as.matrix(scales[rep(1,nrow(fishplot_input$cell.fractions[[i]])),]))
}

pdf("FishPlot.pdf", width = 8, height = 4)
for (i in 1:length(fishplot_input$cell.fractions)){
    fishplot_objects <- createFishObject(fishplot_input$cell.fractions[[i]], fishplot_input$parents[[i]], fix.missing.clones=TRUE)
    fishplot_layout = layoutClones(fishplot_objects)
    fishplot_layout = setCol(fishplot_layout, fishplot_input$clonevol.clone.colors)
    fishPlot(fishplot_layout, shape = "spline", cex.title = 0.7, 
             vlines = seq(1, length(input_names)), vlab = input_names, pad.left = 0.5)
}

dev <- dev.off()