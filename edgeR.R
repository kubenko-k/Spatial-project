BiocManager::install('edgeR')
BiocManager::install('splines', force = TRUE)

install.packages('splines')
library(edgeR)
library(ggplot2)
library(splines)
install.packages("immunarch") 

setwd("/Users/ksenia/Desktop/edgeR")

#Load the data
annotation <- read.csv('annot.csv', sep = ',', header = TRUE) 
annotation <- mutate(annotation, layer = case_when(
  layer == "L1" ~ 1,
  layer == "L2" ~ 2,
  layer == "L3" ~ 3,
  layer == "L4" ~ 4,
  layer == "L5" ~ 5,
  layer == "L6" ~ 6,
  layer == "WM" ~ 7,
))
annotation %>% mutate_at(c('layer'), as.numeric)
annotation$sample_id
s <- as.character(annotation$sample_id)\

annotation <- mutate(annotation, sample_id = case_when(
  sample_id == "151507" ~ "1",
  sample_id == "151508" ~ "1",
  sample_id == "151509" ~ "1",
  sample_id == "151510" ~ "1",
  sample_id == "151669" ~ "2",
  sample_id == "151670" ~ "2",
  sample_id == "151671" ~ "2",
  sample_id == "151672" ~ "2",
  sample_id == "151673" ~ "3",
  sample_id == "151674" ~ "3",
  sample_id == "151675" ~ "3",
  sample_id == "151676" ~ "3", 
  TRUE ~ as.character(sample_id)
))

annotation

expression <- read.csv('expr.csv', sep = ',', header = TRUE) 

#Filter low expressed genes
y <- DGEList(counts=expression, samples=annotation)
y
keep.genes <- filterByExpr(y, group=y$samples$condition, min.count=0.1, min.total.count=1)
table(keep.genes)

#Prepare data
comparison_ann <- annotation
comparison_expr <- expression
rownames(comparison_expr) <- comparison_expr$X
comparison_expr <- comparison_expr[-c(1)]

comparison_ann$layer <- factor(comparison_ann$layer)
comparison_ann$sample_id <- factor(comparison_ann$sample_id)
comparison_ann$condition <- factor(comparison_ann$condition)
head(comparison_ann)

#Filtering and normalization
y <- DGEList(counts=comparison_expr, samples=comparison_ann)
y <- y[keep.genes, , keep=FALSE]
y <- normLibSizes(y)
head(y$samples, n=10L)
summary(y$samples$norm.factors)


#Data exploration
condition <- y$samples$condition
plotMDS(y, pch=c(1:3)[condition], col=c(2:8)[layer])
legend("bottomright", legend=c(levels(layer), levels(condition)), pch=c(rep(16, 7), c(1:5)), col=c(c(2:8), rep(1, 3)), cex=0.8, ncol=2)


#Design Matrix
layer <- as.numeric(y$samples$layer)
sample_id <- y$samples$sample_id
design <- model.matrix(~ sample_id + ns(layer, df=3) + condition:ns(layer, df=3))
colnames(design) <- gsub("condition", "", colnames(design))
colnames(design) <- gsub("sample_id", "", colnames(design))
colnames(design) <- gsub("ns\\(layer_c, df = 3)", "spline_", colnames(design))
head(design)

design <- model.matrix(~ ns(layer, df=3))
layer
design <- model.matrix(~ layer)

layer
sample_id
condition
y$samples

df %>% mutate_at(c('col1', 'col2'), as.numeric)

layer <- as.factor(y$samples$layer)
sample_id <- as.factor(y$samples$sample_id)
condition <- as.factor(y$samples$condition)

#Dispersion estimation

y <- estimateDisp(y, design, robust=TRUE)
y$common.dispersion
plotBCV(y)
fit <- glmQLFit(y, design, robust=TRUE, legacy=FALSE)
colnames(fit)
plotQLDisp(fit)

#DE genes identification

qlf <- glmQLFTest(fit, coef=12:14)

write.csv(qlf$table, 'edgeR_age_sampleid.csv')
summary(decideTestsDGE(qlf))

genes <- qlf$genes

df <- qlf$table
df_all <- cbind(df, genes)
write.csv(df_all, 'edgeR_age.csv')

df_sig <- subset(df_all, df_all$PValue < 0.05)
