# KEGG Gene Set Enrichment Analysis and Visualization of RNA-Seq data

The first goal of this exercise is to test for KEGG Pathways that are enriched for genes with divergent expression patterns between two experimental groups. Secondly, we will visually represent high dimensional RNA-seq data in meaningful ways, including ordination (non-metric multidimensional scaling) and heatmaps organized by gene- and sample-wise clustering. We will also perform a simple statistical test for multivariate differences in transcriptome composition between two experimental groups.
# Instructions

Don’t forget to answer specific questions (those in bold italics in this document), and include all code you are asked to generate, in a single, well-commented .RMD file. Label plots and axes appropriately! Also turn in an .html or .pdf file (properly knitted) from your .RMD file.

# The Data

We will mostly be working with RNA-seq data generated from male brood pouch tissues of 6 pregnant and 6 non-pregnant pipefish. The goal is to understand whether the pouch transcriptome changes significantly during pregnancy and if so, which gene regulatory pathways are affected the most.

We will also take a short look again at the stickleback gut RNA-seq data, this time to investigate gene expression differences between males and females.

All data files may be found on Talapas in: /projects/bgmp/shared/Bi623/assign8
1.	Load the libraries gage, gageData, and pathview, and read in the edgeR differential gene expression results for pipefish.

Note that this file contains a lot of other information. For our purposes now, we really only care about the “ko_ID” and “logFC” fields, which are KEGG gene orthology assignments for pipefish genes and log2 fold change values (fold change here is pregnant mean / non-pregnant mean), respectively.
```
kegg_pouch <- read.delim("pouch_RNAseq.tsv", sep = "\t",    
    stringsAsFactors=FALSE)

head(kegg_pouch)
```
2.	Set up and perform GSEA for all KEGG pathways using gage

First we need to load the KEGG orthology data from the gageData database
```
data(kegg.sets.ko)
data(sigmet.idx.ko)
kegg.sets.ko <- kegg.sets.ko[sigmet.idx.ko]
head(kegg.sets.ko, 3)
```

#### What class of object is kegg.sets.ko, and what kind of information does it contain?

Now we construct a vector of log2 fold change values and use the names() function to name each component with a ko ID (if it has one). This will allow us to do GSEA for all of the KEGG pathways, based on log2 fold change as our per-gene measure of differential expression.

```pouch_foldchanges <- kegg_pouch$logFC
names(pouch_foldchanges) <- kegg_pouch$ko_ID
head(pouch_foldchanges)
```

We now run gage() on our vector of log2 fold change values to test, for each KEGG pathway, enrichment for genes with extreme values.
```
pouch_test <- gage(pouch_foldchanges, gsets=kegg.sets.ko, 
                   same.dir=FALSE)
```

Note that we perform a “two-tailed” test, meaning a significantly enriched KEGG pathway can be defined by extreme fold changes in either direction.

Note that pouch_test is a list. We can call lapply() on the head() function to look at the top entries of the two elements in our list.
```
lapply(pouch_test, head)
head(pouch_test$greater, 30)
```

Note that the columns of pouch_test$greater provide us our hypothesis test results, including raw p-values, FDR-corrected p-values (“q-values”), sizes of the gene sets, etc.

Use subset() with a logical operator to return only those pathways such that our false discovery rate (FDR) is controlled at 0.1

### Which KEGG pathways are enriched for genes with exceptional pregnancy-specific gene expression? Which two of these are related to the immune system?

3.	Visualize pregnancy fold change magnitudes for all genes in the “coagulation and complement” KEGG pathway

First, get this pathway and its ko ID from the gage() output by simple indexing and the substr() function.
```
pouch_pathways <- rownames(pouch_test$greater)[4]
pouch_ids <- substr(pouch_pathways, start=1, stop=7)
pouch_ids
```
This seems unnecessary for a single pathway, but imagine if you wanted to look at many more.

Draw the pathway with a color scale that reflects log2 fold change for each gene. The pathview() function will write 3 files to your working directory. “ko04610.pathview.png” is the one you want to have a look at.
```
pathview(gene.data=pouch_foldchanges, 
         pathway.id=pouch_ids, species="ko", new.signature=FALSE,
         trans.fun = list(gene = NULL, cpd = NULL), 
         low = list(gene = "cornflowerblue", cpd = "blue"), 
         mid = list(gene = "gray", cpd = "gray"), 
         high = list(gene = "darkorange3", cpd = "yellow"), na.col = "transparent")
```
You can change the color scale to your liking by setting the “gene” arguments accordingly. In this case we don’t have compounds (it’s just RNA-seq data), so the “cpd” arguments are irrelevant. For the plot you turn in, color the log2 fold change values on a scale from green to yellow to red.

### What does the plot tell you about components of this pathway and male pregnancy?

4.	Visualize multi-genic expression patterns in pregnant and non-pregnant pouch tissues using non-metric multidimensional scaling

load the packages vegan and MASS

Read in the normalized pipefish pouch expression values, a comma-separated file.
```
pipe_TMMvals <- read.csv("pouch_TMM_values.csv", head = TRUE,   
                           row.names = 1)
dim(pipe_TMMvals)
pipe_TMMvals <- t(pipe_TMMvals)
head(pipe_TMMvals)
```
### What happened to our data frame after using the t() function?

This orientation of data tables is the standard for multivariate analysis.

The first thing we do with our data is compute a dissimilarity matrix (distances between pipefish samples in “transcript space”) using the default for vegdist(), which is a metric called Bray-Curtis dissimilarity.

```pipe.dis <- vegdist(pipe_TMMvals)
```
We then perform multidimensional scaling to find 2 “latent” dimensions for which the distance between samples agrees well with the original distances based on 15,252 dimensions. What do you think is meant by “convergence”? Think about how nmds works, and consult the literature, if necessary.

The fit of the new, 2-D distances to the original, 15,252-D distances can be visualized in a stress plot (or Shepard plot), which plots the 2-D versus 15,252-D distance ranks.
```
pipe.mds0 <- isoMDS(pipe.dis, k=2)
stressplot(pipe.mds0, pipe.dis)
```

Now, given our 2 “new” nMDS variables, we can plot our 12 pouch samples in this 2-D space to assess whether there are broad transcriptional differences between pregnant and non-pregnant males.

First, we need to construct a data frame that links the pouch sample IDs with pregnancy status.
```
Targets <- as.data.frame(rownames(pipe_TMMvals))
Targets$PregStat <- factor(c(rep("preg",6),rep("nonpreg",6)))
colnames(Targets) <- c("ID","PregStat")
```

Then, we define the ordination plotting parameters. The ordiplot() function, as called, will produce an empty space…
```
par(mgp = c(2.5, 1, 0))
preg=as.character(Targets$PregStat)

fig <- ordiplot(pipe.mds0, main="Brood Pouches in Transcript Space",
                ylab="nMDS Dimension 2", xlab="nMDS Dimension 1", 
                font.lab=2, font.axis=2, cex.axis=.7, type="none", cex.main=1,
                xlim=(c(-.2,0.2)))
```
… and then we will add “confidence ellipses” and the individual samples as points using the ordiellipse() and points() functions.
```
ordiellipse(pipe.mds0,groups=preg,label=FALSE, lwd=2,  
              show.groups=preg[1:6], col="darkorange4",
            draw="lines")
ordiellipse(pipe.mds0,groups=preg,label=FALSE, lwd=2, 
              show.groups=preg[7:12], col="dodgerblue4",
            draw="lines")
points(fig, "sites", pch=c(rep(19,6),rep(15,6)), 
         col=c(rep("darkorange3",6),rep("cornflowerblue",6)), cex=1.5)
```
In this case the blue squares are non-pregnant males and the orange circles are pregnant males.

### What does this plot tell you about the brood pouch transcriptomes profiled in this study?

Repeat the ordination (starting with the isoMDS() call), but this time generate 3 nMDS dimensions as opposed to 2. Produce the same type of plot as above, but plot Dimension 2 vs. Dimension 3. Note: You’ll have to index dimensions 2 and 3 of the isoMDS object directly when you use the ordiellipse() function.
5.	Permutational Multivariate Analysis of Variance (perMANOVA) to test for multivariate transcriptional differences between pregnant and non-pregnant males
This is a statistical hypothesis test that asks whether dissimilarity between samples that differ in pregnancy status is large, on average, relative to dissimilarity between samples in general. Another way of thinking about it is, “What proportion of total transcriptome dissimilarity among all 12 brood pouches is explained by pregnancy status?”

The adonis() function, which runs the perMANOVA, requires otu.env as a parameter, so we define it with our Targets data frame.
```
otu.env <- Targets
adonis(pipe.dis ~ PregStat, otu.env, perm=999)
```
### Based on the output of adonis(), do we see a significant effect of pregnancy status?
The R2 value for PregStat in our perMANOVA table corresponds to the proportion of total dissimilarity explained by pregnancy status.
6.	Constructing a heatmap with clustering dendrograms for Coagulation and Complement Cascade KEGG pathway genes.

load the packages gplots, RColorBrewer, and dendextend

Read in the following file, which contains CPM expression data for just those pipefish genes that were mapped to the KEGG “Coagulation and Complement Cascade” pathway.
```
pouch_compcoag <- read.delim("CompCoag_pouch_multivar.tsv", sep = "\t", 
                             row.names = 1, header = F)

colnames(pouch_compcoag) <- c(“KO”,"name","P9","P8","P7","P6","P5","P3",
                                        "NP11","NP10","NP4","NP3","NP2","NP1")
```
Define a vector of gene names we will use later in the heatmap
```
names_compcoag <- pouch_compcoag$name
```
Now define a new data frame consisting of just the CPM columns and call it “pouch_compcoag”, but log-transform (base 2) each CPM value + 0.01.

### Why do we add 0.01 to all values?

We need to mean-center and range data (mean=0, sd=1), but must first transpose for row-wise scaling.
```
pouch_compcoag.n <- scale(t(pouch_compcoag))
Put data back in original orientation
pouch_compcoag.tn <- t(pouch_compcoag.n)
```
### What class of object are we dealing with now?

Now we calculate multivariate dissimilarity for all sample pairs, using Euclidean Distance
```
compcoag.d1 <- dist(pouch_compcoag.n, method = "euclidean", diag = FALSE, 
                    upper = FALSE)
```
Let’s look at the distances
```
round(compcoag.d1,3)
```

### Which two samples are the most dissimilar based on Euclidean Distance?

Similarly, we calculate multivariate dissimilarity for all gene pairs
```
compcoag.d2 <- dist(pouch_compcoag.tn,method = "euclidean", diag = FALSE, 
                             upper = TRUE)
```

Now we cluster samples, then genes, using Ward linkage clustering
```
compcoag.c1 <- hclust(compcoag.d1, method = "ward.D2", members = NULL)
compcoag.c2 <- hclust(compcoag.d2, method = "ward.D2", members = NULL)
```
Let’s take a look at dendrograms based on the clustering
```
par(mfrow=c(2,1),cex=0.5) # Make 2 rows, 1 col plot frame and shrink labels
plot(compcoag.c1)
plot(compcoag.c2)
```
What does the sample dendrogram tell us about pregnant (P) and non-pregnant (NP) pouch transcriptomes?

We can print the order of samples (left to right) in the tree, in case we want to rotate the dendrogram about nodes and re-order them later.
```
compcoag.c1$order
```
Finally, we set the color scale for the heatmap. 299 increments is plenty of resolution!
```
compcoag_pal <- colorRampPalette(c(“green”,"yellow","red")) 
                                (n = 299)
```
We plot the heatmap with heatmap.2 par(cex.main=0.5, cex.axis=0.5, font=2, font.axis=2) # Shrink title fonts on plot
```
heatmap.2(pouch_compcoag.tn,                       
          Colv=rotate(as.dendrogram(compcoag.c1),  
          order=c(12,11,10,9,8,7,6,5,4,3,2,1)),
          Rowv=as.dendrogram(compcoag.c2),
          labRow=names_compcoag,
          density.info="none",
          trace="none",
          scale="none",            
          col = compcoag_pal,                 
          cexRow=0.5,cexCol=0.75,
          margins=c(3,13), lwid=c(.8,3), lhei=c(.8,3),
          srtCol=45, adjCol=c(1,1),
          keysize=1.3)
```
heatmap.2() is a fairly complex plotting function. Take a little bit of time to explore what various arguments mean, and change the values to adjust the appearance of the image. You could also re-set the color gradient, if you don’t like the one above.

### Are there any groups of genes that differentiate between pregnant and non-pregnant males particularly well? If so, name those genes.
Repeat everything above, starting with
```
> pouch_compcoag.n <- scale(t(pouch_compcoag))”
```
but do not call the scale() function.

### Does your heatmap look any different? If so, explain why this might be.

7.	Constructing heatmaps with clustering dendrograms for stickleback gene expression data.
Now that you know how to generate a heatmap from gene expression data, produce 2 heatmaps for the stickleback data set stickleback_CPM.tsv.

You’ll notice that the sample names start with either an “M” (for male) or “F” (for female).

For the first heatmap, subset only genes that have the Genome_Loc value “groupXIX” AND have a Gene_Start_Pos value between 6000000 and 12000000.

(it’s easiest to use data frame indexing with logical operators or the subset() function)

For the second heatmap, take a random sample of genes the same size as your first subset. If your first set contains, say, 100 genes, you could sample this way:
```
subsample2 <- stickle_CPM[sample(nrow(stickle_CPM),100),]
```

