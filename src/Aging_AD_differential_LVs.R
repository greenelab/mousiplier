library(ggplot2)
library(dplyr)
library(ggprism)
library(pheatmap)

### set the working directory
if (rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
} else {
  # If running as a script, finding the file is harder
  # https://stackoverflow.com/a/55322344/10930590
  this_file <- commandArgs() %>%
    tibble::enframe(name = NULL) %>%
    tidyr::separate(col=value, into=c("key", "value"), sep="=", fill='right') %>%
    dplyr::filter(key == "--file") %>%
    dplyr::pull(value)

  setwd(dirname(this_file))
  setwd("../data")
}


# read in data
LVs <- read.table("../output/reformatted_ageing_AD_LVs_day_as_continueous.txt", header = TRUE, colClasses = c("factor", "numeric", "factor", "factor", "numeric"))
# convert day into a numeric variable
#LVs$day <- as.numeric(as.character(LVs$day))

df <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(df) <- c("LV_ID", "treatment", "region", "pvalue")

### perform linear model (linear regression), treating day as numerical variable
for (i in 1: length(levels(LVs$LV_ID))) {
  lv <- paste("LV", as.character(i), sep="")
  for (t in levels(LVs$treatment)) {
    for (r in levels(LVs$region)) {
      temp_lv <- subset(LVs, LVs$LV_ID == lv & LVs$treatment == t & LVs$region == r)
      
      ## one way anova
      #one.way <- aov(lv_value ~ day, data = temp_lv)
      #summary(one.way)
      #pvalue <- summary(one.way)[[1]][[1,"Pr(>F)"]]
        
      model <- lm(lv_value~day, data=temp_lv)
      pvalue <- summary(model)$coefficients[2,4]
      df[nrow(df)+1, ] <- c(lv, t, r, pvalue)
      
    }
  }
}


### BH correction and write the results
df$adjusted_p <- p.adjust(df$pvalue, "BH")
write.csv(df, "../output/Ageing_AD_LVs_pvalues.txt", row.names = FALSE)
sig_df <- subset(df, df$adjusted_p < 0.05)
write.csv(sig_df, "../output/Ageing_AD_sig_LVs_pvalues.txt", row.names = FALSE)



### Venn diagram
wt_mg <-  subset(sig_df, treatment == "WT" & region == "microglia")$LV_ID
wt_as <- subset(sig_df, treatment == "WT" & region == "astrocyte")$LV_ID
ad_mg <- subset(sig_df, treatment == "AD" & region == "microglia")$LV_ID
ad_as <- subset(sig_df, treatment == "AD" & region == "astrocyte")$LV_ID
myV <- plotVenn(list(WT_microglia=wt_mg, WT_astrocyte=wt_as, AD_microglia=ad_mg, AD_astrocyte=ad_as), nCycles = 9000, labelRegions=F)


loading <- read.table("../output/filtered_Z.tsv", header = TRUE)
### WT_microglia and WT_astrocyte
mylv <- intersect(wt_mg, wt_as)
v46 <- loading[order(-loading$V46),]
top20 <- v46[1:20, 46, drop = FALSE]
top20 <- cbind(gene = rownames(top20), top20)
rownames(top20) <- NULL
top20$gene <- factor(top20$gene, levels = top20$gene[order(top20$V46, decreasing = TRUE)])
ggplot(data = top20, aes(x=gene, y=V46)) + 
  geom_bar(stat="identity") +
  ggtitle("WT_microglia and WT_astrocyte")+
  theme_classic()+
  theme(plot.title = element_text(size = 24, hjust = 0.5), line = element_blank(), plot.background = element_blank(), panel.grid.major = element_blank()) + #set the background
  theme(panel.border = element_blank()) +   #set the border
  theme(axis.text.x = element_text()) + 
  theme(axis.title = element_text(size = 24), axis.text.x = element_text(size=24, angle=45, hjust = 1), axis.text.y = element_text(size=24)) + #set the x and y lab
  ylab("V46 loadings") + xlab("") +     #set the name of x-axis and y-axis
  theme( legend.title = element_blank(),  legend.position = "none", legend.text = element_text(size = 12), legend.key.width = unit(0.5, "cm")) + #set legend
  theme(axis.line.x = element_line(color="black", size = .6),  
        axis.line.y = element_line(color="black", size = .6),
        axis.ticks.x = element_line(size = 0.5),
        axis.ticks.y = element_line(size = 0.5),
        axis.ticks.length = unit(0.2, "cm"))


### WT_astrocyte and AD_astrocyte
mylv <- intersect(wt_as, ad_as)
v155 <- loading[order(-loading$V155),]
top20 <- v155[1:20, 155, drop = FALSE]
top20 <- cbind(gene = rownames(top20), top20)
rownames(top20) <- NULL
top20$gene <- factor(top20$gene, levels = top20$gene[order(top20$V155, decreasing = TRUE)])
ggplot(data = top20, aes(x=gene, y=V155)) + 
  geom_bar(stat="identity") +
  ggtitle("WT_astrocyte and AD_astrocyte")+
  theme_classic()+
  theme(plot.title = element_text(size = 24, hjust = 0.5), line = element_blank(), plot.background = element_blank(), panel.grid.major = element_blank()) + #set the background
  theme(panel.border = element_blank()) +   #set the border
  theme(axis.text.x = element_text()) + 
  theme(axis.title = element_text(size = 24), axis.text.x = element_text(size=24, angle=45, hjust = 1), axis.text.y = element_text(size=24)) + #set the x and y lab
  ylab("V155 loadings") + xlab("") +     #set the name of x-axis and y-axis
  theme( legend.title = element_blank(),  legend.position = "none", legend.text = element_text(size = 12), legend.key.width = unit(0.5, "cm")) + #set legend
  theme(axis.line.x = element_line(color="black", size = .6),  
        axis.line.y = element_line(color="black", size = .6),
        axis.ticks.x = element_line(size = 0.5),
        axis.ticks.y = element_line(size = 0.5),
        axis.ticks.length = unit(0.2, "cm"))


### AD_microglia and AD_astrocyte
mylv <- intersect(ad_mg, ad_as)
v40 <- loading[order(-loading$V40),]
top20 <- v40[1:20, 40, drop = FALSE]
top20 <- cbind(gene = rownames(top20), top20)
rownames(top20) <- NULL
top20$gene <- factor(top20$gene, levels = top20$gene[order(top20$V40, decreasing = TRUE)])
ggplot(data = top20, aes(x=gene, y=V40)) + 
  geom_bar(stat="identity") +
  ggtitle("AD_microglia and AD_astrocyte")+
  theme_classic()+
  theme(plot.title = element_text(size = 24, hjust = 0.5), line = element_blank(), plot.background = element_blank(), panel.grid.major = element_blank()) + #set the background
  theme(panel.border = element_blank()) +   #set the border
  theme(axis.text.x = element_text()) + 
  theme(axis.title = element_text(size = 24), axis.text.x = element_text(size=24, angle=45, hjust = 1), axis.text.y = element_text(size=24)) + #set the x and y lab
  ylab("V40 loadings") + xlab("") +     #set the name of x-axis and y-axis
  theme( legend.title = element_blank(),  legend.position = "none", legend.text = element_text(size = 12), legend.key.width = unit(0.5, "cm")) + #set legend
  theme(axis.line.x = element_line(color="black", size = .6),  
        axis.line.y = element_line(color="black", size = .6),
        axis.ticks.x = element_line(size = 0.5),
        axis.ticks.y = element_line(size = 0.5),
        axis.ticks.length = unit(0.2, "cm"))



### plot for significant LVs
for (i in 1:nrow(sig_df)) {
  #print(sig_df[i,])
  lv <- sig_df[i, 1]
  t <- sig_df[i, 2]
  r <- sig_df[i, 3]
  temp_lv <- subset(LVs, LVs$LV_ID == lv & LVs$treatment == t & LVs$region == r)
  model = lm(lv_value~day, data=temp_lv)
  summary(model)
  pvalue <- as.numeric(sig_df[i,5])
  
  mylabel <- paste("adj P = ", as.character(signif(pvalue, digits = 3)), sep="") 
  mytitle <- paste(t, r, sep="_")
  outfile <- paste("/Users/szhang32/Desktop/mousiplier/data/sig_lvs_figures/", lv, t, r, sep="_")
  pdf(file = outfile,   # The directory you want to save the file in
      width = 4, # The width of the plot in inches
      height = 4)
  myplot <- ggplot(temp_lv, aes(x=day, y=lv_value, group=1)) + 
    theme_classic()+
    geom_point() +
    ggtitle(mytitle) +
    theme(plot.title = element_text(size = 20, hjust = 0.5), line = element_blank(), plot.background = element_blank(), panel.grid.major = element_blank()) + #set the background
    theme(panel.border = element_blank()) +   #set the border
    theme(axis.text.x = element_text()) + 
    theme(axis.title = element_text(size = 20), axis.text.x = element_text(size=20, angle=0, hjust = 0.5), axis.text.y = element_text(size=20)) + #set the x and y lab
    ylab(lv) + xlab("Month") +     #set the name of x-axis and y-axis
    theme( legend.title = element_blank(),  legend.position = "none", legend.text = element_text(size = 12), legend.key.width = unit(0.5, "cm")) + #set legend
    theme(axis.line.x = element_line(color="black", size = .6),  
          axis.line.y = element_line(color="black", size = .6),
          axis.ticks.x = element_line(size = 0.5),
          axis.ticks.y = element_line(size = 0.5),
          axis.ticks.length = unit(0.2, "cm")) +
    geom_smooth(formula = y ~ x, method="lm", se=FALSE) +
    annotate("text", x = 7.5, y = max(temp_lv$lv_value)*1.2, label = mylabel, size = 6)
  
  print(myplot)
  dev.off()
}

### U matrix to see the pathways for each significant LV
df_U <- read.table("../output/filtered_U.tsv")

df_sig_U <- df_U[, unique(sig_df$LV_ID)]
df_sig_U <- df_sig_U[apply(df_sig_U, 1, sum) > 0, ]
write.csv(df_sig_U, "../output/Ageing_AD_sig_U.txt", row.names = TRUE)

pheatmap(df_sig_U, treeheight_row = 0, cluster_cols = FALSE, cluster_rows = FALSE, show_rownames = TRUE)

U <- as.matrix(read.table("../output/U.tsv"))
hist(colSums(U != 0), main='Distribution of pathways per LV ', xlab='Pathways per latent variable', cex.lab=1.5, cex.axis=1.5)
dev.off()

#### Differential LV count
LV_count <- read.table("../output/Aging_AD_differential_LVs_linear_regression_count.txt", header = TRUE, sep='\t')
LV_count$sampleID <- factor(LV_count$sampleID, levels = c("WT_microglia", "WT_astrocyte", "AD_microglia", "AD_astrocyte"))
ggplot(LV_count, aes(x=sampleID, y=count)) + 
  geom_bar(stat="identity") +
  ggtitle("")+
  theme_classic()+
  theme(plot.title = element_text(size = 20, hjust = 0.5), line = element_blank(), plot.background = element_blank(), panel.grid.major = element_blank()) + #set the background
  theme(panel.border = element_blank()) +   #set the border
  theme(axis.text.x = element_text()) + 
  theme(axis.title = element_text(size = 20), axis.text.x = element_text(size=20, angle=45, hjust = 1), axis.text.y = element_text(size=20)) + #set the x and y lab
  ylab("Differential LV count") + xlab("") +     #set the name of x-axis and y-axis
  theme( legend.title = element_blank(),  legend.position = "none", legend.text = element_text(size = 12), legend.key.width = unit(0.5, "cm")) + #set legend
  theme(axis.line.x = element_line(color="black", size = .6),  
        axis.line.y = element_line(color="black", size = .6),
        axis.ticks.x = element_line(size = 0.5),
        axis.ticks.y = element_line(size = 0.5),
        axis.ticks.length = unit(0.2, "cm")) 



### read loadings
loading <- read.table("../output/filtered_Z.tsv", header = TRUE)

v142 <- loading[order(-loading$V142),]
write.csv(v142, "../outputv142_sorted.txt", row.names = TRUE)
top20 <- v142[1:20, 142, drop = FALSE]
top20 <- cbind(gene = rownames(top20), top20)
rownames(top20) <- NULL
top20$gene <- factor(top20$gene, levels = top20$gene[order(top20$V142, decreasing = TRUE)])
ggplot(data = top20, aes(x=gene, y=V142)) + 
  geom_bar(stat="identity") +
  ggtitle("V142 loadings")+
  theme_classic()+
  theme(plot.title = element_text(size = 24, hjust = 0.5), line = element_blank(), plot.background = element_blank(), panel.grid.major = element_blank()) + #set the background
  theme(panel.border = element_blank()) +   #set the border
  theme(axis.text.x = element_text()) + 
  theme(axis.title = element_text(size = 24), axis.text.x = element_text(size=24, angle=45, hjust = 1), axis.text.y = element_text(size=24)) + #set the x and y lab
  ylab("V142 loadings") + xlab("") +     #set the name of x-axis and y-axis
  theme( legend.title = element_blank(),  legend.position = "none", legend.text = element_text(size = 12), legend.key.width = unit(0.5, "cm")) + #set legend
  theme(axis.line.x = element_line(color="black", size = .6),  
        axis.line.y = element_line(color="black", size = .6),
        axis.ticks.x = element_line(size = 0.5),
        axis.ticks.y = element_line(size = 0.5),
        axis.ticks.length = unit(0.2, "cm"))


### LV99
v99 <- loading[order(-loading$V99),]
write.csv(v99, "../outputv99_sorted.txt", row.names = TRUE)
top20 <- v99[1:20, 99, drop = FALSE]
top20 <- cbind(gene = rownames(top20), top20)
rownames(top20) <- NULL
top20$gene <- factor(top20$gene, levels = top20$gene[order(top20$V99, decreasing = TRUE)])
ggplot(data = top20, aes(x=gene, y=V99)) + 
  geom_bar(stat="identity") +
  ggtitle("WT_micrglia, AD_microglia, and AD_astrocyte")+
  theme_classic()+
  theme(plot.title = element_text(size = 24, hjust = 0.5), line = element_blank(), plot.background = element_blank(), panel.grid.major = element_blank()) + #set the background
  theme(panel.border = element_blank()) +   #set the border
  theme(axis.text.x = element_text()) + 
  theme(axis.title = element_text(size = 24), axis.text.x = element_text(size=24, angle=45, hjust = 1), axis.text.y = element_text(size=24)) + #set the x and y lab
  ylab("LV99 loadings") + xlab("") +     #set the name of x-axis and y-axis
  theme( legend.title = element_blank(),  legend.position = "none", legend.text = element_text(size = 12), legend.key.width = unit(0.5, "cm")) + #set legend
  theme(axis.line.x = element_line(color="black", size = .6),  
        axis.line.y = element_line(color="black", size = .6),
        axis.ticks.x = element_line(size = 0.5),
        axis.ticks.y = element_line(size = 0.5),
        axis.ticks.length = unit(0.2, "cm"))



### LV105
v105 <- loading[order(-loading$V105),]
write.csv(v105, "../outputv105_sorted.txt", row.names = TRUE)
top20 <- v105[1:20, 105, drop = FALSE]
top20 <- cbind(gene = rownames(top20), top20)
rownames(top20) <- NULL
top20$gene <- factor(top20$gene, levels = top20$gene[order(top20$V105, decreasing = TRUE)])
ggplot(data = top20, aes(x=gene, y=V105)) + 
  geom_bar(stat="identity") +
  ggtitle("LV105 loadings")+
  theme_classic()+
  theme(plot.title = element_text(size = 24, hjust = 0.5), line = element_blank(), plot.background = element_blank(), panel.grid.major = element_blank()) + #set the background
  theme(panel.border = element_blank()) +   #set the border
  theme(axis.text.x = element_text()) + 
  theme(axis.title = element_text(size = 24), axis.text.x = element_text(size=24, angle=45, hjust = 1), axis.text.y = element_text(size=24)) + #set the x and y lab
  ylab("LV105 loadings") + xlab("") +     #set the name of x-axis and y-axis
  theme( legend.title = element_blank(),  legend.position = "none", legend.text = element_text(size = 12), legend.key.width = unit(0.5, "cm")) + #set legend
  theme(axis.line.x = element_line(color="black", size = .6),  
        axis.line.y = element_line(color="black", size = .6),
        axis.ticks.x = element_line(size = 0.5),
        axis.ticks.y = element_line(size = 0.5),
        axis.ticks.length = unit(0.2, "cm"))


### LV74
v74 <- loading[order(-loading$V74),]
top20 <- v74[1:20, 74, drop = FALSE]
top20 <- cbind(gene = rownames(top20), top20)
rownames(top20) <- NULL
top20$gene <- factor(top20$gene, levels = top20$gene[order(top20$V74, decreasing = TRUE)])
ggplot(data = top20, aes(x=gene, y=V74)) + 
  geom_bar(stat="identity") +
  ggtitle("LV74 loadings")+
  theme_classic()+
  theme(plot.title = element_text(size = 24, hjust = 0.5), line = element_blank(), plot.background = element_blank(), panel.grid.major = element_blank()) + #set the background
  theme(panel.border = element_blank()) +   #set the border
  theme(axis.text.x = element_text()) + 
  theme(axis.title = element_text(size = 24), axis.text.x = element_text(size=24, angle=45, hjust = 1), axis.text.y = element_text(size=24)) + #set the x and y lab
  ylab("LV74 loadings") + xlab("") +     #set the name of x-axis and y-axis
  theme( legend.title = element_blank(),  legend.position = "none", legend.text = element_text(size = 12), legend.key.width = unit(0.5, "cm")) + #set legend
  theme(axis.line.x = element_line(color="black", size = .6),  
        axis.line.y = element_line(color="black", size = .6),
        axis.ticks.x = element_line(size = 0.5),
        axis.ticks.y = element_line(size = 0.5),
        axis.ticks.length = unit(0.2, "cm"))

### LV41
v41 <- loading[order(-loading$V41),]
top20 <- v41[1:20, 41, drop = FALSE]
top20 <- cbind(gene = rownames(top20), top20)
rownames(top20) <- NULL
top20$gene <- factor(top20$gene, levels = top20$gene[order(top20$V41, decreasing = TRUE)])
ggplot(data = top20, aes(x=gene, y=V41)) + 
  geom_bar(stat="identity") +
  ggtitle("LV41 loadings")+
  theme_classic()+
  theme(plot.title = element_text(size = 24, hjust = 0.5), line = element_blank(), plot.background = element_blank(), panel.grid.major = element_blank()) + #set the background
  theme(panel.border = element_blank()) +   #set the border
  theme(axis.text.x = element_text()) + 
  theme(axis.title = element_text(size = 24), axis.text.x = element_text(size=24, angle=45, hjust = 1), axis.text.y = element_text(size=24)) + #set the x and y lab
  ylab("LV41 loadings") + xlab("") +     #set the name of x-axis and y-axis
  theme( legend.title = element_blank(),  legend.position = "none", legend.text = element_text(size = 12), legend.key.width = unit(0.5, "cm")) + #set legend
  theme(axis.line.x = element_line(color="black", size = .6),  
        axis.line.y = element_line(color="black", size = .6),
        axis.ticks.x = element_line(size = 0.5),
        axis.ticks.y = element_line(size = 0.5),
        axis.ticks.length = unit(0.2, "cm"))


