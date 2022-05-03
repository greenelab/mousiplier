library(ggplot2)
library(dplyr)

if (rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
} else{
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

LVs <- read.table("reformated_NAc_PFC_VTA_Lvs.txt", header = TRUE, colClasses = c("factor", "factor", "factor", "factor", "numeric"))
#summary(LVs)

df <- data.frame(matrix(ncol = 4, nrow = 0))

colnames(df) <- c("LV_ID", "day", "region", "pvalue")


### perform one way anova
for (i in 1: length(levels(LVs$LV_ID))) {
  lv <- paste("LV", as.character(i), sep="")

  for (d in levels(LVs$day)) {

    for (r in levels(LVs$region)) {

      temp_lv <- subset(LVs, LVs$LV_ID == lv & LVs$day == d & LVs$region == r)

      ## one way anova
      one.way <- aov(lv_value ~ treatment, data = temp_lv)

      #summary(one.way)
      pvalue <- summary(one.way)[[1]][[1,"Pr(>F)"]]


      df[nrow(df)+1, ] <- c(lv, d, r, pvalue)

      ### plot all LVs whose pvalues < 0.05
      if (pvalue < 0.05) {
        outfile <- paste("./sig_lvs_figures/", lv, d, r, sep="_")
        pdf(file = outfile,   # The directory you want to save the file in
            width = 4, # The width of the plot in inches
            height = 4)
        myplot <- ggplot(temp_lv, aes(x=treatment, y=lv_value, fill=treatment)) +
          geom_boxplot()
        print(myplot)
        dev.off()
      }
    }
  }
}

### BH correction and write the results
df$adjusted_p <- p.adjust(df$pvalue, "BH")
write.csv(df, "LVs_pvalues.txt", row.names = FALSE)
sig_df <- subset(df, df$adjusted_p < 0.05)

### LV135
LV135_day1_NAc <- subset(LVs, LVs$LV_ID == "LV135" & LVs$day == "day1" & LVs$region == "NAc")
LV135_day1_NAc$treatment <- factor(LV135_day1_NAc$treatment, levels = c("saline", "food", "cocaine"))
ggplot(LV135_day1_NAc, aes(x=treatment, y=lv_value, fill=treatment)) +
  geom_boxplot() +
  ggtitle("Nucleus accumbens day 1 abstinence")+
  theme_classic()+
  theme(plot.title = element_text(size = 20, hjust = 0.5), line = element_blank(), plot.background = element_blank(), panel.grid.major = element_blank()) + #set the background
  theme(panel.border = element_blank()) +   #set the border
  theme(axis.text.x = element_text()) +
  theme(axis.title = element_text(size = 20), axis.text.x = element_text(size=20, angle=0, hjust = 0.5), axis.text.y = element_text(size=20)) + #set the x and y lab
  ylab("LV135") + xlab("") +     #set the name of x-axis and y-axis
  theme( legend.title = element_blank(),  legend.position = "none", legend.text = element_text(size = 12), legend.key.width = unit(0.5, "cm")) + #set legend
  theme(axis.line.x = element_line(color="black", size = .6),
        axis.line.y = element_line(color="black", size = .6),
        axis.ticks.x = element_line(size = 0.5),
        axis.ticks.y = element_line(size = 0.5),
        axis.ticks.length = unit(0.2, "cm"))


wilcox.test(subset(LV135_day1_NAc, treatment == "cocaine")$lv_value, subset(LV135_day1_NAc, treatment == "saline")$lv_value)

### End of LV135


### read loadings
loading <- read.table("../output/Z.tsv", header = TRUE)

v135 <- loading[order(-loading$V135),]
write.csv(v135, "v135_sorted.txt", row.names = TRUE)

v9 <- loading[order(-loading$V9),]
write.csv(v9, "v9_sorted.txt", row.names = TRUE)

v142 <- loading[order(-loading$V142),]
write.csv(v142, "v142_sorted.txt", row.names = TRUE)


LV9_day28_NAc <- subset(LVs, LVs$LV_ID == "LV9" & LVs$day == "day28" & LVs$region == "NAc")
LV9_day28_NAc$treatment <- factor(LV9_day28_NAc$treatment, levels = c("saline", "cocaine"))
ggplot(LV9_day28_NAc, aes(x=treatment, y=lv_value, fill=treatment)) +
  geom_boxplot() +
  ggtitle("Nucleus accumbens")+
  theme_classic()+
  theme(plot.title = element_text(size = 24, hjust = 0.5), line = element_blank(), plot.background = element_blank(), panel.grid.major = element_blank()) + #set the background
  theme(panel.border = element_blank()) +   #set the border
  theme(axis.text.x = element_text()) +
  theme(axis.title = element_text(size = 24), axis.text.x = element_text(size=24, angle=0, hjust = 0.5), axis.text.y = element_text(size=24)) + #set the x and y lab
  ylab("LV9") + xlab("") +     #set the name of x-axis and y-axis
  theme( legend.title = element_blank(),  legend.position = "none", legend.text = element_text(size = 12), legend.key.width = unit(0.5, "cm")) + #set legend
  theme(axis.line.x = element_line(color="black", size = .6),
        axis.line.y = element_line(color="black", size = .6),
        axis.ticks.x = element_line(size = 0.5),
        axis.ticks.y = element_line(size = 0.5),
        axis.ticks.length = unit(0.2, "cm"))

LV9_day28_PFC <- subset(LVs, LVs$LV_ID == "LV9" & LVs$day == "day28" & LVs$region == "PFC")
LV9_day28_PFC$treatment <- factor(LV9_day28_PFC$treatment, levels = c("saline", "cocaine"))
ggplot(LV9_day28_PFC, aes(x=treatment, y=lv_value, fill=treatment)) +
  geom_boxplot() +
  ggtitle("Prefrontal cortex")+
  theme_classic()+
  theme(plot.title = element_text(size = 24, hjust = 0.5), line = element_blank(), plot.background = element_blank(), panel.grid.major = element_blank()) + #set the background
  theme(panel.border = element_blank()) +   #set the border
  theme(axis.text.x = element_text()) +
  theme(axis.title = element_text(size = 24), axis.text.x = element_text(size=24, angle=0, hjust = 0.5), axis.text.y = element_text(size=24)) + #set the x and y lab
  ylab("LV9") + xlab("") +     #set the name of x-axis and y-axis
  theme( legend.title = element_blank(),  legend.position = "none", legend.text = element_text(size = 12), legend.key.width = unit(0.5, "cm")) + #set legend
  theme(axis.line.x = element_line(color="black", size = .6),
        axis.line.y = element_line(color="black", size = .6),
        axis.ticks.x = element_line(size = 0.5),
        axis.ticks.y = element_line(size = 0.5),
        axis.ticks.length = unit(0.2, "cm"))

LV9_day28_VTA <- subset(LVs, LVs$LV_ID == "LV9" & LVs$day == "day28" & LVs$region == "VTA")
LV9_day28_VTA$treatment <- factor(LV9_day28_VTA$treatment, levels = c("saline", "cocaine"))
ggplot(LV9_day28_VTA, aes(x=treatment, y=lv_value, fill=treatment)) +
  geom_boxplot() +
  ggtitle("Ventral tegmental area")+
  theme_classic()+
  theme(plot.title = element_text(size = 24, hjust = 0.5), line = element_blank(), plot.background = element_blank(), panel.grid.major = element_blank()) + #set the background
  theme(panel.border = element_blank()) +   #set the border
  theme(axis.text.x = element_text()) +
  theme(axis.title = element_text(size = 24), axis.text.x = element_text(size=24, angle=0, hjust = 0.5), axis.text.y = element_text(size=24)) + #set the x and y lab
  ylab("LV9") + xlab("") +     #set the name of x-axis and y-axis
  theme( legend.title = element_blank(),  legend.position = "none", legend.text = element_text(size = 12), legend.key.width = unit(0.5, "cm")) + #set legend
  theme(axis.line.x = element_line(color="black", size = .6),
        axis.line.y = element_line(color="black", size = .6),
        axis.ticks.x = element_line(size = 0.5),
        axis.ticks.y = element_line(size = 0.5),
        axis.ticks.length = unit(0.2, "cm"))


### plot highest top ranked genes from latent variable 9
v9_top <- read.table("v9_topgene.txt", sep=',', header = TRUE)
v9_top$geneID <- factor(v9_top$geneID, levels = c('Eif1ad', 'Khsrp', 'Cdc23', 'Pbx2', 'Rnps1', 'Ppp6r3', 'Sart1', 'Sf1', 'Taf11', 'Ercc3', 'Dpp9', 'Ppp2r5d', 'Traf7', 'Ganab', 'Rbm22', 'Ranbp3', 'Prpf19', 'Men1', 'Ddx39b', 'Sf3b2'))
ggplot(data = v9_top, aes(x=geneID, y=V9)) +
  geom_bar(stat="identity") +
  ggtitle("V9 loadings")+
  theme_classic()+
  theme(plot.title = element_text(size = 24, hjust = 0.5), line = element_blank(), plot.background = element_blank(), panel.grid.major = element_blank()) + #set the background
  theme(panel.border = element_blank()) +   #set the border
  theme(axis.text.x = element_text()) +
  theme(axis.title = element_text(size = 24), axis.text.x = element_text(size=24, angle=45, hjust = 1), axis.text.y = element_text(size=24)) + #set the x and y lab
  ylab("V9 loadings") + xlab("") +     #set the name of x-axis and y-axis
  theme( legend.title = element_blank(),  legend.position = "none", legend.text = element_text(size = 12), legend.key.width = unit(0.5, "cm")) + #set legend
  theme(axis.line.x = element_line(color="black", size = .6),
        axis.line.y = element_line(color="black", size = .6),
        axis.ticks.x = element_line(size = 0.5),
        axis.ticks.y = element_line(size = 0.5),
        axis.ticks.length = unit(0.2, "cm"))
