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

LVs <- read.table("../data/original_cocaine_lvs_reformatted.tsv", header = TRUE, colClasses = c("factor", "factor", "factor", "factor", "numeric"))
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
        outfile <- paste("../output/", lv, d, r, sep="_")
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
write.csv(df, "original_cocaine_LVs_pvalues.txt", row.names = FALSE)
sig_df <- subset(df, df$adjusted_p < 0.05)
print(sig_df)
