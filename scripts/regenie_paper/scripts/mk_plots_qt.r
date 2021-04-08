#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  if(!require(optparse)){ install.packages("optparse"); library(optparse) }
  if(!require(data.table)){ install.packages("data.table"); library(data.table) }
  if(!require(tidyverse)){ install.packages("tidyverse"); library(tidyverse) }
  if(!require(ggplot2)){ install.packages("ggplot2"); library(ggplot2) }
  if(!require(gridExtra)){ install.packages("gridExtra"); library(gridExtra) }
  if(!require(cowplot)){ install.packages("cowplot"); library(cowplot) }
  if(!require(extrafont)){ install.packages("extrafont"); library(extrafont) }
})
#########################################
##
## Script used to run REGENIE/BOLT/fastGWA/SAIGE for GWAS
## for the analyses in the REGENIE 2020 paper
## For more details, visit: https://rgcgithub.github.io/regenie/
##  
#########################################
option_list = list(
  make_option("--loadFuns", type="character", default="",
    help="script with functions for plotting"),
  make_option("--manColors", type="character", default="",
    help="file with colors for Manhattan plot"),
  make_option("--phenoNames", type="character", default="",
    help="phenotype file"),
  make_option("--figfolder", type="character", default="",
    help="prefix of output files"),
  make_option("--prefix", type="character", default="",
    help="prefix of output files")
  );
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

source(opt$loadFuns)
alpha <- .05
p.thr <- 324
log10P.Y <- expression(paste("-log"[10],"P"))
plot.folder <- opt$figfolder
manP.colors <- fread( opt$manColors, header=FALSE) %>% pull
pheno.names <- fread(opt$phenoNames, header = FALSE) %>% pull
################################################


 ## Manhatthan plots
## 3x3 plot
pdf(NULL)
dev.control(displaylist="enable")
layout(matrix(1:9,3, byrow=T))
par(mar = c(6,5,6,0.1), oma = c(1.1, 1.1, 1.1, 0.1))
par(tcl = -0.25)
par(mgp = c(3, 0.6, 0))

dummy <- lapply(1:length(pheno.names), function(ipheno){
  panel.lab <- letters[ipheno]
  pn <- pheno.names[ipheno]
  pn.title <- c("LDL", "BMI", "Bilirubin")
  pn.title <- pn.title[ sapply(pn.title, function(x) grepl(x, pn, ignore.case=TRUE)) ]

  ## BOLT results
  fns <- system(paste0("ls ", opt$prefix, "_bolt_phenoCol", ipheno,"\\.test" ), intern = TRUE)
  bolt.df <- rbindlist(lapply(fns, function(fn){
      fread(fn, colClasses="character", showProgress = FALSE) %>%
        mutate(pval.BOLT = -log10(as.numeric(P_BOLT_LMM)) ,
          A1FREQ = as.numeric(A1FREQ),
          MAF = pmin(A1FREQ, 1 - A1FREQ) ) %>%
      dplyr::filter( MAF >= 0.01 ) %>%
      select( SNP, pval.BOLT )
}))

  ## fastgwa results
  fns <- system(paste0("ls ", opt$prefix, "_fastgwa_phenoCol", ipheno,"\\.test.fastGWA" ), intern = TRUE)
  fastGWA.df <- rbindlist(lapply(fns, function(fn){
      m1 <- fread(fn, showProgress = FALSE) %>%
        mutate( pval.fastGWA = -log10(P))
      # set the p-values of 0 to minimum p-value for plitting
      m1$pval.fastGWA[ m1$ P == 0 ]  <-  p.thr
      m1 %>% select( SNP, pval.fastGWA )
}))

  ## regenie results
  fns <- system(paste0("ls ", opt$prefix, "_regenie_phenoCol", ipheno,"_chr*\\.regenie" ), intern = TRUE)
  regenie.df <- rbindlist(lapply(fns, function(fn){
      fread(fn) %>%
        select( CHROM, GENPOS, ID, LOG10P ) %>%
        setNames(c("CHROM","GENPOS","ID","pval.regenie") )
}))

  ## combine info
  comb.df <- left_join(bolt.df, fastGWA.df, by = c( "SNP" = "SNP") ) %>%
    left_join(., regenie.df, by = c( "SNP" = "ID") ) %>%
    as.data.table

  # convert infinite values to p.thr
  comb.df[ is.infinite(pval.BOLT), `:=`(pval.BOLT=p.thr)]
  comb.df[ pval.regenie > p.thr, `:=`(pval.regenie = p.thr)]
  comb.df <- comb.df %>% drop_na

  man.data <- comb.df  %>%
    rename(POS = GENPOS, CHR = CHROM) %>%
    prep.data
  axisdf <- man.data %>% group_by(CHR) %>% summarize(center=( max(POScum) + min(POScum) ) / 2 )
  chrlims <- man.data %>% group_by(CHR) %>% summarize(top= max(POScum))%>% head(-1)

  # using base R
  all.methods <- c("REGENIE" ,"fastGWA","BOLT-LMM")
  df0 <- man.data %>%
    gather(Method, pval, pval.BOLT, pval.fastGWA, pval.regenie) %>%
    mutate( Method = factor(Method,
        levels = paste0("pval.",c("regenie","fastGWA","BOLT")),
        labels = all.methods ),
      color = manP.colors[CHR]
      ) %>%
    select(-SNP,-POS,-tot,-CHR)
  chrlims.min <- man.data %>% group_by(CHR) %>% summarize(bot= min(POScum))%>%pull
  chrlims.max <- man.data %>% group_by(CHR) %>% summarize(top= max(POScum))%>%pull
  chrlims <- rowMeans(cbind(head(chrlims.max,-1), chrlims.min[-1]))
  myaxis <- c(0, ifelse(max(df0$pval < 300) , 300, 350))

  # compress plot for pvalues >= 20 by 1/20 => x/20 + 19
  df1 <- df0  %>%
    mutate(
      pval = pval * (pval <= 20) + (pval/20+19) * (pval > 20 )
    )
  tot.vals <- ifelse(max(df0$pval) < 300 , 9, 10)
  myaxis <- c(seq(0, 20, by = 4), seq(24, 60, by = 5))[1:tot.vals]
  myaxis.lab <- c(seq(0,20, by=4), seq(100,400,by=100))[1:tot.vals]

  plot.d <- lapply(all.methods, function(my.method){

    df1.meth <- df1 %>% dplyr::filter( Method == my.method ) %>%
      mutate(x=POScum,y=pval) %>%
      select(x,y,color)

    plot(df1.meth$x, df1.meth$y, axes = FALSE,
      col=alpha(df1.meth$color, .6), cex=.8, pch=16, cex.lab=1.7,
      xlab="", ylab="",
      main="", ylim=c(0,max(myaxis))
    )
    abline(h=-log10(5e-8),lty=2);abline(h=20,col="gray45",lty=2)

    # method name
    title(my.method, cex.main=3.2, font.main=1, line=1)

    # add trait title, x/y-axis label
    if(my.method == all.methods[1]) {
      fig_label(panel.lab, cex=4, font=2)
      title(main=pn.title, adj=0, line = 2.5, cex.main=3.5, font.main=1)
      # y-axis
      mtext(log10P.Y, side = 2, outer = F, cex = 1.5, line = 3.5)
      axis(2, at = myaxis, labels=myaxis.lab, las=1, cex.axis=2)
      plotrix::axis.break(axis=2,breakpos=20,bgcol="white",breakcol="black", style="zigzag",brw=0.02)
    } 

    # x-axis labels (jittered)
    abline(v=chrlims, lty=2, col="gray65") #chr breaks
    axis(1, at = c(0,chrlims, tail(chrlims.max,1)), labels=FALSE) #chr ticks
    x.c <- seq(1,22,2)
    text(x=axisdf$center[x.c],  par("usr")[3],
      labels = axisdf$CHR[x.c], pos = 1, xpd = TRUE, cex=2.1, offset=.7)
    text(x=axisdf$center[-x.c],  par("usr")[3],
      labels = axisdf$CHR[-x.c], pos = 1, xpd = TRUE, cex=2.1, offset=2.5)

    return(NULL)
    })

  return(NULL)
  })

p1.base <- recordPlot()
invisible(dev.off())

## TIFF format
tiff(paste0(plot.folder, "Figure1.tiff"), width= 300*6*1.5*3 , height = 300*6*3, res = 300, compression ="lzw")
p1.base
dev.off()

## convert to pdf using imagemagick
system(paste0("convert -density 300 -units PixelsPerInch ", plot.folder, "Figure1.tiff ", plot.folder, "Figure1.pdf"))
