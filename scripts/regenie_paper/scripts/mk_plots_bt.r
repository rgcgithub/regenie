#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  if(!require(optparse)){ install.packages("optparse"); library(optparse) }
  if(!require(data.table)){ install.packages("data.table"); library(data.table) }
  if(!require(dplyr)){ install.packages("dplyr"); library(dplyr) }
  if(!require(tidyr)){ install.packages("tidyr"); library(tidyr) }
  if(!require(ggplot2)){ install.packages("ggplot2"); library(ggplot2) }
  if(!require(scales)){ install.packages("scales"); library(scales) }
  if(!require(viridis)){ install.packages("viridis"); library(viridis) }
  if(!require(gridExtra)){ install.packages("gridExtra"); library(gridExtra) }
  if(!require(cowplot)){ install.packages("cowplot"); library(cowplot) }
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
plot.folder <- opt$figfolder
manP.colors <- fread( opt$manColors, header=FALSE) %>% pull
pheno.names <- fread(opt$phenoNames, header = FALSE) %>% pull
log10P.X <- expression(paste("SAIGE -log"[10],"P"))
log10P.Y <- expression(paste("-log"[10],"P"))
################################################

## Manhatthan plots
pdf(NULL, width = 36, height=24)
dev.control(displaylist="enable")

layout(matrix(1:16,4, byrow=T))
par(mar = c(6,5,8,0.1), oma = c(1.1, 1.1, 1.1, 0.1))
par(tcl = -0.25)
par(mgp = c(3, 0.6, 0))

dummy <- lapply(1:length(pheno.names), function(ipheno){
  panel.lab <- letters[ipheno]
  pn <- pheno.names[ipheno]

  ## SAIGE results
  fns <- system(paste0("ls ", opt$prefix, "_saige_phenoCol", ipheno,"_chr*\\.test" ), intern = TRUE)
  saige.df <- rbindlist(lapply(fns, function(fn){
      fread(fn, fill=TRUE, showProgress = FALSE) %>%
        mutate( pval.SAIGE = -log10(Tstat) ,
          SE.SAIGE = BETA, BETA.SAIGE = N,
          A1FREQ = as.numeric(AC_Allele2),
          MAF = pmin(A1FREQ, 1 - A1FREQ) ) %>%
      select( rsid, MAF, BETA.SAIGE, SE.SAIGE, pval.SAIGE )
}))

  ## BOLT results
  fns <- system(paste0("ls ", opt$prefix, "_bolt_phenoCol", ipheno,"\\.test" ), intern = TRUE)
  bolt.df <- rbindlist(lapply(fns, function(fn){
      fread(fn, colClasses="character", showProgress = FALSE) %>%
        mutate( pval.BOLT = -log10(as.numeric(P_BOLT_LMM)) ) %>%
        select( SNP, pval.BOLT )
}))

  # regenie-Firth
  fns <- system(paste0("ls ", opt$prefix, "_regenie_phenoCol", ipheno,"_Firth_chr*\\.regenie" ), intern = TRUE)
  regenie.firth.df <- rbindlist(lapply(fns, function(fn){
      m1 <-  fread(fn, showProgress = FALSE) %>%
        select( CHROM, GENPOS, ID, BETA, SE, CHISQ,LOG10P ) %>%
        rename_(.dots = setNames(c("BETA","SE","CHISQ","LOG10P"), c("beta.rg.firth","se.rg.firth","test.firth","pval.regenie.firth")) )
      m1 %>%
        select( CHROM, GENPOS, ID, beta.rg.firth, se.rg.firth, pval.regenie.firth)
}))

  # regenie-SPA
  fns <- system(paste0("ls ", opt$prefix, "_regenie_phenoCol", ipheno,"_SPA_chr*\\.regenie" ), intern = TRUE)
  regenie.spa.df <- rbindlist(lapply(fns, function(fn){
      m1 <-  fread(fn, showProgress = FALSE) %>%
        select( ID, BETA, SE, LOG10P ) %>%
        rename_(.dots = setNames(c("BETA","SE","LOG10P"), c("beta.rg.spa","se.rg.spa","pval.regenie.spa")) )
}))

  #     # combine info
  comb.df <- left_join(saige.df, regenie.firth.df, by = c( "rsid" = "ID") )  %>%
    left_join(., bolt.df, by = c( "rsid" = "SNP") )  %>%
    left_join(., regenie.spa.df, by = c( "rsid" = "ID") ) %>%
    as.data.table
  cat("# variants tested for", pn, "=", scales::number(comb.df %>% drop_na %>% nrow, big.mark=",") ,"\n")

  # convert infinite values to p.thr
  comb.df[ is.infinite(pval.BOLT), `:=`(pval.BOLT=p.thr)]
  comb.df[ is.infinite(pval.SAIGE), `:=`(pval.SAIGE=p.thr)]
  comb.df[ pval.regenie.firth > p.thr, `:=`(pval.regenie.firth = p.thr)]
  comb.df[ pval.regenie.spa > p.thr, `:=`(pval.regenie.spa = p.thr)]
  comb.df <- comb.df %>% drop_na

  man.data <- comb.df  %>%
    rename(POS = GENPOS, CHR = CHROM) %>%
    prep.data
  axisdf <- man.data %>% group_by(CHR) %>% summarize(center=( max(POScum) + min(POScum) ) / 2 )

  ## using base R
  all.methods <- c("BOLT-LMM","REGENIE-FIRTH","REGENIE-SPA","SAIGE")
  df0 <- man.data %>%
    gather(Method, pval, pval.BOLT, pval.SAIGE, pval.regenie.firth, pval.regenie.spa) %>%
    mutate( Method = factor(Method,
        levels = paste0("pval.",c("BOLT","regenie.firth", "regenie.spa","SAIGE")),
        labels= all.methods) ,
      pval = as.numeric(pval),
      pval = pval * (pval < p.thr) + p.thr * (pval >= p.thr ),
      color = manP.colors[CHR]
      ) %>%
    select(POScum,Method,pval,color)
  chrlims.min <- man.data %>% group_by(CHR) %>% summarize(bot= min(POScum))%>%pull
  chrlims.max <- man.data %>% group_by(CHR) %>% summarize(top= max(POScum))%>%pull
  chrlims <- rowMeans(cbind(head(chrlims.max,-1), chrlims.min[-1]))
  myaxis <- c(0, ifelse(max(df0$pval) < 70 , 70, 120))


  # compress plot for pvalues >= 20 by 8/100 => x/12.5 + 20 - 1.6
  df1 <- df0  %>%
    mutate(
      pval = pval * (pval <= 20) + (pval/12.5+20-1.6) * (pval > 20 )
    )
  tot.vals <- ifelse(max(df0$pval) < 70 , 7, 8)
  myaxis <- c(seq(0, 20, by = 4), seq(24, 60, by = 4))[1:tot.vals]
  myaxis.lab <- c(seq(0,20,by=4),70,120)[1:tot.vals]

  ddd <- lapply(all.methods, function(my.method){
    cat(my.method,"\n")

    df1.meth <- df1 %>% filter( Method == my.method ) %>%
      mutate(x=POScum,y=pval) %>%
      select(x,y,color)

    plot(df1.meth$x, df1.meth$y, axes = FALSE,
      col=df1.meth$color, cex=1, pch=16, cex.lab=1.7,
      xlab="", ylab="",
      main="", ylim=c(0,max(myaxis)))
    abline(h=-log10(5e-8),lty=2); abline(h=20,col="gray45",lty=2)

    # method name
    title(my.method, cex.main=3.5, font.main=1, line=0)

    if(my.method==all.methods[1]) {
      # add trait title, x/y-axis label
      fig_label(panel.lab, cex=4, font=2)
      title(main=pn, adj=0, line = 4.5, cex.main=3, font.main=1)
      # y-axis
      mtext(log10P.Y, side = 2, outer = F, cex = 1.6, line = 3.5)
      axis(2, at = myaxis, labels=myaxis.lab, las=1,, cex.axis=2)
      plotrix::axis.break(axis=2,breakpos=20,bgcol="white",breakcol="black", style="zigzag",brw=0.02)
    }

    # x-axis labels (jittered)
    abline(v=chrlims, lty=2, col="gray65") #chr breaks
    axis(1, at = c(0,chrlims, tail(chrlims.max,1)), labels=FALSE) #chr ticks
    # x-axes jittered
    x.c <- seq(1,22,2)
    text(x=axisdf$center[x.c],  par("usr")[3],
      labels = axisdf$CHR[x.c], pos = 1, xpd = TRUE, cex=2.1, offset=.7)
    text(x=axisdf$center[-x.c],  par("usr")[3],
      labels = axisdf$CHR[-x.c], pos = 1, xpd = TRUE, cex=2.1, offset=2.7)
    return(NULL)
    })
  return(NULL)
})

p1.base <- recordPlot()
invisible(dev.off())

## TIFF format
tiff(paste0(plot.folder, "Figure2.tiff"), width= 300*6*1.5*4 , height = 300*6*4, res = 300, compression ="lzw")
p1.base
dev.off()
## convert to pdf using imagemagick
system(paste0("convert -density 300 -units PixelsPerInch ", plot.folder, "Figure2.tiff ", plot.folder, "Figure2.pdf"))

