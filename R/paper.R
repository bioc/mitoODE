loadFittedData <- function() {
  ## to avoid warnings during R CMD check
  if (FALSE) {
    pheno <- NULL
    zqc <- NULL
    tab2 <- NULL
    g.kim <- NULL
    g.kmi <- NULL
    g.mit0 <- NULL
    p.lambda <- NULL
  }

  ## load Mitocheck tab and anno
  mitoODEdata:::loadMitocheck()
  
  ## use pconst
  assign("g.kim", 0.025, 1)
  assign("g.kmi", 0.57, 1)
  assign("g.mit0", 0.05, 1)
  assign("p.lambda", 4, 1)
  
  ## load assay fitted phenotypes and compute interphase and mitosis duration
  mitoODEdata:::loadPheno()
  pheno <- cbind(pheno, mitod=0.6/(abs(g.kmi + pheno[,"hmi"] + pheno[,"hmp"] +  pheno[,"ha"])+0.02))
  pheno <- cbind(pheno, interd=1/(abs(g.kim + pheno[,"him"] + pheno[,"ha"])+0.02))
  assign("pheno", pheno, 1)

  ## zqc and tab2
  assign("zqc", buildQC(), 1)
  assign("tab2", buildSuppTab(), 1)
}

buildQC <- function() {
  ## to avoid warnings in R CMD check
  if (FALSE) pheno <- NULL
  
  qscore <- quantile(pheno[,"score"], 0.95, na.rm=TRUE)
  zqc <- tab$qc==TRUE & pheno[,"score"]<qscore & pheno[,"mitod"]>(20/60)
  zqc
}

compute.mtvt <- function(parameter, z, tab2, max.mean=50, max.sd=4) {
  ## to avoid warnings in R CMD check
  if (FALSE) {
    pheno <- NULL
    zqc <- NULL
  }
  
  ssirna <- split(which(z), tab$sirna[z])
  mc.cores <- detectCores()
  mt <- unlist(mclapply(ssirna, function(z) mean((pheno[z, parameter])), mc.cores=mc.cores))
  st <- unlist(mclapply(ssirna, function(z) if (length(z)<=1) NA else sd((pheno[z, parameter])), mc.cores=mc.cores))
  w <- names(which(mt<max.mean & st<max.sd))
  res <- setNames(rep(NA, nrow(tab2)), rownames(tab2))
  res[w] <- mt[w]
  res
}

buildSuppTab <- function() {
  ## to avoid warnings in R CMD check
  if (FALSE) {
    pheno <- NULL
    zqc <- NULL
  }
  
  ## init tab2
  usirna <- sort(unique(tab$sirna))
  tab2 <- data.frame(target.hgnc=getanno(sirna=usirna), stringsAsFactors=FALSE)
  rownames(tab2) <- usirna
  mc.cores <- detectCores()
  
  ## time.quiescence
  z <- zqc & tab$type=="experiment" & pheno[,"him"]<(-0.025) & pheno[,"mu"]<0.5
  tab2$time.quiescence <- compute.mtvt("tim", z, tab2)

  ## time.mitotic.arrest
  z <- zqc & tab$type=="experiment" & pheno[,"hmi"]<(-0.5)  & pheno[,"mu"]<0.5
  tab2$time.mitotic.arrest <- compute.mtvt("tmi", z, tab2)

  ## time.polynucleation
  z <- zqc & tab$type=="experiment" & pheno[,"hmp"]>(0.25)  & pheno[,"mu"]<0.5
  tab2$time.polynucleation <- compute.mtvt("tmp", z, tab2)
 
  ## time.death
  z <- zqc & tab$type=="experiment" & pheno[,"ha"]>0.005  & pheno[,"mu"]<0.5
  tab2$time.death <- compute.mtvt("ta", z, tab2)

  ## duration.mitosis
  z <- zqc & tab$type=="experiment" & pheno[, "tmi"] < 50  & pheno[,"mu"]<0.5
  ssirna <- split(which(z), tab$sirna[z])
  md <- unlist(mclapply(ssirna, function(z) 2^mean(log2(pheno[z, "mitod"])), mc.cores=mc.cores))
  vd <- unlist(mclapply(ssirna, function(z) if (length(z)<=1) NA else 2^sd(log2(pheno[z, "mitod"])), mc.cores=mc.cores))
  sd <- names(which(md>2 & vd<2))
  tab2$duration.mitosis <- NA
  tab2[sd, "duration.mitosis"] <- md[sd]
  
  ## duration.interphase
  z <- zqc & tab$type=="experiment" & pheno[, "tim"] < 50  & pheno[,"mu"]<0.5
  ssirna <- split(which(z), tab$sirna[z])
  md <- unlist(mclapply(ssirna, function(z) 2^mean(log2(pheno[z, "interd"])), mc.cores=mc.cores))
  vd <- unlist(mclapply(ssirna, function(z) if (length(z)<=1) NA else 2^sd(log2(pheno[z, "interd"])), mc.cores=mc.cores))
  sd <- names(which(md>40 & vd<4))
  tab2$duration.interphase <- NA
  tab2[sd, "duration.interphase"] <- md[sd]

  ## keep only siRNAs with at least one non-NA value
  z <- match("target.hgnc", colnames(tab2))
  tab2 <- tab2[rowSums(!is.na(tab2[, -z]))>0,]
  tab2 <- tab2[!is.na(tab2$target.hgnc),]
  sum(tab2[,-z], na.rm=TRUE) ## 51892.28

  ## save formatted df
  df <- cbind(sirna=rownames(tab2), target.hgnc=tab2$target.hgnc, round(tab2[,-z], 1), stringsAsFactors=FALSE)
  df[is.na(df)] <- ""
  write.table(df, "suppTab.tsv", sep="\t", row.names=FALSE, quote=FALSE)

  tab2
}

stats.assay <- function() {
  ## to avoid warnings during R CMD check
  if (FALSE) {
    pheno <- NULL
    zqc <- NULL
    tab2 <- NULL
  }
  
  ## nb spots
  nrow(tab) ## 206592
  
  ## nb target hgnc
  hgncs <- getanno(sirna=unique(tab$sirna))
  length(na.omit(unique(hgncs))) ## nb.hgncs 17293

  ## nb sirnas/hgnc
  z <- tapply(tab$sirna, getanno(sirna=tab$sirna), unique)
  mean(sapply(z, length)>=2) ## 98.7 % of hgnc have at 2 sirnas

  ## nb spot/sirnas
  z <- tapply(1:nrow(tab), getsirna(spot=1:nrow(tab)), unique)
  mean(sapply(z, length)>=3)
  
  ## nb different sirans
  length(unique(tab$sirna)) ## 51766
 
  ## assay layout
  layout <- do.call(rbind, tapply(tab$type, tab$pr, table))
  table(layout[,"scrambled"])
  table(layout[,"bcop"]+layout[,"incenp"]+layout[,"eg5"])

  ## contamination
  z <- tab$type=="experiment"
  z1 <- zqc &pheno[,"mu"]<0.5
  sum(z&z1)/sum(z)  ## fracion of spot experiment that passed QC
  
  ## tab2
  dim(tab2)
  colSums(!is.na(tab2))
}

stats.fitting <- function() {
  ## to avoid warnings during R CMD check
  if (FALSE) {
    pheno <- NULL
  }
  
  ## mre
  mc.cores <- detectCores()
  mres <- parallel::mclapply(sample(1:nrow(pheno), 5000), mc.cores=mc.cores, function(id) {
    y <-  readspot(id)
    yf <- odevaluate(pheno[id,], nt=nrow(y))
    mre <- abs(y-yf$y)/max(y)
    mean(mre)
  })
  quantile(unlist(mres), 0.95) ## 95 % of the spots have a MRE lower than 3.2 %

  ## r2
  r2s <- parallel::mclapply(sample(1:nrow(pheno), 5000), mc.cores=mc.cores, function(id) {
    y <-  readspot(id)
    yf <- odevaluate(pheno[id,], nt=nrow(y))
    sstot <- mean((y-mean(y))^2)
    sserr <- mean((y-yf$y)^2)
    1 - sserr/sstot
  })
  quantile(unlist(r2s), 1-0.95) ## 95 % of the spots have a R2 higher than 0.97
}

figure1 <- function() {
  ## graphics parameters
  width <- 7
  height <- 7
  lwd <- 2
  cex <- 1.5
  ide <- c(scrambled=699L, scrambled=52192L, eg5=156205L, eg5=66722L, bcop=92880)

  ## plot panels
  pdf("raws1.pdf", width=width, height=height)
  plotfit(ide[1], showfit=FALSE, lwd=lwd, cex.axis=cex, cex=cex, xlab="", ylim=c(0, 205), cex.lab=cex)
  dev.off()
  pdf("raws2.pdf", width=width, height=height)
  plotfit(ide[2], showfit=FALSE, lwd=lwd, cex.axis=cex, legend=NULL, xlab="", ylab="", ylim=c(0, 205), cex.lab=cex)
  dev.off()
  pdf("rawk1.pdf", width=width, height=height)
  plotfit(ide[3], showfit=FALSE, lwd=lwd, cex.axis=cex, legend=NULL, ylim=c(0, 55), cex.lab=cex)
  dev.off()
  pdf("rawc1.pdf", width=width, height=height)
  plotfit(ide[5], showfit=FALSE, lwd=lwd, cex.axis=cex, legend=NULL, ylab="", ylim=c(0, 205), cex.lab=cex)
  invisible(dev.off())
}

figure2 <- function() {
  ## graphics parameters
  width <- 7
  height <- 7
  heightk <- 3.5
  lwd <- 2
  cex <- 1.5
  kcol <- c("kim"="#555555", "ka"="#3333ff")
  ide <- c(scrambled=699L, scrambled=52192L, eg5=156205L, eg5=66722L, bcop=92880)
  
  ## penetrances
  pdf("pk1s1.pdf", width=width, height=heightk)
  plotk(ide[1], kk="kim", cex.axis=cex, cex.lab=cex, lwd=lwd, kcol=kcol, xlab="", ylim=c(0, 0.04), ylab="")
  dev.off()
  pdf("pk2s1.pdf", width=width, height=heightk)
  plotk(ide[1], kk="ka", cex.axis=cex, cex.lab=cex, lwd=lwd, kcol=kcol, ylim=c(0, 0.04), ylab="")
  dev.off()
  pdf("pk1c1.pdf", width=width, height=heightk)
  plotk(ide[5], kk="kim", cex.axis=cex, cex.lab=cex, lwd=lwd, kcol=kcol, xlab="", ylim=c(0, 0.04), ylab="")
  dev.off()
  pdf("pk2c1.pdf", width=width, height=heightk)
  plotk(ide[5], kk="ka", cex.axis=cex, cex.lab=cex, lwd=lwd, kcol=kcol, ylim=c(0, 0.04), ylab="")
  dev.off()

  ## fitted examples
  pdf("fits1.pdf")
  plotfit(ide[1], lwd=lwd, kk="kim", cex.axis=cex, kcol=kcol, cex=cex, cex.lab=cex, ylim=c(0, 205))
  dev.off()
  pdf("fitc1.pdf")
  plotfit(ide[5], lwd=lwd, kk=c("kim", "ka"), cex.axis=cex, cex.lab=cex, legend=NULL, kcol=kcol, ylab="", ylim=c(0, 205))
  invisible(dev.off())
}

figure3a <- function() {
  ## to avoid warnings during R CMD check
  if (FALSE) {
    pheno <- NULL
    zqc <- NULL
  }
  
  ## graphics parameters
  cex <- 1.5
  lwd <- 2
  kcol <- c("kim"="#555555", "ka"="#3333ff")
  dc <- data.frame(row.names=c('scrambled', 'eg5', 'incenp', 'bcop', 'experiment'),
                   color=c('#77dd77', '#ffaa77', '#ff7777', '#aaaaff', '#aaaaaa'),
                   name=c( 'siScrambled', 'siKIF11', 'siINCENP', 'siCOPB1', 'sample'), stringsAsFactors=FALSE)
  
  ## cell death penetrance
  pdf("boxplotha.pdf", width=6, height=5)
  boxplot(split(pheno[zqc, "ha"], tab$type[zqc])[rownames(dc)],
          ylim=c(0, 0.02), names=dc$name, main="", ylab="Death penetrance (1/h)", col=dc$color, outline=FALSE)
  dev.off()
  
  ## stats on ha
  sp <- split(pheno[zqc, "ha"], tab$type[zqc])
  sapply(sp[c("scrambled", "eg5", "incenp", "bcop")], mean)
  wilcox.test(sp[["eg5"]], sp[["scrambled"]])
  wilcox.test(sp[["bcop"]], sp[["scrambled"]])

  ##
  sp <- split(pheno[zqc, "hmi"], tab$type[zqc])
  sapply(sp[c("scrambled", "eg5", "incenp", "bcop")], mean)
  wilcox.test(sp[["eg5"]], sp[["scrambled"]])
 
  ## plot MCO_0026491
  w <- getspot(sirna="MCO_0026491")
  pdf("expta1.pdf")
  plotfit(w[1], lwd=lwd, kk=c("ka"), cex.axis=cex, cex.lab=cex, kcol=kcol, cex=cex, ylim=c(0, 105))
  dev.off()
  pdf("expta2.pdf")
  plotfit(w[3], lwd=lwd, kk=c("ka"), cex.axis=cex, cex.lab=cex, legend=NULL, kcol=kcol, ylab="", ylim=c(0, 105))
  invisible(dev.off())
}

figure3b <- function() {
  ## to avoid warnings during R CMD check
  if (FALSE) {
    pheno <- NULL
    zqc <- NULL
    tab2 <- NULL
  }
  
  ## graphics parameters
  cex <- 1.5
  lwd <- 2
  kcol <- c("kim"="#555555", "ka"="#3333ff")
  dc <- data.frame(row.names=c('scrambled', 'eg5', 'incenp', 'bcop', 'experiment'),
                   color=c('#77dd77', '#ffaa77', '#ff7777', '#aaaaff', '#aaaaaa'),
                   name=c( 'siScrambled', 'siKIF11', 'siINCENP', 'siCOPB1', 'sample'), stringsAsFactors=FALSE)
  
  ## mitotic duration
  smitod <- split(pheno[zqc, "mitod"], tab$type[zqc])[rownames(dc)]
   
  pdf("boxplotmitod.pdf", width=6, height=5)
  boxplot(smitod, ylim=c(0.3, 15), names=dc$name, ylab="Mitosis duration (h)", col=dc$color, log="y", outline=FALSE)
  dev.off()

  ## stats
  sapply(smitod[c("scrambled", "eg5", "incenp")], median)
  wilcox.test(smitod[["eg5"]], smitod[["scrambled"]])
  wilcox.test(smitod[["incenp"]], smitod[["scrambled"]])
  length(unique(na.omit(tab2$target.hgnc[!is.na(tab2$duration.mitosis)])))

  ## CCDC9
  w <- getspot(sirna="MCO_0020444")
  pdf("expmitod1.pdf")
  plotfit(w[1], lwd=lwd, cex.axis=cex, cex.lab=cex, legend=NULL, kcol=kcol, ylim=c(0, 30))
  dev.off()
  pdf("expmitod2.pdf")
  plotfit(w[3], lwd=lwd, cex.axis=cex, cex.lab=cex, legend=NULL, kcol=kcol, ylab="", ylim=c(0, 30))
  invisible(dev.off())
}

suppFig1 <- function() {
  ## to avoid warnings during R CMD check
  if (FALSE) {
    pheno <- NULL
    zqc <- NULL
    tab2 <- NULL
  }
  
  pdf("suppFig1.pdf")
  plot(tab2$time.mitotic.arrest, tab2$time.death, xlab="Time of mitotic arrest (h)", ylab="Time of cell death (h)")
  text(tab2$time.mitotic.arrest, tab2$time.death, getanno(sirna=rownames(tab2)), col=2, pos=1, cex=0.8)
  abline(a=0, b=1)
  dev.off()

  ## correlation
  cor(tab2$time.mitotic.arrest, tab2$time.death, use="pair")
}

figure4 <- function() {
  ## to avoid warnings during R CMD check
  if (FALSE) {
    pheno <- NULL
    zqc <- NULL
    tab2 <- NULL
  }
   
  ## graphic parameters
  cex <- 1.5
  dc <- data.frame(row.names=c('scrambled', 'eg5', 'incenp', 'bcop', 'experiment'),
                   color=c('#77dd77', '#ffaa77', '#ff7777', '#aaaaff', '#aaaaaa'),
                   name=c( 'siScrambled', 'siKIF11', 'siINCENP', 'siCOPB1', 'sample'), stringsAsFactors=FALSE)
  
  controls <- split(which(zqc), tab$type[zqc])[c("scrambled", "eg5", "bcop")]
  kk1 <- c("him", "hmi", "hmp", "ha")
  kk2 <- c("tim", "tmi", "tmp", "ta")
  dat <- cbind(log2(abs(pheno[,kk1]+0.001)), pheno[,kk2])
  
  ## controls
  xc <- dat[unlist(controls),]
  yc <- rep(names(controls), sapply(controls, length))
  ldao <- lda(xc, yc)
  lxc <- predict(ldao, xc)

  ## replicated experiments
  z <- zqc
  ssirna <- split(which(z), tab$sirna[z])
  mc.cores <- detectCores()
  xr <- mclapply(ssirna, mc.cores=mc.cores, function(w) if (length(w)>1) apply(dat[w,], 2, median) else NULL)
  xr <- do.call(rbind, xr)
  lxr <- predict(ldao, xr)
  
  ## open plot
  pdf("fig4.pdf", width=8, height=8)
  par(mar=c(2, 2, 2, 2) + 0.1)
  plot(lxc$x, col=NA, xlim=c(-5.2, 4), ylim=c(-2.5, 4), xlab="", ylab="")
  
  ## contours
  for (z in names(controls)) {
    zd <- KernSmooth::bkde2D(lxc$x[yc==z,], bandwidth=c(0.5, 0.5), gridsize=c(128, 128))
    zdf <- zd$fhat/sum(zd$fhat)
    qzdf <- c(0.25, 0.5, 0.75)
    lqzdf <- sapply(qzdf, sapply, function(q) optimize(function(level, q) abs(sum(zdf[zdf>level])-q), interval=range(zdf), q=q)$minimum)
    contour(x=zd$x1, y=zd$x2, zdf, col=dc[z, "color"], levels=lqzdf, labels=paste0(100*qzdf, "%"), drawlabels=TRUE,
            add=TRUE, lwd=2, cex=cex, labcex=0.8, vfont=NULL, method="edge")
  }

  ## highlight gene-of-interest
  roi <- c(COPA="MCO_0021034", COPB2="MCO_0027498", TP53AIP1="MCO_0041314", RAN="MCO_0012676", RAB25="MCO_0012570", C3orf26="MCO_0026491",
           PLK1="MCO_0001748", MAP2K4="MCO_0001007", ANAPC1="MCO_0047613", 
           NEK9="MCO_0001725", NEK10="MCO_0016363",
           ULK4="MCO_0001660", PALM="MCO_0032347", CENPK="MCO_0026292",
           CCDC9="MCO_0020444", CABYR="MCO_0029224", CTRL="MCO_0004203",
           GNAS="MCO_0013034", DYDC2="MCO_0006428", GDPD3="MCO_0013839",
           ROR1="MCO_0001314", BEND7="MCO_0031251", BZW1="MCO_0043373", ITGA3="MCO_0009017", ASB12="MCO_0019734",
           DDX39A="MCO_0014521",
           "MCO_0000079", NEK2="MCO_0001705", "MCO_0034140",
           COPB1="MCO_0029665", KIF11="MCO_0016402", scrambled="MCO_0016403"
           )
  pp <- data.frame(mco=roi, hgnc=getanno(sirna=roi), bg="#ffffff", stringsAsFactors=FALSE)
  pp$bg[pp$mco=="MCO_0029665"] <- dc["bcop", "color"]
  pp$bg[pp$mco=="MCO_0016402"] <- dc["eg5", "color"]
  pp$bg[pp$mco=="MCO_0016403"] <- dc["scrambled", "color"]
  pp$hgnc[pp$mco=="MCO_0016403"] <- "scrambled"
  
  points(lxr$x[pp$mco,], bg=pp$bg, col="#000000", pch=21)
  points(lxr$x[rownames(tab2),], bg=pp$bg, col="#00000077", pch=21, cex=0.1)
  text(lxr$x[pp$mco,], pp$hgnc, pos=1, cex=0.7, xpd=NA)
  
  ## close plot
  legend("topleft", legend=c(dc[names(controls), "name"], "sample"), col=c(dc[names(controls), "color"], "#000000"), pch=1)
  invisible(dev.off())
}

lambda.justification <- function() {
  ## to avoid warnings during R CMD check
  if (FALSE) {
    pheno <- NULL
    zqc <- NULL
  }
  
  ## init dat
  controls <- split(which(zqc), tab$type[zqc])[c("scrambled", "eg5", "bcop")]
  kk <- c("him", "hmi", "hmp", "ha", "tim", "tmi", "tmp", "ta")
  dat <- pheno[,kk]

  ## controls
  xc <- dat[unlist(controls),]
  yc <- rep(names(controls), sapply(controls, length))
  ldao <- lda(xc, yc)
  lxc <- predict(ldao, xc)

  mean(lxc$class!=yc) ## lambda=4: 0.058

  ##
  mean(pheno[zqc, "pen"]/pheno[zqc, "score"])
}
