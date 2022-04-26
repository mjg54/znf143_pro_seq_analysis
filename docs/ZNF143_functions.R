plot.ma.lattice <- function(ma.df, filename = 'file.name', p = 0.01, title.main = "Differential PRO Expression",log2fold = 0.5)
  {;
  pdf(paste("MA_plot", filename, ".pdf", sep=''), width=3.83, height=3.83);
  print(xyplot(ma.df$log2FoldChange ~ log(ma.df$baseMean, base=10),
               groups=(ma.df$padj < p & abs(ma.df$log2FoldChange) > log2fold & !is.na(ma.df$padj)),
               col=c("black","red"), main=title.main, scales="free", aspect=1, pch=20, cex=0.5,
               ylab=expression("log"[2]~"PRO fold change"), xlab=expression("log"[10]~"Mean of Normalized Counts"),
               par.settings=list(par.xlab.text=list(cex=1.1,font=2), par.ylab.text=list(cex=1.1,font=2))));
  dev.off()
}


plot.ma.lattice.mods <- function(ma.df, filename = 'file.name', 
                                 p = 0.01, title.main = "Differential PRO Expression",log2fold = 0.5)
  {;
  pdf(paste("MA_plot", filename, ".pdf", sep=''), 
      useDingbats = FALSE, width=3.83, height=3.83);
  print(xyplot(ma.df$log2FoldChange ~ log(ma.df$baseMean, base=10),
               groups=ma.df$r1881,
               col=c("grey90",   "#ce228e", "#2290cf" , "grey60"),
                main=title.main, scales="free", aspect=1, pch=20, cex=0.5,
               ylab=expression("log"[2]~"PRO fold change"), 
               xlab=expression("log"[10]~"Mean of Normalized Counts"),
               par.settings=list(par.xlab.text=list(cex=1.1,font=2), 
                                 par.ylab.text=list(cex=1.1,font=2))));
  dev.off()
  }


plot.ma.lattice.unadjusted <- function(ma.df, filename = 'file.name', p = 0.01, title.main = "Differential PRO Expression",log2fold =0.5)
  {;
  pdf(paste("MA_plot", filename, ".pdf", sep=''), width=3.83, height=3.83);
  print(xyplot(ma.df$log2FoldChange ~ log(ma.df$baseMean, base=10),
               groups=(ma.df$pvalue < p & abs(ma.df$log2FoldChange) > log2fold & !is.na(ma.df$pvalue)),
               col=c("black","red"), main=title.main, scales="free", aspect=1, pch=20, cex=0.5,
               ylab=expression("log"[2]~"PRO fold change"), xlab=expression("log"[10]~"Mean of Normalized Counts"),
               par.settings=list(par.xlab.text=list(cex=1.1,font=2), par.ylab.text=list(cex=1.1,font=2))));
  dev.off()
}

get.raw.counts.interval <- function(df, path.to.bigWig, file.prefix = 'H') {
    vec.names = c()
    inten.df=data.frame(matrix(ncol = 0, nrow = nrow(df)))
    
    for (mod.bigWig in Sys.glob(file.path(path.to.bigWig, paste(file.prefix, "*plus_body_0-mer.bigWig", sep ='')))) {
        factor.name = strsplit(strsplit(mod.bigWig, "/")[[1]][length(strsplit(mod.bigWig, "/")[[1]])], '_plus')[[1]][1]
        print(factor.name)
        vec.names = c(vec.names, factor.name)
        loaded.bw.plus = load.bigWig(mod.bigWig)
        print(mod.bigWig)
        print(paste(path.to.bigWig,'/',factor.name, '_minus_body_0-mer.bigWig', sep=''))
        loaded.bw.minus = load.bigWig(paste(path.to.bigWig,'/',factor.name, '_minus_body_0-mer.bigWig', sep=''))
        mod.inten = bed6.region.bpQuery.bigWig(loaded.bw.plus, loaded.bw.minus, df)
        inten.df = cbind(inten.df, mod.inten)
    }
    colnames(inten.df) = vec.names
    r.names = paste(df[,1], ':', df[,2], '-', df[,3],'_', df[,4], sep='')
    row.names(inten.df) = r.names
    return(inten.df)
}

run.deseq.list <- function(mat) {
  sample.conditions = factor(c("untreated","untreated","treated","treated"), levels=c("untreated","treated"))        
  deseq.counts.table = DESeqDataSetFromMatrix(mat, DataFrame(sample.conditions), ~ sample.conditions);
  colData(deseq.counts.table)$condition<-factor(colData(deseq.counts.table)$sample.conditions, levels=c('untreated','treated'));
  dds = DESeq(deseq.counts.table);
  res = results(dds);
  #res = res[order(res$padj),];
  return(res)
}

run.deseq.table <- function(mat) {
  sample.conditions = factor(colnames(mat))        
  deseq.counts.table = DESeqDataSetFromMatrix(mat, DataFrame(sample.conditions), ~ sample.conditions);
  colData(deseq.counts.table)$condition<-factor(colData(deseq.counts.table)$sample.conditions, levels=colnames(mat));
  dds = DESeq(deseq.counts.table);
  return(dds)
}

plotPCAlattice <- function(df, file = 'PCA_lattice.pdf') {  
  perVar = round(100 * attr(df, "percentVar"))
  df = data.frame(cbind(df, sapply(strsplit(as.character(df$name), '_rep'), '[', 1)))
  colnames(df) = c(colnames(df)[1:(ncol(df)-1)], 'unique_condition')
  print(df)
  #get colors and take away the hex transparency
  color.x = substring(rainbow(length(unique(df$unique_condition))), 1,7) 
  
  df$color = NA
  df$alpha.x = NA
  df$alpha.y = NA
  df$colpal = NA
  
  for (i in 1:length(unique(df$unique_condition))) {
    
    df[df$unique_condition == unique(df$unique_condition)[[i]],]$color = color.x[i]   
    #gives replicates for unique condition
    reps_col<- df[df$unique_condition == unique(df$unique_condition)[[i]],]
    #gives number of replicates in unique condition
    replicates.x = nrow(reps_col)
    alx <- rev(seq(0.2, 1, length.out = replicates.x))
    
    #count transparency(alx), convert alx to hex(aly), combain color and transparency(cp)
    for(rep in 1:replicates.x) {
    
      na <- reps_col[rep, ]$name
      df[df$name == na, ]$alpha.x = alx[rep]
      aly = as.hexmode(round(alx * 255))
      df[df$name == na, ]$alpha.y = aly[rep]
      cp = paste0(color.x[i], aly)
      df[df$name == na, ]$colpal = cp[rep]
      #print(df)
    }
  }
  colpal = df$colpal
  df$name = gsub('_', ' ', df$name)
  pdf(file, width=6, height=6, useDingbats=FALSE)
  print(xyplot(PC2 ~ PC1, groups = name, data=df,
               xlab = paste('PC1: ', perVar[1], '% variance', sep = ''),
               ylab = paste('PC2: ', perVar[2], '% variance', sep = ''),
               par.settings = list(superpose.symbol = list(pch = c(20), col=colpal)),
               pch = 20, cex = 1.7,
               auto.key = TRUE,
               col = colpal))
  dev.off()
}


run.deseq.list.2.4 <- function(mat) {
  sample.conditions = factor(c("untreated","untreated","untreated","untreated","treated","treated"), levels=c("untreated","treated"))        
  deseq.counts.table = DESeqDataSetFromMatrix(mat, DataFrame(sample.conditions), ~ sample.conditions);
  colData(deseq.counts.table)$condition<-factor(colData(deseq.counts.table)$sample.conditions, levels=c('untreated','treated'));
  dds = DESeq(deseq.counts.table);
  res = results(dds);
  #res = res[order(res$padj),];
  return(res)
}

run.deseq.list.4.4 <- function(mat) {
  sample.conditions = factor(c("untreated","untreated","untreated","untreated","treated","treated", "treated","treated"), levels=c("untreated","treated"))        
  deseq.counts.table = DESeqDataSetFromMatrix(mat, DataFrame(sample.conditions), ~ sample.conditions);
  colData(deseq.counts.table)$condition<-factor(colData(deseq.counts.table)$sample.conditions, levels=c('untreated','treated'));
  dds = DESeq(deseq.counts.table);
  res = results(dds);
  #res = res[order(res$padj),];
  return(res)
}


#if you ave two replicates this is a reasonable way to look at all pairwise correlations--alternatively I could use lattice and xyplot with a panel function for the correlation
pairwise.scatters.each.time.point.reps.MJG <- function(df1, df2, filename = "~/Desktop/PLOS_Gen_rebuttal/Fig_pairwise_scatter_Guertin_reps.pdf") {
  pdf(filename, width=6.83, height=3.4*(ncol(df1)/2))
  grid.newpage()
  playout = grid.layout((ncol(df1)/2), 2, widths=unit(c(3.4, 3.4),c("inches","inches")),heights=unit(c(3.4,3.4),c("inches", "inches")),respect = matrix(data=1, nrow=(ncol(df1)/2), ncol=2), just='center')
  pushViewport(viewport(layout=playout))
  
  for (i in 1:ncol(df1)) {
    if (i%%2 != 0) {
      pushViewport(viewport(layout.pos.col=1, layout.pos.row=ceiling(i/2)))
    } else {
      pushViewport(viewport(layout.pos.col=2, layout.pos.row=ceiling(i/2)))
    }
    pan = xyplot(log(df1[,i], base=10) ~ log(df2[,i], base=10), #groups=padj < pval,
      col=c("#000000"),
      #main=list(paste(strsplit(colnames(df1)[i], 'cit')[[1]][2],' time point', sep=''), cex=0.85, font =1),
      scales="free",
      aspect=1,
      pch=20,
      cex=0.5,
      ylim = c(0,5.2),
      xlim = c(0,5.2), 
      ylab=expression("log"[10]~"Replicate 1 Intensity"),
      xlab=expression("log"[10]~"Replicate 2 Intensity"),
      par.settings=list(par.xlab.text=list(cex=0.85,font=2),
        par.ylab.text=list(cex=0.85,font=2)),
      )
    print(pan, newpage=FALSE)
    grid.text(LETTERS[i], x = unit(0.03, "npc"), y=unit(0.97,"npc"), gp = gpar(fontsize = 12, fontface = 2))
    grid.text(colnames(df1)[i], y= unit(0.965, 'npc'), x = unit(0.55, 'native'), gp = gpar(fontsize = 12, fontface = 1))
    grid.text(paste("rho =",round(cor(df2[,i],df1[,i], use = "complete.obs", method = c('spearman')), digits=4)), y= unit(0.25, 'npc'), x = unit(0.7, 'npc'), just='left')
    upViewport()
  }
  dev.off()
}

                                       #in lattice
                                        #get data in correct format
lattice.pairwise.scatter <- function(rep1, rep2, file = 'xyplot_correlations.pdf', r.1 = '1', r.2 = '2') {
    data.cat = data.frame(matrix(nrow = 0, ncol =3))
    colnames(data.cat) = c('rep1', 'rep2', 'cond')
    for (i in 1:ncol(rep1.corr)) {
        lattice.form = data.frame(rep1[,i], rep2[,i], colnames(rep1[i]))
        colnames(lattice.form) = c('rep1', 'rep2', 'cond')
        data.cat = rbind(data.cat, lattice.form)
    }
    
    pdf(file, width=10, height=7, useDingbats = FALSE)
    print(xyplot(log(rep1, base=10) ~ log(rep2, base=10)|cond, 
                 data = data.cat,
                 scales=list(x=list(cex=1.0, relation = "free", alternating=c(rep(1, 20))),
                             y =list(cex=1.0, alternating=c(rep(1, 20)))),
                 colors = 'black', 
                                        #layout=c(2,43),
                 aspect=1,
                 cex.axis=1.2,
                                        #ylab=expression("log"[10]~paste0("Replicate ", r.1," Intensity")),
                 ylab = bquote("log"[10]~"Replicate"~.(r.1)~"Intensity"),
                 xlab = bquote("log"[10]~"Replicate"~.(r.2)~"Intensity"),
                 par.settings=list(par.xlab.text=list(cex=1.0,font=1),
                     par.ylab.text=list(cex=1.0,font=1),
                     strip.background=list(col="grey85"),
                     par.main.text=list(cex=1.2, font=1)),
                       ylim = c(0,5.2),
                 xlim = c(0,5.2),
                 between=list(y=1.0, x = 1.0),

                 panel=function(x,y, subscripts, ...)
                     {panel.text(3.8,0.5,paste("rho =",round(cor(x,y, use = "complete.obs", method = c("spearman")),digits=2)))
                      panel.xyplot(x,y, pch=20, cex=0.4, col = "grey25")}))
    dev.off()
    #return(data.cat)
}




pro.composites <- function(bed, dir) {
    bedp = cbind(bed[,c(1,2)], bed[,2] + 1, bed[,c(4,5,6)])[bed[,6] == '+',]
    bedm = cbind(bed[,c(1)], bed[,3] - 1, bed[,c(3,4,5,6)])[bed[,6] == '-',]
    colnames(bedp) <- c('chr', 'start', 'end', 'gene', 'xy', 'strand')
    colnames(bedm) <- c('chr', 'start', 'end', 'gene', 'xy', 'strand')
    bed2 = rbind(bedp, bedm)
    return(bed2)
}

composites.func.panels.pro <- function(dat, fact = 'RNA polymerase', summit = 'Summit', class= '', num=90, 
                                   col.lines = c(rgb(0,0,1,1/2), rgb(1,0,0,1/2),  rgb(0.1,0.5,0.05,1/2), rgb(0,0,0,1/2),  
                                       rgb(1/2,0,1/2,1/2), rgb(0,1/2,1/2,1/2), rgb(1/2,1/2,0,1/2)), fill.poly = c(rgb(0,0,1,1/4), 
              rgb(1,0,0,1/4), rgb(0.1,0.5,0.05,1/4),rgb(0,0,0,1/4), rgb(1/2,0,1/2,1/4))) {
    count = length(unique(dat$grp))
    ct.cons = 0
    lst.cons = list()
    unique(dat$grp)[order(unique(dat$grp))]
    #for (i in unique(dat$grp)[order(unique(dat$grp))]) {
    #    ct.cons= ct.cons + 1
    #  lst.cons[[ct.cons]] = c(min(dat[dat$grp == i,]$lower), max(dat[dat$grp == i,]$upper))
    #}
    pdf(paste('composite_', fact, '_signals_', summit, '_peaks', class, '.pdf', sep=''), width=6.83, 
        height=ceiling((count)) * 3.00) 
    print(xyplot(est ~ x|grp, group = cond, data = dat,
                 type = 'l',
               scales=list(x=list(cex=0.8,relation = "free"), y =list(cex=0.8, relation="free")),
                 xlim=c(-(num),(num)),
                 #ylim = lst.cons,
                 col = col.lines,
                 main=list(label=class, cex=0.6),
                 auto.key = list(points=F, lines=T, cex=0.8),
               par.settings = list(superpose.symbol = list(pch = c(16), col=col.lines, cex =0.5), 
                   superpose.line = list(col = col.lines, lwd=c(2,2,2,2,2,2), 
                       lty = c(1,1,1,1,1,1,1,1,1))),
                 cex.axis=1.0,
                 par.strip.text=list(cex=0.9, font=1, col='black'),
                 aspect=1.0,
               between=list(y=0.5, x=0.5),
                                        #lwd=2,
                 ylab = list(label = paste(fact," Density", sep=''), cex =0.8),
                 xlab = list(label = paste("Distance from ", summit, sep=''), cex =0.8),
                 upper = dat$upper,
                 fill = fill.poly,
                 lower = dat$lower,
                 strip = function(..., which.panel, bg) {
                     bg.col = c("grey85")
                 strip.default(..., which.panel = which.panel, bg = rep(bg.col, length = which.panel)[which.panel])
                 },
                 panel = function(x, y, ...){
                     #panel.superpose(x, y, panel.groups = 'my.panel.bands', ...)
                     panel.xyplot(x, y, ...)
       }
                 ))
    dev.off()
}

window.step.matrix<- function(bed, bigWig.plus, bigWig.minus, halfWindow, step) {
    bed[,2] = bed[,2] - halfWindow
    bed[,3] = bed[,3] + halfWindow
    matrix.comp.two = bed6.step.bpQuery.bigWig(bigWig.plus, bigWig.minus, bed, step, op = "avg", follow.strand = TRUE)
    res = do.call(rbind, matrix.comp.two)
  return(list(res, matrix.comp.two))
}


quantiles.metaprofile.adjust <- function(mat, quantiles = c(0.875, 0.5, 0.125)) {
  stopifnot(length(quantiles) == 3)
  stopifnot(all(quantiles < 1 & quantiles > 0))
  
  qTop = quantiles[1]
  qMid = quantiles[2]
  qBottom = quantiles[3]
  
  N = dim(mat)[2]
  
  cTop = sapply(1:N, function(idx) quantile(mat[, idx], qTop, na.rm = TRUE))
  cMid = sapply(1:N, function(idx) quantile(mat[, idx], qMid, na.rm = TRUE))
  cBottom = sapply(1:N, function(idx) quantile(mat[, idx], qBottom, na.rm = TRUE))
  
  # add step if present
  if (!is.null(attr(mat, "step")))
    return(list(step = attr(mat, "step"), top = cTop, middle = cMid, bottom = cBottom))
  else
    return(list(top = cTop, middle = cMid, bottom = cBottom))
}


subsampled.quantiles.metaprofile.adjust <- function(mat, quantiles = c(0.875, 0.5, 0.125), fraction = 0.10, n.samples = 1000) {
  stopifnot(length(quantiles) == 3)
  stopifnot(all(quantiles < 1 & quantiles > 0))
  
  # create sub-samples
  N = dim(mat)[1]
  M = dim(mat)[2]
  
  
  result = matrix(nrow=n.samples, ncol=M)
  K = as.integer(round(N * fraction, 0))


  for (i in 1:n.samples) {
      idx <- sample(N, size=K, replace=FALSE)
    
    result[i, ] = colMeans(as.matrix(mat[idx, ]), na.rm = TRUE)
  }
  
  
  # propagate step attribute
  attr(result, "step") <- attr(mat, "step")

  return(quantiles.metaprofile.adjust(result, quantiles))
}


create.composites.heatmaps <- function(path.dir, composite.input, region=20, step=1, grp = 'PRO', first=FALSE, file.prefix = 'sh') { 
    vec.names = c('chr','start','end')
    hmap.data = list()
    composite.df=data.frame(matrix(ncol = 6, nrow = 0))
    for (mod.bigWig.plus in Sys.glob(file.path(path.dir, paste(file.prefix, "*plus_combined_normalized.bigWig", sep='')))) {
        mod.bigWig.minus=(paste(strsplit(mod.bigWig.plus, 'plus')[[1]][1], 'minus_combined_normalized.bigWig', sep=''))
        factor.name = strsplit(strsplit(mod.bigWig.plus, "/")[[1]][length(strsplit(mod.bigWig.plus, "/")[[1]])], '_plus')[[1]][1]
        print(factor.name)
        vec.names = c(vec.names, factor.name)
        wiggle.plus = load.bigWig(mod.bigWig.plus)
        wiggle.minus = load.bigWig(mod.bigWig.minus)
        bpqueryPM = window.step.matrix(composite.input, wiggle.plus, wiggle.minus, region, step)
        subsample = subsampled.quantiles.metaprofile.adjust((bpqueryPM[[1]]))
        mult.row = ncol(bpqueryPM[[1]])
        hmap.data[[paste(factor.name,'_plus', sep='')]] = list(colMeans(bpqueryPM[[1]]), subsample$top, subsample$bottom, colMeans(bpqueryPM[[1]]), paste(factor.name,'_plus', sep=''), bpqueryPM[[1]])
        df.up <- data.frame(matrix(ncol = 6, nrow = mult.row))
        df.up[, 1] <- colMeans(bpqueryPM[[1]])
        df.up[, 2] <- seq((-1 * region) + 0.5 * step, region - 0.5 * step, by = step)
        df.up[, 3] <- matrix(data = paste(factor.name,'', sep=''), nrow=mult.row, ncol=1)
        df.up[, 4] <- subsample$top
        df.up[, 5] <- subsample$bottom
        df.up[, 6] <- matrix(data = 'Plus', nrow=mult.row, ncol=1)
        composite.df = rbind(composite.df, df.up)
        bpqueryMP = window.step.matrix(composite.input, wiggle.minus, wiggle.plus, region, step)
        subsample = subsampled.quantiles.metaprofile.adjust((bpqueryMP[[1]]))
        mult.row = ncol(bpqueryMP[[1]])
        hmap.data[[paste(factor.name,'_minus', sep='')]] = list(colMeans(bpqueryMP[[1]]), subsample$top, subsample$bottom, colMeans(bpqueryMP[[1]]), paste(factor.name,'_minus', sep=''), bpqueryMP[[1]])
        df.up <- data.frame(matrix(ncol = 6, nrow = mult.row))
        df.up[, 1] <- colMeans(bpqueryMP[[1]])
        df.up[, 2] <- seq((-1 * region) + 0.5 * step, region - 0.5 * step, by = step)
        df.up[, 3] <- matrix(data = paste(factor.name,'', sep=''), nrow=mult.row, ncol=1)
        df.up[, 4] <- subsample$top
        df.up[, 5] <- subsample$bottom
        df.up[, 6] <- matrix(data = 'Minus', nrow=mult.row, ncol=1)
        composite.df = rbind(composite.df, df.up)
        unload.bigWig(wiggle.plus)
        unload.bigWig(wiggle.minus)

    }
    colnames(composite.df) <- c('est', 'x', 'cond', 'upper', 'lower', 'grp')
    for (cond in (1:length(hmap.data))) {
    rownames(hmap.data[[cond]][[6]]) = paste(composite.input[,1], ':',
                composite.input[,2], '-', composite.input[,3], sep='')
    colnames(hmap.data[[cond]][[6]]) = seq((-1 * region) + 0.5 * step, region - 0.5 * step, by = step)
  }
    return(list(composite.df, hmap.data))
}

convert.deseq.to.bed <- function(deseq.df, gene.file) {
    x = gene.file[gene.file[,4] %in% sapply(strsplit(rownames(deseq.df), '_'), '[', 2),]
    y = pro.composites(x)
    return(y)
}


streamplot <- function(df, directory, reg = 400, stp = 20, genelist, class = 'no treatment', file.prefix = 'sh', order.indx = 1, order.layout = FALSE) {
    x = create.composites.heatmaps(directory, convert.deseq.to.bed(df, genelist), region = reg, step = stp, file.prefix = file.prefix)
    composites.func.panels.pro(x[[1]],'RNA Polymerase', summit = 'TSS', num = reg - (reg *0.2), class = class)
    if (order.layout == FALSE) {order.layout = c(1:length(x[[2]]))}
    draw.heatmap.both(x[[2]], file = paste(class, '_heatmap.pdf', sep=''), upstream = -1*reg, downstream = reg, sub.set = 80000, arrangement = 'left', order.indx = order.indx, order.layout=order.layout)
    return(x)
}


composites.func.panels.pro.log <- function(dat, fact = 'RNA polymerase', summit = 'Summit', class= '', num=90, 
                                   col.lines = c(rgb(0,0,1,1/2), rgb(1,0,0,1/2),  rgb(0.1,0.5,0.05,1/2), rgb(0,0,0,1/2),  
                                       rgb(1/2,0,1/2,1/2), rgb(0,1/2,1/2,1/2), rgb(1/2,1/2,0,1/2)), fill.poly = c(rgb(0,0,1,1/4), 
              rgb(1,0,0,1/4), rgb(0.1,0.5,0.05,1/4),rgb(0,0,0,1/4), rgb(1/2,0,1/2,1/4))) {
    count = length(unique(dat$grp))
    ct.cons = 0
    lst.cons = list()
    unique(dat$grp)[order(unique(dat$grp))]
    #for (i in unique(dat$grp)[order(unique(dat$grp))]) {
    #    ct.cons= ct.cons + 1
    #  lst.cons[[ct.cons]] = c(min(dat[dat$grp == i,]$lower), max(dat[dat$grp == i,]$upper))
    #}
    pdf(paste('composite_', fact, '_signals_', summit, '_peaks', class, '.pdf', sep=''), width=6.83, 
        height=ceiling((count)) * 3.00) 
    print(xyplot(log(est, base = 2) ~ x|grp, group = cond, data = dat,
                 type = 'l',
               scales=list(x=list(cex=0.8,relation = "free"), y =list(cex=0.8, relation="free")),
                 xlim=c(-(num),(num)),
                 #ylim = lst.cons,
                 col = col.lines,
                 main=list(label=class, cex=0.6),
                 auto.key = list(points=F, lines=T, cex=0.8),
               par.settings = list(superpose.symbol = list(pch = c(16), col=col.lines, cex =0.5), 
                   superpose.line = list(col = col.lines, lwd=c(2,2,2,2,2,2), 
                       lty = c(1,1,1,1,1,1,1,1,1))),
                 cex.axis=1.0,
                 par.strip.text=list(cex=0.9, font=1, col='black'),
                 aspect=1.0,
               between=list(y=0.5, x=0.5),
                                        #lwd=2,
                 ylab = list(label = paste(fact," Density", sep=''), cex =0.8),
                 xlab = list(label = paste("Distance from ", summit, sep=''), cex =0.8),
                 upper = dat$upper,
                 fill = fill.poly,
                 lower = dat$lower,
                 strip = function(..., which.panel, bg) {
                     bg.col = c("grey85")
                 strip.default(..., which.panel = which.panel, bg = rep(bg.col, length = which.panel)[which.panel])
                 },
                 panel = function(x, y, ...){
                     #panel.superpose(x, y, panel.groups = 'my.panel.bands', ...)
                     panel.xyplot(x, y, ...)
       }
                 ))
    dev.off()
}
streamplot.all.log <- function(df, directory, reg = 400, stp = 20, genelist, class = 'no treatment') {
    x = create.composites.heatmaps(directory, convert.deseq.to.bed(df, genelist), region = reg, step = stp)
    composites.func.panels.pro.log(x[[1]], 'RNA Polymerase', summit = 'TSS', num = reg - (reg *0.2), class = class)
    return(x)
}

draw.heatmap.both <- function(var.name, ord = var.name, filename = "hmap.data.pdf",
                              width.vp = 3, width.space = 0.8, height.vp = 6,
                              height.legend = 1, upstream = -1200, downstream = 1200,
                              sub.set = 1000, order.layout = c(1:length(var.name)),
                              order.indx = 1, arrangement = c("left", "center") ,
                              colors.ramp = c("white", "#FCDBDB", "#FCC0C0",
                                  "#F98F8F", "#F95A5A", "red", "red3", "red4","#3F0202"),
                              convert.zeros = TRUE) {
    arrangement <- match.arg(arrangement)
    plot.num = length(var.name)
    pdf(paste('', filename, sep=''),w=((width.vp+width.space)*plot.num)+(width.space),
        h=(height.vp+height.legend*2))
    grid.newpage()
    step = abs(as.numeric(as.character(colnames(var.name[[1]][[6]])[1])) -
        as.numeric(as.character(colnames(var.name[[1]][[6]])[2])))
    probes = ((downstream-upstream)/step)/2
    playout = grid.layout(3,(2*plot.num)+1, widths=unit(c(width.space,
                         rep(c(width.vp, width.space), plot.num)),
                         c("inches", rep(c("inches", "inches"), plot.num))),
                         heights=unit(c(height.legend, height.vp,height.legend),
                         c("inches","inches", "inches")), respect =
        matrix(data =1, nrow=3, ncol=(plot.num*2)+1))
    pushViewport(viewport(layout=playout))
    plot.list = list()
    x.ord = apply(var.name[[order.indx]][[6]],1,which.max)
    for (cond in (1:length(var.name))) {
        if (arrangement == 'left') {
            var.name[[cond]][[6]] = var.name[[cond]][[6]][order(x.ord),]
        }
        if (arrangement == 'center') {
            var.name[[cond]][[6]] = var.name[[cond]][[6]][
                                order(-ord[[order.indx]][[6]][,ncol(
                                    var.name[[order.indx]][[6]])/2]), ]
        }
        if (nrow(var.name[[cond]][[6]]) > sub.set) {
            first.k = var.name[[cond]][[6]][1:sub.set,
                ((ncol(var.name[[cond]][[6]])/2)-probes):
                ((ncol(var.name[[cond]][[6]])/2)+probes)]
        } else {
        first.k = var.name[[cond]][[6]][,
            ((ncol(var.name[[cond]][[6]])/2)-probes):
            ((ncol(var.name[[cond]][[6]])/2)+probes)]
    }
        plot.list[[cond]] = first.k
    }
    min.range = Inf
    max.range = -Inf
    for (cond in (1:length(var.name))) {
        min.cond = min(plot.list[[cond]][plot.list[[cond]] != 0], na.rm=TRUE)
        max.cond = max(plot.list[[cond]])
        if (min.cond < min.range) {
            min.range = min.cond
        }
        if (max.cond > max.range) {
            max.range = max.cond
        }
    }
    print(min.range)
    print(max.range)
        
  #plot it!
    count = 0 
    for (cond in (1:length(var.name))) {
        count = count + 1
        print(cond)
        pushViewport(viewport(layout.pos.col=(order.layout[count]*2), layout.pos.row=2))
    
        first.k = plot.list[[cond]][nrow(plot.list[[cond]]):1,]
        if (convert.zeros) {
            first.k[first.k==0] <- min.range
        }
        print(levelplot(t(log(first.k, base = 2)),
                        aspect = height.vp/width.vp,
                        col.regions = colorRampPalette(colors.ramp, bias=1)(150),
                        at = seq(log(min.range, base=2), log(max.range, base = 2), length=150), 
                        xlab="",
                        axes = FALSE,
                        ylab="",
                        main='',
                        sub="",
                        colorkey = FALSE,
                        region = TRUE,
                        scales = list(draw = FALSE)),
              panel.width = list(width.vp, "inches"),
              panel.height=list(height.vp, "inches"),
              newpage = FALSE)
        upViewport()
        
                                        #scale
        pushViewport(viewport(layout.pos.col=(order.layout[count]*2), layout.pos.row=3))   
        grid.lines(x = unit(c(0,1), "npc"),
                   y = unit(c(0.9,0.9), "npc"), gp=gpar(lwd=5))
        grid.lines(x = unit(c(0,0), "npc"),
                   y = unit(c(0.8, 0.9), "npc"), gp=gpar(lwd=5))
        grid.lines(x = unit(c(1,1), "npc"),
                   y = unit(c(0.8, 0.9), "npc"), gp=gpar(lwd=5))
        grid.lines(x = unit(c(0.5,0.5), "npc"),
                   y = unit(c(0.8, 0.9), "npc"), gp=gpar(lwd=5))
        grid.text(as.character(upstream), x = unit(0, "npc"),
                  y = unit(0.6, "npc"), gp = gpar(fontsize = 18, fontface = "bold"))
        grid.text(as.character(downstream), x = unit(1, "npc"),
                  y = unit(0.6, "npc"), gp = gpar(fontsize = 18, fontface = "bold"))
        grid.text("0", x = unit(0.5, "npc"),
                  y = unit(0.6, "npc"), gp = gpar(fontsize = 18, fontface = "bold"))
        
        upViewport()
    
        pushViewport(viewport(layout.pos.col=(order.layout[count]*2), layout.pos.row=1))

        grid.text(var.name[[cond]][[5]][1], x = unit(0.5, "npc"),
                  y = unit(0.5, "npc"), gp = gpar(fontsize = 20, fontface = "bold"))

        upViewport()
    
    }
  
    dev.off()

    pdf(paste("scale_", filename, sep=''), width=6.83, height=4.83)
    grid.newpage()
    playout = grid.layout(1,1, widths=unit(c(6.83),c("inches")),
        heights=unit(c(4.35),c("inches")), respect = matrix(data =1, nrow=1, ncol=1))
    pushViewport(viewport(layout=playout))
    pushViewport(viewport(layout.pos.col=1, layout.pos.row=1))
    pushViewport(viewport(x=unit(0.5, "npc"), width=unit(0.5, "npc"), 
                          y=unit(0.5, "npc"), height=unit(0.25, "npc")))
    hmcols = colorRampPalette(colors.ramp)
    pal <- hmcols
    grid.imageFun(1, length(colors.ramp), pal(length(colors.ramp)))
    
#ths needs to be changed if "Error in validDetails.axis(x) :  'labels' and 'at' locations must have same length"
    if (length(seq(round(log(min.range, base = 2)), round(log(max.range, base =2)),by=round((round(log(max.range, base =2))-round(log(min.range, base =2)))/(length(colors.ramp)),1))) == length(colors.ramp)+1) {
        grid.xaxis(seq(0,1,by=1/(length(colors.ramp))), label=round(seq(round(log(min.range, base = 2)), round(log(max.range, base =2)),by=round((round(log(max.range, base =2))-round(log(min.range, base =2)))/(length(colors.ramp)),1)), 2), gp= gpar(cex = 1))
    } else {
        if (length(seq(round(log(min.range, base = 2)), round(log(max.range, base =2)),by=floor((round(log(max.range, base =2))-round(log(min.range, base =2)))/(length(colors.ramp))/0.1)*0.1)) == length(colors.ramp)+1) {
            grid.xaxis(seq(0,1,by=1/(length(colors.ramp))), label=round(seq(round(log(min.range, base = 2)), round(log(max.range, base =2)),by=floor((round(log(max.range, base =2))-round(log(min.range, base =2)))/(length(colors.ramp))/0.1)*0.1), 2), gp = gpar(cex = 1))
        }
    }
    
    grid.text(expression("log"[2]~"PRO Intensity"), y = unit(-3, "lines"))
    grid.rect()
    dev.off()
    
}




image.scale <- function(z, zlim, col = heat.colors(12), breaks, horiz=TRUE, ...) {
  if(!missing(breaks)) {
    if(length(breaks) != (length(col)+1)){stop("must have one more break than colour")}
  }
  if(missing(breaks) & !missing(zlim)) {
    breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1)) 
  }
  if(missing(breaks) & missing(zlim)) {
    zlim <- range(z, na.rm=TRUE)
    zlim[2] <- zlim[2]+c(zlim[2]-zlim[1])*(1E-3)
    zlim[1] <- zlim[1]-c(zlim[2]-zlim[1])*(1E-3)
    breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1))
  }
  poly <- vector(mode="list", length(col))
  for(i in seq(poly)) {
    poly[[i]] <- c(breaks[i], breaks[i+1], breaks[i+1], breaks[i])
  }
  xaxt <- ifelse(horiz, "s", "n")
  yaxt <- ifelse(horiz, "n", "s")
  if(horiz) {
    ylim<-c(0,1)
    xlim<-range(breaks)
  }
  if(!horiz) {
    ylim<-range(breaks)
    xlim<-c(0,1)
  }
  plot(1,1,t="n",ylim=ylim, xlim=xlim, xaxt=xaxt, yaxt=yaxt, xaxs="i", yaxs="i", ...)  
  for(i in seq(poly)) {
    if(horiz) {
      polygon(poly[[i]], c(0,0,1,1), col=col[i], border=NA)
    }
    if(!horiz) { 
      polygon(c(0,0,1,1), poly[[i]], col=col[i], border=NA)
    }
  }
}

grid.imageFun <- function(nrow, ncol, cols, byrow=TRUE) {
  x <- (1:ncol)/ncol
  y <- (1:nrow)/nrow
  if (byrow) {
    right <- rep(x, nrow)
    top <- rep(y, each=ncol)
  } else {
    right <- rep(x, each=nrow)
    top <- rep(y, ncol)
  }
  grid.rect(x=right, y=top,  
            width=1/ncol, height=1/nrow, 
            just=c("right", "top"),
            gp=gpar(col=NA, fill=cols),
            name="image") 
}


get.avg.counts.interval <- function(df, path.to.bigWig, file.prefix = 'H') {
    vec.names = c()
    inten.df=data.frame(matrix(ncol = 0, nrow = nrow(df)))
    
    for (mod.bigWig in Sys.glob(file.path(path.to.bigWig, paste(file.prefix, "*plus_combined_no_scale.bigWig", sep ='')))) {
        factor.name = strsplit(strsplit(mod.bigWig, "/")[[1]][length(strsplit(mod.bigWig, "/")[[1]])], '_plus')[[1]][1]
        print(factor.name)
        vec.names = c(vec.names, factor.name)
        loaded.bw.plus = load.bigWig(mod.bigWig)
        print(mod.bigWig)
        print(paste(path.to.bigWig,'/',factor.name, '_minus_combined_no_scale.bigWig', sep=''))
        loaded.bw.minus = load.bigWig(paste(path.to.bigWig,'/',factor.name, '_minus_combined_no_scale.bigWig', sep=''))
        mod.inten = bed6.region.bpQuery.bigWig(loaded.bw.plus, loaded.bw.minus, df, op = 'avg')
        inten.df = cbind(inten.df, mod.inten)
    }
    colnames(inten.df) = vec.names
    r.names = paste(df[,1], ':', df[,2], '-', df[,3],'_', df[,4], '_', df[,7], sep='')
    row.names(inten.df) = r.names
    return(inten.df)
}

filter.deseq <- function(deseq.df, filter.df) {
    z = deseq.df[sapply(strsplit(rownames(deseq.df), '_'), '[', 2) %in% sapply(strsplit(rownames(filter.df), '_'), '[', 2),]
    return(z)
}


get.raw.counts.interval.pause.body <- function(df, path.to.bigWig, file.prefix = 'H') {
    vec.names = c()
    inten.df=data.frame(matrix(ncol = 0, nrow = nrow(df)))
    
    for (mod.bigWig in Sys.glob(file.path(path.to.bigWig, paste(file.prefix, "*plus_body_0-mer.bigWig", sep ='')))) {
        factor.name = strsplit(strsplit(mod.bigWig, "/")[[1]][length(strsplit(mod.bigWig, "/")[[1]])], '_plus')[[1]][1]
        print(factor.name)
        vec.names = c(vec.names, factor.name)
        loaded.bw.plus = load.bigWig(mod.bigWig)
        print(mod.bigWig)
        print(paste(path.to.bigWig,'/',factor.name, '_minus_body_0-mer.bigWig', sep=''))
        loaded.bw.minus = load.bigWig(paste(path.to.bigWig,'/',factor.name, '_minus_body_0-mer.bigWig', sep=''))
        mod.inten = bed6.region.bpQuery.bigWig(loaded.bw.plus, loaded.bw.minus, df)
        inten.df = cbind(inten.df, mod.inten)
    }
    colnames(inten.df) = vec.names
    r.names = paste(df[,1], ':', df[,2], '-', df[,3],'_', df[,4], '_', df[,7], sep='')
    row.names(inten.df) = r.names
    return(inten.df)
}


                                        #CDF

cdf <- function(genelist){
  number.in.category = dim(genelist)[1]
  #max distance of any AR binding site from the gene
  max.distance = 2251992

  x = genelist[,2]

  diff.vec = vector(mode = 'integer', length = max.distance)
  frac.vec = vector(mode = 'integer', length = max.distance)

  for (i in 1:max.distance) {
    max.meeting = sum(x <= i)
    #print(max.meeting)
    diff.vec[i] = abs(max.meeting - number.in.category)
  }

  frac.vec = 1-(diff.vec/number.in.category)

  return(frac.vec)
}




get.tss <- function(bedfile) {
    if (ncol(bedfile) > 6) {
        bedfile = bedfile[,c(1:6)]
    }
    for (i in 1:nrow(bedfile)) {
        if (bedfile[i,6] == '+') {
            bedfile[i,3] = bedfile[i,2] + 1
        } else {
            bedfile[i,2] = bedfile[i,3] - 1
        }
    }
    return(bedfile)
}



take.relevant.gene.columns <- function(category.overlap.peaks, distance) {
  #category.overlap.peaks = category.overlap.peaks[abs(category.overlap.peaks$dis) <= distance,]
  combined = category.overlap.peaks[,c(T,T,T, F, F,F,F,T,F,F,T,T,T,F,F,T,T, F,F, T,F,F,T,T,T,F,F,T,T, rep(T,14))]
  colnames(combined) = sapply(strsplit(colnames(combined), "_"), function(x) {do.call("paste", as.list(x[1])) })
  combined = combined[,c(1:4,7,9, 5, 8, 6, 10, 13, 15, 11, 14, 12, 16:29)]
  #heatmap <- data.matrix(combined[, c(rep(c(F),3),rep(c(T),6))])
  #rownames(heatmap) <- paste(combined[,1],":",combined[,2],"-",combined[,3], sep="")
  return(combined)
}


filter.deseq.into.bed <- function(deseq.df, gene.file, cat = 'R1881 Activated') {
    deseq.df = deseq.df[deseq.df$response == cat,] 
    x = paste(gene.file[,1], ':', gene.file[,2], '-', gene.file[,3],'_', gene.file[,4], sep='')
    y = gene.file[x %in% rownames(deseq.df),]
    z = get.tss(y)
    return(z)
}



get.tss <- function(bedfile) {
    if (ncol(bedfile) > 6) {
        bedfile = bedfile[,c(1:6)]
    }
    for (i in 1:nrow(bedfile)) {
        if (bedfile[i,6] == '+') {
            bedfile[i,3] = bedfile[i,2] + 1
        } else {
            bedfile[i,2] = bedfile[i,3] - 1
        }
    }
    return(bedfile)
}



bedTools.closest<-function(functionstring="/usr/local/bin/bedtools/closestBed",bed1,bed2,opt.string="") {
  #create temp files
  #a.file=tempfile()
  #b.file=tempfile()
                                        #out   =tempfile()
    #bed1[,1] = as.character(bed1[,1])
    #bed2[,1] = as.character(bed2[,1])
    #bed1[,2] = as.numeric(as.character(bed1[,2]))
    #bed2[,2] = as.numeric(as.character(bed2[,2]))
    #bed1[,3] = as.numeric(as.character(bed1[,3]))
    #bed2[,3] = as.numeric(as.character(bed2[,3]))
    
    #bed1 =  bed1[do.call(order, bed1[,c(colnames(bed1)[1], colnames(bed1)[2])]),]
    #bed2 =  bed2[do.call(order, bed2[,c(colnames(bed2)[1], colnames(bed2)[2])]),]
    options(scipen =99) # not to use scientific notation when writing out
    
  #write bed formatted dataframes to tempfile
    write.table(bed1,file= 'a.file.bed', quote=F,sep="\t",col.names=F,row.names=F)
    write.table(bed2,file= 'b.file.bed', quote=F,sep="\t",col.names=F,row.names=F)
    
                                        # create the command string and call the command using system()
    command1=paste('sort -k1,1 -k2,2n', 'a.file.bed', '> a.file.sorted.bed')
    cat(command1,"\n")
    try(system(command1))
    command2=paste('sort -k1,1 -k2,2n', 'b.file.bed', '> b.file.sorted.bed')
    cat(command2,"\n")
    try(system(command2))
    
    command=paste(functionstring, opt.string,"-a",'a.file.sorted.bed',"-b",'b.file.sorted.bed',">",'out.file.bed',sep=" ")
    cat(command,"\n")
    try(system(command))
    
    res=read.table('out.file.bed',header=F, comment.char='')
                                        #unlink(a.file);unlink(b.file);unlink(out)
    colnames(res) = c(colnames(bed1), colnames(bed2)[1:ncol(bed2)], 'dis' )
    return(res)
}

bed.100.interval <- function(bed) {
    #take the interval +20 - +120 to look for highest density
    bedp = cbind(bed[,c(1)], bed[,2] + 20,bed[,2] + 120, bed[,c(4,5,6)])[bed[,6] == '+',]
    bedm = cbind(bed[,c(1)], bed[,3] - 120, bed[,3] - 20, bed[,c(4,5,6)])[bed[,6] == '-',]
    colnames(bedp) <- c('chr', 'start', 'end', 'gene', 'xy', 'strand')
    colnames(bedm) <- c('chr', 'start', 'end', 'gene', 'xy', 'strand')
    bed2 = rbind(bedp, bedm)
    return(bed2)
}


composites.func.panels.pro.lattice <- function(dat, fact = 'RNA polymerase', 
                                               summit = 'Summit', class= '', num.m = -200, 
                                               num.p =90, y.low =0, y.high = 0.2, 
                                               col.lines = c(rgb(0,0,1,1/2), rgb(1,0,0,1/2),  
                                               rgb(0.1,0.5,0.05,1/2), rgb(0,0,0,1/2),  
                                               rgb(1/2,0,1/2,1/2), rgb(0,1/2,1/2,1/2), 
                                               rgb(1/2,1/2,0,1/2)), fill.poly = c(rgb(0,0,1,1/4), 
              rgb(1,0,0,1/4), rgb(0.1,0.5,0.05,1/4),rgb(0,0,0,1/4), rgb(1/2,0,1/2,1/4))) {
    count = length(unique(dat$grp))
    ct.cons = 0
    lst.cons = list()
    unique(dat$grp)[order(unique(dat$grp))]
    pdf(paste('composite_', fact, '_signals_', summit, '_peaks', class, '.pdf', sep=''), width=8.83, 
        height=4) 
    print(xyplot(est ~ x|auxin, group = cond, data = dat,
                 type = 'l',
               scales=list(x=list(cex=0.8,relation = "free"), y =list(cex=0.8, relation="free")),
                 xlim=c(num.m,num.p),
                 ylim = c(y.low, y.high),
                 col = col.lines,
                 #main=list(label=class, cex=0.6),
                 auto.key = list(points=F, lines=T, cex=0.8),
               par.settings = list(superpose.symbol = list(pch = c(16), col=col.lines, cex =0.5), 
                   superpose.line = list(col = col.lines, lwd=c(2,2,2,2,2,2), 
                       lty = c(1,1,1,1,1,1,1,1,1))),
                 cex.axis=1.0,
                 par.strip.text=list(cex=0.9, font=1, col='black'),
                 aspect=1.0,
                 between=list(y=0.5, x=0.5),
                 index.cond = list(c(1,3,4,2)),
                 lwd=2,
                 ylab = list(label = paste(fact," Density", sep=''), cex =0.8),
                 xlab = list(label = paste("Distance from ", summit, sep=''), cex =0.8),
                 upper = dat$upper,
                 fill = fill.poly,
                 lower = dat$lower,
                 strip = function(..., which.panel, bg) {
                     bg.col = c("#521e77", "grey90","#36771e", "grey60")
                 strip.default(..., which.panel = which.panel, bg = rep(bg.col, length = which.panel)[which.panel])
                 },
               panel = function(x, y, ...){
                   panel.superpose(x, y, panel.groups = 'my.panel.bands', ...)
                   panel.xyplot(x, y, ...)
       }
                 ))
    dev.off()
}

composites.func.Activated <- function(dat, fact = 'RNA polymerase', summit = 'Summit', class= '', num.m = -200, 
                                      num.p =90, y.low =0, y.high = 0.2,
                                   col.lines = c(rgb(0,0,1,1/2), rgb(1,0,0,1/2),  
                                                 rgb(0.1,0.5,0.05,1/2), rgb(0,0,0,1/2),  
                                       rgb(1/2,0,1/2,1/2), rgb(0,1/2,1/2,1/2), rgb(1/2,1/2,0,1/2)), 
                                   fill.poly = c(rgb(0,0,1,1/4), 
              rgb(1,0,0,1/4), rgb(0.1,0.5,0.05,1/4),rgb(0,0,0,1/4), rgb(1/2,0,1/2,1/4))) {
    count = length(unique(dat$grp))
    ct.cons = 0
    lst.cons = list()
    unique(dat$grp)[order(unique(dat$grp))]
    pdf(paste('composite_', fact, '_signals_', summit, '_peaks', class, '.pdf', sep=''), width=3.83, 
        height=3.83) 
    print(xyplot(est ~ x, group = cond, data = dat,
                 type = 'l',
               scales=list(x=list(cex=0.8,relation = "free"), y =list(cex=0.8, relation="free")),
                 xlim=c(num.m,num.p),
                 ylim = c(y.low, y.high),
                 col = col.lines,
                 #main=list(label=class, cex=0.6),
                 auto.key = list(points=F, lines=T, cex=0.8),
               par.settings = list(superpose.symbol = list(pch = c(16), col=col.lines, cex =0.5), 
                   superpose.line = list(col = col.lines, lwd=c(2,2,2,2,2,2), 
                       lty = c(1,1,1,1,1,1,1,1,1))),
                 cex.axis=1.0,
                 par.strip.text=list(cex=0.9, font=1, col='black'),
                 aspect=1.0,
                 between=list(y=0.5, x=0.5),
                 #index.cond = list(c(2, 4, 3, 1)),
                 lwd=2,
                 ylab = list(label = paste(fact," Density", sep=''), cex =0.8),
                 xlab = list(label = paste("Distance from ", summit, sep=''), cex =0.8),
                 upper = dat$upper,
                 fill = fill.poly,
                 lower = dat$lower,
                 strip = function(..., which.panel, bg) {
                     bg.col = c("grey90","#ce228e" , "#2290cf","grey60")
                 strip.default(..., which.panel = which.panel, bg = rep(bg.col, length = which.panel)[which.panel])
                 },
                 panel = function(x, y, ...){
                     panel.xyplot(x, y, ...)
       }
                 ))
    dev.off()
}

composites.func.Repressed <- function(dat, fact = 'RNA polymerase', summit = 'Summit', class= '', num.m = -200, 
                                      num.p =90, y.low =0, y.high = 0.2,
                                   col.lines = c(rgb(0,0,1,1/2), rgb(1,0,0,1/2),  
                                                 rgb(0.1,0.5,0.05,1/2), rgb(0,0,0,1/2),  
                                       rgb(1/2,0,1/2,1/2), rgb(0,1/2,1/2,1/2), rgb(1/2,1/2,0,1/2)), 
                                   fill.poly = c(rgb(0,0,1,1/4), 
              rgb(1,0,0,1/4), rgb(0.1,0.5,0.05,1/4),rgb(0,0,0,1/4), rgb(1/2,0,1/2,1/4))) {
    count = length(unique(dat$grp))
    ct.cons = 0
    lst.cons = list()
    unique(dat$grp)[order(unique(dat$grp))]
    pdf(paste('composite_', fact, '_signals_', summit, '_peaks', class, '.pdf', sep=''), width=3.83, 
        height=3.83) 
    print(xyplot(est ~ x, group = cond, data = dat,
                 type = 'l',
               scales=list(x=list(cex=0.8,relation = "free"), y =list(cex=0.8, relation="free")),
                 xlim=c(num.m,num.p),
                 ylim = c(y.low, y.high),
                 col = col.lines,
                 #main=list(label=class, cex=0.6),
                 auto.key = list(points=F, lines=T, cex=0.8),
               par.settings = list(superpose.symbol = list(pch = c(16), col=col.lines, cex =0.5), 
                   superpose.line = list(col = col.lines, lwd=c(2,2,2,2,2,2), 
                       lty = c(1,1,1,1,1,1,1,1,1))),
                 cex.axis=1.0,
                 par.strip.text=list(cex=0.9, font=1, col='black'),
                 aspect=1.0,
                 between=list(y=0.5, x=0.5),
                 #index.cond = list(c(2, 4, 3, 1)),
                 lwd=2,
                 ylab = list(label = paste(fact," Density", sep=''), cex =0.8),
                 xlab = list(label = paste("Distance from ", summit, sep=''), cex =0.8),
                 upper = dat$upper,
                 fill = fill.poly,
                 lower = dat$lower,
                 strip = function(..., which.panel, bg) {
                     bg.col = c("grey90","#ce228e" , "#2290cf","grey60")
                 strip.default(..., which.panel = which.panel, bg = rep(bg.col, length = which.panel)[which.panel])
                 },
                 panel = function(x, y, ...){
                     panel.xyplot(x, y, ...)
       }
                 ))
    dev.off()
}
#define the pause region as the 50-bp interval with the highest read density within −500 to +500 bp of the TSS. 

raw.pause.interval.slide <- function(df, path.to.bigWig, combined.plus.bw, combined.minus.bw, file.prefix = 'shV') {
    vec.names = c()
    inten.df=data.frame(matrix(ncol = 0, nrow = nrow(df)))
    df.x = df
    df.x[7] = NA
    df.x[7][df.x[,6] == '+',] = df.x[3][df.x[,6] == '+',] 
    df.x[7][df.x[,6] == '-',] = df.x[2][df.x[,6] == '-',] 
    df.x[3][df.x[,6] == '+',] = df.x[2][df.x[,6] == '+',] + 500
    df.x[2][df.x[,6] == '-',] = df.x[3][df.x[,6] == '-',] - 500
    df.x[2][df.x[,6] == '+',] = df.x[3][df.x[,6] == '+',] - 1000
    df.x[3][df.x[,6] == '-',] = df.x[2][df.x[,6] == '-',] + 1000
    dat.x = bed6.step.probeQuery.bigWig(combined.plus.bw, combined.minus.bw, df.x, step = 1, as.matrix = TRUE, op = "avg", follow.strand = FALSE)
    dat.x[is.na(dat.x)] <- 0
                                        #need to get new matrix with rollapply
    dat.x.win = t(apply(dat.x, 1, function(x){rollapply(x, width = 50, FUN = mean, by = 1, by.column = FALSE,align = "left")}))
    
    index.max = apply(dat.x.win, 1, which.max)
    the.max = apply(dat.x.win, 1, max)
    
    index.max[the.max < .0201 & df.x[,6] == '+'] = 525
    index.max[the.max < .0201 & df.x[,6] == '-'] = 425
    df.x.pause = df.x
    df.x.pause[2][df.x.pause[,6] == '-',] = df.x.pause[2][df.x[,6] == '-',] + index.max[df.x.pause[,6] == '-']
    df.x.pause[2][df.x.pause[,6] == '+',] = df.x.pause[2][df.x[,6] == '+',] + index.max[df.x.pause[,6] == '+']
    df.x.pause[3] = df.x.pause[2] + 50
    for (mod.bigWig in Sys.glob(file.path(path.to.bigWig, 
                                          paste(file.prefix, "*plus_body_0-mer.bigWig", sep ='')))) {
        factor.name = strsplit(strsplit(mod.bigWig, "/")[[1]][length(strsplit(mod.bigWig, "/")[[1]])], 
                               '_plus')[[1]][1]
        print(factor.name)
        vec.names = c(vec.names, factor.name)
        loaded.bw.plus = load.bigWig(mod.bigWig)
        print(mod.bigWig)
        print(paste(path.to.bigWig,'/',factor.name, '_minus_body_0-mer.bigWig', sep=''))
        loaded.bw.minus = load.bigWig(paste(path.to.bigWig,'/',factor.name, 
                                            '_minus_body_0-mer.bigWig', sep=''))
        mod.inten = bed6.region.bpQuery.bigWig(loaded.bw.plus, loaded.bw.minus, df.x.pause)
        inten.df = cbind(inten.df, mod.inten)
    }
    colnames(inten.df) = vec.names
    r.names = paste(df[,1], ':', df[,2], '-', df[,3],'_', df[,4], sep='')
    row.names(inten.df) = r.names
    return(inten.df)
}




raw.body.interval <- function(df, path.to.bigWig, file.prefix = 'shV') {
    vec.names = c()
    inten.df=data.frame(matrix(ncol = 0, nrow = nrow(df)))
    df.x = df
    print(head(df.x))
    df.x[3][df.x[,6] == '-',] = df.x[3][df.x[,6] == '-',] - 200
    df.x[2][df.x[,6] == '+',] = df.x[2][df.x[,6] == '+',] + 200
    print(head(df.x))
    for (mod.bigWig in Sys.glob(file.path(path.to.bigWig, 
                                          paste(file.prefix, "*plus_body_0-mer.bigWig", sep ='')))) {
        factor.name = strsplit(strsplit(mod.bigWig, "/")[[1]][length(strsplit(mod.bigWig, "/")[[1]])], 
                               '_plus')[[1]][1]
        print(factor.name)
        vec.names = c(vec.names, factor.name)
        loaded.bw.plus = load.bigWig(mod.bigWig)
        print(mod.bigWig)
        print(paste(path.to.bigWig,'/',factor.name, '_minus_body_0-mer.bigWig', sep=''))
        loaded.bw.minus = load.bigWig(paste(path.to.bigWig,'/',
                                            factor.name, '_minus_body_0-mer.bigWig', sep=''))
        mod.inten = bed6.region.bpQuery.bigWig(loaded.bw.plus, loaded.bw.minus, df.x)
        inten.df = cbind(inten.df, mod.inten)
    }
    colnames(inten.df) = vec.names
    r.names = paste(df[,1], ':', df[,2], '-', df[,3],'_', df[,4], sep='')
    row.names(inten.df) = r.names
    return(inten.df)
}


bwplot.input.pause <- function(path.dir, bed.input, combined.plus.bw, combined.minus.bw, file.prefix = 'sh') { 
    vec.names = c('chr', 'start', 'end', 'gene', 'class', 'strand')
    inten.df=data.frame(matrix(ncol = 0, nrow = nrow(bed.input)))
    df.x = bed.input
    df.x[3][df.x[,6] == '+',] = df.x[2][df.x[,6] == '+',] + 500
    df.x[2][df.x[,6] == '-',] = df.x[3][df.x[,6] == '-',] - 500
    df.x[2][df.x[,6] == '+',] = df.x[3][df.x[,6] == '+',] - 1000
    df.x[3][df.x[,6] == '-',] = df.x[2][df.x[,6] == '-',] + 1000
    dat.x = bed6.step.probeQuery.bigWig(combined.plus.bw, combined.minus.bw, df.x, step = 1, as.matrix = TRUE, op = "avg", follow.strand = FALSE)
    dat.x[is.na(dat.x)] <- 0
                                        #need to get new matrix with rollapply
    dat.x.win = t(apply(dat.x, 1, function(x){rollapply(x, width = 50, FUN = mean, by = 1, by.column = FALSE,align = "left")}))
    
    index.max = apply(dat.x.win, 1, which.max)
    the.max = apply(dat.x.win, 1, max)
    
    #index.max[the.max < .0201 & df.x[,6] == '+'] = 525
    #index.max[the.max < .0201 & df.x[,6] == '-'] = 425
    df.x.pause = df.x
    df.x.pause[2][df.x.pause[,6] == '-',] = df.x.pause[2][df.x[,6] == '-',] + index.max[df.x.pause[,6] == '-']
    df.x.pause[2][df.x.pause[,6] == '+',] = df.x.pause[2][df.x[,6] == '+',] + index.max[df.x.pause[,6] == '+']
    df.x.pause[3] = df.x.pause[2] + 50
    bed.input = df.x.pause
    y = bed.input
    for (mod.bigWig.plus in Sys.glob(file.path(path.dir, 
                                               paste(file.prefix, "*plus_combined_normalized.bigWig", sep='')))) {
        mod.bigWig.minus=(paste(strsplit(mod.bigWig.plus, 'plus')[[1]][1], 
                                'minus_combined_normalized.bigWig', sep=''))
        factor.name = strsplit(strsplit(mod.bigWig.plus, "/")[[1]][length(strsplit(mod.bigWig.plus, "/")[[1]])], 
                               '_plus')[[1]][1]
        print(factor.name)
        vec.names = c(vec.names, factor.name)
        wiggle.plus = load.bigWig(mod.bigWig.plus)
        wiggle.minus = load.bigWig(mod.bigWig.minus)
        x = bed6.region.bpQuery.bigWig(wiggle.plus, wiggle.minus, bed.input, op = "avg")
        y = cbind(y, x)
                #print(y[1:10,])
    }
    colnames(y) = vec.names
    return(y)
}



bwplot.input.body <- function(path.dir, bed.input, gene.annotations, file.prefix = 'sh') {
    reord = gene.annotations[order(match(gene.annotations[,4],bed.input[,4])),]
    vec.names = c('chr', 'start', 'end', 'gene', 'class', 'strand')
    print(head(reord))
    reord[2][reord[,6] == '+',] = reord[2][reord[,6] == '+',] + 200
    reord[3][reord[,6] == '-',] = reord[3][reord[,6] == '-',] - 200
    print(head(reord))
    y = cbind(reord[,1:3], bed.input[,4:6])
    for (mod.bigWig.plus in Sys.glob(file.path(path.dir, 
                                               paste(file.prefix, "*plus_combined_normalized.bigWig", sep='')))) {
        mod.bigWig.minus=(paste(strsplit(mod.bigWig.plus, 'plus')[[1]][1], 'minus_combined_normalized.bigWig', sep=''))
        factor.name = strsplit(strsplit(mod.bigWig.plus, "/")[[1]][length(strsplit(mod.bigWig.plus, "/")[[1]])], 
                               '_plus')[[1]][1]
        print(factor.name)
        vec.names = c(vec.names, factor.name)
        wiggle.plus = load.bigWig(mod.bigWig.plus)
        wiggle.minus = load.bigWig(mod.bigWig.minus)
        x = bed6.region.bpQuery.bigWig(wiggle.plus, wiggle.minus, reord, op = "avg")
        y = cbind(y, x)
                #print(y[1:10,])
    }
    colnames(y) = vec.names
    return(y)
}


run.deseq.list.any <- function(mat, unt = 2, trt =2) {
  sample.conditions = factor(c(rep("untreated",unt), rep("treated",trt)), levels=c("untreated","treated"))        
  deseq.counts.table = DESeqDataSetFromMatrix(mat, DataFrame(sample.conditions), ~ sample.conditions);
  colData(deseq.counts.table)$condition<-factor(colData(deseq.counts.table)$sample.conditions, levels=c('untreated','treated'));
  dds = DESeq(deseq.counts.table);
  res = results(dds);
  #res = res[order(res$padj),];
  return(res)
}

categorize.deseq.df <- function(df, fdr = 0.1, log2fold = 0.0, treat
= 'Auxin') {

     df.activated = data.frame(matrix(nrow = 0, ncol = 0))
     df.repressed = data.frame(matrix(nrow = 0, ncol = 0))

     if (nrow(df[df$padj < fdr & !is.na(df$padj) & df$log2FoldChange > log2fold,]) != 0) {
     	df.activated = df[df$padj < fdr & !is.na(df$padj) & df$log2FoldChange > log2fold,]
	df.activated$response = paste(treat, 'Activated')
	}

     if (nrow(df[df$padj < fdr & !is.na(df$padj) & df$log2FoldChange < -log2fold,]) != 0) {
     	df.repressed = df[df$padj < fdr & !is.na(df$padj) & df$log2FoldChange < -log2fold,]
	df.repressed$response = paste(treat, 'Repressed')
	}
    
    df.unchanged = df[df$padj > 0.5 & !is.na(df$padj) & abs(df$log2FoldChange) < 0.25,]
    df.unchanged$response = paste(treat, 'Unchanged')

    df.dregs = df[!(df$padj < fdr & !is.na(df$padj) & df$log2FoldChange > log2fold) &
                  !(df$padj < fdr & !is.na(df$padj) & df$log2FoldChange < -log2fold) &
                  !(df$padj > 0.5 & !is.na(df$padj) &
    		  abs(df$log2FoldChange) < 0.25), ]
    df.dregs$response = paste(treat, 'All Other Genes')
    
    df.effects.lattice = 
    rbind(df.activated, 
          df.unchanged, 
          df.repressed, 
          df.dregs)
    df.effects.lattice$response = factor(df.effects.lattice$response)
    df.effects.lattice$response = relevel(df.effects.lattice$response, ref = paste(treat, 'Unchanged'))
    df.effects.lattice$response = relevel(df.effects.lattice$response, ref = paste(treat, 'All Other Genes'))
    return(df.effects.lattice)
}


filter.deseq.into.bed <- function(deseq.df, gene.file, cat = 'R1881 Activated') {
    deseq.df = deseq.df[deseq.df$arfauxin == cat,]
                                        #print(dim(deseq.df))
    #scientific notation was messign this up occasionally
    x = paste(gene.file[,1], ':', gene.file[,2], '-', gene.file[,3],'_', gene.file[,4], sep='')
    #print(length(x))
    y = gene.file[x %in% rownames(deseq.df),]
    #print(dim(y))
    z = get.tss(y)
    #print(dim(z))
    return(z)
}


cdf.deseq.df <- function(df, genes = gene.file,
                                chip.peaks = 'Znf143_K562_GM12878_peaks_merged.sorted.bed',
                                treat = 'Auxin', tf.name = 'ZNF143') {
    bed.tss.activated = filter.deseq.into.bed(df, genes, cat = paste(treat, 'Activated'))
    paste(head(bed.tss.activated))
    bed.tss.repressed = filter.deseq.into.bed(df, genes, cat = paste(treat, 'Repressed'))
    bed.tss.unchanged = filter.deseq.into.bed(df, genes, cat = paste(treat, 'Unchanged'))
    bed.tss.dregs = filter.deseq.into.bed(df, genes, cat = paste(treat, 'All Other Genes'))
    peaks = read.table(chip.peaks, sep = '\t')
    print(head(peaks))
    act.distance = bedTools.closest(bed1 = bed.tss.activated, bed2 = peaks[,c(1:3)], opt.string = '-D a')
    unreg.distance = bedTools.closest(bed1 = bed.tss.unchanged, bed2 = peaks[,c(1:3)], opt.string = '-D a')
    repress.distance = bedTools.closest(bed1 = bed.tss.repressed, bed2 = peaks[,c(1:3)], opt.string = '-D a')
    dregs.distance = bedTools.closest(bed1 = bed.tss.dregs, bed2 = peaks[,c(1:3)], opt.string = '-D a')

    df.up.can = cbind(act.distance[,c(4, 10)], paste(treat, 'Activated'), "Canonical")
    df.un.can = cbind(unreg.distance[,c(4, 10)], paste(treat, 'Unchanged'), "Canonical")
    df.down.can = cbind(repress.distance[,c(4, 10)], paste(treat, 'Repressed'), "Canonical")
    df.dregs.can = cbind(dregs.distance[,c(4, 10)], cat = paste(treat, 'All Other Genes'), "Canonical")
    
    colnames(df.up.can) = c(colnames(df.up.can)[1:2], 'status', 'er')
    colnames(df.un.can) = c(colnames(df.up.can)[1:2], 'status', 'er')
    colnames(df.down.can) = c(colnames(df.up.can)[1:2], 'status', 'er')
    colnames(df.dregs.can) = c(colnames(df.up.can)[1:2], 'status', 'er')

    df.all = rbind(df.up.can, df.un.can, df.down.can, df.dregs.can)
    return(df.all)
}
