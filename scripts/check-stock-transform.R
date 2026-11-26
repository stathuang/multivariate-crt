# Deterministic replay from saved prices; no download, new data, or test simulation.
script<-sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1])
root<-normalizePath(file.path(dirname(script),'..'))
args<-commandArgs(TRUE);stopifnot(length(args)==1L)
d<-normalizePath(args[1L],mustWork=TRUE)
out<-file.path(root,'outputs/stock-preprocessing');dir.create(out,recursive=TRUE,showWarnings=FALSE)
pre <- new.env(); load(file.path(d, "Pre-X.RData"), envir = pre)
dat <- new.env(); load(file.path(d, "DataForExp.RData"), envir = dat)
source(file.path(root,'scripts/stock-transform-functions.R'))
returns <- exp(apply(pre$X_array, 2, function(v) diff(log(v)))) - 1
avg <- apply(returns, 2, function(v) tapply(v, (seq_along(v)-1) %/% 5 + 1, mean))
stopifnot(identical(dim(avg), dim(pre$X)), max(abs(pre$X-dat$mat_X)) == 0)
piece <- box <- avg
rows <- vector("list", ncol(avg))
for (j in seq_len(ncol(avg))) {
  z <- avg[,j]; centered <- z-median(z)
  par <- find_best_parameters(centered, c(.1,2), quantile(abs(centered),c(.8,1)))
  piece[,j] <- transform_piece_data(centered,par[1],par[2])
  shift <- -min(z) + .01*diff(range(z))
  positive <- z+shift
  fit <- MASS::boxcox(lm(positive~1), plotit=FALSE)
  lam <- fit$x[which.max(fit$y)]
  # Preserve the executed historical formula, including its lambda convention.
  box[,j] <- (positive^lam-1)/lam
  rows[[j]] <- data.frame(column=j, stock=colnames(avg)[j], median=median(z),
    piece_lambda=par[1], piece_cutoff=par[2],
    piece_shapiro_p=shapiro.test(piece[,j])$p.value,
    box_shift=shift, box_lambda=lam, box_shapiro_p=shapiro.test(box[,j])$p.value)
}
tab <- do.call(rbind,rows)
failed <- which(tab$piece_shapiro_p<.05)
final <- piece; final[,failed] <- box[,failed]
replayed <- scale(final)
tab$adopted_transform <- ifelse(seq_len(ncol(avg)) %in% failed,"Box-Cox","piecewise")
tab$max_abs_saved_difference <- apply(abs(replayed-pre$X),2,max)
write.csv(tab,file.path(out,"stock_transform_parameters.csv"),row.names=FALSE)
cat("R:",R.version.string,"\nMASS:",as.character(packageVersion("MASS")),"\n")
cat("Saved dimensions:",paste(dim(pre$X),collapse=" x "),"\n")
cat("Saved X and DataForExp mat_X difference:",max(abs(pre$X-dat$mat_X)),"\n")
cat("Replayed standardized X maximum absolute difference:",max(abs(replayed-pre$X)),"\n")
print(tab[failed,],digits=16,row.names=FALSE)
stopifnot(length(failed)==2L,all(is.finite(replayed)),max(abs(replayed-pre$X)) < 1e-10)
stopifnot(identical(sort(as.character(tab$stock[failed])),c('MMM','MO')),
  max(abs(tab$box_lambda[match(c('MMM','MO'),tab$stock)]-c(1.4,1.5)))<1e-12)
