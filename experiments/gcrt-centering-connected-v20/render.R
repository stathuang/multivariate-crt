# Render all prespecified cells from verified saved results; no sampling.
options(warn=2)
stopifnot(file.exists('run/RESULTS_CHECK.json'))
s <- read.csv('run/summary.csv'); d <- read.csv('run/paired_differences.csv')
stopifnot(nrow(s)==24,nrow(d)==24,all(s$replicates==2500))
methods <- c('original','centered_pool','centered_identity')
cols <- c('#404040','#0072B2','#D55E00')
pdf('run/connected_comparison.pdf',width=7.2,height=3.6,family='Helvetica',useDingbats=FALSE)
par(mfrow=c(1,2),mar=c(3.6,3.8,2.1,.8),oma=c(0,0,2.6,0),mgp=c(2.2,.6,0),tcl=-.25,cex=.85)
for (b in c(-1,1)) {
  plot(NA,xlim=c(-.02,.62),ylim=c(0,.7),axes=FALSE,xlab='Conditional signal r',
    ylab='Rejection probability',main=if(b<0)'Negative background (b = -1)' else 'Positive background (b = 1)')
  axis(1,at=c(0,.15,.30,.60),labels=c('0','0.15','0.30','0.60'));axis(2,las=1);box(bty='l')
  abline(h=.05,lty=3,col='grey65')
  for (j in 1:3) {
    a <- s[s$b==b & s$method==methods[j],]; a <- a[order(a$r),]; xx <- a$r+(j-2)*.004
    segments(xx,a$ci_lower,xx,a$ci_upper,col=cols[j])
    segments(xx-.003,a$ci_lower,xx+.003,a$ci_lower,col=cols[j])
    segments(xx-.003,a$ci_upper,xx+.003,a$ci_upper,col=cols[j])
    lines(xx,a$rate,lty=j,col=cols[j],lwd=1.2);points(xx,a$rate,pch=c(16,17,15)[j],col=cols[j])
  }
}
par(fig=c(0,1,0,1),new=TRUE,mar=c(0,0,0,0),oma=c(0,0,0,0));plot.new()
legend('top',c('Uncentered ridge','Centered: pooled','Centered: unit'),col=cols,lty=1:3,
  pch=c(16,17,15),ncol=3,bty='n',inset=.003,cex=.8,x.intersp=.8)
dev.off()
con <- file('run/connected_rates.tex','w')
writeLines(c('\\begin{tabular}{ccccc}','\\hline',
  '$b$ & $r$ & Uncentered & Centered: pooled & Centered: unit \\\\', '\\hline'),con)
for (b in c(-1,1)) for (r in c(0,.15,.30,.60)) {
  z <- s[s$b==b & s$r==r,]; z <- z[match(methods,z$method),]
  txt <- sprintf('$%.2f\\;[%.2f,%.2f]$',100*z$rate,100*z$ci_lower,100*z$ci_upper)
  writeLines(paste0(sprintf('$%+d$ & %.2f & ',b,r),paste(txt,collapse=' & '),' \\\\'),con)
}
writeLines(c('\\hline','\\end{tabular}'),con);close(con)
con <- file('run/connected_paired.tex','w')
writeLines(c('\\begin{tabular}{ccccc}','\\hline',
  '$b$ & $r$ & Pooled $-$ uncentered & Unit $-$ uncentered & Pooled $-$ unit \\\\', '\\hline'),con)
for (b in c(-1,1)) for (r in c(0,.15,.30,.60)) {
  z <- d[d$b==b & d$r==r,]
  txt <- sprintf('\\shortstack{$%.2f\\;(%.2f)$\\\\$[%.2f,%.2f]$}',
    100*z$difference,100*z$paired_mcse,100*z$ci_lower,100*z$ci_upper)
  writeLines(paste0(sprintf('$%+d$ & %.2f & ',b,r),paste(txt,collapse=' & '),' \\\\'),con)
}
writeLines(c('\\hline','\\end{tabular}'),con);close(con)
cat('Rendered all 24 method cells and 24 paired comparisons.\n')
