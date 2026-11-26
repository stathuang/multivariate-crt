options(warn=2)
s <- read.csv("centering-summary.csv")
methods <- c("original","centered_pool","centered_identity","scalar_twosided")
labels <- c("Uncentered ridge","Vector centering (pooled Gram)",
            "Vector centering (unit weight)","Centered and squared scalar gains")
cases <- c("null","positive","negative","orthogonal")
titles <- c("Null","Positive direction","Negative direction","Weak orthogonal alternative")
colors <- c("#242424","#1765B0","#198451","#AC3D55")
pdf("centering_comparison_M999.pdf",width=10,height=7)
par(mfrow=c(2,2),mar=c(4,4,2,1),oma=c(5,0,3,0))
for(j in seq_along(cases)) {
  z <- s[s$case==cases[j],]
  top <- if(cases[j]=="null") max(z$ci_upper)*110 else 100
  plot(NA,xlim=c(.7,3.3),ylim=c(0,top),xaxt="n",
       xlab="Sample size",ylab="Rejection rate (%)",main=titles[j])
  axis(1,at=1:3,labels=c(40,120,360))
  abline(h=5,lty=3,col="grey55")
  for(m in seq_along(methods)) {
    q <- z[z$method==methods[m],];q <- q[order(q$n),]
    x <- 1:3+(m-2.5)*.07
    arrows(x,q$ci_lower*100,x,q$ci_upper*100,angle=90,code=3,length=.03,
           col=colors[m])
    lines(x,q$rate*100,type="b",pch=c(16,17,15,18)[m],col=colors[m])
  }
}
mtext("Paired comparison with M=999, L=1; 100 datasets per sample size",
      outer=TRUE,side=3,line=1,cex=1.1)
mtext("Intervals: pointwise 95% exact binomial confidence intervals",
      outer=TRUE,side=1,line=3.8,cex=.8)
par(fig=c(0,1,.03,.13),new=TRUE,mar=c(0,0,0,0),oma=c(0,0,0,0))
plot.new()
legend("center",legend=labels,col=colors,pch=c(16,17,15,18),
       lty=1,ncol=2,bty="n",cex=.85)
dev.off()
