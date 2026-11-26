# Shared plot styles for the paper.
legacy_methods <- function(setting) {
  switch(setting,
    l_l_w=c("LM-SST","LM-SSR","F-test","dCRT (L1)"),
    h_l_w=c("LM-L1-R-SST","LM-L1-R-SSR","De-Lasso","dCRT (L1)"),
    l_l_m=c("RF","RF-D","RF-RR","F-test","dCRT (RF)"),
    h_l_m=c("RF","RF-D","RF-RR","De-Lasso","dCRT (RF)"),
    l_gl_w=c("GLM-Dev","RF","Chi-square","dCRT (L1)"),
    h_gl_w=c("GLM-L1-D","GLM-L1-R-SST","RF","CGM","dCRT (L1)"),
    l_gl_m=c("RF","RF-D","RF-RR","Chi-square","dCRT (RF)"),
    h_gl_m=c("RF","RF-D","RF-RR","CGM","dCRT (RF)"), stop("Unknown setting"))
}

report_curve <- function(summary, methods, mode="simulation", provenance="reconstructed", allow_incomplete=FALSE) {
  if(!requireNamespace("ggplot2",quietly=TRUE)) stop("ggplot2 is required for figures")
  z <- summary[summary$method %in% methods,,drop=FALSE]
  if(!nrow(z)) stop("No cells selected")
  incomplete <- any(!z$complete)
  if(incomplete && !allow_incomplete) stop("Incomplete results: use allow_incomplete=TRUE for a labelled preview")
  z$Test <- factor(z$method,levels=methods)
  group <- mode %in% c("power","FDP")
  noise <- mode=="noise"
  labels <- methods
  if(mode=="simulation") labels[labels=="dCRT (L1)"] <- "dCRT"
  n <- length(methods); sizes <- if(noise) rep(c(.6,.6,1.2,.6),length.out=n) else rep(if(group) 1 else .6,n)
  lines <- rep(if(group||noise) c("solid","dashed","dotted","twodash","longdash","dotdash") else c("dotdash","dashed","solid","twodash","longdash","dotted"),length.out=n)
  caption <- paste(c(if(provenance!="reconstructed") paste("Provenance:",provenance),if(incomplete) "INCOMPLETE: estimates use finite recorded outputs only.\nSee counts.csv for failures and missing repetitions."),collapse="\n")
  p <- ggplot2::ggplot(z,ggplot2::aes(x=theta,y=estimate,color=Test,shape=Test,linetype=Test,linewidth=Test))+
    ggplot2::geom_line(na.rm=TRUE)+ggplot2::geom_point(size=if(group)2 else 3,na.rm=TRUE)+
    ggplot2::scale_color_manual(values=scales::hue_pal()(n),breaks=methods,labels=labels,drop=FALSE)+
    ggplot2::scale_shape_manual(values=rep(c(1,2,3,4,5,6,8,15,16,17,18),length.out=n),breaks=methods,labels=labels,drop=FALSE)+
    ggplot2::scale_linetype_manual(values=lines,breaks=methods,labels=labels,drop=FALSE)+
    ggplot2::scale_linewidth_manual(values=sizes,breaks=methods,labels=labels,drop=FALSE)+
    ggplot2::theme_bw()+ggplot2::theme(text=ggplot2::element_text(size=if(noise)18 else 20),axis.text=ggplot2::element_text(size=if(noise)18 else 20),legend.position="right",legend.key.height=grid::unit(if(group)1 else 1.2,"cm"),legend.key.width=grid::unit(1.5,"cm"),plot.caption=ggplot2::element_text(size=10))+
    ggplot2::labs(x=if(noise)expression(sigma) else expression(theta),y=if(mode=="FDP")"FDR" else if(noise)"Rejection rate" else "Power",color=if(group)"Method" else if(mode=="graphs")"Graph" else "Test",shape=if(group)"Method" else if(mode=="graphs")"Graph" else "Test",linetype=if(group)"Method" else if(mode=="graphs")"Graph" else "Test",linewidth=if(group)"Method" else if(mode=="graphs")"Graph" else "Test",caption=caption)+
    ggplot2::coord_cartesian(ylim=c(0,if(group||noise)1 else min(1,max(c(.05,z$estimate+ifelse(is.finite(z$se),z$se,0)),na.rm=TRUE)*1.2)))
  if(group||noise) p <- p+ggplot2::scale_y_continuous(breaks=seq(0,1,.25),labels=sprintf("%.2f",seq(0,1,.25)))
  if(!noise) p <- p+ggplot2::geom_errorbar(ggplot2::aes(ymin=estimate-se,ymax=estimate+se),linetype="solid",width=.1,linewidth=if(group).5 else .4,show.legend=FALSE,na.rm=TRUE)
  if(mode=="FDP") p <- p+ggplot2::geom_hline(yintercept=.1,linetype="dotted",linewidth=.5)
  if(!group) p <- p+ggplot2::geom_hline(ggplot2::aes(yintercept=.05,size="Level = .05"),linetype="dotted")+ggplot2::scale_size_manual(values=c("Level = .05"=.5))+ggplot2::guides(size=ggplot2::guide_legend(title="\n "))
  p
}

write_report_csv <- function(x,path) {
  # 17 significant digits round-trip IEEE doubles; human rounding is table-only.
  y <- x
  for(name in names(y)) if(is.double(y[[name]])) {
    y[[name]] <- vapply(y[[name]],function(v)if(is.na(v))NA_character_ else sprintf("%.17g",v),character(1))
  }
  write.csv(y,path,row.names=FALSE,na="")
}

