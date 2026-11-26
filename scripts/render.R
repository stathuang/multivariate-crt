# Rebuild figures from the independently aggregated replication-level records.
args<-commandArgs(FALSE);script<-sub('^--file=','',grep('^--file=',args,value=TRUE)[1])
root<-normalizePath(file.path(dirname(script),'..'));setwd(root)
ref<-Sys.getenv('REPLICATION_REFERENCE','reference')
input<-'outputs/rebuilt';out<-'outputs/figures';dir.create(out,recursive=TRUE,showWarnings=FALSE)
stopifnot(file.exists(file.path(input,'verification.json')))
source('experiments/reproducible-study/R/reporting.R')
manifest<-list()
record<-function(asset,data,note) {
  path<-file.path(out,asset);stopifnot(file.exists(path))
  dp<-paste0(path,'.csv');write_report_csv(data,dp)
  manifest[[length(manifest)+1L]]<<-data.frame(asset=asset,rows=nrow(data),basis=note,
    sha256=digest::digest(file=path,algo='sha256',serialize=FALSE))
}
path<-function(asset){p<-file.path(out,asset);dir.create(dirname(p),recursive=TRUE,showWarnings=FALSE);p}
save<-function(p,asset,z,w=8,h=6) {
  ggplot2::ggsave(path(asset),p,width=w,height=h,device='pdf')
  record(asset,z,'Recomputed values; layout may differ from the manuscript PDF.')
}
x<-read.csv(file.path(input,'graphical.csv'));x$complete<-x$failed==0
configs<-jsonlite::fromJSON('experiments/reproducible-study/config/simulations.json',simplifyVector=FALSE)$scenarios
for(sid in unique(x$scenario)) {
  cfg<-Filter(function(c)c$id==sid,configs)[[1L]];z<-x[x$scenario==sid,]
  if(sid=='supergraph_h_l_w') {
    z<-z[z$method=='LM-L1-R-SSR',];z$method<-c('G0','G1','G2')[match(z$graph,c('K=6','K=12','K=18'))]
    p<-report_curve(z,c('G0','G1','G2'),mode='graphs');asset<-'Simulation/Appendix/CRT-h_l_w-LM-R-SSR-OversizedGraphs.pdf'
  } else {
    methods<-legacy_methods(cfg$setting);methods[methods=='Chi-square']<-'Chi-squared test'
    z$method[z$method=='Chi-square']<-'Chi-squared test'
    p<-report_curve(z,methods,allow_incomplete=TRUE)
    asset<-paste0(if(cfg$family=='main')'Simulation/CRT/' else 'Simulation/Appendix/',
      'CRT-',cfg$setting,if(cfg$family=='main')'' else paste0('-',cfg$family),'.pdf')
  }
  p<-p+ggplot2::labs(y='Rejection rate')+ggplot2::theme(legend.text=ggplot2::element_text(size=14))
  if(sid=='main_h_gl_w')p<-p+ggplot2::labs(caption='CGM at theta = 3.5: 398 valid outputs from 400 attempts.')
  b<-ggplot2::ggplot_build(p)$data;zplot<-p$data
  key<-function(g,t)paste(g,sprintf('%.17g',t),sep=':')
  j<-match(key(b[[2]]$group,b[[2]]$x),key(as.integer(zplot$Test),zplot$theta))
  k<-match(key(b[[3]]$group,b[[3]]$x),key(as.integer(zplot$Test),zplot$theta))
  stopifnot(!anyNA(j),!anyNA(k),max(abs(b[[2]]$y-zplot$estimate[j]))<1e-14,
    max(abs(b[[3]]$ymin-(zplot$estimate-zplot$se)[k]))<1e-14,
    max(abs(b[[3]]$ymax-(zplot$estimate+zplot$se)[k]))<1e-14)
  save(p,asset,z)
}
z<-read.csv(file.path(input,'gaussian-power.csv'))
z$setting<-factor(z$setting,levels=c('LD-tight','LD-moderate','p>n target'),
  labels=c('n=60, p=45, t=10','n=80, p=50, t=10','n=80, p=60, t=40'))
z$method[z$method=='LW plug-in CRT (D.10)']<-'LW plug-in CRT'
p<-ggplot2::ggplot(z,ggplot2::aes(theta,estimate,color=method,shape=method,group=method))+
 ggplot2::geom_line()+ggplot2::geom_point()+ggplot2::geom_errorbar(ggplot2::aes(ymin=estimate-se,ymax=estimate+se),width=.02)+
 ggplot2::facet_wrap(~setting,nrow=1)+ggplot2::theme_bw()+ggplot2::labs(x=expression(theta),y='Rejection rate',color='Method',shape='Method')+
 ggplot2::theme(legend.position='bottom')
save(p,'Simulation/CRT/oracle_plugin_power_curve.pdf',z,11,4)
z<-read.csv(file.path(input,'stock-noise.csv'));z$complete<-TRUE
p<-report_curve(z,c('LM-L1-R-SSR','RF-RR','dCRT (L1)','dCRT (RF)'),mode='noise')
save(p,'Application/Stock/Stock-Noise.pdf',z)
z<-read.csv(file.path(input,'nongaussian.csv'));z$distribution<-factor(z$distribution,levels=c('Gaussian','t5','t3'))
p<-ggplot2::ggplot(z,ggplot2::aes(100*alpha,100*estimate))+
 ggplot2::geom_ribbon(ggplot2::aes(ymin=100*lower,ymax=100*upper),fill='#0072B2',alpha=.2)+
 ggplot2::geom_abline(slope=1,intercept=0,linetype=2)+ggplot2::geom_line(color='#0072B2')+
 ggplot2::facet_grid(n~distribution)+ggplot2::theme_bw()+ggplot2::labs(x='Nominal level (%)',y='Rejection rate (%)')
save(p,'Simulation/Appendix/MVN-CRT-nongaussian.pdf',z,8,10)
# Use the manuscript asset name from the frozen inventory, avoiding guessed paths.
inv<-jsonlite::fromJSON('provenance/manuscript-inventory.json',simplifyVector=FALSE)
assets<-vapply(inv$assets,`[[`,'','file')
nonasset<-assets[grepl('nongaussian',assets,ignore.case=TRUE)&grepl('pdf$',assets)]
stopifnot(length(nonasset)==1L)
old<-'Simulation/Appendix/MVN-CRT-nongaussian.pdf'
if(nonasset!=old){file.rename(file.path(out,old),path(nonasset));file.rename(paste0(file.path(out,old),'.csv'),paste0(file.path(out,nonasset),'.csv'));manifest[[length(manifest)]]$asset<-nonasset}
center_asset<-assets[grepl('centering_comparison',assets)&grepl('pdf$',assets)]
stopifnot(length(center_asset)==1L)
z<-read.csv(file.path(input,'centering.csv'));z$case<-factor(z$case,levels=c('null','positive','negative','orthogonal'))
z$method<-factor(z$method,levels=c('original','centered_pool','centered_identity','scalar_twosided'),
  labels=c('Uncentred ridge','Centred: pooled','Centred: unit','Scalar two-sided'))
p<-ggplot2::ggplot(z,ggplot2::aes(n,estimate,color=method,shape=method,group=method))+
 ggplot2::geom_line()+ggplot2::geom_point(position=ggplot2::position_dodge(8))+
 ggplot2::geom_errorbar(ggplot2::aes(ymin=ci_lower,ymax=ci_upper),width=5,position=ggplot2::position_dodge(8))+
 ggplot2::geom_hline(yintercept=.05,linetype=3)+ggplot2::facet_wrap(~case,ncol=2)+ggplot2::theme_bw()+
 ggplot2::labs(x='Sample size n',y='Rejection probability')+ggplot2::theme(legend.position='bottom')
save(p,center_asset,z,9,7)
z<-read.csv(file.path(input,'connected.csv'));z<-z[z$method%in%c('original','centered_pool','centered_identity'),]
z$method<-factor(z$method,levels=c('original','centered_pool','centered_identity'),
  labels=c('Uncentred ridge','Centred: pooled','Centred: unit'))
p<-ggplot2::ggplot(z,ggplot2::aes(r,estimate,color=method,shape=method,group=method))+
 ggplot2::geom_line()+ggplot2::geom_point()+ggplot2::geom_errorbar(ggplot2::aes(ymin=ci_lower,ymax=ci_upper),width=.008)+
 ggplot2::geom_hline(yintercept=.05,linetype=3)+ggplot2::facet_wrap(~b,nrow=1)+ggplot2::theme_bw()+
 ggplot2::labs(x='Conditional signal r',y='Rejection probability')+ggplot2::theme(legend.position='bottom')
ca<-assets[grepl('connected_comparison',assets)&grepl('pdf$',assets)];stopifnot(length(ca)==1L);save(p,ca,z,8,4)
z<-read.csv(file.path(ref,'bcr/UNRESTRICTED_RESULTS.csv'));z<-z[z$method=='LRT',]
p<-ggplot2::ggplot(z,ggplot2::aes(pvalue))+ggplot2::geom_histogram(bins=20,fill='grey',color='white')+
 ggplot2::theme_bw()+ggplot2::labs(x='LRT p-value',y='Frequency')
save(p,'Application/BCR/P-LRT_histogram.pdf',z)
# The graph uses only edges, ticker labels, and sectors; these are public summary data.
source('experiments/readable-experiments/R/application_diagnostics.R')
graph<-readRDS(file.path(ref,'stock/graph.rds'));set.seed(1L)
degree<-rowSums(graph$graph)
graph_summary<-data.frame(vertices=length(degree),maximum_degree=max(degree),median_degree=median(degree),
  edges=sum(graph$graph)/2,edge_density=sum(graph$graph)/(length(degree)*(length(degree)-1)))
stopifnot(graph_summary$vertices==103,graph_summary$maximum_degree==30,graph_summary$median_degree==19,
          abs(graph_summary$edge_density-.184)<.0005)
write.csv(graph_summary,file.path(input,'stock-graph-summary.csv'),row.names=FALSE)
plot_stock_graph(graph,path('Application/Stock/stocks_G.pdf'))
record('Application/Stock/stocks_G.pdf',data.frame(ticker=graph$stocks,sector=graph$sectors),'Fixed graph and labels; seed 1; igraph layout can vary by version.')
gof<-readRDS(file.path(ref,'stock/gof.rds'))
stopifnot(abs((1+sum(gof$randomized>=gof$observed))/(gof$M+1)-gof$pvalues[['F_sum']])<1e-14)
write.csv(do.call(rbind,manifest),file.path(out,'MANIFEST.csv'),row.names=FALSE)
expected<-assets[grepl('pdf$',assets)]
stopifnot(setequal(expected,vapply(manifest,function(x)x$asset,'')))
cat('All',length(manifest),'manuscript figure assets rebuilt; graphical points and SE endpoints verified.\n')
