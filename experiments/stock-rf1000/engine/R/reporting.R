# Named, completeness-aware reporting. No simulation or data preparation occurs here.
report_cell_keys <- c("scenario", "theta", "graph", "method")
report_key <- function(x, cols=report_cell_keys) do.call(paste, c(x[cols], sep="\034"))

expected_report_cells <- function(configs) {
  if (!length(configs)) stop("At least one configuration is required")
  do.call(rbind, lapply(configs, function(cfg) {
    if(identical(cfg$design,"observed") && cfg$repetitions!=1L) stop("Observed analyses require exactly one repetition")
    methods <- if (identical(cfg$family,"group_selection")) c(cfg$adjustments,cfg$baselines) else unique(c(cfg$statistics, cfg$baselines, cfg$methods))
    graphs <- if (!is.null(cfg$report_graphs)) cfg$report_graphs else if (!is.null(cfg$graph_ids)) cfg$graph_ids else if (!is.null(cfg$graphs)) cfg$graphs else if (!is.null(cfg$graph_bandwidths)) paste0("K=",cfg$graph_bandwidths) else "fitted"
    if (is.list(graphs)) stop("Supply character report_graphs for graph specifications")
    z <- expand.grid(scenario=cfg$id, theta=cfg$theta, graph=as.character(graphs), method=methods, stringsAsFactors=FALSE)
    z$expected <- cfg$repetitions
    z
  }))
}

summarize_results <- function(results, expected, metric="rejection", alpha=.05) {
  required <- c(report_cell_keys, "replicate", "status", if(metric=="rejection") "pvalue" else metric)
  if (!all(required %in% names(results))) stop("Missing result columns: ", paste(setdiff(required,names(results)), collapse=", "))
  if (!all(c(report_cell_keys,"expected") %in% names(expected))) stop("Expected-cell schema is incomplete")
  if (anyNA(results[c(report_cell_keys,"replicate","status")]) || anyNA(expected)) stop("Missing identifiers/expected counts")
  if (any(expected$expected < 1 | expected$expected != as.integer(expected$expected))) stop("Expected counts must be positive integers")
  rk <- report_key(results); ek <- report_key(expected)
  if(anyDuplicated(ek)) stop("Duplicate expected cells")
  if(anyDuplicated(report_key(results,c(report_cell_keys,"replicate")))) stop("Duplicate replicate/method records")
  if(any(!rk %in% ek)) stop("Unexpected result cells")
  if(any(results$replicate < 1 | results$replicate != as.integer(results$replicate))) stop("Invalid replicate identifiers")
  rows <- lapply(seq_len(nrow(expected)), function(i) {
    x <- results[rk==ek[i],,drop=FALSE]; n <- expected$expected[i]
    if(any(x$replicate > n)) stop("Replicate exceeds expected count")
    raw <- if(metric=="rejection") x$pvalue else x[[metric]]
    valid <- x$status=="ok" & is.finite(raw) & raw>=0 & raw<=1
    values <- if(metric=="rejection") as.numeric(raw[valid]<alpha) else raw[valid]
    ans <- expected[i,,drop=FALSE]
    ans$metric <- metric; ans$attempts <- nrow(x); ans$valid <- length(values)
    ans$failed <- sum(!valid); ans$missing <- n-nrow(x)
    warnings <- if("warnings" %in% names(x))ifelse(is.na(x$warnings),"",x$warnings) else rep("",nrow(x))
    ans$warning_records <- sum(nzchar(warnings))
    ans$severe_convergence_warnings <- sum(grepl("convergence not obtained|did not converge|failed to converge|non.?convergen",warnings,ignore.case=TRUE))
    ans$rejections <- if(metric=="rejection") sum(values) else NA_integer_
    ans$estimate <- if(length(values)) mean(values) else NA_real_
    ans$se <- if(length(values)>1) sd(values)/sqrt(length(values)) else NA_real_
    ans$complete <- ans$valid==n && ans$attempts==n
    ans
  })
  do.call(rbind,rows)
}

format_pvalue <- function(p) {
  ifelse(!is.finite(p), "--", ifelse(p < .0005, "$<0.001$", sprintf("%.3f",p)))
}
format_rate <- function(mean,se) ifelse(!is.finite(mean),"--",paste0(sprintf("%.3f",mean)," (",ifelse(is.finite(se),sprintf("%.3f",se),"--"),")"))

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

write_named_table <- function(rows, row_ids, methods, file, kind="rate", row_labels=row_ids, method_labels=methods, provenance="reconstructed") {
  # rows need row_id, method, pvalue or estimate/se, and complete.
  key <- paste(rows$row_id,rows$method,sep="\034")
  if(anyDuplicated(key)) stop("Table contains duplicate named cells")
  cells <- lapply(seq_along(row_ids),function(i) {
    x <- rows[match(paste(row_ids[i],methods,sep="\034"),key),,drop=FALSE]
    val <- if(kind=="pvalue") format_pvalue(x$pvalue) else format_rate(x$estimate,x$se)
    val[is.na(x$complete) | !x$complete] <- paste0(val[is.na(x$complete) | !x$complete],"$^{\\dagger}$")
    paste(c(row_labels[i],val),collapse=" & ")
  })
  stock <- all(row_ids %in% c("T1","T2"))
  breast <- grepl("BCR_",basename(file),fixed=TRUE)
  bs <- intToUtf8(92); endrow <- paste0(" ",bs,bs)
  header <- paste0(paste(c(if(stock)paste0("$",bs,"mc{T}$") else "",method_labels),collapse=" & "),endrow)
  lines <- c(paste0("% Provenance: ",provenance),"% Dagger: incomplete cell; see counts.csv. -- means unavailable, never zero.",paste0(bs,"begin{tabular}{l",if(stock)"|" else "",paste(rep("c",length(methods)),collapse=""),"}"),paste0(bs,"hline"),if(breast)paste0(" & ",bs,"multicolumn{",length(methods),"}{c}{Method}",endrow) else NULL,header,paste0(bs,"hline"),paste0(unlist(cells),endrow),paste0(bs,"hline"),paste0(bs,"end{tabular}"))
  writeLines(lines,file)
}

write_report_csv <- function(x,path) {
  # 17 significant digits round-trip IEEE doubles; human rounding is table-only.
  y <- x
  for(name in names(y)) if(is.double(y[[name]])) {
    y[[name]] <- vapply(y[[name]],function(v)if(is.na(v))NA_character_ else sprintf("%.17g",v),character(1))
  }
  write.csv(y,path,row.names=FALSE,na="")
}

render_study <- function(results, configs, out_dir, allow_incomplete=FALSE, provenance="reconstructed") {
  if(is.character(results)) results <- read.csv(results,check.names=FALSE,stringsAsFactors=FALSE)
  dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)
  expected <- expected_report_cells(configs)
  group_ids <- vapply(configs,function(z) identical(z$family,"group_selection"),logical(1))
  summaries <- list()
  for(i in seq_along(configs)) {
    cfg <- configs[[i]]; r <- results[results$scenario==cfg$id,,drop=FALSE]
    e <- expected[expected$scenario==cfg$id,,drop=FALSE]
    metrics <- if(group_ids[i]) c("FDP","power") else if(identical(cfg$design,"observed")) "pvalue" else "rejection"
    for(metric in metrics) summaries[[length(summaries)+1L]] <- summarize_results(r,e,metric)
  }
  counts <- do.call(rbind,summaries)
  # Preserve full double precision in machine-readable exports.
  write_report_csv(counts,file.path(out_dir,"counts.csv"))
  write_report_csv(results,file.path(out_dir,"results.csv"))
  diagnostic_fields <- intersect(c(report_cell_keys,"replicate","status","objective_status",
    "inference_status","numerical_status_policy","error_code","raw_pvalue",
    "full_objective_gap","reduced_objective_gap"),names(results))
  if ("objective_status" %in% names(results))
    write_report_csv(results[,diagnostic_fields,drop=FALSE],file.path(out_dir,"numerical_status.csv"))
  writeLines(c(paste("Provenance:",provenance),paste("Configured scope:",paste(vapply(configs,`[[`,"","id"),collapse=", ")),"Rejection summaries use p < 0.05 and sample-SD/sqrt(valid N) indicator SE. Group summaries use direct FDP/power with the same SE convention. Observed pvalue rows report a single p-value and no Monte Carlo SE.", "Failed and missing repetitions are excluded from the estimate and reported separately.", "Numerical objective certification is distinct from reference-law calibration; see numerical_status.csv when available.", if(all(counts$complete))"All specified cells have the expected valid counts." else "INCOMPLETE: one or more specified cells lack the expected valid counts."),file.path(out_dir,"RENDER_STATUS.txt"))
  if(any(!counts$complete) && !allow_incomplete) stop("Incomplete counts written; figures withheld. Use allow_incomplete=TRUE for labelled previews.")
  save_plot <- function(p,path) {dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE); ggplot2::ggsave(path,p,width=8,height=6,device="pdf")}
  for(cfg in configs) {
    s <- counts[counts$scenario==cfg$id,,drop=FALSE]
    d <- file.path(out_dir,cfg$id); dir.create(d,showWarnings=FALSE)
    if(cfg$family %in% c("main","univariate","t3")) {
      suffix <- if(cfg$family=="main")"" else paste0("-",cfg$family)
      save_plot(report_curve(s,legacy_methods(cfg$setting),provenance=provenance,allow_incomplete=allow_incomplete),file.path(d,paste0("CRT-",cfg$setting,suffix,".pdf")))
    } else if(cfg$family=="supergraph") {
      for(method in cfg$statistics) {
        z <- s[s$method==method,,drop=FALSE]
        graph_ids <- paste0("K=",cfg$graph_bandwidths)
        z$method <- c("Faithful Graph","Oversized Graph 1","Oversized Graph 2")[match(z$graph,graph_ids)]
        stem <- if(method=="LM-L1-R-SSR") "LM-R-SSR" else method
        save_plot(report_curve(z,c("Faithful Graph","Oversized Graph 1","Oversized Graph 2"),"graphs",provenance,allow_incomplete),file.path(d,paste0("CRT-",cfg$setting,"-",stem,"-OversizedGraphs.pdf")))
      }
    } else if(cfg$family=="group_selection") {
      order <- c("Second-Order-Group","BH","BY","e-BH")
      labels <- c("Group Knockoff","G-CRT (BH)","G-CRT (BY)","G-CRT (eBH)")
      s$method <- labels[match(s$method,order)]
      for(metric in c("FDP","power")) save_plot(report_curve(s[s$metric==metric,,drop=FALSE],labels,metric,provenance,allow_incomplete),file.path(d,paste0(if(metric=="FDP")"fdr" else "power","_groupselection_linear.pdf")))
    } else {
      render_application(results[results$scenario==cfg$id,,drop=FALSE],s,cfg,d,allow_incomplete,provenance)
    }
  }
  render_application_combined(results,counts,configs,out_dir,allow_incomplete,provenance)
  invisible(counts)
}

render_application <- function(results,summary,cfg,out_dir,allow_incomplete,provenance) {
  bcr <- identical(cfg$dataset,"bcr")
  if(!bcr) return(invisible(NULL)) # Stock tables require named T1/T2 rows together.
  methods <- c("LRT","RF","dCRT (L1)","dCRT (RF)")
  labels <- c("LRT","$G$-CRT (RF)","dCRT (L1)","dCRT (RF)")
  observed <- cfg$design=="observed"
  z <- summary; z$row_id <- cfg$id
  if(observed) z$pvalue <- z$estimate
  write_named_table(z,cfg$id,methods,file.path(out_dir,if(observed)"BCR_pvalue.tex" else "BCR_Type-I-error.tex"),if(observed)"pvalue" else "rate",row_labels="",method_labels=labels,provenance=paste(provenance,cfg$panel))
  if(!observed) {
    x <- results[results$method=="LRT" & results$status=="ok" & is.finite(results$pvalue),,drop=FALSE]
    caption <- paste(c(paste("Panel:",cfg$panel),if(provenance!="reconstructed")provenance,if(any(!summary$complete))"INCOMPLETE: see counts.csv."),collapse="\n")
    p <- ggplot2::ggplot(x,ggplot2::aes(x=pvalue))+ggplot2::geom_histogram(bins=20,fill="grey",color="white")+ggplot2::theme_minimal()+ggplot2::theme(text=ggplot2::element_text(size=20),axis.title=ggplot2::element_text(size=18),plot.title=ggplot2::element_text(size=20,hjust=.5),plot.caption=ggplot2::element_text(size=10))+ggplot2::labs(title="Histogram",x="Simulated p-value",y="Frequency",caption=caption)
    ggplot2::ggsave(file.path(out_dir,"P-LRT_histogram.pdf"),p,width=8,height=6)
  }
}

render_application_combined <- function(results,counts,configs,out_dir,allow_incomplete,provenance) {
  stock <- Filter(function(x)identical(x$dataset,"stock"),configs)
  if(!length(stock)) return(invisible(NULL))
  methods <- c("LM-L1-R-SSR","RF-RR","dCRT (L1)","dCRT (RF)")
  for(observed in c(TRUE,FALSE)) {
    cfgs <- Filter(function(x) if(observed)x$id %in% c("stock_observed_T1","stock_observed_T2") else identical(x$design,"stock_fixed_start_null"),stock)
    if(!length(cfgs)) next
    ids <- vapply(cfgs,`[[`,"","id"); z <- counts[counts$scenario %in% ids,,drop=FALSE]
    z$row_id <- vapply(cfgs,`[[`,"","group")[match(z$scenario,ids)]
    if(observed) z$pvalue <- z$estimate
    ms <- if(observed)methods else append(methods,"De-Lasso",after=2)
    if(!all(c("T1","T2") %in% z$row_id)) cat("\nINCOMPLETE combined stock table: the other target group was not supplied; daggered cells are unavailable.\n",file=file.path(out_dir,"RENDER_STATUS.txt"),append=TRUE)
    write_named_table(z,c("T1","T2"),ms,file.path(out_dir,if(observed)"Stock_pvalue.tex" else "Stock_Type-I-error.tex"),if(observed)"pvalue" else "rate",row_labels=c("$\\mc{T}_1$","$\\mc{T}_2$"),provenance=provenance)
  }
  noise <- Filter(function(x)identical(x$design,"stock_added_noise"),stock)
  if(length(noise)) {
    ids <- vapply(noise,`[[`,"","id"); z <- counts[counts$scenario %in% ids,,drop=FALSE]
    # A per-noise-scenario completed input is still an incomplete full noise curve.
    full_grid <- c(.05,.10,.15,.20,.25)
    if(!all(vapply(full_grid,function(a)any(abs(z$theta-a)<1e-10),logical(1)))) {
      if(!allow_incomplete) return(invisible(NULL))
      z$complete <- FALSE
    }
    p <- report_curve(z,methods,"noise",provenance,allow_incomplete)
    ggplot2::ggsave(file.path(out_dir,"Stock-Noise.pdf"),p,width=8,height=6)
  }
}
