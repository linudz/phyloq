#!/usr/bin/env Rscript

parse_args<-function(x){o<-list();for(i in seq(1L,length(x),2L))o[[substring(x[[i]],3L)]]<-x[[i+1L]];o}
write_tsv<-function(x,p)utils::write.table(x,p,sep="\t",quote=FALSE,row.names=FALSE,na="NA")
write_gz<-function(x,p){z<-gzfile(p,"wt");on.exit(close(z));write_tsv(x,z)}
curve<-function(d,g)vapply(g,function(x)mean(d<=x),numeric(1))
area_weights<-function(g){d<-diff(g);c(d[1]/2,(head(d,-1)+tail(d,-1))/2,tail(d,1)/2)}
empirical_two<-function(null,obs)(1+sum(abs(null)>=abs(obs)))/(length(null)+1)

a<-parse_args(commandArgs(trailingOnly=TRUE)); expected_by_type<-c(conditional=as.integer(a$`conditional-replicates`),full=as.integer(a$`full-replicates`));types<-names(expected_by_type); grid<-seq(0,1,length.out=as.integer(a$`grid-points`)); weights<-area_weights(grid)
files<-list.files(a$`search-dir`,pattern="[.]bootstrap[.]rds$",recursive=TRUE,full.names=TRUE);if(!length(files))stop("No bootstrap chunks found.")
objects<-lapply(files,readRDS); if(any(vapply(objects,function(x)!identical(x$schema_version,1L),logical(1))))stop("Invalid schema.")
replicates<-do.call(rbind,lapply(objects,`[[`,"replicates")); key<-interaction(replicates$trait,replicates$simulation,replicates$bootstrap_type,replicates$tail,drop=TRUE);if(anyDuplicated(key))stop("Duplicated simulation rows.")
traits<-sort(unique(replicates$trait)); for(t in traits)for(bt in types){ids<-sort(unique(replicates$simulation[replicates$trait==t&replicates$bootstrap_type==bt]));if(!identical(ids,seq_len(expected_by_type[[bt]])))stop("Incomplete simulations for ",t,"/",bt)}

# Each rows-list element has two depth vectors (top, bottom), aligned to two replicate rows.
depth_rows<-list();z<-1L
for(obj in objects){ rr<-obj$replicates; for(i in seq_len(nrow(rr))){ pair_index<-ceiling(i/2); tail_index<-if(rr$tail[i]=="top")1L else 2L; d<-obj$depths[[pair_index]][[tail_index]]; depth_rows[[z]]<-list(meta=rr[i,],depth=d);z<-z+1L }}
if(length(depth_rows)!=nrow(replicates))stop("Depth alignment failure.")

observed_by_trait<-lapply(traits,function(t){o<-objects[[which(vapply(objects,function(x)x$trait==t,logical(1)))[1]]];list(summary=o$observed_summary,depths=o$observed_depths)});names(observed_by_trait)<-traits
curve_rows<-list();summary_rows<-list();z<-1L;s<-1L
for(t in traits)for(bt in types)for(tail_name in c("top","bottom")){
  selected<-depth_rows[vapply(depth_rows,function(x)x$meta$trait==t&&x$meta$bootstrap_type==bt&&x$meta$tail==tail_name,logical(1))]
  mat<-do.call(rbind,lapply(selected,function(x)curve(x$depth,grid))); obs_depth<-observed_by_trait[[t]]$depths[[tail_name]];obs<-curve(obs_depth,grid);null_med<-apply(mat,2,median);lower<-apply(mat,2,quantile,0.025);upper<-apply(mat,2,quantile,0.975);dev<-obs-null_med
  curve_rows[[z]]<-data.frame(trait=t,bootstrap_type=bt,tail=tail_name,opportunity_depth=grid,observed_cdf=obs,null_median=null_med,null_lower=lower,null_upper=upper,deviation=dev);z<-z+1L
  null_areas<-as.numeric((mat-matrix(null_med,nrow(mat),length(grid),byrow=TRUE))%*%weights);obs_area<-sum(dev*weights)
  summary_rows[[s]]<-data.frame(trait=t,bootstrap_type=bt,tail=tail_name,n_replicates=nrow(mat),n_tail=length(obs_depth),observed_signed_area=obs_area,empirical_p_signed_area=empirical_two(null_areas,obs_area),observed_max_abs_deviation=max(abs(dev)),maximum_deviation_depth=grid[which.max(abs(dev))],model_stability=if(bt=="full")mean(vapply(selected,function(x)x$meta$selected_model==x$meta$generating_model,logical(1)))else 1);s<-s+1L
}
summary<-do.call(rbind,summary_rows);summary$BH_p_signed_area<-ave(summary$empirical_p_signed_area,interaction(summary$bootstrap_type,summary$tail),FUN=function(x)p.adjust(x,"BH"))
write_gz(replicates,a$`replicates-output`);write_tsv(do.call(rbind,curve_rows),a$`curves-output`);write_tsv(summary,a$`summary-output`)
