#!/usr/bin/env Rscript

parse_args <- function(x) { out<-list(); for(i in seq(1L,length(x),2L)) out[[substring(x[[i]],3L)]]<-x[[i+1L]]; out }
stable_hash <- function(x) { h<-0; for(b in utf8ToInt(x)) h<-(h*131+b)%%2147483000; as.integer(h) }
fit_models <- function(tree,v) { f<-list(BM=suppressWarnings(geiger::fitContinuous(tree,v,model="BM",ncores=1)),OU=suppressWarnings(geiger::fitContinuous(tree,v,model="OU",ncores=1))); validate_fit(f$BM,"BM"); validate_fit(f$OU,"OU"); f }
covs <- function(tree,f) list(BM=bm_covariance(tree,f$BM$opt$sigsq),OU=ou_covariance(tree,f$OU$opt$alpha,f$OU$opt$sigsq))
gls_mean <- function(v,cov) { one<-rep(1,length(v)); w<-tryCatch(solve(cov,one),error=function(e) qr.solve(cov,one)); z<-sum(w*v)/sum(w); if(is.finite(z)) z else mean(v) }
simulate_trait <- function(mu,cov) { e<-eigen((cov+t(cov))/2,symmetric=TRUE); tol<-1e-9*max(1,max(abs(e$values))); if (min(e$values) < (-tol)) stop("Covariance is not PSD."); setNames(as.numeric(mu+e$vectors%*%(sqrt(pmax(e$values,0))*rnorm(nrow(cov)))),rownames(cov)) }
tail_profile <- function(scores,model,p) { n<-nrow(scores); k<-max(1L,ceiling(n*p)); r<-order(scores$FinalScore,decreasing=TRUE,method="radix"); depth<-(rank(scores$PatristicDistance,ties.method="average")-0.5)/n; idx<-list(top=r[seq_len(k)],bottom=tail(r,k)); list(
  stats=do.call(rbind,lapply(names(idx),function(tn){ii<-idx[[tn]];data.frame(tail=tn,n_pairs=n,n_tail=k,median_depth=median(depth[ii]),mean_depth=mean(depth[ii]),median_score=median(scores$FinalScore[ii]))})),
  depths=lapply(idx,function(ii)as.numeric(depth[ii]))) }

a<-parse_args(commandArgs(trailingOnly=TRUE)); source(a$`pss-core`); require_package("ape"); require_package("geiger")
state<-readRDS(a$state); tree<-state$tree; observed_fits<-state$fits; observed_cov<-covs(tree,observed_fits)
generating<-state$selected_model; mu<-gls_mean(state$values,observed_cov[[generating]])
types<-trimws(strsplit(a$`bootstrap-types`,",",fixed=TRUE)[[1]]); if(!all(types%in%c("conditional","full"))) stop("Invalid bootstrap_types.")
start<-as.integer(a$`simulation-start`); reps<-as.integer(a$replicates); max_attempts<-as.integer(a$`max-attempts`); base_seed<-as.integer(a$seed)
rows<-list(); depths<-list(); z<-1L; RNGkind("L'Ecuyer-CMRG")
for(sim in seq.int(start,length.out=reps)) {
  success<-FALSE; last_error<-NULL
  for(attempt in seq_len(max_attempts)) {
    set.seed(as.integer((as.double(base_seed)+stable_hash(state$trait)+sim*1009+attempt*1000003)%%2147483000))
    candidate<-tryCatch({
      v<-simulate_trait(mu,observed_cov[[generating]]); out<-list()
      if("conditional"%in%types) out$conditional<-list(scores=calculate_pairwise_scores(v,tree,observed_cov$BM,observed_cov$OU,generating),selected=generating,aic_bm=NA,aic_ou=NA)
      if("full"%in%types) { f<-fit_models(tree,v); cv<-covs(tree,f); sm<-select_model(f); out$full<-list(scores=calculate_pairwise_scores(v,tree,cv$BM,cv$OU,sm),selected=sm,aic_bm=fit_aic(f$BM),aic_ou=fit_aic(f$OU)) }
      out
    },error=function(e){last_error<<-conditionMessage(e);NULL})
    if(!is.null(candidate)){success<-TRUE;break}
  }
  if(!success) stop("Simulation ",sim," failed: ",last_error)
  for(type in names(candidate)) { p<-tail_profile(candidate[[type]]$scores,candidate[[type]]$selected,state$tail_proportion); rr<-p$stats; rr$trait<-state$trait; rr$simulation<-sim; rr$bootstrap_type<-type; rr$generating_model<-generating; rr$selected_model<-candidate[[type]]$selected; rr$AIC_BM<-candidate[[type]]$aic_bm; rr$AIC_OU<-candidate[[type]]$aic_ou; rr$attempt<-attempt; rr<-rr[c("trait","simulation","bootstrap_type","generating_model","selected_model","AIC_BM","AIC_OU","attempt","tail","n_pairs","n_tail","median_depth","mean_depth","median_score")]; rows[[z]]<-rr; depths[[z]]<-p$depths; z<-z+1L }
}
object<-list(schema_version=1L,trait=state$trait,simulation_start=start,simulation_end=start+reps-1L,bootstrap_types=types,tail_proportion=state$tail_proportion,replicates=do.call(rbind,rows),depths=depths,observed_summary=state$observed_summary,observed_depths=state$observed_depths,fit_table=state$fit_table)
saveRDS(object,a$output,compress="xz")
