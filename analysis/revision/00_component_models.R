if (!nzchar(Sys.getenv("SUN_BEEM_PROJECT_ROOT"))) stop("Set SUN_BEEM_PROJECT_ROOT to the project directory before running this script.")
suppressPackageStartupMessages({library(data.table);library(survival)})
root <- Sys.getenv("SUN_BEEM_PROJECT_ROOT")
out <- file.path(root,'elife_revision','results'); dir.create(out,recursive=TRUE,showWarnings=FALSE)
setDTthreads(1); options(mc.cores=1); set.seed(20261012)
d <- fread(file.path(root,'washed_data','4_2_2_ana_for_mortality.csv'))
registered <- fread(file.path(root,'washed_data/3_1_1_death_data.csv'),select=c('eid','death_date'))
stopifnot(!anyDuplicated(registered$eid));d[,death_date:=registered$death_date[match(eid,registered$eid)]]
d[,all_cause_mortality:=ifelse(!is.na(death_date),'yes','no')]
d[,death_time:=as.numeric(fifelse(!is.na(death_date),as.Date(death_date),as.Date('2025-08-04'))-as.Date(date_assess))]
saveRDS(d[,.(eid,death_time,all_cause_mortality)],file.path(root,'elife_revision/corrected_death_times.rds'))
covs <- c('sex','age','bmi','deprivation_quintile','highest_edu','phys_activity','sleep_score','smoking_status','drinking')
comps <- c('bi_outdoor','bi_uv_protection','bi_swr2008','bi_solarium_use')
outs <- c('all_cause_mortality','cvd_mortality_narrow','non_skin_cancer_mortality')
for(v in covs) d[[v]] <- factor(d[[v]])
d$deprivation_quintile <- ordered(d$deprivation_quintile,levels=1:5)
d$BEEM <- factor(d$BEEM,levels=c('low','medium','high'));d$BEEM_num <- as.numeric(d$BEEM)-1
allouts <- c(outs,'cancer_mortality','melanoma_skin_cancer_mortality','non_melanoma_mortality')
for(v in allouts) d[[v]] <- if(is.character(d[[v]])) as.integer(d[[v]]=='yes') else as.integer(d[[v]])
stopifnot(!anyDuplicated(d$eid),all(d$death_time>0),all(d$bi_uv_score==rowSums(d[,..comps])))
fwrite(data.table(metric=c('N','followup_median_years','followup_min_years','followup_max_years','missing_covariate_cells','score_mismatches'),value=c(nrow(d),median(d$death_time)/365.25,min(d$death_time)/365.25,max(d$death_time)/365.25,sum(is.na(d[,..covs])),sum(d$bi_uv_score!=rowSums(d[,..comps])))),file.path(out,'cohort_audit.csv'))
fwrite(d[,.(N=.N,events=sum(all_cause_mortality)),by=BEEM],file.path(out,'category_counts.csv'))
fwrite(d[,.(N=.N),by=bi_uv_score],file.path(out,'score_counts.csv'))
delta <- as.numeric(as.Date(d$date_assess))+d$death_time
fwrite(data.table(status=c('alive','dead'),min_end_date=as.character(as.Date(c(min(delta[d$all_cause_mortality==0]),min(delta[d$all_cause_mortality==1])),origin='1970-01-01')),max_end_date=as.character(as.Date(c(max(delta[d$all_cause_mortality==0]),max(delta[d$all_cause_mortality==1])),origin='1970-01-01'))),file.path(out,'followup_dates.csv'))
fwrite(data.table(component=comps,N_high=sapply(comps,function(v)sum(d[[v]]==1)),N_missing=sapply(comps,function(v)sum(is.na(d[[v]])))),file.path(out,'component_counts.csv'))
getfit <- function(y,x,data=d) coxph(as.formula(paste('Surv(death_time,',y,') ~',x,'+',paste(covs,collapse='+'))),data=data,ties='efron',x=TRUE)
rows <- list(); phrows <- list();fits <- list()
extract <- function(f,y,label,terms){s<-summary(f); k<-rownames(s$coefficients)%in%terms;data.table(outcome=y,model=label,term=rownames(s$coefficients)[k],beta=s$coefficients[k,'coef'],se=s$coefficients[k,'se(coef)'],HR=s$conf.int[k,'exp(coef)'],low=s$conf.int[k,'lower .95'],high=s$conf.int[k,'upper .95'],p=s$coefficients[k,'Pr(>|z|)'],N=f$n,events=f$nevent)}
for(y in allouts){
 f<-getfit(y,'BEEM');fits[[y]]<-f;rows[[length(rows)+1]]<-extract(f,y,'original_categories',c('BEEMmedium','BEEMhigh'))
 if(y%in%outs){ph<-cox.zph(f);phrows[[length(phrows)+1]]<-data.table(outcome=y,term=rownames(ph$table),ph$table)
  for(v in comps){f1<-getfit(y,v);rows[[length(rows)+1]]<-extract(f1,y,paste0('component_',v),v)}
  f1<-getfit(y,paste(comps,collapse='+'));rows[[length(rows)+1]]<-extract(f1,y,'mutually_adjusted_components',comps)
  for(v in comps){d$reduced_score<-d$bi_uv_score-d[[v]];f1<-getfit(y,'reduced_score');rows[[length(rows)+1]]<-extract(f1,y,paste0('leave_out_',v),'reduced_score')}
  f1<-getfit(y,'bi_uv_score');rows[[length(rows)+1]]<-extract(f1,y,'original_score_per_point','bi_uv_score')
  d$score_factor<-factor(d$bi_uv_score);f1<-getfit(y,'score_factor');rows[[length(rows)+1]]<-extract(f1,y,'original_score_categorical',grep('score_factor',names(coef(f1)),value=TRUE))
 }
 cat('completed',y,'\n');flush.console()
}
fwrite(rbindlist(rows),file.path(out,'targeted_cox_results.csv'));fwrite(rbindlist(phrows),file.path(out,'ph_tests.csv'))
