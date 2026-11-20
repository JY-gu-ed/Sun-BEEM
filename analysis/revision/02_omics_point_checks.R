if (!nzchar(Sys.getenv("SUN_BEEM_PROJECT_ROOT"))) stop("Set SUN_BEEM_PROJECT_ROOT to the project directory before running this script.")
suppressPackageStartupMessages({library(data.table);library(survival)})
root<-Sys.getenv("SUN_BEEM_PROJECT_ROOT");out<-file.path(root,'elife_revision','results')
dir.create(out,recursive=TRUE,showWarnings=FALSE);setDTthreads(1)
covs<-c('age','sex','bmi','deprivation_quintile','highest_edu','smoking_status','phys_activity','sleep_score','drinking')
meta<-list();checks<-list();summ<-list()
for(layer in c('step8','step10')){
 d<-readRDS(file.path(root,layer,'mea_set_merged.rds'));bms<-readRDS(file.path(root,layer,if(layer=='step8')'biomarkers_olink.rds' else 'biomarkers_nmr.rds'))
 corrected<-readRDS(file.path(root,'elife_revision/corrected_death_times.rds'));d$death_time<-corrected$death_time[match(d$eid,corrected$eid)];stopifnot(!anyNA(d$death_time))
 d$BEEM_num<-match(d$BEEM,c('low','medium','high'))-1
 meta[[layer]]<-data.table(layer=layer,N=nrow(d),biomarkers=length(bms),unique_ids=length(unique(d$eid)),missing_covariates=sum(is.na(d[,covs])),cvd_events=sum(d$cvd_mortality_narrow==1,na.rm=TRUE),cancer_events=sum(d$cancer_mortality==1,na.rm=TRUE))
 selected<-if(layer=='step8')intersect(c('acp1','a1bg','il6','fgf23','pigr','cxcl5','ccl5'),bms) else intersect(c('p23400_i0','p23470_i0','p23480_i0','p23471_i0'),bms)
 for(v in selected){
  summ[[length(summ)+1]]<-data.table(layer=layer,biomarker=v,N_nonmissing=sum(!is.na(d[[v]])),mean=mean(d[[v]],na.rm=TRUE),sd=sd(d[[v]],na.rm=TRUE),min=min(d[[v]],na.rm=TRUE),max=max(d[[v]],na.rm=TRUE))
  for(y in c('cvd_mortality_narrow','cancer_mortality')){
   cols<-c('death_time',y,'BEEM_num',covs,v);dd<-d[complete.cases(d[,cols]),cols]
   if(is.character(dd[[y]]))dd[[y]]<-as.integer(dd[[y]]=='yes')
   f0<-as.formula(paste('Surv(death_time,',y,')~BEEM_num+',paste(covs,collapse='+')))
   f1<-update(f0,paste('~.+',v));m0<-coxph(f0,dd);m1<-coxph(f1,dd)
   a<-coef(m0)['BEEM_num'];b<-coef(m1)['BEEM_num']
   checks[[length(checks)+1]]<-data.table(layer=layer,outcome=y,biomarker=v,N=nrow(dd),events=m0$nevent,beta_without=a,beta_with=b,coefficient_change=a-b,relative_coefficient_change=(a-b)/a)
  }
 }
 cat('completed',layer,'\n');flush.console();rm(d);gc()
}
fwrite(rbindlist(meta),file.path(out,'omics_object_audit.csv'));fwrite(rbindlist(summ),file.path(out,'omics_scale_checks.csv'));fwrite(rbindlist(checks),file.path(out,'omics_point_reproduction.csv'))
cat('OMICS_AUDIT_COMPLETE\n')
