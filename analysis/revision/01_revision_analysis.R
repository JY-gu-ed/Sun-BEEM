if (!nzchar(Sys.getenv("SUN_BEEM_PROJECT_ROOT"))) stop("Set SUN_BEEM_PROJECT_ROOT to the authorised project directory before running this script.")
suppressPackageStartupMessages({library(data.table); library(survival)})
setDTthreads(1); options(mc.cores=1); set.seed(20261104)
if (!nzchar(Sys.getenv("SUN_BEEM_CODE_ROOT"))) stop("Set SUN_BEEM_CODE_ROOT or use tools/run_analysis.R.")
source(file.path(Sys.getenv("SUN_BEEM_CODE_ROOT"), "analysis/shared/followup.R"))

root <- file.path(Sys.getenv("SUN_BEEM_PROJECT_ROOT"), "")
rev <- file.path(root,'JID_R1_20261104'); out <- file.path(rev,'results')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
covs <- c('sex','age','bmi','deprivation_quintile','highest_edu','phys_activity','sleep_score','smoking_status','drinking')
comp <- c('bi_outdoor','bi_uv_protection','bi_swr2008','bi_solarium_use')
ys <- c('all_cause_mortality','cvd_mortality_narrow','non_skin_cancer_mortality')
d <- fread(file.path(root,'washed_data/4_2_2_ana_for_mortality.csv'))
reg <- fread(file.path(root,'washed_data/3_1_1_death_data.csv'),select=c('eid','death_date'))
stopifnot(!anyDuplicated(d$eid),!anyDuplicated(reg$eid))
d[,death_date:=reg$death_date[match(eid,reg$eid)]]
d[,all_cause_mortality:=as.integer(!is.na(death_date))]
d[,death_time:=as.numeric(fifelse(!is.na(death_date),as.Date(death_date),as.Date('2025-08-04'))-as.Date(date_assess))]
prep <- function(z){for(v in covs)z[[v]]<-factor(z[[v]]);z$deprivation_quintile<-ordered(z$deprivation_quintile,levels=1:5);z$BEEM<-factor(z$BEEM,levels=c('low','medium','high'));z}
d <- prep(d)
for(y in ys[-1]) d[[y]]<-as.integer(d[[y]]=='yes')
stopifnot(all(d$death_time>0),all(d$bi_uv_score==rowSums(d[,..comp])))
cox <- function(z,y,x,time='death_time',extra=character())coxph(as.formula(paste0('Surv(',time,',',y,')~',paste(c(x,covs,extra),collapse='+'))),z,ties='efron')
extract <- function(f,y,model,terms){s<-summary(f);k<-which(rownames(s$coefficients)%in%terms);data.table(outcome=y,model=model,term=rownames(s$coefficients)[k],beta=s$coefficients[k,'coef'],se=s$coefficients[k,'se(coef)'],HR=s$conf.int[k,'exp(coef)'],low=s$conf.int[k,'lower .95'],high=s$conf.int[k,'upper .95'],p=s$coefficients[k,'Pr(>|z|)'],N=f$n,events=f$nevent)}
write <- function(z,n)fwrite(z,file.path(out,paste0(n,'.csv')))
old <- fread(file.path(root,'elife_revision/results/targeted_cox_results.csv'))
checks <- list()
for(y in ys){f<-cox(d,y,'BEEM');a<-extract(f,y,'original_categories',c('BEEMmedium','BEEMhigh'));b<-old[outcome==y&model=='original_categories'];stopifnot(all(a$N==b$N),all(a$events==b$events),max(abs(a$HR-b$HR))<1e-8);checks[[y]]<-a;cat('BASELINE VERIFIED',y,'\n');flush.console()}
write(rbindlist(checks),'baseline_verification')
write(old[outcome%in%ys & (grepl('component|leave_out|original_score',model))],'reused_component_robustness')
write(d[,.(N=.N,all_cause_events=sum(all_cause_mortality),CVD_events=sum(cvd_mortality_narrow),non_skin_cancer_events=sum(non_skin_cancer_mortality)),by=BEEM],'cohort_counts')
write(rbindlist(lapply(c('swr2008_kJ','outdoor_avg','vitd'),function(v)d[,.(variable=v,N=sum(is.finite(get(v))),mean=mean(get(v),na.rm=TRUE),SD=sd(get(v),na.rm=TRUE),minimum=min(get(v),na.rm=TRUE),q05=quantile(get(v),.05,na.rm=TRUE),q25=quantile(get(v),.25,na.rm=TRUE),median=median(get(v),na.rm=TRUE),q75=quantile(get(v),.75,na.rm=TRUE),q95=quantile(get(v),.95,na.rm=TRUE),maximum=max(get(v),na.rm=TRUE)),by=BEEM])),'exposure_distributions')
cont <- list(); cuts <- list()
for(v in c('swr2008_kJ','outdoor_avg')){
  other<-setdiff(comp,if(v=='swr2008_kJ')'bi_swr2008' else 'bi_outdoor')
  q<-unique(quantile(d[[v]],c(0,.25,.5,.75,1),na.rm=TRUE));stopifnot(length(q)==5)
  cuts[[v]]<-data.table(variable=v,probability=c(0,.25,.5,.75,1),cutpoint=q)
  d$continuous_exposure<-d[[v]]/if(v=='swr2008_kJ')1000 else 1
  d$exposure_quartile<-cut(d[[v]],breaks=q,include.lowest=TRUE,labels=paste0('Q',1:4))
  for(y in ys){f<-cox(d,y,c('continuous_exposure',other));a<-extract(f,y,paste0(v,'_linear_mutually_adjusted'),'continuous_exposure');cont[[length(cont)+1]]<-a
    f<-cox(d,y,c('exposure_quartile',other));cont[[length(cont)+1]]<-extract(f,y,paste0(v,'_quartiles_mutually_adjusted'),grep('exposure_quartile',names(coef(f)),value=TRUE))}
}
write(rbindlist(cont),'continuous_quartile_sensitivity');write(rbindlist(cuts),'prespecified_empirical_quartiles')
cat('CONTINUOUS COMPLETE\n');flush.console()
spec<-data.table(file=c('4_2_3_ana_for_cvd_inpatient.csv','4_2_8_ana_for_non_skin_cancer_reg.csv','4_2_6_ana_for_melanoma_reg.csv','4_2_7_ana_for_other_skin_cancer_reg.csv'),outcome=c('Cardiovascular incidence','Non-skin cancer incidence','Melanoma incidence','C44 skin cancer incidence'),status=c('cvd_mortality_narrow','non_skin_cancer_status','melanoma_status','other_skin_cancer_status'),prior=c('diagnosed_after_assessment','non_skin_cancer_after_assessment','melanoma_after_assessment','other_skin_cancer_after_assessment'))
inc<-list();ia<-list()
for(i in seq_len(nrow(spec))){
  z<-fread(file.path(root,'washed_data',spec$file[i]));z<-z[is.na(get(spec$prior[i]))|get(spec$prior[i])!='no']
  z$death_time<-d$death_time[match(z$eid,d$eid)];stopifnot(!anyNA(z$death_time))
  dt<-if(i==1)'diag_date' else 'cancer_date';z[,diagnosis_days:=as.numeric(as.Date(get(dt))-as.Date(date_assess))]
  z[,legacy_event:=!is.na(get(spec$status[i]))&get(spec$status[i])=='yes']
  z<-as.data.table(apply_incidence_followup(as.data.frame(z),spec$status[i],dt,registry_end_date(if(i==1)'hospital' else 'cancer')))
  z[,event:=as.integer(get(spec$status[i])=='yes')];z[,time:=diag_time];z<-prep(z)
  ia[[i]]<-data.table(outcome=spec$outcome[i],N=nrow(z),events=sum(z$event),original_events=sum(z$legacy_event),postdeath_events=sum(z$legacy_event&z$diagnosis_days>z$death_time,na.rm=TRUE))
  for(v in comp){f<-cox(z,'event',v,'time');inc[[length(inc)+1]]<-extract(f,spec$outcome[i],paste0('component_',v),v)}
  f<-cox(z,'event',comp,'time');inc[[length(inc)+1]]<-extract(f,spec$outcome[i],'mutually_adjusted_components',comp)
  cat('INCIDENCE COMPLETE',spec$outcome[i],'\n');flush.console()
}
write(rbindlist(inc),'component_incidence');write(rbindlist(ia),'incidence_followup_audit')
r<-fread(file.path(root,'raw_data_3_8/request_an_columns.csv'),select=c('eid','p738_i0','p6142_i0','p20277_i0','p816_i0'))
stopifnot(!anyDuplicated(r$eid));z<-merge(d,r,by='eid',all.x=TRUE)
z[,income:=factor(fifelse(p738_i0%in%1:5,p738_i0,NA_real_),levels=1:5)]
valid_emp<-!is.na(z$p6142_i0)&!grepl('-3',z$p6142_i0,fixed=TRUE)
for(j in 1:7){z[[paste0('employment_',j)]]<-ifelse(valid_emp,as.integer(grepl(paste0('(^|[^0-9-])',j,'([^0-9]|$)'),z$p6142_i0)),NA_integer_)}
stopifnot(all(z$employment_7[z$p6142_i0=='[-7]' & !is.na(z$p6142_i0)]==0))
z[,occupation_SOC_major:=factor(fifelse(p20277_i0>=1000&p20277_i0<=9999,floor(p20277_i0/1000),NA_real_),levels=1:9)]
ev<-paste0('employment_',1:7)
dist<-list()
for(v in c('income',ev,'occupation_SOC_major','p816_i0')){a<-z[,.(N=.N),by=c('BEEM',v)];setnames(a,v,'level');a[,level:=as.character(level)];a[,variable:=v];a[,percent:=100*N/sum(N),by=BEEM];dist[[v]]<-a}
write(rbindlist(dist,fill=TRUE),'extra_confounder_distribution')
ext<-list();avail<-list()
for(label in c('income_employment','income_employment_occupation')){
  extra<-c('income',ev,if(label=='income_employment_occupation')'occupation_SOC_major')
  zz<-z[complete.cases(z[,c(covs,extra),with=FALSE])]
  avail[[label]]<-data.table(model=label,N=nrow(zz),excluded=nrow(z)-nrow(zz),occupation_available=sum(!is.na(z$occupation_SOC_major)))
  for(y in ys){f0<-cox(zz,y,'BEEM');f1<-cox(zz,y,'BEEM',extra=extra);ext[[length(ext)+1]]<-extract(f0,y,paste0(label,'_matched_original'),c('BEEMmedium','BEEMhigh'));ext[[length(ext)+1]]<-extract(f1,y,paste0(label,'_extended'),c('BEEMmedium','BEEMhigh'))}
  cat('EXTENDED COMPLETE',label,nrow(zz),'\n');flush.console()
}
write(rbindlist(ext),'extended_adjustment');write(rbindlist(avail),'extra_confounder_availability')
rm(z,r,zz);gc()
bmrows<-list();flow<-list()
sets<-list(vitamin_D=list(data=d,targets='vitd'),protein=list(path='step8/mea_set_merged.rds',targets=c('il6','fgf23','cxcl5','ccl5','pigr')),metabolite=list(path='step10/mea_set_merged.rds',targets=c('p23470_i0','p23480_i0')))
for(kind in names(sets)){
  item<-sets[[kind]];z<-if(kind=='vitamin_D')copy(d) else as.data.table(readRDS(file.path(root,item$path)));z<-prep(z)
  stopifnot(all(item$targets%in%names(z)),all(z$eid%in%d$eid))
  flow[[kind]]<-data.table(subset=kind,N=nrow(z),target=item$targets,nonmissing=sapply(item$targets,function(v)sum(is.finite(z[[v]]))),assessment_min=as.character(min(as.Date(z$date_assess),na.rm=TRUE)),assessment_max=as.character(max(as.Date(z$date_assess),na.rm=TRUE)))
  for(b in item$targets){
    stopifnot(is.numeric(z[[b]]))
    for(model in c('single_component','mutually_adjusted')){
      for(v in if(model=='single_component')comp else 'joint'){
        xx<-if(model=='single_component')v else comp
        f<-lm(as.formula(paste(b,'~',paste(c(xx,covs),collapse='+'))),data=z)
        s<-summary(f)$coefficients;ci<-confint(f);k<-which(rownames(s)%in%xx)
        bmrows[[length(bmrows)+1]]<-data.table(family=kind,biomarker=b,model=model,component=rownames(s)[k],beta=s[k,1],se=s[k,2],low=ci[k,1],high=ci[k,2],p=s[k,4],N=nobs(f),response_SD=sd(model.response(model.frame(f))))
      }
    }
  }
  cat('BIOMARKER COMPLETE',kind,nrow(z),'\n');flush.console();rm(z);gc()
}
bm<-rbindlist(bmrows);bm[,correction_family:=ifelse(family=='vitamin_D','vitamin_D_4_tests','selected_molecular_28_tests')]
bm[,p_BH:=p.adjust(p,method='BH'),by=.(model,correction_family)]
write(bm,'selected_biomarkers');write(rbindlist(flow),'subset_flow')
capture.output(sessionInfo(),file=file.path(out,'sessionInfo.txt'))
write(data.table(metric=c('main_N','allcause_events','CVD_events','non_skin_cancer_events','assessment_min','assessment_max'),value=c(nrow(d),sum(d$all_cause_mortality),sum(d$cvd_mortality_narrow),sum(d$non_skin_cancer_mortality),as.character(min(as.Date(d$date_assess))),as.character(max(as.Date(d$date_assess))))),'analysis_manifest')
cat('JID_REVISION_ANALYSES_COMPLETE\n')
