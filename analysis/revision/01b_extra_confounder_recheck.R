if (!nzchar(Sys.getenv("SUN_BEEM_PROJECT_ROOT"))) stop("Set SUN_BEEM_PROJECT_ROOT to the authorised project directory before running this script.")
suppressPackageStartupMessages({library(data.table); library(survival)})
setDTthreads(1); options(mc.cores=1); set.seed(20261104)
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

cat('EXTRA_CONFOUNDER_RECHECK_COMPLETE\n')
