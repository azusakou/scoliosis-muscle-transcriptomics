#!/usr/bin/env Rscript
# Post hoc expression-proxy/rank and available-covariate diagnostics only.
# No BMI/fat measurements, authenticated clinical-ID merge or gate rescue.
source("scripts/lib_common.R")
out <- "analysis"
rd <- function(path) read.delim(path, check.names=FALSE)
mat <- function(path) { x<-rd(gzfile(path)); m<-as.matrix(x[,-1]); rownames(m)<-x[[1]]; m }

# Within-patient expression differences preserve the actual unit of replication.
A <- mat("results/tables/bulk/GSE254300_AIS_TMM_logCPM.tsv.gz")
B <- mat("results/tables/bulk_gse311180_requant/GSE311180_TMM_logCPM.tsv.gz")
J <- mat("results/tables/disease_control_gse254300_cs/GSE254300_AIS_CS_joint_TMM_logCPM.tsv.gz")
ma <- rd("metadata/audit/bulk_samples.tsv"); ma$sample_id<-ma$processed_sample_id
mb <- rd("metadata/audit/gse311180_requant_manifest.tsv")
effects <- list(); models <- list()
for (cohort in c("GSE254300_AIS","GSE311180_AIS","GSE254300_CS")) {
  if(cohort=="GSE311180_AIS") { m<-B; info<-mb } else {
    m<-if(cohort=="GSE254300_AIS") A else J
    info<-ma[ma$cohort==if(cohort=="GSE254300_AIS") "AIS" else "congenital_scoliosis",]
  }
  stopifnot(nrow(info)==if(cohort=="GSE311180_AIS") 14L else 10L)
  stopifnot(all(c("ADIPOQ","ANGPT1","LEP") %in% rownames(m)))
  for(pid in sort(unique(info$patient_id))) {
    c<-info$sample_id[info$patient_id==pid & info$side=="concave"]
    v<-info$sample_id[info$patient_id==pid & info$side=="convex"]
    stopifnot(length(c)==1,length(v)==1)
    effects[[length(effects)+1]]<-data.frame(cohort=cohort,patient_id=pid,
      ADIPOQ_difference=m["ADIPOQ",c]-m["ADIPOQ",v],
      ANGPT1_difference=m["ANGPT1",c]-m["ANGPT1",v],
      LEP_difference=m["LEP",c]-m["LEP",v])
  }
}
effects<-do.call(rbind,effects);stopifnot(nrow(effects)==17L,length(unique(effects$cohort))==3L)
for(cohort in unique(effects$cohort)) for(target in c("ANGPT1","LEP")) {
  d<-effects[effects$cohort==cohort,];d$target_difference<-d[[paste0(target,"_difference")]]
  fit<-lm(target_difference~ADIPOQ_difference,d);co<-summary(fit)$coefficients;df<-df.residual(fit)
  ci<-confint(fit)["(Intercept)",]
  models[[length(models)+1]]<-data.frame(cohort=cohort,target=target,n_pairs=nrow(d),
    unadjusted_mean_logCPM_difference=mean(d$target_difference),
    intercept_at_zero_ADIPOQ_difference=co[1,1],SE=co[1,2],residual_df=df,
    CI95_lower=ci[1],CI95_upper=ci[2],nominal_t_p=co[1,4],
    ADIPOQ_slope=co[2,1],pearson_paired_difference_correlation=cor(d$target_difference,d$ADIPOQ_difference),
    scope="Post hoc endogenous expression-proxy sensitivity; not BMI/fat adjustment, causal separation or original edgeR estimator")
}
models<-do.call(rbind,models);stopifnot(nrow(models)==6L)
models$BH_six_proxy_intercept_tests<-p.adjust(models$nominal_t_p,"BH")
write_tsv(effects,file.path(out,"bulk_followup_ADIPOQ_proxy_patient_effects.tsv"))
write_tsv(models,file.path(out,"bulk_followup_ADIPOQ_proxy_models.tsv"))

# Full composition-space rank check: no favorable subset selection/truncation.
props<-rd("results/tables/bulk/deconvolution_two_methods.tsv")
rankrows<-list();pairedrows<-list()
for(method in unique(props$method)) {
  z<-props[props$method==method,];types<-sort(unique(z$cell_type));patients<-sort(unique(z$patient_id))
  delta<-sapply(types,function(ct) sapply(patients,function(pid) {
    q<-z[z$patient_id==pid & z$cell_type==ct,];q$proportion[q$side=="concave"]-q$proportion[q$side=="convex"]
  }))
  colnames(delta)<-types;rownames(delta)<-patients
  variable<-apply(delta,2,sd)>1e-8
  nuis<-delta[,variable,drop=FALSE];rank<-qr(cbind(1,nuis),tol=1e-7)$rank
  nuisance_rank<-qr(nuis,tol=1e-7)$rank
  # More linearly independent nuisance directions than patients cannot provide
  # a residual variance or a non-arbitrary identifiable conditional side effect.
  rankrows[[length(rankrows)+1]]<-data.frame(cohort="GSE254300_AIS",method=method,n_pairs=length(patients),
    reference_lineages=length(types),variable_paired_components=sum(variable),
    augmented_full_composition_design_rank=rank,residual_df=length(patients)-rank,
    nuisance_design_rank=nuisance_rank,conditional_intercept_identifiable=rank>nuisance_rank,
    full_space_inference_available=length(patients)>rank,
    scope="Original atlas-derived logCPM NNLS fractions; shared genes and imperfect reference; rank diagnostic not tissue validation")
  pairedrows[[length(pairedrows)+1]]<-data.frame(method=method,patient_id=patients,delta,check.names=FALSE)
}
write_tsv(do.call(rbind,rankrows),file.path(out,"bulk_followup_full_composition_rank.tsv"))
write_tsv(do.call(rbind,pairedrows),file.path(out,"bulk_followup_composition_patient_differences.tsv"))

# Report known GSE311180 age/Cobb modifiers, without implying within-pair control
# of unrecorded laterality, tissue adiposity or handling order.
clinical<-rd("analysis/tables/audit_patient_clinical_metadata_19_pairs.tsv")
clinical<-clinical[clinical$dataset=="GSE311180",]
eff<-rd("results/tables/bulk_gse311180_requant/frozen_axis_paired_effects.tsv")
stopifnot(all(clinical$analysis_patient_id %in% eff$patient_id))
clinical$score_difference<-eff$paired_difference[match(clinical$analysis_patient_id,eff$patient_id)]
clinical$age_years<-as.numeric(clinical$age_years);clinical$cobb_angle_degrees<-as.numeric(clinical$cobb_angle_degrees)
modifier<-list()
for(label in c("all_seven_metadata_pairs","five_RNA_gate_passing_pairs")) {
  d<-if(label=="all_seven_metadata_pairs") clinical else clinical[!clinical$analysis_patient_id %in% c("STB-084","STB-100"),]
  d$age_centered<-d$age_years-mean(d$age_years);d$cobb_centered<-d$cobb_angle_degrees-mean(d$cobb_angle_degrees)
  fit<-lm(score_difference~age_centered+cobb_centered,d);co<-summary(fit)$coefficients;ci<-confint(fit)
  for(term in rownames(co)) modifier[[length(modifier)+1]]<-data.frame(subset=label,n_pairs=nrow(d),term=term,
    estimate=co[term,1],SE=co[term,2],residual_df=df.residual(fit),CI95_lower=ci[term,1],CI95_upper=ci[term,2],nominal_t_p=co[term,4],
    scope="Post hoc paired-effect modifier model; centered intercept is mean effect; no missing-covariate control or primary replacement")
}
write_tsv(clinical[,c("analysis_patient_id","age_years","cobb_angle_degrees","score_difference")],file.path(out,"bulk_followup_available_clinical_patient_effects.tsv"))
modifier<-do.call(rbind,modifier)
slopes<-modifier$term!="(Intercept)";modifier$BH_four_modifier_tests<-NA_real_
modifier$BH_four_modifier_tests[slopes]<-p.adjust(modifier$nominal_t_p[slopes],"BH")
write_tsv(modifier,file.path(out,"bulk_followup_available_clinical_modifiers.tsv"))
cat("Available-data follow-up sensitivities completed; missing clinical/anatomical evidence remains missing.\n")
