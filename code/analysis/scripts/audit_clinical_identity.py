#!/usr/bin/env python3
"""Retained public metadata and expression plausibility audit; no imputation."""
import csv, gzip, hashlib, json, math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "analysis/tables"
OUT.mkdir(parents=True, exist_ok=True)
def read(path):
    with (gzip.open(path,"rt") if str(path).endswith(".gz") else open(path)) as f:
        return list(csv.DictReader(f,delimiter="\t"))
def write(name, rows):
    with open(OUT/name,"w") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter="\t",lineterminator="\n")
        w.writeheader();w.writerows(rows)
NR="not reported in retained patient-linked public metadata"
fields=["dataset","analysis_patient_id","diagnosis","age_years","sex","cobb_angle_degrees",
        "curve_classification","curve_direction_anatomical","sampled_vertebral_level_concave",
        "sampled_vertebral_level_convex","BMI","fatty_infiltration","sampling_order",
        "operation_side","handedness","clinical_link_status","source"]
rows=[]
sc=read(ROOT/"metadata/audit/scrna_samples.tsv")
for p in dict.fromkeys(x["patient_id"] for x in sc):
    x=next(x for x in sc if x["patient_id"]==p)
    z=dict.fromkeys(fields,NR);z.update(dataset="PRJNA722100",analysis_patient_id=p,diagnosis="AIS",
    age_years=x["age"],sex=x["sex"],clinical_link_status="patient ID and side linked in retained SRA metadata",
    source="https://www.ncbi.nlm.nih.gov/bioproject/PRJNA722100")
    rows.append(z)
bulk=read(ROOT/"metadata/audit/bulk_samples.tsv")
for p in dict.fromkeys(x["patient_id"] for x in bulk):
    x=next(x for x in bulk if x["patient_id"]==p)
    z=dict.fromkeys(fields,NR);z.update(dataset="GSE254300",analysis_patient_id=p,
    diagnosis="AIS" if x["cohort"]=="AIS" else "congenital scoliosis",
    clinical_link_status="clinical tables exist in source article but published case IDs versus GEO submitted aliases are unresolved; no patient-level merge",
    source="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE254300; https://pubmed.ncbi.nlm.nih.gov/38440359/")
    rows.append(z)
re=read(ROOT/"metadata/audit/gse311180_samples.tsv")
paired=[x for x in re if x["complete_pair"]=="TRUE" and x["tissue"]=="Muscle" and x["diagnosis"]=="AIS"]
for p in dict.fromkeys(x["patient_id"] for x in paired):
    x=next(x for x in paired if x["patient_id"]==p and x["side"]=="concave")
    y=next(x for x in paired if x["patient_id"]==p and x["side"]=="convex")
    z=dict.fromkeys(fields,NR);z.update(dataset="GSE311180",analysis_patient_id=p,diagnosis="AIS",
    age_years=x["age"],sex=x["sex"],cobb_angle_degrees=x["cobb_angle"],
    sampled_vertebral_level_concave=x["spinal_level"],sampled_vertebral_level_convex=y["spinal_level"],
    clinical_link_status="retained supplement-linked metadata; separate RNA consistency assessment; anatomical/clinical identity not independently authenticated",
    source="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE311180")
    rows.append(z)
assert len(rows)==19
write("audit_patient_clinical_metadata_19_pairs.tsv",rows)

aliases=[]
for p in dict.fromkeys(x["patient_id"] for x in bulk):
    xx=[x for x in bulk if x["patient_id"]==p]
    aliases.append({"analysis_patient_id":p,"diagnosis":xx[0]["cohort"],
    "GEO_concave_title":next(x["submitted_title"] for x in xx if x["side"]=="concave"),
    "GEO_convex_title":next(x["submitted_title"] for x in xx if x["side"]=="convex"),
    "published_RNAseq_case_ids":"IS1–IS5 or CS1–CS5 (source article Tables 1 and 2)",
    "patient_clinical_merge":"unresolved; matching processed index alone is insufficient to reconcile submitted titles"})
write("audit_gse254300_alias_clinical_link.tsv",aliases)

expr=read(ROOT/"results/tables/bulk_gse311180_requant/GSE311180_TMM_logCPM.tsv.gz")
markers={r["gene"]:r for r in expr if r["gene"] in ("XIST","RPS4Y1","DDX3Y","KDM5D","UTY","EIF1AY")}
sexrows=[]
for x in paired:
    z={"sample_id":x["sample_id"],"patient_id":x["patient_id"],"metadata_side":x["side"],"metadata_sex":x["sex"]}
    for gene in ("XIST","RPS4Y1","DDX3Y","KDM5D","UTY","EIF1AY"):
        z[gene+"_TMM_log2CPM"]=markers.get(gene,{}).get(x["sample_id"],"not retained/detectable")
    z["interpretation_scope"]="sex-associated expression plausibility only; cannot prove unique patient identity or concave/convex anatomical side"
    sexrows.append(z)
write("audit_gse311180_sex_marker_expression.tsv",sexrows)

# Pair-level concordance of an orthogonal sex marker is qualitative; no inferred genotype.
pairrows=[]
for p in dict.fromkeys(x["patient_id"] for x in paired):
    xx=[x for x in sexrows if x["patient_id"]==p]
    z={"patient_id":p,"metadata_sex":xx[0]["metadata_sex"],"concave_sample":next(x["sample_id"] for x in xx if x["metadata_side"]=="concave"),"convex_sample":next(x["sample_id"] for x in xx if x["metadata_side"]=="convex")}
    for gene in ("XIST","RPS4Y1","DDX3Y","KDM5D","UTY","EIF1AY"):
        cc=next(x[gene+"_TMM_log2CPM"] for x in xx if x["metadata_side"]=="concave")
        cv=next(x[gene+"_TMM_log2CPM"] for x in xx if x["metadata_side"]=="convex")
        try:z[gene+"_paired_log2CPM_difference"]=float(cc)-float(cv)
        except ValueError:z[gene+"_paired_log2CPM_difference"]="not retained/detectable"
    pairrows.append(z)
write("audit_gse311180_pair_sex_marker_differences.tsv",pairrows)

paths=["scripts/32_freeze_multicellular_signature.R","run.py","results/tables/bulk/frozen_multicellular_signature.tsv",
       "results/tables/bulk/frozen_multicellular_signature_manifest.tsv","results/tables/bulk/frozen_multicellular_signature_direction_conflicts.tsv",
       "scripts/60_bulk_validation.R","scripts/64_gse311180_requant_validation.R"]
hashes=[]
for p in paths:
    b=(ROOT/p).read_bytes()
    hashes.append({"path":p,"sha256_at_import":hashlib.sha256(b).hexdigest(),
                   "provenance_scope":"current file identity; does not timestamp prospective analytic decisions"})
write("audit_current_file_sha256.tsv",hashes)
print(json.dumps({"clinical_pairs":len(rows),"sex_marker_samples":len(sexrows),"retained_sex_markers":list(markers)},indent=2))
