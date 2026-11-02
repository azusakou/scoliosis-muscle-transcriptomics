#!/usr/bin/env python3
"""Resolution-aware patient-side tests for the two unchanged candidate genes."""
import csv
import itertools
import math
from pathlib import Path
from verify_results import bh, read_table

ROOT = Path(__file__).resolve().parents[1]


def generate(root=ROOT):
    rev = root / "analysis"
    records = read_table(rev / "bulk_gate_fixed_candidate_patient_effects.tsv")
    results = []
    for norm in ("locked_original_TMM", "retained_library_reTMM"):
        for gene in ("ANGPT1", "LEP"):
            rows = [r for r in records if r["normalization"] == norm and r["gene"] == gene]
            if len(rows) != 5 or len({r["patient_id"] for r in rows}) != 5: raise ValueError("Candidate pairs")
            values = [float(r["paired_logCPM_difference"]) for r in rows]
            sums = [math.fsum(x*s for x,s in zip(values,signs)) for signs in itertools.product((-1,1),repeat=5)]
            p = sum(abs(v) >= abs(math.fsum(values))-1e-12 for v in sums)/32
            results.append(dict(normalization=norm,gene=gene,n_pairs=5,positive_pairs=sum(v>0 for v in values),
                                paired_logCPM_mean=math.fsum(values)/5,exact_two_sided_p=p,
                                minimum_attainable_two_sided_p=2/32,
                                role="retrospective_patient_level_resolution_check_not_QL_FDR"))
    for row, adj in zip(results,bh([r["exact_two_sided_p"] for r in results])):
        row["BH_four_candidate_normalization_family"] = adj
    with (rev/"bulk_gate_candidate_patient_exact.tsv").open("w",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=list(results[0]),delimiter="\t");writer.writeheader();writer.writerows(results)
    return results


if __name__ == "__main__":
    generate()
