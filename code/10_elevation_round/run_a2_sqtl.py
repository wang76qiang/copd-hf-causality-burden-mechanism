#!/usr/bin/env python3
"""A2 extension v2: screen SERPINE1 (ENSG00000106366) in GTEx v8 exon/tx/txrev/leafcutter
permuted top-variant files (only these are distributed for non-ge quantifications)."""
import csv, os, subprocess, gzip, io

BASE = "https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000015"
DATASETS = {
    "lung":            {"exon": "QTD000272", "tx": "QTD000273", "txrev": "QTD000274", "leafcutter": "QTD000275"},
    "heart_LV":        {"exon": "QTD000257", "tx": "QTD000258", "txrev": "QTD000259", "leafcutter": "QTD000260"},
    "heart_AA":        {"exon": "QTD000252", "tx": "QTD000253", "txrev": "QTD000254", "leafcutter": "QTD000255"},
    "artery_aorta":    {"exon": "QTD000132", "tx": "QTD000133", "txrev": "QTD000134", "leafcutter": "QTD000135"},
    "artery_coronary": {"exon": "QTD000137", "tx": "QTD000138", "txrev": "QTD000139", "leafcutter": "QTD000140"},
}
os.makedirs("eqtl_extended", exist_ok=True)
os.makedirs("results_new", exist_ok=True)
summary = []
for tissue, quants in DATASETS.items():
    for qm, qtd in quants.items():
        url = f"{BASE}/{qtd}/{qtd}.permuted.tsv.gz"
        out = f"eqtl_extended/{tissue}_{qm}_{qtd}_permuted.tsv"
        if not os.path.exists(out):
            subprocess.run(["curl", "-s", "--retry", "3", "-o", out + ".gz", url], check=True, timeout=300)
            try:
                with gzip.open(out + ".gz", "rt", encoding="utf-8", errors="replace") as fh:
                    txt = fh.read()
            except Exception as e:
                print(f"{tissue}/{qm} {qtd}: download/gunzip failed: {e}")
                os.path.exists(out + ".gz") and os.remove(out + ".gz")
                continue
            open(out, "w", encoding="utf-8").write(txt)
            os.remove(out + ".gz")
        rows = list(csv.DictReader(open(out, encoding="utf-8"), delimiter="\t"))
        if not rows:
            print(f"{tissue}/{qm} {qtd}: empty"); continue
        cols = list(rows[0].keys())
        serp = [r for r in rows if "ENSG00000106366" in (r.get("gene_id","") or "") + (r.get("molecular_trait_id","") or "")]
        if tissue == "lung" and qm == "exon":
            print("cols:", cols)
        pcol = "qval" if "qval" in cols else ("pval_beta" if "pval_beta" in cols else "pvalue")
        best = min(serp, key=lambda r: float(r[pcol])) if serp else None
        summary.append({"tissue": tissue, "quant": qm, "qtd": qtd,
                        "n_traits_total": len(rows), "n_serpine1_traits": len(serp),
                        "pcol": pcol,
                        "best_trait": best.get("molecular_trait_id") if best else None,
                        "best_variant": best.get("rsid") or best.get("variant") if best else None,
                        "best_p": best.get(pcol) if best else None,
                        "best_beta": best.get("beta") if best else None,
                        "significant_005": (float(best[pcol]) < 0.05) if best else None})
        if best:
            print(f"{tissue}/{qm} ({qtd}): SERPINE1 traits={len(serp)}, best {pcol}={best[pcol]} ({best.get('rsid') or best.get('variant')})")
with open("results_new/a2_sqtl_screen_summary.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(summary[0].keys())); w.writeheader(); w.writerows(summary)
print("DONE")
