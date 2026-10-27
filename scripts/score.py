#!/usr/bin/env python3
"""Score SQANTI3 rules-filter runs against an answer key.

The answer key says, for every isoform, whether it matches a reference
transcript exactly (gffcompare class "="). You will never have this for your
own species: here it lets us see whether a filtering decision is good.

Usage:
    python scripts/score.py results/01_strict [results/02_rescue ...] --key data/celegans_chrI_answer_key.tsv
    python scripts/score.py results/01_strict --key data/celegans_chrI_answer_key.tsv --losses
"""
import argparse
import glob
import os
import re

import pandas as pd


def read_key(key_file):
    key = pd.read_csv(key_file, sep="\t")
    return set(key.loc[key["matches_reference"].astype(str) == "True", "isoform"])


def find_file(run_dir, suffix):
    hits = glob.glob(os.path.join(run_dir, f"*{suffix}"))
    if not hits:
        raise SystemExit(f"No *{suffix} file in {run_dir}. Did the SQANTI3 run finish?")
    return hits[0]


def score_run(run_dir, real):
    classif = pd.read_csv(find_file(run_dir, "_RulesFilter_classification.txt"), sep="\t",
                          usecols=["isoform", "filter_result"])
    kept = set(classif.loc[classif["filter_result"] == "Isoform", "isoform"])
    real_here = real & set(classif["isoform"])
    real_kept = len(kept & real_here)
    return {
        "Run": os.path.basename(os.path.normpath(run_dir)),
        "Isoforms": len(classif),
        "Kept": len(kept),
        "Real kept": real_kept,
        "Recovered (%)": round(100 * real_kept / len(real_here), 1),
        "Precision (%)": round(100 * real_kept / len(kept), 1) if kept else 0.0,
    }


def score(run_dirs, key_file):
    """Table with one row per run: how many real transcripts we kept, and how clean the result is."""
    real = read_key(key_file)
    return pd.DataFrame([score_run(d, real) for d in run_dirs])


def condition_name(reason):
    # "perc_A_downstream_TTS: 70.0 > 59" -> perc_A_downstream_TTS ; "NA value in CDS_length" -> CDS_length
    reason = reason.strip()
    na = re.match(r"NA value in (\S+)", reason)
    return na.group(1) if na else reason.split(":")[0]


def losses(run_dir, key_file):
    """For each condition: how many isoforms fail it, how many of those are real,
    and how many real ones fail ONLY that condition (the ones a rescue rule could bring back)."""
    real = read_key(key_file)
    reasons = pd.read_csv(find_file(run_dir, "_filtering_reasons.txt"), sep="\t")
    reasons["conditions"] = reasons["filter_reason"].fillna("").apply(
        lambda r: sorted({condition_name(x) for x in r.split(";") if x.strip()}))
    reasons["real"] = reasons["isoform"].isin(real)
    rows = []
    for cond in sorted({c for cs in reasons["conditions"] for c in cs}):
        fails = reasons[reasons["conditions"].apply(lambda cs: cond in cs)]
        only = fails[fails["conditions"].apply(len) == 1]
        rows.append({"Condition": cond,
                     "Isoforms failing": len(fails),
                     "Real among them": int(fails["real"].sum()),
                     "Real failing ONLY this": int(only["real"].sum())})
    return pd.DataFrame(rows).sort_values("Real failing ONLY this", ascending=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("runs", nargs="+", help="SQANTI3 filter output directories")
    parser.add_argument("--key", required=True, help="Answer key (isoform, matches_reference, ...)")
    parser.add_argument("--losses", action="store_true", help="Which conditions remove real transcripts (first run only)")
    args = parser.parse_args()

    pd.set_option("display.width", 200)
    print(score(args.runs, args.key).to_string(index=False))
    if args.losses:
        print(f"\nWhy were isoforms removed in {args.runs[0]}?")
        print(losses(args.runs[0], args.key).to_string(index=False))
