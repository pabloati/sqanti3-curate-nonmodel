"""Builds rules_filter_activity.ipynb (the Colab hands-on for the 29 Oct 2026 seminar).

Edit the cell texts here and re-run `python3 build/build_notebook_A.py`; don't edit the .ipynb by hand.
"""
import json
from pathlib import Path

DATA_URL = "https://raw.githubusercontent.com/pabloati/sqanti3-curate-nonmodel/main/data"
COLAB_URL = "https://colab.research.google.com/github/pabloati/sqanti3-curate-nonmodel/blob/main/rules_filter_activity.ipynb"

cells = []


def md(text):
    cells.append({"cell_type": "markdown", "metadata": {}, "source": text.strip("\n")})


def code(text, form=False):
    meta = {"cellView": "form"} if form else {}
    cells.append({"cell_type": "code", "metadata": meta, "execution_count": None,
                  "outputs": [], "source": text.strip("\n")})


# ---------------------------------------------------------------- intro
md("""
# 🧹 Curating a long-read transcriptome with (almost) no extra evidence
**Hands-on: the SQANTI3 rules filter** · Pablo Atienza López, ConesaLab (I2SysBio, CSIC)

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](__COLAB_URL__)

You sequenced the transcriptome of your favourite non-model species with long reads, built transcript models (here with IsoQuant) and ran SQANTI3 quality control. Now you have a table with thousands of isoforms, and many of them are artifacts. **Which ones do you keep?**

In a model organism you would lean on the reference annotation, short reads, CAGE or polyA data. In most biodiversity projects you have **none of these**. In this activity we see how far the SQANTI3 **rules filter** gets with only what the long reads tell you, and how **combining several rules** recovers real transcripts that a single strict rule throws away.

### Today's trick: an answer key
We use *Caenorhabditis elegans* (PacBio HiFi, 20,417 isoforms). SQANTI3 was run against a *placebo* annotation, so it knows nothing about the worm's genes, just like in a newly sequenced species. But *C. elegans* has one of the best annotations there is (WormBase), so we keep an **answer key**: for each isoform, whether its coding sequence matches a WormBase transcript exactly. You will never have this key for your own species; here it lets us *see* whether our filtering decisions are good.

### How to use this notebook
Run the cells in order with ▶ (or `Shift + Enter`). You only need to edit code in **Step 8**. If something breaks, `Runtime → Restart and run all`.
""".replace("__COLAB_URL__", COLAB_URL))

# ---------------------------------------------------------------- setup
md("## Step 0 · Set-up (≈ 20 s)")
code(f'''
# @title ▶ Run me: load the data and the helper functions
import json
import pandas as pd
import matplotlib.pyplot as plt
from IPython.display import display, Markdown

DATA_URL = "{DATA_URL}"

classif = pd.read_csv(f"{{DATA_URL}}/celegans_hifi_classification.txt.gz", sep="\\t", dtype={{"chrom": str}})
answer_key = pd.read_csv(f"{{DATA_URL}}/celegans_answer_key.tsv.gz", sep="\\t")
classif = classif.merge(answer_key, on="isoform", how="left")
real = classif["matches_WormBase"]
N_ISO, N_REAL = len(classif), int(real.sum())

BLUE, ORANGE, GREY, INK = "#2a78d6", "#eb6834", "#a3a29c", "#52514e"
plt.rcParams.update({{"axes.spines.top": False, "axes.spines.right": False, "axes.edgecolor": GREY,
                     "axes.labelcolor": INK, "xtick.color": INK, "ytick.color": INK,
                     "axes.titleweight": "bold", "axes.titlesize": 12, "font.size": 10}})

# The SQANTI3 rules filter, rewritten to work on the whole table at once.
# It makes the same decisions as SQANTI3 (src/utilities/filter/sqanti3_rules_filter.py),
# but takes a fraction of a second instead of a minute. Step 9 lets you check this.
JUNCTION_COLS = {{"RTS_stage", "all_canonical", "min_cov", "min_cov_pos", "sd_cov",
                 "min_intron_length", "n_indels_junc", "bite", "predicted_NMD"}}

def _condition(df, col, rule):
    v = df[col]
    if isinstance(rule, list) and all(isinstance(x, (int, float)) for x in rule):
        ok = v.between(min(rule), max(rule))                           # [min, max] range
    elif isinstance(rule, list):
        ok = v.astype(str).str.lower().isin([str(x).lower() for x in rule])  # one of several values
    elif isinstance(rule, (int, float)):
        ok = v >= rule                                                 # single number = minimum
    else:
        ok = v.astype(str).str.lower() == str(rule).lower()            # text must match
    ok = ok & v.notna()                                                # NA fails the condition
    if col in JUNCTION_COLS:
        ok = ok | (df["exons"] == 1)                                   # mono-exons skip junction rules
    return ok

def sqanti_rules(df, rules):
    """True for isoforms kept by a SQANTI3 rules JSON (as a Python dict)."""
    cat = df["structural_category"].where(df["structural_category"].isin(rules.keys()), "rest")
    keep = pd.Series(False, index=df.index)
    for category, blocks in rules.items():
        for block in blocks:                       # OR between blocks
            ok = cat == category
            for col, rule in block.items():        # AND inside a block
                ok &= _condition(df, col, rule)
            keep |= ok
    return keep

scoreboard = []

def score(rules, name, show=True):
    keep = sqanti_rules(classif, rules)
    kept, good = int(keep.sum()), int((keep & real).sum())
    row = dict(rule_set=name, kept=kept, wormbase_kept=good, other_kept=kept - good,
               recovered=100 * good / N_REAL, precision=100 * good / kept if kept else 0)
    scoreboard[:] = [r for r in scoreboard if r["rule_set"] != name] + [row]
    if show:
        display(Markdown(
            f"**{{name}}** keeps **{{kept:,}}** of {{N_ISO:,}} isoforms. "
            f"**{{good:,}} of the {{N_REAL:,}} WormBase transcripts survive ({{row['recovered']:.0f}}% recovered)**, "
            f"and {{kept - good:,}} kept isoforms don't match WormBase. "
            f"**Precision: {{row['precision']:.0f}}%** of what you keep is a WormBase transcript."))
    return keep

def plot_tradeoff(curve=None, highlight=None):
    nofilter = [r for r in scoreboard if r["rule_set"] == "No filter"]
    sb = pd.DataFrame([r for r in scoreboard if r["rule_set"] != "No filter"])
    highlight = highlight or sb["rule_set"].iloc[-1]
    fig, ax = plt.subplots(figsize=(7, 4))
    if curve is not None:
        ax.plot(curve["recovered"], curve["precision"], color=BLUE, lw=2, alpha=0.35, zorder=1)
    key = []
    for i, r in sb.iterrows():
        on = r["rule_set"] == highlight
        ax.scatter(r["recovered"], r["precision"], s=230, color=BLUE if on else GREY,
                   edgecolor="white", linewidth=2, zorder=5 if on else 3)
        ax.text(r["recovered"], r["precision"], str(i + 1), ha="center", va="center",
                fontsize=8, color="white", fontweight="bold", zorder=6 if on else 4)
        key.append((f"{{i + 1}}  {{r['rule_set']}}:  {{r['recovered']:.1f}}% · {{r['precision']:.1f}}%", on))
    for j, (line, on) in enumerate(key):
        ax.text(1.03, 0.97 - j * 0.085, line, transform=ax.transAxes, va="top", fontsize=9,
                color="black" if on else INK, fontweight="bold" if on else "normal")
    note = "(recovered · precision)"
    if nofilter:
        note += f"\\nNot shown: no filter, {{nofilter[0]['recovered']:.0f}}% · {{nofilter[0]['precision']:.0f}}%"
    ax.text(1.03, 0.97 - len(key) * 0.085 - 0.03, note, transform=ax.transAxes,
            va="top", fontsize=8, color=INK, style="italic", linespacing=1.6)
    lo_x, lo_y, hi_y = sb["recovered"].min(), sb["precision"].min(), sb["precision"].max()
    if curve is not None:
        lo_x, lo_y, hi_y = min(lo_x, curve["recovered"].min()), min(lo_y, curve["precision"].min()), max(hi_y, curve["precision"].max())
    ax.set_xlim(lo_x - 2, 100.5); ax.set_ylim(lo_y - 1.5, hi_y + 1.5)
    ax.set_xlabel("WormBase transcripts recovered (%)  →  more complete")
    ax.set_ylabel("Precision (%)  →  cleaner")
    ax.set_title("Each dot is a rule set: up and to the right is better", loc="left")
    ax.grid(color="#e6e5e0", lw=0.8); ax.set_axisbelow(True)
    plt.show()

def explain_losses(block):
    """For each condition of one rule block: isoforms that fail ONLY that condition."""
    full = sqanti_rules(classif, {{"rest": [block]}})
    rows = []
    for col in block:
        if col == "exons":
            continue
        without = {{k: v for k, v in block.items() if k != col}}
        lost = sqanti_rules(classif, {{"rest": [without]}}) & ~full
        rows.append(dict(condition=col, removed_only_by_this=int(lost.sum()),
                         of_which_WormBase=int((lost & real).sum())))
    t = pd.DataFrame(rows).sort_values("removed_only_by_this")
    fig, ax = plt.subplots(figsize=(7, 3.2))
    other = t["removed_only_by_this"] - t["of_which_WormBase"]
    ax.barh(t["condition"], t["of_which_WormBase"], color=BLUE, label="matches WormBase (lost real transcripts)")
    ax.barh(t["condition"], other, left=t["of_which_WormBase"] + 30, color=ORANGE, label="doesn't match WormBase")
    for y, (n, w) in enumerate(zip(t["removed_only_by_this"], t["of_which_WormBase"])):
        ax.text(n + 60, y, f"{{n:,}} removed · {{w:,}} real", va="center", fontsize=9, color=INK)
    ax.set_xlim(0, t["removed_only_by_this"].max() * 1.45)
    ax.set_xlabel("Isoforms removed by this condition alone")
    ax.set_title("Which condition throws away real transcripts?", loc="left")
    ax.legend(frameon=False, loc="lower right", fontsize=9)
    plt.show()
    return t.sort_values("removed_only_by_this", ascending=False).reset_index(drop=True)

print(f"Loaded {{N_ISO:,}} isoforms; {{N_REAL:,}} of them match a WormBase transcript (our answer key).")
''', form=True)

# ---------------------------------------------------------------- step 1
md("""
## Step 1 · What do we actually know about these isoforms?
SQANTI3 describes each isoform with ~50 columns, but many of them need extra data you may not have. Let's see which ones are filled in.
""")
code('''
evidence = {
    "Short-read junction coverage (min_cov)": "min_cov",
    "Short-read 5′ coverage ratio (ratio_TSS)": "ratio_TSS",
    "CAGE peaks (dist_to_CAGE_peak)": "dist_to_CAGE_peak",
    "PolyA site data (dist_to_polyA_site)": "dist_to_polyA_site",
    "PolyA motif (polyA_motif)": "polyA_motif",
    "Short-read expression (iso_exp)": "iso_exp",
    "Long-read count (FL)": "FL",
    "% A after the 3′ end (perc_A_downstream_TTS)": "perc_A_downstream_TTS",
    "RT-switching (RTS_stage)": "RTS_stage",
    "Coding potential (coding)": "coding",
    "PSAURON protein score (psauron_score)": "psauron_score",
    "Canonical junctions (all_canonical)": "all_canonical",
    "NMD prediction (predicted_NMD)": "predicted_NMD",
}
filled = pd.Series({name: 100 * classif[col].notna().mean() for name, col in evidence.items()}).sort_values()
fig, ax = plt.subplots(figsize=(7, 4))
ax.barh(filled.index, filled.values, color=[BLUE if v > 0 else GREY for v in filled.values])
for y, v in enumerate(filled.values):
    ax.text(v + 1 if v > 0 else 1, y, f"{v:.0f}%" if v > 0 else "empty", va="center", fontsize=9, color=INK)
ax.set_xlim(0, 112); ax.set_xlabel("Isoforms with a value (%)")
ax.set_title("Evidence available for this transcriptome", loc="left")
plt.show()

print("Structural categories:", classif["structural_category"].value_counts().to_dict())
''')
md("""
**Everything that needs extra data is empty**, and since the reference was a placebo, SQANTI3 calls almost every isoform *intergenic*: the structural categories carry no information either. (The 25% of isoforms without a value in `all_canonical`, `psauron_score` or `predicted_NMD` are mono-exonic or have no ORF, so those columns don't apply to them.)

What we have left are the features SQANTI3 computes from the genome sequence and the reads themselves:

| Column | What it measures | Artifact it catches |
|---|---|---|
| `perc_A_downstream_TTS` | % of A in the 20 genomic bases right after the 3′ end | **Intrapriming**: the oligo-dT primer bound an A-rich stretch *inside* the gene instead of the real polyA tail |
| `RTS_stage` | TRUE if a junction sits between short direct repeats | **RT switching**: the reverse transcriptase jumped between repeats and created a fake intron |
| `all_canonical` | whether all splice junctions use GT-AG, GC-AG or AT-AC | Alignment errors that create **non-canonical junctions** |
| `coding`, `CDS_length`, `psauron_score` | whether there is an ORF, its length, and how protein-like it looks (PSAURON, 0–1) | **Fragments** and noise with no real coding sequence |
| `predicted_NMD` | TRUE if the stop codon falls well before the last junction, so the mRNA would be degraded by nonsense-mediated decay | **Mis-spliced** or broken models |
| `FL` | number of full-length reads behind the isoform | Not an artifact flag but **support**: our only "extra" evidence |
""")

# ---------------------------------------------------------------- step 2
md("""
## Step 2 · How the rules filter thinks
The rules live in a JSON file. This example says *keep an isoform if it has no intrapriming AND canonical junctions, OR if it has at least 10 reads*:

```json
{
  "rest": [
    { "perc_A_downstream_TTS": [0, 59], "all_canonical": "canonical" },
    { "FL": 10 }
  ]
}
```

- **Top-level keys are structural categories** (`full-splice_match`, `novel_not_in_catalog`…). `"rest"` covers every category without its own rules. Today everything is intergenic, so we only use `"rest"`.<br>⚠️ Always include a `"rest"` key: in SQANTI3 v6.1, isoforms in categories without rules are discarded.
- **Inside a block, every condition must hold** (AND).
- **Between blocks, passing any one is enough** (OR). This is what lets us *rescue* isoforms.
- **Values:** `[min, max]` is a range, a single number is a minimum, and text must match exactly.
- **A missing value (NA) fails the condition**, except that junction conditions (`all_canonical`, `RTS_stage`, `predicted_NMD`, `min_cov`…) are skipped for mono-exonic isoforms.
""")

# ---------------------------------------------------------------- step 3
md("""
## Step 3 · One strict rule
Our first attempt: one block that asks for everything at once.
""")
code('''
strict = {
    "rest": [
        {
            "exons": 1,                         # any number of exons (minimum 1)
            "perc_A_downstream_TTS": [0, 59],   # no intrapriming
            "all_canonical": "canonical",       # canonical splice junctions
            "coding": "coding",                 # has an ORF...
            "CDS_length": 100,                  # ...of at least 100 nt...
            "psauron_score": 0.7,               # ...that looks like a real protein
            "predicted_NMD": "FALSE"            # not an NMD target
        }
    ]
}
score({"rest": [{"exons": 1}]}, "No filter")
score(strict, "1 strict rule")
plot_tradeoff()
''')
md("""
Precision jumps from 46% to 67%, which is great. But **we threw away about 1,200 real WormBase transcripts** (13%). In an annotation project, those are genes or isoforms missing from your final gene models. Which condition is to blame?
""")

# ---------------------------------------------------------------- step 4
md("## Step 4 · Who did we lose, and why?")
code('''
explain_losses(strict["rest"][0])
''')
md("""
**The intrapriming condition alone removes ~2,160 isoforms, and 900 of them are real transcripts.** The worm genome is AT-rich, so many genuine 3′ ends happen to sit next to an A-rich stretch. Can we tell the real ones from the artifacts without WormBase? Let's look at how many reads support them.
""")
code('''
flagged = classif[~classif["perc_A_downstream_TTS"].between(0, 59)]
bins = pd.cut(flagged["FL"], [0, 2, 4, 9, 19, 49, float("inf")], labels=["2", "3–4", "5–9", "10–19", "20–49", "≥ 50"])
pct = 100 * flagged.groupby(bins, observed=True)["matches_WormBase"].mean()
n = flagged.groupby(bins, observed=True).size()

fig, ax = plt.subplots(figsize=(7, 3.4))
ax.bar(pct.index.astype(str), pct.values, color=BLUE, width=0.6)
for x, (p, k) in enumerate(zip(pct.values, n.values)):
    ax.text(x, p + 1.5, f"{p:.0f}%", ha="center", fontsize=10, color=INK)
    ax.text(x, -9, f"n = {k:,}", ha="center", fontsize=8, color=INK)
ax.set_ylim(-12, 80); ax.set_yticks([0, 20, 40, 60, 80])
ax.set_xlabel("Full-length reads supporting the isoform (FL)", labelpad=14)
ax.set_ylabel("Match WormBase (%)")
ax.set_title("Isoforms flagged for intrapriming: how many are real?", loc="left")
plt.show()
''')
md("""
With few reads, a flagged isoform is almost always an artifact. **With 50 reads or more, two out of three are real.** Read support is the evidence that lets us forgive the intrapriming flag, and the rules filter can express exactly that with a second block.
""")

# ---------------------------------------------------------------- step 5
md("## Step 5 · A second rule to rescue them")
code('''
rescue = {
    "rest": [
        strict["rest"][0],                  # block 1: the strict rule, unchanged
        {                                   # block 2: everything except the intrapriming check...
            "exons": 1,
            "all_canonical": "canonical",
            "coding": "coding",
            "CDS_length": 100,
            "psauron_score": 0.7,
            "predicted_NMD": "FALSE",
            "FL": 10                        # ...as long as at least 10 reads support it
        }
    ]
}
score(rescue, "Strict + rescue (FL ≥ 10)")
plot_tradeoff()
''')
md("""
**About 680 real transcripts come back, and precision doesn't move.** One extra block, and the filter got more complete without getting dirtier.
""")

# ---------------------------------------------------------------- step 6
md("""
## Step 6 · Fine-tune the rescue
Is 10 reads the right threshold? Move the slider. The faint line shows every threshold from 2 to 200 reads.
""")
code('''
import ipywidgets as widgets

curve = []
for fl in [2, 3, 4, 5, 7, 10, 15, 20, 30, 50, 75, 100, 150, 200]:
    block2 = dict(rescue["rest"][1], FL=fl)
    keep = sqanti_rules(classif, {"rest": [strict["rest"][0], block2]})
    good = (keep & real).sum()
    curve.append(dict(FL=fl, recovered=100 * good / N_REAL, precision=100 * good / keep.sum()))
curve = pd.DataFrame(curve)

@widgets.interact(min_reads=widgets.IntSlider(value=10, min=2, max=200, step=1, continuous_update=False,
                                              description="FL ≥"))
def tune(min_reads):
    block2 = dict(rescue["rest"][1], FL=min_reads)
    score({"rest": [strict["rest"][0], block2]}, "Your FL threshold")
    plot_tradeoff(curve=curve, highlight="Your FL threshold")
''')
md("""
There is no single right answer. If the transcripts become **hints for a gene predictor**, you may want completeness; if they are your **final high-confidence set**, you may want precision. The rules make that choice explicit and reproducible.
""")

# ---------------------------------------------------------------- step 7
md("""
## Step 7 · One rescue per risky condition
The same idea, applied to every condition that cost us real transcripts in Step 4: each extra block drops **one** condition and asks for ≥ 10 reads instead. This is the rule set we use in SQANTI-evidence to annotate genomes from long reads alone.
""")
code('''
designed = {
    "rest": [
        strict["rest"][0],                                                       # strict
        {k: v for k, v in strict["rest"][0].items() if k != "perc_A_downstream_TTS"} | {"FL": 10},  # forgive intrapriming
        {k: v for k, v in strict["rest"][0].items() if k != "psauron_score"} | {"FL": 10},          # forgive a low PSAURON score
        {k: v for k, v in strict["rest"][0].items() if k != "all_canonical"} | {"FL": 10},          # forgive a non-canonical junction
    ]
}
print(json.dumps(designed, indent=2))
score(designed, "Designed (4 blocks)")
plot_tradeoff()
''')

# ---------------------------------------------------------------- step 8
md("""
## Step 8 · Your turn 🏆
Edit `my_rules` below and run the cell. Can you beat the designed set? **Post your *recovered* and *precision* numbers in the chat.**

Some ideas:
- **Drop mono-exonic isoforms:** set `"exons": 2` in every block. Only 45 of the 5,124 mono-exonic isoforms are in WormBase.
- **Add `"RTS_stage": "FALSE"`** to block 1. Does RT switching matter in this dataset?
- **Make the protein check stricter or looser** (`psauron_score`, `CDS_length`).
- **Add a rescue block for `predicted_NMD`**, or ask for more reads in a block.
""")
code('''
my_rules = {
    "rest": [
        {"exons": 1, "perc_A_downstream_TTS": [0, 59], "all_canonical": "canonical",
         "coding": "coding", "CDS_length": 100, "psauron_score": 0.7, "predicted_NMD": "FALSE"},
        {"exons": 1, "all_canonical": "canonical",
         "coding": "coding", "CDS_length": 100, "psauron_score": 0.7, "predicted_NMD": "FALSE", "FL": 10},
        {"exons": 1, "perc_A_downstream_TTS": [0, 59], "all_canonical": "canonical",
         "coding": "coding", "CDS_length": 100, "predicted_NMD": "FALSE", "FL": 10},
        {"exons": 1, "perc_A_downstream_TTS": [0, 59],
         "coding": "coding", "CDS_length": 100, "psauron_score": 0.7, "predicted_NMD": "FALSE", "FL": 10},
    ]
}
score(my_rules, "My rules")
plot_tradeoff()
display(pd.DataFrame(scoreboard).round(1))
''')

# ---------------------------------------------------------------- step 9
md("""
## Step 9 · Take it home
Save your rules as a JSON file. On your own data, the same file goes straight into SQANTI3:

```bash
sqanti3_filter.py rules \\
    --sqanti_class my_species_classification.txt \\
    -j my_rules.json \\
    --filter_gtf my_species_corrected.gtf \\
    -d filter_output -o my_species
```
""")
code('''
with open("my_rules.json", "w") as fh:
    json.dump(my_rules, fh, indent=2)
print(open("my_rules.json").read())

try:
    from google.colab import files
    files.download("my_rules.json")
except ImportError:
    pass
''')
md("""
**Optional (≈ 2–3 min): check against the real SQANTI3 code.** This cell downloads SQANTI3 and runs its own rules filter on the classification file with your rules. It should keep exactly the same isoforms as our fast version.
""")
code(f'''
import os, sys, subprocess, logging
if not os.path.exists("SQANTI3"):
    subprocess.run(["git", "clone", "-q", "--depth", "1", "https://github.com/ConesaLab/SQANTI3.git"], check=True)
sys.path.insert(0, "SQANTI3")
from src.module_logging import filter_logger, update_logger
from src.utilities.filter.sqanti3_rules_filter import rules_filter

os.makedirs("filter_output/logs", exist_ok=True)
update_logger(filter_logger, "filter_output", "filter", "INFO")
classif.drop(columns=answer_key.columns[1:]).to_csv("filter_output/celegans_classification.txt", sep="\\t", index=False)
rules_filter("filter_output/celegans_classification.txt", "my_rules.json", False,
             "filter_output/celegans", logging.getLogger("seminar"))

sqanti_pass = set(open("filter_output/celegans_pass_isoforms.txt").read().split())
ours = set(classif.loc[sqanti_rules(classif, my_rules), "isoform"])
print(f"SQANTI3 kept {{len(sqanti_pass):,}} isoforms; our version kept {{len(ours):,}}. Identical: {{sqanti_pass == ours}}")
display(pd.read_csv("filter_output/celegans_filtering_reasons.txt", sep="\\t").head())
''')

# ---------------------------------------------------------------- wrap-up
md("""
## Wrap-up
- **A single strict rule buys precision at the cost of real transcripts.** Several blocks let you say *"I accept this risk only when there's other support"*.
- **In a low-evidence setting, read support (FL) is your best rescuer.** Check what it rescues before trusting it, as we did in Step 4.
- **"Not in WormBase" doesn't mean "artifact".** WormBase is incomplete and some of those isoforms are real novel ones, so the precision you saw here is a lower bound.

**When you get more evidence, add it to the rules:**

| Evidence | Column | Example condition |
|---|---|---|
| Short reads | `min_cov` | `"min_cov": 3` (every junction covered by ≥ 3 short reads) |
| CAGE | `within_CAGE_peak` | `"within_CAGE_peak": "TRUE"` |
| PolyA sites / motifs | `within_polyA_site`, `polyA_motif_found` | `"polyA_motif_found": "TRUE"` |
| A good reference annotation | structural categories | separate blocks: loose for `full-splice_match`, strict for `novel_not_in_catalog` |

With a good reference you can also try the **SQANTI3 machine-learning filter**, which learns from known transcripts. Without one, the rules filter is the way to go.

📚 SQANTI3: <https://github.com/ConesaLab/SQANTI3> · Wiki on filtering: <https://github.com/ConesaLab/SQANTI3/wiki/Running-SQANTI3-filter>
""")

nb = {
    "cells": cells,
    "metadata": {
        "colab": {"provenance": [], "toc_visible": True},
        "kernelspec": {"display_name": "Python 3", "name": "python3"},
        "language_info": {"name": "python"},
    },
    "nbformat": 4,
    "nbformat_minor": 0,
}
out = Path(__file__).resolve().parent.parent / "rules_filter_activity.ipynb"
out.write_text(json.dumps(nb, indent=1, ensure_ascii=False))
print(f"Wrote {out} ({len(cells)} cells)")
