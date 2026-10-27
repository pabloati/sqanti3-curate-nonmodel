#!/usr/bin/env python3
"""Build rules_filter_activity_B.ipynb, the VALT-style version of the seminar practical.

Edit this script, not the notebook:  python3 build/build_notebook_B.py
"""
import json
from pathlib import Path

# Where the data and scripts are hosted
REPO_RAW = "https://raw.githubusercontent.com/pabloati/sqanti3-curate-nonmodel/main"
COLAB_URL = "https://colab.research.google.com/github/pabloati/sqanti3-curate-nonmodel/blob/main/rules_filter_activity_B.ipynb"
IMAGES = "https://raw.githubusercontent.com/ConesaLab/courses-SQANTI_verse/main/data"

cells = []


def md(text):
    cells.append({"cell_type": "markdown", "metadata": {}, "source": text.strip("\n")})


def code(text):
    cells.append({"cell_type": "code", "metadata": {}, "execution_count": None, "outputs": [],
                  "source": text.strip("\n")})


# ---------------------------------------------------------------------------------------------
md("""
# Curating a long-read transcriptome when you have (almost) nothing: the SQANTI3 rules filter

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](__COLAB_URL__)

This tutorial will guide you through the curation of a long-read transcriptome in the situation that most biodiversity projects face: **a genome, a long-read RNA-seq experiment, and nothing else**. No reference annotation, no short reads, no CAGE peaks, no polyA sites. We will use the **SQANTI3 rules filter** to remove the artifacts, and we will see how **combining several rules** recovers real transcripts that a single strict rule throws away.

> ⏳ **Before anything else, run the first code cell (Software installation).** Installing SQANTI3 takes 5-10 minutes, so let it work while you read.

<details>
<summary><strong> 📖 Theoretical background</strong></summary><br>

Long-read transcriptomes are not perfect. Next to the real isoforms, every long-read pipeline produces **artifacts**: fragments of degraded RNA, molecules primed inside the gene instead of at the polyA tail (**intrapriming**), fake introns created when the reverse transcriptase jumps between repeats (**RT switching**), or junctions broken by sequencing errors. **SQANTI3** characterizes every isoform with dozens of descriptors (the *classification* file) and lets you remove the artifacts with its **filter** module.

In a model organism, the most informative descriptors come from the reference annotation (structural categories, distance to the annotated TSS and TTS) and from orthogonal data (short-read junction support, CAGE peaks, polyA sites). In a newly sequenced species we have none of them. That is exactly the situation of the first step of **[SQANTI-evidence](https://github.com/ConesaLab/SQANTI_evidence)**, our pipeline to annotate genomes from long reads: the transcriptome is curated with only what the reads and the genome sequence tell us, and the curated isoforms are then given to AUGUSTUS as hints.

---
</details><br>

We will use two datasets:

1. **_Caenorhabditis elegans_** chromosome I (PacBio HiFi, 3,325 isoforms reconstructed with IsoQuant). The worm has one of the best annotations there is (WormBase), so we can check our decisions against it.
2. **Zebrafish** chromosome 5 (PacBio, 5,103 isoforms) as a bonus: much noisier data, to see how the same rules behave.

## Table of contents

- [0. Pre-requisites](#0-pre-requisites)
  - [Software installation](#software-installation)
  - [Data download](#data-download)
- [1. The data: a transcriptome with (almost) no information](#1-the-data-a-transcriptome-with-almost-no-information)
  - [1.1 The answer key](#11-the-answer-key)
- [2. SQANTI3 rules filter](#2-sqanti3-rules-filter)
  - [2.1 Rules](#21-rules)
  - [2.2 One strict rule](#22-one-strict-rule)
  - [2.3 Who did we lose?](#23-who-did-we-lose)
  - [2.4 A rescue rule](#24-a-rescue-rule)
- [3. Fine-tuning: the SQANTI-evidence rules](#3-fine-tuning-the-sqanti-evidence-rules)
  - [3.1 Your turn!](#31-your-turn)
- [4. BONUS: a noisier dataset](#4-bonus-a-noisier-dataset)
- [5. Conclusion](#5-conclusion)

⚠️ This tutorial uses the latest release of SQANTI3 (v6.1).
""".replace("__COLAB_URL__", COLAB_URL))

md("""
# 0. Pre-requisites

### Software installation

SQANTI3 can be directly installed from **bioconda**. On your own computer, that is one single command (see below). Google Colab, however, does not come with conda. For that reason, we will use **micromamba**, a tiny and fast version of conda that understands the same commands, and we will call it `conda`. This way, every command in this tutorial is exactly the one you would type on your own computer or cluster.

SQANTI3 will live in its own environment, called `sqanti3`, and we will call its scripts with `conda run -n sqanti3 <command>`, which runs a command inside that environment without activating it.

<details>
<summary><strong> 💻 Installing SQANTI3 on your own computer</strong></summary><br>

```bash
conda create -n sqanti3 -c conda-forge -c bioconda sqanti3
conda activate sqanti3
sqanti3_qc.py --version
```

If you prefer to install it from GitHub, or you run into problems, check the [installation page of the SQANTI3 wiki](https://github.com/ConesaLab/SQANTI3/wiki/Dependencies-and-installation).

---
</details><br>
""")

code("""
%%bash
# ⏳ Installs SQANTI3 from bioconda (5-10 min). Run it now and keep reading!
cd /content
curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xj bin/micromamba

# A small "conda" command that calls micromamba, so the commands look like on your computer
cat > /usr/local/bin/conda << 'EOF'
#!/bin/bash
export MAMBA_ROOT_PREFIX=/content/micromamba
exec /content/bin/micromamba "$@"
EOF
chmod +x /usr/local/bin/conda

conda create -y -q -n sqanti3 -c conda-forge -c bioconda sqanti3=6.1 > install.log 2>&1 \\
    && echo "SQANTI3 installed 🥳" || { echo "Something went wrong, check install.log"; tail install.log; }
conda run -n sqanti3 sqanti3_filter.py --version
""")

md("""
### Data download

Before starting the tutorial, it is key to have a clean and organized working environment. In bioinformatics, an organized workspace is vital, so when you come back after some time to your project, you can find and understand what you were doing, rather than spend hours searching through weirdly named directories. It is important to always create three directories:

- *scripts*: all the scripts will be stored here, with meaningful names
- *data*: raw data will go in here and, if you want and need, databases
- *results*: create a sub directory for every different process you do. If you run a process multiple times with different parameters, include them in the directory name, so you will differentiate them in the future.

Today we add a fourth one, *rules*, where our filter rules will live. The data and a small helper script are hosted in the [tutorial repository](https://github.com/pabloati/sqanti3-curate-nonmodel):
""")

code("""
%%bash
mkdir -p scripts data results rules
REPO=__REPO_RAW__

wget -q -P scripts $REPO/scripts/score.py
for f in celegans_chrI_classification.txt celegans_chrI_corrected.gtf celegans_chrI_answer_key.tsv \\
         zebrafish_chr5_classification.txt zebrafish_chr5_corrected.gtf zebrafish_chr5_answer_key.tsv; do
    wget -q -O - $REPO/data/$f.gz | gunzip > data/$f
done
ls -lh data scripts
""".replace("__REPO_RAW__", REPO_RAW))

# ---------------------------------------------------------------------------------------------
md("""
# 1. The data: a transcriptome with (almost) no information

The worm reads were reconstructed into a transcriptome with [IsoQuant](https://github.com/ablab/IsoQuant), without using any annotation (we want to annotate the genome, so it would be a bit of cheating to use the annotation to build the transcriptome 😉). Then, the transcriptome went through **SQANTI3 QC**. However, SQANTI3 needs a reference annotation to run. So... how can we use an annotation to generate an annotation? The solution is simple: we create a **"placebo" annotation**, with one fake transcript in each chromosome, so SQANTI3 can run. The only thing to keep in mind is that many of the SQANTI3 classification columns will be useless, as no true reference annotation is given.

<details>
<summary><strong> 🛠️ How was this classification made?</strong></summary><br>

```bash
## DO NOT RUN, this data has been generated for you :)
conda run -n sqanti3 \\
sqanti3_qc.py \\
    --isoforms isoquant/c_elegans.transcript_models.gtf \\
    --refGTF reference/sqanti_placebo.gtf \\
    --refFasta reference/c_elegans.WBcel235.dna.toplevel.fa \\
    --fl_count isoquant/c_elegans.discovered_transcript_counts.tsv \\
    --include_ORF --dir sqanti --output c_elegans
```

`--include_ORF` predicts the coding sequences with [TD2](https://github.com/Markusjsommer/TD2), and `--fl_count` adds the number of reads supporting each isoform (the `FL` column). Only chromosome I is used today, to keep everything fast.

---
</details><br>

Let's take a look at the classification file:
""")

code("""
import pandas as pd
pd.set_option("display.max_columns", 60)

classif = pd.read_csv("data/celegans_chrI_classification.txt", sep="\\t")
print(classif.shape)
classif["structural_category"].value_counts()
""")

md("""
**❓Trivia:** In which structural category are all the isoforms? Why?

<details><summary>Answer</summary>

All 3,325 isoforms are **intergenic**. The placebo annotation has no real genes, so SQANTI3 cannot match any isoform to a known transcript. The structural categories (and everything that depends on the reference, like `diff_to_TSS` or `associated_gene`) carry no information here.
</details><br>

So, what is left? The descriptors that come from **the reads themselves and from the genome sequence**. These are the ones we can use for filtering:

<details>
<summary><strong> 📊 SQANTI3 concepts: the columns we can still use</strong></summary><br>

| Column | What it measures | Artifact it catches |
|---|---|---|
| `FL` | Number of long reads supporting the isoform | Low-support noise |
| `perc_A_downstream_TTS` | % of A in the 20 genomic bases right after the 3′ end | **Intrapriming**: the oligo-dT primer bound an A-rich stretch *inside* the gene instead of the real polyA tail |
| `RTS_stage` | TRUE if a junction sits between short direct repeats | **RT switching**: the reverse transcriptase jumped between repeats and created a fake intron |
| `all_canonical` | Whether all splice junctions use GT-AG, GC-AG or AT-AC | Alignment or sequencing errors at junctions |
| `coding`, `CDS_length`, `psauron_score` | Whether there is an ORF, its length, and how protein-like it looks ([PSAURON](https://github.com/salzberg-lab/PSAURON), 0 to 1) | Fragments and noise with no real coding sequence |
| `predicted_NMD` | TRUE if the stop codon falls more than 50 bp before the last junction | **Mis-spliced** or broken models (the mRNA would be degraded by nonsense-mediated decay) |

* **Canonical junctions**: GT-AG, GC-AG and AT-AC are found in over 99.9% of introns, so a non-canonical junction is often a sign of an error.
* **RT switching** events are enriched in low-abundance transcripts from highly expressed genes.
* **Intrapriming** produces isoforms that end too early, in an A-rich region of the genome. SQANTI3 flags them through the A content after the 3′ end: above 60%, the isoform is suspicious.

---
</details><br>
""")

code("""
cols = ["isoform", "exons", "length", "FL", "perc_A_downstream_TTS", "RTS_stage", "all_canonical",
        "coding", "CDS_length", "psauron_score", "predicted_NMD"]
classif[cols].head(10)
""")

md("""
## 1.1 The answer key

Here comes a little bit of cheating. *C. elegans* has a very good annotation, so we compared our isoforms against it with [gffcompare](https://github.com/gpertea/gffcompare). Isoforms whose intron chain matches a WormBase transcript exactly (class `=`) are considered **real**. You will **never** have this for your own species, but today it lets us *see* whether our filtering decisions are good. For each filter we will compute two numbers:

- **Recovered**: the percentage of the real isoforms that the filter keeps (we want it high).
- **Precision**: the percentage of the kept isoforms that are real (we want it high too).
""")

code("""
key = pd.read_csv("data/celegans_chrI_answer_key.tsv", sep="\\t")
key["matches_reference"].value_counts()
""")

md("""
Only 1,469 of the 3,325 isoforms (44%) match WormBase exactly. That is our starting point: **without any filter, more than half of the transcriptome is not a known transcript**.
""")

# ---------------------------------------------------------------------------------------------
md("""
# 2. SQANTI3 rules filter

The **SQANTI3 filter** module provides two approaches for curating long-read transcriptomes: a **rules-based filter** and a **machine learning filter**. The machine learning filter trains a random forest on sets of true and false isoforms, which are usually defined from the reference annotation... which we don't have. So, today we will focus on the rules filter, which allows for more flexibility and customization from the user.

## 2.1 Rules

The rules filter is based on a set of rules (surprise!) that are defined in a JSON file. The way the rules work is the following:

1. Create a set of rules for a specific structural category, or for all of them under the tag `"rest"`.
2. Within a rule, the user can set multiple conditions, which are based on the columns of the SQANTI3 classification file. These conditions are treated as **AND** conditions, meaning that all conditions within a rule have to be met for the isoform to be kept.
3. A category can have several rules, which are treated as **OR** conditions, meaning that if an isoform passes any of the rules, it will be kept.

The conditions must have the same name as the columns in the classification file, and you can use 3 possible value types:

1. A single number to denote the minimum value (`"FL": 10` means at least 10 reads).
2. A range, using the format `[min,max]` (`"perc_A_downstream_TTS": [0,59]`).
3. A text or a list of texts to denote a specific value or a list of values (`"coding": "coding"`).

>Note: In case of a condition being added but the column being an NA (not having information in the classification file), the rule check will automatically fail. The only exception are the junction-related columns (`all_canonical`, `RTS_stage`, `predicted_NMD`...), which are skipped for mono-exonic isoforms, since they have no junctions.

<details>
<summary><strong> ⚠️ Two traps when writing rules</strong></summary><br>

1. **Always include a `"rest"` key.** Categories with no rules of their own are filtered with `"rest"`. If there is no `"rest"` key, all of those isoforms are discarded (even if the log says otherwise).
2. **Write TRUE and FALSE in quotes**: `"RTS_stage": "FALSE"`. Without the quotes, JSON reads `false` as a boolean, and SQANTI3 v6.1 turns it into the number 0, so the condition keeps everything. We found both while preparing this tutorial, and they will be fixed in the next release 🐛.

---
</details><br>

Since all our isoforms are intergenic, **all of today's rules go under `"rest"`**.

## 2.2 One strict rule

Let's start with the obvious approach: one rule asking for everything we would like a real isoform to have. Keep an isoform if:

* It is **not** a potential intrapriming product (`perc_A_downstream_TTS` between 0 and 59) **AND**
* All its junctions are canonical **AND**
* It is coding, with a CDS of at least 100 nt and a PSAURON score of at least 0.7 **AND**
* It is not predicted to be degraded by NMD.

The `%%writefile` command at the top of the cell saves the content of the cell into a file:
""")

code("""
%%writefile rules/01_strict.json
{
    "rest": [
        {
            "perc_A_downstream_TTS": [0, 59],
            "all_canonical": "canonical",
            "coding": "coding",
            "CDS_length": 100,
            "psauron_score": 0.7,
            "predicted_NMD": "FALSE"
        }
    ]
}
""")

md("""
Now, let's dive into running it. The filter module has two subcommands, `rules` and `ml`. The inputs are the classification file, the JSON with the rules and, optionally, the files that we want to be filtered (here, the GTF of the transcriptome):
""")

code("""
%%bash
conda run -n sqanti3 \\
sqanti3_filter.py rules \\
    --sqanti_class data/celegans_chrI_classification.txt \\
    --filter_gtf data/celegans_chrI_corrected.gtf \\
    --json_filter rules/01_strict.json \\
    --skip_report \\
    --dir results/01_strict --output worm
ls results/01_strict
""")

md("""
<details>
<summary><strong> 📤 Output files</strong></summary><br>

- `worm_RulesFilter_classification.txt`: the classification file with a new column, `filter_result`, which says whether each isoform is an `Isoform` (kept) or an `Artifact`.
- `worm_filtering_reasons.txt`: for every artifact, the conditions it failed. We will love this file in a minute.
- `worm_pass_isoforms.txt`: the list of isoforms that passed the filter.
- `worm.filtered.gtf`: the GTF with only the isoforms that passed. This is your curated transcriptome!
- `worm_params.txt`: the parameters of the run, to keep track of what you did.

`--skip_report` skips the HTML report, which needs R and takes a while. On your own computer, leave it out. SQANTI3 also tells you a lot about what it is doing (the log). From now on, we will add `--log_level WARNING`, so it only talks when something goes wrong and the results are easier to see.

---
</details><br>

Let's see how well we did, using the `scripts/score.py` helper and the answer key:
""")

code("""
%%bash
python scripts/score.py results/01_strict --key data/celegans_chrI_answer_key.tsv
""")

md("""
**❓Trivia:** How much did the precision improve? And how many real transcripts did we lose on the way?

<details><summary>Answer</summary>

The precision went from **44%** (no filter) to **67.2%**: two out of three isoforms we keep are now real. But we kept only 1,274 of the 1,469 real isoforms, so **we lost 195 real transcripts (13%)**. A single strict rule cleans a lot, but it also throws away the baby with the bath water.
</details><br>

## 2.3 Who did we lose?

To understand what went wrong, SQANTI3 tells us, for every discarded isoform, which conditions it failed:
""")

code("""
%%bash
head -5 results/01_strict/worm_filtering_reasons.txt
""")

md("""
An isoform can fail several conditions at once (they are separated by `;`). Let's count how many isoforms fail each condition, using some classic bash:
""")

code("""
%%bash
cut -f3 results/01_strict/worm_filtering_reasons.txt | tail -n +2 | \\
    tr ';' '\\n' | sed -E 's/^ +//; s/^NA value in //; s/:.*//' | \\
    sort | uniq -c | sort -nr
""")

md("""
This counts all the isoforms. But we care about the **real** ones that we threw away. Let's cross the reasons with the answer key (`--losses`). The last column is the important one: real isoforms that failed **only** that condition, so they would come back if we were more flexible with it.
""")

code("""
%%bash
python scripts/score.py results/01_strict --key data/celegans_chrI_answer_key.tsv --losses
""")

md("""
**❓Trivia:** Which condition throws away the most real transcripts? Can you think of a biological reason?

<details><summary>Answer</summary>

**Intrapriming** (`perc_A_downstream_TTS`): 1,051 isoforms fail it, 162 of them are real, and **157 real isoforms fail only that condition**. Many real genes simply have an A-rich stretch right after their 3′ end, so the intrapriming flag is a probability, not a certainty. The good news: a real isoform from a well-expressed gene is supported by many reads, while an intrapriming artifact usually is not.
</details><br>

## 2.4 A rescue rule

Remember that the rules in a category are treated as **OR**? This is where they shine. We keep the strict rule, and we add a second one: **forgive the intrapriming flag, as long as at least 10 reads support the isoform**. Everything else stays the same.
""")

code("""
%%writefile rules/02_rescue_intrapriming.json
{
    "rest": [
        {
            "perc_A_downstream_TTS": [0, 59],
            "all_canonical": "canonical",
            "coding": "coding",
            "CDS_length": 100,
            "psauron_score": 0.7,
            "predicted_NMD": "FALSE"
        },
        {
            "all_canonical": "canonical",
            "coding": "coding",
            "CDS_length": 100,
            "psauron_score": 0.7,
            "predicted_NMD": "FALSE",
            "FL": 10
        }
    ]
}
""")

code("""
%%bash
conda run -n sqanti3 \\
sqanti3_filter.py rules \\
    --sqanti_class data/celegans_chrI_classification.txt \\
    --filter_gtf data/celegans_chrI_corrected.gtf \\
    --json_filter rules/02_rescue_intrapriming.json \\
    --skip_report --log_level WARNING \\
    --dir results/02_rescue_intrapriming --output worm

python scripts/score.py results/01_strict results/02_rescue_intrapriming --key data/celegans_chrI_answer_key.tsv
""")

md("""
**❓Trivia:** How many real transcripts came back? What did it cost?

<details><summary>Answer</summary>

The second rule brought back **191 isoforms, and 118 of them are real** (62%, far above the 44% of the raw transcriptome). We now recover **94.8%** of the real transcripts instead of 86.7%, and the precision barely moves (67.2% → 66.7%). One extra rule, almost for free!
</details><br>

🥳 You just wrote your first rescue rule. Congrats!
""")

# ---------------------------------------------------------------------------------------------
md("""
# 3. Fine-tuning: the SQANTI-evidence rules

We can apply the same trick to the other conditions that removed real transcripts. In SQANTI-evidence, the rules are built this way:

1. **One strict rule**, with everything we would like a real isoform to have.
2. **One rescue rule per "forgivable" condition**: the strict rule without that condition, but asking for at least 10 reads. Each rescue rule forgives **one** thing (intrapriming, a low PSAURON score, or a non-canonical junction), never two at the same time.
""")

code("""
%%writefile rules/03_sqanti_evidence.json
{
    "rest": [
        {
            "perc_A_downstream_TTS": [0, 59], "all_canonical": "canonical",
            "coding": "coding", "CDS_length": 100, "psauron_score": 0.7, "predicted_NMD": "FALSE"
        },
        {
            "all_canonical": "canonical",
            "coding": "coding", "CDS_length": 100, "psauron_score": 0.7, "predicted_NMD": "FALSE", "FL": 10
        },
        {
            "perc_A_downstream_TTS": [0, 59], "all_canonical": "canonical",
            "coding": "coding", "CDS_length": 100, "predicted_NMD": "FALSE", "FL": 10
        },
        {
            "perc_A_downstream_TTS": [0, 59],
            "coding": "coding", "CDS_length": 100, "psauron_score": 0.7, "predicted_NMD": "FALSE", "FL": 10
        }
    ]
}
""")

code("""
%%bash
conda run -n sqanti3 \\
sqanti3_filter.py rules \\
    --sqanti_class data/celegans_chrI_classification.txt \\
    --filter_gtf data/celegans_chrI_corrected.gtf \\
    --json_filter rules/03_sqanti_evidence.json \\
    --skip_report --log_level WARNING \\
    --dir results/03_sqanti_evidence --output worm

python scripts/score.py results/01_strict results/02_rescue_intrapriming results/03_sqanti_evidence \\
    --key data/celegans_chrI_answer_key.tsv
""")

md("""
The SQANTI-evidence rules recover **95.5%** of the real transcripts with a precision of **66.4%**. On the whole worm genome, these rules keep 13,376 isoforms, which become the curated layer of the SQANTI-evidence annotation.

Let's put all the runs in one picture. The best filters are in the top-right corner: they keep almost all the real transcripts, and little else.
""")

code("""
import glob
import sys
import matplotlib.pyplot as plt
sys.path.append("scripts")
import score

runs = sorted(glob.glob("results/0*"))
table = score.score([r for r in runs if "zebrafish" not in r], "data/celegans_chrI_answer_key.tsv")
display(table)

fig, ax = plt.subplots(figsize=(8, 4.5))
fig.subplots_adjust(right=0.62)
ax.scatter(table["Recovered (%)"], table["Precision (%)"], s=80, color="#1f77b4")
seen = {}
for _, r in table.iterrows():
    xy = (r["Recovered (%)"], r["Precision (%)"])
    n = seen.get(xy, 0)          # runs with the same result get their labels stacked
    seen[xy] = n + 1
    ax.annotate(r["Run"][:2], xy, textcoords="offset points", xytext=(8, 4 - 13 * n), fontsize=10)
ax.set_xlim(table["Recovered (%)"].min() - 2, 100)
ax.set_ylim(table["Precision (%)"].min() - 3, table["Precision (%)"].max() + 3)
ax.set_xlabel("Real transcripts recovered (%)")
ax.set_ylabel("Kept isoforms that are real (%)")
ax.set_title("Top right is best", fontsize=11)
key_text = "\\n".join(f"{r['Run'][:2]}  {r['Run'][3:]}" for _, r in table.iterrows())
ax.text(1.04, 1, key_text + "\\n\\nNo filter: 100% recovered,\\n44.2% real", transform=ax.transAxes,
        va="top", fontsize=9, family="monospace")
plt.show()
""")

md("""
## 3.1 Your turn!

Now it is your time to show all you have learned. Edit the rules below, run the two cells and **post your *Recovered* and *Precision* numbers in the chat**. Can you beat the SQANTI-evidence rules? Some ideas:

- **Drop the mono-exonic isoforms** by adding `"exons": 2` to every rule. How many real transcripts have a single exon?
- **Add `"RTS_stage": "FALSE"`** to every rule. Does RT switching matter in this dataset?
- **Make the protein check stricter or looser** (`psauron_score`, `CDS_length`).
- **Ask for more (or fewer) reads** in the rescue rules.

<details>
<summary>Click here if desperate (or lazy)</summary>

- Adding `"exons": 2` to every rule gives **95.3% recovered and 68.2% precision**: only 3 of the 906 mono-exonic isoforms are real, so we lose almost nothing and the transcriptome gets cleaner.
- Adding `"RTS_stage": "FALSE"` to every rule drops the recovery to **84.6%** with no gain in precision. Of the 289 multi-exon isoforms flagged for RT switching, 165 are real: in this dataset the flag removes more real transcripts than artifacts.
</details><br>
""")

code("""
%%writefile rules/04_my_rules.json
{
    "rest": [
        {
            "perc_A_downstream_TTS": [0, 59], "all_canonical": "canonical",
            "coding": "coding", "CDS_length": 100, "psauron_score": 0.7, "predicted_NMD": "FALSE"
        },
        {
            "all_canonical": "canonical",
            "coding": "coding", "CDS_length": 100, "psauron_score": 0.7, "predicted_NMD": "FALSE", "FL": 10
        },
        {
            "perc_A_downstream_TTS": [0, 59], "all_canonical": "canonical",
            "coding": "coding", "CDS_length": 100, "predicted_NMD": "FALSE", "FL": 10
        },
        {
            "perc_A_downstream_TTS": [0, 59],
            "coding": "coding", "CDS_length": 100, "psauron_score": 0.7, "predicted_NMD": "FALSE", "FL": 10
        }
    ]
}
""")

code("""
%%bash
conda run -n sqanti3 \\
sqanti3_filter.py rules \\
    --sqanti_class data/celegans_chrI_classification.txt \\
    --filter_gtf data/celegans_chrI_corrected.gtf \\
    --json_filter rules/04_my_rules.json \\
    --skip_report --log_level WARNING \\
    --dir results/04_my_rules --output worm

python scripts/score.py results/03_sqanti_evidence results/04_my_rules --key data/celegans_chrI_answer_key.tsv
""")

# ---------------------------------------------------------------------------------------------
md("""
# 4. BONUS: a noisier dataset

The worm data is pretty clean. What happens with the same rules on a noisier dataset? Let's try **zebrafish chromosome 5**, processed exactly in the same way (IsoQuant, SQANTI3 against a placebo annotation). Here, the answer key comes from the Ensembl annotation.
""")

code("""
%%bash
conda run -n sqanti3 \\
sqanti3_filter.py rules \\
    --sqanti_class data/zebrafish_chr5_classification.txt \\
    --filter_gtf data/zebrafish_chr5_corrected.gtf \\
    --json_filter rules/03_sqanti_evidence.json \\
    --skip_report --log_level WARNING \\
    --dir results/05_zebrafish --output zebrafish

python scripts/score.py results/05_zebrafish --key data/zebrafish_chr5_answer_key.tsv
""")

md("""
Without any filter, only **9.2%** of the zebrafish isoforms match Ensembl. The SQANTI-evidence rules take the precision to **35.9%** (four times higher!) and keep **88.5%** of the real transcripts. Much better, but far from the worm.

**❓Trivia:** Why is the precision so much lower than in the worm? Use the cell below to compare both datasets (*hint: look at the number of exons and at the reads*).
""")

code("""
zf = pd.read_csv("data/zebrafish_chr5_classification.txt", sep="\\t")
zf_key = pd.read_csv("data/zebrafish_chr5_answer_key.tsv", sep="\\t")
worm = classif.merge(key, on="isoform")
zf = zf.merge(zf_key, on="isoform")

for name, df in [("worm", worm), ("zebrafish", zf)]:
    mono = df["exons"] == 1
    print(f"{name:10s} isoforms: {len(df):5d} | mono-exonic: {mono.mean():.0%} "
          f"({df.loc[mono, 'matches_reference'].sum()} real) | median reads per isoform: {df['FL'].median():.0f}")
""")

md("""
<details><summary>Answer</summary>

- **60% of the zebrafish isoforms are mono-exonic** (3,060), and only 3 of them match Ensembl. In the worm they are 27%.
- **The zebrafish isoforms have fewer reads**: a median of 4 reads per isoform, against 13 in the worm. Less depth means less evidence for each isoform, and the `FL` rescue rules can rescue less.
- **Ensembl is not WormBase**: a real zebrafish transcript missing from Ensembl counts as "not real" here, so the true precision is higher than 35.9%.

Adding `"exons": 2` to every rule raises the precision to 38.6% while keeping 88.1% of the real transcripts. Try it! The lesson: **the rules have to be tuned to your data**, and the noisier the data, the more the filter matters.
</details><br>
""")

md("""
# 5. Conclusion

With this, you have reached the end of this tutorial 🥳. You have learned how to curate a long-read transcriptome when the only information you have comes from the reads and the genome: how the SQANTI3 rules are written (AND inside a rule, OR between rules), how to find out *why* isoforms are discarded with the filtering reasons, and how **rescue rules recover real transcripts** that a single strict rule throws away, at almost no cost in precision.

In your own species, you will not have an answer key. To judge your filters, look at how many genes (loci) you keep, at the fraction of coding isoforms, and at the completeness of the resulting proteome with [BUSCO](https://busco.ezlab.org/) or [OMArk](https://github.com/DessimozLab/OMArk). And remember: the filtered GTF that you just produced (`results/03_sqanti_evidence/worm.filtered.gtf`) is exactly what SQANTI-evidence gives to AUGUSTUS as hints to annotate the genome.

Want to know more? Check the [SQANTI3 wiki](https://github.com/ConesaLab/SQANTI3/wiki) and the [SQANTI-evidence repository](https://github.com/ConesaLab/SQANTI_evidence).

![Picture of what happens if you use the SQANTI-colors right](__IMAGES__/ConesaColors_happy.jpg "Conesa Color")
""".replace("__IMAGES__", IMAGES))

notebook = {
    "cells": cells,
    "metadata": {
        "colab": {"provenance": [], "toc_visible": True},
        "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
        "language_info": {"name": "python"},
    },
    "nbformat": 4,
    "nbformat_minor": 0,
}
out = Path(__file__).resolve().parent.parent / "rules_filter_activity_B.ipynb"
out.write_text(json.dumps(notebook, indent=1, ensure_ascii=False) + "\n")
print(f"Wrote {out} ({len(cells)} cells)")
