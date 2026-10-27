# sqanti3-curate-nonmodel

Practical session materials on how to curate a long-read transcriptome when you have no good reference annotation, using the **SQANTI3 rules filter**.

Both notebooks run in Google Colab: click a badge, sign in with a Google account and run the cells in order. To keep your changes, use **File → Save a copy in Drive**.

| Notebook | What it covers | Open |
|---|---|---|
| **Quick hands-on** (`rules_filter_activity.ipynb`) | Whole *C. elegans* transcriptome. Rules are scored directly in Python, no installation needed (≈ 20 s set-up). | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/pabloati/sqanti3-curate-nonmodel/blob/main/rules_filter_activity.ipynb) |
| **Full tutorial** (`rules_filter_activity_B.ipynb`) | Installs SQANTI3 and runs `sqanti3_filter.py rules` on *C. elegans* chromosome I, with a bonus on zebrafish chromosome 5. Installation takes 5-10 min. | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/pabloati/sqanti3-curate-nonmodel/blob/main/rules_filter_activity_B.ipynb) |

## Repository layout

```
├── rules_filter_activity.ipynb      quick hands-on
├── rules_filter_activity_B.ipynb    full tutorial
├── build/                           scripts that generate the notebooks
├── scripts/score.py                 scores a filter run against the answer key
└── data/                            SQANTI3 classifications, GTFs and answer keys (gzipped)
```

## Editing the notebooks

The notebooks are generated: edit the scripts in `build/`, not the `.ipynb` files, and rebuild from the repository root:

```bash
python3 build/build_notebook_A.py
python3 build/build_notebook_B.py
```

The notebooks download their data from this repository's `main` branch, so changes to `data/` or `scripts/` only reach Colab once they are pushed.
