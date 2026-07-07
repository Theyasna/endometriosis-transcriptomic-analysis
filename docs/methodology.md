# Methodological Notes

## Why Per-Dataset Analysis (Not Cross-Study Pooling or ComBat)

The two datasets use different Affymetrix platforms (GPL571 for GSE25628, GPL570 for GSE7305) and were generated in different labs. Their biological groups also do not overlap except for "Normal":

```
              Diseased  Ectopic  Eutopic  Normal
GSE25628             0        8        8       6
GSE7305             10        0        0      10
```

Pooling samples across the two datasets would confound platform and lab with biology. An earlier version of this pipeline attempted cross-study batch correction (ComBat), but with this group structure ComBat cannot separate batch from biology for any group present in only one dataset — it would either fail to remove the batch effect or remove real biological signal, with no way to tell which.

**The approach used here:** run each dataset's differential expression independently, then compare the two *lesion-vs-control* results at the gene-list level. Both datasets happen to contain a lesion-vs-control comparison:
- GSE25628: ectopic lesion vs. healthy-donor endometrium
- GSE7305: ovarian lesion vs. matched patient control endometrium

Overlapping these two independently-derived DEG lists (608 shared genes, 37.9%) provides **cross-cohort replication** across two platforms and two control definitions — a stronger validation than any single-dataset result, and one that is robust precisely because the two controls differ.

## The GSE7305 Paired Design (and Why the Current Analysis Treats It as Unpaired)

In GSE7305, lesion and control samples come from the **same patients**, collected at the same time. A statistically ideal analysis would use a paired/blocked limma model (`~ patient + group`) to account for within-patient correlation and gain power.

The current pipeline uses an unpaired two-group model (`~ 0 + group`). This is a deliberate, conservative choice: it does not exploit the pairing for extra power, but it keeps the GSE7305 analysis structurally parallel to GSE25628 (which is not paired). The unpaired result is therefore a lower bound on statistical power, not an inflated one. A paired reanalysis (`~ patient + group`) is a reasonable future refinement.

## Why C7 ImmuneSigDB (Not C8)

An earlier version used MSigDB C8 (a pan-tissue single-cell reference atlas). Checking the returned gene sets, roughly 82% corresponded to non-immune cell types (fetal organ, neural, smooth muscle, etc.), producing significant enrichment scores with no immunological interpretation. C7 ImmuneSigDB is MSigDB's curated immune-specific collection (T/B cell, macrophage, dendritic, NK subsets, etc.); switching to it replaced ~707 largely non-immune signatures with 4,872 genuine immune-biology gene sets.

## KEGG Infection-Pathway Artifacts

Several KEGG pathways named for specific infections (*S. aureus* infection, HTLV-1 infection, Malaria, Leishmaniasis, etc.) appear enriched across contrasts. This is a known property of KEGG: its "infectious disease" pathways are built partly from host complement, coagulation, and immune-signaling genes. Real complement/inflammatory signal in the data can therefore make these pathways appear significant via gene overlap, with no implication of infection. The same pathways recur across both independent datasets, consistent with a reproducible statistical artifact. Each KEGG result table includes a `likely_gene_overlap_artifact` column flagging them.

## Extended Limitations

- **No cross-dataset group pairing beyond replication:** GSE25628-only and GSE7305-only groups cannot be directly compared; only same-comparison-type replication is valid.
- **Two different control definitions:** healthy-donor endometrium (GSE25628) vs. matched patient control (GSE7305). These differ biologically, but the difference is buried under the larger platform effect and cannot be isolated here.
- **Undefined `-M`/`-G` sample-ID suffix** in GSE7305 (not defined in the source paper; appears to be a within-patient sampling/processing label, symmetric across groups, so it does not confound the main comparison and is not used in the analysis).
- **Immune axis interpretation:** the dominant PC1 axis in GSE7305 immune scores near-perfectly separates lesion from control. Consistent with genuine immune biology, but a technical difference correlated with lesion status cannot be excluded without GEO-level processing metadata. Report as "a dominant, disease-associated immune axis," not "N independent significant signatures."
