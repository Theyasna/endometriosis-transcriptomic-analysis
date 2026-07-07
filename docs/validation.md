# Validation Notes

## Cross-Cohort Replication (the headline result)

Both datasets contain a lesion-vs-control comparison. Overlapping the two independently-derived significant-DEG lists:

- GSE25628 Ectopic lesion vs Normal: 3,604 DEGs
- GSE7305 Ovarian lesion vs Control: 1,603 DEGs
- **Overlap: 608 genes (37.9% of the smaller set)**

This is a clean replication: two independent patient cohorts, two different Affymetrix platforms, and two different control definitions (healthy-donor endometrium vs. matched patient control). The overlap is computed from independently-generated gene lists, never from pooled samples. Robustness to the control-definition difference makes the shared signal more credible, not less.

## The Two "Normal" Groups Are Different — and Why It Doesn't Break Anything

The two control groups are not the same kind of tissue:
- GSE25628 "Normal" = healthy-donor endometrium
- GSE7305 "Normal" = matched control endometrium from endometriosis patients

Testing the 16 normal samples directly: they separate almost perfectly by dataset along PC1 (≈62% of variance; Mann-Whitney p ≈ 0.00025), and 31.7% of probes differ by more than 2-fold between the two groups. This separation is far larger and more global than the expected biological difference between healthy-donor and patient-control endometrium — it is dominated by the platform/lab difference (GPL571 vs GPL570). The biological difference is real but cannot be isolated from the technical one here.

**Why this is not a problem for the results:** the per-dataset design never pools the two normal groups. Each dataset's lesion comparison uses only its own control. The two-normal difference is a documented limitation, not a source of error in any reported comparison.

## PCA Diagnostic: Immune Signature Axes

PCA on the ssGSEA enrichment score matrix within each dataset:

- **GSE25628** (ectopic/eutopic/normal, n=22): PC1 = **60.9%** of variance; PC1 differs across groups (Kruskal-Wallis p = 0.0165).
- **GSE7305** (lesion/control, n=20): PC1 = **75.6%** of variance; PC1 near-perfectly separates lesion from control (Mann-Whitney p = 0.000246, near the theoretical minimum for n=10 per group).

**Interpretation:** the immune-signature hit rates (35.2% and 72.7%) are driven by one dominant, correlated axis per dataset, not thousands of independent findings. C7 ImmuneSigDB gene sets are highly redundant, so a single strong biological signal pushes many correlated sets past significance at once.

**Unresolved caveat:** PC1 correlating with group status is necessary but not sufficient to prove the axis is purely biological. A technical difference within GSE7305 correlated with lesion status would produce an identical signature. GEO-level processing metadata was not available to test this. Report accordingly.

## KEGG Artifact Cross-Dataset Reproducibility

Pathways flagged as likely gene-overlap artifacts recurred across independent datasets:

| Pathway | GSE25628 ectopic vs normal | GSE25628 eutopic vs normal | GSE7305 lesion vs control |
|:---|:---:|:---:|:---:|
| Human papillomavirus infection | flagged | flagged | flagged |
| Human T-cell leukemia virus 1 infection | flagged | — | flagged |
| Malaria | flagged | — | flagged |
| Leishmaniasis | flagged | flagged | — |
| Rheumatoid arthritis | flagged | flagged | flagged |
| Staphylococcus aureus infection | flagged | — | flagged |
| Viral myocarditis | flagged | — | flagged |

Recurrence across two independent cohorts and platforms is consistent with these being reproducible statistical artifacts (shared host complement/immune gene membership), not infection biology.
