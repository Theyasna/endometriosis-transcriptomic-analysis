# Dataset Notes

## GSE25628

- **Platform:** Affymetrix Human Genome U133A 2.0 Array (**GPL571**) — verified directly from the GEO series record
- **Samples:** 22 total — Ectopic lesion (8, samples g–p), Eutopic (8, samples q–z), Normal (6, samples a–f)
- **Tissue types:** Ectopic = endometriotic lesion tissue (ovary or pelvic peritoneum); Eutopic = endometrium from affected women; Normal = endometrium from healthy donors.
- **Menstrual cycle control:** All samples collected in the **proliferative phase** — confirmed from the GEO overall design. Cycle phase is therefore controlled within this dataset.
- **Study aim (from GEO):** "To identify genes potentially involved in growth and maintenance of the ectopic endometrium." Published: Crispi S et al., *J Cell Physiol* 2013;228(9):1927–34. PMID: 23460397.
- **Group assignment:** Extracted from the sample title field (source_name_ch1 was uninformative). Title-based assignment (g–p = ectopic, q–z = eutopic, a–f = normal) matches the GEO overall design exactly.
- **Probes:** 22,277.

## GSE7305

- **Platform:** Affymetrix Human Genome U133 Plus 2.0 Array (GPL570).
- **Samples:** 20 total — "Diseased" (10, E- prefix), "Normal"/control (10, N- prefix).
- **CRITICAL CORRECTION — what these groups actually are:** In the source publication (Hever et al., *PNAS* 2007), the "Diseased" (E-) samples are **ovarian endometriosis lesions** (ectopic tissue), and the "Normal" (N-) samples are **matched control endometrium from the same patients, collected at the same time**. This is a **paired within-patient, lesion-vs-control** design. Neither group is "eutopic endometrium from patients vs. eutopic from healthy controls" — an earlier version of this project's documentation described it that way in error; that description has been corrected throughout.
- **Why the pairing matters:** the original study explicitly compared lesion to same-patient control endometrium *specifically to remove menstrual-cycle differences* as a confounder. This makes GSE7305 a lesion-vs-control comparison directly analogous to GSE25628's ectopic-vs-normal — which is the basis for the cross-cohort replication in this project.
- **Sample-ID suffix (`-M` / `-G`):** each patient's samples carry an `-M` or `-G` suffix. The source publication (Hever et al. 2007) does not define these labels — the design is described only as ovarian endometriosis lesions paired with matched control endometrium from the same patients, with no myometrium or tissue-compartment split mentioned for the endometriosis comparison. The suffix pattern is symmetric between the Diseased and Normal groups for each patient (the same patient numbers, 01–08, appear in both groups), indicating it is a within-patient sampling or processing label rather than a disease variable. It is not used in this analysis and cannot confound the lesion-vs-control comparison.
- **Log transformation:** Phase 0 detected raw (non-log) values and applied log2(x + 0.5). GSE25628 was already log-transformed.
- **Probes:** 54,675.

## Cross-Dataset Design

After merging, the batch × group cross-tabulation is:

```
              Diseased  Ectopic  Eutopic  Normal
GSE25628             0        8        8       6
GSE7305             10        0        0      10
```

Only "Normal" appears in both datasets. Because the two datasets use different platforms and were processed in different labs, samples are **not pooled across datasets**. Instead, each dataset's lesion-vs-control comparison is computed independently and the resulting DEG lists are overlapped for replication.

### A note on the two "Normal" groups

The two control groups are **not identical in kind**:
- GSE25628 "Normal" = endometrium from **healthy donors** (no endometriosis).
- GSE7305 "Normal" = **matched control endometrium from the endometriosis patients themselves**.

Empirically, the 16 normal samples separate almost perfectly by dataset (PC1 ≈ 62% variance, Mann-Whitney p ≈ 0.00025). This separation is dominated by the platform/processing difference (GPL571 vs GPL570, different labs) rather than the subtle biological difference between healthy-donor and patient-control endometrium — the two effects cannot be disentangled here. Because the per-dataset design never pools these two normals, this difference does not contaminate any reported comparison; it is documented as a limitation, not a source of error.

## Common Probes

22,277 probes were shared between GPL571 and GPL570 and used for all downstream analyses.
