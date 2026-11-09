# MG_191–MG_192: a matched Pedro–experimental interface case study

This is **one distinct protein pair and one AF3 prediction**, compared with repeated experimental copies and structures. The primary comparison uses the four main corresponding interfaces of 6RUT (2.65 Å): assembly 1 A:B and G-2:H-2; assembly 2 C:D and E-2:F-2. The copies are not independent targets, and no significance test is performed.

## Exact structures and definitions

Pedro observation `a063_j0007:r425`, job `250310_mgen_allbyall_1488`, model 0: MG_192 is chain C (1,053 residues), MG_191 is chain E (1,444 residues). The selected model was recovered from `AF3_completed_zips/folds_2025_04_06_16_31.zip`, member `250310_mgen_allbyall_1488/fold_250310_mgen_allbyall_1488_model_0.cif`; its SHA256 matches the previous BSA calculation. Pair ipTM is 0.76. It is a chain pair within a multichain prediction, not an independently predicted binary complex.

A contact is one residue pair whose minimum non-H/non-D atom distance is ≤8 Å. A close contact has **3 < d_min ≤4 Å**. BSA is (SASA_A + SASA_B − SASA_AB)/2, using only the selected two chains, the existing 1.4 Å probe and 100-point Shrake–Rupley implementation. All full-structure metrics were reused; the selected AF3 result was independently reconciled. Newly computed values are restricted to this target.

## Main result: separate sequence coverage from local geometry

The full AF3 chain pair has BSA **5,474.1 Å²**, **1,437** residue contacts, and **151** close contacts (**10.51%**). This cannot be compared without qualification with the truncated/partially observed 6RUT chains. After retaining the same reliably mapped, observed residue positions on both sides, **636–638 AF3 contacts** and **61–61 AF3 close contacts** lie outside the comparable residue scope and are no longer included. This is a scope sensitivity analysis, not deletion from the original results.

| 6RUT interface | Scope | BSA Exp / AF3 (Å²) | 8 Å count Exp / AF3 | 3–4 Å count Exp / AF3 | 3–4 Å fraction Exp / AF3 |
| --- | --- | --- | --- | --- | --- |
| 6RUT asm 1 A:B | Full | 2665.0 / 5474.1 | 758 / 1437 | 100 / 151 | 13.19% / 10.51% |
| 6RUT asm 1 A:B | Common residues | 2665.0 / 2790.9 | 758 / 801 | 100 / 90 | 13.19% / 11.24% |
| 6RUT asm 1 G-2:H-2 | Full | 2651.0 / 5474.1 | 771 / 1437 | 105 / 151 | 13.62% / 10.51% |
| 6RUT asm 1 G-2:H-2 | Common residues | 2651.0 / 2788.1 | 771 / 799 | 105 / 90 | 13.62% / 11.26% |
| 6RUT asm 2 C:D | Full | 2606.5 / 5474.1 | 758 / 1437 | 100 / 151 | 13.19% / 10.51% |
| 6RUT asm 2 C:D | Common residues | 2606.5 / 2790.9 | 758 / 801 | 100 / 90 | 13.19% / 11.24% |
| 6RUT asm 2 E-2:F-2 | Full | 2684.3 / 5474.1 | 760 / 1437 | 101 / 151 | 13.29% / 10.51% |
| 6RUT asm 2 E-2:F-2 | Common residues | 2684.3 / 2788.1 | 760 / 799 | 101 / 90 | 13.29% / 11.26% |

On common residues, the AF3-minus-experimental close-contact fraction difference ranges from **-2.35 to -1.96 percentage points** across the four copies, and the close-contact count difference ranges from **-15 to -10**. The common-residue BSA ratio (AF3/experimental) ranges from **1.039 to 1.071**. These are within-target descriptive ranges, not uncertainty intervals.

## What happens to experimental close contacts?

All experimental contact residues in the 13 investigated chain pairs have stable sequence correspondences. An absent experimental contact in AF3 is therefore not a missing-residue mapping failure. Exact AF3 distances were calculated for the corresponding residue pairs, including those beyond 8 Å; they were not censored or assigned an arbitrary distance.

| 6RUT interface | Exp close | Still (3,4] Å | AF3 ≤3 Å | AF3 >4 Å | Of these: >8 Å | New AF3 close pairs |
| --- | --- | --- | --- | --- | --- | --- |
| 6RUT asm 1 A:B | 100 | 71 | 8 | 21 | 0 | 19 |
| 6RUT asm 1 G-2:H-2 | 105 | 68 | 10 | 27 | 0 | 22 |
| 6RUT asm 2 C:D | 100 | 68 | 11 | 21 | 0 | 22 |
| 6RUT asm 2 E-2:F-2 | 101 | 68 | 10 | 23 | 1 | 22 |

For each copy, leaving the 3–4 Å bin includes both shortening and lengthening; only the >8 Å subset loses its 8 Å interface contact. New close contacts can partly offset these losses. Among experimental close contacts, 0–0 residue pairs per copy involve missing expected backbone/side-chain atoms. Repeating the minimum-distance comparison with atom names present in both structures retains 29–37 lost close contacts per copy. See the detailed transition table for atom identities, missing-atom flags, common-atom distances and inter-residue Cα distance changes.

These controls distinguish incomplete sequence coverage from changes in distances between corresponding atomic groups. They do not uniquely attribute each change to a side-chain rotamer, a backbone displacement, refinement uncertainty or assembly motion. Cα fits and inter-residue Cα distance changes are provided as diagnostics, not a causal classification.

## Which other experimental interfaces are comparable?

6YRK A:B and 8PBX/8PBY A:B provide additional corresponding views of the same protein pair. The common-residue fraction difference is −2.24 percentage points for 6YRK, but reverses direction for 8PBX (+0.65 points) and 8PBY (+1.20 points); AF3 has 5 and 10 more close contacts, respectively. Therefore, even within this protein pair, the direction depends on the experimental reference. The stable-mapping rule excludes ambiguous positions in 6YRK, while retaining all of its observed interface contact residues. 8PBZ A:B/C:D and 8PC0/8PC1 A:B are retained as exploratory comparisons with low-resolution coordinates (11/17/18 Å). 8PBZ A:D/B:C have partial overlap with the AF3 contact surface and are labelled alternative-interface diagnostics, not primary equivalent-interface matches. These labels distinguish contact surfaces and existing resolution limitations; no new numerical QC exclusion threshold was introduced.

Four 6RUT cross-copy chain pairs with **zero shared residue contacts** are excluded from genuine 1:1 interface comparisons. Six homomeric pairs are also outside the MG_191–MG_192 target. All 23 reference records and the reasons for their status remain in `interface_comparability.csv`. None of the 13 displayed comparisons represents another independent protein pair.

## Interpretation for discussion

The matched case can establish whether the same interface region and residue contacts are reproduced; it cannot estimate a general AF3 prediction bias. The high 6RUT 8 Å contact recovery (91.6–93.3%) indicates good recovery of the main contact surface, while the distance-bin transitions show differences in local atomic arrangements. Full-chain counts/BSA additionally reflect AF3 contacts outside the experimental coordinate scope. The common-residue analysis is essential to separate these effects. Sequence coverage is not a physical filter, and the restricted BSA values describe trimmed coordinate subsets: removing residues can also change solvent occlusion.

The existing severe-overlap flag is False, but that does not certify the absence of local clashes. Full and common-scope counts below 2 Å and 0.5 Å are retained in `interface_metrics.csv`. No structures were automatically excluded on a new distance threshold. This study does not measure binding affinity or demonstrate that every predicted interaction is biological.

## Figures

1. `01_6rut_full_and_common_residue_metrics`: full versus common-residue BSA, total contacts, close-contact counts and fractions. The repeated full AF3 points refer to the same model.
2. `02_6rut_experimental_close_contact_distances`: one point per experimental close residue pair, with its exact corresponding AF3 distance; separate panels for the four experimental copies.
3. `03_6rut_close_contact_turnover`: destination distance bins and gains/losses within the common-residue scope. Contact observations and copies are correlated.
4. `04_same_target_reference_sensitivity`: common-residue fractions and contact recovery across corresponding and exploratory structures of this single target.
5. `05_6rut_contact_correspondence_map`: shared and unmatched contacts in sequence space for the preselected 6RUT assembly 1 A:B copy, before and after correspondence restriction.

The selected PNG figures are included in this repository. The complete contact-level table is `experimental_close_contact_distances.csv`; `experimental_close_contacts_not_reproduced.csv` contains the subset outside the AF3 (3,4] Å bin. `6rut_assembly1_AB_largest_close_contact_shifts.csv` lists explicit atom-level examples without choosing the best experimental copy retrospectively.

## Local geometry and reference dependence: concrete observations

The four main 6RUT two-chain Cα RMSDs are 1.459–1.483 Å, consistent with overall structural correspondence. In contrast, 8PBZ A:D/B:C joint RMSDs are approximately 28.5 Å: these are alternative partial-contact surfaces, not equivalent primary 1:1 interfaces. 8PBX/8PBY joint RMSDs are approximately 1.34/1.51 Å; the opposite direction of their close-fraction difference therefore cannot simply be dismissed as an unrelated protein pair.

For 6RUT assembly 1 A:B, experimental A:ASN806 maps to MG_192/AF3 C:ASN807, and experimental B:GLY705 maps to MG_191/AF3 E:GLY705. The minimum changes from ND2–CA at 3.464 Å to N–CA at 5.302 Å, while the Cα separation changes from 5.949 to 5.890 Å. Expected heavy atoms are present on both sides; common-atom restriction leaves this result unchanged. This is evidence of a local atomic-group arrangement difference, not a missing sequence position. It does not uniquely identify a physical cause or energetic penalty.

The sole primary experimental close contact that moves beyond 8 Å is MG_192 position 529–MG_191 position 1189 in 6RUT assembly 2 E-2:F-2: experimental 3.987 Å, AF3 10.077 Å. Other experimental copies need not place this residue pair in the same close-contact bin. All transition values and atom identities are retained in the CSV, including shorter contacts.

Full AF3 has one residue contact below 2 Å (minimum 1.017 Å), but none of the four common-residue AF3 scopes contains it. It is outside the experimental coordinate overlap; this observation is retained as a QC caveat, not silently removed from the full analysis.
