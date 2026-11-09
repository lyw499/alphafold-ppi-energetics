# Within-type 3–4 Å contact fractions

Numerator: unique residue-pair contacts with `3 < d_min ≤ 4 Å`. Denominator: all unique residue-pair contacts of that contact type with `d_min ≤ 8 Å`. Values are pooled across contacts, not per-interface averages. Mgen includes only the 23 nonzero interfaces in the existing reference table; the sensitivity subset excludes every interface from 8PBZ, 8PC0, and 8PC1 (15 interfaces remain).

| Dataset | Contact type | 3–4 Å contacts (numerator) | ≤8 Å contacts (denominator) | Within-type fraction |
|---|---|---:|---:|---:|
| Pedro AF3 | BB-BB | 26,908 | 695,115 | 3.87% |
| Pedro AF3 | BB-SC | 193,034 | 3,504,759 | 5.51% |
| Pedro AF3 | SC-SC | 332,667 | 2,471,790 | 13.46% |
| FoldBench Experimental | BB-BB | 784 | 6,554 | 11.96% |
| FoldBench Experimental | BB-SC | 1,843 | 18,468 | 9.98% |
| FoldBench Experimental | SC-SC | 2,477 | 10,517 | 23.55% |
| SAFES Experimental | BB-BB | 133 | 1,502 | 8.85% |
| SAFES Experimental | BB-SC | 478 | 4,714 | 10.14% |
| SAFES Experimental | SC-SC | 604 | 2,706 | 22.32% |
| Mgen Experimental (all nonzero interfaces) | BB-BB | 28 | 1,087 | 2.58% |
| Mgen Experimental (all nonzero interfaces) | BB-SC | 525 | 5,857 | 8.96% |
| Mgen Experimental (all nonzero interfaces) | SC-SC | 579 | 3,584 | 16.16% |
| Mgen Experimental (excluding 8PBZ, 8PC0, 8PC1) | BB-BB | 13 | 723 | 1.80% |
| Mgen Experimental (excluding 8PBZ, 8PC0, 8PC1) | BB-SC | 386 | 3,356 | 11.50% |
| Mgen Experimental (excluding 8PBZ, 8PC0, 8PC1) | SC-SC | 355 | 1,975 | 17.97% |

The two Mgen rows are overlapping sensitivity populations and must not be interpreted as independent datasets. All bars show pooled-contact fractions; they do not give equal weight to each interface.
