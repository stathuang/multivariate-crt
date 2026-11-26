# Validation record

Validation date: 20 November 2026 (Asia/Singapore). The original sections below describe the first local candidate; the revision-2 addendum records the separate review and GitHub packages. No manuscript source, remote repository, or cloud manuscript was changed. No Pro review or independent-agent review is claimed.

## Environment and isolation

The package was copied to fresh directories outside the research project. Commands used only package-relative code and explicitly supplied application inputs. A new Python 3.12 virtual environment was created and its pinned dependencies installed successfully from the package index. The first network attempt in the restricted environment failed; the authorized host-network installation succeeded. Python numerical checks ran in that new environment.

R 4.3.1 used the existing installed library. All 149 direct/transitive versions in `environment/R-dependencies.csv` matched. A fresh R library installation, container build, and complete historical source-download test were not performed. The runtime inventory is therefore evidence of the tested environment, not a guarantee for every platform.

The initial clean copy was used for computation trials. A second fresh copy of the final source set was used to repeat reconstruction, plotting, population calculations, raw-data preparation, and the modified MaxCoefT dispatcher. All application inputs and patient-derived validation objects stayed outside the deliverable archive.

## Completed checks

| Check | Result and scope |
|---|---|
| Active manuscript inventory | 42 loaded files, 27 figure/table environments, 31 external assets. Every loaded source hash was rechecked before packaging. The unloaded group-selection appendix was excluded. |
| Numerical table reconstruction | All 12 numerical tables recomputed from saved replication-level records. Each printed value agrees at its displayed precision. The 13th table defines statistics and is mapped to the implementations. |
| Graphical simulation curves | All plotted estimates and sample SEs agree with the adopted counts. Both CGM failures remain recorded; the affected cell uses 398 valid results out of 400. |
| Figure generation | All 20 expected PDF assets regenerated. Plot input CSVs accompany them. All PDF pages were rendered and visually inspected; no clipping or broken labels was observed. Layout differences are documented, and byte-for-byte PDF equality is not claimed. |
| Non-Gaussian MVN-CRT | All 12,000 formal p-values aggregated; 240 level-specific Wilson intervals match the adopted diagnostics. The independent implementation check and all 60 pilot datasets ran successfully. The 12,000 formal datasets were not rerun. |
| E.5 population approximations | Regenerated the original 500 backgrounds per design, with data hashes checked; leading differences 0.0407600449497 and 0.0944721688355, conditional-effect ranges, and variance 0.5625 reproduced. |
| E.8 scan norm | Precision-induced norm 0.0598150221680 reproduced. |
| Direction example | Complete 1,500-dataset, 999-copy protocol rerun; all 4,500 test rows match the adopted scores, p-values, rejection indicators, and seeds. |
| Independent-target centering | Complete 300-dataset, 19-copy protocol rerun; all 4,800 test rows match the adopted scores, p-values, rejection indicators, and seeds. |
| Connected-target centering | Independent one-dataset implementation check passed. The 100,000 adopted test rows, five-method rates, fixed-weight differences, and background-sign differences were reconstructed. The full 2,500-dataset experiment was not rerun. |
| Matched baseline | Eight independent smoke datasets ran with 999 copies; no failures. These smoke IDs are separate from formal evidence. |
| MaxCoefT | Pilot evaluated 28 formal-grid datasets with 999 copies. All data hashes, new p-values, and reused baseline p-values match the adopted records exactly. The final dispatcher also produces a combined raw JSONL file. The 7,000-dataset full grid was not rerun. |
| Gaussian benchmark | E01 and E02 trial modes passed, including the Ledoit–Wolf method. Formal saved counts were independently recomputed. |
| Main/extension entry trials | E03, E05, E07, E08, E09, E12, and E14 passed; all method p-values in these trials are finite and have status `ok`. These cover linear, nonlinear continuous, logistic, nonlinear binary, univariate, and t3 entry paths. Other full graphical experiments were not rerun. |
| Raw input preparation | Exact input hashes checked. Stock and ER-positive BCR scientific fields reproduce the adopted prepared objects. Only the stock provenance-description string differs. Small-class fitting warnings are retained. |
| Stock transformation | Reconstructed all 103 standardized columns from the saved price arrays; maximum difference from the accepted X matrix is zero. MMM/MO Box-Cox powers are 1.4/1.5. The fitted graph is independently reproduced during preparation. |
| Stock observed analysis | Full T1 observed protocol, including 1,000 copies and 1,000-tree forests, rerun. All six engine p-values match the adopted record exactly, including the four displayed methods. |
| Stock fitted null | First T1 fitted-null replication rerun at its full copy/tree counts. All five p-values match the adopted record exactly. The 2,802-analysis complete stock campaign was not rerun. |
| Stock graph/goodness of fit | Graph renderer ran. Maximum/median degrees, density, and all saved GoF values checked; F-sum rank p-value recalculated from 1,000 saved randomized statistics. A fresh 1,000-copy GoF calculation was not run. |
| BCR observed | Adopted RF/dCRT random stream replayed and certified LRT recalculated. All four p-values match exactly. |
| BCR fitted null and RF repair | Replicates 1 and 74 rerun, including an original RF empty-level failure. All four adopted p-values match within 3.4e-16. The full 400-replication experiment was not rerun. |
| BCR likelihood evidence | Extracted scalar summaries for the observed fit and 400 adopted null fits. All full/reduced likelihood gaps are at most 1e-5, exact-face/dual flags pass, and the corresponding LRT p-values agree. Patient-derived certificate arrays were not included. These summaries record existing certification; they are not 401 fresh likelihood fits. |
| Fresh-run import | A completed direction run was imported into a separate reference tree, then all numerical summaries rebuilt successfully. The bundled reference remained unchanged. |
| Source and archive checks | Source adaptations recorded as diffs; Python/R source syntax checked; local absolute paths and application input matrices excluded from the clean distribution; final checksum manifest and ZIP integrity checked. |

Machine-readable evidence is in `validation/`, `outputs/rebuilt/verification.json`, `outputs/rebuilt/population.json`, and `outputs/figures/MANIFEST.csv`. The author retains the larger temporary computation outputs separately. Rejected trial preparations and transient adapter errors were repaired before the final checks; no failed method record was converted to a successful result or dropped from the adopted data.

## Changes made for this distribution

The package reuses the established implementations, formal outputs, and repair records. New adapters provide explicit package-relative paths, a single dispatcher, the preserved BCR RNG/repair sequence, exact snapshot preparation, result import, independent aggregation, and figure rendering. The shared stock transformation functions were recovered from the preserved source. Duplicate unused engine files were removed. The original experiment files in the research project were not edited.

Development checks caught and fixed a missing helper import in raw preparation, table-parser conventions, null-string parsing in pandas, the output collection step in the MaxCoefT dispatcher, and selection of the adopted observed-class degrees of freedom when extracting historical likelihood summaries. The population check also caught an initially incorrect background-penalty setting before it was changed to the actual 0.2n protocol. The final results above use the corrected code and verified source definitions.

## Limits and remaining actions

Full simulations were reconstructed from adopted records, with selected full reruns and entry trials as listed above. This is not a claim that every original experiment was rerun in a new operating system. The historical stock campaign alone cost about 49.5 process-hours; this validation did not repeat that workload. Fresh R installation remains unverified.

The first delivery left source and data-access questions open. The revision-2 addendum and DCRT_ADAPTATIONS.md supersede that status: the historical code source and six public input routes have been verified. Public redistribution of the dCRT adaptation and the public DIA archive location remain separate from the complete review package. Author code reuse terms and the final repository/archive address are author choices. No publication or manuscript-link change has been performed.


## Revision 2: review and GitHub distributions

The review archive now includes all seven raw application inputs and the historical dCRT script for optional inspection. The separate GitHub candidate excludes these raw files and the full historical script, while retaining the implemented comparator and its documented adaptation patch. Both use the same adopted numerical code. The input exclusion described in the original validation above applies to the first distribution and the GitHub candidate, not to the new review archive.

New checks on 20 November 2026:

- Downloaded all six public input files from the fixed GGMtest revision and NCBI. Their SHA-256 values and sizes agree with `data/required-inputs.json`. The seventh file is the exact locally preserved DIA snapshot, included in the review archive and separately verified.
- Ran the review README's offline input check and raw preparation in a new directory outside the research workspace. All 28 fields of the resulting stock and breast-cancer objects agree, at zero comparison tolerance, with the previously validated prepared objects. R reported the existing preparation warnings; no failed preparation was treated as a pass.
- Extracted the two dCRT functions from the pinned 2021 source, applied the supplied patch, and obtained exactly the same extracted source as the historical adopted file.
- Compared all three restricted dCRT implementations, with Gaussian, binomial, and RF response fits, against both the historical file and the patched functions. All nine p-values agree exactly in this check, and all nine final random-number states agree. The direct historical check was also repeated in the clean review directory. The package includes `validation/revision-2/dcrt-parity.csv` and the preparation comparison.
- Verified that every pre-existing R and Python computation file is unchanged from the first validated package. The new scripts provide explicit downloading, function extraction, and the adaptation check; they do not replace the manuscript's simulation or analysis functions.

The first draft of the adaptation checker lacked two helper imports, and an initial preparation invocation used the wrong working directory. Those setup attempts failed and were corrected before the successful clean-directory checks. This revision does not repeat the full simulations, figure rendering, or fresh R installation; the earlier validation and its limits remain applicable because the numerical implementation and adopted result records are unchanged.
