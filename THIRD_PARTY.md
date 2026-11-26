# Attribution and distribution scope

## Review package

The separate `jrssb-review` archive includes the adopted comparator implementations, all seven exact application inputs, and a historical dCRT source copy for inspection. Reviewers do not need to download application data or dCRT. Installed R and Python dependencies are still required. This packaging decision addresses review reproducibility; it does not declare a new licence over third-party code or data.

## Private GitHub repository

This repository contains the same analysis code, saved computational records, source documentation, and input checksums. Raw application inputs and the full historical dCRT script are excluded. Restricted dCRT implementations and the documented adaptation patch are still present: this is not an implementation based solely on an external dependency. The repository is private and has not been assigned a blanket open-source licence.

The dCRT source is attributed to Molei Liu's Distillation-CRT repository and the paper by Liu, Katsevich, Janson, and Ramdas. Its fixed historical source, upstream history, local changes, and reproducible adaptation commands are in DCRT_ADAPTATIONS.md. No explicit upstream software licence was found in the inspected sources. Public redistribution of this implementation remains a separate decision; source acknowledgement or a patch does not itself establish permission.

The six historical/public inputs have verified external acquisition routes in README.md. The public location for the exact DIA.json remains to be filled after the author chooses an archive. The exact snapshot is already included for review. No additional long data-source statement is required in this code package.

## Author code and installed packages

The author can separately choose reuse terms for their own code. That choice is distinct from the final repository URL and must not purport to grant third-party rights. R and Python libraries are installed dependencies, not bundled source distributions. Their recorded versions and R licence fields are in environment/. The eventual public URL and DOI belong in the submission/publication statements when available.
