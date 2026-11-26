# dCRT source record

The historical script attributes Molei Liu and https://github.com/moleibobliu/Distillation-CRT. Its local original SHA-256 is `c8df5c7717dc9ea525e42979a3df179b8ee51168b5910108cc40f996790bb864`. The public GGMtest copy differs only in its final newline and is fixed at https://github.com/xiaotong0201/GGMtest/blob/35a3647910792d39e27986b05e07b04b0a6f6461/GGM-Test-Version0/CRT/dCRT_function.R. The same historical content is already present in the GGMtest revision dated 21 December 2023.

This package reconstructs the no-screening, no-interaction, estimated Gaussian conditional, analytic-tail branch. The complete implementation history and executable adaptation recipe are in ../DCRT_ADAPTATIONS.md. The current upstream update was dated 2 January 2024; no 2026 update was found in the inspected history. The RF difference predates that update.

The original adjacent group wrapper `CRT-Stat.R` has SHA-256 `3d3dff7f41a661e398e03f6e1369b9c69e32c16d54ea9989d681c1f9bdebae5d`. Group combination is Bonferroni. External packages are called through installed namespaces. No explicit source-distribution licence was verified; the review archive and the public candidate are distinguished in THIRD_PARTY.md.
