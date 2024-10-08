# Changelog
All notable changes to the StableLift pipeline.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

This project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [1.0.0] - 2024-09-10

### Added

- Add workflow for SNV callers (Mutect2, HaplotypeCaller, Strelka2, Muse2, SomaticSniper)
- Add workflow for SV caller (Delly2)
- Add pipeline diagram
- Add reverse liftover (GRCh38 -> GRCh37) for SNV branch
- Add reverse liftover (GRCh38 -> GRCh37) for SV branch
- Add optional `target_threshold` and `target_specificity` parameters
- Add NFTest cases for all variant callers and directions

### Changed

- Sort VCF after liftover in SV branch
- Hide parametric complexity behind simpler user inputs
- Downgrade to score version 1.16
