# Changelog
All notable changes to the StableLift pipeline.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

This project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [1.1.0] - 2024-10-17

### Added

- Add Delly2-sSV support
- Add preprocessing script for Strelka2 VCFs to add GT field for Funcotator
- Add process to split large VCFs to parse in chunks

### Changed

- Standardize annotation workflow to always annotate variants after LiftOver
- Rename Delly2 to Delly2-gSV
- Upgrade to score version 1.20

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
