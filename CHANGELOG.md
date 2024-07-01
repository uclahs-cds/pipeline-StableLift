# Changelog
All notable changes to the pipeline-name pipeline.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

This project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [Unreleased]
### Added
- Add `sample_id` extraction from BAM
- Add template input YAMLs
- Add pipeline-Nexflow-config as submodule and redirect set_resources_allocation
- Add pipeline-Nextflow-module as submodule
- Additional out of memory exit code
- Pipeline release action
- Template for NFTest testing results in PR template
- Enable dependabot
- Add example PlantUML image to README
- Add workflow to build documentation
- Add workflows to run Nextflow configuration tests

### Changed
- Switch resource limit checks to external scripts
- Update links in on-prem Confluence to point to cloud-based Confluence
- Fix `CODEOWNERS` file
- Use `schema.check_path` for `workDir` validation
- Add `Discussions` and `Contributors` to the Table of Contents in `README.md`
- Update from DSL1 to DSL2
- Standardize config structure
- Restructure repo so main script is main.nf
- Reorganize contributors and metadata
- Reorganize PR template so description is at top
- Update automatic node detection to allow for F2 detection
- Update Issue Template
- Standardize input/output/parameter structure in README
- Avoid modification of input parameter `output_dir`
- Create default docker container registry parameter for tools
- Use `methods.setup_process_afterscript()` to capture log files

---

## [1.0.0] - YYYY-MM-DD
### Added
- For new features.
- Added item 1.

### Changed
- For changes in existing functionality.
- Changed item 1.

### Deprecated
- For soon-to-be removed features.

### Removed
- For now removed features.
- Removed item 1.

### Fixed
- For any bug fixes.
- Fixed item 1.

### Security
- In case of vulnerabilities.
