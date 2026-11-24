# Analysis data

The tracked files here are the Hallmark gene-set definition and a manifest for all source snapshots. Raw data and participant-level derived outputs are rebuilt locally from public data and are ignored by Git. The source accession, mapping rules and survival fields are documented in [data provenance](../docs/data_provenance.md).

`Rscript scripts/fetch_inputs.R` downloads the named source snapshots. `Rscript scripts/run_pipeline.R` verifies the hashes, runs the analysis and generates publication tables and figures. Copying existing source snapshots into the listed paths is equivalent if their size and checksum match.
