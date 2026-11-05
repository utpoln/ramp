# Gene Ontology file governed by the frozen protocol

| Field | Value |
|---|---|
| File in this directory | `go-basic_2026-05-19.obo` |
| Source URL | https://release.geneontology.org/2026-05-19/ontology/go-basic.obo |
| GO release (file header) | `data-version: releases/2026-05-19` |
| Size | 32,212,518 bytes |
| SHA-256 | `923645000afcb1df0e72e252b64950b1c2bcc98c9393fec8c1dd1a0f7c9d4d77` |
| Frozen record | `frozen_experiment_config.json` → `go_obo_sha256_16` = `923645000afcb1df` (16-character prefix match; no full SHA-256 was recorded at the time) |
| `[Term]` stanzas | 48,321 (10,058 obsolete) |
| Server `Last-Modified` | Thu, 28 May 2026 21:45:18 GMT |
| Downloaded and hashed | 3 November 2026 |

The archived file is byte-identical (full SHA-256) to the `go-basic.obo` copy retained in the
study's working directory. When the experiments were run, only a 16-character prefix of the
SHA-256 was recorded (one of six 16-character prefixes in the frozen configuration). The file read
at that time is therefore matched by that prefix, not by a full-hash comparison. The frozen score and the
primary endpoint do not read this file. It is needed only by the post-hoc identifier-validity and
ancestor-credit analyses (`scripts/posthoc/`), which refuse to run if the file's full SHA-256 differs from the value above.

The release label `2026-03-25` written in the frozen files is incorrect; see `GO_RELEASE_CORRECTION.md`.

## Licence and attribution

The Gene Ontology is released under the Creative Commons Attribution 4.0 International licence
(CC BY 4.0, https://creativecommons.org/licenses/by/4.0/), as stated in the file's own header
(`property_value: terms:license http://creativecommons.org/licenses/by/4.0/`). The file is
redistributed here unmodified, with attribution to the Gene Ontology Consortium
(https://geneontology.org; release archive https://release.geneontology.org/2026-05-19/).
Please cite the Gene Ontology as described at https://geneontology.org/docs/go-citation-policy/.

## Repository constraints

32.2 MB is below GitHub's 50 MB warning threshold and 100 MB hard limit for a single file, so
Git LFS is not required. The scripts verify the SHA-256 above before reading the file
(`scripts/posthoc/paths.py`, `require_obo()`) and stop if it differs or is missing.
