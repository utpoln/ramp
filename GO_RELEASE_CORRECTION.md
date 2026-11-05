# Correction: Gene Ontology release label

**Date of correction:** 31 October 2026
**Affects:** the GO release label in `FROZEN_PROTOCOL.md` (Sec. 6, "GO ontology" row) and `frozen_experiment_config.json` (`go_obo_release`)
**Does not affect:** any input file, any recorded 16-character SHA-256 hash prefix, score, cohort, analysis, or result

## What is wrong

Both frozen files record the GO ontology as `go-basic.obo`, **release 2026-03-25**. That label is incorrect. The file actually governed by the frozen protocol is **release 2026-05-19**.

## Evidence

- **Recorded hash prefix.** The frozen configuration identifies the ontology file by a 16-character SHA-256 hash prefix, `"go_obo_sha256_16": "923645000afcb1df"`. This is one of six 16-character prefixes recorded in that file: the homology reference database (FASTA), the GO ontology file, and four scripts (`draw_heldout_sample.py`, `run_blast_retrieval.py`, `llm_prompting_rag.py`, `evaluate_rag.py`). No full SHA-256 was recorded for any of them, and no hash was recorded for the other governed inputs.
- **The file whose SHA-256 begins with that prefix** has the following properties:
  - SHA-256: `923645000afcb1df0e72e252b64950b1c2bcc98c9393fec8c1dd1a0f7c9d4d77`
  - Header: `format-version: 1.2` and `data-version: releases/2026-05-19`
  - 48,321 `[Term]` stanzas, of which 10,058 are marked obsolete. This matches the corrected term count already recorded in both frozen files.
  - Size: 32,212,518 bytes.
- **The archived release matches.** The file at <https://release.geneontology.org/2026-05-19/ontology/go-basic.obo> was downloaded on 3 Nov 2026. Its full SHA-256 is the value above, and it has the same size and header as the copy of `go-basic.obo` kept in the study's working directory (file date 3 June 2026). The two files are byte-identical. A copy is included as `reference_data/go-basic_2026-05-19.obo`; see `reference_data/GO_FILE_VERIFICATION.md`.
- **What this does and does not establish.** When the experiments were run, only the first 16 hexadecimal characters of the ontology file's SHA-256 were recorded (`go_obo_sha256_16`). No full SHA-256 of the file read at that time was recorded. The full-hash comparison above is therefore between the archived release and the retained working-directory copy. The link to the file read during the experiments rests on the 16-character (64-bit) prefix match.
- **The 2026-03-25 release** (`data-version: releases/2026-03-25`) has a different SHA-256, beginning `a77e356737dab39a`, and 48,291 `[Term]` stanzas. It does not match the recorded hash.
- **How the error arose.** When the configuration was generated, the hash prefix was computed from the ontology file itself. The release string and the original term count of 48,291 were entered as fixed text. They match, and appear to have been carried over from, documentation describing the PROBE repository's own copy of the ontology, which is the 2026-03-25 release. The term count was later corrected to 48,321 (`FROZEN_PROTOCOL.md`, Sec. 15). The release label was not corrected at that time.

The determination rests on the 16-character prefix match together with the file's own `data-version` header, not on the term count.

## Why the frozen files are not corrected in place

`frozen_experiment_config.json` is unmodified since it was frozen. `FROZEN_PROTOCOL.md` has received one later edit, unrelated to this correction: on 18 October 2026, after the held-out run, three lines were changed to rename the project from RAGO to RAMP (the title, Sec. 3 and Sec. 5; commit `c657d0c`). The edit changed no rule, number, or threshold. The last version before it is in commit `44fcdff`. That protocol file's own statement that it is "retained unmodified" predates the rename.

Neither file is edited to correct the GO label, so that both remain the record of the prespecified design. This note supersedes only their GO release label.

## Impact

None on results. The frozen configuration identifies the ontology file by a 16-character SHA-256 hash prefix. The configuration records the same prefix before and after the term-count correction. The frozen score and primary endpoint do not read the ontology file: claims are matched to ground truth by exact GO identifier within a namespace. Only the post-hoc identifier-validity and ancestor-credit analyses (`scripts/posthoc/`) read it, and they verify its hash before running.

## Access date

The date on which the ontology file was originally downloaded was not recorded and could not be recovered from the available project records, so no access date is given for it.
