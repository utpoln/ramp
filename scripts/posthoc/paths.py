"""Input locations for the post-hoc analyses (no machine-specific paths).

REPO            the RAMP repository root (two levels above this file)
RAMP_GO_OBO     go-basic.obo, GO release 2026-05-19 (SHA-256 923645000afcb1df...);
                default: <REPO>/reference_data/go-basic_2026-05-19.obo
RAMP_PROBE_RESULTS  the PROBE development responses (probe_results/ of github.com/utpoln/PROBE);
                default: <REPO>/../PROBE/probe_results  (a sibling clone of the PROBE repository)
"""
import hashlib, os, sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
RAGO = REPO  # name kept from the original analysis code
OBO = Path(os.environ.get("RAMP_GO_OBO", REPO / "reference_data" / "go-basic_2026-05-19.obo"))
PR = Path(os.environ.get("RAMP_PROBE_RESULTS", REPO.parent / "PROBE" / "probe_results"))
# The frozen configuration recorded only a 16-character SHA-256 prefix (go_obo_sha256_16); the full
# SHA-256 below is that of the archived release file distributed in reference_data/.
FROZEN_GO_PREFIX16 = "923645000afcb1df"
GO_SHA256 = "923645000afcb1df0e72e252b64950b1c2bcc98c9393fec8c1dd1a0f7c9d4d77"
assert GO_SHA256.startswith(FROZEN_GO_PREFIX16)


def require_obo():
    """Stop unless the GO file is the distributed release file (whose SHA-256 matches the frozen prefix)."""
    if not OBO.exists():
        sys.exit(f"GO file not found at {OBO}; see scripts/posthoc/README.md (set RAMP_GO_OBO)")
    h = hashlib.sha256(OBO.read_bytes()).hexdigest()
    if h != GO_SHA256:
        sys.exit(f"GO file {OBO} has SHA-256 {h}, not {GO_SHA256} (the distributed release file, matching the frozen 16-character prefix {FROZEN_GO_PREFIX16})")


def require_probe():
    if not (PR / "mistral-7b.jsonl").exists():
        sys.exit(f"PROBE responses not found at {PR}; clone github.com/utpoln/PROBE next to this "
                 "repository or set RAMP_PROBE_RESULTS")
