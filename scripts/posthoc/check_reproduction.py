"""Compares every output in the given directory with scripts/posthoc/expected/ (CSV: byte-exact; JSON: equal after parsing)."""
import json, sys
from pathlib import Path

exp = Path(__file__).parent / "expected"; out = Path(sys.argv[1]); bad = 0
for e in sorted(exp.iterdir()):
    o = out / e.name
    if not o.exists():
        print("MISSING   ", e.name); bad += 1; continue
    same = (o.read_text() == e.read_text()) if e.suffix == ".csv" else (json.load(open(o)) == json.load(open(e)))
    print(("IDENTICAL " if same else "DIFFERENT ") + e.name); bad += (not same)
print("ALL OUTPUTS REPRODUCED" if not bad else f"{bad} OUTPUT(S) DIFFER"); sys.exit(1 if bad else 0)
