"""Markdown tables for the README from the jsonl logs of lp_bound.py,
exact_certificate.py and lift_certificate.py."""

import glob
import json
from fractions import Fraction

ORDER = "1CRN,1SHG,1UBQ,1QYS,3CHY,1AKI,1A6M,2LZM,8ABP,1ADE,1GAI,1QOP,3PGK,4AKE".split(",")


def rows(pattern):
    out = {}
    for f in sorted(glob.glob(pattern)):
        for line in open(f):
            r = json.loads(line)
            if "iters" in r or "red_bound" in r or "soplex_s" in r:
                out[(r["protein"], r["chi"])] = r
    return out


lp = rows("experiments/lp_bound.[0-9].jsonl")
tight = rows("experiments/tighten.jsonl")
red = rows("experiments/exact_certificate.R.jsonl")
lift = rows("experiments/lift_certificate.*.jsonl")

print("| instance | rotamers | LP gap (a) | triplets | exact, reduced | exact, unpruned (lifted) | rounds | final set |")
print("|---|---|---|---|---|---|---|---|")
n_red = n_lift = 0
for chi in (1, 2, 3):
    for p in ORDER:
        k = (p, chi)
        a = lp.get(k)
        r = red.get(k)
        L = lift.get(k)
        gap = a["red_bound"]["gap"] if a else None
        gtxt = "≤1e-11" if gap is not None and gap < 1e-9 else (f"**{gap:.3g}**" if gap is not None else "–")
        rtxt = "✔ gap 0" if r and r["optimal"] and Fraction(r["gap"]) == 0 else ("✘" if r else "–")
        ltxt = "✔ gap 0" if L and L["optimal"] and Fraction(L["gap"]) == 0 else ("✘" if L else "not run")
        n_red += rtxt.startswith("✔")
        n_lift += ltxt.startswith("✔")
        print(f"| {p} χ≤{chi} | {a['rot'] if a else '–'} | {gtxt} | {r['triplets'] if r else '–'} | "
              f"{rtxt} | {ltxt} | {L['iters'] if L else '–'} | "
              f"{str(L['final_kept']) + ' / ' + str(L['labels']) if L else '–'} |")
print(f"\nreduced exact: {n_red}/42   unpruned exact: {n_lift}/42")
