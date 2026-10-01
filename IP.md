# IP registrations

Evidence of authorship registered through [iptrust.org](https://iptrust.org) with
`iptrust2`. Each record commits, hash-only, to the repository's tracked files at
the tagged commit; it is anchored on BSV (NotaryHash) and Bitcoin
(OpenTimestamps). A record is evidence that the content existed and that the
signer made the statement. It is not a registered right.

| date | covers | type | record hash | signer | tag |
|---|---|---|---|---|---|
| 2026-09-30 | Exact LP-dual optimality certificates: local-polytope LP + triplet tightening, SoPlex exact duals, Fraction checker, certificate lifting to unpruned instances (42/42 reduced, 38/42 unpruned) | research | [`5bf36e96f735fea98535c18d93db32375668a33b04bb1abb53147c76d296acd9`](https://iptrust.org/v/5bf36e96f735fea98535c18d93db32375668a33b04bb1abb53147c76d296acd9) | Gregory J. Ward | `ip-2026-09-30` |

## Verifying

The private evidence (`.iptrust/v2/<recordHash>/`, with `opening.json` and
`manifest.json`) is kept off this public repository, on the authors' machine.
With it, verify in a separate worktree so the main working copy is untouched:

```bash
git worktree add /tmp/at-tag ip-2026-09-30 && cp -r .iptrust /tmp/at-tag/
(cd /tmp/at-tag && iptrust2 verify 5bf36e96f735fea98535c18d93db32375668a33b04bb1abb53147c76d296acd9 --fetch-headers)   # expect: valid, matched
git worktree remove --force /tmp/at-tag
```
