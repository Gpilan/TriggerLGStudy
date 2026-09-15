# Code versions and simulation records

- `main`: integrated code after relevant build and regression checks.
- `fix/<topic>`, `study/<topic>`, `docs/<topic>`: one active task per branch; remove the branch after its changes have been incorporated or deliberately archived.
- `snapshot/<original-version>`: historical branch-tip commits. These are records, not certifications of correctness. Do not move existing tags.
- `baseline/YYYY-MM-DD`: a tested code reference, accompanied by its validation scope. Add a suffix if more than one baseline is needed on the same day.

LG/noLG, beam positions, event counts and random seeds are configuration values, not branch names. Use a new run directory and manifest for each condition. Record the commit and any dirty patch, configuration, source/binary hashes, dependencies, analysis settings, event denominator and extraction failures. A tag alone does not archive external data or guarantee reproducibility.

Before retiring a branch, preserve its exact tip with a remote tag and verify the SHA. A branch with commits absent from main may be archived without merging obsolete physics back into the current model. Keep the reason for archival in the snapshot index.

Keep code, run configuration, meaningful tests, documentation and selected small reference results in Git. Keep automatic caches, raw ROOT data and bulk run output in their designated data storage with a manifest. Adding an ignore rule does not remove previously tracked files. Review existing result dependencies before changing their tracking.

Personal assistant context and local maintenance backups are not publication artifacts. They are excluded from Git. Git version management is not a backup of external or ignored research data.
