---
description: Ingest every .md paper under a directory into references/, one isolated subagent per paper.
argument-hint: <dir-or-glob> [--force]
allowed-tools: Glob, Read, Bash, Task
---

# Ingest papers

Ingest papers into the `references/` database.

Target: `$ARGUMENTS` (a directory, a glob, or a single .md path). If a `--force` flag is present, pass it through so existing papers are reprocessed.

Do this:

1. Enumerate the target with Glob (`**/*.md`), excluding anything already under `references/papers/` and excluding `references/CONVENTIONS.md`. List what you found and how many.
2. For **each** file, delegate to the `paper-ingestor` subagent with that single file path (and `--force` if requested). Run them as independent tasks - each paper gets its own context window, so batch size doesn't dilute attention. Prefer parallel/background delegation for throughput; keep going if one fails.
3. Collect the one-line result from each subagent. Do NOT pull their intermediate work into this context - only their status lines.
4. Print a summary table: file -> status (created / updated / routing-only / skipped-exists / error) -> key. Then counts per status.
5. If any errored, list those paths so they can be retried; don't stop the batch for a single failure.

Do not clean, extract, or edit papers yourself - that is the subagent's job. Your role is orchestration and the final summary. Read the final `references/index.yaml` file and check for any inconsistencies/inaccuracies/contradictions.
