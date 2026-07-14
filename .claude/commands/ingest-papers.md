---
description: Ingest every .md paper under a directory into references/, one isolated subagent per paper.
argument-hint: <dir-or-glob> [--force]
allowed-tools: Glob, Read, Bash, Task
---

# Ingest papers

Ingest papers into the `references/` database. Your role is orchestration and the final summary; you
MUST NOT clean, extract, or edit papers yourself - that is the `paper-ingestor` subagent's job.

Target: `$ARGUMENTS` (a directory, a glob, or a single .md path). If a `--force` flag is present, pass
it through so existing papers are reprocessed.

Steps SHALL run in this order:

1. Enumerate the target with Glob (`**/*.md`), excluding anything already under `references/papers/`
   and excluding `references/CONVENTIONS.md`. List what you found and how many.
2. For **each** file, delegate to the `paper-ingestor` subagent with that single file path (and
   `--force` if requested). Each paper SHALL run as its own task with its own context window, so batch
   size does not dilute attention. You SHOULD delegate in parallel for throughput. A single failure
   MUST NOT stop the batch; record the error and continue.
3. Collect the one-line result from each subagent. You MUST NOT pull their intermediate work into this
   context - only their status lines.
4. Read the final `references/index.yaml` and check it for inconsistencies, inaccuracies, or
   contradictions introduced by the batch (duplicate keys, entries whose file is missing, unsorted
   order).
5. Print a summary table: file -> status (created / updated / routing-only / skipped-exists / error)
   -> key, then counts per status, then any index findings from step 4. If any paper errored, list
   those paths so they can be retried.
