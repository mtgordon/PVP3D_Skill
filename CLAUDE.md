# PVP3D pelvic-floor model in FEBio: instructions for Claude Code

The Abaqus PVP3D model (`PVP3DModel_job.inp`, the source and ground truth) is being converted to FEBio 4.13 and made to
converge. Two people (Mark and a colleague) take turns. Each has their own Claude Code, and this repo carries the work between
them. The person whose turn it is steers. **Ask them before starting new runs, and before any commit or push.**

## Step 1, every session
1. `git pull`, then read the top entry of `HANDOFF.md` (whose turn it is, what is running, what is next).
2. Use the project skill `abaqus-febio-fea-pipeline` (in `.claude/skills/`; Claude Code loads it from there). Read its
   `SKILL.md` and the reference files it points to as needed.
3. Read the newest `NEXT_SESSION_PROMPT_*.md` (the detailed hand-off), `claude_diag/README_2026-10-09.md` section 8 (the run
   journal, newest at the end), and `notes/MEMORY_NOTES.md` (the rules both sides follow, especially the run rules and the
   "NOT IN SOURCE" rule: say plainly when a change is ours, not the source's).

## Setup (once per machine)
- FEBio 4.13 (FEBioStudio). The tools expect `C:\Program Files\FEBioStudio\bin\febio4.exe`; set `FEBIO4` if it is elsewhere.
- Python 3.10 with numpy, scipy, matplotlib, openpyxl. The tools are run as `py -3.10 claude_diag/tools/<tool>.py`.
- Git Bash (the run scripts are bash and use PowerShell's `Get-CimInstance` to see febio4.exe, so Windows is assumed).
- Before cloning: `git config --global core.longpaths true` (some base run paths exceed Windows' 260 characters).
- Paths: the tools find everything relative to `claude_diag/tools/paths.py`, so the repo can live anywhere; a short path
  (e.g. `F:\PVP3D`) helps, because Windows' 260-character limit bites on the long run names.
- The skill's scripts are found at `$PVP3D_SKILL_SCRIPTS`, else `.claude/skills/abaqus-febio-fea-pipeline/scripts` here, else
  `~/.claude/skills/...`.

## Runs
- A model variant is built from a base run by a `claude_diag/tools/build_batch<N>.py SRC_RUN NEW_RUN [options]` script. It
  writes `claude_diag/runs/NEW_RUN/NEW_RUN.feb` plus `.feb.changes.txt` (the base and every change, NOT IN SOURCE marked).
  One change per variant.
- Start runs only with `bash claude_diag/tools/queue_threads.sh LIMIT 8 RUN...`, where LIMIT = this machine's PHYSICAL cores
  (Mark: 32). Runs never share cores. After a run starts: `bash claude_diag/tools/autostop_runs.sh 30 900 RUN...` (stops a
  run that makes no progress for 30 min).
- Watch them: `bash claude_diag/tools/wait_any_end.sh 15 RUN...` (in the background), then `py -3.10 claude_diag/tools/run_eta.py RUN...`
  and `bash claude_diag/tools/progress.sh RUN...`. Report a table with the full run names, t, and the ETAs to t 1 and t 1.5.
- Judge with `curve_cross.py` (curve self-crossing), `strap_cross.py`, `pm_cross.py`, `pop_q_bump.py` (POP-Q as Bump 1996
  defines it) at t 0.5, t 1 and the end of the hold (t 1.5).
- Keep run names short. Run numbers: Mark uses `L<number>`; the colleague uses `C<number>` (e.g. `C001_...`), so the two never
  clash. A base you build on is named in `.feb.changes.txt`; if it is the other person's run, it must be in this repo.
- `.xplt` / `.log` files are never committed (`.gitignore`). Results go into words: `py -3.10 claude_diag/tools/run_notes.py RUN
  --why ... --result ...` (writes `claude_diag/run_log_notes.csv`) and the README section 8 journal.
  `py -3.10 claude_diag/tools/run_log.py` rebuilds `PVP3D_run_log.xlsx` (~20 min) from the runs on THIS machine; the other
  person's rows come only from the notes CSV.

## End of a turn (before handing over)
1. Update README section 8, `run_notes`, `claude_diag/PAPER_ADDENDUM_since_2026-10-19.md` (paper decisions), and
   `notes/MEMORY_NOTES.md` if a new rule was set.
2. Write a new `NEXT_SESSION_PROMPT_<date>.md`, and a new top entry in `HANDOFF.md`.
3. `git add -f` the `.feb` + `.feb.changes.txt` of every run the next person needs as a base (not every run).
4. Show the person `git status`; commit and push only when they say so. Work on `main` in small commits; never force-push.
