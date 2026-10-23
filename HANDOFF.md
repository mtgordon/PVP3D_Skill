# Hand-over log (newest first)

New here? Read `notes/PROJECT_BRIEF.md` first.

Each entry: date, who, what ran (full run names) and what it showed, what is still running and where, what is next, and who
has the turn. Only the person who has the turn starts runs.

## 2026-10-23, Mark: the repo becomes the hand-over repo
- **Turn:** Mark (runs still going on Mark's machine).
- **Bases in this repo** (`claude_diag/runs/`): L403 (fine strap), L423b (canal at the source width), L431 (+ strap x1/3:
  first clean t 1.5 finish, old curve material), L439 (fibre curves + strap-side contacts; stalled t 0.64), L440g (+ fine
  distal curve + strap merged onto it), L444 (L440g line + 4-fibre AVW / PVW wall, source k), L445 (L444 + PM_conn ends on
  the AVW only), L446 (L444 + strap x1/3), L447 (L444 + stiffer wall matrix c1 0.0005, E 0.011), L448 (L446 + PM_conn on every AVW distal-row node).
- **Running on Mark's machine:** L444, L445, L447, L448 (L440g stopped at t 0.53: looked bad; L446 stopped for L448: it kept
  the old PM_conn layout). Results go into README section 8 and `run_log_notes.csv`.
- **Learned today:** a 2-fibre (0 / 90 deg) wall with a low bulk modulus blows up (L442c); use 4 in-plane fibres + the source
  k (`build_batch444.py`). The curves' fibre material is still 2-fibre: probably the same weakness (next to test).
- **Next:** judge L444 / L445; then the softer strap (x1/3, as L431) on top; 4-fibre curves.
- Details: `NEXT_SESSION_PROMPT_2026-10-23b.md` and the next one.
