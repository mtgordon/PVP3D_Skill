# Working notes (from Mark's Claude Code memory, copied 2026-10-23)

A tidied copy of the auto-memory the Claude sessions on this project kept, one section per note. Dates are when the rule
was set; later notes and `HANDOFF.md` win where they disagree. Left out: machine-specific notes (Box folder, disk space,
a Claude Code outage). Add new rules here (and in your own memory) so both sides see them.

## Rules and preferences (how to work)

### `feedback-avw-damping-with-rest`
*From 2026-10-22: AVW mass damping (against the late surge / fold-back) is OK only with a rest (hold) period at the end or evidence it does not change the final state: the answer must be quasi-static*

The user (2026-10-22 ~08:00): "I'm not opposed to adding damping to the avw as long as we either have a rest time at the end
or we don't think it's affecting the final results so we're actually getting a quasi-static solution."

**Why:** AVW mass damping C 2000 (L392 / L393) removed the anterior-wall surge / fold-back but the wall was still creeping at
t 1 (p95 ~0.1 mm per 0.001 t), so the t 1 state is held back (lag), not an equilibrium.

**How to apply:** every damped candidate gets a hold (time_steps 100 -> 150 or more, full load held past t 1) and is judged
at the end of the hold, with node speeds ~settled (node_speed.py p95 ~0.05 or less, falling); report the change from t 1 to the
hold end. Mass damping does not change a static equilibrium once velocities are ~0, so a settled hold is quasi-static.
Links: `feedback-jitter-mass-damping`, `project-curve-pressure-fix`.

### `feedback-cervix-stiffer`
*From 2026-10-09 the cervix may differ from the vaginal walls: start at 3x the wall stiffness, up to ~10x; the AVW and PVW stay one tissue*

The user, 2026-10-09: "I'm fine with you changing the cervix tissue now. It should start at 3x the stiffness of the vaginal
wall and can maybe go up to 10x."

Earlier (2026-10-07) the rule was "keep the entire vaginal wall (and cervix) the same tissue" (one Yeoh, as the PVP3D source).
Now only the cervix (Vagina_Cervix, _PickedSet346) may be scaled, as a whole tissue (c1, c2, k together:
build_batch145.py `cervix=S`, bounds x3-x10 per the user); the AVW and PVW stay the source's one wall tissue.

**Why:** the user's call (likely: a real cervix is much stiffer than the vaginal wall; it may hold the apex so less CL is needed).

**How to apply:** new fits may start from cervix x3; say it is a whole-tissue change of the cervix only; NOT FAITHFUL (the source
gives the cervix the wall's Marlow data). Related: `feedback-lips-minimal-influence`, `project-default-bases`.

2026-10-09 (later): the cervix x3 of the wall is now the NORMAL setting for new models on both lines (`project-range-of-motion-goal`).

### `feedback-curve-first-no-popq`
*From 2026-10-20 afternoon: don't report POP-Q / range measurements until the curve interaction (curve self-collapse, strap-curve) is solved*

The user (2026-10-20 ~16:40): "don't worry about any of the measurements until we figure out the curve interaction. I think that will greatly change the measurement results anyway."

**Why:** the curve collapse / self-penetration and the strap-curve interaction change the wall motion, so POP-Q (Aa/Ba/Bp/C) and range comparisons are not meaningful yet.

**How to apply:** while working the curve problem, judge runs on the curve itself (self-penetration depth and timing, tip-edge kinks, strap-curve contact behaviour, pace) and leave POP-Q / range out of reports and run decisions until the user says the curve interaction is settled. Related: `project-curve-selfpen-rules`, `project-range-of-motion-goal`.

### `feedback-equal-end-attachments`
*From 2026-10-06 every new tube model has all curve-end attachments (PM_avw_bottom L/R + curve_end_PM_conn) at the same stiffness, PM_conn_mat*

From 2026-10-06 (the user: "all models going forward should have the same stiffness on the end attachments"): in every new
tube model the springs at the curves' distal ends -- the source's PM_avw_bottom_left/right_conn (2 a side; source law PM-LA
x1e-5, i.e. off) and the added curve_end_PM_conn (3 a side) -- all use PM_conn_mat's force points (tools/build_batch187.py
does it). Mark NOT IN SOURCE / NOT FAITHFUL.

**Why:** the source's bottom springs were ~1e5 x softer than the added end springs; the user wants the end attachments uniform.
L187 (the change alone on L169) did not fix the curve fold (19-24 nodes inside vs L169's 2-26).

**Scope (the user, 2026-10-06):** TUBE models only -- not the seam-band (non-tube) models. "Just the ones that attach to (or
are fixed on) the PM arch": PM_conn (26), PM-LA-x%stiff (1), curve_end_PM_conn (6), PM_avw_bottom L/R (2 + 2), all anchored on
the PM's outer arc, all at PM_conn_mat's curve (PM-LA's is the same curve). NOT the PM_PeB springs (tied to the PM shell, not
the arch), AVW-Para, sphincter, PeB-constrin or Parcus.

**Position (the user, 2026-10-09):** when CL / USL / AVW-Para / Parcus wall ends are moved onto the curves' outer surface
(tools/build_batch242.py attach), the PM-end springs (PM_conn, curve_end_PM_conn, PM_avw_bottom) stay where they are.

**How to apply:** start new tube builds from L187-style end springs (or apply build_batch187's change); L187 / L190 already
satisfy the rule.
Related: `feedback-pm-conn-inner-arc`, `feedback-three-base-versions`.

### `feedback-extra-runs-slow-first`
*When the user asks about running other things alongside a run they are waiting on, remind them it slows that run (shared physical cores / memory bandwidth)*

When the user asks whether to run other things while waiting on a run, remind them up front that extra runs will slow the
run they are waiting for, and say roughly by how much.

**Why:** 2026-10-15: filling the machine from 40 to 64 threads cut L286_around2_MAX from ~6 to ~2.5 steps/min and pushed its
t 0.7 estimate from ~11:30 to ~13:00-14:00; the user wants that trade-off stated before deciding. The machine is 2 x Xeon
Silver 4216 = 32 physical cores / 64 logical (hyperthreads): 64 threads means two threads per core, plus shared memory
bandwidth for FEBio's sparse solver, so each run slows even though none is oversubscribed.

**How to apply:** with any "run more while we wait?" recommendation, give the expected slowdown of the priority run (e.g.
"<= 32 threads keeps one thread per core; beyond that the waiting run slows ~2x"), then let them choose. Related:
`feedback-febio-run-workflow`.

### `feedback-febio-docs-and-notes`
*Consult the FEBio manuals and forum, not only run results; keep what is learned in the Claude Doc of notes*

When investigating FEBio behaviour (shells, contact, convergence, how the model compares with Abaqus), also read the FEBio
manuals (https://febio.org/knowledgebase/, e.g. help.febio.org user and theory manuals) and the forum
(https://forums.febio.org/), not only the run results. Keep what is learned in a document: the Claude Doc
"PVP3D FEBio notes: shells, contact and the LA" (https://claude.ai/code/artifact/18b2a765-6850-4206-9660-fa08fbb42fd6),
mainly for what the manuals and forum say, plus other findings as useful. Cite the pages opened.

**Why:** the user, 2026-09-27: "don't just look through results of runs, also look at the manuals ... and the forum" and
"keep notes of what you discover in a document" ("that last comment about the document was for things learned from the
manuals and the forums. Although put other things in the document as you see fit").

**How to apply:** before proposing a solver/contact/shell setting, check the manual's parameter table and search the forum;
add the finding (with its link) to the Doc. Related: `feedback-febio-run-workflow`, `project-pvp3d-handoff`.

### `feedback-febio-run-workflow`
*How the user wants FEBio variant runs proposed, started, stopped and reported in the PVP3D project*

Ask before starting any FEBio run and say what each run is for; build and check variants freely. One change per
variant, reported as a comparison table (variant, change, t reached, failed attempts, failure location). Runs may be
stopped without asking once they make no new converged step for 30 min (claude_diag/tools/autostop_runs.sh). The user
pauses to understand results before choosing: give the findings, the choices and a recommendation, then ask. Screen on
the fast base (arcus chain mass x10); confirm like-for-like later. Runs stay on F:; never overwrite the user's files.

**Why:** the user said "ask me before starting", chose "stop after 30 min stalled", and repeatedly paused runs to
understand what had been learned before picking the next ones (2026-09-24).

**How to apply:** at every batch boundary summarize results and options with a recommendation, then use a question
to get approval; label anything not in the Abaqus source as NOT IN SOURCE. See `project-pvp3d-handoff`.

File names, 2026-09-26: in every results table, summary and review (chat, SESSION_SUMMARY, the README, the next
NEXT_SESSION_PROMPT) name each run by its full file name, never shortened with "…" or a nickname ("the recipe", "L53",
"the hold run"), and say where it is: claude_diag/runs/<run>/<run>.feb, with the results in the .xplt beside it. In
chat, make the names clickable links. The user asked for this after a table that shortened the run names ("add the
filenames so I know what files they are"). Carry this rule into the "How I like to work" part of each new
NEXT_SESSION_PROMPT.

Waiting, 2026-09-27: when working on my own, cap every background wait at ~30 min (an until-loop that also exits on a
deadline), then re-check every run. A waiter on "run X ends" sat 6.5 h on one crawling run while six crawling runs went
unchecked and the machine did no useful work. The autostop script only stops runs that make NO progress; a crawl in tiny
steps is never stopped, so check the rate (t gained per 30 min) myself.

Fixes that change the answer, 2026-09-28: the user rejected a seam tie that made runs converge but moved the anterior wall
14 mm ("I don't like that the ties are affecting the results that much"). For every convergence fix, show how far it moves
the results (and their trends across load cases) next to its convergence gain, and prefer fixes that leave the physics
alone; label physics changes NOT IN SOURCE.

Exception, 2026-09-25: when the user says to keep working while away ("I'll check back in about 8 hours"), start runs
without asking. Keep to one change per variant, stop runs that have answered their question (log each stop in
runs/_autostop.txt), and leave their choices (lofts vs connectors, Box copies they said to ask about) for their return.
2026-10-05 (unattended night): do not wind runs down or end the session's runs at a set time ("8-10 h" is not a deadline);
keep going until the user checks back, then they decide what to stop. Crawl / answered-question stops still apply.
2026-10-06 morning: the user allowed healthy APICAL (CL / USL) scales above x5 "if needed" (to hold Ba); other levers keep
x0.2-x5 (build_batch145.OVER5).
2026-10-06: never oversubscribe the machine: the threads of all running febio4.exe together <= 64 (logical cores); e.g. 6 runs
x 10 or 8 x 8, counting runs already going. On 10-06 night 9 x 8-10 (~78 threads) slowed the tube runs ~2x (L172 8 h vs
L169 4.5 h). queue_runs.sh now refuses MAXPROC x THREADS > 64.

2026-10-17: the user: "if I don't respond for 15 minutes go ahead and keep doing what you think is necessary." -> when a question to the
user (e.g. start a run) gets no answer for ~15 min, proceed with the recommended option and keep working; say so in the next message.

### `feedback-fibre-walls`
*from 2026-10-23: future versions put the AVW and PVW (not only the curves) in the tension-only fibre material so the wall can buckle/fold*

From 2026-10-23 the user: the vaginal wall should not be stiff in compression; it folds. First the curves were switched to a fibre material (tools/build_batch434.py, Curve_fibre):
- a soft Mooney-Rivlin matrix plus two tension-only fiber-pow-linear-uncoupled fibres, one along the canal and one around;
- fitted in a one-hex strip test (shell_tests/curve_fibre/gen.py) to match the old material in tension;
- 5-14x softer in compression.

Then the user decided: "include the avw and pvw in the fiber modification as well. Let's do that for future versions."

**Why:** today's Yeoh/Marlow-type laws are stiffer in compression than in tension (the curve material: 6.9 vs 3.0 kPa at +/-30 %). So the tissue resists like a block and is forced through itself instead of buckling.

**How to apply:**
- Every new model gets the curves AND the AVW / PVW in fibre form, each fitted to its own source law in tension by the strip test.
- AVW and PVW stay one tissue (`feedback-cervix-stiffer`).
- WALLS (10-23 L442c failure): two fibres at 0/90 leave no in-plane shear stiffness (45 deg pull 6x softer) and k scaled down with the matrix lost volume control: the AVW stretched diagonally, V/V0 0.5-1.6 by t 0.25, then exploded. Use FOUR in-plane fibres at 0/45/90/135 deg (build_batch444: E 0.014, beta 2.5, lam0 1.04, c1 0.000208) and keep the SOURCE k (0.35): k does not stiffen uniaxial compression, it only stops volume drift. Strip-test pulls at 0 AND 22.5/45 deg, not just along a fibre. A spherical fibre distribution (cfd) is stiffer in compression than tension (lateral fibres engage) - not usable.
- Curves (L434): scaling k with the matrix was used there; recheck the curves with the same off-axis test.
- Keep the curve self-contact penalty at 5: x3 and x4 crawl (L417, L437b, L438).
- Related: `project-canal-width-source`, `project-curve-selfpen-rules`.

### `feedback-finer-strap-mesh`
*from 2026-10-22: every future fan strap (PM_avw_bottom sheet solid) gets the 4x finer in-plane mesh (build_batch404)*

From 2026-10-22 (the user, after L403): all future straps get the finer mesh, the 4x in-plane refinement of build_batch404.py (strap edges ~0.7 mm instead of ~2.9 mm, children keep the fibre directions, mid-edge ties). Apply it to every new strap model, on both lines.

**Why:** in L403 the finer strap halved how far the strap pushed into the curve at the tip (0.5-0.9 vs 1.1-1.7 mm). The user said it "kept it together better".

**How to apply:** when you build a strap variant, run build_batch404 on it (or start from an sref base) unless the user says not to. Say so in the run's changes.txt. The user also asked for the curve's corner material at the back to be meshed as finely as the front. Related: `project-curve-selfpen-rules`, `project-pvp3d-handoff`.

### `feedback-jitter-mass-damping`
*When a run jitters or crawls on a few chattering nodes, try local mass/damping on the jittering part (its own variant), and check the lag*

The user (2026-10-18): whenever a run jitters (nodes chattering, steps crawling), try adding some mass / damping, and keep
this in mind as a standard thing to try.

**Why:** it worked before: L309 (near-zero lofts + loft density x100 + mass damping C 200) ran at normal pace where L305
(near-zero lofts alone) crawled; distal curve damping delayed the curve snap (L289, L314).

**How to apply:** when a crawl is traced to a few fast-moving nodes (velocity between converged states), make a variant with
density x10-100 and `mass damping` (body_load, elem_set = that part only) on the JITTERING PART ONLY, one change, next to
whatever else is tried. Then check the lag: damping holds the motion back (L314 lost 16 mm of Bp at t 1, L309 ended at
Ba +9 / Bp +8), so compare the range at t 1 and, for a candidate, hold the loads past t 1 to see it settle (L292-type hold).
Large or global damping crawls (L299 x10 on the chains line). Related: `feedback-febio-run-workflow`.

### `feedback-lips-minimal-influence`
*Seam lips (tools/seam_extend.py) are a convergence aid only: never tune their stiffness/size to steer Ba/Bp; they should affect the result as little as possible*

**Name (the user, 2026-10-09): call the seam extensions "ledges" from now on** (and from later on 2026-10-09 the lips line is the "ledge line"; the
tube's wrap parts are the "curves").

The user, 2026-10-07 evening: "I don't want the lip stiffness to be adjusted to try and control the movement. They should
affect the run as little as possible as they are not in the original design." (Said after I proposed stiff AVW-side lips,
L207_ext5_avwlipx3/x10 — built, never run.)

**Why:** the lips (unloaded extrusions of the walls past the contact-only seam) are NOT IN SOURCE; their only job is to stop
seam edge nodes flipping off the other wall's last facet.

Measured 2026-10-07 (P1, 5 mm lips): lips at the wall material hold the canal up 3-4 mm; x0.1 vs x0.01 differ by
<= 0.6 mm (x0.01 needed 28 retries; the late movers in these runs are PM_Plane nodes, not lips: last_moves.py prints plot
indices; tools/lip_motion.py shows lip steps ~0.004 mm). **Use lips at c1, c2, k x0.1 of the wall (tools/build_batch208.py).**

**How to apply:** keep lips on the wall's own (or a softer) material and short; check that the answer barely depends on them
(lip length / stiffness sensitivity) rather than using them as a fit lever. Fit P1 / P2 with the paper's tissues (supports,
walls) instead. Related: `feedback-say-whole-vs-local`, `project-default-bases`.

### `feedback-no-shared-cores`
*From 2026-10-19 the user: start new runs only if they won't share threads/cores with running ones (32 physical cores) - queue with limit 32*

The user (2026-10-19): "run them only if they won't share threads and slow down the other runs ... or share cores".
The machine has 32 physical cores / 64 logical; 4 runs x 8 threads already fill the physical cores.

**Why:** extra runs on hyperthreads slow the runs the user is waiting on (`feedback-extra-runs-slow-first`).
**How to apply:** queue new runs with `tools/queue_threads.sh 32 8 RUN ...` (limit 32, not 64) so each starts only when a
running one ends; CPU-heavy side jobs (run_log.py rebuild) wait for a free core too. This tightens the older "<= 64 threads"
rule (`feedback-febio-run-workflow`) unless the user says otherwise.

### `feedback-no-spring-compression`
*From 2026-10-10 do not give the connector springs (e.g. AVW-Para) compression to hold the healthy AVW; use lofts + springs (chains in sheets) instead*

The user (2026-10-10): "for ledge, don't add compression on the springs yet as I don't think that's physiologically accurate.
You can try it with the loft + springs to see if that changes anything." The source's connectors are tension-only; keep them so.

**Why:** compression in a ligament spring is not physiological; the user prefers the lofts + springs representation (sheets carry
in-plane load like tissue).

**How to apply:** to hold the healthy AVW, use the springs + lofts version (tools/build_batch262.py, sheets as a fraction of
the family's fibre E1; 30 % worked: L271_ledge_cx3_chains_sheet30_*) or the cervix 3-10x, not build_batch217 compression.
Related: `project-range-of-motion-goal`, `feedback-three-base-versions`.

### `feedback-pm-conn-inner-arc`
*PVP3D FEBio - every new model anchors the 27 PM_conn-family springs (AVW to PM) on the PM_Plane's OUTER arc (the clamped BC), not the inner arc*

**From 2026-10-02 on (supersedes the inner-arc rule of 2026-10-01):** every new PVP3D model has the PM as a structure
(build_batch113 step 2) and anchors the 27 PM_conn-family springs (the "top" springs from the AVW to the PM arch: the 26
PM_conn + the middle one, set `PM-LA-x%stiff_mat`; source law `PM-LA-x%stiff`) on the PM_Plane's **outer arc**, the clamped
edge (boundary condition BC-PM-outer-arc): each anchor moved along its spring to where the spring crosses the outer arc
(springs shorten from 14.04 to 11.6-13.9 mm), the anchor fixed there (its PM_anchor_ties entry dropped), the law's
elongation axis scaled by the mean length change (as build_batch129b). Builder: claude_diag/tools/build_batch149.py.
The near-zero families (PM-PeB, PM-AVW-bottom, PeB-constrain: 0.001 % of PM-LA in the source) stay as in the source (the
user's choice, 2026-10-01).

**Why:** 2026-10-02 the user: "I actually want the top springs that attach the AVW to the PM arch to be moved to the top of
the arch which is a boundary condition. Make a note of that for future builds (both here and going forward). That could
maybe help with the speed as well." (The inner-arc anchoring doubled run time and pulled the free PM edge 5 mm out of plane.)
Earlier: "Make all future (don't redo these) models attach the PM_Plane"; the middle spring is included ("keep it fixed").

**How to apply:** include it in every new builder (use build_batch149's function, not build_batch129). Carry it into each
NEXT_SESSION_PROMPT. Related: `feedback-three-base-versions`.

### `feedback-popq-bump-definitions`
*Measure POP-Q points (Aa, Ba, Ap, Bp, C) as Bump et al. 1996 define them, with pop_q_bump.py; the user checks definitions against primary sources*

Report POP-Q points as Bump et al. 1996 (Am J Obstet Gynecol 175:10-17) define them, with
claude_diag/tools/pop_q_bump.py (summary_table.py uses the same): Aa / Ap are points fixed 3 cm in from the urethral meatus /
hymen, tracked; Ba / Bp are the most dependent point of the upper wall (fornix down to Aa / Ap), found anew at each state;
C is the most distal edge of the cervix; all signed from the plane of the hymen, negative above. Report Aa, Ba, Ap, Bp, C.

**Why:** on 2026-09-28 the user asked how Ba / Bp were measured, supplied the Bump paper, and was right that Ba / Bp are
not tracked points: the old tool (pop_q.py) reported a tracked node 3 cm above the hymen (Bump's Aa / Ap), following the
source paper's one-line paraphrase. The numbers happened to agree (within 0.1 mm), but the definition was wrong.

**How to apply:** when a result is compared with a clinical or published measure, read its primary definition and implement
it as defined; state plainly what was measured before, and correct it without softening. Related: `project-pvp3d-handoff`,
`feedback-febio-run-workflow`.

### `feedback-run-log-spreadsheet`
*Keep PVP3D_run_log.xlsx (one row per FEBio run, every run folder) updated for every run; layout approved 2026-09-26*

The user wants a spreadsheet of every run, kept updated as runs are done: F:\Developing PVP3D Model for FEBio\PVP3D_run_log.xlsx,
one row per run folder (all 277 on 2026-09-26, in build order). It is rebuilt by claude_diag/tools/run_log.py from the run
folders. The hand-written parts live in claude_diag/run_log_notes.csv, with columns run, why, result, label, focus and source.
- When a run is set up, add its row with "why" and "focus" (rotation / pvw / la, or blank for convergence-only runs).
- When it has answered its question, fill "result" (and source = written <date>), then rebuild.
- If the xlsx is open in Excel, the script writes a dated copy instead.
- 265 rows were drafted after the fact from the READMEs, SESSION_SUMMARY and the builders; their source says so.

Layout the user asked for (2026-09-26 evening):
- one column per property of every material (tissues, lofts, connector spring materials, the arcus chain);
- "connectors" or "springs", never "cables";
- rotation and other measurements only for runs where they were the focus;
- what was changed highlighted; yellow = changed in this run, light orange = changed earlier (differs from the line's
  reference, only for its descendants), grey = not in the model;
- row order does not matter to the user.
The user approved this layout with "yes, do it for all the runs".

Base models (the user's request, 2026-09-28): a "Base models" sheet lists every adopted base per line, in order, with the
decision and every property that differs from the previous base (read from the two .feb files, contacts setting by
setting); a green "Base model" column marks those rows on the Runs sheet. The list lives in claude_diag/tools/run_log.py
(BASES): add a row whenever the user adopts a new base.

**Why:** the user wants to look at one sheet and know what every run was, why it was run, what changed and what came of
it, and when a new base was adopted and how it differs from the one before.

**How to apply:** update the notes and rebuild after every batch, as part of reporting results. Name runs by their full file
names, as in `feedback-febio-run-workflow`.

### `feedback-say-whole-vs-local`
*When proposing a stiffness change, say plainly whether it changes a whole tissue or only a local region (the seam band); the user prefers tightening support tissues over stiffening whole walls*

When proposing or reporting a material / stiffness change, state explicitly whether it applies to a **whole tissue**
(every element of the AVW, cervix, PVW) or only a **local region** (e.g. the seam band within 3 mm of the seam). On
2026-10-05 the user approved "AVW + cervix fully x3 / x5" thinking it was a wider seam band and was surprised it stiffened the
entire tissue ("I didn't realize the stiffness of the entire tissue was being changed on those ones"). Their direction
instead: keep the seam-only band and **tighten the support (connective) tissues** so things move less.

**Why:** the user reasons about the anatomy; a whole-wall stiffness change is a much larger departure from the source than a
local band, and they want to choose it knowingly.

**How to apply:** in every proposal table, a column or phrase "whole tissue" vs "seam only / region"; never shorten "whole
AVW" to "AVW x3" in a question; prefer support-tissue levers (apical, AVW-Para, PM_conn, Parcus) for fits. Related:
`feedback-febio-run-workflow`, `feedback-sphincter-side-softer`.

### `feedback-sphincter-side-softer`
*Every new model: LA_sphincter_side_conn_mat at 25 % of the source (scale x0.25; was x0.5 earlier on 2026-10-05) so the perineal body rotates and Bp can increase*

Every new PVP3D model: LA_sphincter_side_conn_mat at **25 % of the source** (its `<scale>` x0.25), the 16
LA_sphincter_side_conn springs between the levator (LA_PCMPRM) and the perineal body (_PickedSet66). The user, 2026-10-05:
first "in all future models the LA_sphincter_side_conn_mat need to be less stiff. Start with 50% less stiff", then the same
evening "decrease the connector strength even further. It's still not rotating as much as I would expect. Maybe decrease
it to 25% of the original".

**Why:** the user expects these springs hold the perineal body / distal PVW too firmly: it does not rotate as much as they
expect, and softening should let Bp increase (all fits so far have Bp short where Ba is right).

**How to apply:** apply in every new build with the other standing rules (`feedback-pm-conn-inner-arc` PM_conn on the outer
arc, `feedback-three-base-versions`); mark it NOT IN SOURCE; LA_sphincter_post_conn unchanged unless asked. The x0.5
runs of 2026-10-05 (batch 174) are the earlier step; existing runs are not rebuilt without asking.

### `feedback-three-base-versions`
*PVP3D FEBio - every new base model must be made in three versions (springs, springs+lofts chains, lofts)*

When a new base model is adopted for the PVP3D FEBio project, keep three versions of it:
1. springs only (every connector a spring; the springs line),
2. springs + lofts (each connector a chain of springs stitched into a weak sheet, as batch 111 `L111_chains_all_m05`;
   claude_diag/tools/build_batch132.py `chains`),
3. lofts only (the lofts line, from `L87_lofts_newline_vwyeoh_rhoi0`; build_batch132.py `lofts`).

**Why:** the user (2026-10-01) wants each base change carried across all three connector representations to compare them.

**How to apply:** whenever a base is adopted (added to `BASES` in claude_diag/tools/run_log.py), build and run the other two
versions with the same change. Related: `feedback-pm-conn-inner-arc`, `project-pvp3d-handoff`.

## Project state and decisions

### `project-canal-width-source`
*2026-10-23: narrowing walls + curves + cervix to the SOURCE width (build_batch422) DELAYED the hold collapse of the curves to t ~1.45 (L423b), not removed; less prolapse*

On 2026-10-22 the user asked for the AVW, PVW and curves (and the cervix) to be the same width as the source's AVW/PVW/cervix (PVP3DModel_job.inp). Our curves had made the canal 8-28 % wider than the source. tools/build_batch422.py narrows the canal in x by level, shrinking only, and spreads the 26 PM_conn springs in order along the whole distal edge, curve corners included, so none cross.

Result on L403's smaller strap, run L423b (MAXl_ant1):
- **Hold:** CORRECTED 10-23 07:20 (the user): the collapse is only delayed, from t ~1.03 to ~1.45-1.48. Settled to t 1.4 with a steady ~1.1 mm distal-curve self-crossing (rest z -31, just above the strap / opening), then at constant load it jumps to 2.4 mm and the walls follow.
- **Softer strap (L431, strap material x1/3) on the same model: NO late collapse.** It ran to t 1.5, settled, with the curve self-crossing steady at ~1.24 mm. Ba/Bp -0.0/-1.4 at t 1.5.
- **Prolapse:** much less (Ba/Bp -2.7/-3.9 vs +7.7/+9.4).
- **Healthy (L423h):** -18.8/-21.5.

L423 itself was slow because build_batch288 --only had added a second spring_patch_ties block; that is fixed.

**Why:** a wider-than-source canal gets dragged through the hiatus and folds; the strap cannot narrow it again (the user's reasoning).

**How to apply:** not settled as the base. On 10-23 the user doubted that narrowing addresses the problem: the failure is still at the narrow opening, so the ratio of opening to canal width pushed through may be what matters. A softer strap may only move the failure further out. Options raised: an in-between width, softer curves / strap / walls, less curve pressure. Ask before choosing the base. Contact tweaks alone did not stop the fold: curve_strap penalty, all strap faces, a weak canal-PM contact, curve mesh refinement. Related: `project-curve-selfpen-rules`, `feedback-finer-strap-mesh`, `project-pvp3d-handoff`.

### `project-curve-pressure-fix`
*From 2026-10-21: OUR curve pressure (not in the source: it has no curves) stops the distal-curve collapse (L359) by jacking the apex up ~22 mm, then the lift releases at t ~0.9 (the slam) and the run crashes*

**The PVP3D source (PVP3DModel_job.inp) has NO curves and NO curve pressure** (the user, 10-21: "The original abaqus
file didn't have curves"). The curves (canal_wrap_hex_*, build_batch135, 10-01), their pressure Load-curves (the wall
pressure 0.014 MPa carried onto their outer skin, L140; removed on the tube line L254 / L258) and Curve_soft (Vagina_AVW
x0.6) are all OURS. Never call them "as in source"; the 10-20 / 10-21 night records did and were corrected 10-21.

2026-10-21: putting the curve pressure back (L359 = L352 + pressure, tools/build_batch364.py) stopped the distal-curve
hairpin (<= 0.46 mm to t 0.88). Mechanism: the pressure squeezes the soft, nearly incompressible curves, they lengthen and
jack the cervix up ~22 mm above rest (L352: 4-13 mm down); at t 0.88-0.93 the lift is released (cervix drops 42 mm; the
user: "the vaginal wall flies up and then slams back down right before it crashes") and the distal PVW / curve fly off.
Milder MAX-lite (L368 la60, L369 P-arcus x0.075) lift the same way. The user wants runs to reach t = 1.

Rejected 10-20 night: more self-contact faces (L361), pre-bend (L362 / L365), sphincter side 0.1 (L360 / L363), strap <->
curve contact (L357), distal curve mass x100 (L366: drags the PeB).

**Why:** the "fix" is our own load acting as a jack, not fidelity; the slam is its price.
**How to apply:** treat curve pressure, curve material and curve geometry as modelling choices to justify; next proposed
test L372 = L359 with Curve_soft x0.6 -> x1 of the wall (less squeeze / lift). Links: `project-max-lite`,
`feedback-curve-first-no-popq`, `project-curve-selfpen-rules`, `reference-normal-generic-model`.

Update 2026-10-21 evening: the pressurised / stiffened curves act as RAILS: the canal moves out as a rigid tube (no mid-canal
bulge), so the length goes into the arch, not the prolapse (L352 bulges 33-37 mm). Best so far: L381 (full candidate + half
curve pressure) t 1 clean but Ba / Bp -9.8 / -9.5. Hand-off NEXT_SESSION_PROMPT_2026-10-21b.md.

Update 2026-10-22 morning (night stretch, NEXT_SESSION_PROMPT_2026-10-22.md): the full candidate is capped (no pressure:
Ba +1 at t 1; quarter pressure L384: Ba -6.1, cleanest). The springs-line L353 (stiff strap + strap-curve contact) + quarter
/ eighth curve pressure reaches the targets with a clean curve (L389 / L391, Ba +10-12 at t 0.88-0.92) and then the AVW
SURGES and folds back. AVW mass damping C 2000 (numerical aid) removes the surge but lags: L393 (1/8 + C 2000) t 1 Ba / Bp
+8.2 / +9.3, curve <= 0.69, no arch, mid-canal bulge (best). Open: its hold L394, C 600 runs L395 / L396, whether the user
accepts the damping.

### `project-curve-selfpen-rules`
*From 2026-10-16 night: curve self-penetration \"works\" = <= 0.3 mm (MAX); allowed levers patches/damping, curve-end mesh, contact - NOT curve material; split effort with loft contacts; combine only after each alone*

The user's rules (2026-10-16 ~20:35, before a 10-12 h unattended stretch; "keep trying things"):
- 2026-10-17 morning, after the night found nothing <= 0.3 mm: the user: "If we have to go up to 1mm allowable penetration I guess
  that's fine" -> up to ~1 mm accepted (L288c-type, peaks 0.9 mm, counts; a collapse past 1 mm does not).
- (10-16 night, superseded) "Works" for the tube curves means NO self-penetration (the user's word; I had been calling it "the fold"), measured by
  claude_diag/tools/curve_selfpen.py: deepest curve node inside a curve hex <= ~0.3 mm counts as none. L288c (peaks 0.9 mm)
  and L294 chains (1.2 mm) do NOT count. Say "self-penetration", not "the fold".
- Effort split equally: (a) a MAX version (chains + lofts or springs) without self-penetration and with the range met,
  (b) the loft-contacts plan (`project-loft-contacts-plan`).
- Allowed levers against self-penetration: spring patches + distal damping + held loads (build_batch288 / 289), curve-end
  mesh changes near the PeB, contact on the curves (self-contact variants, AVW-Para <-> walls). NOT the curve-end material
  (no stiffer tips).
- Combinations (e.g. patches + AVW-Para <-> PM + AVW-Para <-> walls) are allowed once each change has run on its own.

**Why:** the user saw the curves fold through themselves in L283 MAX and wants them not to penetrate, without trading range
for a stiffer tip; one change at a time shows which change causes trouble.

**How to apply:** report every MAX run's deepest self-penetration against 0.3 mm alongside its range; keep one change per
variant until combining. Related: `feedback-febio-run-workflow`, `project-range-of-motion-goal`.

**2026-10-22 correction: measure with tools/curve_cross.py, not curve_selfpen.py.** curve_selfpen counts only surface nodes
ending inside another curve hex and reports their distance to the nearest surface, which for the 1-2 mm curve wall cannot
exceed ~half the wall, and misses nodes that went right through (the user saw L399c's curve through itself while it read
<= 0.93 mm). curve_cross counts surface edges piercing curve faces (with the depth past the face): every run so far crosses
> 1 mm late (L350b ~3 mm distal; pressured runs 1-2 mm distal + ~1-1.4 mm at the mid-curve bend z ~-8..+2). The 10-20..10-22
"curve <= 1 mm" verdicts are withdrawn.

### `project-default-bases`
*Default models: band line (L199_band75av_0pv_pmd_sph025 + apical x6.5, posterior x0.35), tube line (L202_curvepress_apexclosed_P1), ledge line (called the lips line until 2026-10-09) (L228_ext5sl01_pmc_pmd_sph025 + CL x35, Parcus x0.59); build new variants from these*

The user, 2026-10-07: "Set the current tube line and band lines as defaults for their line."

- **Band line default:** base `L199_band75av_0pv_pmd_sph025` (L149_seamspr_pm_outer, seam springs out = contact-only seam,
  7.5 mm stiff seam band c1 x3 on the AVW + cervix only (59 % of the AVW), PVW band-free, wall + PM damping C 200, sphincter
  side x0.25) with healthy supports apical x6.5, posterior x0.35 via `build_batch145.py ... apical=6.5 parcus=0.35 --base
  L199_band75av_0pv_pmd_sph025`. Runs L201_b75a0p_ap65parcus035_{healthy,P1,P2}: -4.7/-16.7, -0.8/+3.1, +7.5/+10.2.
- **Tube line default:** `L202_curvepress_apexclosed_P1` (a P1 file; healthy / P2 not built yet): L169 tube + equal PM-arch
  springs (`feedback-equal-end-attachments`) + 0.014 MPa on the curves' outer faces + curve damping C 1000 + the apex
  closed (tools/apex_close.py). P1 -6.5/-10.5/-39.1.
- **Ledge line default (the user, 2026-10-09, then called the lips line; renamed "the ledge line" later that day: "make CL x35 Parcus x0.59 the default for the lips line"):** base
  `L228_ext5sl01_pmc_pmd_sph025` (band-free L199 + 5 mm unloaded seam lips at x0.1 of the wall + the PM shell remeshed,
  tools/pm_coarsen.py) with `build_batch145.py ... cl=35 parcus=0.59 --base L228_ext5sl01_pmc_pmd_sph025` (USL, para at the
  source). Runs L230_pmc_cl35parcus059_src_{healthy,P1,P2}: -3.5/-16.4, -0.5/+4.9, +3.3/+8.2 (~45 min each). Both Bp within
  ~0.2 mm of their edges. The remeshed PM lowers Bp slightly (P1 ~0.4, P2 ~1.2 mm): never mix old-mesh and remeshed fits.
- **Band line PAUSED (the user, 2026-10-09):** "let these run, but then let's pause on those as I think the lips is a better
  version". The running band-remesh checks (L232_b75a0p_pmc_ap65parcus035_*) finish and get recorded; no new band runs unless
  the user asks. The ledge line is the main line.
- All three are in run_log.py BASES (Base models sheet).
- The user's decisions, 2026-10-07 afternoon: **accept P2's Ba miss** (+7.5; no lever); **keep the 7.5 mm band** (not 9 mm);
  tube: if L203 (elastic curve-PVW join) opens the Bp - Ba gap, **duplicate the 200 shared curve-PVW nodes for a fully
  elastic join, then build the tube default's healthy and P2**.

**Why:** the night of 10-06/07 fit; the 7.5 mm band set meets healthy, P1 and the P2 Bp; the closed apex stopped the tube's
apical inversion.

**How to apply:** start new band, tube or ledge-line variants from these (one change each), not from L175 / L183 / L169 / L192.
Related: `project-pvp3d-handoff`, `feedback-say-whole-vs-local`.

- **Tube main from 2026-10-10 (the user: "adopt your version as the tube main"):** `L258_tubecand_P1` (cases via build_batch244: L258_tubecand_healthy / _MAX) = L202 + remeshed PM + cervix x3 + new apex cap + curves cervix band/cap + springs on the curves + curves x0.6 + sphincter x0.25 + NO curve pressure; replaces L202_curvepress_apexclosed_P1. Curve-pressure variants (full / distal half) L261_tubemain_* run next.

### `project-febio-studio-memory-commit`
*FEBio runs die mid-step when Windows commit memory runs out; a long-open FEBio Studio held 53 GB of it (2026-09-26)*

On 2026-09-26 two hold runs died mid-step with no error in the log and no Windows crash event:
L61_springs_pc10_cu30_antside0_hold at 14:02 and L63_lofts_pc10_cu30_avwparaconn_antside0_hold at 14:56. Their exit code
was 127, the same code a stopped run gets. Windows had 1.7 GB of committed memory left out of about 100 GB, while 22 GB of
RAM was free. The user's FEBio Studio, open since the day before and used to view result files, held 53 GB of commit with
a 1.7 GB working set. The page file is on C:, which is nearly full, so the commit limit cannot grow. Running out of commit
is the most likely cause; it was inferred, not proven. The user closed their files and FEBio Studio dropped to 11.5 GB.

**Why:** a solver that cannot commit memory dies without a message. This looks like a numerical crash, but it isn't one.

**How to apply:** before starting runs, and whenever a run "exited without a termination line", check the commit
(Win32_PerfFormattedData_PerfOS_Memory: CommittedBytes vs CommitLimit) and the top processes by PrivateMemorySize64. If
FEBio Studio holds tens of GB, ask the user to close their files or restart it; never close it yourself.
Both queries work here (checked 2026-09-27): the PerfOS_Memory class and Win32_OperatingSystem FreeVirtualMemory /
TotalVirtualMemorySize (KB). NEXT_SESSION_PROMPT_2026-09-27.md wrongly says the PerfOS class does not exist. Related:
`feedback-disk-space`, `feedback-febio-run-workflow`.

### `project-loft-contacts-plan`
*From 2026-10-16: after the tube springs+lofts fan run (L294b) looks good, add contacts lofts <-> vaginal wall / PM arch / levator, one at a time; AVW-Para must contact the PM arch and all wall parts; CL / USL left out*

The user (2026-10-16): once `L294b_junction_chains_fan_MAX` (tube springs + lofts version of the L283 junction, with the
curve-tip -> PM fan lofted as one sheet per side) looks good, add contacts between ALL the loft sheets and (1) the vaginal
wall parts (AVW, PVW, cervix, curves, PeB), (2) the PM arch (PM_Plane), (3) the levator ani (LA shells) - one at a time,
each its own variant; START with AVW-Para <-> PM arch (the user, 10-16). Leave out CL and USL ("they don't intersect anything"). The AVW-Para sheets (both sides) DEFINITELY
need contact with the PM arch and every part of the vaginal wall.

**Why:** the lofts should not pass through the tissues they wrap; adding one contact family at a time shows which one
causes run trouble.

**How to apply:** shell-to-solid sliding-elastic contact is ONE-SIDED (tested 10-16 on a 2-block file: only the face the
shell normal points to engages; shell_bottom_primary 1 flips it; two_pass does not fix it), so each loft pair needs a
top-face and a bottom-face contact. Shell-to-shell (loft <-> PM arch / LA shells) WORKS (corrected 2026-10-16 evening; the first
"never engages" test fixed the whole target): wind each contact <Surface>'s facets to face the other part (surface winding
decides, not the element normals); tests in claude_diag/shell_tests (gen_ss2.py). If a step causes trouble, try AVW-Para alone first. Related:
`feedback-three-base-versions`, `feedback-febio-run-workflow`.

2026-10-17: shell lofts are never caught when they meet the PM edge-on (shell_tests/gen_edge.py), and a one-pass contact
misses the PM's own nodes; search_radius 3 (not 1) stops lost nodes. A thin SOLID sheet is caught edge-on: the user OK'd
turning the sheets into thin solids "if we need to" (tools/build_batch312.py, runs L312a/b). Solid ropes fold where the
AVW-Para connectors go slack at MAX (not usable for that family).

2026-10-17 evening: the user on L312b (solid lofts + PM contact): a solid loft that dips into the PM and is pushed back out
is "probably okay" (judge a contact by whether it recovers, not by zero transient penetration); the user suggested the PM
arch as a thin solid too (tools/build_batch315.py, runs L315a/b; extruded on the side away from the lofts). No shell
contact setting holds the PM's free edge (L296j two_pass, L296k reversed facets).

### `project-max-lite`
*From 2026-10-18 the user: MAX goes beyond what is needed; test MAX-lite cases (less-impaired supports) that still give a lot of deformation*

The user (2026-10-18): "They should still produce a lot of deformation, but I think the current ones go beyond what is
necessary. That's fine if they don't error out, but we may not need to push that far if it's causing problems."

MAX = levator 40 %, CL / USL x0.2, AVW-Para x0.5, P-arcus x0.05 (tools/build_batch244.py CASES). It overshoots the range
targets (Ba / Bp >= +10, apex motion >= 35 mm; `project-range-of-motion-goal`): L288c +21.2 / +17.7, L294 +18.9 / +25.0.
The chains-line curve folds (snap at t 0.73-0.78) appear only at MAX (milder cases <= 0.24 mm).

**Why:** a MAX that overshoots is what drives the curve-end snap; a lighter case that still reaches the targets avoids it.
**How to apply:** test MAX-lite variants one less-impaired support at a time (e.g. P-arcus x0.2, CL / USL x0.5, levator 60 %),
and per-target cases (anterior-heavy for Ba, posterior-heavy for Bp); judge on targets met + folds < ~1 mm + healthy up.
Keep MAX runs that work; do not drop MAX. Say whole-tissue vs local (`feedback-say-whole-vs-local`).

### `project-paper-sync`
*From 2026-10-19 a paper on the new PVP3D model is being written; log paper-relevant decisions in the addendum as we go, tag runs with run_notes --paper, sync into the methods draft at batch end*

The user (2026-10-19) is writing a paper on the new model (how it differs from the Abaqus source / first conversion, new
features, anatomical accuracy, technical work). A separate session ("PVP3D paper write-up and model-evolution map") wrote a
methods draft (Claude Doc) and the PVP3D Run Timeline artifact (https://claude.ai/artifact/3nakpQE8zrxnutiu5tXS4c).

**The loop the user approved:**
1. Every paper-relevant DECISION (a change adopted, rejected, or still open, with its evidence runs) gets one line in
   claude_diag/PAPER_ADDENDUM_since_2026-10-19.md: time | run / step | what | result | status (adopted / rejected / open).
   Routine progress stays in README section 8. Plain words, NOT IN SOURCE / NOT FAITHFUL marked.
2. Tag runs: `py -3.10 tools/run_notes.py RUN --paper methods|supplement|skip` (the `paper` column of run_log_notes.csv,
   carried into PVP3D_run_log.xlsx by run_log.py).
3. At the end of each batch (or daily): the paper session folds the new addendum lines into the methods draft (adopted ->
   methods, rejected -> "what was tested" supplement), refreshes the timeline from the rebuilt run log, and the addendum gets
   a "synced up to <time>" line so nothing is added twice.

**Why:** the user wants the paper kept current without re-reading every record each time.
**How to apply:** do 1 and 2 as runs are judged; at batch end (or in the hand-off) trigger 3 (message the paper session, or
start one) after rebuilding the run log. Related: `project-pvp3d-handoff`, `feedback-run-log-spreadsheet`.

### `project-pvp3d-handoff`
*Where the PVP3D FEBio project's session hand-off lives, which skill copy to use, and the GitHub share*

Each session ends by writing NEXT_SESSION_PROMPT_<date>[b|c].md in the project root (latest: NEXT_SESSION_PROMPT_2026-10-07.md: the band / tube defaults, the soft PVW-join tube test L203, the open questions; before it NEXT_SESSION_PROMPT_2026-10-05.md: an 8-10 h unattended fit the user authorized on the non-tube seam-band base with tightened healthy supports + sphincter side x0.25; before it NEXT_SESSION_PROMPT_2026-10-03.md, the session of 2026-10-02 to 10-03, which starts with decisions for the user; before it NEXT_SESSION_PROMPT_2026-10-02.md; before it NEXT_SESSION_PROMPT_2026-09-29.md,
the session of 2026-09-28 night to 09-29 midday; before it 2026-09-28.md), plus an "Attempt N" section in SESSION_SUMMARY.md
(latest: Attempt 24) and a claude_diag/README_<date>[b|c].md (latest: README_2026-10-06.md, covering 10-06 and 10-07).
Skill backups from 2026-09-26 on go to F:\Developing PVP3D Model for FEBio\skill_backups (C: is nearly full).
Use the LOCAL skill at ~/.claude/skills/abaqus-febio-fea-pipeline, not the older account-synced anthropic-skills copy;
the prompt lists grep checks that prove the copy is current, and skill backups are zips in Downloads.

A private GitHub repo https://github.com/mtgordon/PVP3D_Skill shares the skill, the Abaqus source, the tools and the most
stable model with another researcher; the local repository is F:\Developing PVP3D Model for FEBio\github_PVP3D_Skill.
Pushing from Claude is blocked by the permission guard, so the user runs the push; strip personal paths from anything
added there. The user's terminal is Windows PowerShell 5.1 (no `&&`): give `git -C "<repo>" push`.
Last synced 2026-10-12 (committed locally, 81bb3b0 skill + b360711 models / tools / notes; the user pushes): tube main L258_tubecand_*, junction candidate L283_junction_fine2z18_*, ledge L228 + L271 sheet30, READMEs 10-06 / 10-09, all tools.
Before that, synced and pushed 2026-10-05: skill (24a3f26), then all tools (personal paths made repo-relative), L149_seamspr_pm_outer
(+ _fast), L159_C3full_parcus042_{healthy,P1,P2}, README_2026-10-02.md, updated README.md (3904c32).

Named models (user, 2026-09-25): the fast-base reference is L21_lofts4_pen_r2 (the user's four fitted lofts, soft lofts as
connectors); **"the Abaqus replica" is L19_pen5** (every Abaqus connector a spring, no lofts) — keep track of it so the
work can come back to it; the Abaqus U/V comparison model is L21D_allconn_pen (the replica like for like). Until the
LA's larger-than-expected displacement is explained, the user wants reduced LA force (Load-LA x 1/3) on both lines (the
springs line and the lofts line), and does not want that question worked on yet. Planned soon: add the lofts back one at
a time to the springs line and measure each one's effect on performance. Every new model also gets the source's
display-body part PM_Plane (reference only; added as a fixed FEBio rigid body, build_batch36.py add_display_body).

2026-09-26: the user wants the perineal body (_PickedSet66) to rotate about x, cranial end +y and -z, "up to 90 degrees",
by changing connective tissue, not the vaginal-wall materials (AVW, PVW, cervix, PeB), on both lines; told to keep pushing
toward 90 after ~40-45 deg was reached. The lever is P-arcus (then the side sphincters, CL/USL). The user approved the
lofts line with the two AVW-Para lofts replaced by their Abaqus connectors (the lofts crawl in compression once P-arcus is
weakened); CL, USL, PM stay lofts. Details: claude_diag/README_2026-09-26.md.
Later on 2026-09-26 the user said they can't run the Abaqus file, so no Abaqus U/V output is coming: the comparison
against L21D_allconn_pen is off the task list unless they say otherwise. At the end of that day: "I think I've learned
enough of how to make it rotate now". No change set was adopted; the lines stay L26_springs_la3_pm / L26_lofts_la3_pm. Next
session (NEXT_SESSION_PROMPT_2026-09-26b.md): (A) which fitted loft makes the lofts line faster than the springs line; (B)
why runs stall when the PVW first touches the LA (the PVW_LA contact onset).

2026-09-27 (session after NEXT_SESSION_PROMPT_2026-09-26b.md): the user asked for ~8 h of autonomous work and **reopened
the LA question**: find why the LA moves much more in FEBio than in Abaqus by comparing how the two models are built; the
user's guess is the LA's thickness (in Abaqus "a surface that has a thickness", in FEBio a shell with a thickness) is not
treated the same. Runs allowed to test it. The user also asked to read the FEBio knowledgebase/forum and keep notes in a
document: the Claude Doc "PVP3D FEBio notes: shells, contact and the LA"
(https://claude.ai/code/artifact/18b2a765-6850-4206-9660-fa08fbb42fd6). Findings (details: claude_diag/README_2026-09-27.md):
thickness treatment is not the cause (panel test within 4 %; full model mid-surface shells -8 %, back face held -5 %,
LA x2 thick -20 %); the LA hangs on the end-pinned, initially slack truss chains (converted exactly). The ~10 mm is a
recollection, not in the paper; the paper's POP-Q cases give numbers instead (tools/pop_q.py): FEBio's Bp tracks the paper
(~5 mm high), Ba is 13-20 mm too low. The vaginal-wall and perineal-body Ogden fits do not match the source's Marlow data
(walls 3x stiff at small strain, 2.3x soft equibiaxial at 1.5; PeB 3-4x soft): refit runs L70_* started 03:00.
Test A answer: the AVW-Para loft (it also holds the AVW ~9 mm higher: the lines are not like for like).
The user (08:00): their concern is the LA's **stretching/ballooning**, not its motion (area x2.16 at the source loads, which
the source's own LA law requires at 0.014 MPa; measure stretch, not displacement, for the LA question). They do not know
which version of the Abaqus .inp they were given, so the paper's 0 % connector baseline is unknown. Their reference is the
paper's Figure 3 and "a similar model I've seen run". Decisions: try things that might reduce the ballooning; materials
"stay true to the material properties from Abaqus, but adapt them for how the FEBio model is different (lofts vs trusses,
etc.)" = the Yeoh refits for walls and PeB (as the LA), lofts keep their connector-fitted materials, chains exact.
Still open for them: PVW_LA contact option (penalty 0.5 vs the 30 pinch facets out; recommended the latter), which line.
Decided 2026-09-27 08:40 (session after NEXT_SESSION_PROMPT_2026-09-27.md): **remove the 30 pinch facets** (the PVW_LA
primary facets holding a P-arcus connector end; NOT IN SOURCE; penalty 5 kept) "for now"; **the springs line is the
comparison line**, the lofts line may be used for faster screening runs. Then ~12 h of autonomous work (to ~20:40),
the user may check in during it. Outcome (claude_diag/README_2026-09-27b.md, SESSION_SUMMARY Attempt 17): option (b) alone
crawls wherever the PVW meets the LA earlier (Load-LA 1/3, softer walls): the PVW edge row touches at other facets; PVW_LA
seg_up 2 fixes it (new springs line L77_springs_la3_nopinch_pvwsegup2, same solution as L26); the faithful PeB runs; the
faithful (Yeoh) walls still crawl (canal seam + PVW edge), best with PVW_LA penalty 0.5; ballooning: healthy LA x1.60,
chains held x1.37 (PVW edge then pinches), Load-LA 1/3 x1.22, both x1.10.
The user's decisions (2026-09-27 ~20:00): (1) PVW_LA on both lines = the 30 facets out + penalty 0.5 + seg_up 2 (bases
L81_springs_la3_nopinch_pvwsegup2_pen05, L82_lofts_la3_nopinch_pvwpen05_segup2); (3) ballooning: the healthy LA from the
source's own PCM-LA_Yamada100% table (Yeoh c1 0.02347128 c2 0.03283498 k 2) AND Load-LA 1/3; (4) keep pushing the
faithful Yeoh walls (FEBio 4.13 has "contact potential": mini model first); (5) no runs left overnight. PeB refit
adopted. **The new lines: L83_springs_la3_la100_pebyeoh and L83_lofts_la3_la100_pebyeoh** (t = 1, LA area x1.09).
Contact potential failed on a mini model once surfaces slid (parked); batch 88 tests canal settings on the walls.
2026-09-28 (README_2026-09-27c.md): the faithful walls stall FEBio's implicit solver at the canal's contact-only side seam
(AVW/cervix and PVW share no nodes; frictionless contact, as in the Abaqus source, which runs *Dynamic, Explicit). No
setting fixes it; a seam tie or springs >= 0.01 N/mm converge but lift the walls (Ogden line Ba -3.3 -> -14.6 / -17.7) and
in the paper's cases stop the posterior wall responding (P1->P2 Bp +0.6 vs paper +5); the user disliked that ("I don't like
that the ties are affecting the results that much"). Softer springs and a tension (sliding) seam do not converge. FEBio's
explicit-solid runs every model feature in mini tests but its dyn_damping over-damps the heavy chain; ~7-14 h a run. Lines
unchanged (L83_*). The user decided: Ogden walls on the lines "for now when testing other things"; next context window
(started by the user): "keep investigating how to get softer vaginal wall materials without it causing problems", first the
Yeoh walls with c1 x1.5 and x2 (NOT FAITHFUL, "most faithful that still converges"), then maybe fibre (directional)
stiffness. The user then adopted solver rhoi 0 on both lines ("go ahead and change rhoi"): **the lines are now
L87_springs_newline_rhoi0 and L87_lofts_newline_rhoi0** (Ogden walls; bases in run_log.py BASES). The run log has a Base
models sheet (run_log.py BASES: add a row whenever the user adopts a base).

**Why:** the user continues across context windows via these files and asked for the GitHub share on 2026-09-24.

**How to apply:** start from the newest prompt file; before editing shared files, check for other sessions and
modification times. Related: `feedback-febio-run-workflow`.

2026-09-28 evening (README_2026-09-28.md, Attempt 19): the walls converge on the springs line down to Yeoh c1 x1.5
(within 10 % of the source's data at 5-10 % strain) and with a 3 mm c1 x3 band along the canal seam (answer within ~1 mm of
the faithful walls); the lofts line stalls with any softer wall. Paper cases on the springs line: Bp close (+1.6 / +7.8 vs
+4 / +9), Ba 13-32 mm below the hymen in every run (the paper: above); the gap is the anterior wall's support, not its
material. Lines unchanged (Ogden walls); no softer-wall option adopted (the user's decision, on hold). Batch 108 (the
distal anterior support: PM_avw_bottom connectors x1e5, NOT IN SOURCE as given) was handed over running. **The user's next
main task: a third line, "fibre lofts"**: the springs line (no lofts), the lofts line, and now every family ever lofted
(the 14 *_fan domains of L9_fitall_edge) as a weak matrix (a few % of the fibres) + fibres carrying the source connector
laws; Ogden walls; benchmarks: the springs line, the paper's P1 / P2, PeB rotation. The user was away ~8 h from 22:15.

2026-09-29 (README_2026-09-29.md, Attempt 20): **the user adopted a new springs-line base, L91_newline_vwyeoh_rhoi0_seamspr001**
(faithful Yeoh walls + 0.01 N/mm canal seam springs; "for now"). Built, not adopted: fibre lofts (too compliant in situ; a
fibre attaches only through its element's nodes), connector chains in a weak sheet (L111_*: matches the springs line; the
user's "option 3"), the PM (PM_Plane) as a clamped Yeoh shell with the PM-family anchors tied to it (L113/L114; the user's
wish: it should be a real structure, not a reference), the canal as a tube (L118_tube_fold, L119_tube_wrap: converges 3x
faster, Ba -14.0 vs -9.0), Jing et al. 2012 levator / PeB laws (L120_*). Next (the user): make the tube open more.

2026-09-29 to 10-02 (README_2026-09-30.md, README_2026-10-01.md, Attempt 21): the tube with a 3 mm rounded lumen end
(L123_tube_round_r30) opens the seam; SlidingElastic1 seg_up 2 fixed its slow P1 (same answer); the user adopted
L127_tube_r30s, then **L128_tube_r30s_pm** (+ the PM as a structure) as the springs base, with chains / lofts versions
L132_chains_tube_pm / L132_lofts_tube_pm (the lofts line now reaches t = 1 with faithful walls). Wall-stiffness scaling
cannot bring Bp to the paper (the gap is posterior support). Then, on the user's requests: 2.75 mm walls with the
channel opened, narrowed 3 mm (L137), curves to the PeB end embedded in the PVW / PeB (L143_thin275_narrow3_wrapend_noflare,
P1 / P2 = L130_thin275_wrapend_paperP1/P2: Ba above the hymen, Bp -12 / -10). Not adopted as base yet (the user: "not
yet"). Next (NEXT_SESSION_PROMPT_2026-10-02.md): a new base = L91 seam springs + PM structure + PM_conn inner arc + all
contacts (Yeoh walls), then fit the paper's P1 / P2 by changing healthy connective tissue and wall stiffness.

2026-10-02 (README_2026-10-02.md): new springs-line base **L149_seamspr_pm_outer** (L91 seam springs + Yeoh walls + seg_up 2 +
PM structure + the 27 PM_conn springs on the PM's clamped OUTER arc (the user's rule from now on) + walls_LA / walls_PM); fast
learning copy L149_seamspr_pm_outer_fast (no wall contacts). Not yet formally adopted into run_log.py BASES, and no chains /
lofts versions built (ask). The user chose one shared wall scale (AVW = cervix = PVW). Fit result within x0.2-x5: Bp OR Ba,
not both. C3 = healthy Parcus x0.42: Bp +5.0 / +9.2 (P1 / P2, full base), Ba 9-15 mm below the hymen. C2 = Parcus x0.42 +
walls x5: Ba at the hymen, Bp 10 mm short. All levers move Ba and Bp together (seam springs couple the walls). Runs with
dozens of retries can converge to spurious upward jumps: trace before using an end state.
2026-10-03: seam tests on C3: seam springs x0.5 and removing them (with the edges held in x) both drop Ba (the springs hold the AVW up on the PVW); the x constraint makes the cervix tip crawl. Session closed with NEXT_SESSION_PROMPT_2026-10-03.md (decisions first: walls that differ / accept C3 / anterior support / bounds; adopt L149 + its three versions).

2026-10-05 (README_2026-10-05.md, Attempt 23): tube line (L162-L173): wall-only mass damping C 200/s needed for walls x0.35;
right curve fixed to reach the end; 5 springs a side on the curves' distal end ring to the PM arc; the curves' distal ends fold
through themselves from t ~0.8 (not yet functional; refined ends + self-contact penalty 1 = L172 is the candidate fix). Non-tube
contact-only seam + seam band (L167/L168): Bp on target, Ba 28-30 mm low; the user wants tightened healthy supports, not
whole-wall stiffening. New rule: LA_sphincter_side_conn_mat x0.25 in every new model.

2026-10-07 afternoon (README_2026-10-07.md): the user's seam lips (tools/seam_extend.py: unloaded 3 / 5 mm extrusions of the
walls past the contact-only seam, no band) fix the seam stall (control crawled and was stopped at t 0.775; lip runs t 1 in
~30 min) but without the band the AVW drops (P1 gap -3.3 vs the band's +3.9); posterior x0.25 only gets -1.2. Decision
pending: adopt lips + refit, or keep the band default. queue_runs.sh now locks (runs/.queue_lock). Another session runs
the tube's fully elastic join (L204_apexclosed_pvwsoftjoin_full_P1, its own build_batch204.py).

2026-10-08 night audit (README_2026-10-08.md): with lips x0.1 the SOURCE supports run (L209_src_sl01_*): Bp meets P1 / P2
(+4.5 / +9.2) unfitted, Ba 25-30 mm out. Audit: Load-LA, damping, LA 50 % change nothing; source chain mass (line x10) moves
Bp only (P1 +7.9 at t 1, settling ~+6). Cause of the Ba miss: the AVW's paravaginal springs (35 per side to ONE pinned RP)
shorten as the AVW descends and, tension-only like the source, carry nothing; with compression (diagnostic L215_paracomp)
Ba -5.7. Apical alone saturates the gap at ~-15. Decision pending for the user (README section 6).

2026-10-08/09 (README_2026-10-08.md sections 7-11): CL holds the AVW, USL the PVW (build_batch145 cl= / usl=); Parcus drops
Bp. The PM shell's 0.36-0.68 mm-wide quads caused ~88 % of the negative Jacobians: tools/pm_coarsen.py (rows 0/5/10) -> 2.6x
faster, answer within 0.5 mm; new base L228_ext5sl01_pmc_pmd_sph025. On it, CL x35 + Parcus x0.59 (L230_pmc_cl35parcus059_src_*)
MEETS healthy -3.5/-16.4, P1 -0.5/+4.9, P2 Bp +8.2 (Ba +3.3, accepted) -- not yet adopted as a default (user to decide).
last_moves.py prints plot index AND node id (the "lip nodes" blamed for late steps were PM_Plane nodes).

2026-10-09: hand-off NEXT_SESSION_PROMPT_2026-10-09.md (the user asked for a fresh context). Open task: close the tube's apical
corners (10 user descriptions quoted there; 7 rejected attempts; screenshots claude_diag/apex_screenshots/). Show the user a picture
of the intended SHAPE and get a yes before building. Cervix x3 lets CL drop to x20 on the lips line (P2 Bp met, Ba at the hymen).

2026-10-09/10 (README_2026-10-09.md): new goal = range of motion (`project-range-of-motion-goal`). Candidates, NOT adopted: ledge L246_ledge_cx3_paracomp03 (cervix x3 + AVW-Para compression x0.3); tube L258_tubecand (curves x0.6, sphincter x0.25, no curve pressure). Hand-off NEXT_SESSION_PROMPT_2026-10-10.md. queue_threads.sh counts threads (mixed 12 / 8-thread runs).

2026-10-10 afternoon: hand-off NEXT_SESSION_PROMPT_2026-10-10b.md (the user away 3-10+ h; runs L274 MAX, L276 noembed / spreadembed, L277b / L277c self-contact going). The user suspects the PVW-end / curves / PeB junction needs remeshing or re-joining (conforming junction proposed, picture first).
2026-10-10 evening: hand-off NEXT_SESSION_PROMPT_2026-10-10c.md (all runs done: the curve-end point ties cause the tearing / fold, no self-contact setting helps; ledge 25 % full set; junction rebuild and 25 vs 30 % for the user).
2026-10-11 night: hand-off NEXT_SESSION_PROMPT_2026-10-11.md (junction rebuilt: no PVW flare, curves tied to the PeB, curve rings split; best L283_junction_fine2z18_*, fold-free, limits met; its healthy / F8A / P1 running).
2026-10-12: hand-off NEXT_SESSION_PROMPT_2026-10-12.md (session set up after the junction work; skill updated with item 14 / gotcha 30 / element_volume_report.py; GitHub commits waiting for the user's push).
2026-10-14/15 (user away ~10 h; README_2026-10-09.md section 8): the curves of L283 MAX (and every junction since L258) fold through themselves ~1 mm at the distal (PeB) ends from t 0.65 (tools/curve_selfpen.py). Changing the WHOLE curve self-contact stalls at the curve-top / apex_contour junction; build_batch285.py gives the distal facets (rest z < -15) their own self-contact. Radius, penalty, AUGLAG, node-on-facet, seg_up 0: no good; two_pass 1 (L285_distal_tp_*) holds it to ~0.6 mm at t 0.68 (best; MAX + healthy/F8A/P1 running at hand-back).
2026-10-16: hand-off NEXT_SESSION_PROMPT_2026-10-16.md (curve fold-through work 10-14..16 in README_2026-10-09 section 8; spring patches build_batch288, damping build_batch289, springs+lofts L294 / fan L294b running; next: loft contacts, memory project-loft-contacts-plan; skill updated: gotchas 31-32, 4 lessons).
2026-10-16 night / 10-17 (the user away ~11 h; README_2026-10-09.md section 8 end, 'Night summary'): no MAX run reached <= 0.3 mm curve self-penetration; on chains every lever failed (patches, damping crawls, two_pass, wall contact, ring split to 30 mm); best springs L288c / L301 0.6-0.9 mm but L301 blew up at t 0.981. Loft contacts: shell-to-shell works (surfaces must face each other); AVW-Para <-> PM holds only at penalty 0.5 (crawl at the PM free edge; L296f rerun mostly past it); AVW-Para <-> walls cannot hold where curves collapse. Tools: build_batch296/297/302, loft_contacts/*.py, curve_selfpen_study/*.py.
2026-10-17: hand-off NEXT_SESSION_PROMPT_2026-10-17.md (curves wrinkle on the compressed side of the ~90 deg bend; <= 1 mm accepted; loft <-> PM: shells are never caught edge-on, one-pass misses PM nodes, sr 3 helps -> solid lofts L312a/b running; ropes fold when slack; L309 heavy near-zero lofts, L310 fold contact running; skill gotchas 32-33 updated, 72 lessons rows).
2026-10-18: hand-off NEXT_SESSION_PROMPT_2026-10-18.md (chains-line curve collapse = one-step snap; L309 near-zero lofts never snap -> L323 curve-end loft band; solid-loft contacts: PM (L312b, crawls at the PM free edge; L324 curled lip, L326 mass), walls L318, P-arcus-LA L321b, fan-PM L322; PM as solid fails; skill gotcha 34 contact potential).
2026-10-20 morning: hand-off NEXT_SESSION_PROMPT_2026-10-20.md (overnight: the springs-line AVW-Para fibre lofts crawled because two near-zero-stiffness fibre lofts (strap + AVW-Para) meet at the curve tip; candidate L348b = AVW-Para loft matrix x10, strap matrix x3, tip corners on patch anchors, AVW-Para <-> PM contact (build_batch343) + 12 mm band + curled lip with lip damping; chains-line lip + band L343a holds to t 0.96; tools/tip_kinks.py).
2026-10-20 evening: hand-off NEXT_SESSION_PROMPT_2026-10-20b.md (direct strap-edge tie fixes the tip dimple; the curve collapse = the distal curve buckling into a hairpin between the descending wall and the PeB, whose rotation stalls at ~28 deg from t 0.6; the user's three tests on L352: more self-contact faces, pre-bent distal curve 4 mm, LA_sphincter_side 0.1).
2026-10-21 night / 10-22: hand-off NEXT_SESSION_PROMPT_2026-10-22.md (the user away 22:10 -> morning, run freely: candidate line capped at Ba ~+1; springs-line L353 + 1/8-1/4 curve pressure reaches the targets then the AVW surges / folds back; AVW damping C 2000 = L393 Ba / Bp +8.2 / +9.3 at t 1, clean; holds L388b (user's, on L384) and L394 (on L393) running; new tool midline_bulge.py; keep run paths <= 255 chars).
2026-10-22 afternoon: hand-off NEXT_SESSION_PROMPT_2026-10-22b.md (curve_selfpen under-reports -> tools/curve_cross.py; every run crosses > 1 mm late; strap / fillet / curve-refine variants L402-L411; the user's hold L388c running).

### `project-range-of-motion-goal`
*From 2026-10-09 the goal is range of motion, not the paper's numbers - Ba/Bp/apex must reach at least the paper's values, healthy unprolapsed; remeshed PM + cervix x3 + source stiffnesses are the new standing settings*

From 2026-10-09 (the user): "We're no longer trying to match the numbers from the paper. What we do want to be able to do is get
about the same range of motion ... make sure the movement of our model isn't too restrictive." "Range is just that the Bp and
Ba can get at least to the points that they got to in the paper. If the paper has apical location check for that too. Healthy
case should stay unprolapsed." Both lines (ledge and tube).

Standing settings for new models: the remeshed PM (tools/pm_coarsen.py) is THE mesh; the cervix x3 of the wall is normal (see
`feedback-cervix-stiffer`); every support and wall back at the source value (drop the CL / Parcus fit scales); a material may
differ from the source if the range of motion needs it (label NOT FAITHFUL). Kept: Load-LA x1/3, sphincter side x0.25
(`feedback-sphincter-side-softer`), arcus chain mass x10, ledges x0.1, PM structure + PM_conn outer arc. Tube walls stay
x0.35 "for now, revisit later".

Targets (Luo et al. 2015, 140 cm H2O; tools/build_batch243.py cases): healthy Ba, Bp < 0; P1 Bp >= +4; P2 Bp >= +9; Fig 7A
rectocele Bp >= +10; Fig 8A cystocele Ba >= +10; Fig 6 apex motion >= 25 mm (no apical impairment) / 35 mm (apical 60 %),
anterior 24 / 29 mm, posterior 41 / 50 mm of motion. First batch: L243_ledge_cx3_src_* (7 cases).

**Correction (the user, 2026-10-09 evening): "we're not trying to match P1 and P2 anymore. We're just trying to be able to
get those ranges, but maybe with different impairments."** The paper's cases are NOT tests to pass one by one: the model passes
if SOME impairment (any levels, any combination) reaches each position (Bp >= +10, Ba >= +10, apex motion >= 35 mm, the most
the paper shows), while healthy stays unprolapsed. Do not tune supports to make P1 / P2 hit +4 / +9.

**Why:** the user wants a model whose motion is not over-constrained, rather than a fit to single published values.
**How to apply:** judge variants on reaching these limits (and healthy staying up), not on matching P1 / P2 exactly. Related:
`project-default-bases`.

## References

### `reference-normal-generic-model`
*PVP3DModel_job.inp is THE source; Normal_Generic_copy (1).inp (Downloads) is only a reference for another model's CL/US tissue properties*

The user, 2026-10-08: "I'd view the pvp3dmodel as the source material that we are working from. I just gave you the other
one to look at a different model with CL tissue properties."

- **Source (judge faithfulness against this):** `F:\Developing PVP3D Model for FEBio\PVP3DModel_job.inp`.
- **Reference only:** `C:\Users\mgordon\Downloads\Normal_Generic_copy (1).inp` ("Normal_2_withoutPMLAtie", Abaqus/CAE 2016,
  anterior-compartment OPAL325 model: AVW on a GI-filler solid, no PVW). CL / US / paravaginal / ATFP are T3D2 truss chains,
  1 mm2, 0.5 mm segments, welded (CONN3D2 Weld) to the AVW and pinned at the far end; sections point to undefined
  `CL_Analog` / `US_Analog` (placeholders for CLUS_H1 / H5 / H10 / H50 / H100). CLUS_H100 (the user's pasted table): nearly
  linear, E 0.027 -> 0.046 MPa over 0-100 %, no compression. CL 6 chains a side (50-83 mm), US 3 a side (43-50 mm): per side
  26-140x SOFTER than the PVP3D source's CL / USL connectors. Paravaginal: 6 chains a side (9-22 mm) fanning to ONE pinned
  point per side (the ATFP's anterior end), the same topology as the PVP3D source's BC-RP-ParaSupport.

Say "the PVP3D source" and "the Normal_Generic file", never "the original", to avoid ambiguity.
Related: `project-default-bases`, `feedback-lips-minimal-influence`.
