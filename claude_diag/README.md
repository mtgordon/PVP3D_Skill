# LA convergence diagnosis — 2026-09-22

Question: the no-LA model runs to t≈0.60, but no LA version got past ≈0.25, and
`PVP3DModel_v21_v3_HYBRID_withLA_TRUSS.feb` dies at t≈0.0002. Why?

Everything below was tested by running FEBio 4.13 (4 threads per run), one change
at a time. The variant builder is in `tools/variants.py`, the runs are in `runs/<name>/`,
and each `<name>.feb.changes.txt` lists that variant's exact edits.

## 1. Why the TRUSS file fails at t = 0.0002

| Cause | Evidence |
|---|---|
| Each **phantom node** is held by a single linear spring (one direction only). The two other directions are held only by the tied contact, whose penalty is **ramped from 0** (`lc="2"` TieRamp), so at t≈0 they have almost no stiffness and the matrix is singular. | Residual 9e8, displacement 7.7e6 in the first step. |
| The **Posterior Arcus ties use LA_ICM only** as primary. In Abaqus the slave surface is **all 1170 LA elements**. The upper arcus (truss nodes 1–5, 15; 3–24 mm from ICM) never engages, so those chain nodes and their phantoms are unsupported. | `tools/analyze_truss.py` distance table. |
| **`tied-node-on-facet` closes initial gaps.** It pulls secondary nodes onto the primary, unlike Abaqus `*Tie, adjust=no`. Any engaged node that sits off the surface gets snapped (the old 0.22 vs 0.25 `max_distance` knife-edge). `node_reloc=1` removes the snap. | `mini/tie_gap05*.feb` (20-line test models). |

**Fix used (T-series):** delete the phantom nodes/links and the 4 tied contacts. Add the Abaqus
ties exactly as FEBio linear constraints (`<Constraints><constraint type="linear constraint">`,
u_LA = u_truss in x, y, z). Each LA slave node is tied to its nearest truss node within 4 mm:
18 + 18 LA nodes on ATLA, 82 + 81 on Posterior Arcus (from `_PickedSet473/475` and
`_PickedSurf737/739` in the .inp). Verified in `mini/lc*.feb`. The constraints converge in
2–3 augmentations per step. In the fixed model at t = 0.56, |u_LA − u_truss| is at most
1e-6 mm over all 199 ties, and the truss chains move 0.23 mm (ATLA) and 0.39 mm (Arcus),
symmetric left/right.

## 2. Why LA versions stall at ≈0.13–0.27

1. **The 4 mm LA shell is pinched at rest.** FEBio shells integrate a 3-D volume through
   the thickness along *averaged nodal normals*. At the tight folds of the posterior LA
   sling, a 4 mm thickness is more than the local fold radius allows. **6 LA elements are
   inverted on one face before any load is applied, and 12 more are nearly pinched**
   (face area < 25 % of mid-surface). The elements that invert in every LA run are these:
   23695/23694 (T3 diagnostic), 23617/23769 (T2), 23727 (E2). Abaqus S4R conventional
   shells don't form that 3-D Jacobian, so the source model never noticed.
   **Fix: `<shell_normal_nodal>0</shell_normal_nodal>`** in the four LA ShellDomains
   (element normals). The thickness stays 4 mm and no element is pinched.
2. **LA-to-perineal-body coupling.** Abaqus uses 40 soft axial connectors
   (Sphincter-L-1..8, R-1..8, P-1..24; `LA-Y-sphincter-side-*`, `LA-Y-spincter-p`). The
   FEBio model merged 16 LA nodes into the stiff E = 21 Sphincter-L/R fans and left the
   Sphincter-P fan (24 connectors' worth) unattached. All 80 connector endpoints exist
   exactly (0.0000 mm) in the FEBio mesh, so the connectors are modelled directly.
3. The LA is **not** a problem as a body: un-merged and unloaded (E1), the model is identical
   to no-LA (0.6043).

## 3. Pre-existing conversion gaps found along the way (they affect the no-LA model too)

- **Remesh-lost pins (since v3_HYBRID):** the MMG remesh gave the fan anchor nodes new IDs,
  and the BC node sets still point at the old (now orphan) nodes.
  `BC-RP-ParaSupport` went from 2/2 to 0/2 active, so the **AVW-Para (paravaginal) fan apexes are free**.
  `BC-RP-PeB-constrain-origins` went from 26/26 to 0/26 active, so the PeB-constrain fan dangles. Re-pinning by
  coordinate (apex nodes 16333, 22091; edge nodes 18502…18518) makes the no-LA model
  run **retry-free to 0.60** (N1), instead of needing 79 retries.
- **Missing PVW–LA contact:** Abaqus `Int-PVW-LA` is frictionless contact between the
  PVW back face and all LA elements. It is absent in FEBio. The PVW back face starts
  3.03 mm from the LA mid-surface (1.03 mm from its face). In the fixed model the gap
  was 2.77 mm at t = 0.38 (`tools/pvw_la_gap.py`), so the contact would engage around
  t ≈ 0.5–0.6. It does not explain the early stalls, but the model needs it before
  results past t ≈ 0.5 mean anything.
- **Analysis type:** the Abaqus step is `*Dynamic, Explicit` (1 s, bulk viscosity), which
  tolerates floppy parts that a static implicit solve cannot. Implicit DYNAMIC in FEBio
  with the Abaqus densities was tried and was *worse* (stalls at ≈0.07; D0/D1).
- **LA bulk modulus:** k = 1.13 came from the wrong μ₀ formula. The Ogden fit is correct
  (c1 = 2μ, μ₀ = c1/2 = 0.0214 MPa). For Abaqus `poisson=0.47`, k = 0.349.

## 4. Run table (t reached = last converged time)

| Run | Change (single step from the row it builds on) | t reached | retries | end |
|---|---|---|---|---|
| TRUSS (yours) | your TRUSS file | 0 | 21 | failed |
| TRUSS_NOTIE (yours) | your TRUSS, ties off | 0.0190871 | 35 | failed |
| CORRECTED_testing (yours) | your best LA run (ties inert) | 0.219372 | 101 | failed |
| E3 | no-LA v32 + LA solver settings | 0.604337 | 79 | failed |
| E0 | CORRECTED, unchanged (control) | 0.219706 | 80 | failed |
| E4 | E0 with v32 solver settings | 0.101359 | 21 | failed |
| E1 | E0, LA un-merged from sphincter fans (LA unloaded) | 0.604337 | 79 | failed |
| E2 | E1 + 8 Abaqus side connectors | 0.263072 | 43 | failed |
| E5 | E1 + all 40 Abaqus connectors | 0.234193 | 32 | failed |
| T1 | TRUSS, ties -> exact linear constraints (fan merge kept) | 0.01 | 13 | stopped (by me) |
| T2 | T1, un-merged + 16 side connectors | 0.267143 | 33 | failed |
| T3 | T2 + 24 posterior connectors | 0.201429 | 31 | stopped (by me) |
| T3s | T3 + remesh-lost pins restored | 0.204286 | 19 | stopped (by me) |
| T3f | T3s - redundant sphincter fans | 0.205 | 15 | stopped (by me) |
| T3f_k | T3f + LA k=0.349 | 0.202876 | 5 | stopped (by me) |
| T3f_3f | T3f + three-field LA shells | 0.205 | 6 | stopped (by me) |
| T3f_t15 | T3f, LA t=1.5 membrane-equivalent | 0.139286 | 3 | stopped (by me) |
| D1 | T3f, implicit DYNAMIC + Abaqus densities | 0.0758762 | 16 | stopped (by me) |
| **T3f_snn0** | **T3f + LA element normals (shell_normal_nodal=0)** | 0.564902 | 48 | stopped (by me) |
| **T3f_k_snn0** | **T3f_snn0 + LA k=0.349 -> TRUSS_fixed** | 0.564902 | 46 | stopped (by me) at the wall |
| E0n | E0 + LA element normals only (fan merge kept) | 0.19914 | 55 | stopped (by me) |
| N1 | no-LA v32 + remesh-lost pins restored | 0.600898 | 5 | stopped (by me) |

The **T3f_k_snn0** row is `PVP3DModel_v21_v3_HYBRID_withLA_TRUSS_fixed.feb`. It runs retry-free to 0.564, then crawls (steps < 0.001) in the same soft-mode stall as the no-LA wall at 0.604. The most-strained elements are the same in both (USL-L_fan, `_PickedSet346` elem 7115).

## 5. Tools

- `tools/variants.py` builds variants: `remove_truss_ties`, `add_abaqus_ties_lc`, `unmerge_la`,
  `add_all_sphincter_connectors`, `restore_remesh_lost_pins`, `remove_domains`,
  `set_la_shell_normal_nodal`, `set_la_material`, `set_dynamic`.
- `tools/xplt.py` reads the .xplt plot file. `tools/postmortem.py` finds where a failed run
  deforms or blows up. `tools/conn_watch.py` tracks connector loads.
- `tools/abq_ties.py` resolves the Abaqus LA–truss `*Tie` pairs. `tools/progress.sh` gives
  a one-line status per run.
- Per-element negative-Jacobian IDs:
  `printf 'set output_negative_jacobians 1\nrun -i model.feb\nquit\n' | febio4.exe`
