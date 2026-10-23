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
| **T3f_k_snn0** | **T3f_snn0 + LA k=0.349 -> TRUSS_fixed** | 0.564902 | 46 | running |
| E0n | E0 + LA element normals only (fan merge kept) | 0.19914 | 55 | stopped (by me) |
| N1 | no-LA v32 + remesh-lost pins restored | 0.600898 | 5 | stopped (by me) |
