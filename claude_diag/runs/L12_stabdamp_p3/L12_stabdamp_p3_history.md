# How L12_stabdamp_p3.feb was built: every change, oldest first

Each step lists the model it started from and its changes (from the `.feb.changes.txt` written by the build tools). "NOT IN SOURCE" marks anything that is not in the Abaqus model `PVP3DModel_job.inp`.

## PVP3DModel_v21_v3_HYBRID_withLA_TRUSS_fixed2.feb

Started from: (the earliest file in the chain; no change list)


## LP2DM_m10.feb

Started from: `PVP3DModel_v21_v3_HYBRID_withLA_TRUSS_fixed2.feb`

- added Abaqus Load-LA: pressure 0.014 (lc 1) on the top face of all LA elements
- LA material -> uncoupled Yeoh (I1-only, like Abaqus Marlow LA_Yamada50%): c=(0.0100393, 0.0172014), k=0.32795 (was Ogden density=1.06e-09, c1=0.042712, m1=6.488300, k=0.3488)
- LA bulk modulus k 0.32795 -> 1
- added 58 mass-only truss2 elements along the chain springs ('ATLA_Left_springs', 'ATLA_Right_springs', 'Posterior_Arcus_Left_springs', 'Posterior_Arcus_Right_springs') (linear truss, rho=0.00011, E=1e-06, area=1; Abaqus truss density)
- DYNAMIC analysis, rhoi=0.5, Abaqus densities x1.0
- chain mass x10 (mass scaling on the chains only)

## LP2DM_m10_settle.feb

Started from: `LP2DM_m10.feb`

- settle phase: run to t=2.0 with loads held, mass damping C=20/s ramped in over t=1.0-1.05 (lc 3)

## LPFmsTk_m10s.feb

Started from: `LP2DM_m10_settle.feb`

- format only: truss2 -> line2 + BeamDomain linear-truss for ['chain_mass'] (FEBio Studio)
- P-arcus-L_fan merged onto the arcus chain: 45 element-node references replaced (fan node -> chain node) [(13335, 23701), (13343, 23715), (13351, 23702), (13359, 23703), (13367, 23704), (13375, 23705), (13383, 23706), (13391, 23707), (13399, 23708), (13407, 23709), (13415, 23710), (13423, 23711), (13431, 23716)]
- P-arcus-R_fan merged onto the arcus chain: 45 element-node references replaced (fan node -> chain node) [(13487, 23720), (13495, 23734), (13503, 23721), (13511, 23722), (13519, 23723), (13527, 23724), (13535, 23725), (13543, 23726), (13551, 23727), (13559, 23728), (13567, 23729), (13575, 23730), (13583, 23735)]
- fan-edge nodes moved onto the chain line (distance mm, left unattached): [(14320, 1.0), (14325, 1.0), (14340, 1.0), (14343, 1.0), (14345, 1.0), (14349, 0.99), (14351, 0.99), (14372, 0.99), (14415, 1.0), (14420, 1.0), (14435, 1.0), (14438, 1.0), (14440, 1.0), (14444, 0.99), (14446, 0.99), (14467, 0.99)]
- P-arcus fans only: material Beam-CL-USL (isotropic elastic, density=7.8e-07, E=21, v=0.3) -> P-arcus_fan_soft: Ogden c1=0.25908, m1=3.5, k=32.38 (mu0 0.1295 MPa), density 7.8e-07; fitted to the 13 LA-Y-parcus connectors per side (fit_fan_to_connectors.py)
- P-arcus lofts: shell thickness 0.125 -> 0.49 mm (as the other lofts), c1, k and density x0.2551 (membrane pull and loft mass unchanged, bending x15.4)

## L5X_all.feb

Started from: `LPFmsTk_m10s.feb`

- run length only: time_steps 200 -> 120 (ends at t = 1.2); the model is unchanged
- NOT IN SOURCE: Load-LA pressure 0.014 -> 0.00466667 MPa (x0.3333); the other pressure loads unchanged
- tied-node-on-facet AVW_Para_L_fan__PickedSet346_tnof: penalty 0.0005 -> 100
- tied-node-on-facet AVW_Para_L_fan__PickedSet347_tnof: penalty 0.0005 -> 100
- tied-node-on-facet USL_L_fan__PickedSet64_tnof: penalty 0.0005 -> 100
- tied-node-on-facet USL_L_fan__PickedSet346_tnof: penalty 0.0005 -> 100
- P-arcus-L loft edge nodes on the chain line made chain nodes (8; node, distance mm): [(14320, 0.0), (14325, 0.0), (14340, 0.0), (14343, 0.0), (14345, 0.0), (14349, 0.0), (14351, 0.0), (14372, 0.0)]; Posterior_Arcus_Left_springs 18 -> 26 springs
- P-arcus-R loft edge nodes on the chain line made chain nodes (8; node, distance mm): [(14415, 0.0), (14420, 0.0), (14435, 0.0), (14438, 0.0), (14440, 0.0), (14444, 0.0), (14446, 0.0), (14467, 0.0)]; Posterior_Arcus_Right_springs 18 -> 26 springs
- removed 328 stab_* ground springs on BC-VW-mid (discrete_material ids renumbered 1..N: FEBio reads dmat as a list position)

## L6_ls.feb

Started from: `L5X_all.feb`

- run length only: time_steps 120 -> 100 (ends at t = 1); the model is unchanged
- solver: lstol 0 -> 0.9

## L8_fitall.feb

Started from: `L6_ls.feb`

- AVW-Para-L_fan only: material Beam-CL-USL (isotropic elastic, density=7.8e-07, E=21, v=0.3) -> AVW-Para-L_fan_fit: Ogden c1=8.73413, m1=4.75, k=1091.77 (250 mu0), density 7.8e-07; loft_fit_all.py strip fit to its 35 Abaqus connectors, rms 5.3 % over u = 2-47 mm
- AVW-Para-R_fan only: material Beam-CL-USL (isotropic elastic, density=7.8e-07, E=21, v=0.3) -> AVW-Para-R_fan_fit: Ogden c1=8.74029, m1=4.75, k=1092.54 (250 mu0), density 7.8e-07; loft_fit_all.py strip fit to its 35 Abaqus connectors, rms 5.3 % over u = 2-47 mm
- CL-L_fan only: material Beam-CL-USL (isotropic elastic, density=7.8e-07, E=21, v=0.3) -> CL-L_fan_fit: Ogden c1=0.291678, m1=6.5, k=36.4597 (250 mu0), density 7.8e-07; loft_fit_all.py strip fit to its 28 Abaqus connectors, rms 22.4 % over u = 2-47 mm
- CL-R_fan only: material Beam-CL-USL (isotropic elastic, density=7.8e-07, E=21, v=0.3) -> CL-R_fan_fit: Ogden c1=0.291797, m1=6.5, k=36.4747 (250 mu0), density 7.8e-07; loft_fit_all.py strip fit to its 28 Abaqus connectors, rms 22.4 % over u = 2-47 mm
- USL-L_fan only: material Beam-CL-USL (isotropic elastic, density=7.8e-07, E=21, v=0.3) -> USL-L_fan_fit: Ogden c1=0.0719942, m1=5, k=8.99927 (250 mu0), density 7.8e-07; loft_fit_all.py strip fit to its 13 Abaqus connectors, rms 24.3 % over u = 2-47 mm
- USL-R_fan only: material Beam-CL-USL (isotropic elastic, density=7.8e-07, E=21, v=0.3) -> USL-R_fan_fit: Ogden c1=0.0703702, m1=5, k=8.79628 (250 mu0), density 7.8e-07; loft_fit_all.py strip fit to its 13 Abaqus connectors, rms 24.3 % over u = 2-47 mm
- PM_fan only: material Beam-CL-USL (isotropic elastic, density=7.8e-07, E=21, v=0.3) -> PM_fan_fit: Ogden c1=42.7027, m1=4, k=5337.84 (250 mu0), density 7.8e-07; loft_fit_all.py strip fit to its 26 Abaqus connectors, rms 4.5 % over u = 2-47 mm
- PM_PeB_Left_fan only: material Beam-CL-USL (isotropic elastic, density=7.8e-07, E=21, v=0.3) -> PM_PeB_Left_fan_fit: Ogden c1=4.00225e-05, m1=3, k=0.00500281 (250 mu0), density 1.06e-09; loft_fit_all.py strip fit to its 8 Abaqus connectors, rms 8.1 % over u = 2-47 mm
- PM_PeB_Right_fan only: material Beam-CL-USL (isotropic elastic, density=7.8e-07, E=21, v=0.3) -> PM_PeB_Right_fan_fit: Ogden c1=3.98006e-05, m1=3, k=0.00497508 (250 mu0), density 1.06e-09; loft_fit_all.py strip fit to its 8 Abaqus connectors, rms 8.1 % over u = 2-47 mm
- PM_avw_bottom_left_fan only: material Beam-CL-USL (isotropic elastic, density=7.8e-07, E=21, v=0.3) -> PM_avw_bottom_left_fan_fit: Ogden c1=0.00013905, m1=3.5, k=0.0173813 (250 mu0), density 1.06e-09; loft_fit_all.py strip fit to its 2 Abaqus connectors, rms 4.8 % over u = 2-47 mm
- PM_avw_bottom_right_fan only: material Beam-CL-USL (isotropic elastic, density=7.8e-07, E=21, v=0.3) -> PM_avw_bottom_right_fan_fit: Ogden c1=0.000136846, m1=3.5, k=0.0171057 (250 mu0), density 1.06e-09; loft_fit_all.py strip fit to its 2 Abaqus connectors, rms 4.7 % over u = 2-47 mm
- PeB-constrin_fan only: material Beam-CL-USL (isotropic elastic, density=7.8e-07, E=21, v=0.3) -> PeB-constrin_fan_fit: Ogden c1=0.000681394, m1=4.75, k=0.0851742 (250 mu0), density 1.06e-09; loft_fit_all.py strip fit to its 26 Abaqus connectors, rms 4.5 % over u = 2-47 mm

## L9_fitall_edge.feb

Started from: `L8_fitall.feb`

- NOT IN SOURCE (loft convention, user): PM_PeB_Left_fan whole anchor edge under BC-RP-PM-PeB-origins: + 7 free edge nodes between the 8 fixed ones [14554, 14551, 14693, 14692, 14691, 14690, 14686] (off the straight segments between them: max 8.04e-08 mm)
- NOT IN SOURCE (loft convention, user): PM_PeB_Right_fan whole anchor edge under BC-RP-PM-PeB-origins: + 9 free edge nodes between the 8 fixed ones [14752, 14846, 14879, 14765, 14855, 14857, 14858, 14859, 14861] (off the straight segments between them: max 6.38e-08 mm)

## L11_stab.feb

Started from: `L9_fitall_edge.feb`

- removed 145 stab_* ground springs (discrete_material ids renumbered 1..N: FEBio reads dmat as a list position)

## L12_stabdamp_p3.feb

Started from: `L11_stab.feb`

- NOT IN SOURCE (numerical): mass damping settle_damping (C = 20/s) on from t = 0: load curve settle_on ['1.0,0', '1.05,1'] -> [0,1; 1.05,1] (it drives only the damping); Abaqus explicit has no such damping
- NOT IN SOURCE: Load-AVW pressure 0.014 -> 0.00466667 MPa (x0.3333); the other pressure loads unchanged
- NOT IN SOURCE: Load-PVW pressure 0.014 -> 0.00466667 MPa (x0.3333); the other pressure loads unchanged
- NOT IN SOURCE: Load-PeB-top pressure 0.014 -> 0.00466667 MPa (x0.3333); the other pressure loads unchanged
- NOT IN SOURCE: Load-top pressure 0.014 -> 0.00466667 MPa (x0.3333); the other pressure loads unchanged
