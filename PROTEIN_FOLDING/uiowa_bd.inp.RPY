---------------------------------------------------------------------------------
 teprint  ttprint  tmprint  num_lst_stp  num_fmd_stp  num_hyd_stp
   100.0   1000.0 1000000.0         100           10     100
---------------------------------------------------------------------------------
 num_threads  bond_dev_quit  i_do_lincs?  i_do_YHL?
    6           10.0           no         no
---------------------------------------------------------------------------------
 f_typs  f_mols  i_debug  q_desired  mol_Q1  mol_Q2  go_eps_low
      1       1       no      0.900       1       1        -1.0
---------------------------------------------------------------------------------
     xmin     xmax      ymin     ymax      zmin     zmax   i_pbc
-250.000 250.000 -250.000 250.000 -250.000 250.000       1
---------------------------------------------------------------------------------
 i_look_for_crashes?  periodic_bonds_okay?  i_use_high_mem?
                 yes                    no               no
---------------------------------------------------------------------------------
 i_append_movie  i_limit_verbosity?  i_use_v_typ  steepest_descent?
              0                 yes            1                 no
---------------------------------------------------------------------------------
 temperature  ionic_strength(mM)  r_ion  dielectric  viscosity(cP)  r_f_st
      298.15               150.0    0.0       78.40           0.89     0.0
---------------------------------------------------------------------------------
parameter_file    no_elec  wrap_molecules  RPY/none  scale_nb     BD/LD
parameter.file.01     yes               0       RPY       2.0  brownian
---------------------------------------------------------------------------------
go_potentials_file          i_use_go_pairs?  go_primacy  i_skip_intermol_go?
001.go.parameters                       yes           2                   no
---------------------------------------------------------------------------------
nomove_file   i_dont_move_some?
nomove_file                  no
---------------------------------------------------------------------------------
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
001.charge.parameters             ! charges/hydrodynamic radii in this file
001.internal.parameters           ! bonds/angles/dihedrals
1                                 ! # of copies of this molecule type
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
time_step  totsimtime  vdw_s  vdw_m  goe_s  goe_m  ele_s  ele_m  f_f_cell_size
    0.025  10000000.0  25.00  25.00  25.00  25.00   0.00   0.00         100.00
---------------------------------------------------------------------------------
position_restraint_file  i_do_pos_restraints?  mission_creep?
position_restraint_file                    no              no
---------------------------------------------------------------------------------
r_size   l_size   r_size_fac   n_size   f_size
  -1.0     -1.0          1.0       -1     0.10
--------------------------------------------------------------------------------
fixman?  fixman_tol  fixman_order  fixman_override?    lmin    lmax
     no       0.010           100                no  0.0010    10.0
--------------------------------------------------------------------------------
treecode?   theta   order   shrink   maxatm
       no   0.100       4        1       25
--------------------------------------------------------------------------------
walls?   num_walls   wall_file
    no           0   wall_file
--------------------------------------------------------------------------------
