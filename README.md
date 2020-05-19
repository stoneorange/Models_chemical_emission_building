# Models_chemical_emission_building

This repository contains matlabs codes for the different model solutions for chemical emissions from building materials.

The model deals with a single diffuisonal source of the chemical with a single layer and one convective surface.

Model solutions include:
1) Analytical solution. Files: compute_qn_vector.m, compute_An.m, compute_y.m, compute_c.m, compute_me.m
2) Numercial solutions.
  2a) Huang's even method. File: MOL_single_even_par_modal.m
  2b) Yan's even method. File: SS_single_even_par_modal.m
  2c) Guo's uneven method. File: SS_single_uneven_Guo_par_modal.m
  2d) Huang's uneven method. File: MOL_single_uneven_q_2_par_modal.m
3) Simplified solutions.
  3a) D-limited model. Files: sim_D_limited.m, compute_q1_proxy_1b.m
  3b) SVOC constant concentration model. File: Sim_Little_model_no_sorption_par.m
  3c) K-limited model. File: Sim_K_limited_no_sorption_par.m
