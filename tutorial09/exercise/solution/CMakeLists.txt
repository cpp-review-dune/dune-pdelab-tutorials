# dune_add_generated_executable(TARGET solution_navier_stokes
#                               UFLFILE navier_stokes.ufl
#                               INIFILE navier_stokes.ini
#                               )

dune_add_generated_executable(TARGET solution_nonlinear_poisson_dg
                              UFLFILE nonlinear_poisson_dg.ufl
                              INIFILE solution_nonlinear_poisson_dg.ini
                              )

# dune_symlink_to_source_files(FILES navier_stokes.ini nonlinear_poisson_dg.ini)

dune_symlink_to_source_files(FILES solution_nonlinear_poisson_dg.ini)
