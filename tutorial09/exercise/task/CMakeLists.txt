dune_add_generated_executable(TARGET nonlinear_poisson_dg
                              UFLFILE nonlinear_poisson_dg.ufl
                              INIFILE nonlinear_poisson_dg.ini
                              )

dune_symlink_to_source_files(FILES nonlinear_poisson_dg.ini)
