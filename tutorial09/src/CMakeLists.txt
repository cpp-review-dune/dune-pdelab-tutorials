dune_add_generated_executable(TARGET poisson
                              UFLFILE poisson.ufl
                              INIFILE poisson.ini
                              )

dune_add_generated_executable(TARGET nonlinear_poisson
                              UFLFILE nonlinear_poisson.ufl
                              INIFILE nonlinear_poisson.ini
                              )

dune_add_generated_executable(TARGET nonlinear_poisson_nitschebc
                              UFLFILE nonlinear_poisson_nitschebc.ufl
                              INIFILE nonlinear_poisson_nitschebc.ini
                              )

dune_symlink_to_source_files(FILES poisson.ini nonlinear_poisson.ini nonlinear_poisson_nitschebc.ini)
