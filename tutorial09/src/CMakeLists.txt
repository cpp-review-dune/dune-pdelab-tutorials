dune_add_formcompiler_system_test(UFLFILE poisson.ufl
                                  BASENAME poisson
                                  INIFILE poisson.mini
                                  NO_TESTS
                                  )

dune_add_formcompiler_system_test(UFLFILE nonlinear_poisson.ufl
                                  BASENAME nonlinear_poisson
                                  INIFILE nonlinear_poisson.mini
                                  NO_TESTS
                                  )

dune_add_formcompiler_system_test(UFLFILE nonlinear_poisson_nietschebc.ufl
                                  BASENAME nonlinear_poisson_nietschebc
                                  INIFILE nonlinear_poisson_nietschebc.mini
                                  NO_TESTS
                                  )
