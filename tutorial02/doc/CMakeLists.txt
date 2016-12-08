dune_symlink_to_source_files(FILES ../src/driver.hh
                                   ../src/nonlinearpoissonfv.hh
                                   ../src/problem.hh
                                   ../src/tutorial02.cc
                                   ../src/tutorial02.ini)
add_latex_document(tutorial02.tex
                   IMAGE_DIRS .
                   FORCE_PDF
                  )
