dune_symlink_to_source_files(FILES ../src/tutorial00.cc
                                   ../src/tutorial00.ini
                                   ../src/bctype.hh
                                   ../src/driver.hh
                                   ../src/poissonp1.hh
                                   ../src/gfunction.hh
                                   ../src/ffunction.hh)

add_latex_document(tutorial00.tex TARGET_NAME tutorial00_pdf
                   BIBFILES tutorial00.bib
                   IMAGE_DIRS .
                   FORCE_PDF
                  )

