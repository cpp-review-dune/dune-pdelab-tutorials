dune_symlink_to_source_files(FILES ../src/tutorial00.cc
                                   ../src/tutorial00.ini
                                   ../src/driver.hh
                                   ../src/poissonp1.hh)

add_latex_document(tutorial00.tex
                   BIBFILES tutorial00.bib
                   IMAGE_DIRS .
                   DEFAULT_PDF
                   MANGLE_TARGET_NAMES
                  )

