dune_symlink_to_source_files(FILES ../src/tutorial02.cc
                             ../src/driver.hh ../src/tutorial02.ini
                             ../src/nonlinearpoissonfv.hh)

add_latex_document(tutorial02.tex
                   IMAGE_DIRS .
                   DEFAULT_PDF
                   MANGLE_TARGET_NAMES
                  )
