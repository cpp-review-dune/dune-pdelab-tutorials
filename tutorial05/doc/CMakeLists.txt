dune_symlink_to_source_files(FILES ../src/tutorial05.cc
../src/driver.hh
../src/tutorial05.ini
../src/nonlinearpoissonfemestimator.hh
../src/problem.hh
)

add_latex_document(tutorial05.tex
#                   BIBFILES tutorial05.bib
                   IMAGE_DIRS .
                   DEFAULT_PDF
                   MANGLE_TARGET_NAMES
                  )
