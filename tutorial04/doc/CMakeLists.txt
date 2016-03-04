dune_symlink_to_source_files(FILES ../src/tutorial04.cc
../src/driver.hh
../src/wavefem.hh
../src/problem.hh
../src/tutorial04.ini
)

add_latex_document(tutorial04.tex
                   BIBFILES tutorial04.bib
                   IMAGE_DIRS .
                   DEFAULT_PDF
                   MANGLE_TARGET_NAMES
                  )
