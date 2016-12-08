dune_symlink_to_source_files(FILES ../src/tutorial03.cc
../src/driver.hh
../src/nonlinearheatfem.hh
../src/problem.hh
../src/tutorial03.ini
)

add_latex_document(tutorial03.tex
                   BIBFILES tutorial03.bib
                   IMAGE_DIRS .
                   FORCE_PDF
                  )
