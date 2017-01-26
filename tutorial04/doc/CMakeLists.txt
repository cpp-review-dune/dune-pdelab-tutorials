dune_symlink_to_source_files(FILES ../src/tutorial04.cc
../src/driver.hh
../src/wavefem.hh
../src/tutorial04.ini
)

add_latex_document(tutorial04.tex TARGET_NAME tutorial04_pdf
                   BIBFILES tutorial04.bib
                   IMAGE_DIRS .
                   FORCE_PDF
                  )
