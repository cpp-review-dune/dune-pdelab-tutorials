dune_symlink_to_source_files(FILES ../src/tutorial01.cc
                             ../src/driver.hh
                             ../src/problem.hh
                             ../src/nonlinearpoissonfem.hh
                             ../src/tutorial01.ini )

add_latex_document(tutorial01.tex TARGET_NAME tutorial01_pdf
                   BIBFILES tutorial01.bib
                   IMAGE_DIRS .
                   FORCE_PDF
                  )
