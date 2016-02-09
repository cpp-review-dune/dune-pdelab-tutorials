dune_symlink_to_source_files(FILES ../src/tutorial06.cc)

add_latex_document(tutorial06.tex
                   BIBFILES tutorial06.bib
                   IMAGE_DIRS .
                   DEFAULT_PDF
                   MANGLE_TARGET_NAMES
                  )
