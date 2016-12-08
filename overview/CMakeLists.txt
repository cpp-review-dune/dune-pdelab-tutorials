add_latex_document(abstractions.tex
                   BIBFILES abstractions.bib
                   FORCE_PDF
                  )

add_latex_document(overview.tex
                   FORCE_PDF
                  )

add_latex_document(exercise-workflow.tex
                   FORCE_PDF
                  )
dune_symlink_to_source_files(FILES debug.opts release.opts)
