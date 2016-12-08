add_latex_document(abstractions.tex TARGET_NAME abstractions_pdf
                   BIBFILES abstractions.bib
                   FORCE_PDF
                  )

add_latex_document(overview.tex TARGET_NAME overview_pdf
                   FORCE_PDF
                  )

add_latex_document(exercise-workflow.tex TARGET_NAME exercise-workflow_pdf
                   FORCE_PDF
                  )
dune_symlink_to_source_files(FILES debug.opts release.opts)
