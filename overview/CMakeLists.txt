add_latex_document(abstractions.tex
                   BIBFILES abstractions.bib
                   DEFAULT_PDF
                   MANGLE_TARGET_NAMES
                  )

add_latex_document(overview.tex
                   DEFAULT_PDF
                   MANGLE_TARGET_NAMES
                  )

add_latex_document(exercise-workflow.tex
                   DEFAULT_PDF
                   MANGLE_TARGET_NAMES
                  )
dune_symlink_to_source_files(FILES debug.opts release.opts)
