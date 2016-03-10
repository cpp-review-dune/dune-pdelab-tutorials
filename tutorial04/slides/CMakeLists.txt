file(GLOB_RECURSE _slides04_tex_inputs /tikz/*.tex)
add_latex_document(slides04.tex
                   IMAGE_DIRS .
                   INPUTS ${_slides04_tex_inputs}
                   BIBFILES slides04.bib
                   DEFAULT_PDF
                   MANGLE_TARGET_NAMES
                  )
