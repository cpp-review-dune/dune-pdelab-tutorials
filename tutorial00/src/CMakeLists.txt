add_executable(tutorial00 tutorial00.cc)
add_dune_alberta_flags(tutorial00)

set(CLEANFILES  dgfparser.log q1*.vtu *.vtu *.png)
dune_symlink_to_source_files(FILES tutorial00.ini unitsquare.geo unitsquare.msh unitcube.geo unitcube.msh)
