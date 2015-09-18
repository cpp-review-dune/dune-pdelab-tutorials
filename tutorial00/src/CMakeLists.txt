add_executable(tutorial00 tutorial00.cc)

dune_symlink_to_source_files(FILES tutorial00.ini unitsquare.geo unitsquare.msh unitcube.geo unitcube.msh)
