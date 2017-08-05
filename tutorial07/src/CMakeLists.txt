add_executable("tutorial07-maxwell" tutorial07-maxwell.cc)
add_executable("tutorial07-acoustics" tutorial07-acoustics.cc)

dune_symlink_to_source_files(FILES tutorial07-maxwell.ini)
dune_symlink_to_source_files(FILES tutorial07-acoustics.ini)
