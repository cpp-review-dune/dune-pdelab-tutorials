add_executable("navier-stokes-cylinder" navier-stokes-cylinder.cc)

dune_symlink_to_source_files(FILES navier-stokes-cylinder.ini)
dune_symlink_to_source_files(FILES cylinder2d.geo)
dune_symlink_to_source_files(FILES cylinder2dmesh1.msh)
dune_symlink_to_source_files(FILES cylinder2dmesh2.msh)
dune_symlink_to_source_files(FILES cylinder2dmesh3.msh)
dune_symlink_to_source_files(FILES cylinder2dmesh4.msh)

