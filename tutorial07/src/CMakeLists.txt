foreach(dim 2 3)
  add_executable("tutorial07-${dim}d" tutorial07.cc)
  target_compile_definitions("tutorial07-${dim}d" PUBLIC "GRIDDIM=${dim}")
endforeach()

dune_symlink_to_source_files(FILES tutorial07.ini)
