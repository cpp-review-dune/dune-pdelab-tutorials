int main(int argc, char** argv) {
  try {
    //Maybe initialize Mpi
    Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);
    if(Dune::MPIHelper::isFake)
      std::cout<< "This is a sequential program." << '\n';
    else {
      if(helper.rank()==0)
        std::cout << "parallel run on " << helper.size() << " process(es)" << '\n';
    }

    if (argc!=2) {
      if(helper.rank()==0)
        std::cout << "usage: ./example01 <level>" << '\n';
      return 1;
    }

    int level;
    sscanf(argv[1],"%d",&level);

    // sequential version
    if (helper.size()==1) {
      Dune::FieldVector<double,2> L(1.0);
      auto N = Dune::filledArray<2, int>(1);
      std::bitset<2> periodic(false);
      int overlap=0;
      Dune::YaspGrid<2> grid(L,N,periodic,overlap);
      grid.globalRefine(level);
      typedef Dune::YaspGrid<2>::LeafGridView GV;
      const GV& gv=grid.leafGridView();
      example01a_Qk<1>(gv); // Q1
      example01a_Qk<2>(gv); // Q2
      example01a_RT(gv);
      example01b_Q2(gv);
    }
  }
  catch (Dune::Exception &e) {
    std::cerr << "Dune reported error: " << e << '\n';
    return 1;
  }
  catch (...) {
    std::cerr << "Unknown exception thrown!" << '\n';
    return 1;
  }
}
