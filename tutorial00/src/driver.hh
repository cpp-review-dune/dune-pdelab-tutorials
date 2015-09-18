/******************************************************/
//! \brief Driver function to set up and solve the problem
/******************************************************/
template<class GV>
int driver (GV& gv, Dune::ParameterTree& ptree)
{
  // dimension and important types
  const int dim = GV::dimension;
  typedef typename GV::Grid::ctype Coord; // type for ccordinates
  typedef double Real;                    // type for computations

  // make an object of all user defined classes
  typedef FFunction<GV,Real> F;
  F f(gv);
  BCType bct;
  typedef GFunction<GV,Real> G;
  G g(gv);

  // Make grid function space
  typedef Dune::PDELab::PkLocalFiniteElementMap<GV,Coord,Real,1> FEM;
  FEM fem(gv);
  typedef Dune::PDELab::ConformingDirichletConstraints CON; // constraints class
  typedef Dune::PDELab::ISTLVectorBackend<> VBE;
  typedef Dune::PDELab::GridFunctionSpace<GV,FEM,CON,VBE> GFS;
  GFS gfs(gv,fem);
  gfs.name("solution");
  typedef typename GFS::template ConstraintsContainer<Real>::Type CC;
  CC cc;
  Dune::PDELab::constraints(bct,gfs,cc); // assemble constraints
  std::cout << "constrained dofs=" << cc.size() << " of " << gfs.globalSize() << std::endl;

  // Initialize the solution vector
  using U = Dune::PDELab::Backend::Vector<GFS,Real>;
  U u(gfs); // initial value
  Dune::PDELab::interpolate(g,gfs,u);

  // Make a local operator
  typedef PoissonP1<G,F,FEM> LOP;
  LOP lop(g,f,fem.find(*gv.template begin<0>()));

  // Make a global operator
  typedef Dune::PDELab::istl::BCRSMatrixBackend<> MBE;
  MBE mbe(1<<(dim+1));
  typedef Dune::PDELab::GridOperator<
    GFS,GFS,        /* ansatz and test space */
    LOP,            /* local operator */
    MBE,            /* matrix backend */
    Real,Real,Real, /* field types for domain, range and jacobian */
    CC,CC           /* constraints transformation  for ansatz and test space */
    > GO;
  GO go(gfs,cc,gfs,cc,lop,mbe);

  // How well did we estimate the number of entries per matrix row?
  // => print Jacobian pattern statistics
  typename GO::Traits::Jacobian jac(go);
  std::cout << jac.patternStatistics() << std::endl;

  // Select a linear solver backend
  typedef Dune::PDELab::ISTLBackend_SEQ_CG_AMG_SSOR<GO> LS;
  LS ls(100,true);

  // Assemble and solve linear problem
  typedef Dune::PDELab::StationaryLinearProblemSolver<GO,LS,U> SLP;
  SLP slp(go,ls,u,1e-10);
  slp.apply();

  // Write VTK output file
  U uexact(gfs); // initial value
  Dune::PDELab::interpolate(g,gfs,uexact);
  Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTK::conforming);
  //  Dune::PDELab::addSolutionToVTKWriter(vtkwriter,gfs,uexact);
  typedef Dune::PDELab::DiscreteGridFunction<GFS,U> UDGF;
  UDGF udgf(gfs,u);
  vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<UDGF>(udgf,"computed solution"));
  UDGF uexactdgf(gfs,uexact);
  vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<UDGF>(uexactdgf,"interpolated solution"));
  vtkwriter.write(ptree.get("output.filename","output"),Dune::VTK::appendedraw);
}
