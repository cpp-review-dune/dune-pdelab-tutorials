/********************************************************/
// Beware of line number changes, they may corrupt docu!
//! \brief Driver function to set up and solve the problem
/********************************************************/

template<typename GV, typename FEM>
void driver (const GV& gv, const FEM& fem, Dune::ParameterTree& ptree)
{
  // dimension and important types
  const int dim = GV::dimension;
  typedef typename GV::Grid::ctype DF; // type for ccordinates
  typedef double RF;                   // type for computations

  // make user functions
  RF eta = ptree.get("problem.eta",(RF)1.0);
  Problem<RF> problem(eta);
  auto glambda = [&](const auto& e, const auto& x){return problem.g(e,x);};
  auto g = Dune::PDELab::makeGridFunctionFromLocalLambda(gv,glambda);
  auto blambda = [&](const auto& i, const auto& x){return problem.b(i,x);};
  auto b = Dune::PDELab::makeBoundaryConditionFromLocalLambda(blambda);

  // Make grid function space
  typedef Dune::PDELab::ConformingDirichletConstraints CON;
  typedef Dune::PDELab::istl::VectorBackend<> VBE;
  typedef Dune::PDELab::GridFunctionSpace<GV,FEM,CON,VBE> GFS;
  GFS gfs(gv,fem);
  gfs.name("Vh");

  // Assemble constraints
  typedef typename GFS::template
    ConstraintsContainer<RF>::Type CC;
  CC cc;
  Dune::PDELab::constraints(b,gfs,cc); // assemble constraints
  std::cout << "constrained dofs=" << cc.size() << " of "
            << gfs.globalSize() << std::endl;

  // A coefficient vector
  using Z = Dune::PDELab::Backend::Vector<GFS,RF>;
  Z z(gfs); // initial value

  // Make a grid function out of it
  typedef Dune::PDELab::DiscreteGridFunction<GFS,Z> ZDGF;
  ZDGF zdgf(gfs,z);

  // Fill the coefficient vector
  Dune::PDELab::interpolate(g,gfs,z);

  // Make a local operator
  typedef NonlinearPoissonFEM<Problem<RF>,FEM> LOP;
  LOP lop(problem);

  // Make a global operator
  typedef Dune::PDELab::istl::BCRSMatrixBackend<> MBE;
  int degree = ptree.get("fem.degree",(int)1);
  MBE mbe((int)pow(1+2*degree,dim));
  typedef Dune::PDELab::GridOperator<
    GFS,GFS,  /* ansatz and test space */
    LOP,      /* local operator */
    MBE,      /* matrix backend */
    RF,RF,RF, /* domain, range, jacobian field type*/
    CC,CC     /* constraints for ansatz and test space */
    > GO;
  GO go(gfs,cc,gfs,cc,lop,mbe);

  // Select a linear solver backend
  typedef Dune::PDELab::ISTLBackend_SEQ_CG_AMG_SSOR<GO> LS;
  LS ls(100,0);

  // solve nonlinear problem
  Dune::PDELab::Newton<GO,LS,Z> newton(go,z,ls);
  newton.setReassembleThreshold(0.0);
  newton.setVerbosityLevel(2);
  newton.setReduction(1e-10);
  newton.setMinLinearReduction(1e-4);
  newton.setMaxIterations(25);
  newton.setLineSearchMaxIterations(10);
  newton.apply();

  // Write VTK output file
  Dune::SubsamplingVTKWriter<GV> vtkwriter(gv,ptree.get("output.subsampling",(int)0));
  typedef Dune::PDELab::VTKGridFunctionAdapter<ZDGF> VTKF;
  vtkwriter.addVertexData(std::shared_ptr<VTKF>(new
                                         VTKF(zdgf,"fesol")));
  vtkwriter.write(ptree.get("output.filename","output"),
                  Dune::VTK::appendedraw);
}
