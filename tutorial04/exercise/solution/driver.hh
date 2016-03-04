/********************************************************/
// Beware of line number changes, they may corrupt docu!
//! \brief Driver function to set up and solve the problem
/********************************************************/

template<typename GV, typename FEM>
void driver (const GV& gv, const FEM& fem, Dune::ParameterTree& ptree)
{
  // dimension and important types
  const int dim = GV::dimension;
  typedef double RF;                   // type for computations

  // Make grid function space used per component
  typedef Dune::PDELab::ConformingDirichletConstraints CON;
  typedef Dune::PDELab::istl::VectorBackend<> VBE0;
  typedef Dune::PDELab::GridFunctionSpace<GV,FEM,CON,VBE0> GFS0;
  GFS0 gfs0(gv,fem);

  // Make grid function space for the system
  typedef Dune::PDELab::istl::VectorBackend<Dune::PDELab::istl::Blocking::fixed> VBE;
  typedef Dune::PDELab::EntityBlockedOrderingTag OrderingTag;
  typedef Dune::PDELab::PowerGridFunctionSpace<GFS0,2,VBE,OrderingTag> GFS;
  GFS gfs(gfs0);

  // subspaces
  typedef Dune::PDELab::GridFunctionSubSpace
    <GFS,Dune::TypeTree::TreePath<0> > U0SUB;
  U0SUB u0sub(gfs);
  typedef Dune::PDELab::GridFunctionSubSpace
    <GFS,Dune::TypeTree::TreePath<1> > U1SUB;
  U1SUB u1sub(gfs);

  // define the initial condition

  auto ulambda = [dim](const auto& x){
    Dune::FieldVector<RF,2> rv(0.0);
    for (int i=0; i<dim; i++) rv[0] += (x[i]-0.375)*(x[i]-0.375);
    rv[0] = std::max(0.0,1.0-8.0*sqrt(rv[0]));
    return rv;
  };
  /*
  //minak: define stationary sine initial condition
  auto ulambda = [dim](const auto& x){
    Dune::FieldVector<RF,2> rv(0.0);
    for (int i=0; i<dim; i++) rv[0] += sin(2*x[i]);
    return rv;
  };
  */
  auto u = Dune::PDELab::makeGridFunctionFromCallable(gv,ulambda);

  // set up coefficient vector filled with initial condition
  using Z = Dune::PDELab::Backend::Vector<GFS,RF>;
  Z z(gfs); // initial value
  Dune::PDELab::interpolate(u,gfs,z);

  // Make discrete grid functions for components
  typedef Dune::PDELab::DiscreteGridFunction<U0SUB,Z> U0DGF;
  U0DGF u0dgf(u0sub,z);
  typedef Dune::PDELab::DiscreteGridFunction<U1SUB,Z> U1DGF;
  U1DGF u1dgf(u1sub,z);

  // assemble constraints
  auto b0lambda = [](const auto& x){return true;};
  auto b0 = Dune::PDELab::
    makeBoundaryConditionFromCallable(gv,b0lambda);
  auto b1lambda = [](const auto& x){return true;};
  auto b1 = Dune::PDELab::
    makeBoundaryConditionFromCallable(gv,b1lambda);
  typedef Dune::PDELab::CompositeConstraintsParameters<decltype(b0),decltype(b1)> B;
  B b(b0,b1);
  typedef typename GFS::template ConstraintsContainer<RF>::Type CC;
  CC cc;
  Dune::PDELab::constraints(b,gfs,cc);
  std::cout << "constrained dofs=" << cc.size() << " of "
            << gfs.globalSize() << std::endl;
  set_constrained_dofs(cc,0.0,z); // set zero Dirichlet boundary conditions

  // prepare VTK writer and write first file
  int subsampling=ptree.get("output.subsampling",(int)0);
  typedef Dune::SubsamplingVTKWriter<GV> VTKWRITER;
  VTKWRITER vtkwriter(gv,subsampling);
  //std::string filename=ptree.get("output.filename","output");
  std::string filename=ptree.get("output.filename","output") +ptree.get("grid.dim","grid") + "d";

  struct stat st;
  if( stat( filename.c_str(), &st ) != 0 )
    {
      int stat = 0;
      stat = mkdir( filename.c_str(), S_IRWXU | S_IRWXG | S_IRWXO );
      if( stat != 0 && stat != -1)
        std::cout << "Error: Cannot create directory "
                  << filename << std::endl;
    }
  typedef Dune::VTKSequenceWriter<GV> VTKSEQUENCEWRITER;
  VTKSEQUENCEWRITER vtkSequenceWriter(std::make_shared<VTKWRITER>(vtkwriter),filename,filename,"");
  typedef Dune::PDELab::VTKGridFunctionAdapter<U0DGF> VTKF0;
  vtkSequenceWriter.addVertexData(std::shared_ptr<VTKF0>(new VTKF0(u0dgf,"u0")));
  typedef Dune::PDELab::VTKGridFunctionAdapter<U1DGF> VTKF1;
  vtkSequenceWriter.addVertexData(std::shared_ptr<VTKF1>(new VTKF1(u1dgf,"u1")));
  vtkSequenceWriter.write(0.0,Dune::VTK::appendedraw);
  



  // Make instationary grid operator
  double speedofsound=ptree.get("problem.speedofsound",(double)1.0);
  typedef WaveFEM<FEM> LOP;
  //typedef WaveFEMElip<FEM> LOP;
  LOP lop(speedofsound);
  typedef WaveL2<FEM> TLOP;
  //typedef WaveElip<FEM> TLOP;
  TLOP tlop;
  typedef Dune::PDELab::istl::BCRSMatrixBackend<> MBE;
  int degree = ptree.get("fem.degree",(int)1);
  MBE mbe((int)pow(1+2*degree,dim));
  typedef Dune::PDELab::GridOperator<GFS,GFS,LOP,MBE,RF,RF,RF,CC,CC> GO0;
  GO0 go0(gfs,cc,gfs,cc,lop,mbe);
  typedef Dune::PDELab::GridOperator<GFS,GFS,TLOP,MBE,RF,RF,RF,CC,CC> GO1;
  GO1 go1(gfs,cc,gfs,cc,tlop,mbe);
  typedef Dune::PDELab::OneStepGridOperator<GO0,GO1> IGO;
  IGO igo(go0,go1);

  // Linear problem solver
  typedef Dune::PDELab::ISTLBackend_SEQ_BCGS_SSOR LS;
  LS ls(5000,false);
  typedef Dune::PDELab::
    StationaryLinearProblemSolver<IGO,LS,Z> SLP;
  SLP slp(igo,ls,1e-8);

  // select and prepare time-stepping scheme
  int torder = ptree.get("fem.torder",(int)1);
  Dune::PDELab::OneStepThetaParameter<RF> method1(0.5);
  Dune::PDELab::Alexander2Parameter<RF> method2;
  Dune::PDELab::Alexander3Parameter<RF> method3;
  Dune::PDELab::TimeSteppingParameterInterface<RF>* pmethod=&method1;
  if (torder==1) pmethod = &method1;
  if (torder==2) pmethod = &method2;
  if (torder==3) pmethod = &method3;
  if (torder<1||torder>3) std::cout<<"torder should be in [1,3]"<<std::endl;
  Dune::PDELab::OneStepMethod<RF,IGO,SLP,Z,Z> osm(*pmethod,igo,slp);
  osm.setVerbosityLevel(2);

  // initialize simulation time
  RF time = 0.0;

  // time loop
  RF T = ptree.get("problem.T",(RF)1.0);
  RF dt = ptree.get("fem.dt",(RF)0.1);
  while (time<T-1e-8)
    {
      // do time step
      Z znew(z);
      osm.apply(time,dt,z,znew);

      // accept time step
      z = znew;
      time+=dt;

      // Energy - integrate u square
	  typedef Dune::PDELab::SqrGridFunctionAdapter<U1DGF> SQRU1DGF;
	  SQRU1DGF u1sqrdgf(u1dgf);

	  typename SQRU1DGF::Traits::RangeType sqru1(0.0);
      Dune::PDELab::integrateGridFunction(u1sqrdgf, sqru1, 2*degree);

	  //typedef Dune::PDELab::GridFunctionSubSpace<U0DGF,Dune::TypeTree::TreePath<0> > SUBU0;
      typedef Dune::PDELab::DiscreteGridFunctionGradient<U0SUB, Z> DGFU0G;
      DGFU0G dgfu0g(u0sub,z);
      //gradient square
      typedef Dune::PDELab::SqrGridFunctionAdapter<DGFU0G> SQRDGFU0G;
	  SQRDGFU0G sqrdgfu0g(dgfu0g);
      //integrate gradient square
      typename SQRDGFU0G::Traits::RangeType sqru0g(0.0);
      Dune::PDELab::integrateGridFunction(sqrdgfu0g, sqru0g, 2*degree);

      std::cout << "::: energy          " << sqru1 + sqru0g << std::endl;

      // output to VTK file
      vtkSequenceWriter.write(time,Dune::VTK::appendedraw);
    }
}
