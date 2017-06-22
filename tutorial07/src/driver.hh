// -*- tab-width: 4; indent-tabs-mode: nil -*-

//===============================================================
// driver for general pouropse hyperbolic solver
//===============================================================

template<typename GV, typename FEMDG>
void driver (const GV& gv, const FEMDG& femdg, Dune::ParameterTree& ptree)
{
  //std::cout << "using degree " << degree << std::endl;

  // Choose domain and range field type
  using RF = double;                   // type for computations
  const int dim = GV::dimension;

  //Minak it woudl be better to extract it from fem
  //int degree = ptree.get("fem.degree",(int)1);

  // make PDE parameter class

  typedef RiemannProblem<GV,RF> Param;
  Param param;

  //minak: no need for dirichlet BC
  //auto glambda = [&](const auto& e, const auto& x)
  //  {return param.g(e,x);};
  //auto g = Dune::PDELab::
  //  makeBaundaryConditionsFromCallable(gv,glambda);

  //initial condition
  auto u0lambda = [&](const auto& i, const auto& x)
    {return param.u0(i,x);};
  auto u0 = Dune::PDELab::
    makeGridFunctionFromCallable(gv,u0lambda);


  //TODO figure out proper blocksize
  // <<<2>>> Make grid function space
  //const int blocksize = Dune::PB::PkSize<degree,dim>::value;


  typedef Dune::PDELab::NoConstraints CON;
  //typedef Dune::PDELab::ISTLVectorBackend
  //  <Dune::PDELab::istl::Parameters::static_blocking,blocksize> VBE;

  // blocking fixed is not allowed???
  using VBE0 = Dune::PDELab::istl::VectorBackend<>;

  using VBE = Dune::PDELab::istl::VectorBackend<Dune::PDELab::istl::Blocking::fixed>;
  using OrderingTag = Dune::PDELab::EntityBlockedOrderingTag;
  typedef Dune::PDELab::GridFunctionSpace<GV,FEMDG,CON,VBE0> GFSDG;
  GFSDG gfsdg(gv,femdg);

  //TODO use vector grid function space
  //typedef Dune::PDELab::PowerGridFunctionSpace
  //  <GFSDG,dim+1,Dune::PDELab::ISTLVectorBackend<> > GFS;

  typedef Dune::PDELab::PowerGridFunctionSpace<GFSDG,dim+1,VBE,OrderingTag> GFS;
  GFS gfs(gfsdg);

  /* for the future
  typedef Dune::PDELab::VectorGridFunctionSpace
    <GV,FEMDG,dim,VBE0,VBE,CON> GFS;
  GFS gfs(gv,femdg);
  gfs.name("u");
  */

  // Add names to the components for VTK output
  using namespace Dune::TypeTree::Indices;
  gfs.child(_0).name("u0");
  gfs.child(_1).name("u1");
  gfs.child(_2).name("u2");
  

  typedef typename GFS::template ConstraintsContainer<RF>::Type C;
  C cg;
  gfs.update(); // initializing the gfs
  std::cout << "degrees of freedom: " << gfs.globalSize() << std::endl;



  // Make instationary grid operator
  using LOP = Dune::PDELab::DGLinearAcousticsSpatialOperator<Param,FEMDG>;
  LOP lop(param);
  using TLOP = Dune::PDELab::DGLinearAcousticsTemporalOperator<Param,FEMDG>;
  TLOP tlop(param);

  using MBE = Dune::PDELab::istl::BCRSMatrixBackend<>;
  MBE mbe(5); // Maximal number of nonzeroes per row can be cross-checked by patternStatistics().

  using GO0 = Dune::PDELab::GridOperator<GFS,GFS,LOP,MBE,RF,RF,RF,C,C>;
  GO0 go0(gfs,cg,gfs,cg,lop,mbe);
  using GO1 = Dune::PDELab::GridOperator<GFS,GFS,TLOP,MBE,RF,RF,RF,C,C>;
  GO1 go1(gfs,cg,gfs,cg,tlop,mbe);
  using IGO = Dune::PDELab::OneStepGridOperator<GO0,GO1,false>;
  IGO igo(go0,go1);


  // select and prepare time-stepping scheme
  int torder = ptree.get("fem.torder",(int)1);

  Dune::PDELab::ExplicitEulerParameter<RF> method1;
  Dune::PDELab::HeunParameter<RF> method2;
  Dune::PDELab::Shu3Parameter<RF> method3;
  Dune::PDELab::RK4Parameter<RF> method4;
  Dune::PDELab::TimeSteppingParameterInterface<RF> *method;

  if (torder==0) {method=&method1; std::cout << "setting explicit Euler" << std::endl;}
  if (torder==1) {method=&method2; std::cout << "setting Heun" << std::endl;}
  if (torder==2) {method=&method3; std::cout << "setting Shu 3" << std::endl;}
  if (torder==3) {method=&method4; std::cout << "setting RK4" << std::endl;}
  if (torder<1||torder>3) std::cout<<"torder should be in [1,3]"<<std::endl;


  igo.setMethod(*method);

  // <<<5>>> set initial values
  typedef typename IGO::Traits::Domain V;
  V xold(gfs,0.0);
  //Dune::PDELab::LinearAcousticsInitialValueAdapter<Param> u0(gv,param);

  Dune::PDELab::interpolate(u0,gfs,xold);
 

  //only interpolate components TODO interpolation of vector function space


  // <<<6>>> Make a linear solver backend
  //typedef Dune::PDELab::ISTLBackend_SEQ_CG_SSOR LS;
  //LS ls(10000,1);
  typedef Dune::PDELab::ISTLBackend_OVLP_ExplicitDiagonal<GFS> LS;
  LS ls(gfs);

  // <<<8>>> time-stepper
  typedef Dune::PDELab::CFLTimeController<RF,IGO> TC;
  TC tc(0.999,igo);
  Dune::PDELab::ExplicitOneStepMethod<RF,IGO,LS,V,V,TC> osm(*method,igo,ls,tc);
  osm.setVerbosityLevel(2);

  // prepare VTK writer and write first file
  int subsampling=ptree.get("output.subsampling",(int)0);
  using VTKWRITER = Dune::SubsamplingVTKWriter<GV>;
  VTKWRITER vtkwriter(gv,subsampling);

  std::string filename=ptree.get("output.filename","output"); //+ STRINGIZE(GRIDDIM) "d";

  struct stat st;
  if( stat( filename.c_str(), &st ) != 0 )
    {
      int stat = 0;
      stat = mkdir(filename.c_str(),S_IRWXU|S_IRWXG|S_IRWXO);
      if( stat != 0 && stat != -1)
        std::cout << "Error: Cannot create directory "
                  << filename << std::endl;
    }
  using VTKSEQUENCEWRITER = Dune::VTKSequenceWriter<GV>;
  VTKSEQUENCEWRITER vtkSequenceWriter(
    std::make_shared<VTKWRITER>(vtkwriter),filename,filename,"");
  // add data field for all components of the space to the VTK writer
  Dune::PDELab::addSolutionToVTKWriter(vtkSequenceWriter,gfs,xold);
  vtkSequenceWriter.write(0.0,Dune::VTK::appendedraw);

  // initialize simulation time
  RF time = 0.0;

  // time loop
  RF T = ptree.get("problem.T",(RF)1.0);
  RF dt = ptree.get("fem.dt",(RF)0.1);
  const int every = ptree.get<int>("output.every");

  V x(gfs,0.0);

  while (time < T)
    {
      // do time step
      osm.apply(time,dt,xold,x);

      // TODO output to VTK file every n-th timestep 
      //if(osm.result().total.timesteps % every == 0)
        vtkSequenceWriter.write(time,Dune::VTK::appendedraw);

      xold = x;
      time += dt;
    }
}

