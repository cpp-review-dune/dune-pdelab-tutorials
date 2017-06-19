// -*- tab-width: 4; indent-tabs-mode: nil -*-
/** \file
    \brief Acoustic wave propagation in a simple 2D cavity

    Example 4 from "L. Krivodonova, J. Xin, J.-F. Remacle, N. Chevaugeon, J.E. Flaherty:
    Shock detection and limiting with discontinuous Galerkin methods for hyperbolic
    conservation laws. Applied Numerical Mathematics, 48, 323-338, 2004.
*/

//==============================================================================
// Parameter class for the linear acoustics problem
//==============================================================================

template<typename GV, typename RF>
class Krivodonova4Problem
{
public:
  typedef Dune::PDELab::LinearAcousticsParameterTraits<GV,RF> Traits;

  Krivodonova4Problem ()
    : pi(3.141592653589793238462643), time(0.0)
  {
  }

  //! speed of sound
  typename Traits::RangeFieldType
  c (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    return 340.0;
  }

  //! Dirichlet boundary condition value
  typename Traits::StateType
  g (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x, const typename Traits::StateType& s) const
  {
    typename Traits::DomainType xglobal = is.geometry().global(x);
    if (xglobal[0]<1e-6)
      {
        typename Traits::StateType u(0.0);
        u[1] = 1.224*(1+0.5*sin(2*pi*1500.0*time));
        return u;
      }
    if (xglobal[0]>1.0-1e-6)
      {
        typename Traits::StateType u(0.0);
        return u;
      }
    typename Traits::StateType u(0.0);
    u[2] = 0.0;
    return u;
  }

  //! right hand side
  typename Traits::StateType
  q (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    typename Traits::StateType rhs(0.0);
    return rhs;
  }

  //! initial value
  typename Traits::StateType
  u0 (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    typename Traits::StateType u(0.0);
    return u;
  }

  //! set time for subsequent evaluation
  void setTime (RF t)
  {
    time = t;
  }

private:
  double pi;
  RF time;
};


template<typename GV, typename RF>
class RiemannProblem
{
public:
  typedef Dune::PDELab::LinearAcousticsParameterTraits<GV,RF> Traits;

  RiemannProblem ()
    : time(0.0),  pi(3.141592653589793238462643)
  {
  }

  //! speed of sound
  typename Traits::RangeFieldType
  c (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    typename Traits::DomainType xglobal = e.geometry().global(x);
    if ( xglobal[1] < 1-(0.6/0.9)*(xglobal[0]-0.1) ) return 1.0;
    //    if (xglobal[0]>0.5) return 2.0;
    return 0.5;
  }

  //! Dirichlet boundary condition value
  typename Traits::StateType
  g (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x, const typename Traits::StateType& s) const
  {
    typename Traits::StateType u(0.0);
    u[0] = s[0];
    u[1] = -s[1];
    u[2] = -s[2];
    return u;
  }

  //! right hand side
  typename Traits::StateType
  q (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    typename Traits::StateType rhs(0.0);
    return rhs;
  }

  //! initial value
  typename Traits::StateType
  u0 (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    typename Traits::DomainType xglobal = e.geometry().global(x);
    typename Traits::StateType u(0.0);
    if (xglobal[0]>0.45 && xglobal[0]<0.55 && xglobal[1]>0.3 && xglobal[1]<0.4)
      {
        u[0] = sin(pi*(xglobal[0]-0.45)/0.1)*sin(pi*(xglobal[0]-0.45)/0.1)*sin(pi*(xglobal[1]-0.3)/0.1)*sin(pi*(xglobal[1]-0.3)/0.1);
        u[1] = 0;
        u[2] = 0;
      }
    return u;
  }

  //! set time for subsequent evaluation
  void setTime (RF t)
  {
    time = t;
  }

private:
  RF time;
  RF pi;
};


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

  //it woudl be better to extract it from fem
  //int degree = ptree.get("fem.degree",(int)1);


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

  // <<<2b>>> define problem parameters
  typedef RiemannProblem<GV,RF> Param;
  Param param;

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
  Dune::PDELab::LinearAcousticsInitialValueAdapter<Param> u0(gv,param);

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

