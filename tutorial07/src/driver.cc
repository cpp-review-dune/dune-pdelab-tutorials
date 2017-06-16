// -*- tab-width: 4; indent-tabs-mode: nil -*-
/** \file
    \brief Acoustic wave propagation in a simple 2D cavity

    Example 4 from "L. Krivodonova, J. Xin, J.-F. Remacle, N. Chevaugeon, J.E. Flaherty:
    Shock detection and limiting with discontinuous Galerkin methods for hyperbolic
    conservation laws. Applied Numerical Mathematics, 48, 323-338, 2004.
*/
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include<math.h>
#include<iostream>
#include<vector>
#include<map>
#include<string>

#include<sys/stat.h>

#include<dune/common/parallel/mpihelper.hh>
#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>
#include<dune/common/typetraits.hh>
#include<dune/common/timer.hh>
#include<dune/grid/io/file/vtk.hh> //added 
#include<dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include<dune/grid/io/file/gmshreader.hh>
#include<dune/grid/yaspgrid.hh>
#if HAVE_UG
#include<dune/grid/uggrid.hh>
#endif
#include<dune/istl/bvector.hh>
#include<dune/istl/operators.hh>
#include<dune/istl/solvers.hh>
#include<dune/istl/preconditioners.hh>
#include<dune/istl/io.hh>
#include<dune/istl/superlu.hh>

#include<dune/pdelab/finiteelementmap/opbfem.hh>

#include<dune/pdelab/gridfunctionspace/subspace.hh>
#include <dune/pdelab/gridfunctionspace/vectorgridfunctionspace.hh> //added by minak
#include<dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include<dune/pdelab/gridfunctionspace/genericdatahandle.hh>
#include<dune/pdelab/gridfunctionspace/interpolate.hh>
#include<dune/pdelab/gridfunctionspace/vtk.hh> //added by minak
#include<dune/pdelab/constraints/common/constraints.hh>
#include<dune/pdelab/common/function.hh>
#include<dune/pdelab/common/vtkexport.hh>
#include<dune/pdelab/gridoperator/gridoperator.hh>
#include<dune/pdelab/gridoperator/onestep.hh>
#include<dune/pdelab/backend/istl.hh>
#include<dune/pdelab/instationary/onestep.hh>
#include<dune/pdelab/common/instationaryfilenamehelper.hh>

#include<dune/pdelab/gridoperator/onestep.hh>
#include<dune/pdelab/instationary/onestep.hh>
#include<dune/pdelab/gridoperator/gridoperator.hh>

#include"linearacousticsdg.hh"


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
// Some variants to solve the nonlinear diffusion problem
//===============================================================


// example using explicit time-stepping
template<class GV, class FEMDG, int degree>
void explicit_scheme (const GV& gv, const FEMDG& femdg, double Tend, double timestep, std::string name, int modulo)
{
  std::cout << "using degree " << degree << std::endl;
  // <<<1>>> Choose domain and range field type
  typedef double Real;
  const int dim = GV::dimension;



  // <<<2>>> Make grid function space
  const int blocksize = Dune::PB::PkSize<degree,dim>::value;
  typedef Dune::PDELab::NoConstraints CON;
  //typedef Dune::PDELab::ISTLVectorBackend
  //  <Dune::PDELab::istl::Parameters::static_blocking,blocksize> VBE;

  // blocking fixed is not allowed???
  using VBE0 = Dune::PDELab::istl::VectorBackend<>;

  using VBE = Dune::PDELab::istl::VectorBackend<Dune::PDELab::istl::Blocking::fixed>;
  using OrderingTag = Dune::PDELab::EntityBlockedOrderingTag;
  typedef Dune::PDELab::GridFunctionSpace<GV,FEMDG,CON,VBE0> GFSDG;
  GFSDG gfsdg(gv,femdg);

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
  

  typedef typename GFS::template ConstraintsContainer<Real>::Type C;
  C cg;
  gfs.update(); // initializing the gfs
  std::cout << "degrees of freedom: " << gfs.globalSize() << std::endl;

  // <<<2b>>> define problem parameters
  typedef RiemannProblem<GV,Real> Param;
  Param param;

  // <<<4>>> Make grid operator
  typedef Dune::PDELab::DGLinearAcousticsSpatialOperator<Param,FEMDG> LOP;
  LOP lop(param);
  typedef Dune::PDELab::DGLinearAcousticsTemporalOperator<Param,FEMDG> TLOP;
  TLOP tlop(param);
  Dune::PDELab::ExplicitEulerParameter<Real> method1;
  Dune::PDELab::HeunParameter<Real> method2;
  Dune::PDELab::Shu3Parameter<Real> method3;
  Dune::PDELab::RK4Parameter<Real> method4;
  Dune::PDELab::TimeSteppingParameterInterface<Real> *method;
  if (degree==0) {method=&method1; std::cout << "setting explicit Euler" << std::endl;}
  if (degree==1) {method=&method2; std::cout << "setting Heun" << std::endl;}
  if (degree==2) {method=&method3; std::cout << "setting Shu 3" << std::endl;}
  if (degree==3) {method=&method4; std::cout << "setting RK4" << std::endl;}

  typedef Dune::PDELab::istl::BCRSMatrixBackend<> MBE;
  MBE mbe(5); // Maximal number of nonzeroes per row can be cross-checked by patternStatistics().

  typedef Dune::PDELab::GridOperator
    <GFS,GFS,LOP,MBE,Real,Real,Real,C,C> GO0;
  GO0 go0(gfs,cg,gfs,cg,lop,mbe);

  typedef Dune::PDELab::GridOperator
    <GFS,GFS,TLOP,MBE,Real,Real,Real,C,C> GO1;
  GO1 go1(gfs,cg,gfs,cg,tlop,mbe);

  typedef Dune::PDELab::OneStepGridOperator<GO0,GO1,false> IGO;
  IGO igo(go0,go1);
  igo.setMethod(*method);

  // <<<5>>> set initial values
  typedef typename IGO::Traits::Domain V;
  V xold(gfs,0.0);
  Dune::PDELab::LinearAcousticsInitialValueAdapter<Param> u0(gv,param);

  Dune::PDELab::interpolate(u0,gfs,xold);
 

  //only interpolate components
  


  // <<<6>>> Make a linear solver backend
  //typedef Dune::PDELab::ISTLBackend_SEQ_CG_SSOR LS;
  //LS ls(10000,1);
  typedef Dune::PDELab::ISTLBackend_OVLP_ExplicitDiagonal<GFS> LS;
  LS ls(gfs);

  // <<<8>>> time-stepper
  typedef Dune::PDELab::CFLTimeController<Real,IGO> TC;
  TC tc(0.999,igo);
  Dune::PDELab::ExplicitOneStepMethod<Real,IGO,LS,V,V,TC> osm(*method,igo,ls,tc);
  osm.setVerbosityLevel(2);

  // <<<10>>> graphics for initial guess
  /*
  Dune::PDELab::FilenameHelper fn(name);
  int counter=0;
  {
    typedef Dune::PDELab::VectorDiscreteGridFunction<GFS,V> DGF;
    DGF xdgf(gfs,xold);
    int refinement = std::max(degree-1,0);
    if (degree>=2) refinement+=2;
    Dune::SubsamplingVTKWriter<GV> vtkwriter(gv,refinement);
    vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<DGF>(xdgf,"u"));
    vtkwriter.pwrite(fn.getName(),"vtk","",Dune::VTK::appendedraw);
    fn.increment();
  }
  */

  // prepare VTK writer and write first file
  //int subsampling=ptree.get("output.subsampling",(int)0);
  int subsampling = 0;
  using VTKWRITER = Dune::SubsamplingVTKWriter<GV>;
  VTKWRITER vtkwriter(gv,subsampling);

  //std::string filename=ptree.get("output.filename","output") + STRINGIZE(GRIDDIM) "d";
  //Dune::PDELab::FilenameHelper filename(name);
  std::string filename(name);

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




  // <<<11>>> time loop
  Real time = 0.0;
  Real dt = timestep;
  V x(gfs,0.0);
  while (time < Tend)
    {
      // do time step
      osm.apply(time,dt,xold,x);

      /*
      // graphics
      counter++;
      if (counter%modulo==0)
        {
          typedef Dune::PDELab::VectorDiscreteGridFunction<GFS,V> DGF;
          DGF xdgf(gfs,x);
          int refinement = std::max(degree-1,0);
          if (degree>=2) refinement+=2;
          Dune::SubsamplingVTKWriter<GV> vtkwriter(gv,refinement);
          vtkwriter.addVertexData(new Dune::PDELab::VTKGridFunctionAdapter<DGF>(xdgf,"u"));
          vtkwriter.pwrite(fn.getName(),"vtk","",Dune::VTK::appendedraw);
          fn.increment();
        }
      */

      // output to VTK file
      vtkSequenceWriter.write(time,Dune::VTK::appendedraw);

      xold = x;
      time += dt;
    }
}

//===============================================================
// Main program with grid setup
//===============================================================

int main(int argc, char** argv)
{
  try{
    //Maybe initialize Mpi
    Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);
    if(Dune::MPIHelper::isFake)
      std::cout<< "This is a sequential program." << std::endl;
    else
	  {
		if(helper.rank()==0)
		  std::cout << "parallel run on " << helper.size() << " process(es)" << std::endl;
	  }
    /*
    if (argc!=7)
      {
        if(helper.rank()==0)
          {
            std::cout << "usage: " << argv[0] << " <end time> <time step> <grid file> <refinement> <degree> <modulo>" << std::endl;
            std::cout << "         <grid file> = 'yaspgrid' || <a gmsh file>"  << std::endl;
            std::cout << "         <refinement> = nonnegative integer, initial mesh has h=1/20" << std::endl;
            std::cout << "         <modulo> = write vtk file every modulo'th time step" << std::endl;
            std::cout << "coarse example:" << std::endl;
            std::cout << "./heterogeneoussquare 0.1 0.001 'yaspgrid' 1 1 5" << std::endl;
          }
        return 1;
      }
    */
    double Tend(1.0);
    //sscanf(argv[1],"%lg",&Tend);
    double timestep(0.001);
    //sscanf(argv[2],"%lg",&timestep);
    std::string grid_file = "yaspgrid"; //(argv[3]);
    int max_level(1);// sscanf(argv[4],"%d",&max_level);
    int p(0);// sscanf(argv[5],"%d",&p);
    int modulo(1);// sscanf(argv[6],"%d",&modulo);

    // parallel overlapping yaspgrid version
    if (grid_file=="yaspgrid")
      {
        const int dim=2;
        Dune::FieldVector<double,dim> L(1.0);
        Dune::array<int,dim> N(Dune::fill_array<int,dim>(20));
        std::bitset<dim> periodic(false);
        int overlap=1;
        //Dune::YaspGrid<2> grid(helper.getCommunicator(),L,N,periodic,overlap);
        Dune::YaspGrid<2> grid(L,N,periodic,overlap,helper.getCommunicator());
 
        for (int i=0; i<max_level; i++) grid.globalRefine(1);
        typedef Dune::YaspGrid<2>::LeafGridView GV;
        const GV& gv=grid.leafGridView();
        if (p==0)
          {
            const int degree=0;
            typedef Dune::PDELab::OPBLocalFiniteElementMap<GV::Grid::ctype,double,degree,dim,Dune::GeometryType::cube> FEM;
            FEM fem;
            std::stringstream fullname;
            fullname << grid_file << "_l" << max_level << "_k" << p;
            explicit_scheme<GV,FEM,degree>(gv,fem,Tend,timestep,fullname.str(),modulo);
          }
        if (p==1)
          {
            const int degree=1;
            typedef Dune::PDELab::OPBLocalFiniteElementMap<GV::Grid::ctype,double,degree,dim,Dune::GeometryType::cube> FEM;
            FEM fem;
            std::stringstream fullname;
            fullname << grid_file << "_l" << max_level << "_k" << p;
            explicit_scheme<GV,FEM,degree>(gv,fem,Tend,timestep,fullname.str(),modulo);
          }
        if (p==2)
          {
            const int degree=2;
            typedef Dune::PDELab::OPBLocalFiniteElementMap<GV::Grid::ctype,double,degree,dim,Dune::GeometryType::cube> FEM;
            FEM fem;
            std::stringstream fullname;
            fullname << grid_file << "_l" << max_level << "_k" << p;
            explicit_scheme<GV,FEM,degree>(gv,fem,Tend,timestep,fullname.str(),modulo);
          }
        if (p==3)
          {
            const int degree=3;
            typedef Dune::PDELab::OPBLocalFiniteElementMap<GV::Grid::ctype,double,degree,dim,Dune::GeometryType::cube> FEM;
            FEM fem;
            std::stringstream fullname;
            fullname << grid_file << "_l" << max_level << "_k" << p;
            explicit_scheme<GV,FEM,degree>(gv,fem,Tend,timestep,fullname.str(),modulo);
          }
        return 0;
      }

#if HAVE_UG
    if (true)
      {
        // make uggrid
        const int dim=2;
        typedef Dune::UGGrid<dim> GridType;
        GridType grid;
        typedef std::vector<int> GmshIndexMap;
        GmshIndexMap boundary_index_map;
        GmshIndexMap element_index_map;
        Dune::GridFactory<GridType> factory(&grid);
        Dune::GmshReader<GridType> gmsh_reader;
        gmsh_reader.read(factory,grid_file,boundary_index_map,element_index_map,true,false);
        for (int i=0; i<max_level; i++) grid.globalRefine(1);
        typedef GridType::LeafGridView GV;
        const GV& gv=grid.leafGridView();
        if (p==0)
          {
            const int degree=0;
            typedef Dune::PDELab::OPBLocalFiniteElementMap<GV::Grid::ctype,double,degree,dim,Dune::GeometryType::simplex> FEM;
            FEM fem;
            std::stringstream fullname;
            fullname << grid_file << "_l" << max_level << "_k" << p;
            explicit_scheme<GV,FEM,degree>(gv,fem,Tend,timestep,fullname.str(),modulo);
          }
        if (p==1)
          {
            const int degree=1;
            typedef Dune::PDELab::OPBLocalFiniteElementMap<GV::Grid::ctype,double,degree,dim,Dune::GeometryType::simplex> FEM;
            FEM fem;
            std::stringstream fullname;
            fullname << grid_file << "_l" << max_level << "_k" << p;
            explicit_scheme<GV,FEM,degree>(gv,fem,Tend,timestep,fullname.str(),modulo);
          }
        if (p==2)
          {
            const int degree=2;
            typedef Dune::PDELab::OPBLocalFiniteElementMap<GV::Grid::ctype,double,degree,dim,Dune::GeometryType::simplex> FEM;
            FEM fem;
            std::stringstream fullname;
            fullname << grid_file << "_l" << max_level << "_k" << p;
            explicit_scheme<GV,FEM,degree>(gv,fem,Tend,timestep,fullname.str(),modulo);
          }
        if (p==3)
          {
            const int degree=3;
            typedef Dune::PDELab::OPBLocalFiniteElementMap<GV::Grid::ctype,double,degree,dim,Dune::GeometryType::simplex> FEM;
            FEM fem;
            std::stringstream fullname;
            fullname << grid_file << "_l" << max_level << "_k" << p;
            explicit_scheme<GV,FEM,degree>(gv,fem,Tend,timestep,fullname.str(),modulo);
          }
      }
#endif

  }
  catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
	return 1;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
	return 1;
  }
}

