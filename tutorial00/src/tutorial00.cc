// -*- tab-width: 4; indent-tabs-mode: nil -*-
/** \file

    \brief Solve elliptic problem in unconstrained spaces with conforming finite elements
*/
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include<math.h>
#include<iostream>
#include<vector>
#include<map>
#include<string>

#include<dune/common/parallel/mpihelper.hh>
#include<dune/common/parametertreeparser.hh>
#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>
#include<dune/common/static_assert.hh>
#include<dune/common/timer.hh>
#include<dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include<dune/grid/io/file/gmshreader.hh>
#include<dune/grid/yaspgrid.hh>
#include<dune/grid/onedgrid.hh>
#if HAVE_ALBERTA
#include<dune/grid/albertagrid.hh>
#include <dune/grid/albertagrid/dgfparser.hh>
#endif
#if HAVE_UG
#include<dune/grid/uggrid.hh>
#endif
#if HAVE_DUNE_ALUGRID
#include<dune/alugrid/grid.hh>
#include<dune/grid/io/file/dgfparser/dgfalu.hh>
#include<dune/grid/io/file/dgfparser/dgfparser.hh>
#endif
#include<dune/istl/bvector.hh>
#include<dune/istl/operators.hh>
#include<dune/istl/solvers.hh>
#include<dune/istl/preconditioners.hh>
#include<dune/istl/io.hh>
#include<dune/istl/superlu.hh>
#include<dune/pdelab/newton/newton.hh>
#include<dune/pdelab/finiteelementmap/p0fem.hh>
#include<dune/pdelab/finiteelementmap/pkfem.hh>
#include<dune/pdelab/finiteelementmap/qkfem.hh>
#include<dune/pdelab/finiteelementmap/rannacherturekfem.hh>
#include<dune/pdelab/constraints/common/constraints.hh>
#include<dune/pdelab/constraints/common/constraintsparameters.hh>
#include<dune/pdelab/constraints/conforming.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include<dune/pdelab/gridfunctionspace/genericdatahandle.hh>
#include<dune/pdelab/gridfunctionspace/interpolate.hh>
#include<dune/pdelab/gridfunctionspace/vtk.hh>
#include<dune/pdelab/common/function.hh>
#include<dune/pdelab/common/vtkexport.hh>
#include<dune/pdelab/backend/istl.hh>
#include<dune/pdelab/stationary/linearproblem.hh>
#include<dune/pdelab/gridoperator/gridoperator.hh>
#include<dune/geometry/referenceelements.hh>
#include<dune/geometry/quadraturerules.hh>
#include<dune/pdelab/common/geometrywrapper.hh>
#include<dune/pdelab/localoperator/defaultimp.hh>
#include<dune/pdelab/localoperator/pattern.hh>
#include<dune/pdelab/localoperator/flags.hh>

/******************************************************/
/** \brief Class defining the right hand side
 */
/******************************************************/
template<typename GV, typename RF>
class FFunction
  : public Dune::PDELab::GridFunctionBase<Dune::PDELab::
           GridFunctionTraits<GV,RF,1,Dune::FieldVector<RF,1> >, FFunction<GV,RF> >
{
  const GV& gv;
  RF time;
public:
  typedef Dune::PDELab::GridFunctionTraits<GV,RF,1,Dune::FieldVector<RF,1> > Traits;

  //! construct from grid view
  FFunction (const GV& gv_) : gv(gv_) {}

  //! evaluate extended function on element
  inline void evaluate (const typename Traits::ElementType& e, 
                        const typename Traits::DomainType& xlocal,
                        typename Traits::RangeType& y) const
  {  
    const int dim = Traits::GridViewType::Grid::dimension;
    typedef typename Traits::GridViewType::Grid::ctype ctype;
    Dune::FieldVector<ctype,dim> x = e.geometry().global(xlocal);
    for (int i=0; i<dim; i++) y -= 2.0;
    return;
  }
  
  //! get a reference to the grid view
  inline const GV& getGridView () {return gv;}
};

/******************************************************/
/** \brief Class defining Dirichlet boundary conditions and initial guess
 */
/******************************************************/
template<typename GV, typename RF>
class GFunction
  : public Dune::PDELab::GridFunctionBase<Dune::PDELab::
           GridFunctionTraits<GV,RF,1,Dune::FieldVector<RF,1> >, GFunction<GV,RF> >
{
  const GV& gv;
  RF time;
public:
  typedef Dune::PDELab::GridFunctionTraits<GV,RF,1,Dune::FieldVector<RF,1> > Traits;

  //! construct from grid view
  GFunction (const GV& gv_) : gv(gv_) {}

  //! evaluate extended function on element
  inline void evaluate (const typename Traits::ElementType& e, 
                        const typename Traits::DomainType& xlocal,
                        typename Traits::RangeType& y) const
  {  
    const int dim = Traits::GridViewType::Grid::dimension;
    typedef typename Traits::GridViewType::Grid::ctype ctype;
    Dune::FieldVector<ctype,dim> x = e.geometry().global(xlocal);
    y = 0.0;
    for (int i=0; i<dim; i++) y += x[i]*x[i];
    return;
  }
  
  //! get a reference to the grid view
  inline const GV& getGridView () {return gv;}
};

/******************************************************/
//! \brief Parameter class selecting type of boundary conditions
/******************************************************/
class BCType
  : public Dune::PDELab::DirichletConstraintsParameters
{
public:
  //! Test whether boundary is Dirichlet-constrained
  template<typename I>
  bool isDirichlet(const I & intersection,
                   const Dune::FieldVector<typename I::ctype, I::dimension-1> & coord
                   ) const
  {
    return true;  // Dirichlet b.c. all over
  }
};


/******************************************************/
/** a local operator for solving the linear convection-diffusion equation with standard FEM
 *
 * \f{align*}{
 *   \nabla\cdot(-A(x) \nabla u + b(x) u) + c(x)u &=& f \mbox{ in } \Omega,  \\
 *                                              u &=& g \mbox{ on } \partial\Omega_D \\
 *                (b(x,u) - A(x)\nabla u) \cdot n &=& j \mbox{ on } \partial\Omega_N \\
 *                        -(A(x)\nabla u) \cdot n &=& j \mbox{ on } \partial\Omega_O
 * \f}
 * Note:
 *  - This formulation is valid for velocity fields which are non-divergence free.
 *  - Outflow boundary conditions should only be set on the outflow boundary
 *
 * \tparam T model of ConvectionDiffusionParameterInterface
 */
/******************************************************/
template<typename G, typename F, typename FiniteElementMap>
class PoissonP1 :
  public Dune::PDELab::FullVolumePattern,
  public Dune::PDELab::LocalOperatorDefaultFlags
{
private:
  // define useful types
  typedef typename FiniteElementMap::Traits::FiniteElementType FiniteElementType;
  typedef typename FiniteElementType::Traits::LocalBasisType LocalBasisType;
  typedef typename LocalBasisType::Traits::JacobianType JacobianType;
  typedef typename LocalBasisType::Traits::RangeType RangeType;
  typedef typename LocalBasisType::Traits::DomainType DomainType;
  typedef typename LocalBasisType::Traits::RangeFieldType RF;

  // data members
  const G& g;
  const F& f;
  enum {dim=LocalBasisType::Traits::dimDomain}; 
  enum {n=dim+1};
  DomainType qp;
  double weight;  // quadrature weight on reference element
  double phihat[n]; // basis functions at qp
  double gradhat[dim][n]; // coordinate x #basisfct 
  
public:
  // pattern assembly flags
  enum { doPatternVolume = true };

  // residual assembly flags
  enum { doAlphaVolume = true };
  enum { doLambdaVolume = true };

  PoissonP1 (const G& g_, const F& f_, const FiniteElementType& fel)
    : g(g_), f(f_)
  {
    // determine the single quadrature point
    // and evaluate basis functions on the reference element

    // select quadrature rule
    Dune::GeometryType gt = fel.type();
    const Dune::QuadratureRule<RF,dim>& rule = Dune::QuadratureRules<RF,dim>::rule(gt,1);

    // loop over quadrature points 
    int q=0;
    for (const auto& ip : rule) {
      if (q>0) break; // there should only be one quadrature point !

      // position and weight of the quadrature point
      weight = ip.weight();
      qp = ip.position();

      // print quadrature point data
      std::cout << "Quadrature: " << ip.position() << " " << ip.weight() << std::endl;
      ++q;
    }

    // evaluate basis functions
    std::vector<RangeType> phi(n);
    fel.localBasis().evaluateFunction(qp,phi);
    for (int i=0; i<n; i++) phihat[i] = phi[i];

    // evaluate gradients of basis functions
    std::vector<JacobianType> js(n);
    fel.localBasis().evaluateJacobian(qp,js);
    for (int i=0; i<n; i++)
      for (int j=0; j<dim; j++)
        gradhat[j][i] = js[i][0][j];
  }

  // volume integral depending on test and ansatz functions
  template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
  void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
  {
    // get Jacobian and determinant
    // assume the transformation is linear
    const auto geo = eg.geometry();
    const auto S = geo.jacobianInverseTransposed(qp);
    RF factor = weight*geo.integrationElement(qp);

    // compute gradients of basis functions in transformed element
    double grad[dim][n] = {0.0};  // coordinate x #basisfct
    for (int i=0; i<dim; i++) // rows of S
      for (int k=0; k<dim; k++) // columns of S
        for (int j=0; j<n; j++) // columns of gradhat
          grad[i][j] += S[i][k] * gradhat[k][j];

    // compute gradient u_h 
    double graduh[dim] = {0.0};
    double z[n];
    for (int j=0; j<n; j++) z[j] = x(lfsu,j); // read data into array
    for (int k=0; k<dim; k++) // rows of grad
      for (int j=0; j<n; j++) // columns of grad
        graduh[k] += grad[k][j]*z[j];

    // scalar products
    double s[n] = {0.0};
    for (int k=0; k<dim; k++) // rows of grad
      for (int j=0; j<n; j++) 
        s[j] += graduh[k]*grad[k][j];

    // store in result
    for (int i=0; i<n; i++) 
      r.accumulate(lfsv,i,s[i]*factor);
  }

  //! apply local jacobian of the volume term
  template<typename EG, typename LFSU, typename X, typename LFSV,
           typename Y>
  void jacobian_apply_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, Y& y) const
  {
    alpha_volume(eg,lfsu,x,lfsv,y);
  }
    
  // jacobian of volume term
  template<typename EG, typename LFSU, typename X, typename LFSV, typename M>
  void jacobian_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv,
                        M& mat) const
  {
    // get Jacobian and determinant
    // assume the transformation is linear
    const auto geo = eg.geometry();
    const auto S = geo.jacobianInverseTransposed(qp);
    RF factor = weight*geo.integrationElement(qp);

    // compute gradients of basis functions in transformed element
    double grad[dim][n] = {0.0}; // coordinate x #basisfct
    for (int i=0; i<dim; i++) // rows of S
      for (int k=0; k<dim; k++) // columns of S
        for (int j=0; j<n; j++) // columns of gradhat
          grad[i][j] += S[i][k] * gradhat[k][j];

    // compute grad^T * grad
    double A[n][n] = {0.0};
    for (int i=0; i<n; i++)
      for (int k=0; k<dim; k++)
        for (int j=0; j<n; j++)
          A[i][j] += grad[k][i]*grad[k][j];
    
    // store in result
    for (int i=0; i<n; i++)
      for (int j=0; j<n; j++)
        mat.accumulate(lfsu,i,lfsu,j,A[i][j]*factor);
  }

  // volume integral depending only on test functions
  template<typename EG, typename LFSV, typename R>
  void lambda_volume (const EG& eg, const LFSV& lfsv, R& r) const
  {
    typename F::Traits::RangeType fval;
    f.evaluate(eg.entity(),qp,fval);
    RF factor = fval*weight*eg.geometry().integrationElement(qp);
    for (int i=0; i<n; i++)
      r.accumulate(lfsv,i,-factor*phihat[i]);
  }
};


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
  typedef Dune::PDELab::ISTLBackend_SEQ_CG_SSOR LS;
  LS ls(5000,true);

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

//===============================================================
// Main program with grid setup
//===============================================================
int main(int argc, char** argv)
{
  try{
    // Maybe initialize Mpi
    Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);
    if(Dune::MPIHelper::isFake)
      std::cout<< "This is a sequential program." << std::endl;
    else
      std::cout << "Parallel code run on " << helper.size() << " process(es)" << std::endl;

    // Read parameters from ini file
    Dune::ParameterTree ptree;
    Dune::ParameterTreeParser ptreeparser;
    ptreeparser.readINITree("tutorial00.ini",ptree);
    ptreeparser.readOptions(argc,argv,ptree);

    // make grid
    const unsigned int dim = ptree.get<int>("grid.dim");
    const unsigned int refinement = ptree.get<int>("grid.refinement");
    if (dim==1)
      {
        // read grid parameters from input file
        typedef Dune::OneDGrid::ctype ctype;
        const ctype a = ptree.get<ctype>("grid.oned.a");
        const ctype b = ptree.get<ctype>("grid.oned.b");
        const unsigned int N = ptree.get<int>("grid.oned.elements");

        // create equidistant intervals
        typedef std::vector<ctype> Intervals;
        Intervals intervals(N+1);
        for(unsigned int i=0; i<N+1; ++i)
          intervals[i] = a + ctype(i)*(b-a)/ctype(N);
        
        // Construct grid
        typedef Dune::OneDGrid Grid;
        Grid grid(intervals);
        grid.globalRefine(refinement);

        // call generic function
        typedef Grid::LeafGridView GV;
        GV gv = grid.leafGridView();
        driver(gv,ptree);
      }
#if HAVE_UG
    if (dim==2)
      {
        std::string filename = ptree.get("grid.twod.filename","unitsquare.msh");
        typedef Dune::UGGrid<2> Grid;
        Dune::GridFactory<Grid> factory;
        Dune::GmshReader<Grid>::read(factory,filename,true,true);
        std::shared_ptr<Grid> gridp = std::shared_ptr<Grid>(factory.createGrid());
        Dune::Timer timer;
        gridp->globalRefine(refinement);
        std::cout << "Time for mesh refinement " << timer.elapsed() << " seconds" << std::endl;
        typedef Grid::LeafGridView GV;
        GV gv = gridp->leafGridView();
        driver(gv,ptree);
      }
    if (dim==3)
      {
        std::string filename = ptree.get("grid.threed.filename","unitcube.msh");
        typedef Dune::UGGrid<3> Grid;
        Dune::GridFactory<Grid> factory;
        Dune::GmshReader<Grid>::read(factory,filename,true,true);
        std::shared_ptr<Grid> gridp = std::shared_ptr<Grid>(factory.createGrid());
        Dune::Timer timer;
        gridp->globalRefine(refinement);
        std::cout << "Time for mesh refinement " << timer.elapsed() << " seconds" << std::endl;
        typedef Grid::LeafGridView GV;
        GV gv = gridp->leafGridView();
        driver(gv,ptree);
      }
#else
      std::cout << "This example requires UG!" << std::endl;
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
