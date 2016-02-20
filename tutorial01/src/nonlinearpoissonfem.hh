#include<vector>

#include<dune/common/exceptions.hh>
#include<dune/common/fvector.hh>
#include<dune/geometry/type.hh>

#include<dune/pdelab/common/referenceelements.hh>
#include<dune/pdelab/common/quadraturerules.hh>
#include<dune/pdelab/localoperator/pattern.hh>
#include<dune/pdelab/localoperator/flags.hh>
#include<dune/pdelab/localoperator/defaultimp.hh>
#include<dune/pdelab/finiteelement/localbasiscache.hh>

/** a local operator for solving the nonlinear Poisson equation with conforming FEM
 *
 * \f{align*}{
 *   -\Delta u(x) + q(u(x)) &=& f(x) x\in\Omega,  \\
 *                     u(x) &=& g(x) x\in\partial\Omega_D \\
 *  -\nabla u(x) \cdot n(x) &=& j(x) x\in\partial\Omega_N \\
 * \f}
 *
 */
template<typename Param, typename FEM>
class NonlinearPoissonFEM :
  public Dune::PDELab::NumericalJacobianVolume<NonlinearPoissonFEM<Param,FEM> >,
  public Dune::PDELab::NumericalJacobianApplyVolume<NonlinearPoissonFEM<Param,FEM> >,
  public Dune::PDELab::FullVolumePattern,
  public Dune::PDELab::LocalOperatorDefaultFlags
{
  typedef typename FEM::Traits::FiniteElementType::Traits::LocalBasisType LocalBasis;
  Dune::PDELab::LocalBasisCache<LocalBasis> cache; // a cache for local basis evaluations
  Param& param; // parameter functions
  int incrementorder; // additional increment for integration order

public:
  // pattern assembly flags
  enum { doPatternVolume = true };

  // residual assembly flags
  enum { doLambdaVolume = true };
  enum { doLambdaBoundary = true };
  enum { doAlphaVolume = true };

  //! constructor stores a copy of the parameter object
  NonlinearPoissonFEM (Param& param_, int incrementorder_=0)
    : param(param_), incrementorder(incrementorder_)
  {}

  //! right hand side integral
  template<typename EG, typename LFSV, typename R>
  void lambda_volume (const EG& eg, const LFSV& lfsv,
                      R& r) const
  {
    // types & dimension
    const int dim = EG::Entity::dimension;

    // select quadrature rule
    auto geo = eg.geometry();
    const int order = incrementorder+2*lfsv.finiteElement().localBasis().order();
    auto rule = Dune::PDELab::quadratureRule(geo,order);

    // loop over quadrature points
    for (const auto& ip : rule)
      {
        // evaluate basis functions
        auto& phihat = cache.evaluateFunction(ip.position(),lfsv.finiteElement().localBasis());

        // integrate -f*phi_i
        decltype(ip.weight()) factor = ip.weight() * geo.integrationElement(ip.position());
        for (size_t i=0; i<lfsv.size(); i++)
          r.accumulate(lfsv,i,-param.f(eg.entity(),ip.position())*phihat[i]*factor);
      }
  }

    // boundary integral
  template<typename IG, typename LFSV, typename R>
  void lambda_boundary (const IG& ig, const LFSV& lfsv, R& r) const
  {
    // types
    const int dim = IG::Entity::dimension;

    // evaluate boundary condition type
    auto localgeo = ig.geometryInInside();
    auto facecenterlocal = Dune::PDELab::referenceElement(localgeo).position(0,0);
    bool isdirichlet = param.b(ig.intersection(),facecenterlocal);

    // skip rest if we are on Dirichlet boundary
    if (isdirichlet) return;

    // select quadrature rule
    auto globalgeo = ig.geometry();
    const int order = incrementorder+2*lfsv.finiteElement().localBasis().order();
    auto rule = Dune::PDELab::quadratureRule(globalgeo,order);

    // loop over quadrature points and integrate normal flux
    for (const auto& ip : rule)
      {
        // position of quadrature point in local coordinates of element
        auto local = localgeo.global(ip.position());

        // evaluate shape functions (assume Galerkin method)
        auto& phihat = cache.evaluateFunction(local,lfsv.finiteElement().localBasis());

        // integrate j
        decltype(ip.weight()) factor = ip.weight()*globalgeo.integrationElement(ip.position());
        for (size_t i=0; i<lfsv.size(); i++)
          r.accumulate(lfsv,i,param.j(ig.intersection(),ip.position())*phihat[i]*factor);
      }
  }

  //! volume integral depending on test and ansatz functions
  template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
  void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
  {
    // types & dimension
    const int dim = EG::Entity::dimension;
    typedef decltype(makeZeroBasisFieldValue(lfsu)) RF;

    // select quadrature rule
    auto geo = eg.geometry();
    const int order = incrementorder+2*lfsu.finiteElement().localBasis().order();
    auto rule = Dune::PDELab::quadratureRule(geo,order);

    // loop over quadrature points
    for (const auto& ip : rule)
      {
        // evaluate basis functions
        auto& phihat = cache.evaluateFunction(ip.position(),lfsu.finiteElement().localBasis());

        // evaluate u
        RF u=0.0;
        for (size_t i=0; i<lfsu.size(); i++) u += x(lfsu,i)*phihat[i];

        // evaluate gradient of shape functions
        auto& gradphihat = cache.evaluateJacobian(ip.position(),lfsu.finiteElement().localBasis());

        // transform gradients of shape functions to real element
        const auto S = geo.jacobianInverseTransposed(ip.position());
        auto gradphi = makeJacobianContainer(lfsu);
        for (size_t i=0; i<lfsu.size(); i++)
          S.mv(gradphihat[i][0],gradphi[i][0]);

        // compute gradient of u
        Dune::FieldVector<RF,dim> gradu(0.0);
        for (size_t i=0; i<lfsu.size(); i++)
          gradu.axpy(x(lfsu,i),gradphi[i][0]);

        // integrate (grad u)*grad phi_i + q(u)*phi_i
        RF factor = ip.weight() * geo.integrationElement(ip.position());
        for (size_t i=0; i<lfsu.size(); i++)
          r.accumulate(lfsu,i,(gradu*gradphi[i][0] + param.q(u)*phihat[i])*factor);
      }
  }
};