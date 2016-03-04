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


/** a local operator for residual-based error estimation for the problem
 *
 * \f{align*}{
 *   -\Delta u(x) + q(u(x)) &=& f(x) x\in\Omega,  \\
 *                     u(x) &=& g(x) x\in\partial\Omega_D \\
 *  -\nabla u(x) \cdot n(x) &=& j(x) x\in\partial\Omega_N \\
 * \f}
 *
 * A call to residual() of a grid operator space will assemble
 * the quantity \f$\eta_T^2\f$ for each cell. Note that the squares
 * of the cell indicator \f$\eta_T\f$ is stored. To compute the global
 * error estimate sum up all values and take the square root.
 *
 * Assumptions and limitations:
 * - Assumes that LFSU is \f$P_k\f$/\f$Q_k\f$ finite element space
 *   and LFSV is a \f$P_0\f$ finite element space (one value per cell).
 * - However, the second order derivatives are ignored!
 *
 */
template<typename Param, typename FEM>
class NonlinearPoissonFEMEstimator
  : public Dune::PDELab::LocalOperatorDefaultFlags
{
  typedef typename FEM::Traits::FiniteElementType::Traits::LocalBasisType LocalBasis;
  Dune::PDELab::LocalBasisCache<LocalBasis> cache; // a cache for local basis evaluations
  Param& param; // parameter functions
  int incrementorder; // additional increment for integration order

  // a function to compute the diameter of an entity
  template<class GEO>
  typename GEO::ctype diameter (const GEO& geo) const
  {
    typedef typename GEO::ctype DF;
    DF hmax = -1.0E00;
    for (int i=0; i<geo.corners(); i++)
      {
        auto xi = geo.corner(i);
        for (int j=i+1; j<geo.corners(); j++)
          {
            auto xj = geo.corner(j);
            xj -= xi;
            hmax = std::max(hmax,xj.two_norm());
          }
      }
    return hmax;
  }

public:
  // pattern assembly flags
  enum { doPatternVolume = false };
  enum { doPatternSkeleton = false };

  // residual assembly flags
  enum { doAlphaVolume  = true };
  enum { doAlphaSkeleton  = true };
  enum { doAlphaBoundary  = true };

  //! constructor: pass parameter object
  NonlinearPoissonFEMEstimator (Param& param_, int incrementorder_=0)
    : param(param_), incrementorder(incrementorder_)
  {}

  // volume integral depending on test and ansatz functions
  template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
  void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
  {
    // types & dimension
    typedef decltype(makeZeroBasisFieldValue(lfsu)) RF;

    // select quadrature rule
    auto geo = eg.geometry();
    const int order = incrementorder+2*lfsu.finiteElement().localBasis().order();
    auto rule = Dune::PDELab::quadratureRule(geo,order);

    // loop over quadrature points
    RF sum(0.0);
    for (const auto& ip : rule)
      {
        // evaluate basis functions
        auto& phihat = cache.evaluateFunction(ip.position(),lfsu.finiteElement().localBasis());

        // evaluate u
        RF u=0.0;
        for (size_t i=0; i<lfsu.size(); i++) u += x(lfsu,i)*phihat[i];

        // evaluate reaction term
        auto q = param.q(u);

        // evaluate right hand side parameter function
        auto f = param.f(eg.entity(),ip.position());

        // integrate f^2
        RF factor = ip.weight() * geo.integrationElement(ip.position());
        sum += (f-q)*(f-q)*factor;
      }

    // accumulate cell indicator
    auto h_T = diameter(eg.geometry());
    r.accumulate(lfsv,0,h_T*h_T*sum);
  }

  // skeleton integral depending on test and ansatz functions
  // each face is only visited ONCE!
  template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
  void alpha_skeleton (const IG& ig,
                       const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                       const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n,
                       R& r_s, R& r_n) const
  {
    // geometries in local coordinates of the elements
    auto insidegeo = ig.geometryInInside();
    auto outsidegeo = ig.geometryInOutside();

    // inside and outside cells
    auto cell_inside = ig.inside();
    auto cell_outside = ig.outside();

    // geometries from local to global in elements
    auto geo_s = cell_inside.geometry();
    auto geo_n = cell_outside.geometry();

    // dimensions
    const int dim = IG::Entity::dimension;

    // select quadrature rule
    auto globalgeo = ig.geometry();
    const int order = incrementorder+2*lfsu_s.finiteElement().localBasis().order();
    auto rule = Dune::PDELab::quadratureRule(globalgeo,order);

    // loop over quadrature points and integrate normal flux
    typedef decltype(makeZeroBasisFieldValue(lfsu_s)) RF;
    RF sum(0.0);
    for (const auto& ip : rule)
      {
        // position of quadrature point in local coordinates of elements
        auto iplocal_s = insidegeo.global(ip.position());
        auto iplocal_n = outsidegeo.global(ip.position());

        // unit outer normal direction
        auto n_F = ig.unitOuterNormal(ip.position());

        // gradient in normal direction in self
        auto& gradphihat_s = cache.evaluateJacobian(iplocal_s,lfsu_s.finiteElement().localBasis());
        const auto S_s = geo_s.jacobianInverseTransposed(iplocal_s);
        RF gradun_s = 0.0;
        for (size_t i=0; i<lfsu_s.size(); i++)
          {
            Dune::FieldVector<RF,dim> v;
            S_s.mv(gradphihat_s[i][0],v);
            gradun_s += x_s(lfsu_s,i)*(v*n_F);
          }

        // gradient in normal direction in neighbor
        auto& gradphihat_n = cache.evaluateJacobian(iplocal_n,lfsu_n.finiteElement().localBasis());
        const auto S_n = geo_n.jacobianInverseTransposed(iplocal_n);
        RF gradun_n = 0.0;
        for (size_t i=0; i<lfsu_n.size(); i++)
          {
            Dune::FieldVector<RF,dim> v;
            S_n.mv(gradphihat_n[i][0],v);
            gradun_n += x_n(lfsu_n,i)*(v*n_F);
          }

        // integrate
        RF factor = ip.weight()*globalgeo.integrationElement(ip.position());
        RF jump = gradun_s-gradun_n;
        sum += jump*jump*factor;
      }

    // accumulate indicator
    auto h_T = diameter(globalgeo);
    r_s.accumulate(lfsv_s,0,0.5*h_T*sum);
    r_n.accumulate(lfsv_n,0,0.5*h_T*sum);
  }

  // boundary integral depending on test and ansatz functions
  // We put the Dirchlet evaluation also in the alpha term to save some geometry evaluations
  template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
  void alpha_boundary (const IG& ig,
                       const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                       R& r_s) const
  {
    // geometries in local coordinates of the elements
    auto insidegeo = ig.geometryInInside();

    // inside and outside cells
    auto cell_inside = ig.inside();

    // geometries from local to global in elements
    auto geo_s = cell_inside.geometry();

    // dimensions
    const int dim = IG::Entity::dimension;

    // select quadrature rule
    auto globalgeo = ig.geometry();
    const int order = incrementorder+2*lfsu_s.finiteElement().localBasis().order();
    auto rule = Dune::PDELab::quadratureRule(globalgeo,order);

    // loop over quadrature points and integrate normal flux
    typedef decltype(makeZeroBasisFieldValue(lfsu_s)) RF;
    RF sum(0.0);
    for (const auto& ip : rule)
      {
        // skip body if we are on Dirichlet boundary
        if (param.b(ig.intersection(),ip.position())) continue;

        // position of quadrature point in local coordinates of elements
        auto iplocal_s = insidegeo.global(ip.position());

        // unit outer normal direction
        auto n_F = ig.unitOuterNormal(ip.position());

        // gradient in normal direction in self
        auto& gradphihat_s = cache.evaluateJacobian(iplocal_s,lfsu_s.finiteElement().localBasis());
        const auto S_s = geo_s.jacobianInverseTransposed(iplocal_s);
        RF gradun_s = 0.0;
        for (size_t i=0; i<lfsu_s.size(); i++)
          {
            Dune::FieldVector<RF,dim> v;
            S_s.mv(gradphihat_s[i][0],v);
            gradun_s += x_s(lfsu_s,i)*(v*n_F);
          }

        // Neumann boundary condition value
        auto j = param.j(ig.intersection(),ip.position());

        // integrate
        RF factor = ip.weight()*globalgeo.integrationElement(ip.position());
        RF jump = gradun_s+j;
        sum += jump*jump*factor;
      }

    // accumulate indicator
    auto h_T = diameter(globalgeo);
    r_s.accumulate(lfsv_s,0,h_T*sum);
  }
};
