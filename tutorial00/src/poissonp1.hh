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
  DomainType qp;          // center of mass of refelem
  double weight;          // quadrature weight on refelem
  double phihat[n];       // basis functions at qp
  double gradhat[dim][n]; // coordinate x #basisfct 
  
public:
  // pattern assembly flags
  enum { doPatternVolume = true };

  // residual assembly flags
  enum { doAlphaVolume = true };
  enum { doLambdaVolume = true };

  // Constructor precomputes element independent data
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
      if (q>0) {
        std::cout << "This is the wrong quadrature rule!" << std::endl;
        exit(1);
      }
      
      // position and weight of the quadrature point
      weight = ip.weight();
      qp = ip.position();

      // print quadrature point data
      std::cout << "Quadrature: " << ip.position() << " " << ip.weight() << std::endl;
      ++q;
    }

    // evaluate basis functions on refelem
    std::vector<RangeType> phi(n);
    fel.localBasis().evaluateFunction(qp,phi);
    for (int i=0; i<n; i++) phihat[i] = phi[i];

    // evaluate gradients of basis functions on refelem
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
