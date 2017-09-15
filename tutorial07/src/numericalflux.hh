/*template<typename MODEL>
class NumericalFlux
{
public:

  using Model = MODEL;

  using RF = typename MODEL::RangeField;
  static constexpr int m = MODEL::Model::m;

  using NumericalFlux = Dune::FieldVector<RF,m>;

  const Model& model() const;

  //NumericalFlux numericalFlux(ig,x, u_s, u_n, gradu_s, gradu_n) const;
};
*/

//Flux Vector splitting
template<typename MODEL>
class FluxVectorSplitting
{

public:

  static constexpr int dim = MODEL::dim;
  static constexpr int m = MODEL::Model::m;

  using RF = typename MODEL::RangeField; // type for computations
  using DF = RF;

  FluxVectorSplitting (const MODEL& model_)
    : model(model_)
  {
  }

  //using NumericalFlux = Dune::FieldVector<RF,m>;

  template<typename IG, typename X>
  static void numericalFlux(const IG& ig, const X& x,
                              const Dune::FieldVector<RF,m>& u_s,
                              const Dune::FieldVector<RF,m>& u_n,Dune::FieldVector<RF,m>& f)
                              //const std::vector<Dune::FieldVector<RF,dim> >& gradu_s,
                              //const std::vector<Dune::FieldVector<RF,dim> >& gradu_n) const
  {
    // Normal: assume faces are planar

    const Dune::FieldVector<DF,dim> n_F = ig.centerUnitOuterNormal();

    Dune::FieldMatrix<DF,m,m> D(0.0);
    // fetch eigenvalues
    MODEL::diagonal(D);

    Dune::FieldMatrix<DF,m,m> Dplus(0.0);
    Dune::FieldMatrix<DF,m,m> Dminus(0.0);

    for (size_t i =0 ; i<m;i++)
      (D[i][i] > 0) ? Dplus[i][i] = D[i][i] : Dminus[i][i] = D[i][i];

    // fetch eigenvectors
    Dune::FieldMatrix<DF,m,m> Rot;
    MODEL::eigenvectors(n_F,Rot);

    // compute B+ = RD+R^-1 and B- = RD-R^-1
    Dune::FieldMatrix<DF,m,m> Bplus(Rot);
    Dune::FieldMatrix<DF,m,m> Bminus(Rot);

    //multiply by D+-
    Bplus.rightmultiply(Dplus);
    Bminus.rightmultiply(Dminus);

    //multiply by R^-1
    Rot.invert();
    Bplus.rightmultiply(Rot);
    Bminus.rightmultiply(Rot);

    // Initialize vectors outside for loop
    //Dune::FieldVector<RF,m> u_s(0.0);
    //Dune::FieldVector<RF,m> u_n(0.0);
    //Dune::FieldVector<RF,m> f(0.0);

    // Compute numerical flux at  the integration point
    f = 0.0;
    // f = Bplus*u_s + Bminus*u_n
    Bplus.umv(u_s,f);
    Bminus.umv(u_n,f);

  }

  const MODEL& model;

};// FVS

//local Lax-Friedrichs Flux
template<typename MODEL>
class LLFflux
{

public:

  static constexpr int dim = MODEL::dim;
  static constexpr int m = MODEL::Model::m;

  using RF = typename MODEL::RangeField; // type for computations
  using DF = RF;

  LLFflux (const MODEL& model_)
    : model(model_)
  {
  }

  //using NumericalFlux = Dune::FieldVector<RF,m>;

  template<typename IG, typename X>
  static void numericalFlux(const IG& ig, const X& x,
                              const Dune::FieldVector<RF,m>& u_s,
                              const Dune::FieldVector<RF,m>& u_n,Dune::FieldVector<RF,m>& f)
  {
    // Normal: assume faces are planar
    const Dune::FieldVector<DF,dim> n_F = ig.centerUnitOuterNormal();

    Dune::FieldMatrix<DF,m,m> D(0.0);
    // fetch eigenvalues
    MODEL::diagonal(D);

    Dune::FieldMatrix<DF,m,m> Dplus(0.0);
    Dune::FieldMatrix<DF,m,m> Dminus(0.0);

    for (size_t i =0 ; i<m;i++)
      (D[i][i] > 0) ? Dplus[i][i] = D[i][i] : Dminus[i][i] = D[i][i];

    // fetch eigenvectors
    Dune::FieldMatrix<DF,m,m> Rot;
    MODEL::eigenvectors(n_F,Rot);

    // compute B+ = RD+R^-1 and B- = RD-R^-1
    Dune::FieldMatrix<DF,m,m> Bplus(Rot);
    Dune::FieldMatrix<DF,m,m> Bminus(Rot);

    //multiply by D+-
    Bplus.rightmultiply(Dplus);
    Bminus.rightmultiply(Dminus);

    //multiply by R^-1
    Rot.invert();
    Bplus.rightmultiply(Rot);
    Bminus.rightmultiply(Rot);

    // Compute numerical flux at  the integration point
    f = 0.0;
    // f = Bplus*u_s + Bminus*u_n
    Bplus.umv(u_s,f);
    Bminus.umv(u_n,f);

  }

  const MODEL& model;

};// FVS
