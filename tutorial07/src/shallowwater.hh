#ifndef SHALLOWWATER_MODEL
#define SHALLOWWATER_MODEL
/*
 Shallow Water model class
*/

template<int dim, typename PROBLEM>
class Model ;

template<typename PROBLEM>
class Model<1,PROBLEM>
{

public:
  static constexpr int dim = 1;
  static constexpr int m = 2;

  using RangeField = typename PROBLEM::RangeField;

  Model (PROBLEM& p)
  : problem(p)
  {
  }

  template<typename RF>
  static void max_eigenvalue (const Dune::FieldVector<RF,m>& u_s,
                              const Dune::FieldVector<RF,m>& u_n,
                              const Dune::FieldVector<RF,dim>& n_F,
                              RF& alpha)
  {
    int g = 1;

    RF alpha_s(0.0);
    RF alpha_n(0.0);

      alpha_s = std::abs(u_s[1]/u_s[0]) + sqrt(g*u_s[0]);
      alpha_n = std::abs(u_n[1]/u_n[0]) + sqrt(g*u_n[0]);

      alpha = std::max(alpha_s, alpha_n);
  }

  //Flux function
  template<typename E, typename X, typename RF>
  void flux (const E& e, const X& x, const Dune::FieldVector<RF,m>& u, Dune::FieldMatrix<RF,m,dim>& F) const
  {
    int g = 1;

    F[0][0] = u[1]                          ;
    F[1][0] = u[1]*u[1]/u[0] + 0.5*u[0]*u[0];
  }

  const PROBLEM& problem;
};


template<typename PROBLEM>
class Model<2,PROBLEM>
{

public:
  static constexpr int dim = 2;
  static constexpr int m = 3;

  using RangeField = typename PROBLEM::RangeField;

  Model (PROBLEM& p)
  : problem(p)
  {
  }

  template<typename RF>
  static void max_eigenvalue (const Dune::FieldVector<RF,m>& u_s,
                              const Dune::FieldVector<RF,m>& u_n,
                              const Dune::FieldVector<RF,dim>& n_F,
                              RF& alpha)
  {
    int g = 1;

    RF alpha_s(0.0);
    RF alpha_n(0.0);

    for(size_t k=0;k<dim;++k)
      {
        alpha_s +=  n_F[k]*u_s[k+1]/u_s[0];
        alpha_n += -n_F[k]*u_n[k+1]/u_n[0];
      }
      alpha_s = std::abs(alpha_s) + sqrt(g*u_s[0]);
      alpha_n = std::abs(alpha_n) + sqrt(g*u_n[0]);

      alpha = std::max(alpha_s, alpha_n);
  }


  //Flux function
  template<typename E, typename X, typename RF>
  void flux (const E& e, const X& x, const Dune::FieldVector<RF,m>& u, Dune::FieldMatrix<RF,m,dim>& F) const
  {
    int g = 1;

    F[0][0] = u[1]                          ; F[0][1] = u[2];
    F[1][0] = u[1]*u[1]/u[0] + 0.5*u[0]*u[0]; F[1][1] = u[1]*u[2]/u[0];
    F[2][0] = u[1]*u[2]/u[0]                ; F[2][1] = u[2]*u[2]/u[0] + 0.5*u[0]*u[0];
  }

  const PROBLEM& problem;
};

#endif //SHALLOWWATER_MODEL
