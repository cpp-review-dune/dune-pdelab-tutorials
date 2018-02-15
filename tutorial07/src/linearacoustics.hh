#ifndef ACOUSTICS_MODEL
#define ACOUSTICS_MODEL
/*
 Linear Acoustics model class
*/
/** \brief provide matrix which contains rowwise the eigenvectors of linear acoustics problem
    \tparam dim the space dimension
    \param c speed of sound
    \param n unit outer normal vector
    \param RT matrix to be filled
*/

template<int dim, typename PROBLEM>
class Model ;

template<typename PROBLEM>
class Model<2,PROBLEM>
{
public:
  static constexpr int dim = 2;   // space dimension
  static constexpr int m = dim+1; // system size
  static constexpr int mplus = 1; // number of positive eigenvalues
  static constexpr int mminus = 1; // number of negative eigenvalues
  static constexpr int mstar = mplus+mminus; // number of nonzero eigenvalues

  using RangeField = typename PROBLEM::RangeField;

  Model (PROBLEM& p)
  : problem(p)
  {
  }

  template<typename E, typename X, typename T2, typename T3>
  void eigenvectors (const E& e, const X& x, 
                     const Dune::FieldVector<T2,dim>& n,
                     Dune::FieldMatrix<T3,m,m>& RT) const
  {
    auto c = problem.c(e,x);

    RT[0][0] =  1.0/c;  RT[0][1] = -1.0/c;  RT[0][2] = 0.0;
    RT[1][0] =  n[0]; RT[1][1] = n[0];  RT[1][2] = -n[1];
    RT[2][0] =  n[1]; RT[2][1] = n[1];  RT[2][2] = n[0];
  }

  template<typename E, typename X, typename RF>
  void diagonal (const E& e, const X& x, Dune::FieldMatrix<RF,m,m>& D) const
  {
    auto c = problem.c(e,x);

    D[0][0] = c;    D[0][1] = 0.0; D[0][2] = 0.0;
    D[1][0] = 0.0;  D[1][1] = -c ; D[1][2] = 0.0;
    D[2][0] = 0.0;  D[2][1] = 0.0; D[2][2] = 0.0;
  }

  template<typename RF>
  static void max_eigenvalue (const Dune::FieldVector<RF,m>& u_s,
                              const Dune::FieldVector<RF,m>& u_n,
                              const Dune::FieldVector<RF,dim>& n_F,
                              RF& alpha)
  {
    alpha = 1.0;
  }


  //Flux function
  template<typename E, typename X, typename RF>
  void flux (const E& e, const X& x, 
             const Dune::FieldVector<RF,m>& u,
             Dune::FieldMatrix<RF,m,dim>& F) const
  {
    auto c = problem.c(e,x);

    F[0][0] = u[1]    ; F[0][1] = u[2];
    F[1][0] = c*c*u[0]; F[1][1] = 0.0;
    F[2][0] = 0.0     ; F[2][1] = c*c*u[0];
  }

  const PROBLEM& problem;
};

#endif //ACOUSTICS_MODEL
