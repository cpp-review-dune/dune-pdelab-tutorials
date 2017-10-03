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
  static constexpr int mplus = 1; // number of positive eigenvectors
  static constexpr int mminus = 1; // number of negative eigenvectors
  static constexpr int mstar = mplus+mminus; // number of nonzero eigenvalues

  using RangeField = typename PROBLEM::RangeField;

  Model (PROBLEM& p)
  : problem(p)
  {
  }

  template<typename E, typename X, typename T2, typename T3>
  void eigenvectors (const E& e, const X& x, const Dune::FieldVector<T2,dim>& n, Dune::FieldMatrix<T3,m,m>& RT) const
  {
    auto c = problem.c(e,x);

    RT[0][0] =  1.0/c;  RT[0][1] = -1.0/c;  RT[0][2] = 0.0;
    RT[1][0] =  n[0]; RT[1][1] = n[0];  RT[1][2] = -n[1];
    RT[2][0] =  n[1]; RT[2][1] = n[1];  RT[2][2] = n[0];
  }

  // where do we need this matrix for?
  // template<typename RF>
  // static void coefficients (Dune::FieldMatrix<RF,m,m>& A)
  // {
  //   int c2 = 1;

  //   A[0][0] = 0.0; A[0][1] = 1.0; A[0][2] = 1.0;
  //   A[1][0] = c2;  A[1][1] = 0.0; A[1][2] = 0.0;
  //   A[2][0] = c2;  A[2][1] = 0.0; A[2][2] = 0.0;
  // }

  template<typename E, typename X, typename RF>
  void diagonal (const E& e, const X& x, Dune::FieldMatrix<RF,m,m>& D) const
  {
    auto c = problem.c(e,x);

    D[0][0] = c;    D[0][1] = 0.0; D[0][2] = 0.0;
    D[1][0] = 0.0;  D[1][1] = -c ; D[1][2] = 0.0;
    D[2][0] = 0.0;  D[2][1] = 0.0; D[2][2] = 0.0;
  }

  //Flux function
  template<typename E, typename X, typename RF>
  void flux (const E& e, const X& x, const Dune::FieldVector<RF,m>& u, Dune::FieldMatrix<RF,m,dim>& F) const
  {
    auto c = problem.c(e,x);

    F[0][0] = u[1]    ; F[0][1] = u[2];
    F[1][0] = c*c*u[0]; F[1][1] = 0.0;
    F[2][0] = 0.0     ; F[2][1] = c*c*u[0];
  }

  const PROBLEM& problem;
};

/*
not tested
template<>
class Model<1>
{

public:
  static constexpr int dim = 1;
  static constexpr int m = 2;


  template<typename T1, typename T2, typename T3>
  static void eigenvectors_transposed (T1 c, const Dune::FieldVector<T2,dim>& n, Dune::FieldMatrix<T3,m,m>& RT)
  {
    RT[0][0] =  1; RT[0][1] = c*n[0];
    RT[1][0] = -1; RT[1][1] = c*n[0];
  }

  template<typename T1, typename T2, typename T3>
  static void eigenvectors (T1 c, const Dune::FieldVector<T2,dim>& n, Dune::FieldMatrix<T3,m,m>& RT)
  {
    RT[0][0] =  1; RT[1][0] = c*n[0];
    RT[0][1] = -1; RT[1][1] = c*n[0];
  }

};

template<>
class Model<3>
{

public:
  static constexpr int dim = 3;
  static constexpr int m = 4;

  template<typename T1, typename T2, typename T3>
  static void eigenvectors_transposed (T1 c, const Dune::FieldVector<T2,dim>& n, Dune::FieldMatrix<T3,m,m>& RT)
  {
    Dune::FieldVector<T2,dim> s;
    s[0] = n[1]-n[2]; s[1] = n[2]-n[0]; s[2] = n[0]-n[1];
    if (s.two_norm()<1e-5)
      {
        s[0] = n[1]+n[2]; s[1] = n[2]-n[0]; s[2] = -(n[0]+n[1]);
      }

    Dune::FieldVector<T2,dim> t; // compute cross product s * n
    t[0] = n[1]*s[2] - n[2]*s[1];
    t[1] = n[2]*s[0] - n[0]*s[2];
    t[2] = n[0]*s[1] - n[1]*s[0];

    RT[0][0] =  0;  RT[0][1] =   s[0];  RT[0][2] =   s[1];  RT[0][3] =   s[2];
    RT[1][0] =  0;  RT[1][1] =   t[0];  RT[1][2] =   t[1];  RT[1][3] =   t[2];
    RT[2][0] =  1;  RT[2][1] = c*n[0];  RT[2][2] = c*n[1];  RT[2][3] = c*n[2];
    RT[3][0] = -1;  RT[3][1] = c*n[0];  RT[3][2] = c*n[1];  RT[3][3] = c*n[2];
  }

  template<typename T1, typename T2, typename T3>
  static void eigenvectors (T1 c, const Dune::FieldVector<T2,dim>& n, Dune::FieldMatrix<T3,m,m>& RT)
  {
    Dune::FieldVector<T2,dim> s;
    s[0] = n[1]-n[2]; s[1] = n[2]-n[0]; s[2] = n[0]-n[1];
    if (s.two_norm()<1e-5)
      {
        s[0] = n[1]+n[2]; s[1] = n[2]-n[0]; s[2] = -(n[0]+n[1]);
      }

    Dune::FieldVector<T2,dim> t; // compute cross product s * n
    t[0] = n[1]*s[2] - n[2]*s[1];
    t[1] = n[2]*s[0] - n[0]*s[2];
    t[2] = n[0]*s[1] - n[1]*s[0];

    RT[0][0] =  0;  RT[1][0] =   s[0];  RT[2][0] =   s[1];  RT[3][0] =   s[2];
    RT[0][1] =  0;  RT[1][1] =   t[0];  RT[2][1] =   t[1];  RT[3][1] =   t[2];
    RT[0][2] =  1;  RT[1][2] = c*n[0];  RT[2][2] = c*n[1];  RT[3][2] = c*n[2];
    RT[0][3] = -1;  RT[1][3] = c*n[0];  RT[2][3] = c*n[1];  RT[3][3] = c*n[2];
  }
};
*/
#endif //ACOUSTICS_MODEL
