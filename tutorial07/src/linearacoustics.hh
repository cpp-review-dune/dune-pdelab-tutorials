/*
 Linear Acoustics model class
*/
/** \brief provide matrix which contains rowwise the eigenvectors of linear acoustics problem
	  \tparam dim the space dimension
    \param c speed of sound
    \param n unit outer normal vector
    \param RT matrix to be filled
*/

template<int dim>
class Model ;


template<>
class Model<2>
{

public:
  static constexpr int dim = 2;
  static constexpr int m = 3;

	template<typename T1, typename T2, typename T3>
	static void eigenvectors (T1 c, const Dune::FieldVector<T2,dim>& n, Dune::FieldMatrix<T3,m,m>& RT)
	{
	  RT[0][0] =  0; RT[1][0] =  -n[1];  RT[2][0] = n[0];
	  RT[0][1] =  1; RT[1][1] = c*n[0];  RT[2][1] = c*n[1];
	  RT[0][2] = -1; RT[1][2] = c*n[0];  RT[2][2] = c*n[1];
	}

  //one can also provide eigenvectors inverse

  template<typename RF>
  static void coefficients (RF c2, Dune::FieldMatrix<RF,m,m>& A) 
  {
    A[0][0] = 0.0; A[0][1] = 1.0; A[0][2] = 1.0;
    A[1][0] = c2;  A[1][1] = 0.0; A[1][2] = 0.0;
    A[2][0] = c2;  A[2][1] = 0.0; A[2][2] = 0.0;
	}


  template<typename RF>
  static void diagonal (RF c, Dune::FieldMatrix<RF,m,m>& D) 
  {
    D[0][0] = 0.0;  D[0][1] = 0.0; D[0][2] = 0.0;
    D[1][0] = 0.0;  D[1][1] = c  ; D[1][2] = 0.0;
    D[2][0] = 0.0;  D[2][1] = 0.0; D[2][2] = -c;
	}

  //Flux function
  template<typename RF>
  static void flux (RF c, Dune::FieldVector<RF,m>& u, Dune::FieldMatrix<RF,dim,m>& F) 
	{
    F[0][0] = u[1]    ; F[0][1] = u[2];
    F[1][0] = c*c*u[0]; F[1][1] = 0.0;
    F[2][0] = 0.0     ; F[2][1] = c*c*u[0];
	}


};

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

