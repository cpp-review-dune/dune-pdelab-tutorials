/*
 Maxwell model class
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
class Model<3>
{

public:
  static constexpr int dim = 3;
  static constexpr int m = 6;


  //minak
  /*
	template<typename T1, typename T2, typename T3>
	static void eigenvectors (T1 c, const Dune::FieldVector<T2,dim>& n, Dune::FieldMatrix<T3,m,m>& R)
  {
    R[0][0] = 0.0;  R[0][1] = 0.0;  R[0][2] = 0.0;  R[0][3] = 0.0;  R[0][4] = n[2]; R[0][5] =-n[1];
    R[1][0] = 0.0;  R[1][1] = 0.0;  R[1][2] = 0.0;  R[1][3] =-n[2]; R[1][4] = 0.0;  R[1][5] = n[0];
    R[2][0] = 0.0;  R[2][1] = 0.0;  R[2][2] = 0.0;  R[2][3] = n[1]; R[2][4] =-n[0]; R[2][5] = 1.0;
		R[3][0] = 0.0;  R[3][1] =-n[2]; R[3][2] = n[1]; R[3][3] = 0.0;  R[3][4] = 1.0;  R[3][5] = 1.0;
		R[4][0] = n[2]; R[4][1] = 0.0;  R[4][2] =-n[0]; R[4][3] = 0.0;  R[4][4] = 1.0;  R[4][5] = 1.0;
		R[5][0] =-n[1]; R[5][1] = n[0]; R[5][2] = 0.0;  R[5][3] = 0.0;  R[5][4] = 1.0;  R[5][5] = 1.0;

    //c = sqrtmueps
		R *= 1/c;
  }
  */

  //peter
	template<typename T1, typename T2, typename T3>
	static void eigenvectors (T1 cc, const Dune::FieldVector<T2,dim>& n, Dune::FieldMatrix<T3,m,m>& R)
  {

    auto a=n[0], b=n[1], c=n[2];

    R[0][0] = a*b;  R[0][1] =-c;    R[0][2] = c;    R[0][3] = a*b;  R[0][4] = n[0]; R[0][5] = 0.0;
    R[1][0] = b*b-1;R[1][1] = 0.0;  R[1][2] = 0.0;  R[1][3] = b*b-1;R[1][4] = n[1]; R[1][5] = 0.0;
    R[2][0] = b*c;  R[2][1] = a;    R[2][2] = -a;   R[2][3] = b*c;  R[2][4] = n[2]; R[2][5] = 0.0;
		R[3][0] = c;    R[3][1] = a*b;  R[3][2] = a*b;  R[3][3] =-c;    R[3][4] = 0.0;  R[3][5] = n[0];
		R[4][0] = 0.0;  R[4][1] = b*b-1;R[4][2] = b*b-1;R[4][3] = 0.0;  R[4][4] = 0.0;  R[4][5] = n[1];
		R[5][0] = -a;   R[5][1] = b*c;  R[5][2] = b*c;  R[5][3] = a;    R[5][4] = 0.0;  R[5][5] = n[2];

    //c = sqrtmueps
		R *= 1/cc;
  }


  //one can also provide eigenvectors inverse

  template<typename RF>
  static void coefficients (RF mu, Dune::FieldMatrix<RF,m,m>& A) 
  {
    A[0][0] = 0.0;   A[0][1] = 0.0;   A[0][2] = 0.0;   A[0][3] = 0.0;   A[0][4] = 1./mu; A[0][5] =-1./mu;
    A[1][0] = 0.0;   A[1][1] = 0.0;   A[1][2] = 0.0;   A[1][3] =-1./mu; A[1][4] = 0.0;   A[1][5] = 1./mu;
    A[2][0] = 0.0;   A[2][1] = 0.0;   A[2][2] = 0.0;   A[2][3] = 1./mu; A[2][4] =-1./mu; A[2][5] = 1.0;
		A[3][0] = 0.0;   A[3][1] =-1./mu; A[3][2] = 1./mu; A[3][3] = 0.0;   A[3][4] = 1.0;   A[3][5] = 1.0;
		A[4][0] = 1./mu; A[4][1] = 0.0;   A[4][2] =-1./mu; A[4][3] = 0.0;   A[4][4] = 1.0;   A[4][5] = 1.0;
		A[5][0] =-1./mu; A[5][1] = 1./mu; A[5][2] = 0.0;   A[5][3] = 0.0;   A[5][4] = 1.0;   A[5][5] = 1.0;
	}

  template<typename RF>
  static void diagonal (RF c, Dune::FieldMatrix<RF,m,m>& D) 
  {
    D[0][0] = c;   D[0][1] = 0.0; D[0][2] = 0.0; D[0][3] = 0.0; D[0][4] = 0.0; D[0][5] = 0.0;
    D[1][0] = 0.0; D[1][1] = c;   D[1][2] = 0.0; D[1][3] = 0.0; D[1][4] = 0.0; D[1][5] = 0.0;
    D[2][0] = 0.0; D[2][1] = 0.0; D[2][2] = -c;  D[2][3] = 0.0; D[2][4] = 0.0; D[2][5] = 0.0;
		D[3][0] = 0.0; D[3][1] = 0.0; D[3][2] = 0.0; D[3][3] = -c;  D[3][4] = 0.0; D[3][5] = 0.0;
		D[4][0] = 0.0; D[4][1] = 0.0; D[4][2] = 0.0; D[4][3] = 0.0; D[4][4] = 0.0; D[4][5] = 0.0;
		D[5][0] = 0.0; D[5][1] = 0.0; D[5][2] = 0.0; D[5][3] = 0.0; D[5][4] = 0.0; D[5][5] = 0.0;
	}

  //Flux function
  template<typename RF>
  static void flux (RF c, Dune::FieldVector<RF,m>& u, Dune::FieldMatrix<RF,dim,m>& F) 
	{
    RF mu(1.0);
    RF ep(1.0);

    F[0][0] = 0.0      ; F[0][1] =-1/mu*u[5]; F[0][2] = 1/mu*u[4];  
    F[1][0] = 1/mu*u[5]; F[1][1] = 0.0;       F[1][2] =-1/mu*u[3]; 
    F[2][0] =-1/mu*u[4]; F[2][1] = 1/mu*u[3]; F[2][2] = 0.0;
    F[3][0] = 0.0      ; F[3][1] = 1/ep*u[2]; F[3][2] =-1/ep*u[1];  
    F[4][0] =-1/ep*u[2]; F[4][1] = 0.0;       F[4][2] = 1/ep*u[0]; 
    F[5][0] = 1/ep*u[1]; F[5][1] =-1/ep*u[0]; F[5][2] = 0.0;
	}


};

