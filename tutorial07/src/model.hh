
/*
 Linear Acoustics model class
*/

template<int dim>
class Model ;


template<>
class Model<2>
{

public:

  static constexpr int dim = 2;
  static constexpr int m = 3;

  //Model ();


	/** \brief provide matrix which contains rowwise the eigenvectors of linear acoustics problem
		  \tparam dim the space dimension
	*/
	//class Eigenvectors
	//{
	//	public:
		/**
		   \param c speed of sound
		   \param n unit outer normal vector
		   \param RT matrix to be filled
		*/
		template<typename T1, typename T2, typename T3>
		static void eigenvectors_transposed (T1 c, const Dune::FieldVector<T2,dim>& n, Dune::FieldMatrix<T3,m,m>& RT)
		{
		  RT[0][0] =  0; RT[0][1] =  -n[1];  RT[0][2] = n[0];
		  RT[1][0] =  1; RT[1][1] = c*n[0];  RT[1][2] = c*n[1];
		  RT[2][0] = -1; RT[2][1] = c*n[0];  RT[2][2] = c*n[1];
		}

		template<typename T1, typename T2, typename T3>
		static void eigenvectors (T1 c, const Dune::FieldVector<T2,dim>& n, Dune::FieldMatrix<T3,m,m>& RT)
		{
		  RT[0][0] =  0; RT[1][0] =  -n[1];  RT[2][0] = n[0];
		  RT[0][1] =  1; RT[1][1] = c*n[0];  RT[2][1] = c*n[1];
		  RT[0][2] = -1; RT[1][2] = c*n[0];  RT[2][2] = c*n[1];
		}
  //};
};



