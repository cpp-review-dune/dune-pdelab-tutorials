

template<typename MODEL>
class NumericalFlux 
{



  //Flux Vector splitting
  template<typename RF>
  static void fluxvectorsplitting (RF c, Dune::FieldVector<RF,m>& u, Dune::FieldMatrix<RF,dim,m>& F) 
	{
    F[0][0] = u[1]    ; F[0][1] = u[2];
    F[1][0] = c*c*u[0]; F[1][1] = 0.0;
    F[2][0] = 0.0     ; F[2][1] = c*c*u[0];
	}

  //local Lax-Friedrichs Flux


}
