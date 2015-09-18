/******************************************************/
//! \brief Parameter class selecting type of boundary conditions
/******************************************************/
class BCType
  : public Dune::PDELab::DirichletConstraintsParameters
{
public:
  //! Test whether boundary is Dirichlet-constrained
  template<typename I>
  bool isDirichlet(const I & intersection,
                   const Dune::FieldVector<typename I::ctype, I::dimension-1> & coord
                   ) const
  {
    return true;  // Dirichlet b.c. all over
  }
};
