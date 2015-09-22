/******************************************************/
/** \brief Class defining Dirichlet boundary conditions 
    and initial guess
 */
/******************************************************/
template<typename GV, typename RF>
class GFunction
 : public Dune::PDELab::GridFunctionBase<Dune::PDELab::
          GridFunctionTraits<GV,RF,1,Dune::FieldVector<RF,1> >,
          GFunction<GV,RF> >
{
  const GV& gv;
public:
  typedef Dune::PDELab::
  GridFunctionTraits<GV,RF,1,Dune::FieldVector<RF,1> > Traits;

  //! construct from grid view
  GFunction (const GV& gv_) : gv(gv_) {}

  //! get a reference to the grid view
  inline const GV& getGridView () {return gv;}

  //! evaluate extended function on element
  inline void evaluate (const typename Traits::ElementType& e, 
                        const typename Traits::DomainType& xl,
                        typename Traits::RangeType& y) const
  {  
    const int dim = Traits::GridViewType::Grid::dimension;
    typename Traits::DomainType xg = e.geometry().global(xl);
    y = 0.0;
    for (int i=0; i<dim; i++) y += xg[i]*xg[i];
    return;
  }
};
