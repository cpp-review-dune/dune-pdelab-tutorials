/******************************************************/
/** \brief Class defining the right hand side
 */
/******************************************************/
template<typename GV, typename RF>
class FFunction
 : public Dune::PDELab::GridFunctionBase<Dune::PDELab::
          GridFunctionTraits<GV,RF,1,Dune::FieldVector<RF,1> >,
          FFunction<GV,RF> >
{
  const GV& gv;
public:
  typedef Dune::PDELab::
  GridFunctionTraits<GV,RF,1,Dune::FieldVector<RF,1> > Traits;

  //! construct from grid view
  FFunction (const GV& gv_) : gv(gv_) {}

  //! evaluate extended function on element
  inline void evaluate (const typename Traits::ElementType& e, 
                        const typename Traits::DomainType& xl,
                        typename Traits::RangeType& y) const
  {  
    const int dim = Traits::GridViewType::Grid::dimension;
    typename Traits::DomainType xg = e.geometry().global(xl);
    y = 0.0;
    for (int i=0; i<dim; i++) y -= 2.0;
    return;
  }
  
  //! get a reference to the grid view
  inline const GV& getGridView () {return gv;}
};
