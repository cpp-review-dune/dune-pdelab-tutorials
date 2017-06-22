#ifndef RIEMANNPROBLEM
#define RIEMANNPROBLEM
//template<typename GV, typename Range, typename M>
template<typename GV, typename Number>
class RiemannProblem
{
public:
  //typedef Dune::PDELab::LinearAcousticsParameterTraits<GV,Number> Traits;

  using RangeField = Number;
  using Range = Dune::FieldVector<Number,3>;


  //RiemannProblem (M& m)
  //  : m_(m), time(0.0),  pi(3.141592653589793238462643)
  //{
  //}
  RiemannProblem ()
    : time(0.0),  pi(3.141592653589793238462643)
  {
    const int dim = 2;
  }

  
  //! speed of sound
  /*
  typename Traits::RangeFieldType
  c (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    typename Traits::DomainType xglobal = e.geometry().global(x);
    if ( xglobal[1] < 1-(0.6/0.9)*(xglobal[0]-0.1) ) return 1.0;
    //    if (xglobal[0]>0.5) return 2.0;
    return 0.5;
  }
  */

  //! speed of sound
  template<typename E, typename X>
  Number c (const E& e, const X& x) const
  {
    X xglobal = e.geometry().global(x);
    if ( xglobal[1] < 1-(0.6/0.9)*(xglobal[0]-0.1) ) return 1.0;
    //    if (xglobal[0]>0.5) return 2.0;
    return 0.5;
  }



  //! Neumann boundary condition
  template<typename I, typename X>
  Number j (const I& i, const X& x) const
  {
    return 0.0;
  }


  //! Boundary condition value - reflecting bc
  template<typename I, typename X, typename R>
  Range g (const I& is, const X& x, const R& s) const
  {
    Range u(0.0);
    u[0] =  s[0];
    u[1] = -s[1];
    u[2] = -s[2];
    return u;
  }

  //! right hand side
  /*
  typename Traits::StateType
  q (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    typename Traits::StateType rhs(0.0);
    return rhs;
  }
  */
  //! right hand side
  template<typename E, typename X>
  Range q (const E& e, const X& x) const
  {
    Range rhs(0.0);
    return rhs;
  }

  //! initial value 
  template<typename E, typename X>
  Range u0 (const E& e, const X& x) const
  {
    X xglobal = e.geometry().global(x);
    Range u(0.0);
    if (xglobal[0]>0.45 && xglobal[0]<0.55 && xglobal[1]>0.3 && xglobal[1]<0.4)
      {
        u[0] = sin(pi*(xglobal[0]-0.45)/0.1)*sin(pi*(xglobal[0]-0.45)/0.1)*sin(pi*(xglobal[1]-0.3)/0.1)*sin(pi*(xglobal[1]-0.3)/0.1);
        u[1] = 0;
        u[2] = 0;
      }
    return u;
  }


  //! set time for subsequent evaluation
  void setTime (Number t)
  {
    time = t;
  }

private:
  //M& m_;
  Number time;
  Number pi;

};
#endif
