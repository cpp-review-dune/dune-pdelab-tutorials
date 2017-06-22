template<typename GV, typename RF>
class Krivodonova4Problem
{
public:
  typedef Dune::PDELab::LinearAcousticsParameterTraits<GV,RF> Traits;

  Krivodonova4Problem ()
    : pi(3.141592653589793238462643), time(0.0)
  {
  }

  //! speed of sound
  typename Traits::RangeFieldType
  c (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    return 340.0;
  }

  //! Dirichlet boundary condition value
  typename Traits::StateType
  g (const typename Traits::IntersectionType& is, const typename Traits::IntersectionDomainType& x, const typename Traits::StateType& s) const
  {
    typename Traits::DomainType xglobal = is.geometry().global(x);
    if (xglobal[0]<1e-6)
      {
        typename Traits::StateType u(0.0);
        u[1] = 1.224*(1+0.5*sin(2*pi*1500.0*time));
        return u;
      }
    if (xglobal[0]>1.0-1e-6)
      {
        typename Traits::StateType u(0.0);
        return u;
      }
    typename Traits::StateType u(0.0);
    u[2] = 0.0;
    return u;
  }

  //! right hand side
  typename Traits::StateType
  q (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    typename Traits::StateType rhs(0.0);
    return rhs;
  }

  //! initial value
  typename Traits::StateType
  u0 (const typename Traits::ElementType& e, const typename Traits::DomainType& x) const
  {
    typename Traits::StateType u(0.0);
    return u;
  }

  //! set time for subsequent evaluation
  void setTime (RF t)
  {
    time = t;
  }

private:
  double pi;
  RF time;
};

