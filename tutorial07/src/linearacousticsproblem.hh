#ifndef ACOUSTICS_RIEMANNPROBLEM
#define ACOUSTICS_RIEMANNPROBLEM
template<typename GV, typename NUMBER>
class Problem
{
public:

  using RangeField = NUMBER;

  //problem specification depends on dimension
  static constexpr int dim = 2;
  static constexpr int m = 3;

  using Range = Dune::FieldVector<NUMBER,m>;

  Problem ()
    : time(0.0),  pi(3.141592653589793238462643)
  {
  }

  //! speed of sound
  template<typename E, typename X>
  NUMBER c (const E& e, const X& x) const
  {
    X xglobal = e.geometry().global(x);
    if ( xglobal[1] < 1-(0.6/0.9)*(xglobal[0]-0.1) ) return 1.0;

    return 0.5;
  }

  //! Neumann boundary condition
  template<typename I, typename X>
  NUMBER j (const I& i, const X& x) const
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
  template<typename E, typename X>
  Range q (const E& e, const X& x) const
  {
    Range rhs(0.0);
    return rhs;
  }

  //! initial value -> the same as tutorial04
  template<typename E, typename X>
  Range u0 (const E& e, const X& x) const
  {
    X xglobal = e.geometry().global(x);
    Range u(0.0);
    for (int i=0; i<dim; i++)
      u[0] += (xglobal[i]-0.375)*(xglobal[i]-0.375);
    u[0] = std::max(0.0,1.0-8.0*sqrt(u[0]));

    return u;
  }

  //! set time for subsequent evaluation
  void setTime (NUMBER t)
  {
    time = t;
  }

private:

  NUMBER time;
  NUMBER pi;

};
#endif //ACOUSTICS_RIEMANNPROBLEM
