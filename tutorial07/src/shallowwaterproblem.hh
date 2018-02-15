#ifndef SHALLOWWATER_RIEMANNPROBLEM
#define SHALLOWWATER_RIEMANNPROBLEM
template<const int dim, typename GV, typename NUMBER>
class Problem
{
public:

  using RangeField = NUMBER;

  //problem specification depends on dimension
  static constexpr int m = dim+1;

  using Range = Dune::FieldVector<NUMBER,m>;

  Problem ()
    : time(0.0),  pi(3.141592653589793238462643)
  {
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
    for (int i=0; i<dim; i++)
       u[i+1] = -s[i+1];

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

    //hump
    auto tmp = 0.0;
    for (int i=0; i<dim; i++){
      tmp += std::pow((xglobal[i]-.5)/.2,2)/2;
    }

    u[0] += exp(-tmp) + .5;

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
#endif //SHALLOWWATER_RIEMANNPROBLEM
