#ifndef MAXWELL_RIEMANNPROBLEM
#define MAXWELL_RIEMANNPROBLEM
template<typename GV, typename NUMBER>
class Problem
{
public:

  using RangeField = NUMBER;

  //problem specification depends on dimension
  static constexpr int dim = 3;
  static constexpr int m = 6;

  using Range = Dune::FieldVector<NUMBER,m>;

  Problem ()
    : time(0.0),  pi(3.141592653589793238462643)
  {
  }

  //! permittivity
  template<typename E, typename X>
  NUMBER eps (const E& e, const X& x) const
  {
    return 1.0;
  }

  //! permeability
  template<typename E, typename X>
  NUMBER mu (const E& e, const X& x) const
  {
    return 1.0;
  }

  //! speed of sound
  template<typename E, typename X>
  NUMBER c (const E& e, const X& x) const
  {
    return 1.0;
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
    // u[0] = -s[0];
    // u[1] = -s[1];
    // u[2] = -s[2];
    // u[3] = -s[3];
    // u[4] = -s[4];
    // u[5] = -s[5];
    return u;
  }

  //! right hand side
  template<typename E, typename X>
  Range q (const E& e, const X& x) const
  {
    Range rhs(0.0);
    return rhs;
  }

  //! initial value
  template<typename E, typename X>
  Range u0 (const E& e, const X& p) const
  {
    X xglobal = e.geometry().global(p);


    Range u(0.0);

    auto x=xglobal[0], y=xglobal[1], z=xglobal[2];
    /*
    if (x>0.4 && x<0.6 && y>0.4 && y<0.6 && z>0.4 && z<0.6)
      {
        auto p1 = sin(pi*(x-0.4)/0.2)*sin(pi*(x-0.4)/0.2);
        auto p1prime = 2.0*sin(pi*(x-0.4)/0.2)*cos(pi*(x-0.4)/0.2)*pi/0.2;
        auto q2 = sin(pi*(y-0.4)/0.2)*sin(pi*(y-0.4)/0.2);
        auto q2prime = 2.0*sin(pi*(y-0.4)/0.2)*cos(pi*(y-0.4)/0.2)*pi/0.2;
        auto r3 = sin(pi*(z-0.4)/0.2)*sin(pi*(z-0.4)/0.2);
        auto r3prime = 2.0*sin(pi*(z-0.4)/0.2)*cos(pi*(z-0.4)/0.2)*pi/0.2;
        u[0] = -p1*q2prime*r3prime;
        u[1] = -0.5*p1prime*q2*r3prime;
        u[2] = -0.5*p1prime*q2prime*r3;
      }
    */

    auto alpha = 1.0;

    auto s1 = std::sin(alpha*x);
    auto c1 = std::cos(alpha*x);
    auto st = std::sin(alpha*0.0);
    auto ct = std::cos(alpha*0.0);
    //u[0] += 0;     // E_x
    u[1] += -c1*st;  // E_y
    u[2] += s1*ct;   // E_z
    //u[3] += 0;     // H_x
    u[4] += c1*st;   // H_y
    u[5] += s1*ct;   // H_z

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
#endif //MAXWELL_RIEMANNPROBLEM
