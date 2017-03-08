const auto &is = *iit;

// evaluate fluxes
Dune::FieldVector<ctype, dim> center = is.geometry().center();
if (is.neighbor())
{
    // mean flux
    flux = ( myshapefkt.gradient(center)
           - nbshapefkt.gradient(center) )
           * is.centerUnitOuterNormal()
           * is.geometry().volume();
}
// boundary intersection
else if (is.boundary())
{
    // neumann boundary condition
    flux = j(center);
}
