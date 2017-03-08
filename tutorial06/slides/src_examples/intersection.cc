auto it = gridview.template begin<0>();

// iterate over intersection of current entity
for (auto iit = gridview.ibegin(*it); iit != gridview.iend(*it); ++iit)
{
  // neighbor intersection
  if (iit->neighbor())
  {
      // do something ...
  }
  // boundary intersection
  if (iit->boundary())
  {
      // do something else ...
  }
}
