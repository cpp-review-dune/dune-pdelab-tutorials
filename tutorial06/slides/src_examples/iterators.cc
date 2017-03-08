#include <dune/grid/somegrid.hh>

void iterate_the_grid()
{
    // we have a grid
    Dune::SomeGrid grid(parameters);

    // iterate over level 2
    auto levelGV = grid.levelGridView(2);
    for (auto it = levelGV.template begin<0>(); it != levelGV.template end<0>(); ++it) {
        ...
    }

    // iterate over the leaf
    auto leafGV = grid.leafGridView();
    for (auto it = leafGV.template begin<0>(); it != leafGV.template end(); ++it) {
        ...
    }
}
