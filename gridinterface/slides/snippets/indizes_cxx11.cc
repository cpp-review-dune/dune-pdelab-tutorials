#include <dune/grid/common/mcmgmapper.hh>
...

typedef Dune::SomeGrid::LeafGridView GridView;
...

/* create a mapper*/
// Layout description (equivalent to Dune::MCMGElementLayout)
template<int dim>
struct CellData {
    bool contains (Dune::GeometryType gt) {
        return gt.dim() == dim;
    }
};

// mapper for elements (codim=0) on leaf
using Mapper =
    Dune::MultipleCodimMultipleGeomTypeMapper<GridView,CellData>;
Mapper mapper(gridview);

/* setup sparsity pattern */
// iterate over the leaf
for (auto it = gridview.template begin<0>();
     it != gridview.template end<0>(); ++it)
{
    int index = mapper.map(*it);

    for (auto iit = gridview.ibegin(*it);
         iit != gridview.iend(*it);; ++iit) {
        // neighbor intersection
        if (iit->neighbor()) {
            int nindex = mapper.map(*(iit->outside()));
            matrix[index].insert(nindex);
        }
    }
}
