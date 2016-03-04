#ifndef COMMUNICATE_HH
#define COMMUNICATE_HH


template <class IndexSetImp,
    class GlobalIdSetImp,
    class DataVectorType >
class ExampleDataHandle
  : public Dune::CommDataHandleIF< ExampleDataHandle< IndexSetImp, GlobalIdSetImp, DataVectorType >, typename DataVectorType::value_type >
{
  const IndexSetImp & iset_;
  const GlobalIdSetImp & ids_;
  int cdim_;
  DataVectorType & data1_;
public:
  typedef typename DataVectorType :: value_type DataType;
  ExampleDataHandle(const IndexSetImp & iset,
                    const GlobalIdSetImp & ids,
                    int cdim,
                    DataVectorType & d1) :
    iset_(iset), ids_(ids) , cdim_(cdim), data1_(d1)
  {}

  //! returns true if data for this codim should be communicated
  bool contains (int dim, int codim) const
  {
    return (codim==cdim_);
  }

  //! returns true if size per entity of given dim and codim is a constant
  bool fixedsize (int dim, int codim) const
  {
    // this problem is a fixed size problem,
    // but to simulate also non-fixed size problems
    // we set this to false, should work anyway
    return false;
  }

  /*! how many objects of type DataType have to be sent for a given entity
     Note: Only the sender side needs to know this size.
   */
  template<class EntityType>
  size_t size (EntityType& e) const
  {
    typedef typename EntityType::Geometry Geometry;
    return Geometry::coorddimension;
  }

  //! pack data from user to message buffer
  template<class MessageBuffer, class EntityType>
  void gather (MessageBuffer& buff, const EntityType& e) const
  {
    // std::cout << "palpo mydim: " << e.mydimension << std::endl;
    // std::cout << "palpo dim: " << e.dimension << std::endl;

    int idx = iset_.index(e);
    buff.write(data1_[idx]);
  }

  /*! unpack data from message buffer to user
     n is the number of objects sent by the sender
   */
  template<class MessageBuffer, class EntityType>
  void scatter (MessageBuffer& buff, const EntityType& e, size_t n)
  {
    // as this problem is a fixed size problem we can check the sizes
    assert( n == size(e) );

    // else do normal scatter
    int idx = iset_.index(e);
    DataType x=0.0;
    buff.read(x);
    data1_[idx] = x;
  }
};


template <typename GV>
void communicate(const GV& gv){
    typedef typename GV::IndexSet IndexSet;
    typedef typename GV::Traits::Grid::Traits::GlobalIdSet GlobalIdSet;
    typedef std::vector<double> ArrayType;

    const IndexSet& indexSet = gv.indexSet();
    const GlobalIdSet& globalIdSet(gv.grid().globalIdSet());
    int cdim = 0;
    const int dataSize = indexSet.size( cdim );
    ArrayType data(dataSize, 0.0);

    // set initial data
    for(int i=0 ; i<dataSize; ++i)
    {
      data[i] = gv.comm().rank();
    }


    using DH = ExampleDataHandle<IndexSet, GlobalIdSet, ArrayType>;
    DH dh(indexSet, globalIdSet, cdim, data);

    int myrank = gv.comm().rank();
    if( myrank == 0 )
    {
      std::cout << "TEST ";
      std::cout << " communication for codim " << cdim << std::endl;
    }

    gv.communicate(dh,Dune::All_All_Interface,Dune::ForwardCommunication);
    // gv.communicate(dh,Dune::InteriorBorder_All_Interface,Dune::ForwardCommunication);
    // gv.communicate(dh,Dune::InteriorBorder_InteriorBorder_Interface,Dune::ForwardCommunication);
    int sum(0);
    int total_data(0);
    if (myrank==0){
      for(int i=0 ; i<dataSize; ++i)
      {
        total_data += 1;
        sum += data[i];
        std::cout << data[i] << " ";
      }
      std::cout << std::endl;
      std::cout << "sum: " << sum << std::endl;
      std::cout << "total_data: " << total_data << std::endl;
    }
}

#endif
