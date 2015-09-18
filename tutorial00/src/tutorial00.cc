// -*- tab-width: 4; indent-tabs-mode: nil -*-
/** \file
    \brief Solve Poisson equation with P1 conforming finite elements
*/
// always include the config file
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
// C++ includes
#include<math.h>
#include<iostream>
// dune-common includes
#include<dune/common/parallel/mpihelper.hh>
#include<dune/common/parametertreeparser.hh>
#include<dune/common/timer.hh>
// dune-geometry includes
#include<dune/geometry/referenceelements.hh>
#include<dune/geometry/quadraturerules.hh>
// dune-grid includes
#include<dune/grid/onedgrid.hh>
#include<dune/grid/io/file/vtk/vtkwriter.hh>
#include<dune/grid/io/file/gmshreader.hh>
#if HAVE_UG
#include<dune/grid/uggrid.hh>
#endif
#if HAVE_DUNE_ALUGRID
#include<dune/alugrid/grid.hh>
#include<dune/grid/io/file/dgfparser/dgfalu.hh>
#include<dune/grid/io/file/dgfparser/dgfparser.hh>
#endif
// dune-istl included by pdelab
// dune-pdelab includes
#include<dune/pdelab/common/function.hh>
#include<dune/pdelab/common/vtkexport.hh>
#include<dune/pdelab/finiteelementmap/pkfem.hh>
#include<dune/pdelab/constraints/common/constraints.hh>
#include<dune/pdelab/constraints/common/constraintsparameters.hh>
#include<dune/pdelab/constraints/conforming.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include<dune/pdelab/gridfunctionspace/interpolate.hh>
#include<dune/pdelab/gridfunctionspace/vtk.hh>
#include<dune/pdelab/gridoperator/gridoperator.hh>
#include<dune/pdelab/localoperator/defaultimp.hh>
#include<dune/pdelab/localoperator/pattern.hh>
#include<dune/pdelab/localoperator/flags.hh>
#include<dune/pdelab/backend/istl.hh>
#include<dune/pdelab/stationary/linearproblem.hh>

// include all components making up this tutorial
#include"ffunction.hh"
#include"gfunction.hh"
#include"bctype.hh"
#include"poissonp1.hh"
#include"driver.hh"

//===============================================================
// Main program with grid setup
//===============================================================
int main(int argc, char** argv)
{
  try{
    // Maybe initialize Mpi
    Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);
    if(Dune::MPIHelper::isFake)
      std::cout<< "This is a sequential program." << std::endl;
    else
      std::cout << "Parallel code run on "
                << helper.size() << " process(es)" << std::endl;

    // Read parameters from ini file
    Dune::ParameterTree ptree;
    Dune::ParameterTreeParser ptreeparser;
    ptreeparser.readINITree("tutorial00.ini",ptree);
    ptreeparser.readOptions(argc,argv,ptree);

    // make grid
    const unsigned int dim = ptree.get<int>("grid.dim");
    const unsigned int refinement = ptree.get<int>("grid.refinement");
    if (dim==1)
      {
        // read grid parameters from input file
        typedef Dune::OneDGrid::ctype ctype;
        const ctype a = ptree.get<ctype>("grid.oned.a");
        const ctype b = ptree.get<ctype>("grid.oned.b");
        const unsigned int N = ptree.get<int>("grid.oned.elements");

        // create equidistant intervals
        typedef std::vector<ctype> Intervals;
        Intervals intervals(N+1);
        for(unsigned int i=0; i<N+1; ++i)
          intervals[i] = a + ctype(i)*(b-a)/ctype(N);
        
        // Construct grid
        typedef Dune::OneDGrid Grid;
        Grid grid(intervals);
        grid.globalRefine(refinement);

        // call generic function
        typedef Grid::LeafGridView GV;
        GV gv = grid.leafGridView();
        driver(gv,ptree);
      }
    if (dim==2)
      {
        std::string filename = ptree.get("grid.twod.filename","unitsquare.msh");
#if HAVE_DUNE_ALUGRID
        typedef Dune::ALUGrid<2,2,Dune::simplex,Dune::nonconforming> Grid;
#elif HAVE_UG
        typedef Dune::UGGrid<2> Grid;
#else  // ! (HAVE_UG || HAVE_DUNE_ALUGRID)
        std::cout << "This example requires a simplex grid!" << std::endl;
#endif
        Dune::GridFactory<Grid> factory;
        Dune::GmshReader<Grid>::read(factory,filename,true,true);
        std::shared_ptr<Grid> gridp = std::shared_ptr<Grid>(factory.createGrid());
        Dune::Timer timer;
        gridp->globalRefine(refinement);
        std::cout << "Time for mesh refinement " << timer.elapsed() << " seconds" << std::endl;
        typedef Grid::LeafGridView GV;
        GV gv = gridp->leafGridView();
        driver(gv,ptree);
      }
    if (dim==3)
      {
        std::string filename = ptree.get("grid.threed.filename","unitcube.msh");
#if HAVE_DUNE_ALUGRID
        typedef Dune::ALUGrid<3,3,Dune::simplex,Dune::nonconforming> Grid;
#elif HAVE_UG
        typedef Dune::UGGrid<3> Grid;
#else  // ! (HAVE_UG || HAVE_DUNE_ALUGRID)
        std::cout << "This example requires a simplex grid!" << std::endl;
#endif
        Dune::GridFactory<Grid> factory;
        Dune::GmshReader<Grid>::read(factory,filename,true,true);
        std::shared_ptr<Grid> gridp = std::shared_ptr<Grid>(factory.createGrid());
        Dune::Timer timer;
        gridp->globalRefine(refinement);
        std::cout << "Time for mesh refinement " << timer.elapsed() << " seconds" << std::endl;
        typedef Grid::LeafGridView GV;
        GV gv = gridp->leafGridView();
        driver(gv,ptree);
      }
  }
  catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
	return 1;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
	return 1;
  }
}
