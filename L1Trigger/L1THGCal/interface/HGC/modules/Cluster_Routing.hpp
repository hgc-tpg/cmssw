// -----------------------------------------------------------------------------------------------------------------------
// Andrew W. Rose, 2026
// Imperial College London Particles community
// -----------------------------------------------------------------------------------------------------------------------
#pragma once

#include <array>

#include <cstdlib>
#include <ctime>

#include "L1Trigger/L1THGCal/interface/HGC/HgcTypes.hpp"
#include "L1Trigger/L1THGCal/interface/HGC/generated/PkgLutValues.hpp"

namespace HGC
{

void Cluster_Routing( const std::array< std::array< LinkTriggerCell , 192 > , 162 >& CellsIn ,
                            std::array< std::array< LinkTriggerCell , 289 > , 162 >& CellsInt ,
                            std::array< std::array< LinkTriggerCell , 154 > , 162 >& CellsOut )
{
  for ( int i(0); i!= 162; ++i ) { // clock
    for ( int j(0); j!= 289; ++j ) CellsInt.at(i).at(j) = CellsIn .at(i).at( RoutingLut1[j][i] );
    for ( int j(0); j!= 154; ++j ) CellsOut.at(i).at(j) = CellsInt.at(i).at( RoutingLut2[j][i] );
  }
}

}
