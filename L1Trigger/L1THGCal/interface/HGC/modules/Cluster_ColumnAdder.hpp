// -----------------------------------------------------------------------------------------------------------------------
// Andrew W. Rose, 2026
// Imperial College London Particles community
// -----------------------------------------------------------------------------------------------------------------------
#pragma once

#include <array>
#include <cstdlib>

#include "L1Trigger/L1THGCal/interface/HGC/HgcTypes.hpp"
#include "L1Trigger/L1THGCal/interface/HGC/generated/PkgLutValues.hpp"

namespace HGC
{


  void Cluster_ColumnAdder( const std::array< std::array< Cluster , 154 > , 162 >& ClusterIn ,
                                  std::array< std::array< Cluster ,  44 > ,  41 >& ClusterOut )
  {

    ClusterOut = std::array< std::array< Cluster , 44 > , 41 >(); // Start afresh

    for ( int j(0); j!=44 ; ++j ) { // col
      for ( int i(0); i!=41 ; ++i ) { // row
        // Add the debugging-fields to the output cell
        Cluster& lOut = ClusterOut.at( i ).at( j );
        lOut.Cell = (4*i) + (j%4);
        lOut.DebugRow = i;
        lOut.DebugCol = j - 21; // No negative indices in C++
        lOut.DataValid = (i<40);
        lOut.Last = (lOut.Cell==161);
      }
    }

    for ( int j(0); j!= 154; ++j ) { // decoder
      for ( int i(0); i!= 162; ++i ) { // clock
        const Cluster& lIn = ClusterIn.at( i ).at( j );
        Cluster& lOut = ClusterOut.at( lIn.DebugRow ).at( lIn.DebugCol + 21 ); // No negative indices in C++
        lOut += lIn;
      }
    }

  }

}
