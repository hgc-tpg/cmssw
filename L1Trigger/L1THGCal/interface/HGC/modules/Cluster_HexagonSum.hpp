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


  void Cluster_HexagonSum( const std::array< std::array< Cluster , 44 > ,  41 >& ClusterIn ,
                                 std::array< std::array< Cluster , 44 > ,  41 >& ClusterOut )
  {

    ClusterOut = std::array< std::array< Cluster , 44 > ,  41 >(); // Start afresh

    static const bool phase[ 8 ] = { false , true , false , true , true , false , true , false };

    for ( int j(0); j!=44 ; ++j ) { // col
      for ( int i(0); i!=40 ; ++i ) { // row
        // Add the debugging-fields to the output cell
        const Cluster& lIn = ClusterIn.at( i ).at( j );

        Cluster& lOut = ClusterOut.at( i ).at( j );
        lOut.Cell = lIn.Cell;
        lOut.DebugRow = lIn.DebugRow;
        lOut.DebugRow2 = i == 0 ? -999 : (i-1);
        lOut.DebugRow3 = i == 39 ? -999 : (i+1)%41;
        lOut.DebugCol = lIn.DebugCol;
        lOut.DebugCol2 = lIn.DebugCol2; // j == 43 ? -999 : (lIn.DebugCol+1);
        lOut.DataValid = i < 40 ? phase[ lIn.Cell%8 ] : false;
        lOut.Last = lIn.Last;

        for( int J( j ) ; J!=std::min(44,j+2) ; ++J ) {
          for( int I( std::max(0,i-1) ) ; I!=std::min(41,i+2) ; ++I ) {
            lOut += ClusterIn.at( I ).at( J );
          }
        }

      }
    }


  }

}
