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

  void Cluster_Accumulator( const std::array< std::array< Cluster , 154 > , 162 >& ClusterIn ,
                                  std::array< std::array< Cluster , 154 > , 162 >& ClusterOut )
  {
    static const int DebugColPhase = 3;
    ClusterOut = std::array< std::array< Cluster , 154 > , 162 >(); // Start afresh

    for ( int j(0); j!= 154; ++j ) { // decoder
      auto lDecoderCol = DecoderColumns[ j ];
      for ( int i(0); i!= 162; ++i ) { // clock

        {
          // Add the debugging-fields to the output cell
          Cluster& lOut = ClusterOut.at( i ).at( j );
          lOut.Cell = i;
          lOut.DebugRow = i/4;
          lOut.DebugCol = lDecoderCol + ( ( 444 + DebugColPhase + i - lDecoderCol ) % 4 ); // C++ % is actually a remainder, not a modulo - Add constant to make it positive
          lOut.DataValid = i<160;
          lOut.Last = i==161;
        }

        {
          // For each input cell, add it into the histogram
          const Cluster& lIn = ClusterIn.at( i ).at( j );
          if ( lIn.Cell > 160 ) continue;
          Cluster& lOut = ClusterOut.at( lIn.Cell ).at( j );
          lOut.Field0 += lIn.Field0;
          lOut.Field1 += lIn.Field1;
          lOut.Field2 += lIn.Field2;
          lOut.Field3 += lIn.Field3;
          lOut.Field4 += lIn.Field4;
        }

      }
    }

  }


}
