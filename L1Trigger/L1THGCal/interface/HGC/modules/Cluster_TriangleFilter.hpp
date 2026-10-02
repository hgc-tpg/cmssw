// -----------------------------------------------------------------------------------------------------------------------
// Andrew W. Rose, 2026
// Imperial College London Particles community
// -----------------------------------------------------------------------------------------------------------------------
#pragma once

#include <array>
#include <cstdlib>

#include "L1Trigger/L1THGCal/interface/HGC/HgcTypes.hpp"

namespace HGC
{

  void Cluster_TriangleFilter( const std::array< std::array< Cluster , 44 > ,  41 >& HexagonsIn ,
                               const std::array< std::array< Cluster , 44 > ,  41 >& TrianglesIn ,
                                     std::array< std::array< Cluster , 44 > ,  41 >& TrianglesOut )
  {
    TrianglesOut = TrianglesIn;

    for( auto& lIt2 : HexagonsIn )
    {
      for( auto& lIt : lIt2 )
      {
        if( lIt.DataValid ){
          if( lIt.DebugCol >= -21 and lIt.DebugCol <= 22 ) {
            if( lIt.DebugRow  != -999 ) TrianglesOut.at( lIt.DebugRow  ).at( lIt.DebugCol+21  ).DataValid = false;
            if( lIt.DebugRow2 != -999 ) TrianglesOut.at( lIt.DebugRow2 ).at( lIt.DebugCol+21  ).DataValid = false;
            if( lIt.DebugRow3 != -999 ) TrianglesOut.at( lIt.DebugRow3 ).at( lIt.DebugCol+21  ).DataValid = false;
          }
          if( lIt.DebugCol2 >= -21 and lIt.DebugCol2 <= 22 ) {
            if( lIt.DebugRow  != -999 ) TrianglesOut.at( lIt.DebugRow  ).at( lIt.DebugCol2+21 ).DataValid = false;
            if( lIt.DebugRow2 != -999 ) TrianglesOut.at( lIt.DebugRow2 ).at( lIt.DebugCol2+21 ).DataValid = false;
            if( lIt.DebugRow3 != -999 ) TrianglesOut.at( lIt.DebugRow3 ).at( lIt.DebugCol2+21 ).DataValid = false;
          }
        }
      }
    }

    for( auto& lIt2 : TrianglesOut )
    {
      for( auto& lIt : lIt2 )
      {
        if( !lIt.DataValid ){
          lIt.Field0 = Cluster::tField0();
          lIt.Field1 = Cluster::tField1();
          lIt.Field2 = Cluster::tField2();
          lIt.Field3 = Cluster::tField3();
          lIt.Field4 = Cluster::tField4();
        }
      }
    }

  }

}
