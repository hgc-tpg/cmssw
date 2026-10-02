// -----------------------------------------------------------------------------------------------------------------------
// Andrew W. Rose, 2026
// Imperial College London Particles community
// -----------------------------------------------------------------------------------------------------------------------
#pragma once

#include <queue>
#include <array>
#include <cstdlib>

#include "L1Trigger/L1THGCal/interface/HGC/HgcTypes.hpp"
#include "L1Trigger/L1THGCal/interface/HGC/generated/PkgLutValues.hpp"

namespace HGC
{

  void SingleChannelBuffer_1( const std::array< Cluster , 162 >& ClusterIn ,
                              const std::array< int , 4 >& Thresholds ,
                                    std::array< std::vector< Cluster >  , 4 >& Buffers
                      ){

    // Clear the buffers
    for( auto& lIt : Buffers ) lIt.clear();

    // Write the event into the buffers
    for ( int i(0); i!= 162; ++i ) {
      Cluster lOut = ClusterIn.at(i);
      if( !lOut.DataValid ) continue;
      for ( int j(0); j!=4; ++j ) {
        if (lOut.Field0.E >= static_cast<uint32_t>(Thresholds.at(j))) {
          lOut.DebugCol = lOut.DebugCol2 = lOut.DebugRow = lOut.DebugRow2 = lOut.DebugRow3 = -999;
          lOut.DataValid = false;
          Buffers.at(j).push_back( lOut );
          break;
        }
      }
    }

  }


  void SingleChannelBuffer_2( const std::array< std::vector< Cluster >  , 4 >& Buffers ,
                              const std::array< int  , 162 >& Bin ,
                              const std::array< int  , 162 >& Index ,
                                    std::array< Cluster , 162 >& ClustersOut
                      ){

    for ( int i(0); i!= 162; ++i ) {
      auto& lBuf = Buffers.at( Bin.at(i) );
      if( !lBuf.size() ) continue;
      ClustersOut.at(i) = lBuf.at( Index.at(i) );
    }

  }





  void Cluster_Buffer2_1( const std::array< std::vector< Cluster >  , 4 >& Buffers1 ,
                          const std::array< std::vector< Cluster >  , 4 >& Buffers2 ,
                          const int& aEvent ,
                                std::array< int , 162 >& Ram ,
                                std::array< int , 162 >& Event ,
                                std::array< int , 162 >& Bin ,
                                std::array< int , 162 >& Index ,
                                std::array< bool , 162 >& BinValid
                      ){

    int r = 0;
    int b = 0;
    int i = 0;
    bool v = true;

    for ( int idx(0); idx!=162 ; ++idx ){

      Ram.at(idx) = r;
      Event.at(idx) = aEvent;
      Bin.at(idx) = b;
      Index.at(idx) = i;
      BinValid.at(idx) = v;

      const auto& lBuf = r ? Buffers2.at( b ) : Buffers1.at( b );

      if (lBuf.empty() || static_cast<size_t>(i) == lBuf.size() - 1) {
        if ( r == 1 ) {
          r = 0;
          i = 0;

          if ( b == 3 ) v = false;
          b = (b+1)%4;

        } else {
          r = 1;
          i = 0;
        }
      } else {
        i = i + 1;
      }


    }

  }




  void Cluster_Buffer2_2( const std::array< std::vector< Cluster >  , 4 >& Buffers1 ,
                          const std::array< std::vector< Cluster >  , 4 >& Buffers2 ,
                          const std::array< int , 162 >& Ram ,
                          const std::array< int , 162 >& Event ,
                          const std::array< int , 162 >& Bin ,
                          const std::array< int , 162 >& Index ,
                          const std::array< bool , 162 >& BinValid ,
                                std::array< Cluster , 162 >& ClusterOut
                      ){

      for( int i(0); i!=162; ++i ){
        const auto& I = Index.at(i);
        const auto& b = Bin.at(i);
        const auto& r = Ram.at(i);
        const auto& lBuf = r ? Buffers2.at( b ) : Buffers1.at( b ) ;

        if (BinValid.at(i) && !(lBuf.empty() || I < 0 || static_cast<size_t>(I) >= lBuf.size())) {
          ClusterOut.at(i) = lBuf.at( I );
          ClusterOut.at(i).DataValid = true;
        } else {
          ClusterOut.at(i) = Cluster();
        }

      }
      ClusterOut.at(161).Last = true;

  }


  __attribute__((flatten))
  void Cluster_Buffer2( const std::array< Cluster , 162 >& ClustersIn1 ,
                        const std::array< Cluster , 162 >& ClustersIn2 ,
                              std::array< Cluster , 162 >& ClustersOut ){

    const std::array< int , 4 > Thresholds1 = { 3000 , 1000 , 500 , 1 };
    const std::array< int , 4 > Thresholds2 = { 3000 , 1000 , 500 , 1 };

    // CMSSW supplies one complete packet per event. Keep this packet-local so
    // independent framework streams cannot share mutable firmware state.
    uint32_t EventFlag = 0;
    std::array< std::vector< Cluster >  , 4 > Buffers1 , Buffers2;
    std::array< int , 162 > RamTemp , EventTemp , BinTemp , IndexTemp;
    std::array< bool , 162 > BinValidTemp;

    HGC::SingleChannelBuffer_1( ClustersIn1 , Thresholds1 , Buffers1 );
    HGC::SingleChannelBuffer_1( ClustersIn2 , Thresholds2 , Buffers2 );
    HGC::Cluster_Buffer2_1( Buffers1 , Buffers2 , EventFlag , RamTemp , EventTemp , BinTemp , IndexTemp , BinValidTemp );
    HGC::Cluster_Buffer2_2( Buffers1 , Buffers2 , RamTemp , EventTemp , BinTemp , IndexTemp , BinValidTemp , ClustersOut );
  }


}
