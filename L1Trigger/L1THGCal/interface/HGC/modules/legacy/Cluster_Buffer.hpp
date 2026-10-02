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


  void Cluster_Buffer1( const int aResetOffset ,
                        const std::array< std::array< Cluster , 1 > , 162 >& ClusterIn ,
                              std::array< std::array< uint16_t , 8 > , 162 >& Counts ,
                              std::array< uint16_t , 162 >& BufferWaddr ,
                              // std::array< uint16_t , 162 >& BufferRaddr,
                              std::array< Cluster  , 512 >& Buffer
                      ){

    const std::array< uint32_t , 4 > Constants = { 3000 , 1000 , 500 , 1 };

    // -----------------------------------------------------------
    // THIS NEED TO PERSIST BETWEEN EVENTS FOR FIRMWARE EMULATION
    static bool EventFlag = 1;
    // -----------------------------------------------------------

    int ReadBin = 4;
    int ReadCtr = 0;

    for ( int i(0); i!= 162; ++i ) { // clock

      auto& lCluster = ClusterIn.at(i).at(0);

      // ============================================
      Counts.at(i) = Counts.at( (i+161)%162 );
      // ============================================

      // ============================================
      std::array< bool , 4 > FavouredBin, CountOk;
      for ( int j(0); j!=4; ++j ){
        auto index = (4*lCluster.Field0.EventFlag) + j;
        CountOk.at(j)     = ( Counts.at( (i+160)%162 ).at(index) < 61 );
        FavouredBin.at(j) = ( lCluster.Field0.E >= Constants.at(j) );
      }

      int Selected = 0;
      for ( ; Selected!=4; ++Selected ){
        if( FavouredBin.at(Selected) and CountOk.at(Selected) ) break;
      }

      if( Selected==4 ) {
        for ( Selected=3 ; Selected>=0; --Selected ){
          if( CountOk.at(Selected) ) break;
        }
      }
      // ============================================

      // ============================================
      {
        auto index = (4*lCluster.Field0.EventFlag) + Selected;
        BufferWaddr.at(i) = (256*lCluster.Field0.EventFlag) + (64*Selected) + Counts.at(i).at(index);
        Buffer.at( BufferWaddr.at(i) ) = lCluster;

        if( lCluster.DataValid ) Counts.at(i).at(index) += 1;

        if( i == aResetOffset )
        {
          for ( int j(0); j!=4; ++j ){
            auto index = (4*EventFlag) + j;
            Counts.at( i ).at(index) = 0;
          }
        }
      }
      // ============================================
    }

    EventFlag = !EventFlag;
  }




  void Cluster_Buffer2( const std::array< std::array< uint16_t , 8 > , 162 >& Counts ,
                        const std::array< Cluster  , 512 >& Buffer ,
                              std::array< uint16_t , 162 >& BufferRaddr,
                              std::array< std::array< Cluster , 1 > , 162 >& ClusterOut
                      ){


    // -----------------------------------------------------------
    // THESE NEED TO PERSIST BETWEEN EVENTS FOR FIRMWARE EMULATION
    static int ReadPage = 1;
    // -----------------------------------------------------------

    int ReadBin = 0;
    int ReadCtr = 0;
    ReadPage = !ReadPage;

    for ( int i(0); i!= 162; ++i ) { // clock
      BufferRaddr.at(i) = ( (256*ReadPage) + (64*ReadBin) + ReadCtr ) % 512;
      int index = ( 4 * ReadPage ) + ReadBin;

      Cluster& lOut = ClusterOut.at( i ).at( 0 );

      lOut = Buffer.at( BufferRaddr.at(i) );
      lOut.DebugCol = lOut.DebugCol2 = lOut.DebugRow = lOut.DebugRow2 = lOut.DebugRow3 = -999;
      lOut.Last = ( i == 161 );

      if( ReadBin == 4 )
      {
        lOut.DataValid = false;
      }
      else if( ReadCtr < Counts.at(i).at(index) )
      {
        lOut.DataValid = true;
        ReadCtr += 1;
      }
      else
      {
        lOut.DataValid = false;
        ReadBin += 1;
        ReadCtr = 0;
      }
    }

  }


}
