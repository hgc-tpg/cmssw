// -----------------------------------------------------------------------------------------------------------------------
// Andrew W. Rose, 2026
// Imperial College London Particles community
// -----------------------------------------------------------------------------------------------------------------------
#pragma once

#include <vector>
#include <array>
#include <cstdlib>

#include "L1Trigger/L1THGCal/interface/HGC/HgcTypes.hpp"
#include "L1Trigger/L1THGCal/interface/HGC/generated/PkgLutValues.hpp"

namespace HGC
{

  struct ClusterFifo
  {
    static const int size = 8;
    std::array< Cluster , size > RAM;
    uint32_t ReadAddr, WriteAddr;

    ClusterFifo() : ReadAddr(0) , WriteAddr(0)
    {
      for( auto& lIt: RAM ) lIt.Cell = 0;
    }

    void push( const Cluster& aCluster , const uint8_t& aEventFlag )
    {
      auto& lCell = RAM.at( WriteAddr );
      lCell = aCluster;
      lCell.DebugCol = lCell.DebugCol2 = lCell.DebugRow = lCell.DebugRow2 = lCell.DebugRow3 = -999;
      lCell.Last = 0;
      lCell.Field0.EventFlag = aEventFlag;

      WriteAddr = ( WriteAddr + 1 ) % size;
    }

    Cluster& pop()
    {
      Cluster& lRet = RAM.at( ReadAddr );
      if( ReadAddr == WriteAddr ) lRet.DataValid = false;
      else ReadAddr = ( ReadAddr + 1 ) % size;
      return lRet;
    }

  };


  void Cluster_Funnel1( const std::array< std::array< Cluster , 11 > , 162 >& ClusterIn ,
                             std::array< std::array< Cluster , 11 > , 162 >& ClusterReg ,
                             std::array< uint16_t , 162 >& Sel )
  {

    // -----------------------------------------------------------
    // CMSSW supplies one complete packet per event. Keep this packet-local so
    // independent framework streams cannot share mutable firmware state.
    std::array< ClusterFifo , 11 > lFifos;
    std::array< Cluster , 11 > ClusterRegInt;
    uint16_t SelInt(0);
    uint8_t EventFlag = 1;
    // -----------------------------------------------------------

    for ( int i(0); i!= 163; ++i ) { // clock

      if( i<162 ) {

        int j(0);
        for ( ; j!=11; ++j ) {
          SelInt = (SelInt+1)%11;
          if( ClusterRegInt.at( SelInt ).DataValid ) break;
        }
        if( j==11 ) SelInt = 15;

        for ( int j(0); j!= 11; ++j ) { // channel

          Cluster& lOut = ClusterRegInt.at(j);

          // Clock 1
          if( lOut.DataValid ) {
            if( j == SelInt ) lOut.DataValid = false;
          } else {
            lOut = lFifos.at(j).pop();
          }

          // Clock 0
          const Cluster& lIn = ClusterIn.at(i).at(j);
          if( lIn.DataValid ) lFifos.at(j).push( lIn , EventFlag );
        }
      }

      if( i > 0 ){
        ClusterReg.at(i-1) = ClusterRegInt;
        Sel.at(i-1) = SelInt;
      }

    }

    EventFlag = ( EventFlag + 1 ) % 4;
  }



  void Cluster_Funnel2( const std::array< std::array< Cluster , 11 > , 162 >& ClusterReg ,
                        const std::array< uint16_t , 162 >& Sel ,
                              std::array< Cluster  , 162 >& ClusterOut )
  {

    // ClusterOut = std::array< std::array< Cluster , 1 >  , 162 >(); // Start afresh

    for ( int i(0); i!= 162; ++i ) { // clock
      auto lSel = Sel.at( i );

      if( lSel < 11 )
      {
        ClusterOut.at(i) = ClusterReg.at(i).at( lSel );
        ClusterOut.at(i).DataValid = true;
      } else {
        ClusterOut.at(i) = Cluster();

        // std::cout << i << " Sel:" << lSel << " | ";
        // for( auto& x : ClusterReg.at(i) ) std::cout << (int)(x.Cell) << " " << x.DataValid << ", ";
        // std::cout << std::endl;

        // ClusterOut.at(i) = ClusterReg.at(i).at( 7 );
        // ClusterOut.at(i).DataValid = false;
      }

      ClusterOut.at(i).Last = (i==161);
    }

  }



  void Cluster_Funnel( const std::array< std::array< Cluster , 44 > ,  41 >& ClusterIn ,
                             std::array< Cluster  , 162 >& ClusterOut )
  {
    std::array< std::array< Cluster , 11 > , 162 > ClusterReg;
    std::array< uint16_t , 162 > Sel;

    std::array< std::array< Cluster , 11 > , 162 > lFirmware;
    PhysicalToFirmware( ClusterIn , lFirmware );
    Cluster_Funnel1( lFirmware , ClusterReg , Sel );
    Cluster_Funnel2( ClusterReg , Sel , ClusterOut );
  }


}
