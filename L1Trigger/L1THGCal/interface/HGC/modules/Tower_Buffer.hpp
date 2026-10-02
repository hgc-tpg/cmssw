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

void Tower_Buffer( const std::array< std::array< FormattedTower ,  5 > , 162 >& TowersIn ,
                         std::array< std::array< FormattedTower , 24 > ,  20 >& TowersOut ,
                         std::array< std::array< lword          ,  4 > , 162 >& TowerLinksOut )
{

  TowerLinksOut = std::array< std::array< lword          ,  4 > , 162 >(); // Start afresh

  uint32_t k(0);
  for ( int i(0); i!= 162; ++i ) { // clock

    for ( int j(0); j!= 4; ++j ){
      auto& lOut = TowerLinksOut.at(i).at(j);
      lOut.strobe = true;
      lOut.valid  = (i<30);
      lOut.start  = (i==0);
      lOut.last   = (i==161);
    }

    if( (i%9) > 5 ) continue;

    for ( int j(0); j!= 5; ++j ){
      const auto& lTower = TowersIn.at(i).at(j);

      uint32_t channel = lTower.DEBUG_PHI/6;

      if( channel > 3 ) continue;

      uint32_t clk = (lTower.DEBUG_ETA/4) + 5*(lTower.DEBUG_PHI%6);
      uint32_t offset = lTower.DEBUG_ETA%4;

      auto& lOut = TowerLinksOut.at( clk ).at( channel );

      if( lOut.data & ( ( (uint64_t)0xFFFF ) << (16*offset) ) ){ // Should be able to figure out a better way than checking for existing data...
        // std::cout << i << " " << j << " : " << clk << " " << channel << " " << offset << " : " << lTower.DEBUG_PHI << " " << lTower.DEBUG_ETA << std::endl;
        continue;
      }

      lOut.data |= ( ( (((uint64_t)lTower.Saturated)<<13) | (lTower.Ratio<<10) | (lTower.Value<<0) ) << (16*offset) );

      TowersOut.at( lTower.DEBUG_ETA ).at( lTower.DEBUG_PHI ) = lTower;

    }

  }

}

}
