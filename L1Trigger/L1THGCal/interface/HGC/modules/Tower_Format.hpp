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


void Tower_Format( const std::array< std::array<          Tower , 5 > , 162 >& TowersIn ,
                         std::array< std::array< FormattedTower , 5 > , 162 >& TowersOut )
{

  const uint32_t Ratios[7] = { 0x0FFFF , 0x0EFFF , 0x0DFFF , 0x0CFFF , 0x0BFFF , 0x0AFFF , 0x09FFF };
  const uint32_t Scale = 0x1000;


  TowersOut = std::array< std::array< FormattedTower ,  5 > , 162 >(); // Start afresh

  for ( int i(0); i!= 162; ++i ) { // clock
    for ( int j(0); j!= 5; ++j ){
      const    Tower& lTowerIn = TowersIn.at(i).at(j);
      FormattedTower& lTowerOut = TowersOut.at(i).at(j);

      lTowerOut.Last      = lTowerIn.Last;
      lTowerOut.DataValid = lTowerIn.DataValid;
      lTowerOut.DEBUG_ETA = lTowerIn.DEBUG_ETA;
      lTowerOut.DEBUG_PHI = lTowerIn.DEBUG_PHI;

      lTowerOut.Ratio = 7;
      for ( int k(0); k!=7; ++k ){
        if( int( lTowerIn.CEH * Ratios[k] / std::pow( 2.0 , 17.0 ) ) <= lTowerIn.CEE ) break;
        lTowerOut.Ratio = k;
      }

      lTowerOut.Value = ( lTowerIn.CEE + lTowerIn.CEH ) * Scale / std::pow( 2.0 , 17.0 );

      if ( lTowerOut.Value > 0x3FF ) {
        lTowerOut.Value     = 0x3FF;
        lTowerOut.Saturated = true;
      } else {
        lTowerOut.Saturated = false;
      }

    }
  }

}

}
