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

void Tower_Unpack( const std::array< std::array< LinkTower , 96 > , 162 >& TowersIn ,
                         std::array< std::array<     Tower , 96 > , 162 >& TowersOut )
{
  for ( int i(0); i!= 162; ++i ) { // clock
    for ( int j(0); j!= 96; ++j ){

      const LinkTower& lTowerIn  = TowersIn.at(i).at(j);
      Tower&           lTowerOut = TowersOut.at(i).at(j);

      lTowerOut = Tower(); // By default, set null
      lTowerOut.DataValid = lTowerIn.DataValid;
      lTowerOut.Last = lTowerIn.Last;

      // Unpack the energy
      if( lTowerIn.EcalExponent ) lTowerOut.CEE = ( lTowerIn.EcalMantissa | 0x10 ) << ( lTowerIn.EcalExponent - 1 );
      else                        lTowerOut.CEE = lTowerIn.EcalMantissa;

      if( lTowerIn.HcalExponent ) lTowerOut.CEH = ( lTowerIn.HcalMantissa | 0x10 ) << ( lTowerIn.HcalExponent - 1 );
      else                        lTowerOut.CEH = lTowerIn.HcalMantissa;


    }
  }
}

}
