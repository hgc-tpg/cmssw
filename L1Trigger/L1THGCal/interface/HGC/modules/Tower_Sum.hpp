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

  const uint16_t FibreMap[6][16] = { { 0 , 1 , 2 , 3 , 4 , 5 , 6 , 7 , 8 , 9 , 10 , 11 , 12 , 13 , 14 , 15 } , // Fibre 0 {Secondary}
                                     {16 , 17 , 18 , 19 , 20 , 21 , 22 , 23 , 24 , 25 , 26 , 27 , 28 , 29 , 30 , 31 } , // Fibre 1 {Secondary}
                                     {32 , 33 , 34 , 35 , 36 , 37 , 38 , 39 , 40 , 41 , 42 , 43 , 44 , 45 , 46 , 47 } , // Fibre 2 {Primary}
                                     {48 , 49 , 50 , 51 , 52 , 53 , 54 , 55 , 56 , 57 , 58 , 59 , 60 , 61 , 62 , 63 } , // Fibre 3 {Primary}
                                     {64 , 65 , 66 , 67 , 68 , 69 , 70 , 71 , 72 , 73 , 74 , 75 , 76 , 77 , 78 , 79 } , // Fibre 4 {Primary}
                                     {80 , 81 , 82 , 83 , 84 , 85 , 86 , 87 , 88 , 89 , 90 , 91 , 92 , 93 , 94 , 95 } }; // Fibre 5 {Primary}



void Tower_Sum( const std::array< std::array< Tower , 96 > , 162 >& TowersIn ,
                      std::array< std::array< Tower ,  5 > , 162 >& TowersOut )
{

  TowersOut = std::array< std::array< Tower ,  5 > , 162 >(); // Start afresh

  for ( int j(0); j!= 6; ++j ){
    int cnt = 0;
    for ( int i(0); i!= 162; ++i ) { // clock
      Tower& lTowerOut = TowersOut.at(i).at(j%5);

      lTowerOut.DEBUG_ETA = cnt % 20;
      lTowerOut.DEBUG_PHI = ( cnt / 20 ) + ( 5 * j ) + (  j > 1 ? -10 : 15 );
      lTowerOut.Last      = ( i == 161 );

      if( TowersIn.at(i).at( FibreMap[j][0] ).DataValid ){
        lTowerOut.DataValid = ( cnt < 100 );
        cnt += 1;
      }

      for ( int k(0); k!= 16; ++k ){
        const Tower& lTowerIn = TowersIn.at(i).at( FibreMap[j][k] );
        lTowerOut.CEE += lTowerIn.CEE;
        lTowerOut.CEH += lTowerIn.CEH;
      }

    }
  }
}

}
