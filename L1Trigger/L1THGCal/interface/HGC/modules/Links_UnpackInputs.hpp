// -----------------------------------------------------------------------------------------------------------------------
// Andrew W. Rose, 2026
// Imperial College London Particles community
// -----------------------------------------------------------------------------------------------------------------------
#pragma once

#include <array>

#include <cstdlib>
#include <ctime>

#include "L1Trigger/L1THGCal/interface/HGC/HgcTypes.hpp"

namespace HGC
{

void Links_UnpackInput( const std::array< std::array< lword           , 96   > , 162 >& linksIn ,
                              std::array< std::array< LinkTriggerCell , 96*2 > , 162 >& CellsOut ,
                              std::array< std::array< LinkTower       , 96   > , 162 >& TowersOut ,
                              std::array< std::array< LinkFlags       , 96   > , 162 >& FlagsOut )
{
    uint64_t lCell;

    for ( int i(0); i!= 162; ++i ) { // clock
      for ( int j(0); j!= 96; ++j ) { // channel

        // Trigger Cells
        for ( int k(0); k!= 2; ++k ) { // cell
          if( i%9 < 6 )     lCell = linksIn.at(i)  .at(j).data >> (15*k); // We mask below
          else if( k == 0 ) lCell = linksIn.at(i-6).at(j).data >> (15*2); // We mask below
          else /* k == 1 */ lCell = linksIn.at(i-3).at(j).data >> (15*2); // We mask below
          CellsOut.at(i).at((2*j)+k) = { /*Mantissa*/ (lCell >> 0 ) & 0xF , /*Exponent*/ (lCell >> 4) & 0x1F , /*TcId*/ (lCell >> 9) & 0x3F , /*Last*/ i==161 , /*DataValid*/ true };
        }

        const uint64_t& lWord  = linksIn.at(i).at(j).data; // We mask below
        const bool&     lValid = linksIn.at(i).at(j).valid; // <================ This should probably be "strobe"

        // Towers
        TowersOut.at(i).at(j) = { /*EcalExponent*/ (lWord >> 49) & 0xF , /*EcalMantissa*/ (lWord >> 45) & 0xF , /*HcalExponent*/ (lWord >> 57) & 0xF , /*HcalMantissa*/ (lWord >> 53) & 0xF , /*Last*/ i==161 , /*DataValid*/ lValid };

        // Flags
        FlagsOut.at(i).at(j) = { /*Data*/ (lWord >> 61) & 0x7 , /*DataValid*/ lValid , /*Last*/ i==161 };

      }
    }

}

}