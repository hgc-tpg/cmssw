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

void Links_MuxOutput( const std::array< std::array< lword , 4 > , 162 >& TowerLinksIn ,
                      const std::array< std::array< lword , 4 > , 162 >& ClusterLinksIn ,
                            std::array< std::array< lword , 4 > , 162 >& LinksOut )
{

  static uint16_t EventId = 0;
  const uint8_t TmIndex = 0;
  const uint8_t BoardId = 0;

  for ( int j( 1); j!= 31; ++j ) LinksOut.at(j) =   TowerLinksIn.at( j-1 );
  for ( int j(31); j!=162; ++j ) LinksOut.at(j) = ClusterLinksIn.at( j-1 );

  for ( int i( 0); i!=  4; ++i ){
    LinksOut.at(0).at(i).valid = true;
    LinksOut.at(0).at(i).start = true;
    LinksOut.at(0).at(i).data = ((uint64_t)0xABC0  << 48 ) |
                                ((uint64_t)     0  << 40 ) |
                                ((uint64_t)TmIndex << 32 ) |
                                ((uint64_t)BoardId << 24 ) |
                                ((uint64_t)      i << 16 ) |
                                ((uint64_t)EventId <<  0 );

    LinksOut.at(1).at(i).start = false; // We copied this from the Towers, and it is now incorrect...
    LinksOut.at(161).at(i).last = true;

    for ( int j( 0); j!=162; ++j ) LinksOut.at(j).at(i).strobe = true;

  }


  EventId++;
}

}
