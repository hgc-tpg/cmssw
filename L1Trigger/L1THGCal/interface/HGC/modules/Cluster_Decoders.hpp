// -----------------------------------------------------------------------------------------------------------------------
// Andrew W. Rose, 2026
// Imperial College London Particles community
// -----------------------------------------------------------------------------------------------------------------------
#pragma once

#include <array>
#include <cstdlib>

#include "L1Trigger/L1THGCal/interface/HGC/HgcTypes.hpp"
#include "L1Trigger/L1THGCal/interface/HGC/generated/PkgLutValues.hpp"

namespace HGC
{

  void Cluster_Decode(const LinkTriggerCell& lCell, const TcDecoder& lDecoded, Cluster& lCluster) {
    lCluster = Cluster();
    lCluster.Last = lCell.Last;
    lCluster.Cell = lDecoded.Cell;
    if (!lCell.DataValid || !lDecoded.DataValid)
      return;

    lCluster.DataValid = true;
    lCluster.Field0.ShapeQ = 1;
    lCluster.Field1.N_TC = 1;
    lCluster.Field2.N_TC_W = 1;

    uint64_t E = lCell.Mantissa;
    if (lCell.Exponent)
      E = (E | 0x10) << lCell.Exponent;

    lCluster.Field0.E = (E >> 1) & 0x0003FFFF;
    lCluster.Field0.xE = (E >> 19);

    lCluster.Field0.W = (E >> 2) & 0x0000FFFF;
    if (lCluster.Field0.W == 0)
      lCluster.Field0.W = 1;

    lCluster.Field1.WZ = (lCluster.Field0.W * lDecoded.LayerDepth) >> 4;
    lCluster.Field2.Wroz = (lCluster.Field0.W * lDecoded.RoverZ) >> 4;
    lCluster.Field3.Wphi = (lCluster.Field0.W * lDecoded.Phi) >> 4;

    lCluster.Field0.W2 = ((uint64_t)lCluster.Field0.W * lCluster.Field0.W);
    lCluster.Field1.WZ2 = ((uint64_t)lCluster.Field1.WZ * lDecoded.LayerDepth) >> 4;
    lCluster.Field2.Wroz2 = ((uint64_t)lCluster.Field2.Wroz * lDecoded.RoverZ) >> 4;
    lCluster.Field3.Wphi2 = ((uint64_t)lCluster.Field3.Wphi * lDecoded.Phi) >> 4;

    uint64_t Eem = ((uint64_t)lCluster.Field0.E * (uint64_t)lDecoded.LayerWeight) + (uint64_t)131071;
    lCluster.Field4.Eem = (Eem >> 18) & 0x0003FFFF;
    lCluster.Field4.xEem = (Eem >> 36);

    if (lDecoded.TriggerLayer >= 15 and lDecoded.TriggerLayer <= 18) {
      lCluster.Field4.Ehearly = lCluster.Field0.E;
      lCluster.Field4.xEhearly = lCluster.Field0.xE;
    }

    if (lDecoded.TriggerLayer >= 4 and lDecoded.TriggerLayer <= 8) {
      lCluster.Field4.Eemcore = lCluster.Field4.Eem;
      lCluster.Field4.xEemcore = lCluster.Field4.xEem;
    }

    uint64_t Layer = (uint64_t)(0x1) << lDecoded.TriggerLayer;
    lCluster.Field1.LayerBits = (Layer >> 3) & 0x7;
    lCluster.Field2.LayerBits = (Layer >> 6) & 0x7;
    lCluster.Field3.LayerBits = (Layer >> 9) & 0x3FFF;
    lCluster.Field4.LayerBits = (Layer >> 23) & 0x3FFF;
  }

  // ----------------------------------------------------------------
  std::array< TcDecoder , 4096 > InitDummyClusterDecoderROM()
  {

    std::array< TcDecoder , 4096 > lArray;

    const uint32_t cDepths[51] = { 0 , // No zero layer? True.
          0 , 30 , 59 , 89 , 118 , 148 , // CE-E (early)
          178 , 208 , 237 , 267 , 297 , 327 , 356 , 386 , // CE-E (core)
          415 , 445 , 475 , 505 , 534 , 564 , 594 , 624 , 653 , 683 , 712 , 742 , 772 , 802 , // CE-E (back)
          911 , 1020 , 1129 , 1238 , // CE-H (early)
          1347 , 1456 , 1565 , 1674 , 1783 , 1892 , 2001 , 2110 , 2281 , 2452 , 2623 , 2794 , 2965 , 3136 , 3307 , 3478 , 3649 , 3820 }; // CE-H (back)

    const uint32_t cTriggerLayers[51] = { 0 , // No zero layer
          1 , 0 , 2 , 0 , 3 , 0 , // CE-E (early)
          4 , 0 , 5 , 0 , 6 , 0 , 7 , 0 , 8 , 0 , // CE-E (core)
          9 , 0 , 10 , 0 , 11 , 0 , 12 , 0 , 13 , 0 , 14 , 0 , // CE-E (back)
          15 , 16 , 17 , 18 , // CE-H (early)
          19 , 20 , 21 , 22 , 23 , 24 , 25 , 26 , 27 , 28 , 29 , 30 , 31 , 32 , 33 , 34 , 35 , 36 }; // CE-H (back)

    const uint32_t cLayerWeights_E_EM[51] = { 0 , // No zero layer
          252969 , 0 , 254280 , 0 , 255590 , 0 , // CE-E (early)
          256901 , 0 , 258212 , 0 , 259523 , 0 , 260833 , 0 , 262144 , 0 , // CE-E (core)
          263455 , 0 , 264765 , 0 , 266076 , 0 , 267387 , 0 , 268698 , 0 , 270008 , 0 , // CE-E (back)
          0 , 0 , 0 , 0 , // CE-H (early)
          0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 }; // CE-H (back)

    for( int i(1); i!=4096; ++i )
    {
      int Plane                 = ( i % 49 ) + 1;
      lArray.at(i).RoverZ       = ( 4091 * i ) % 4096;
      lArray.at(i).Phi          = ( 4093 * i ) % 4096;
      lArray.at(i).TriggerLayer = cTriggerLayers[Plane];
      lArray.at(i).LayerDepth   = cDepths[Plane];
      lArray.at(i).Cell         = ( 4 * ( i % 54 ) ) + ( i % 4 ) ;
      lArray.at(i).LayerWeight  = cLayerWeights_E_EM[Plane];
    }

    return lArray;
  }
  // ----------------------------------------------------------------



  void Cluster_Decoders(const std::array<std::array<LinkTriggerCell, 154>, 162>& cellsIn,
                        const std::array<TcDecoder, 4096>& clusterDecoderROM,
                        std::array<std::array<Cluster, 154>, 162>& clusterOut) {

    for ( int j(0); j!= 154; ++j ) { // decoder
      for ( int i(0); i!= 162; ++i ) { // clock

        const LinkTriggerCell& lCell = cellsIn.at(i).at(j);
        Cluster& lCluster = clusterOut.at(i).at(j);
        lCluster = Cluster();
        lCluster.Last = lCell.Last;

        // Convert clock-cycle and channel to LUT block-index
        uint16_t BlockIndex = UnpackingLut[j][i];

        // Empty virtual lanes need no ROM lookup.  The complete firmware LUT
        // addresses more blocks than the temporary CMSSW virtual ROM, so make
        // the bounds requirement explicit before forming its index.
        const unsigned decoderAddress = (64 * BlockIndex) + lCell.TcId;
        if (BlockIndex == 0 || !lCell.DataValid || decoderAddress >= clusterDecoderROM.size())
          continue;

        // Look-up the decoder value
        const TcDecoder& lDecoded = clusterDecoderROM[decoderAddress];
        if (!lDecoded.DataValid)
          continue;
        Cluster_Decode(lCell, lDecoded, lCluster);
      }
    }
  }

  void Cluster_Decoders(const std::array<std::array<LinkTriggerCell, 154>, 162>& cellsIn,
                        std::array<std::array<Cluster, 154>, 162>& clusterOut) {
    static const std::array<TcDecoder, 4096> clusterDecoderROM = InitDummyClusterDecoderROM();
    Cluster_Decoders(cellsIn, clusterDecoderROM, clusterOut);
  }


}
