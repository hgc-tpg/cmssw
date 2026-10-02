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

void Cluster_PackLinks( const std::array< ClusterProperty         , 162 >& ClusterIn ,
                              std::array< std::array< lword , 4 > , 162 >& ClusterLinksOut )
{


  for ( int i(0); i!= 162; ++i ) { // clock

    auto& lIn = ClusterIn.at(i);

    {
      auto& lOut = ClusterLinksOut.at(i).at(0);
      lOut.data = 0;
      lOut.data |= ( (uint64_t) lIn.ET << (0 % 64) ); //  13 DOWNTO
      lOut.data |= ( (uint64_t) lIn.e_or_gamma_ET << (14 % 64) ); //  27 DOWNTO
      lOut.data |= ( (uint64_t) lIn.GCT_e_or_gamma_Select_0 << (28 % 64) ); //
      lOut.data |= ( (uint64_t) lIn.GCT_e_or_gamma_Select_1 << (29 % 64) ); //
      lOut.data |= ( (uint64_t) lIn.GCT_e_or_gamma_Select_2 << (30 % 64) ); //
      lOut.data |= ( (uint64_t) lIn.GCT_e_or_gamma_Select_3 << (31 % 64) ); //
      lOut.data |= ( (uint64_t) lIn.Fraction_in_CE_E << (32 % 64) ); //  39 DOWNTO
      lOut.data |= ( (uint64_t) lIn.Fraction_in_core_CE_E << (40 % 64) ); //  47 DOWNTO
      lOut.data |= ( (uint64_t) lIn.Fraction_in_front_CE_H << (48 % 64) ); //  55 DOWNTO
      lOut.data |= ( (uint64_t) lIn.First_Layer << (56 % 64) ); //  61 DOWNTO
      // lOut.data |= ( (uint64_t) lIn.Spare_2 << (62 % 64) ); //  63 DOWNTO

      lOut.strobe = true;
      lOut.valid  = lIn.DataValid;
      lOut.start  = (i==0);
      lOut.last   = (i==161);
    }

    {
      auto& lOut = ClusterLinksOut.at(i).at(1);
      lOut.data = 0;
      lOut.data |= ( (uint64_t) lIn.ET_Weighted_Eta << (64 % 64) ); //  73 DOWNTO
      lOut.data |= ( (uint64_t) lIn.ET_Weighted_Phi << (74 % 64) ); //  82 DOWNTO
      lOut.data |= ( (uint64_t) lIn.ET_Weighted_Z << (83 % 64) ); //  94 DOWNTO
      // lOut.data |= ( (uint64_t) lIn.Spare_3 << (95 % 64) ); //  95 DOWNTO
      lOut.data |= ( (uint64_t) lIn.Number_of_Cells << (96 % 64) ); //  105 DOWNTO
      lOut.data |= ( (uint64_t) lIn.Saturated_Trigger_Cell << (106 % 64) ); //
      lOut.data |= ( (uint64_t) lIn.Quality_of_Fraction_in_CE_E << (107 % 64) ); //
      lOut.data |= ( (uint64_t) lIn.Quality_of_Fraction_in_core_CE_E << (108 % 64) ); //
      lOut.data |= ( (uint64_t) lIn.Quality_of_Fraction_in_front_CE_H << (109 % 64) ); //
      lOut.data |= ( (uint64_t) lIn.Quality_of_Sigmas_and_Means << (110 % 64) ); //
      lOut.data |= ( (uint64_t) lIn.Saturated_Phi << (111 % 64) ); //
      lOut.data |= ( (uint64_t) lIn.Nominal_Phi << (112 % 64) ); //
      // lOut.data |= ( (uint64_t) lIn.Spare_4 << (113 % 64) ); //  127 DOWNTO

      lOut.strobe = true;
      lOut.valid  = lIn.DataValid;
      lOut.start  = (i==0);
      lOut.last   = (i==161);
    }

    {
      auto& lOut = ClusterLinksOut.at(i).at(2);
      lOut.data = 0;
      lOut.data |= ( (uint64_t) lIn.Sigma_E << (128 % 64) ); //  134 DOWNTO
      lOut.data |= ( (uint64_t) lIn.Last_Layer << (135 % 64) ); //  140 DOWNTO
      lOut.data |= ( (uint64_t) lIn.Shower_Length << (141 % 64) ); //  146 DOWNTO
      lOut.data |= ( (uint64_t) lIn.CoreShowerLen << (147 % 64) ); //  152 DOWNTO
      // lOut.data |= ( (uint64_t) lIn.Spare_5 << (153 % 64) ); //  159 DOWNTO
      lOut.data |= ( (uint64_t) lIn.Sigma_ZZ << (160 % 64) ); //  166 DOWNTO
      lOut.data |= ( (uint64_t) lIn.Sigma_PhiPhi << (167 % 64) ); //  173 DOWNTO
      lOut.data |= ( (uint64_t) lIn.Sigma_EtaEta << (174 % 64) ); //  180 DOWNTO
      lOut.data |= ( (uint64_t) lIn.Sigma_RoZRoZ << (181 % 64) ); //  187 DOWNTO
      lOut.data |= ( (uint64_t) lIn.LayerBits_Hi << (188 % 64) ); //  189 DOWNTO
      // lOut.data |= ( (uint64_t) lIn.Spare_6 << (190 % 64) ); //  191 DOWNTO

      lOut.strobe = true;
      lOut.valid  = lIn.DataValid;
      lOut.start  = (i==0);
      lOut.last   = (i==161);
    }

    {
      auto& lOut = ClusterLinksOut.at(i).at(3);
      lOut.data = 0;
      lOut.data |= ( (uint64_t) lIn.LayerBits_Lo << (192 % 64) ); //  223 DOWNTO
      // lOut.data |= ( (uint64_t) lIn.Spare_8 << (224 % 64) ); //  255 DOWNTO

      lOut.strobe = true;
      lOut.valid  = lIn.DataValid;
      lOut.start  = (i==0);
      lOut.last   = (i==161);
    }

  }
}

}
