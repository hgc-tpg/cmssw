// -----------------------------------------------------------------------------------------------------------------------
// Andrew W. Rose, 2026
// Imperial College London Particles community
// -----------------------------------------------------------------------------------------------------------------------

#pragma once

#include <cstdint>
#include <array>
#include <format>
#include "L1Trigger/L1THGCal/interface/Magic.hpp"


// ----------------------------------------------------------------------------------------------------------------------
struct lword : public magic< lword >
{
  lword() = default;
  MAGIC_CONSTRUCTOR( lword );
  MAGIC( data , valid , start , start_of_orbit , strobe , last );

  uint64_t data;
  bool valid;
  bool start;
  bool start_of_orbit;
  bool strobe;
  bool last;
};
// ----------------------------------------------------------------------------------------------------------------------

// ----------------------------------------------------------------------------------------------------------------------
struct LinkTriggerCell : public magic< LinkTriggerCell >
{
  LinkTriggerCell() = default;
  MAGIC_CONSTRUCTOR( LinkTriggerCell );
  MAGIC( Mantissa , Exponent , TcId , Last , DataValid );

  uint8_t Mantissa;
  uint8_t Exponent;
  uint8_t TcId;
  bool Last;
  bool DataValid;
};
// ----------------------------------------------------------------------------------------------------------------------

// ----------------------------------------------------------------------------------------------------------------------
struct LinkTower : public magic< LinkTower >
{
  LinkTower() = default;
  MAGIC_CONSTRUCTOR( LinkTower );
  MAGIC( EcalExponent , EcalMantissa , HcalExponent , HcalMantissa , Last , DataValid );

  uint8_t EcalExponent;
  uint8_t EcalMantissa;
  uint8_t HcalExponent;
  uint8_t HcalMantissa;
  bool Last;
  bool DataValid;
};
// ----------------------------------------------------------------------------------------------------------------------

// ----------------------------------------------------------------------------------------------------------------------
struct LinkFlags : public magic< LinkFlags >
{
  LinkFlags() = default;
  MAGIC_CONSTRUCTOR( LinkFlags );
  MAGIC( Data , DataValid , Last );

  uint64_t Data;
  bool Last;
  bool DataValid;
};
// ----------------------------------------------------------------------------------------------------------------------




// ----------------------------------------------------------------------------------------------------------------------
struct Cluster : public magic< Cluster >
{

  Cluster() : Cell( 0xFF ) , ColumnSet(0) , DataValid(false) , Last(false) ,
              DebugCol( -999 ) , DebugCol2( -999 ) , DebugRow( -999 ) , DebugRow2( -999 ) , DebugRow3( -999 ) ,
              Field0() , Field1() , Field2() , Field3() , Field4()
  {}

  MAGIC_CONSTRUCTOR( Cluster );
  MAGIC( Cell , ColumnSet , DataValid , Last , DebugCol , DebugCol2 , DebugRow , DebugRow2 , DebugRow3 , Field0 , Field1 , Field2 , Field3 , Field4 );


  struct tField0 : public magic< tField0 >
  {
    tField0() = default;
    MAGIC_CONSTRUCTOR( tField0 );
    MAGIC( W2 , xW2 , W , xW , E , xE , EventFlag , ShapeQ );

    uint32_t W2;
    bool xW2;
    uint16_t W;
    bool xW;
    uint32_t E;
    bool xE;
    uint8_t EventFlag;
    uint8_t ShapeQ;

  };

  struct tField1 : public magic< tField1 >
  {
    tField1() = default;
    MAGIC_CONSTRUCTOR( tField1 );
    MAGIC( WZ2 , xWZ2 , WZ , xWZ , N_TC , xN_TC , LayerBits );

    uint32_t WZ2;
    bool xWZ2;
    uint32_t WZ;
    bool xWZ;
    uint16_t N_TC;
    bool xN_TC;
    uint16_t LayerBits;

  };

  struct tField2 : public magic< tField2 >
  {
    tField2() = default;
    MAGIC_CONSTRUCTOR( tField2 );
    MAGIC( Wroz2 , xWroz2 , Wroz , xWroz , N_TC_W , xN_TC_W , LayerBits );

    uint32_t Wroz2;
    bool xWroz2;
    uint32_t Wroz;
    bool xWroz;
    uint16_t N_TC_W;
    bool xN_TC_W;
    uint16_t LayerBits;

  };

  struct tField3 : public magic< tField3 >
  {
    tField3() = default;
    MAGIC_CONSTRUCTOR( tField3 );
    MAGIC( Wphi2 , xWphi2 , Wphi , xWphi , LayerBits );

    uint32_t Wphi2;
    bool xWphi2;
    uint32_t Wphi;
    bool xWphi;
    uint16_t LayerBits;

  };

  struct tField4 : public magic< tField4 >
  {
    tField4() = default;
    MAGIC_CONSTRUCTOR( tField4 );
    MAGIC( Ehearly , xEhearly , Eem , xEem , Eemcore , xEemcore , LayerBits );

    uint32_t Ehearly;
    bool xEhearly;
    uint32_t Eem;
    bool xEem;
    uint32_t Eemcore;
    bool xEemcore;
    uint16_t LayerBits;

  };

  uint8_t Cell;
  int8_t ColumnSet;
  bool DataValid;
  bool Last;
  int32_t DebugCol , DebugCol2;
  int32_t DebugRow , DebugRow2 , DebugRow3;
  tField0 Field0;
  tField1 Field1;
  tField2 Field2;
  tField3 Field3;
  tField4 Field4;

};

template< typename T >
void AddField( const uint32_t& aSize , const uint32_t& aMask , const T& aA , const bool& xA , const T& aB , const bool& xB , T& aC , bool& xC )
{
  auto lSum = (uint64_t)aA + (uint64_t)aB;
  aC = (uint32_t)lSum & aMask;
  xC = xA or xB or (bool)(lSum >> aSize);
}

Cluster::tField0& operator+= ( Cluster::tField0& left , const Cluster::tField0& right )
{
  AddField( 32 , 0xFFFFFFFF , left.W2 , left.xW2 , right.W2 , right.xW2 , left.W2 , left.xW2 );
  AddField( 16 , 0x0000FFFF , left.W  , left.xW  , right.W  , right.xW  , left.W  , left.xW  );
  AddField( 18 , 0x0003FFFF , left.E  , left.xE  , right.E  , right.xE  , left.E  , left.xE  );
  left.ShapeQ = left.ShapeQ | right.ShapeQ;
  return left;
}

Cluster::tField1& operator+= ( Cluster::tField1& left , const Cluster::tField1& right )
{
  AddField( 32 , 0xFFFFFFFF , left.WZ2  , left.xWZ2  , right.WZ2  , right.xWZ2  , left.WZ2  , left.xWZ2  );
  AddField( 24 , 0x00FFFFFF , left.WZ   , left.xWZ   , right.WZ   , right.xWZ   , left.WZ   , left.xWZ   );
  AddField( 10 , 0x000003FF , left.N_TC , left.xN_TC , right.N_TC , right.xN_TC , left.N_TC , left.xN_TC );
  left.LayerBits = left.LayerBits | right.LayerBits;
  return left;
}

Cluster::tField2& operator+= ( Cluster::tField2& left , const Cluster::tField2& right )
{
  AddField( 32 , 0xFFFFFFFF , left.Wroz2  , left.xWroz2  , right.Wroz2  , right.xWroz2  , left.Wroz2  , left.xWroz2  );
  AddField( 24 , 0x00FFFFFF , left.Wroz   , left.xWroz   , right.Wroz   , right.xWroz   , left.Wroz   , left.xWroz   );
  AddField( 10 , 0x000003FF , left.N_TC_W , left.xN_TC_W , right.N_TC_W , right.xN_TC_W , left.N_TC_W , left.xN_TC_W );
  left.LayerBits = left.LayerBits | right.LayerBits;
  return left;
}

Cluster::tField3& operator+= ( Cluster::tField3& left , const Cluster::tField3& right )
{
  AddField( 32 , 0xFFFFFFFF , left.Wphi2 , left.xWphi2 , right.Wphi2 , right.xWphi2 , left.Wphi2 , left.xWphi2 );
  AddField( 24 , 0x00FFFFFF , left.Wphi  , left.xWphi  , right.Wphi  , right.xWphi  , left.Wphi  , left.xWphi  );
  left.LayerBits = left.LayerBits | right.LayerBits;
  return left;
}

Cluster::tField4& operator+= ( Cluster::tField4& left , const Cluster::tField4& right )
{
  AddField( 18 , 0x0003FFFF , left.Ehearly , left.xEhearly , right.Ehearly , right.xEhearly , left.Ehearly , left.xEhearly );
  AddField( 18 , 0x0003FFFF , left.Eem     , left.xEem     , right.Eem     , right.xEem     , left.Eem     , left.xEem     );
  AddField( 18 , 0x0003FFFF , left.Eemcore , left.xEemcore , right.Eemcore , right.xEemcore , left.Eemcore , left.xEemcore );
  left.LayerBits = left.LayerBits | right.LayerBits;
  return left;
}

Cluster& operator+= ( Cluster& left , const Cluster& right )
{
  left.Field0 += right.Field0;
  left.Field1 += right.Field1;
  left.Field2 += right.Field2;
  left.Field3 += right.Field3;
  left.Field4 += right.Field4;
  return left;
}
// ----------------------------------------------------------------------------------------------------------------------

// ----------------------------------------------------------------------------------------------------------------------
struct TcDecoder : public magic< TcDecoder >
{
  TcDecoder() = default;
  MAGIC_CONSTRUCTOR( TcDecoder );
  MAGIC( RoverZ , Phi , TriggerLayer , LayerDepth , Cell , LayerWeight , Index , E , Eem , IsEemcore , IsEhearly , DataValid , Last );

  uint16_t RoverZ;
  uint16_t Phi;
  uint8_t TriggerLayer;
  uint16_t LayerDepth;
  uint8_t Cell;
  uint32_t LayerWeight;
  //UNSIGNED( 21 DOWNTO 19 ) Spare;
// Utility fields
  uint16_t Index;
  uint64_t E , Eem;
  bool IsEemcore , IsEhearly;
  bool DataValid;
  bool Last;

};
// ----------------------------------------------------------------------------------------------------------------------

// ----------------------------------------------------------------------------------------------------------------------
struct ClusterProperty : public magic< ClusterProperty >
{

  ClusterProperty() = default;
  MAGIC_CONSTRUCTOR( ClusterProperty );
  MAGIC( ET , e_or_gamma_ET , GCT_e_or_gamma_Select_0 , GCT_e_or_gamma_Select_1 , GCT_e_or_gamma_Select_2 , GCT_e_or_gamma_Select_3 ,
          Fraction_in_CE_E , Fraction_in_core_CE_E , Fraction_in_front_CE_H , First_Layer ,
          ET_Weighted_Eta , ET_Weighted_Phi , ET_Weighted_Z ,
          Number_of_Cells , Saturated_Trigger_Cell , Quality_of_Fraction_in_CE_E , Quality_of_Fraction_in_core_CE_E , Quality_of_Fraction_in_front_CE_H , Quality_of_Sigmas_and_Means , Saturated_Phi , Nominal_Phi ,
          Sigma_E , Last_Layer , Shower_Length , CoreShowerLen ,
          Sigma_ZZ , Sigma_PhiPhi , Sigma_EtaEta , Sigma_RoZRoZ , LayerBits_Hi , LayerBits_Lo ,
          DataValid , Last );

// First 32b word
  uint16_t ET;
  uint16_t e_or_gamma_ET;
  bool GCT_e_or_gamma_Select_0;
  bool GCT_e_or_gamma_Select_1;
  bool GCT_e_or_gamma_Select_2;
  bool GCT_e_or_gamma_Select_3;

// Second 32b word
  uint16_t Fraction_in_CE_E;
  uint16_t Fraction_in_core_CE_E;
  uint16_t Fraction_in_front_CE_H;
  uint8_t First_Layer;
  // uint8_t Spare_2;

// Third 32b word
  uint16_t ET_Weighted_Eta;
  int16_t ET_Weighted_Phi;
  uint16_t ET_Weighted_Z;
  // uint8_t Spare_3;

// Fourth 32b word
  uint16_t Number_of_Cells;
  bool Saturated_Trigger_Cell;
  bool Quality_of_Fraction_in_CE_E;
  bool Quality_of_Fraction_in_core_CE_E;
  bool Quality_of_Fraction_in_front_CE_H;
  bool Quality_of_Sigmas_and_Means;
  bool Saturated_Phi;
  bool Nominal_Phi;
  // uint16_t Spare_4;

// Fifth 32b word
  uint16_t Sigma_E;
  uint8_t Last_Layer;
  uint8_t Shower_Length;
  uint8_t CoreShowerLen;
  // uint8_t Spare_5;

// Sixth 32b word
  uint16_t Sigma_ZZ;
  uint16_t Sigma_PhiPhi;
  uint16_t Sigma_EtaEta;
  uint16_t Sigma_RoZRoZ;
  uint8_t LayerBits_Hi;
  // uint8_t Spare_6;

// Seventh 32b word
  uint32_t LayerBits_Lo;

// Eighth 32b word
  // uint32_t Spare_8;

//    -- (Other stuff we have available)
//    E_EM_core                             : UNSIGNED( 21 DOWNTO 0 );
//    E_H_early                             : UNSIGNED( 21 DOWNTO 0 );
//    ET_Weighted_RoZ                       : uint16_t;

  bool DataValid;
  bool Last;

};

// ----------------------------------------------------------------------------------------------------------------------
struct DebugBoolean : public magic< DebugBoolean >
{
  DebugBoolean() = default;
  MAGIC_CONSTRUCTOR( DebugBoolean );
  MAGIC( Value , A , B );

  bool Value;
  Cluster A , B;
};
// ----------------------------------------------------------------------------------------------------------------------

// ----------------------------------------------------------------------------------------------------------------------
struct DebugBoolean2 : public magic< DebugBoolean2 >
{
  DebugBoolean2() = default;
  MAGIC_CONSTRUCTOR( DebugBoolean2 );
  MAGIC( Value , Ce , UpLt , Up , UpRt , DnLt , Dn , DnRt );

  bool Value;
  Cluster Ce , UpLt , Up , UpRt , DnLt , Dn , DnRt;
};
// ----------------------------------------------------------------------------------------------------------------------




// ----------------------------------------------------------------------------------------------------------------------
struct Tower : public magic< Tower >
{
  Tower() = default;
  MAGIC_CONSTRUCTOR( Tower );
  MAGIC( CEE , CEH , Index , Index2 , DataValid , Last , DEBUG_ETA , DEBUG_PHI );

  uint32_t CEE;
  uint32_t CEH;
  uint16_t Index;
  uint16_t Index2;
  bool DataValid;
  bool Last;
// Debug info
  uint16_t DEBUG_ETA;
  uint16_t DEBUG_PHI;

};

struct FormattedTower : public magic< FormattedTower >
{
  FormattedTower() = default;
  MAGIC_CONSTRUCTOR( FormattedTower );
  MAGIC( Value , Ratio , Saturated , DataValid , Last , DEBUG_ETA , DEBUG_PHI );

  uint32_t Value;
  uint32_t Ratio;
  bool     Saturated;
  bool DataValid;
  bool Last;
// Debug info
  uint16_t DEBUG_ETA;
  uint16_t DEBUG_PHI;
};


struct PackedTower : public magic< PackedTower >
{
  PackedTower() = default;
  MAGIC_CONSTRUCTOR( PackedTower );
  MAGIC( Value , DEBUG_ETA0 , DEBUG_PHI0 , DEBUG_ETA1 , DEBUG_PHI1 , DEBUG_ETA2 , DEBUG_PHI2 , DEBUG_ETA3 , DEBUG_PHI3 );

  lword Value;

  uint16_t DEBUG_ETA0;
  uint16_t DEBUG_PHI0;
  uint16_t DEBUG_ETA1;
  uint16_t DEBUG_PHI1;
  uint16_t DEBUG_ETA2;
  uint16_t DEBUG_PHI2;
  uint16_t DEBUG_ETA3;
  uint16_t DEBUG_PHI3;
};
// ----------------------------------------------------------------------------------------------------------------------




// ----------------------------------------------------------------------------------------------------------------------
void PhysicalToFirmware( const std::array< std::array< Cluster , 44 > , 41 >& aPhysical , std::array< std::array< Cluster , 11 > , 162 >& aFirmware )
{
  // Convert coordinate-oriented data to firmware-oriented data
  for ( int j(0); j!=44 ; ++j ) { // col
    for ( int i(0); i!=41 ; ++i ) { // row
      const Cluster& lRef = aPhysical.at(i).at(j);
      if( lRef.Cell > 161 ) continue;
      aFirmware.at( lRef.Cell ).at( j/4 ) = lRef;
    }
  }
}

void PhysicalToFirmware( const std::array< std::array< DebugBoolean , 44 > , 41 >& aPhysical , std::array< std::array< DebugBoolean , 13 > , 162 >& aFirmware )
{
  // Convert coordinate-oriented data to firmware-oriented data
  for ( int j(0); j!=44 ; ++j ) { // col
    for ( int i(0); i!=41 ; ++i ) { // row
      const DebugBoolean& lRef = aPhysical.at(i).at(j);
      if( lRef.A.Cell > 161 ) continue;
      aFirmware.at( lRef.A.Cell ).at( (j/4)+1 ) = lRef;
    }
  }
}

void PhysicalToFirmware( const std::array< std::array< DebugBoolean2 , 44 > , 41 >& aPhysical , std::array< std::array< DebugBoolean2 , 11 > , 162 >& aFirmware )
{
  // Convert coordinate-oriented data to firmware-oriented data
  for ( int j(0); j!=44 ; ++j ) { // col
    for ( int i(0); i!=41 ; ++i ) { // row
      const DebugBoolean2& lRef = aPhysical.at(i).at(j);
      if( lRef.Ce.Cell > 161 ) continue;
      aFirmware.at( lRef.Ce.Cell ).at( j/4 ) = lRef;
    }
  }
}

void PhysicalToFirmware( const std::array< std::array< FormattedTower , 24 > , 20 >& aPhysical , std::array< std::array< FormattedTower , 16 > , 162 >& aFirmware )
{
  // Convert coordinate-oriented data to firmware-oriented data
  for ( int j(0); j!=24 ; ++j ) { // phi
    for ( int i(0); i!=20 ; ++i ) { // eta
      const FormattedTower& lRef = aPhysical.at(i).at(j);
      // if( lRef.DataValid ) {
        // if( j != lRef.DEBUG_PHI ) throw std::runtime_error( std::format( "Phi: {} != {}" , j , lRef.DEBUG_PHI ) );
        // if( i != lRef.DEBUG_ETA ) throw std::runtime_error( std::format( "Eta: {} != {}" , i , lRef.DEBUG_ETA ) );
      // }
      uint32_t clk     = (i/4) + 5*(j%6);
      uint32_t channel = (i%4) + 4*(j/6);

      aFirmware.at( clk ).at( channel ) = lRef;
    }
  }
}
// ----------------------------------------------------------------------------------------------------------------------
