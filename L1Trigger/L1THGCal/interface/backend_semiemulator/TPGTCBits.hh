#ifndef TPGTCBits_hh
#define TPGTCBits_hh

#include <iostream>
#include <iomanip>
#include <cstdint>
#include <cstring>
#include <cassert>

#include "ap_int.h"

// Updated following the code from Paul, dated 06/11/2024
// The format below from Andy

//////////////// OLD Format /////////////////////////
/*

 TYPE tData IS RECORD
    Energy    : UNSIGNED( 17 DOWNTO 0 );
    XoverZ    :   SIGNED( 11 DOWNTO 0 );
    YoverZ    :   SIGNED( 11 DOWNTO 0 );
    Z         : UNSIGNED( 14 DOWNTO 0 );
    Plane     : UNSIGNED(  5 DOWNTO 0 );
    U         :   SIGNED(  6 DOWNTO 0 );
    U_bar     : UNSIGNED(  6 DOWNTO 0 );
    DataValid : BOOLEAN;
    Last      : BOOLEAN;
    BypassFlag : BOOLEAN;
    EndOfCluster : BOOLEAN;

 */
/////////////////////////////////////////////////////

////////////// TCBit format //////////////////////////
// https://gitlab.cern.ch/hgcal-tpg/Stage2/-/blob/master/firmware/hdl/types/PkgTriggerCell.vhd?ref_type=heads
/*
    // Bit widths =====================
    RoverZ    : UNSIGNED( 11 DOWNTO 0 );
    Phi       : UNSIGNED( 11 DOWNTO 0 );
    -- Z         : UNSIGNED( 14 DOWNTO 0 );
    Plane     : UNSIGNED(  5 DOWNTO 0 );
    Latent    : STD_LOGIC_VECTOR( 10 DOWNTO 0 );
    U         :   SIGNED(  5 DOWNTO 0 );
    U_bar     : UNSIGNED(  5 DOWNTO 0 );
    Energy    : UNSIGNED( 17 DOWNTO 0 );
    DataValid : BOOLEAN;
    Last      : BOOLEAN;
    BypassFlag : BOOLEAN;
    EndOfCluster : BOOLEAN;
    //=================================

    // Bit packing ====================
    lRet( 11 DOWNTO 0 ) := STD_LOGIC_VECTOR( aData.RoverZ );
    lRet( 23 DOWNTO 12 ) := STD_LOGIC_VECTOR( aData.Phi );
    lRet( 29 DOWNTO 24 ) := STD_LOGIC_VECTOR( aData.Plane );
    lRet( 40 DOWNTO 30 ) := aData.Latent;
    lRet( 46 DOWNTO 41 ) := STD_LOGIC_VECTOR( aData.U );
    lRet( 52 DOWNTO 47 ) := STD_LOGIC_VECTOR( aData.U_bar );
    lRet( 70 DOWNTO 53 ) := STD_LOGIC_VECTOR( aData.Energy );
    lRet( 71 )           := TO_STD_LOGIC( aData.BypassFlag );
    //=================================
 */
///////////////////////////////////////////////////

class TPGTCBits {
public:
  TPGTCBits() : cmsswIndex_(-1) { setZero(); }
  void setZero() { _data=0; }

  ap_uint<72>  getData() const { return _data; }
  void setData(ap_uint<72> d) { _data=d; }

  ap_uint<12>	getROverZ() const { return _data.range(11,0); }
  ap_uint<12>	getPhi() const { return _data.range(23,12); }
  ap_uint<6>	getLayer() const { return _data.range(29,24); }
  //ap_uint<11> getLatent() const { return _data.range(40,30); }
  ap_int<6>	getU() const { return _data.range(46,41); }
  ap_uint<6>	getUbar() const {  return _data.range(52,47); }
  ap_uint<18>	getEnergy() const { return _data.range(70,53); }
  bool          getBypassFlag() { return (_data.range(71,71) == 1)?true:false; }

  void setROverZ(uint16_t roz) { _data(11,0) = (roz>0xfff)?0xfff:roz;}
  void setPhi(uint16_t phi) { _data(23,12) = (phi>0xfff)?0xfff:phi;}
  void setLayer(uint16_t l) { _data(29,24) = (l>0x3f)?0x3f:l;}
  void setU(int16_t u) { _data(46,41) = (u>0x1f or u<-0x20)?-0x20:u;}
  void setUbar(uint16_t u) { _data(52,47) = (u>0x3f)?0x3f:u;}
  void setEnergy(uint32_t e) { _data(70,53) = (e>0x3ffff)?0x3ffff:e;}
  void setBypassFlag(bool flag) { _data(71,71) = (flag)?1:0;}

  int getCMSSWIndex() const { return cmsswIndex_; }
  void setCMSSWIndex(int idx) { cmsswIndex_ = idx; }

  void print() {
    std::cout << "TPGTCBits(" << this << ")::print() = " ;
	      // << std::hex << std::setfill('0')
	      // << std::setw(18) << _data
	      // << std::dec << std::setfill(' ')
	      // << std::endl;

    std::cout << " Energy " << getEnergy() << " *2^-9 GeV"
	      << ", (r/z, phi) = (" << getROverZ() << " *1.4e-4"
	      << ", " << getPhi() << " *0.8 mrad)"
	      << ", layer = " << getLayer()
	      << ", (u, u-bar) = (" << getU()
	      << ", " << getUbar() << ")" << std::endl;
  }

private:
  ap_uint<72> _data;

  int cmsswIndex_;
};

#endif
