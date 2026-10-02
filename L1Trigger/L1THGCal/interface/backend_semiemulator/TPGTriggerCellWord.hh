#ifndef TPGTriggerCellWord_hh
#define TPGTriggerCellWord_hh

#include <iostream>
#include <iomanip>
#include <cstdint>
#include <cstring>
#include <cassert>

// Code from Paul, dated 06/11/2024
// The format below from Andy

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


class TPGTriggerCellWord {
public:
  TPGTriggerCellWord() {
    setZero();
  }

  void setZero() {
    _data=0;
  }

  uint64_t getData() const {
    return _data;
  }

  void setData(uint64_t d) {
    _data=d;
  }

  uint32_t getEnergy() const {
    return _data&0x3ffff;
  }

  int16_t getXOverZ() const {
    uint16_t u((_data>>18)&0xfff);
    return u<0x800?u:int16_t(u)-0x1000;
  }

  int16_t getYOverZ() const {
    uint16_t u((_data>>30)&0xfff);
    return u<0x800?u:int16_t(u)-0x1000;
  }
  /*
  uint16_t getZ() const {
    return (_data>>42)&0x7fff;
  }
  */
  uint16_t getLayer() const {
    return (_data>>42)&0x3f;
  }

  int16_t getU() const {
    uint16_t u((_data>>48)&0x7f);
    return u<0x40?u:int16_t(u)-0x80;
  }

  uint16_t getUbar() const {
    return (_data>>55)&0x7f;
  }

  void setEnergy(uint32_t e) {
    if(e>0x3ffff) e=0x3ffff;
    _data&=0xfffffffffffc0000;
    _data|=e;
  }

  void setXOverZ(int16_t x) {
    uint64_t xd;
    if(x> 0x7ff) xd= 0x7ff;
    else if(x>=0) xd=x;
    else if(x<-0x800) xd=0x800;
    else xd=0x1000+x;
    _data&=0xffffffffc003ffff;
    _data|=(xd&0xfff)<<18;
  }

  void setYOverZ(int16_t y) {
    uint64_t yd;
    if(y> 0x7ff) yd= 0x7ff;
    else if(y>=0) yd=y;
    else if(y<-0x800) yd=0x800;
    else yd=0x1000+y;
    _data&=0xfffffc003fffffff;
    _data|=(yd&0xfff)<<30;
  }

  void setLayer(uint16_t l) {
    if(l>47) l=63;
    _data&=0xffff03ffffffffff;
    _data|=uint64_t(l&0x3f)<<42;
  }

  void setU(int16_t u) {
    uint64_t ud;
    if(u> 0x3f) ud= 0x3f;
    else if(u>=0) ud=u;
    else if(u<-0x4f) ud=0x4f;
    else ud=0x80+u;
    _data&=0xff80ffffffffffff;
    _data|=(ud&0xfff)<<48;
  }

  void setUbar(uint16_t u) {
    if(u>0x7f) u=0x7f;
    _data&=0xc07fffffffffffff;
    _data|=uint64_t(u&0x7f)<<55;
  }

  void print() {
    std::cout << "TPGTriggerCellWord(" << this << ")::print() = 0x"
	      << std::hex << std::setfill('0')
	      << std::setw(16) << _data
	      << std::dec << std::setfill(' ')
	      << std::endl;

    std::cout << " Energy " << getEnergy() << " 2^-8 GeV" << std::endl
	      << " x/z, y/z = " << getXOverZ()
	      << ", " << getYOverZ()
	      << ", layer = " << getLayer() << std::endl
	      << " u, u-bar = " << getU()
	      << ", " << getUbar() << std::endl;
  }

private:
  uint64_t _data;
};

#endif
