#ifndef TPGClusterData_hh
#define TPGClusterData_hh

#include <iostream>
#include <iomanip>
#include <cstdint>
#include <cstring>
#include <cassert>

// Code from Paul, dated 10/11/2024

class TPGClusterData {
public:
  TPGClusterData() {
    setZero();
  }

  void setZero() {
    for(unsigned w(0);w<4;w++) _data[w]=0;
  }

  const uint64_t* getData() const {
    return _data;
  }

  void setData(const uint64_t *d) {
    for(unsigned w(0);w<4;w++) _data[w]=d[w];
  }

  uint16_t getEnergy() const {
    return _data[0]&0x3fff;
  }

  uint8_t getCeeFraction() const {
    return (_data[0]>>32)&0xff;
  }

  uint16_t getEta() const {
    return _data[1]&0x1ff;
  }

  int16_t getPhi() const {
    uint16_t f((_data[1]>>9)&0x1ff);
    return f<0x100?f:int16_t(f)-0x200;
  }

  uint16_t getZ() const {
    return (_data[1]>>18)&0xfff;
  }

  uint16_t getNumberOfTcs() const {
    return (_data[1]>>32)&0x3ff;
  }

  /*
  int16_t getYOverZ() const {
    uint16_t u((_data[0]>>30)&0xfff);
    return u<0x800?u:int16_t(u)-0x1000;
  }

  uint16_t getLayer() const {
    return (_data[0]>>42)&0x3f;
  }

  int16_t getU() const {
    uint16_t u((_data[0]>>48)&0x7f);
    return u<0x40?u:int16_t(u)-0x80;
  }

  uint16_t getUbar() const {
    return (_data[0]>>55)&0x7f;
  }
  */
  void setEnergy(uint16_t e) {
    if(e>0x3fff) e=0x3fff;
    _data[0]&=0xffffffffffffc000;
    _data[0]|=e;
  }

  void setCeeFraction(uint8_t f) {
    _data[0]&=0xffffff00ffffffff;
    _data[0]|=uint64_t(f)<<32;
  }

  void setEta(uint16_t h) {
    if(h>0x1ff) h=0x1ff;
    _data[1]&=0xfffffffffffffe00;
    _data[1]|=h;
  }

  void setPhi(int16_t f) {
    uint64_t fd;
    if(f> 0x0ff) fd= 0x0ff;
    else if(f>=0) fd=f;
    else if(f<-0x100) fd=0x100;
    else fd=0x200+f;
    _data[1]&=0xfffffffffffc01ff;
    _data[1]|=(fd&0x1ff)<<9;
  }

  void setZ(uint16_t z) {
    if(z>0xfff) z=0xfff;
    _data[1]&=0xffffffffc003ffff;
    _data[1]|=uint64_t(z&0xfff)<<18;
  }

  void setNumberOfTcs(uint16_t n) {
    if(n>0x3ff) n=0x3ff;
    _data[1]&=0xfffffc00ffffffff;
    _data[1]|=uint64_t(n)<<32;
  }

  /*
  void setXOverZ(int16_t x) {
    uint64_t xd;
    if(x> 0x7ff) xd= 0x7ff;
    else if(x>=0) xd=x;
    else if(x<-0x800) xd=0x800;
    else xd=0x1000+x;
    _data[0]&=0xffffffffc003ffff;
    _data[0]|=(xd&0xfff)<<18;
  }

  void setYOverZ(int16_t y) {
    uint64_t yd;
    if(y> 0x7ff) yd= 0x7ff;
    else if(y>=0) yd=y;
    else if(y<-0x800) yd=0x800;
    else yd=0x1000+y;
    _data[0]&=0xfffffc003fffffff;
    _data[0]|=(yd&0xfff)<<30;
  }

  void setLayer(uint16_t l) {
    if(l>47) l=63;
    _data[0]&=0xffff03ffffffffff;
    _data[0]|=uint64_t(l&0x3f)<<42;
  }

  void setU(int16_t u) {
    uint64_t ud;
    if(u> 0x3f) ud= 0x3f;
    else if(u>=0) ud=u;
    else if(u<-0x4f) ud=0x4f;
    else ud=0x80+u;
    _data[0]&=0xff80ffffffffffff;
    _data[0]|=(ud&0xfff)<<48;
  }

  void setUbar(uint16_t u) {
    if(u>0x7f) u=0x7f;
    _data[0]&=0xc07fffffffffffff;
    _data[0]|=uint64_t(u&0x7f)<<55;
  }
  */
  void print() {
    std::cout << "TPGClusterData(" << this << ")::print()" << std::endl;

    for(unsigned w(0);w<4;w++) {
      std::cout << " Word " << w << " "
		<< std::hex << std::setfill('0')
		<< std::setw(16) << _data[w]
		<< std::dec << std::setfill(' ')
		<< std::endl;
    }

    std::cout << " Energy " << getEnergy() << " 2^-2 GeV" << std::endl;
    std::cout << " Eta, phi, z " << getEta() << " x pi/720, "
	      << getPhi() << " x pi/720, "
	      << getZ() << " x 0.05 cm"
	      << std::endl;
  }

private:
  uint64_t _data[4];
};

#endif
