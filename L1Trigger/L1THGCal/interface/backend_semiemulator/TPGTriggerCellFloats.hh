#ifndef TPGTriggerCellFloats_hh
#define TPGTriggerCellFloats_hh

#include <iostream>
#include <iomanip>
#include <cstdint>
#include <cstring>
#include <cassert>

#include "TPGTriggerCellWord.hh"

// Code from Paul, dated 06/11/2024 (updated 10/11/2024)

class TPGTriggerCellFloats : public TPGTriggerCellWord {
public:
  TPGTriggerCellFloats() {
    TPGTriggerCellWord::setZero();
  }

  TPGTriggerCellFloats(const TPGTriggerCellWord &w) : TPGTriggerCellWord(w) {
  }

  double getEnergyGeV() const {
    return TPGTriggerCellWord::getEnergy()/_eScale;
  }

  double getXOverZF() const {
    return TPGTriggerCellWord::getXOverZ()/_xyScale;
  }

  double getYOverZF() const {
    return TPGTriggerCellWord::getYOverZ()/_xyScale;
  }

  double getZCm() const {
    uint16_t l(getLayer());
    return l>0 && l<48?_zLayer[l-1]:0.0;
  }

  void setEnergyGeV(double e) {
    setEnergy(uint32_t(_eScale*e+0.5));
  }

  void setXOverZF(double x) {
    setXOverZ(int16_t(x>=0.0?_xyScale*x+0.5:_xyScale*x-0.5));
  }

  void setYOverZF(double y) {
    setYOverZ(int16_t(y>=0.0?_xyScale*y+0.5:_xyScale*y-0.5));
  }

  void setXYOverZF(double x, double y, unsigned s) {
    static double s120(0.5*sqrt(3.0)),c120(-0.5);

    if(     s==0) {setXOverZF( y            );setYOverZF(-x            );}
    else if(s==1) {setXOverZF( y*c120-x*s120);setYOverZF(-x*c120-y*s120);}
    else if(s==2) {setXOverZF( y*c120+x*s120);setYOverZF(-x*c120+y*s120);}
    else if(s==3) {setXOverZF( y            );setYOverZF( x            );}
    else if(s==4) {setXOverZF( y*c120+x*s120);setYOverZF( x*c120-y*s120);}
    else if(s==5) {setXOverZF( y*c120-x*s120);setYOverZF( x*c120+y*s120);}

    else assert(false);
  }

  void print() {
    std::cout << "TPGTriggerCellFloats(" << this << ")::print()"
	      << std::endl << " ";
    TPGTriggerCellWord::print();
    std::cout << " Energy " << getEnergyGeV() << " GeV" << std::endl
	      << " x/z, y/z, z = " << getXOverZF()
	      << ", " << getYOverZF()
	      << ", " << getZCm() << std::endl;
  }

private:
  static double _zLayer[47];
  static double _eScale;
  static double _xyScale;
};

double TPGTriggerCellFloats::_zLayer[]={
  322.155,     302, 325.212,     304, 328.269,     306, 331.326,     308, 334.383,     310,
  337.440,     312, 340.497,     314, 343.554,     316, 346.611,     318, 349.993,     320,
  353.375,     322, 356.757,     324, 360.139,     326, 367.976,     374, 380.586, 386.891,
  393.196, 399.501, 405.806, 412.111, 418.416, 424.721, 431.026, 439.251, 447.476, 455.701,
  463.926, 472.151, 480.376, 488.406, 496.826, 505.051, 513.276
};

double TPGTriggerCellFloats::_eScale=256.0;
double TPGTriggerCellFloats::_xyScale=2500.0;

#endif
