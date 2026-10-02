#ifndef TPGClusterFloats_hh
#define TPGClusterFloats_hh

#include <iostream>
#include <iomanip>
#include <cstdint>
#include <cstring>
#include <cassert>

#include "TPGClusterData.hh"

// Code from Paul, dated 10/11/2024

class TPGClusterFloats : public TPGClusterData {
public:
  TPGClusterFloats() {
    TPGClusterData::setZero();
  }

  TPGClusterFloats(const TPGClusterData &d) : TPGClusterData(d) {
  }

  double getEnergyGeV() const {
    return TPGClusterData::getEnergy()/_eScale;
  }

  double getCeeFractionF() const {
    return (TPGClusterData::getCeeFraction()+0.5)/_ceeScale;
  }

  double getLocalPhiRad() const {
    return TPGClusterData::getPhi()/_phiScale;
  }

  double getGlobalPhiRad(unsigned s) const {
    //if(s<3) return -(0.5+2.0*s/3.0)*acos(-1)+getLocalPhiRad();
    //    return (0.5-2.0*(s-3)/3.0)*acos(-1)-getLocalPhiRad();

    int32_t gPhi;
    if(s<3) gPhi= 4*90+4*120*s    +getPhi();
    else    gPhi= 4*90-4*120*(s-3)-getPhi();

    if(gPhi<-720) gPhi+=2*720;
    if(gPhi>=720) gPhi-=2*720;
    return gPhi/_phiScale;
  }

  double getEtaRad() const {
    return (TPGClusterData::getEta()+_etaOffset)/_etaScale;
  }

  double getGlobalEtaRad(unsigned s) const {
    return (s<3) ? (-getEtaRad()) : getEtaRad() ;
  }

  double getLocalXOverZF() const {
    return cos(getLocalPhiRad())/sinh(getEtaRad());
  }

  double getLocalYOverZF() const {
    return sin(getLocalPhiRad())/sinh(getEtaRad());
  }

  double getGlobalXOverZF(unsigned s) const {
    return cos(getGlobalPhiRad(s))/sinh(getEtaRad());
  }

  double getGlobalYOverZF(unsigned s) const {
    return sin(getGlobalPhiRad(s))/sinh(getEtaRad());
  }

  double getZCm() const {
    return _zOffset+TPGClusterData::getZ()/_zScale;
  }

  double getRhoOverZF() const {
    return 1.0/sinh(getEtaRad());
  }



  /*
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
    static double sn(0.5*sqrt(3.0)),cs(-0.5);

    if(     s==0) {setXOverZF( x        );setYOverZF( y        );}
    else if(s==1) {setXOverZF( x*cs+y*sn);setYOverZF( y*cs-x*sn);}
    else if(s==2) {setXOverZF( x*cs-y*sn);setYOverZF( y*cs+x*sn);}
    else if(s==3) {setXOverZF(-x        );setYOverZF( y        );}
    else if(s==4) {setXOverZF(-x*cs+y*sn);setYOverZF( y*cs+x*sn);}
    else if(s==5) {setXOverZF(-x*cs-y*sn);setYOverZF( y*cs-x*sn);}

    else assert(false);
  }

  void print() {
    std::cout << "TPGClusterFloats(" << this << ")::print()"
	      << std::endl << " ";
    TPGClusterData::print();
    std::cout << " Energy " << getEnergyGeV() << " GeV" << std::endl
	      << " x/z, y/z, z = " << getXOverZF()
	      << ", " << getYOverZF()
	      << ", " << getZCm() << std::endl;
  }
  */
private:
  static double _eScale;
  static double _ceeScale;
  static double _phiScale;
  static double _etaScale;
  static double _etaOffset;
  static double _zScale;
  static double _zOffset;
};

double TPGClusterFloats::_eScale=4.0;
double TPGClusterFloats::_ceeScale=256.0;
double TPGClusterFloats::_phiScale=720.0/acos(-1);
double TPGClusterFloats::_etaScale=720.0/acos(-1);
double TPGClusterFloats::_etaOffset=256.0;
double TPGClusterFloats::_zScale=20.0;
double TPGClusterFloats::_zOffset=322.155;

#endif
