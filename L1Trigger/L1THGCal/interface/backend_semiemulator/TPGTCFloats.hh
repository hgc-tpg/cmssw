#ifndef TPGTCFloats_hh
#define TPGTCFloats_hh

#include <iostream>
#include <iomanip>
#include <cstdint>
#include <cstring>
#include <cassert>

#include "TPGLSBScales.hh"
#include "TPGTCBits.hh"
#include "TMath.h"

// Code from Paul, dated 06/11/2024 (updated 12/04/2025 by Indra)

class TPGTCFloats : public TPGTCBits {
public:
  TPGTCFloats() { setZero(); }
  TPGTCFloats(const TPGTCBits &w) {
    setData(w.getData());
    setCMSSWIndex(w.getCMSSWIndex());
  }

  double getEnergyGeV() const {
    return getEnergy() * lsbScales.LSB_E_TC();
  }

  double getPhiF() const {
    double phi = getPhi()*lsbScales.LSB_phi_TC() ;
    phi -= TMath::Pi()/2.; //shifts the phi
    return phi;
  }

  double getROverZF() const {
    return getROverZ() * lsbScales.LSB_roz_TC();
  }

  double getXOverZF() const {
    return getROverZF() * cos(getPhiF()) ;
  }

  double getYOverZF() const {
    return getROverZF() * sin(getPhiF()) ;
  }

  double getZCm() const {
    uint16_t l(getLayer());
    return lsbScales.getZcm(l);
  }

  uint16_t getZ() const {
    return uint16_t(getZCm()*10.*lsbScales.LSB_z_TC()+0.5);
  }

  void setEnergyGeV(double e) {
    setEnergy(uint32_t(e/lsbScales.LSB_E_TC()+0.5));
  }

  void setROverZF(double roz) {
    setROverZ(uint16_t(roz/lsbScales.LSB_roz_TC()+0.5));
  }

  void setPhiF(double phi) {
    phi += TMath::Pi()/2.;            //to fit the LSB bound of 0 to pi
    if(phi<0.0 or phi>TMath::Pi()) phi = 1.18*TMath::Pi(); //saturate anything else
    setPhi(uint16_t(phi/lsbScales.LSB_phi_TC()+0.5));
  }

  void setROverZPhiF(double x/*xoz*/, double y/*yoz*/, unsigned s) {
    const double s120(0.5*sqrt(3.0)),c120(-0.5);
    double xl = 0.;
    double yl = 0.;  // local reference frame

    switch(s){
    case 0:
      xl = y ; yl = -x;
      break;
    case 1:
      xl = y*c120-x*s120 ; yl = -x*c120-y*s120;
      break;
    case 2:
      xl = y*c120+x*s120 ; yl = -x*c120+y*s120;
      break;
    case 3:
      xl = y ; yl = x;
      break;
    case 4:
      xl = y*c120+x*s120 ; yl = x*c120-y*s120;
      break;
    case 5:
      xl = y*c120-x*s120 ; yl = x*c120+y*s120;
      break;
    default:
      assert(false);
    }

    setROverZF(sqrt(xl*xl+yl*yl));
    setPhiF(atan2(yl,xl));
    //setPhiF(atan(yl/xl));
  }

  void print() {
    std::cout << "TPGTCFloats(" << this << ")::print()"
	      << std::endl << " ";
    TPGTCBits::print();
    std::cout << " Energy " << getEnergyGeV() << " GeV" << std::endl
	      << " x/z, y/z, z = " << getXOverZF()
	      << ", " << getYOverZF()
	      << ", " << getZCm() << std::endl;
  }

private:
  TPGLSBScales::TPGStage2ClusterLSB lsbScales;
};

#endif
