/**********************************************************************
 Created on : 22/04/2025
 Purpose    : HGCAL TPG Cluster information
 Author     : Indranil Das, Research Associate, Imperial
 Email      : indranil.das@cern.ch | indra.ehep@gmail.com
**********************************************************************/
#ifndef TPGCluster_h
#define TPGCluster_h

#include "TPGLSBScales.hh"
#include "TPGBEDataformat.hh"

#ifdef CMSSW_GIT_HASH
#include "L1Trigger/L1THGCal/interface/backend_semiemulator/HGCalCluster_HW.h"
#else
#include "CMSSWCode/DataFormats/L1THGCal/interface/HGCalCluster_HW.h"
#endif

//Class to access cluster information
class TPGCluster{
public:
  TPGCluster() : maxFinderPass_(0), cmsswIndices_() {}
  TPGCluster(const l1thgcfirmware::HGCalCluster_HW *clus) : cldata(*clus), cmsswIndices_() {}
  void setHGCalCluster_HW(const l1thgcfirmware::HGCalCluster_HW *clus) { cldata = *clus ; }
  const l1thgcfirmware::HGCalCluster_HW& getClData() const { return cldata; }

  void setMaxFinderPass(unsigned int pass) { maxFinderPass_ = pass; }
  unsigned int getMaxFinderPass() const { return maxFinderPass_; }

  uint32_t getEnergy() const { return cldata.e.to_uint(); }
  double getEnergyGeV() const { return getEnergy() * lsbScales.LSB_E(); }
  double getCeeFractionF() const { return cldata.fractionInCE_E.to_uint() * lsbScales.LSB_Frac(); }

  int getLocalPhi() const { return cldata.w_phi.to_int(); }
  double getLocalPhiRad() const { return getLocalPhi() * lsbScales.LSB_phi(); }

  double getGlobalPhiRad(unsigned s) const {
    int32_t gPhi;
    if(s<3) gPhi= 4*90+4*120*s    +getLocalPhi();
    else    gPhi= 4*90-4*120*(s-3)-getLocalPhi();

    if(gPhi<-720) gPhi+=2*720;
    if(gPhi>=720) gPhi-=2*720;
    return gPhi*lsbScales.LSB_phi();
  }

  uint32_t getEta() const { return cldata.w_eta.to_uint(); }
  double getEtaRad() const { return getEta() * lsbScales.LSB_eta(); }
  double getGlobalEtaRad(unsigned s) const { return (s<3) ? (-getEtaRad()) : getEtaRad() ;  }

  double getLocalXOverZF() const { return cos(getLocalPhiRad())/sinh(getEtaRad()); }
  double getLocalYOverZF() const { return sin(getLocalPhiRad())/sinh(getEtaRad()); }
  double getGlobalXOverZF(unsigned s) const { return cos(getGlobalPhiRad(s))/sinh(getGlobalEtaRad(s)); }
  double getGlobalYOverZF(unsigned s) const { return sin(getGlobalPhiRad(s))/sinh(getGlobalEtaRad(s)); }

  uint32_t getZ() const { return cldata.w_z.to_uint(); }
  double getZCm() const { return getZ() * lsbScales.LSB_z(); }

  double getRhoOverZF() const { return 1.0/sinh(getEtaRad()); }
  double getGlobalRhoOverZF(unsigned s) const { return 1.0/sinh(getGlobalEtaRad(s)); }

  // Getter and setter for cmsswIndices_
  const std::vector<int>& getCMSSWIndices() const { return cmsswIndices_; }
  void setCMSSWIndices(const std::vector<int>& indices) { cmsswIndices_ = indices; }

private:
  l1thgcfirmware::HGCalCluster_HW cldata;
  TPGLSBScales::TPGStage2ClusterLSB lsbScales;

  unsigned int maxFinderPass_; // Which pass of the local maxima finding this cluster was found in

  std::vector<int> cmsswIndices_;
};

// double TPGClusterFloats::_eScale=4.0;
// double TPGClusterFloats::_ceeScale=256.0;
// double TPGClusterFloats::_phiScale=720.0/acos(-1);
// double TPGClusterFloats::_etaScale=720.0/acos(-1);
// double TPGClusterFloats::_etaOffset=256.0;
// double TPGClusterFloats::_zScale=20.0;
// double TPGClusterFloats::_zOffset=322.155;

#endif
