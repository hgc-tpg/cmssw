#ifndef __L1Trigger_L1THGCal_HGCalVFELinearizationImpl_h__
#define __L1Trigger_L1THGCal_HGCalVFELinearizationImpl_h__

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/HGCDigi/interface/HGCDigiCollections.h"
#include "L1Trigger/L1THGCal/interface/HGCalTriggerTools.h"
#include "SimCalorimetry/HGCalSimAlgos/interface/HGCalSiNoiseMap.h"
#include "SimCalorimetry/HGCalSimAlgos/interface/HGCalSciNoiseMap.h"

#include <vector>
#include <utility>

class HGCalVFELinearizationImpl {
public:
  HGCalVFELinearizationImpl(const edm::ParameterSet& conf, DetId::Detector det);
  void setGeometry(const HGCalTriggerGeometryBase* const geom);

  void linearize(const std::vector<HGCalDataFrame>&, std::vector<std::pair<DetId, uint32_t>>&);

private:
  DetId::Detector detector_;
  //
  double adcLSB_;
  double linLSB_;
  double adcsaturation_;
  uint32_t tdcnBits_;
  double tdcOnset_;
  uint32_t adcnBits_;
  double tdcsaturation_;
  double tdcLSB_;
  //
  uint32_t linMax_;
  uint32_t linnBits_;
  std::vector<double> oot_coefficients_;
  //
  HGCalTriggerTools triggerTools_;
  bool old_digi_ = false;
  mutable HGCalSiNoiseMap<HGCSiliconDetId> noise_map_;
};

#endif
