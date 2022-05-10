#ifndef __L1Trigger_L1THGCal_HGCalHistoClustering_h__
#define __L1Trigger_L1THGCal_HGCalHistoClustering_h__

#include "L1Trigger/L1THGCal/interface/backend_emulator/HGCalTriggerCell_SA.h"
#include "L1Trigger/L1THGCal/interface/backend_emulator/HGCalHistogramCell_SA.h"
#include "L1Trigger/L1THGCal/interface/backend_emulator/CentroidHelper.h"
#include "L1Trigger/L1THGCal/interface/backend_emulator/HGCalHistoClusteringConfig_SA.h"

namespace l1thgcfirmware {

class HGCalHistoClustering {
public:
  HGCalHistoClustering(l1thgcfirmware::ClusterAlgoConfig& config);
  ~HGCalHistoClustering() {}

  void runClustering(const l1thgcfirmware::HGCalTriggerCellSAPtrCollection& triggerCellsIn, const l1thgcfirmware::HGCalHistogramCellSAPtrCollection& histogramIn, l1thgcfirmware::HGCalTriggerCellSAPtrCollection& clusteredTriggerCellsOut, l1thgcfirmware::CentroidHelperPtrCollection& readoutFlagsOut ) const;

private:

  // Clustering
  void clusterizer( const l1thgcfirmware::HGCalTriggerCellSAPtrCollection& triggerCellsIn, const l1thgcfirmware::HGCalHistogramCellSAPtrCollection& histogram, l1thgcfirmware::HGCalTriggerCellSAPtrCollection& clusteredTriggerCells, l1thgcfirmware::HGCalTriggerCellSAPtrCollection& unclusteredTriggerCells, l1thgcfirmware::CentroidHelperPtrCollection& prioritizedMaxima, l1thgcfirmware::CentroidHelperPtrCollection& readoutFlags ) const;

  l1thgcfirmware::ClusterAlgoConfig& config_;
};
}

#endif