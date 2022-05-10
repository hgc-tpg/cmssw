#ifndef __L1Trigger_L1THGCal_HGCalTCDistribution_h__
#define __L1Trigger_L1THGCal_HGCalTCDistribution_h__

#include "L1Trigger/L1THGCal/interface/backend_emulator/HGCalTriggerCell_SA.h"
#include "L1Trigger/L1THGCal/interface/backend_emulator/DistServer.h"
#include "L1Trigger/L1THGCal/interface/backend_emulator/HGCalHistoClusteringConfig_SA.h"

#include <vector>

namespace l1thgcfirmware {

class HGCalTCDistribution {
public:
  HGCalTCDistribution(l1thgcfirmware::ClusterAlgoConfig& config);
  ~HGCalTCDistribution() {}

  void runTriggerCellDistribution(const l1thgcfirmware::HGCalTriggerCellSAPtrCollections& triggerCellsIn, l1thgcfirmware::HGCalTriggerCellSAPtrCollection& triggerCellsOut ) const;

private:
  // TC input step
  l1thgcfirmware::HGCalTriggerCellSAPtrCollection triggerCellInput( const l1thgcfirmware::HGCalTriggerCellSAPtrCollections& inputs ) const;

  // TC distribution steps
  void triggerCellDistribution0( l1thgcfirmware::HGCalTriggerCellSAPtrCollection& triggerCellsIn ) const;
  l1thgcfirmware::HGCalTriggerCellSAPtrCollections triggerCellDistribution1( l1thgcfirmware::HGCalTriggerCellSAPtrCollection& triggerCellsIn ) const;
  l1thgcfirmware::HGCalTriggerCellSAPtrCollections triggerCellDistribution2( l1thgcfirmware::HGCalTriggerCellSAPtrCollection& triggerCellsIn, l1thgcfirmware::HGCalTriggerCellSAPtrCollections& inTriggerCellDistributionGrid ) const;
  l1thgcfirmware::HGCalTriggerCellSAPtrCollections triggerCellDistribution3( l1thgcfirmware::HGCalTriggerCellSAPtrCollection& triggerCellsIn, l1thgcfirmware::HGCalTriggerCellSAPtrCollections& inTriggerCellDistributionGrid ) const;
  void triggerCellDistribution4( l1thgcfirmware::HGCalTriggerCellSAPtrCollection& triggerCellsIn ) const;
  void triggerCellDistribution5( l1thgcfirmware::HGCalTriggerCellSAPtrCollection& triggerCellsIn, l1thgcfirmware::HGCalTriggerCellSAPtrCollections& inTriggerCellDistributionGrid ) const;

  // Useful functions
  void initializeTriggerCellDistGrid( l1thgcfirmware::HGCalTriggerCellSAPtrCollections& grid, unsigned int nX, unsigned int nY ) const;

  void runDistServers( const l1thgcfirmware::HGCalTriggerCellSAPtrCollections& gridIn,
                       l1thgcfirmware::HGCalTriggerCellSAPtrCollections& gridOut,
                       l1thgcfirmware::HGCalTriggerCellSAPtrCollection& tcsOut,
                       unsigned int latency,
                       unsigned int nDistServers,
                       unsigned int nInputs,
                       unsigned int nOutputs,
                       unsigned int nInterleave,
                      bool setOutputGrid ) const;

  l1thgcfirmware::ClusterAlgoConfig& config_;
};
}

#endif