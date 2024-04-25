#ifndef __L1Trigger_L1THGCal_HGCalBackendLayer1ProcessorTruncationFw_h__
#define __L1Trigger_L1THGCal_HGCalBackendLayer1ProcessorTruncationFw_h__

#include "L1Trigger/L1THGCal/interface/HGCalProcessorBase.h"
#include "L1Trigger/L1THGCal/interface/HGCalAlgoWrapperBase.h"

#include "DataFormats/L1THGCal/interface/HGCalTriggerCell.h"
#include "DataFormats/L1THGCal/interface/HGCalCluster.h"

#include "L1Trigger/L1THGCal/interface/backend_emulator/HGCalLayer1TruncationFwImpl.h"
#include "L1Trigger/L1THGCal/interface/backend/HGCalClusteringDummyImpl.h"

class HGCalBackendLayer1ProcessorTruncationFw : public HGCalBackendLayer1ProcessorBase {
public:
  HGCalBackendLayer1ProcessorTruncationFw(const edm::ParameterSet& conf);

  void run(const std::pair<uint32_t, std::vector<edm::Ptr<l1t::HGCalTriggerCell>>>& fpga_id_tcs,
           l1t::HGCalClusterBxCollection& clusters) override;

private:
  std::unique_ptr<HGCalClusteringDummyImpl> clusteringDummy_;
  std::unique_ptr<HGCalLayer1TruncationWrapperBase> truncationWrapper_;
  const edm::ParameterSet conf_;
};

#endif
