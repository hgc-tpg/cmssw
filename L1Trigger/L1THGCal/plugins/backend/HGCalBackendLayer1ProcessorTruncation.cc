#include "L1Trigger/L1THGCal/interface/backend/HGCalBackendLayer1ProcessorTruncation.h"

DEFINE_EDM_PLUGIN(HGCalBackendLayer1Factory, HGCalBackendLayer1ProcessorTruncation, "HGCalBackendLayer1ProcessorTruncation");

HGCalBackendLayer1ProcessorTruncation::HGCalBackendLayer1ProcessorTruncation(const edm::ParameterSet& conf)
    : HGCalBackendLayer1ProcessorBase(conf) {
  clusteringDummy_ = std::make_unique<HGCalClusteringDummyImpl>(conf.getParameterSet("clustering_dummy_parameters"));
  truncation_ = std::make_unique<HGCalLayer1TruncationImpl>(conf.getParameterSet("truncation_parameters"));
}

void HGCalBackendLayer1ProcessorTruncation::run(
    const std::pair<uint32_t, std::vector<edm::Ptr<l1t::HGCalTriggerCell>>>& fpga_id_tcs,
                                      l1t::HGCalClusterBxCollection& clusters) {
  if (clusteringDummy_)
    clusteringDummy_->setGeometry(geometry());
  if (truncation_)
    truncation_->setGeometry(geometry());

  std::vector<edm::Ptr<l1t::HGCalTriggerCell>> truncated_tcs;
  truncation_->run(fpga_id_tcs.first, fpga_id_tcs.second, truncated_tcs);
  clusteringDummy_->clusterizeDummy(truncated_tcs, clusters);
}
