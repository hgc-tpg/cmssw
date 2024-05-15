#include "L1Trigger/L1THGCal/interface/backend/HGCalBackendLayer1ProcessorTruncationFw.h"
#include "DataFormats/ForwardDetId/interface/HGCalTriggerBackendDetId.h"

DEFINE_EDM_PLUGIN(HGCalBackendLayer1Factory,
                  HGCalBackendLayer1ProcessorTruncationFw,
                  "HGCalBackendLayer1ProcessorTruncationFw");

HGCalBackendLayer1ProcessorTruncationFw::HGCalBackendLayer1ProcessorTruncationFw(const edm::ParameterSet& conf)
    : HGCalBackendLayer1ProcessorBase(conf), conf_(conf) {
  const edm::ParameterSet& truncationParamConfig = conf.getParameterSet("truncation_parameters");
  const std::string& truncationWrapperName = truncationParamConfig.getParameter<std::string>("AlgoName");

  clusteringDummy_ = std::make_unique<HGCalClusteringDummyImpl>(conf.getParameterSet("clustering_dummy_parameters"));
  truncationWrapper_ = std::unique_ptr<HGCalLayer1TruncationWrapperBase>{
      HGCalLayer1TruncationWrapperBaseFactory::get()->create(truncationWrapperName, truncationParamConfig)};
}

void HGCalBackendLayer1ProcessorTruncationFw::run(
    const std::pair<uint32_t, std::vector<edm::Ptr<l1t::HGCalTriggerCell>>>& fpga_id_tcs,
    l1t::HGCalClusterBxCollection& clusters) {
  const unsigned sector120 = HGCalTriggerBackendDetId(fpga_id_tcs.first).sector();
  const uint32_t fpga_id = fpga_id_tcs.first;

  // Configuration
  clusteringDummy_->setGeometry(geometry());
  const std::tuple<const HGCalTriggerGeometryBase* const, unsigned, uint32_t> configuration{
      geometry(), sector120, fpga_id};
  truncationWrapper_->configure(configuration);

  //std::vector<edm::Ptr<l1t::HGCalTriggerCell>> truncated_tcs;
  truncationWrapper_->process(fpga_id_tcs.second, clusters);
}
