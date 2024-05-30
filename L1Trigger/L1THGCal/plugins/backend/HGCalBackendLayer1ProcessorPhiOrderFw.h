#include "L1Trigger/L1THGCal/interface/backend/HGCalBackendLayer1ProcessorPhiOrderFw.h"
#include "DataFormats/ForwardDetId/interface/HGCalTriggerBackendDetId.h"

DEFINE_EDM_PLUGIN(HGCalBackendLayer1Factory,
                  HGCalBackendLayer1ProcessorPhiOrderFw,
                  "HGCalBackendLayer1ProcessorPhiOrderFw");

HGCalBackendLayer1ProcessorPhiOrderFw::HGCalBackendLayer1ProcessorPhiOrderFw(const edm::ParameterSet& conf)
    : HGCalBackendLayer1ProcessorBase(conf), conf_(conf) {
  const edm::ParameterSet& PhiOrderParamConfig = conf.getParameterSet("phiorder_parameters");
  const std::string& PhiOrderWrapperName = PhiOrderParamConfig.getParameter<std::string>("AlgoName");

  clusteringDummy_ = std::make_unique<HGCalClusteringDummyImpl>(conf.getParameterSet("clustering_dummy_parameters"));
  PhiOrderWrapper_ = std::unique_ptr<HGCalLayer1PhiOrderWrapperBase>{
      HGCalLayer1PhiOrderWrapperBaseFactory::get()->create(PhiOrderWrapperName, PhiOrderParamConfig)};
}

void HGCalBackendLayer1ProcessorPhiOrderFw::run(
    const std::pair<uint32_t, std::vector<edm::Ptr<l1t::HGCalTriggerCell>>>& fpga_id_tcs,
    l1t::HGCalClusterBxCollection& clusters) {
  const unsigned sector120 = HGCalTriggerBackendDetId(fpga_id_tcs.first).sector();
  const uint32_t fpga_id = fpga_id_tcs.first;

  // Configuration
  clusteringDummy_->setGeometry(geometry());
  const std::tuple<const HGCalTriggerGeometryBase* const, unsigned, uint32_t> configuration{
      geometry(), sector120, fpga_id};
  PhiOrderWrapper_->configure(configuration);

  PhiOrderWrapper_->process(fpga_id_tcs.second, clusters);
}
