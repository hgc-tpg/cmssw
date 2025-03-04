#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/ESGetToken.h"

#include "DataFormats/L1THGCal/interface/HGCalTriggerCell.h"
#include "DataFormats/L1THGCal/interface/HGCalTriggerSums.h"
#include "DataFormats/HGCDigi/interface/HGCDigiCollections.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "L1Trigger/L1THGCal/interface/HGCalTriggerGeometryBase.h"
#include "DataFormats/ForwardDetId/interface/HGCalTriggerBackendDetId.h"

#include "L1Trigger/L1THGCal/interface/HGCalProcessorBase.h"

#include "DataFormats/L1THGCal/interface/HGCalCluster.h"

#include <memory>

class HGCalBackendLayer1Producer : public edm::stream::EDProducer<> {
public:
  HGCalBackendLayer1Producer(const edm::ParameterSet&);
  ~HGCalBackendLayer1Producer() override {}

  void beginRun(const edm::Run&, const edm::EventSetup&) override;
  void produce(edm::Event&, const edm::EventSetup&) override;

private:
  // inputs
  edm::EDGetToken input_cell_, input_sums_;
  edm::ESHandle<HGCalTriggerGeometryBase> triggerGeometry_;
  edm::ESGetToken<HGCalTriggerGeometryBase, CaloGeometryRecord> triggerGeomToken_;

  std::unique_ptr<HGCalBackendLayer1ProcessorBase> backendProcess_;
  bool bypass_be_mapping_ = false;
};

DEFINE_FWK_MODULE(HGCalBackendLayer1Producer);

HGCalBackendLayer1Producer::HGCalBackendLayer1Producer(const edm::ParameterSet& conf)
    : input_cell_(consumes<l1t::HGCalTriggerCellBxCollection>(conf.getParameter<edm::InputTag>("InputTriggerCells"))),
      triggerGeomToken_(esConsumes<HGCalTriggerGeometryBase, CaloGeometryRecord, edm::Transition::BeginRun>()),
      bypass_be_mapping_(conf.getParameter<bool>("BypassBackendMapping")) {
  //setup Backend parameters
  const edm::ParameterSet& beParamConfig = conf.getParameterSet("ProcessorParameters");
  const std::string& beProcessorName = beParamConfig.getParameter<std::string>("ProcessorName");
  backendProcess_ = std::unique_ptr<HGCalBackendLayer1ProcessorBase>{
      HGCalBackendLayer1Factory::get()->create(beProcessorName, beParamConfig)};

  produces<l1t::HGCalClusterBxCollection>(backendProcess_->name());
}

void HGCalBackendLayer1Producer::beginRun(const edm::Run& /*run*/, const edm::EventSetup& es) {
  triggerGeometry_ = es.getHandle(triggerGeomToken_);
  backendProcess_->setGeometry(triggerGeometry_.product());
}

void HGCalBackendLayer1Producer::produce(edm::Event& e, const edm::EventSetup& es) {
  // Output collections
  auto be_cluster_output = std::make_unique<l1t::HGCalClusterBxCollection>();

  // Input collections
  edm::Handle<l1t::HGCalTriggerCellBxCollection> trigCellBxColl;
  e.getByToken(input_cell_, trigCellBxColl);

  // Split trigger cell collection per FPGA
  std::unordered_map<uint32_t, std::vector<edm::Ptr<l1t::HGCalTriggerCell>>> tcs_per_fpga;

  for (unsigned i = 0; i < trigCellBxColl->size(); ++i) {
    edm::Ptr<l1t::HGCalTriggerCell> tc_ptr(trigCellBxColl, i);
    uint32_t fpga = 0;  // null detid
    if (!bypass_be_mapping_) {
      uint32_t module = triggerGeometry_->getModuleFromTriggerCell(tc_ptr->detId());
      fpga = triggerGeometry_->getStage1FpgaFromModule(module);
    }
    // if bypass_be_mapping, TCs from the entire detector will be put in a single collection
    tcs_per_fpga[fpga].push_back(tc_ptr);
  }

  for (auto& fpga_tcs : tcs_per_fpga) {
    backendProcess_->run(fpga_tcs, *be_cluster_output);
  }

  e.put(std::move(be_cluster_output), backendProcess_->name());
}
