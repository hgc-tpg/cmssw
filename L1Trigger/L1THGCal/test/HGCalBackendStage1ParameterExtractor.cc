#include <iostream> // std::cout
#include <fstream>  // std::ofstream

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/stream/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/ForwardDetId/interface/ForwardSubdetector.h"
#include "DataFormats/ForwardDetId/interface/HGCalDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCScintillatorDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCSiliconDetId.h"
#include "DataFormats/ForwardDetId/interface/HFNoseTriggerDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCalTriggerDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCalTriggerModuleDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCalTriggerBackendDetId.h"

#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "L1Trigger/L1THGCal/interface/HGCalTriggerGeometryBase.h"
#include "L1Trigger/L1THGCal/interface/HGCalTriggerTools.h"

#include <nlohmann/json.hpp>
using json = nlohmann::ordered_json;  // using ordered_json for readability

class HGCalBackendStage1ParameterExtractor : public edm::stream::EDAnalyzer<> {
public:
  explicit HGCalBackendStage1ParameterExtractor(const edm::ParameterSet&);
  ~HGCalBackendStage1ParameterExtractor();

  virtual void beginRun(const edm::Run&, const edm::EventSetup&);
  virtual void analyze(const edm::Event&, const edm::EventSetup&);

private:
  void fillTriggerGeometry(json& json_file);
  std::tuple<uint32_t,uint32_t,uint32_t> getTCaddress(int& tc_ueta, int& tc_vphi, bool& isScint);

  HGCalTriggerTools triggerTools_;

  edm::ESHandle<HGCalTriggerGeometryBase> triggerGeometry_;
  edm::ESGetToken<HGCalTriggerGeometryBase, CaloGeometryRecord> triggerGeomToken_;

  // Metadata
  std::string detector_version_;
  std::string cmssw_era_;

  uint32_t the_fpga_;

  // geometry data
  std::vector<uint32_t> disconnected_layers_;
  std::string json_mapping_file_;
  uint32_t scintillator_trigger_cell_size_;
  std::string trigger_geom_;

  // truncation parameters
  double roz_min_;
  double roz_max_;
  uint32_t roz_bins_;
  std::vector<uint32_t> max_tcs_per_bins_;
  std::vector<double> phi_edges_;

  // TC map parameters:
  std::map<std::pair<uint32_t, uint32_t>, uint32_t> tc_coord_uv_;

  // output json name:
  std::string outJSONname_;

  typedef std::unordered_map<uint32_t, std::unordered_set<uint32_t>> trigger_map_set;
};

HGCalBackendStage1ParameterExtractor::HGCalBackendStage1ParameterExtractor(const edm::ParameterSet& conf)
    : triggerGeomToken_(esConsumes<HGCalTriggerGeometryBase, CaloGeometryRecord, edm::Transition::BeginRun>())

{
  // get name of output JSON config file
  outJSONname_ = conf.getParameter<std::string>("outJSONname");

  // get ID of tested FPGA
  the_fpga_ = conf.getParameter<uint32_t>("testedFpga");

  // get meta data
  const edm::ParameterSet& metaData = conf.getParameterSet("MetaData");
  cmssw_era_ = "dummy1";//metaData.getParameter<std::string>("CMSSW_era");
  detector_version_ = "dummy2";//metaData.getParameter<std::string>("detector_version_");

  // get geometry configuration
  const edm::ParameterSet& triggerGeom = conf.getParameterSet("TriggerGeometryParam");
  disconnected_layers_ = triggerGeom.getParameter<std::vector<uint32_t>>("DisconnectedLayers");
  json_mapping_file_ = triggerGeom.getParameter<edm::FileInPath>("JsonMappingFile").relativePath();
  scintillator_trigger_cell_size_ = triggerGeom.getParameter<uint32_t>("ScintillatorTriggerCellSize");
  trigger_geom_ = triggerGeom.getParameter<std::string>("TriggerGeometryName");

  // get TCid -> uv coordinate correspondance
  const edm::ParameterSet& TCcoord_uv = conf.getParameterSet("TCcoord_UV");
  std::vector<uint32_t> tc_coord_u = TCcoord_uv.getParameter<std::vector<uint32_t>>("TCu");
  std::vector<uint32_t> tc_coord_v = TCcoord_uv.getParameter<std::vector<uint32_t>>("TCv");
  if (tc_coord_u.size() != tc_coord_v.size())
    throw cms::Exception("BadParameter") << "TCu and TCv vectors should be of same size";
  for (size_t i = 0; i < tc_coord_u.size(); ++i)
    tc_coord_uv_.emplace(std::make_pair(tc_coord_u.at(i), tc_coord_v.at(i)), i);

  // Get truncation parameters
  const edm::ParameterSet& truncationParamConfig = conf.getParameterSet("BackendStage1Params");
  roz_min_ = truncationParamConfig.getParameter<double>("rozMin");
  roz_max_ = truncationParamConfig.getParameter<double>("rozMax");
  roz_bins_ = truncationParamConfig.getParameter<uint32_t>("rozBins");
  max_tcs_per_bins_ = truncationParamConfig.getParameter<std::vector<uint32_t>>("maxTcsPerBin");
  phi_edges_ = truncationParamConfig.getParameter<std::vector<double>>("phiSectorEdges");
}

HGCalBackendStage1ParameterExtractor::~HGCalBackendStage1ParameterExtractor() {}

void HGCalBackendStage1ParameterExtractor::beginRun(const edm::Run& /*run*/, const edm::EventSetup& es) {
  triggerGeometry_ = es.getHandle(triggerGeomToken_);

  json outJSON;

  // fill MetaData
  outJSON["MetaData"]["CMSSWEra"] = cmssw_era_;
  outJSON["MetaData"]["DetectorVersion"] = detector_version_;
  outJSON["MetaData"]["fpgaId"] = the_fpga_;

  // fill geometry configuration
  outJSON["TriggerGeometryConfig"]["TriggerGeometryName"] = trigger_geom_;
  outJSON["TriggerGeometryConfig"]["DisconnectedLayers"] = disconnected_layers_;
  outJSON["TriggerGeometryConfig"]["JsonMappingFile"] = json_mapping_file_;
  outJSON["TriggerGeometryConfig"]["ScintillatorTriggerCellSize"] = scintillator_trigger_cell_size_;

  // fill truncation parameters
  outJSON["TruncationConfig"]["rozMin"] = roz_min_;
  outJSON["TruncationConfig"]["rozMax"] = roz_max_;
  outJSON["TruncationConfig"]["rozBins"] = roz_bins_;
  outJSON["TruncationConfig"]["maxTcsPerBin"] = max_tcs_per_bins_;
  outJSON["TruncationConfig"]["phiSectorEdges"] = phi_edges_;

  // fill trigger geometry
  fillTriggerGeometry(outJSON);

  // write out JSON file
  std::ofstream outputfile(outJSONname_.c_str());
  outputfile << std::setw(4) << outJSON << std::endl;
}

void HGCalBackendStage1ParameterExtractor::fillTriggerGeometry(json& json_file) {
  std::vector<uint32_t> trigger_cells;

  // retrieve valid trigger cells
  for (const auto& id : triggerGeometry_->eeGeometry()->getValidDetIds()) {
    HGCSiliconDetId detid(id);
    if (triggerGeometry_->eeTopology().valid(id))
      trigger_cells.push_back(triggerGeometry_->getTriggerCellFromCell(id));
  }
  for (const auto& id : triggerGeometry_->hsiGeometry()->getValidDetIds()) {
    HGCSiliconDetId detid(id);
    if (triggerGeometry_->hsiTopology().valid(id))
      trigger_cells.push_back(triggerGeometry_->getTriggerCellFromCell(id));
  }
  for (const auto& id : triggerGeometry_->hscGeometry()->getValidDetIds()) {
    trigger_cells.push_back(triggerGeometry_->getTriggerCellFromCell(id));
  }

  // loop over trigger cells
  edm::LogPrint("JSONFilling") << "Filling JSON map";
  for (const auto& triggercell : trigger_cells) {
    DetId id(triggercell);
    // get module ID and check if relevant module (sector0, zside=1, and connected to FPGA)
    uint32_t moduleId = triggerGeometry_->getModuleFromTriggerCell(id);
    if (moduleId == 0 || triggerGeometry_->disconnectedModule(moduleId))
      continue;
    HGCalTriggerModuleDetId tc_module(moduleId);
    if (!(tc_module.isHGCalModuleDetId()) || (tc_module.zside() < 0) || (tc_module.sector() != 0))
      continue;

    // only retrieve mapping for the tested fpga
    uint32_t fpgaId = triggerGeometry_->getStage1FpgaFromModule(moduleId);
    HGCalTriggerBackendDetId tc_fpga(fpgaId);
    if (!(tc_fpga.isStage1FPGA()) || (tc_fpga.sector() != 0) || (tc_fpga.zside() < 0))
      continue;
    if ((uint32_t)tc_fpga.label() != the_fpga_)
      continue;

    // retrieve information to be saved
    int triggerCellSubdet = 0;
    int triggerCellLayer = 0;
    int triggerCellUEta = 0;
    int triggerCellVPhi = 0;

    bool isScintillatorCell = triggerTools_.isScintillator(id);
    if (isScintillatorCell) {
      HGCScintillatorDetId id_sc(triggercell);
      triggerCellSubdet = id_sc.subdet();
      triggerCellLayer = id_sc.layer();
      triggerCellUEta = id_sc.ietaAbs();
      triggerCellVPhi = id_sc.iphi();
    } else {
      HGCalTriggerDetId id_si_trig(triggercell);
      triggerCellSubdet = id_si_trig.subdet();
      triggerCellLayer = id_si_trig.layer();
      triggerCellUEta = id_si_trig.triggerCellU();
      triggerCellVPhi = id_si_trig.triggerCellV();
    }
    GlobalPoint position = triggerGeometry_->getTriggerCellPosition(triggercell);
    float triggerCellX = position.x();
    float triggerCellY = position.y();
    float triggerCellZ = position.z();
    float triggerCellEta = position.eta();
    float triggerCellPhi = position.phi();
    float triggerCellRoverZ = sqrt(triggerCellX * triggerCellX + triggerCellY * triggerCellY) / triggerCellZ;
    
    std::tuple<uint32_t,uint32_t,uint32_t> tcAddressAndCoord = getTCaddress(triggerCellUEta, triggerCellVPhi, isScintillatorCell);

    // save TC info into JSON
    std::string strModId = std::to_string(moduleId);
    std::string strTCAddr = std::to_string(std::get<0>(tcAddressAndCoord));
    json_file["TriggerCellMap"]["module_hash"][strModId]["tc_id"][strTCAddr]["subdet"] = triggerCellSubdet;
    json_file["TriggerCellMap"]["module_hash"][strModId]["tc_id"][strTCAddr]["layer"] = triggerCellLayer;
    json_file["TriggerCellMap"]["module_hash"][strModId]["tc_id"][strTCAddr]["ueta"] = std::get<1>(tcAddressAndCoord);
    json_file["TriggerCellMap"]["module_hash"][strModId]["tc_id"][strTCAddr]["vphi"] = std::get<2>(tcAddressAndCoord);
    json_file["TriggerCellMap"]["module_hash"][strModId]["tc_id"][strTCAddr]["x"] = triggerCellX;
    json_file["TriggerCellMap"]["module_hash"][strModId]["tc_id"][strTCAddr]["y"] = triggerCellY;
    json_file["TriggerCellMap"]["module_hash"][strModId]["tc_id"][strTCAddr]["z"] = triggerCellZ;
    json_file["TriggerCellMap"]["module_hash"][strModId]["tc_id"][strTCAddr]["roz"] = triggerCellRoverZ;
    json_file["TriggerCellMap"]["module_hash"][strModId]["tc_id"][strTCAddr]["eta"] = triggerCellEta;
    json_file["TriggerCellMap"]["module_hash"][strModId]["tc_id"][strTCAddr]["phi"] = triggerCellPhi;
  }
}

std::tuple<uint32_t,uint32_t,uint32_t> HGCalBackendStage1ParameterExtractor::getTCaddress(int& tc_ueta, int& tc_vphi, bool& isScintillator) {

    // transform HGCalHSc TC coordinates to define a TC address in [0,47]
    if (isScintillator) {  // HGCalHSc
      tc_vphi = (tc_vphi - 1) % 4;
      if (tc_ueta <= 3) {
        tc_ueta = tc_ueta;
      } else if (tc_ueta <= 9) {
        tc_ueta = tc_ueta - 4;
      } else if (tc_ueta <= 13) {
        tc_ueta = tc_ueta - 10;
      } else if (tc_ueta <= 17) {
        tc_ueta = tc_ueta - 14;
      } else {
        tc_ueta = tc_ueta - 18;
      }
    }

    // attribute ID to TC according to subdetector
    uint32_t tcaddress = 99;
    if (isScintillator) {  //HGCalHSc(10)
      tcaddress = (tc_ueta << 2) + tc_vphi;
    } else {  //HGCalHSiTrigger(2) or HGCalEE(1)
      tcaddress = tc_coord_uv_.find(std::make_pair(tc_ueta, tc_vphi))->second;
    }
    return std::make_tuple(tcaddress, tc_ueta, tc_vphi);

}

void HGCalBackendStage1ParameterExtractor::analyze(const edm::Event& e, const edm::EventSetup& es) {}

// define this as a plug-in
DEFINE_FWK_MODULE(HGCalBackendStage1ParameterExtractor);
