#include <iostream>
#include <fstream>  // std::ofstream
#include <vector>
#include <cstdlib>
#include <any>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/ESTransientHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/ForwardDetId/interface/ForwardSubdetector.h"
#include "DataFormats/ForwardDetId/interface/HFNoseTriggerDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCalDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCalTriggerDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCScintillatorDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCSiliconDetId.h"

#include "Geometry/Records/interface/CaloGeometryRecord.h"

#include "L1Trigger/L1THGCal/interface/HGCalTriggerGeometryBase.h"
#include "DataFormats/ForwardDetId/interface/HGCalTriggerModuleDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCalTriggerBackendDetId.h"

#include <nlohmann/json.hpp>
using json = nlohmann::ordered_json;  // using ordered_json for readability

class HGCalBackendStage1ParameterExtractor : public edm::stream::EDAnalyzer<> {
public:
  explicit HGCalBackendStage1ParameterExtractor(const edm::ParameterSet&);
  ~HGCalBackendStage1ParameterExtractor();

  virtual void beginRun(const edm::Run&, const edm::EventSetup&);
  virtual void analyze(const edm::Event&, const edm::EventSetup&);

private:
  bool checkMappingConsistency();
  void fillMetaData(json& json_file);
  void fillMethodConfig(json& json_file);
  void fillTriggerGeometry(json& json_file);

  edm::ESHandle<HGCalTriggerGeometryBase> triggerGeometry_;
  edm::ESGetToken<HGCalTriggerGeometryBase, CaloGeometryRecord> triggerGeomToken_;

  // Metadata
  std::string cmssw_version_;
  int the_fpga_;

  // geometry data
  std::vector<uint32_t> disconnected_layers_;
  std::vector<uint32_t> disconnected_modules_;
  std::string json_mapping_file_;
  std::string L1T_links_mapping_;
  std::string L1T_modules_mapping_;
  uint32_t scintillator_links_per_module_;
  uint32_t scintillator_module_size_;
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

  // mapping consistency check result
  bool no_trigger_ = false;

  typedef std::unordered_map<uint32_t, std::unordered_set<uint32_t>> trigger_map_set;
};

HGCalBackendStage1ParameterExtractor::HGCalBackendStage1ParameterExtractor(const edm::ParameterSet& conf)
    : triggerGeomToken_(esConsumes<HGCalTriggerGeometryBase, CaloGeometryRecord, edm::Transition::BeginRun>())

{
  // get name of output JSON config file
  outJSONname_ = conf.getParameter<std::string>("outJSONname");

  // get ID of tested FPGA
  the_fpga_ = conf.getParameter<int>("testedFpga");

  // get meta data
  const edm::ParameterSet& metaData = conf.getParameterSet("MetaData");
  cmssw_version_ = metaData.getParameter<std::string>("CMSSW_version");

  // get geometry configuration
  const edm::ParameterSet& triggerGeom = conf.getParameterSet("TriggerGeometryParam");
  disconnected_layers_ = triggerGeom.getParameter<std::vector<uint32_t>>("DisconnectedLayers");
  disconnected_modules_ = triggerGeom.getParameter<std::vector<uint32_t>>("DisconnectedModules");
  json_mapping_file_ = triggerGeom.getParameter<edm::FileInPath>("JsonMappingFile").relativePath();
  L1T_links_mapping_ = triggerGeom.getParameter<edm::FileInPath>("L1TLinksMapping").relativePath();
  L1T_modules_mapping_ = triggerGeom.getParameter<edm::FileInPath>("L1TModulesMapping").relativePath();
  scintillator_links_per_module_ = triggerGeom.getParameter<uint32_t>("ScintillatorLinksPerModule");
  scintillator_module_size_ = triggerGeom.getParameter<uint32_t>("ScintillatorModuleSize");
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
  const edm::ParameterSet& truncationParamConfig =
      conf.getParameterSet("BackendStage1Params").getParameterSet("truncation_parameters");
  roz_min_ = truncationParamConfig.getParameter<double>("rozMin");
  roz_max_ = truncationParamConfig.getParameter<double>("rozMax");
  roz_bins_ = truncationParamConfig.getParameter<uint32_t>("rozBins");
  max_tcs_per_bins_ = truncationParamConfig.getParameter<std::vector<uint32_t>>("maxTcsPerBin");
  phi_edges_ = truncationParamConfig.getParameter<std::vector<double>>("phiSectorEdges");
}

HGCalBackendStage1ParameterExtractor::~HGCalBackendStage1ParameterExtractor() {}

void HGCalBackendStage1ParameterExtractor::beginRun(const edm::Run& /*run*/, const edm::EventSetup& es) {
  triggerGeometry_ = es.getHandle(triggerGeomToken_);

  // note: identical function in geometry tester; should it be imported instead?
  no_trigger_ = !checkMappingConsistency();

  json outJSON;

  // fill MetaData
  outJSON["MetaData"]["CMSSWversion"] = cmssw_version_;
  outJSON["MetaData"]["fpgaId"] = the_fpga_;

  // fill geometry configuration
  outJSON["TriggerGeometryConfig"]["TriggerGeometryName"] = trigger_geom_;
  outJSON["TriggerGeometryConfig"]["DisconnectedLayers"] = disconnected_layers_;
  outJSON["TriggerGeometryConfig"]["DisconnectedModules"] = disconnected_modules_;
  outJSON["TriggerGeometryConfig"]["JsonMappingFile"] = json_mapping_file_;
  outJSON["TriggerGeometryConfig"]["L1TLinksMapping"] = L1T_links_mapping_;
  outJSON["TriggerGeometryConfig"]["L1TModulesMapping"] = L1T_modules_mapping_;
  outJSON["TriggerGeometryConfig"]["ScintillatorLinksPerModule"] = scintillator_links_per_module_;
  outJSON["TriggerGeometryConfig"]["ScintillatorModuleSize"] = scintillator_module_size_;
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

bool HGCalBackendStage1ParameterExtractor::checkMappingConsistency() {
  try {
    // Set of (subdet,layer,waferU,waferV) with module mapping errors
    std::set<std::tuple<uint32_t, uint32_t, int, int>> module_errors;
    trigger_map_set modules_to_triggercells;
    trigger_map_set modules_to_cells;
    trigger_map_set triggercells_to_cells;
    // EE
    for (const auto& id : triggerGeometry_->eeGeometry()->getValidDetIds()) {
      HGCSiliconDetId detid(id);
      if (!triggerGeometry_->eeTopology().valid(id))
        continue;
      // fill trigger cells
      uint32_t trigger_cell = triggerGeometry_->getTriggerCellFromCell(id);
      auto itr_insert = triggercells_to_cells.emplace(trigger_cell, std::unordered_set<uint32_t>());
      itr_insert.first->second.emplace(id);
      // fill modules
      uint32_t module = 0;
      try {
        module = triggerGeometry_->getModuleFromCell(id);
        triggerGeometry_->getLinksInModule(module);
      } catch (const std::exception& e) {
        module_errors.emplace(std::make_tuple(HGCalTriggerModuleDetId(module).triggerSubdetId(),
                                              HGCalTriggerModuleDetId(module).layer(),
                                              HGCalTriggerModuleDetId(module).moduleU(),
                                              HGCalTriggerModuleDetId(module).moduleV()));
        continue;
      }
      itr_insert = modules_to_cells.emplace(module, std::unordered_set<uint32_t>());
      itr_insert.first->second.emplace(id);
    }

    // HSi
    for (const auto& id : triggerGeometry_->hsiGeometry()->getValidDetIds()) {
      if (!triggerGeometry_->hsiTopology().valid(id))
        continue;
      // fill trigger cells
      uint32_t trigger_cell = triggerGeometry_->getTriggerCellFromCell(id);
      auto itr_insert = triggercells_to_cells.emplace(trigger_cell, std::unordered_set<uint32_t>());
      itr_insert.first->second.emplace(id);
      // fill modules
      uint32_t module = 0;
      try {
        module = triggerGeometry_->getModuleFromCell(id);
        triggerGeometry_->getLinksInModule(module);
      } catch (const std::exception& e) {
        module_errors.emplace(std::make_tuple(HGCalTriggerModuleDetId(module).triggerSubdetId(),
                                              HGCalTriggerModuleDetId(module).layer(),
                                              HGCalTriggerModuleDetId(module).moduleU(),
                                              HGCalTriggerModuleDetId(module).moduleV()));
        continue;
      }
      itr_insert = modules_to_cells.emplace(module, std::unordered_set<uint32_t>());
      itr_insert.first->second.emplace(id);
    }

    // HSc
    for (const auto& id : triggerGeometry_->hscGeometry()->getValidDetIds()) {
      // fill trigger cells
      uint32_t layer = HGCScintillatorDetId(id).layer();
      if (HGCScintillatorDetId(id).type() != triggerGeometry_->hscTopology().dddConstants().getTypeTrap(layer)) {
        std::cout << "Sci cell type = " << HGCScintillatorDetId(id).type()
                  << " != " << triggerGeometry_->hscTopology().dddConstants().getTypeTrap(layer) << "\n";
      }
      uint32_t trigger_cell = triggerGeometry_->getTriggerCellFromCell(id);
      auto itr_insert = triggercells_to_cells.emplace(trigger_cell, std::unordered_set<uint32_t>());
      itr_insert.first->second.emplace(id);
      // fill modules
      uint32_t module = triggerGeometry_->getModuleFromCell(id);
      if (module != 0) {
        itr_insert = modules_to_cells.emplace(module, std::unordered_set<uint32_t>());
        itr_insert.first->second.emplace(id);
      }
    }

    // NOSE
    if (triggerGeometry_->isWithNoseGeometry()) {
      for (const auto& id : triggerGeometry_->noseGeometry()->getValidDetIds()) {
        if (!triggerGeometry_->noseTopology().valid(id))
          continue;
        // fill trigger cells
        uint32_t trigger_cell = triggerGeometry_->getTriggerCellFromCell(id);
        auto itr_insert = triggercells_to_cells.emplace(trigger_cell, std::unordered_set<uint32_t>());
        itr_insert.first->second.emplace(id);
        // fill modules
        uint32_t module = triggerGeometry_->getModuleFromCell(id);
        if (module != 0) {
          itr_insert = modules_to_cells.emplace(module, std::unordered_set<uint32_t>());
          itr_insert.first->second.emplace(id);
        }
      }
    }

    if (module_errors.size() > 0) {
      throw cms::Exception("BadGeometry") << "HGCalTriggerGeometry: Found  module mapping problems. Check the produced "
                                             "tree to see the list of problematic wafers";
    }

    edm::LogPrint("TriggerCellCheck") << "Checking cell -> trigger cell -> cell consistency";
    // Loop over trigger cells
    for (const auto& triggercell_cells : triggercells_to_cells) {
      DetId id(triggercell_cells.first);

      // fill modules
      uint32_t module = triggerGeometry_->getModuleFromTriggerCell(id);
      if (module != 0) {
        auto itr_insert = modules_to_triggercells.emplace(module, std::unordered_set<uint32_t>());
        itr_insert.first->second.emplace(id);
      }

      // Check consistency of cells included in trigger cell
      HGCalTriggerGeometryBase::geom_set cells_geom = triggerGeometry_->getCellsFromTriggerCell(id);
      const auto& cells = triggercell_cells.second;
      for (auto cell : cells) {
        if (cells_geom.find(cell) == cells_geom.end()) {
          if (id.det() == DetId::HGCalHSc) {
            edm::LogProblem("BadTriggerCell")
                << "Error: \n Cell " << cell << "(" << HGCScintillatorDetId(cell)
                << ")\n has not been found in \n trigger cell " << HGCScintillatorDetId(id);
            std::stringstream output;
            output << " Available cells are:\n";
            for (auto cell_geom : cells_geom)
              output << "     " << HGCScintillatorDetId(cell_geom) << "\n";
            edm::LogProblem("BadTriggerCell") << output.str();
          } else if (HFNoseTriggerDetId(id).subdet() == HGCalTriggerSubdetector::HFNoseTrigger) {
            edm::LogProblem("BadTriggerCell")
                << "Error: \n Cell " << cell << "(" << HFNoseDetId(cell) << ")\n has not been found in \n trigger cell "
                << HFNoseTriggerDetId(triggercell_cells.first);
            std::stringstream output;
            output << " Available cells are:\n";
            for (auto cell_geom : cells_geom)
              output << "     " << HFNoseDetId(cell_geom) << "\n";
            edm::LogProblem("BadTriggerCell") << output.str();
          } else if (HGCalTriggerDetId(id).subdet() == HGCalTriggerSubdetector::HGCalEETrigger ||
                     HGCalTriggerDetId(id).subdet() == HGCalTriggerSubdetector::HGCalHSiTrigger) {
            edm::LogProblem("BadTriggerCell")
                << "Error: \n Cell " << cell << "(" << HGCSiliconDetId(cell)
                << ")\n has not been found in \n trigger cell " << HGCalTriggerDetId(triggercell_cells.first);
            std::stringstream output;
            output << " Available cells are:\n";
            for (auto cell_geom : cells_geom)
              output << "     " << HGCSiliconDetId(cell_geom) << "\n";
            edm::LogProblem("BadTriggerCell") << output.str();
          } else {
            edm::LogProblem("BadTriggerCell")
                << "Unknown detector type " << id.det() << " " << id.subdetId() << " " << id.rawId() << "\n";
            edm::LogProblem("BadTriggerCell") << " cell " << std::hex << cell << std::dec << " "
                                              << "\n";
            edm::LogProblem("BadTriggerCell")
                << "Cell ID " << HGCSiliconDetId(cell) << " or " << HFNoseDetId(cell) << "\n";
          }
          throw cms::Exception("BadGeometry")
              << "HGCalTriggerGeometry: Found inconsistency in cell <-> trigger cell mapping";
        }
      }
    }
    edm::LogPrint("ModuleCheck") << "Checking trigger cell -> module -> trigger cell consistency";
    // Loop over modules
    for (const auto& module_triggercells : modules_to_triggercells) {
      HGCalTriggerModuleDetId id(module_triggercells.first);
      // Check consistency of trigger cells included in module
      HGCalTriggerGeometryBase::geom_set triggercells_geom = triggerGeometry_->getTriggerCellsFromModule(id);
      const auto& triggercells = module_triggercells.second;
      for (auto cell : triggercells) {
        if (triggercells_geom.find(cell) == triggercells_geom.end()) {
          if (id.triggerSubdetId() == HGCalTriggerSubdetector::HGCalHScTrigger) {
            HGCScintillatorDetId cellid(cell);
            edm::LogProblem("BadModule") << "Error: \n Trigger cell " << cell << "(" << cellid
                                         << ")\n has not been found in \n module " << HGCalTriggerModuleDetId(id);
            std::stringstream output;
            output << " Available trigger cells are:\n";
            for (auto cell_geom : triggercells_geom) {
              output << "     " << HGCScintillatorDetId(cell_geom) << "\n";
            }
            edm::LogProblem("BadModule") << output.str();
            throw cms::Exception("BadGeometry")
                << "HGCalTriggerGeometry: Found inconsistency in trigger cell <->  module mapping";
          } else if (id.triggerSubdetId() == HGCalTriggerSubdetector::HFNoseTrigger) {
            HFNoseTriggerDetId cellid(cell);
            edm::LogProblem("BadModule") << "Error : \n Trigger cell " << cell << "(" << cellid
                                         << ")\n has not been found in \n module " << HGCalTriggerModuleDetId(id);
            std::stringstream output;
            output << " Available trigger cells are:\n";
            for (auto cell_geom : triggercells_geom) {
              output << "     " << HFNoseTriggerDetId(cell_geom) << "\n";
            }
            edm::LogProblem("BadModule") << output.str();
            throw cms::Exception("BadGeometry")
                << "HGCalTriggerGeometry: Found inconsistency in trigger cell <->  module mapping";
          } else {
            HGCalTriggerDetId cellid(cell);
            edm::LogProblem("BadModule") << "Error : \n Trigger cell " << cell << "(" << cellid
                                         << ")\n has not been found in \n module " << HGCalTriggerModuleDetId(id);
            std::stringstream output;
            output << " Available trigger cells are:\n";
            for (auto cell_geom : triggercells_geom) {
              output << "     " << HGCalTriggerDetId(cell_geom) << "\n";
            }
            edm::LogProblem("BadModule") << output.str();
            throw cms::Exception("BadGeometry")
                << "HGCalTriggerGeometry: Found inconsistency in trigger cell <->  module mapping";
          }
        }
      }
    }
    edm::LogPrint("ModuleCheck") << "Checking cell -> module -> cell consistency";
    for (const auto& module_cells : modules_to_cells) {
      HGCalTriggerModuleDetId id(module_cells.first);
      // Check consistency of cells included in module
      HGCalTriggerGeometryBase::geom_set cells_geom = triggerGeometry_->getCellsFromModule(id);
      const auto& cells = module_cells.second;
      for (auto cell : cells) {
        if (cells_geom.find(cell) == cells_geom.end()) {
          if (id.triggerSubdetId() == HGCalTriggerSubdetector::HGCalHScTrigger) {
            edm::LogProblem("BadModule") << "Error: \n Cell " << cell << "(" << HGCScintillatorDetId(cell)
                                         << ")\n has not been found in \n module " << HGCalTriggerModuleDetId(id);
          } else if (id.triggerSubdetId() == HGCalTriggerSubdetector::HFNoseTrigger) {
            edm::LogProblem("BadModule") << "Error: \n Cell " << cell << "(" << HFNoseDetId(cell)
                                         << ")\n has not been found in \n module " << HGCalTriggerModuleDetId(id);
          } else {
            edm::LogProblem("BadModule") << "Error: \n Cell " << cell << "(" << HGCSiliconDetId(cell)
                                         << ")\n has not been found in \n module " << HGCalTriggerModuleDetId(id);
          }
          std::stringstream output;
          output << " Available cells are:\n";
          for (auto cell_geom : cells_geom) {
            output << cell_geom << " ";
          }
          edm::LogProblem("BadModule") << output.str();
          throw cms::Exception("BadGeometry") << "HGCalTriggerGeometry: Found inconsistency in cell <-> module mapping";
        }
      }
    }

    // Filling Stage 1 FPGA -> modules

    edm::LogPrint("ModuleCheck") << "Checking module -> stage-1 -> module consistency";
    trigger_map_set stage1_to_modules;
    for (const auto& module_tc : modules_to_triggercells) {
      HGCalTriggerModuleDetId id(module_tc.first);
      HGCalTriggerGeometryBase::geom_set lpgbts = triggerGeometry_->getLpgbtsFromModule(id);
      if (lpgbts.size() == 0)
        continue;  //Module is not connected to an lpGBT and therefore not to a Stage 1 FPGA
      uint32_t stage1 = 0;
      for (const auto& lpgbt : lpgbts) {
        uint32_t stage1_tmp = triggerGeometry_->getStage1FpgaFromLpgbt(lpgbt);
        if (stage1 != 0 && stage1_tmp != stage1) {
          throw cms::Exception("BadGeometry") << "HGCalTriggerGeometry: Module " << HGCalTriggerModuleDetId(id)
                                              << " is split is split into more than one Stage-1 FPGA";
        }
        stage1 = stage1_tmp;
      }
      auto itr_insert = stage1_to_modules.emplace(stage1, std::unordered_set<uint32_t>());
      itr_insert.first->second.emplace(id);
    }
    // checking S1 -> module consistency

    for (const auto& stage1_modules : stage1_to_modules) {
      HGCalTriggerBackendDetId stage1(stage1_modules.first);
      HGCalTriggerGeometryBase::geom_set modules_geom;
      // Check consistency of modules going to Stage-1 FPGA
      HGCalTriggerGeometryBase::geom_set lpgbts = triggerGeometry_->getLpgbtsFromStage1Fpga(stage1);
      for (const auto& lpgbt : lpgbts) {
        HGCalTriggerGeometryBase::geom_set modules = triggerGeometry_->getModulesFromLpgbt(lpgbt);
        modules_geom.insert(modules.begin(), modules.end());
      }
      const auto& modules = stage1_modules.second;
      for (auto module : modules) {
        if (modules_geom.find(module) == modules_geom.end()) {
          edm::LogProblem("BadStage1") << "Error: \n Module " << module << "(" << HGCalTriggerModuleDetId(module)
                                       << ")\n has not been found in \n stage-1 " << HGCalTriggerBackendDetId(stage1);
          std::stringstream output;
          output << "   Available modules are:\n";
          for (auto module_geom : modules_geom) {
            output << module_geom << " ";
          }
          output << "   Connected lpgbts are:\n";
          for (auto lpgbt : lpgbts) {
            output << lpgbt << " ";
          }
          edm::LogProblem("BadStage1") << output.str();
          throw cms::Exception("BadGeometry")
              << "HGCalTriggerGeometry: Found inconsistency in Stage1 <-> module mapping";
        }
      }
    }

    // Filling Stage 2 FPGA -> Stage 1 FPGA

    edm::LogPrint("ModuleCheck") << "Checking Stage 1 -> Stage 2 -> Stage 1 consistency";
    trigger_map_set stage2_to_stage1;
    for (const auto& stage1 : stage1_to_modules) {
      HGCalTriggerBackendDetId id(stage1.first);
      HGCalTriggerGeometryBase::geom_set stage2FPGAs = triggerGeometry_->getStage2FpgasFromStage1Fpga(id);
      for (const auto& stage2 : stage2FPGAs) {
        auto itr_insert = stage2_to_stage1.emplace(stage2, std::unordered_set<uint32_t>());
        itr_insert.first->second.emplace(id);
      }
    }
    // checking S1 -> S2 consistency

    for (const auto& stage2_modules : stage2_to_stage1) {
      HGCalTriggerBackendDetId stage2(stage2_modules.first);

      // Check consistency of Stage-1 FPGA going to Stage 2 FPGA
      HGCalTriggerGeometryBase::geom_set stage1FPGAs = triggerGeometry_->getStage1FpgasFromStage2Fpga(stage2);

      const auto& stage1fpgas = stage2_modules.second;

      for (auto stage1fpga : stage1fpgas) {
        if (stage1FPGAs.find(stage1fpga) == stage1FPGAs.end()) {
          edm::LogProblem("BadStage2") << "Error: \n Stage-1 FPGA " << stage1fpga << "("
                                       << HGCalTriggerBackendDetId(stage1fpga)
                                       << ")\n has not been found in \n Stage-2 " << HGCalTriggerBackendDetId(stage2);
          std::stringstream output;
          output << "\n   Available Stage-1 FPGAs are:\n";
          for (auto stage1FPGA : stage1FPGAs) {
            output << HGCalTriggerBackendDetId(stage1FPGA) << "\n";
          }
          edm::LogProblem("BadStage2") << output.str();
          throw cms::Exception("BadGeometry")
              << "HGCalTriggerGeometry: Found inconsistency in Stage2 <-> Stage1 mapping";
        }
      }
    }

  } catch (const cms::Exception& e) {
    edm::LogWarning("HGCalTriggerGeometryTester")
        << "Problem with the trigger geometry detected. Only the basic cells tree will be filled\n";
    edm::LogWarning("HGCalTriggerGeometryTester") << e.message() << "\n";
    return false;
  }
  return true;
}

void HGCalBackendStage1ParameterExtractor::fillTriggerGeometry(json& json_file) {
  trigger_map_set trigger_cells;

  // check geometry and retrieve trigger cells
  std::cout << "Checking EE geometry\n";
  for (const auto& id : triggerGeometry_->eeGeometry()->getValidDetIds()) {
    HGCSiliconDetId detid(id);
    if (!triggerGeometry_->eeTopology().valid(id))
      continue;
    // fill trigger cells
    if (!no_trigger_) {
      uint32_t trigger_cell = triggerGeometry_->getTriggerCellFromCell(id);
      auto itr_insert = trigger_cells.emplace(trigger_cell, std::unordered_set<uint32_t>());
      itr_insert.first->second.emplace(id);
    }
  }
  std::cout << "Checking HSi geometry\n";
  for (const auto& id : triggerGeometry_->hsiGeometry()->getValidDetIds()) {
    HGCSiliconDetId detid(id);
    if (!triggerGeometry_->hsiTopology().valid(id))
      continue;
    // fill trigger cells
    if (!no_trigger_) {
      uint32_t trigger_cell = triggerGeometry_->getTriggerCellFromCell(id);
      auto itr_insert = trigger_cells.emplace(trigger_cell, std::unordered_set<uint32_t>());
      itr_insert.first->second.emplace(id);
    }
  }
  std::cout << "Checking HSc geometry\n";
  for (const auto& id : triggerGeometry_->hscGeometry()->getValidDetIds()) {
    if (!no_trigger_) {
      uint32_t trigger_cell = triggerGeometry_->getTriggerCellFromCell(id);
      auto itr_insert = trigger_cells.emplace(trigger_cell, std::unordered_set<uint32_t>());
      itr_insert.first->second.emplace(id);
    }
  }
  if (triggerGeometry_->isWithNoseGeometry()) {
    std::cout << "Checking NOSE geometry\n";
    for (const auto& id : triggerGeometry_->noseGeometry()->getValidDetIds()) {
      // fill trigger cells
      if (!no_trigger_) {
        uint32_t trigger_cell = triggerGeometry_->getTriggerCellFromCell(id);
        auto itr_insert = trigger_cells.emplace(trigger_cell, std::unordered_set<uint32_t>());
        itr_insert.first->second.emplace(id);
      }
    }
  }

  // if problem detected in the trigger geometry, don't produce the TC map
  if (no_trigger_)
    return;

  // loop over trigger cells
  edm::LogPrint("JSONFilling") << "Filling JSON map";
  for (const auto& triggercell_cells : trigger_cells) {
    DetId id(triggercell_cells.first);

    // get module ID and check if relevant module (sector0, zside=-1)
    uint32_t moduleId = triggerGeometry_->getModuleFromTriggerCell(id);
    if (moduleId == 0)
      continue;
    HGCalTriggerModuleDetId tc_module(moduleId);
    if (!(tc_module.isHGCalModuleDetId()) || (tc_module.zside() < 0) || (tc_module.sector() != 0))
      continue;

    // check that module is connected to a stage 1 FPGA
    HGCalTriggerGeometryBase::geom_set lpgbts = triggerGeometry_->getLpgbtsFromModule(moduleId);
    if (lpgbts.size() == 0)
      continue;

    // only retrieve mapping for the tested fpga
    uint32_t fpgaId = triggerGeometry_->getStage1FpgaFromModule(moduleId);
    HGCalTriggerBackendDetId tc_fpga(fpgaId);
    if (!(tc_fpga.isStage1FPGA()) || (tc_fpga.sector() != 0) || (tc_fpga.zside() < 0))
      continue;
    int fpgaLabel = tc_fpga.label();
    if (fpgaLabel != the_fpga_)
      continue;

    // retrieve information to be saved
    int triggerCellSubdet = 0;
    int triggerCellLayer = 0;
    int triggerCellUEta = 0;
    int triggerCellVPhi = 0;
    if (id.det() == DetId::HGCalHSc) {
      HGCScintillatorDetId id_sc(triggercell_cells.first);
      triggerCellSubdet = id_sc.subdet();
      triggerCellLayer = id_sc.layer();
      triggerCellUEta = id_sc.ietaAbs();
      triggerCellVPhi = id_sc.iphi();
    } else if (HFNoseTriggerDetId(triggercell_cells.first).det() == DetId::HGCalTrigger &&
               HFNoseTriggerDetId(triggercell_cells.first).subdet() == HGCalTriggerSubdetector::HFNoseTrigger) {
      HFNoseTriggerDetId id_nose_trig(triggercell_cells.first);
      triggerCellSubdet = id_nose_trig.subdet();
      triggerCellLayer = id_nose_trig.layer();
      triggerCellUEta = id_nose_trig.triggerCellU();
      triggerCellVPhi = id_nose_trig.triggerCellV();
    } else {
      HGCalTriggerDetId id_si_trig(triggercell_cells.first);
      triggerCellSubdet = id_si_trig.subdet();
      triggerCellLayer = id_si_trig.layer();
      triggerCellUEta = id_si_trig.triggerCellU();
      triggerCellVPhi = id_si_trig.triggerCellV();
    }
    GlobalPoint position = triggerGeometry_->getTriggerCellPosition(triggercell_cells.first);
    float triggerCellX = position.x();
    float triggerCellY = position.y();
    float triggerCellZ = position.z();
    float triggerCellEta = position.eta();
    float triggerCellPhi = position.phi();
    float triggerCellRoverZ = sqrt(triggerCellX * triggerCellX + triggerCellY * triggerCellY) / triggerCellZ;

    // transform HGCalHSc TC coordinates to define a TC id in [0,47]
    uint32_t tce = triggerCellUEta;
    uint32_t tcp = triggerCellVPhi;
    if (triggerCellSubdet == 10) {  // HGCalHSc
      tcp = (triggerCellVPhi - 1) % 4;
      if (triggerCellUEta <= 3) {
        tce = triggerCellUEta;
      } else if (triggerCellUEta <= 9) {
        tce = triggerCellUEta - 4;
      } else if (triggerCellUEta <= 13) {
        tce = triggerCellUEta - 10;
      } else if (triggerCellUEta <= 17) {
        tce = triggerCellUEta - 14;
      } else {
        tce = triggerCellUEta - 18;
      }
    }

    // attribute ID to TC according to subdetector
    uint32_t tcaddress = 99;
    if (triggerCellSubdet == 10) {  //HGCalHSc(10)
      tcaddress = (int(tce) << 2) + int(tcp);
    } else {  //HGCalHSiTrigger(2) or HGCalEE(1)
      tcaddress = tc_coord_uv_.find(std::make_pair(tce, tcp))->second;
    }

    // save TC info into JSON
    std::string strModId = std::to_string(moduleId);
    std::string strTCAddr = std::to_string(tcaddress);
    json_file["TriggerCellMap"]["module_hash"][strModId]["tc_id"][strTCAddr]["subdet"] = triggerCellSubdet;
    json_file["TriggerCellMap"]["module_hash"][strModId]["tc_id"][strTCAddr]["layer"] = triggerCellLayer;
    json_file["TriggerCellMap"]["module_hash"][strModId]["tc_id"][strTCAddr]["ueta"] = tce;
    json_file["TriggerCellMap"]["module_hash"][strModId]["tc_id"][strTCAddr]["vphi"] = tcp;
    json_file["TriggerCellMap"]["module_hash"][strModId]["tc_id"][strTCAddr]["x"] = triggerCellX;
    json_file["TriggerCellMap"]["module_hash"][strModId]["tc_id"][strTCAddr]["y"] = triggerCellY;
    json_file["TriggerCellMap"]["module_hash"][strModId]["tc_id"][strTCAddr]["z"] = triggerCellZ;
    json_file["TriggerCellMap"]["module_hash"][strModId]["tc_id"][strTCAddr]["roz"] = triggerCellRoverZ;
    json_file["TriggerCellMap"]["module_hash"][strModId]["tc_id"][strTCAddr]["eta"] = triggerCellEta;
    json_file["TriggerCellMap"]["module_hash"][strModId]["tc_id"][strTCAddr]["phi"] = triggerCellPhi;
  }
}

void HGCalBackendStage1ParameterExtractor::analyze(const edm::Event& e, const edm::EventSetup& es) {}

// define this as a plug-in
DEFINE_FWK_MODULE(HGCalBackendStage1ParameterExtractor);
