#include "Geometry/HGCalCommonData/interface/HGCalGeomRotation.h"
#include "DataFormats/ForwardDetId/interface/ForwardSubdetector.h"
#include "DataFormats/ForwardDetId/interface/HFNoseTriggerDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCScintillatorDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCSiliconDetIdToROC.h"
#include "DataFormats/ForwardDetId/interface/HGCalDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCalTriggerDetId.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "L1Trigger/L1THGCal/interface/geometries/HGCalTriggerGeometryV9.h"
#include "DataFormats/ForwardDetId/interface/HGCalTriggerModuleDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCalTriggerBackendDetId.h"
#include "DataFormats/ForwardDetId/interface/HFNoseDetIdToModule.h"

#include <fstream>
#include <vector>

#include <nlohmann/json.hpp>
using json = nlohmann::json;

class HGCalTriggerGeometryV9Imp3 : public HGCalTriggerGeometryV9 {
public:
  HGCalTriggerGeometryV9Imp3(const edm::ParameterSet& conf);
  void reset() final;

  geom_set getStage1LinksFromStage2Fpga(const unsigned) const final;
  unsigned getStage2FpgaFromStage1Link(const unsigned) const final;

private:
  void fillMaps() final;
  edm::FileInPath jsonMappingFile_;

  std::unordered_map<unsigned, bool> stage1links_samesector_;
};

HGCalTriggerGeometryV9Imp3::HGCalTriggerGeometryV9Imp3(const edm::ParameterSet& conf)
    : HGCalTriggerGeometryV9(conf), jsonMappingFile_(conf.getParameter<edm::FileInPath>("JsonMappingFile")) {}

void HGCalTriggerGeometryV9Imp3::reset() {
  HGCalTriggerGeometryV9::reset();
  stage1links_samesector_.clear();
}

HGCalTriggerGeometryBase::geom_set HGCalTriggerGeometryV9Imp3::getStage1LinksFromStage2Fpga(
    const unsigned stage2_id) const {
  geom_set stage1link_ids;
  HGCalTriggerBackendDetId id(stage2_id);
  auto stage2_itrs = stage2_to_stage1links_.equal_range(id.label());
  for (auto stage2_itr = stage2_itrs.first; stage2_itr != stage2_itrs.second; stage2_itr++) {
    unsigned label = stage2_itr->second;
    if (stage1links_samesector_.at(label) == true) {  //link and stage2 FPGA are the same sector
      stage1link_ids.emplace(
          HGCalTriggerBackendDetId(id.zside(), HGCalTriggerBackendDetId::BackendType::Stage1Link, id.sector(), label));
    } else {  //link is from the next sector (anti-clockwise)
      stage1link_ids.emplace(HGCalTriggerBackendDetId(
          id.zside(), HGCalTriggerBackendDetId::BackendType::Stage1Link, getNextSector(id.sector()), label));
    }
  }

  return stage1link_ids;
}

unsigned HGCalTriggerGeometryV9Imp3::getStage2FpgaFromStage1Link(const unsigned link_id) const {
  HGCalTriggerBackendDetId id(link_id);
  bool same_sector = stage1link_to_stage2_.at(id.label());
  unsigned sector = id.sector();

  if (!same_sector) {
    sector = getPreviousSector(sector);
  }

  return HGCalTriggerBackendDetId(id.zside(), HGCalTriggerBackendDetId::BackendType::Stage2FPGA, sector, 0);
}

void HGCalTriggerGeometryV9Imp3::fillMaps() {
  // read json mapping file
  json mapping_config;
  std::ifstream json_input_file(jsonMappingFile_.fullPath());
  if (!json_input_file.is_open()) {
    throw cms::Exception("MissingDataFile") << "Cannot open HGCalTriggerGeometry L1TMapping file\n";
  }
  json_input_file >> mapping_config;

  try {
    //Stage 2 to Stage 1 links mapping
    for (unsigned stage2_id = 0; stage2_id < mapping_config.at("Stage2").size(); stage2_id++) {
      for (unsigned link_id = 0; link_id < mapping_config.at("Stage2").at(stage2_id).at("Stage1Links").size();
           link_id++) {
        stage2_to_stage1links_.emplace(stage2_id, link_id);
        stage1links_samesector_.emplace(
            link_id, mapping_config.at("Stage2").at(stage2_id).at("Stage1Links").at(link_id).at("SameSector"));
      }
    }
  } catch (const json::exception& e) {
    edm::LogError("HGCalTriggerGeometryV9Imp3")
        << "The mapping input json file does not have the expected structure for the Stage2 block";
  }

  try {
    for (unsigned link_id = 0; link_id < mapping_config.at("Stage1Links").size(); link_id++) {
      //Stage 1 links to Stage 1 FPGAs mapping
      stage1link_to_stage1_.emplace(link_id, mapping_config.at("Stage1Links").at(link_id).at("Stage1"));

      //Stage 1 links to Stage 2 mapping
      stage1link_to_stage2_.emplace(link_id, mapping_config.at("Stage1Links").at(link_id).at("Stage2SameSector"));
    }
  } catch (const json::exception& e) {
    edm::LogError("HGCalTriggerGeometryV9Imp3")
        << "The mapping input json file does not have the expected structure for the Stage1Links block";
  }

  try {
    for (unsigned stage1_id = 0; stage1_id < mapping_config.at("Stage1").size(); stage1_id++) {
      //Stage 1 to Stage 1 links mapping
      for (auto& link_id : mapping_config.at("Stage1").at(stage1_id).at("Stage1Links")) {
        stage1_to_stage1links_.emplace(stage1_id, link_id);
      }

      //Stage 1 to lpgbt mapping
      std::vector<unsigned> lpgbt_id_vec;
      for (auto& lpgbt_id : mapping_config.at("Stage1").at(stage1_id).at("lpgbts")) {
        lpgbt_id_vec.push_back(lpgbt_id);
      }
      stage1_to_lpgbts_.emplace(stage1_id, lpgbt_id_vec);
    }

  } catch (const json::exception& e) {
    edm::LogError("HGCalTriggerGeometryV9Imp3")
        << "The mapping input json file does not have the expected structure for the Stage1 block";
  }

  try {
    for (unsigned lpgbt_id = 0; lpgbt_id < mapping_config.at("lpgbt").size(); lpgbt_id++) {
      //lpgbt to Stage 1 mapping
      unsigned stage1_id = mapping_config.at("lpgbt").at(lpgbt_id).at("Stage1");
      lpgbt_to_stage1_.emplace(lpgbt_id, stage1_id);

      //lpgbt to module mapping
      for (auto& modules : mapping_config.at("lpgbt").at(lpgbt_id).at("Modules")) {
        unsigned layer = modules.at("layer");
        int subdetId = 0;
        bool isSilicon = modules.at("isSilicon");
        layerWithoutOffsetAndSubdetId(layer, subdetId, isSilicon);
        unsigned packed_value = packLayerSubdetWaferId(layer, subdetId, modules.at("u"), modules.at("v"));
        lpgbt_to_modules_.emplace(lpgbt_id, packed_value);

        //fill subsiduary module to stage 1 mapping
        auto result = module_to_stage1_.emplace(packed_value, stage1_id);
        if (result.second == false &&
            stage1_id != result.first->second) {  //check that the stage1_id is the same as in the existing map
          edm::LogError("HGCalTriggerGeometryV9Imp3") << "One module is connected to two separate Stage1 FPGAs";
        }
      }
    }

  } catch (const json::exception& e) {
    edm::LogError("HGCalTriggerGeometryV9Imp3")
        << "The mapping input json file does not have the expected structure for the lpGBT block";
  }

  try {
    //module to lpgbt mapping
    for (unsigned module = 0; module < mapping_config.at("Module").size(); module++) {
      unsigned num_elinks = 0;  //Sum number of e-links in each module over lpGBTs
      unsigned layer = mapping_config.at("Module").at(module).at("layer");
      unsigned moduleU = mapping_config.at("Module").at(module).at("u");
      unsigned moduleV = mapping_config.at("Module").at(module).at("v");
      bool isSilicon = mapping_config.at("Module").at(module).at("isSilicon");
      int subdetId = 0;
      layerWithoutOffsetAndSubdetId(layer, subdetId, isSilicon);

      for (auto& lpgbt : mapping_config.at("Module").at(module).at("lpgbts")) {
        module_to_lpgbts_.emplace(packLayerSubdetWaferId(layer, subdetId, moduleU, moduleV), lpgbt.at("id"));
        num_elinks += unsigned(lpgbt.at("nElinks"));
      }
      int packed_module = packLayerSubdetWaferId(layer, subdetId, moduleU, moduleV);
      links_per_module_.emplace(packed_module, num_elinks);
    }
  } catch (const json::exception& e) {
    edm::LogError("HGCalTriggerGeometryV9Imp3")
        << "The mapping input json file does not have the expected structure for the Module block";
  }

  json_input_file.close();
}

DEFINE_EDM_PLUGIN(HGCalTriggerGeometryFactory, HGCalTriggerGeometryV9Imp3, "HGCalTriggerGeometryV9Imp3");
