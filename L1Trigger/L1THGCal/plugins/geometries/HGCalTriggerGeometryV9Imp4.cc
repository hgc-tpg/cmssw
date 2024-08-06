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

#include "L1Trigger/L1THGCal/interface/mappingTools/HgcConfigReader.hpp"
#include "L1Trigger/L1THGCal/interface/mappingTools/GUID.hpp"

class HGCalTriggerGeometryV9Imp4 : public HGCalTriggerGeometryV9 {
public:
  HGCalTriggerGeometryV9Imp4(const edm::ParameterSet& conf);
  void reset() final;

  geom_set getStage1LinksFromStage2Fpga(const unsigned) const final;
  unsigned getStage2FpgaFromStage1Link(const unsigned) const final;

private:
  void fillMaps() final;
  edm::FileInPath xmlMappingFile_;

  std::unordered_map<unsigned, bool> stage1links_whichsector_;
};

HGCalTriggerGeometryV9Imp4::HGCalTriggerGeometryV9Imp4(const edm::ParameterSet& conf)
    : HGCalTriggerGeometryV9(conf), xmlMappingFile_(conf.getParameter<edm::FileInPath>("xmlMappingFile")) {}

void HGCalTriggerGeometryV9Imp4::reset() {
  HGCalTriggerGeometryV9::reset();
  stage1links_whichsector_.clear();
}

HGCalTriggerGeometryBase::geom_set HGCalTriggerGeometryV9Imp4::getStage1LinksFromStage2Fpga(
    const unsigned stage2_id) const {
  geom_set stage1link_ids;
  HGCalTriggerBackendDetId id(stage2_id);

  auto stage2_itrs = stage2_to_stage1links_.equal_range(id.label());

  for (auto stage2_itr = stage2_itrs.first; stage2_itr != stage2_itrs.second; stage2_itr++) {
    unsigned label = stage2_itr->second;
    const int whichSector = stage1links_whichsector_.at(label);
    if (whichSector == 0) {  //link and stage2 FPGA are the same sector
      stage1link_ids.emplace(
          HGCalTriggerBackendDetId(id.zside(), HGCalTriggerBackendDetId::BackendType::Stage1Link, id.sector(), label));
    } else if (whichSector == 1) {  //link is from the next sector (anti-clockwise)
      stage1link_ids.emplace(HGCalTriggerBackendDetId(
          id.zside(), HGCalTriggerBackendDetId::BackendType::Stage1Link, getNextSector(id.sector()), label));
    } else {
      stage1link_ids.emplace(HGCalTriggerBackendDetId(
          id.zside(), HGCalTriggerBackendDetId::BackendType::Stage1Link, getPreviousSector(id.sector()), label));
    }
  }

  return stage1link_ids;
}

unsigned HGCalTriggerGeometryV9Imp4::getStage2FpgaFromStage1Link(const unsigned link_id) const {
  HGCalTriggerBackendDetId id(link_id);
  int whichSector = stage1link_to_stage2_.at(id.label());
  unsigned sector = id.sector();

  if (whichSector == 1) {
    sector = getPreviousSector(sector);
  } else if (whichSector == -1) {
    sector = getNextSector(sector);
  }

  return HGCalTriggerBackendDetId(id.zside(), HGCalTriggerBackendDetId::BackendType::Stage2FPGA, sector, 0);
}

void HGCalTriggerGeometryV9Imp4::fillMaps() {
  l1thgcmapping::Maps lMaps;
  // l1thgcmapping::OpenBackendMapping( "scenarios/CombinedTD.60.MixedTypes.NoSplit/BackendMapping.xml" , lMaps );

  // Extract base path of mapping files from xmlMappingFile_
  std::string searchFor = "mapping/";
  std::string basePath = "";
  unsigned pos = xmlMappingFile_.fullPath().rfind(searchFor);
  if (pos != std::string::npos) {
    pos += searchFor.length();
    basePath = xmlMappingFile_.fullPath().substr(0, pos);
  }

  // Parse xml files and fill maps
  l1thgcmapping::OpenBackendMapping(xmlMappingFile_.fullPath(), basePath, lMaps);

  unsigned int nBELinksToS2In = 0;

  // Stage 2 board to Stage 1 links mapping
  for (const auto& [belink, s2in] : lMaps.belink_to_s2in) {
    if (l1thgcmapping::get_endcap(s2in) == 0 && l1thgcmapping::get_sector(s2in) == 0 &&
        l1thgcmapping::get_s2tm(s2in) == 0) {
      // Get Stage 1 sector
      const auto s1out = lMaps.belink_to_s1out.at(belink);
      const auto s1outSector = l1thgcmapping::get_sector(s1out);
      // Is the S1 sector the same, the next sector (clockwise) or previous sector (anticlockwise)
      int whichSector = (s1outSector == 2) ? -1 : s1outSector;

      // Equivalent to existing mapping structures
      const auto link_id = l1thgcmapping::get_io_index(s2in);
      stage2_to_stage1links_.emplace(0, link_id);
      stage1links_whichsector_.emplace(link_id, whichSector);
      ++nBELinksToS2In;
    }
  }

  // Stage 1 links to Stage 1 board mapping
  // Loop links between S1 and S2.  Select those in 1 endcap, s1 sector, and s2 tm period
  // Number the links from 0 to N (84 in one scenario being considered), and map to s1 board (stage1link_to_stage1_)
  // And vice versa (stage1_to_stage1links_)
  // and also map to which (relative) s2 sector the links is going to (-1, 0, 1) (stage1link_to_stage2_)
  unsigned int iLink = 0;
  for (const auto& [belink, s1out] : lMaps.belink_to_s1out) {
    const auto s2in = lMaps.belink_to_s2in.at(belink);
    const auto s2tm = l1thgcmapping::get_s2tm(s2in);
    const auto endcap = l1thgcmapping::get_endcap(s1out);
    const auto s1Sector = l1thgcmapping::get_sector(s1out);

    if (endcap == 0 && s2tm == 0 && s1Sector == 0) {
      const auto s1Index = l1thgcmapping::get_s1index(belink);
      const auto link_id = l1thgcmapping::get_io_index(s2in);
      stage1link_to_stage1_.emplace(link_id, s1Index);
      stage1_to_stage1links_.emplace(s1Index, link_id);

      // Is the S2 sector the same, the next sector (clockwise) or previous sector (anticlockwise)
      const auto s2Sector = l1thgcmapping::get_sector(s2in);
      int whichSector = (s2Sector == 2) ? -1 : s2Sector;
      stage1link_to_stage2_.emplace(link_id, whichSector);
      ++iLink;
    }
  }

  // Stage 1 boards to input lpgbt mapping
  // Loop S1 boards, which correspond to those within one 120 S1 sector
  for (const auto& [s1, region] : lMaps.stage1_to_region) {
    const auto s1Index = l1thgcmapping::get_s1index(s1);
    std::vector<unsigned> lpgbts;

    // Get motherboards within this region
    auto motherboardRange = lMaps.region_to_motherboard.equal_range(region);
    for (auto it = motherboardRange.first; it != motherboardRange.second; ++it) {
      const auto lpgbtIds = lMaps.motherboard_to_trigLPGBTs.at(it->second);
      lpgbts.reserve(lpgbts.size() + lpgbtIds.size());
      lpgbts.insert(lpgbts.end(), lpgbtIds.begin(), lpgbtIds.end());

      for (const auto lpgbt : lpgbtIds) {
        lpgbt_to_stage1_.emplace(lpgbt, s1Index);
      }
    }

    if (stage1_to_lpgbts_.find(s1Index) != stage1_to_lpgbts_.end()) {
      stage1_to_lpgbts_[s1Index].insert(stage1_to_lpgbts_[s1Index].end(), lpgbts.begin(), lpgbts.end());
    } else {
      stage1_to_lpgbts_[s1Index] = lpgbts;
    }
  }

  // Loop motherboards->lpGBTs, map to S1 boards
  unsigned int maxMBID = 0;
  for (const auto& [motherboard, region] : lMaps.motherboard_to_region) {
    if (uint16_t(motherboard & 0xFFFF) > maxMBID)
      maxMBID = uint16_t(motherboard & 0xFFFF);

    unsigned nLPGBTs = lMaps.motherboard_to_nTrigLPGBT.at(motherboard);
    if (nLPGBTs == 0)
      continue;  // Disconnected CE-E layers for trigger
    const auto s1 = lMaps.region_to_stage1.at(region);
    const auto s1Index = l1thgcmapping::get_s1index(s1);
    const auto lpgbtIds = lMaps.motherboard_to_trigLPGBTs.at(motherboard);

    // Modules in this motherboard
    auto moduleRange = lMaps.motherboard_to_module.equal_range(motherboard);
    for (auto it = moduleRange.first; it != moduleRange.second; ++it) {
      const auto plane = l1thgcmapping::get_plane(it->second);
      const auto type = l1thgcmapping::get_object_type(it->second);

      bool isSilicon = (type == 0);
      int subdetId = 0;
      unsigned layer = plane;
      layerWithoutOffsetAndSubdetId(layer, subdetId, isSilicon);
      const auto module_uv = l1thgcmapping::get_module_uv(it->second);
      unsigned packed_value = packLayerSubdetWaferId(layer, subdetId, module_uv.first, module_uv.second);

      for (const auto lpgbt : lpgbtIds) {
        lpgbt_to_modules_.emplace(lpgbt, packed_value);
        module_to_lpgbts_.emplace(packed_value, lpgbt);
      }
      module_to_stage1_.emplace(packed_value, s1Index);
      links_per_module_.emplace(packed_value, 0);  // dummy
    }
  }
}

DEFINE_EDM_PLUGIN(HGCalTriggerGeometryFactory, HGCalTriggerGeometryV9Imp4, "HGCalTriggerGeometryV9Imp4");
