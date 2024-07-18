#include "Geometry/HGCalCommonData/interface/HGCalGeomRotation.h"
#include "DataFormats/ForwardDetId/interface/ForwardSubdetector.h"
#include "DataFormats/ForwardDetId/interface/HFNoseTriggerDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCScintillatorDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCSiliconDetIdToROC.h"
#include "DataFormats/ForwardDetId/interface/HGCalDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCalTriggerDetId.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "L1Trigger/L1THGCal/interface/HGCalTriggerGeometryBase.h"
#include "DataFormats/ForwardDetId/interface/HGCalTriggerModuleDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCalTriggerBackendDetId.h"
#include "DataFormats/ForwardDetId/interface/HFNoseDetIdToModule.h"

#include <fstream>
#include <vector>

#include "L1Trigger/L1THGCal/interface/mappingTools/HgcConfigReader.hpp"
#include "L1Trigger/L1THGCal/interface/mappingTools/GUID.hpp"

class HGCalTriggerGeometryV9Imp4 : public HGCalTriggerGeometryBase {
public:
  HGCalTriggerGeometryV9Imp4(const edm::ParameterSet& conf);

  void initialize(const HGCalGeometry*, const HGCalGeometry*, const HGCalGeometry*) final;
  void initialize(const HGCalGeometry*, const HGCalGeometry*, const HGCalGeometry*, const HGCalGeometry*) final;
  void reset() final;

  unsigned getTriggerCellFromCell(const unsigned) const final;
  unsigned getModuleFromCell(const unsigned) const final;
  unsigned getModuleFromTriggerCell(const unsigned) const final;

  geom_set getCellsFromTriggerCell(const unsigned) const final;
  geom_set getCellsFromModule(const unsigned) const final;
  geom_set getTriggerCellsFromModule(const unsigned) const final;

  geom_ordered_set getOrderedCellsFromModule(const unsigned) const final;
  geom_ordered_set getOrderedTriggerCellsFromModule(const unsigned) const final;

  geom_set getNeighborsFromTriggerCell(const unsigned) const final;

  geom_set getStage1FpgasFromStage2Fpga(const unsigned) const final;
  geom_set getStage2FpgasFromStage1Fpga(const unsigned) const final;

  geom_set getStage1LinksFromStage2Fpga(const unsigned) const final;
  unsigned getStage1FpgaFromStage1Link(const unsigned) const final;
  unsigned getStage2FpgaFromStage1Link(const unsigned) const final;
  geom_set getStage1LinksFromStage1Fpga(const unsigned) const final;
  std::vector<unsigned> getLpgbtsFromStage1Fpga(const unsigned) const final;
  unsigned getStage1FpgaFromLpgbt(const unsigned) const final;
  geom_set getModulesFromLpgbt(const unsigned) const final;
  geom_set getLpgbtsFromModule(const unsigned) const final;
  unsigned getStage1FpgaFromModule(const unsigned module_id) const final;

  unsigned getLinksInModule(const unsigned module_id) const final;
  unsigned getModuleSize(const unsigned module_id) const final;

  GlobalPoint getTriggerCellPosition(const unsigned) const final;
  GlobalPoint getModulePosition(const unsigned) const final;

  bool validCell(const unsigned) const final;
  bool validTriggerCell(const unsigned) const final;
  bool disconnectedModule(const unsigned) const final;
  unsigned lastTriggerLayer() const final { return last_trigger_layer_; }
  unsigned triggerLayer(const unsigned) const final;
  const std::vector<unsigned>& triggerLayers() const final { return trigger_layers_; }

private:
  // HSc trigger cell grouping
  unsigned hSc_triggercell_size_ = 2;
  static constexpr unsigned hSc_num_panels_per_sector_ = 12;
  static constexpr unsigned hSc_tcs_per_module_phi_ = 4;
  static constexpr unsigned hSc_front_layers_split_ = 12;
  static constexpr unsigned hSc_back_layers_split_ = 8;
  static constexpr unsigned hSc_layer_for_split_ = 40;
  static constexpr int hSc_tc_layer0_min_ = 24;
  static constexpr int ntc_per_wafer_ = 48;
  static constexpr int nSectors_ = 3;

  edm::FileInPath xmlMappingFile_;

  // rotation class
  HGCalGeomRotation geom_rotation_120_ = {HGCalGeomRotation::SectorType::Sector120Degrees};

  // module related maps
  std::unordered_map<unsigned, unsigned> links_per_module_;

  std::unordered_multimap<unsigned, unsigned> stage2_to_stage1links_;
  std::unordered_map<unsigned, bool> stage1links_samesector_;
  std::unordered_map<unsigned, int> stage1links_whichsector_;
  std::unordered_map<unsigned, int> stage1link_to_stage2_;
  std::unordered_map<unsigned, unsigned> stage1link_to_stage1_;
  std::unordered_multimap<unsigned, unsigned> stage1_to_stage1links_;
  std::unordered_map<unsigned, std::vector<unsigned>> stage1_to_lpgbts_;
  std::unordered_map<unsigned, unsigned> stage1_to_nLpgbts_;
  std::unordered_map<unsigned, unsigned> lpgbt_to_stage1_;
  std::unordered_multimap<unsigned, unsigned> lpgbt_to_modules_;
  std::unordered_multimap<unsigned, unsigned> module_to_lpgbts_;
  std::unordered_map<unsigned, unsigned> module_to_stage1_;

  // Disconnected modules and layers
  std::unordered_set<unsigned> disconnected_layers_;
  std::vector<unsigned> trigger_layers_;
  std::vector<unsigned> trigger_nose_layers_;
  unsigned last_trigger_layer_ = 0;

  // layer offsets
  unsigned heOffset_ = 0;
  unsigned noseLayers_ = 0;
  unsigned totalLayers_ = 0;

  void fillMaps();

  bool validCellId(unsigned det, unsigned cell_id) const;
  bool validTriggerCellFromCells(const unsigned) const;

  int detIdWaferType(unsigned det, unsigned layer, short waferU, short waferV) const;
  void layerWithoutOffsetAndSubdetId(unsigned& layer, int& subdetId, bool isSilicon) const;
  unsigned packLayerSubdetWaferId(unsigned layer, int subdet, int waferU, int waferV) const;
  void unpackLayerSubdetWaferId(unsigned wafer, unsigned& layer, int& subdet, int& waferU, int& waferV) const;
  HGCalGeomRotation::WaferCentring getWaferCentring(unsigned layer, int subdet) const;
  void etaphiMappingFromSector0(int& ieta, int& iphi, unsigned sector) const;
  unsigned tcEtaphiMappingToSector0(int& tc_ieta, int& tc_iphi) const;
  void getScintillatoriEtaiPhi(int& ieta, int& iphi, int tc_eta, int tc_phi, unsigned layer) const;
  unsigned layerWithOffset(unsigned) const;
  unsigned getNextSector(const unsigned sector) const;
  unsigned getPreviousSector(const unsigned sector) const;
};

HGCalTriggerGeometryV9Imp4::HGCalTriggerGeometryV9Imp4(const edm::ParameterSet& conf)
    : HGCalTriggerGeometryBase(conf),
      hSc_triggercell_size_(conf.getParameter<unsigned>("ScintillatorTriggerCellSize")),
      xmlMappingFile_(conf.getParameter<edm::FileInPath>("xmlMappingFile")) {
  std::vector<unsigned> tmp_vector = conf.getParameter<std::vector<unsigned>>("DisconnectedLayers");
  std::move(tmp_vector.begin(), tmp_vector.end(), std::inserter(disconnected_layers_, disconnected_layers_.end()));
}

void HGCalTriggerGeometryV9Imp4::reset() {
  stage2_to_stage1links_.clear();
  stage1links_samesector_.clear();
  stage1links_whichsector_.clear();
  stage1link_to_stage2_.clear();
  stage1link_to_stage1_.clear();
  stage1_to_stage1links_.clear();
  stage1_to_lpgbts_.clear();
  stage1_to_nLpgbts_.clear();
  lpgbt_to_stage1_.clear();
  lpgbt_to_modules_.clear();
  module_to_lpgbts_.clear();
  module_to_stage1_.clear();
}

void HGCalTriggerGeometryV9Imp4::initialize(const HGCalGeometry* hgc_ee_geometry,
                                            const HGCalGeometry* hgc_hsi_geometry,
                                            const HGCalGeometry* hgc_hsc_geometry) {
  setEEGeometry(hgc_ee_geometry);
  setHSiGeometry(hgc_hsi_geometry);
  setHScGeometry(hgc_hsc_geometry);
  heOffset_ = eeTopology().dddConstants().layers(true);
  totalLayers_ = heOffset_ + hsiTopology().dddConstants().layers(true);
  trigger_layers_.resize(totalLayers_ + 1);
  trigger_layers_[0] = 0;  // layer number 0 doesn't exist
  unsigned trigger_layer = 1;
  for (unsigned layer = 1; layer < trigger_layers_.size(); layer++) {
    if (disconnected_layers_.find(layer) == disconnected_layers_.end()) {
      // Increase trigger layer number if the layer is not disconnected
      trigger_layers_[layer] = trigger_layer;
      trigger_layer++;
    } else {
      trigger_layers_[layer] = 0;
    }
  }
  last_trigger_layer_ = trigger_layer - 1;
  fillMaps();
}

void HGCalTriggerGeometryV9Imp4::initialize(const HGCalGeometry* hgc_ee_geometry,
                                            const HGCalGeometry* hgc_hsi_geometry,
                                            const HGCalGeometry* hgc_hsc_geometry,
                                            const HGCalGeometry* hgc_nose_geometry) {
  setEEGeometry(hgc_ee_geometry);
  setHSiGeometry(hgc_hsi_geometry);
  setHScGeometry(hgc_hsc_geometry);
  setNoseGeometry(hgc_nose_geometry);

  heOffset_ = eeTopology().dddConstants().layers(true);
  totalLayers_ = heOffset_ + hsiTopology().dddConstants().layers(true);

  trigger_layers_.resize(totalLayers_ + 1);
  trigger_layers_[0] = 0;  // layer number 0 doesn't exist
  unsigned trigger_layer = 1;
  for (unsigned layer = 1; layer < trigger_layers_.size(); layer++) {
    if (disconnected_layers_.find(layer) == disconnected_layers_.end()) {
      // Increase trigger layer number if the layer is not disconnected
      trigger_layers_[layer] = trigger_layer;
      trigger_layer++;
    } else {
      trigger_layers_[layer] = 0;
    }
  }
  last_trigger_layer_ = trigger_layer - 1;
  fillMaps();

  noseLayers_ = noseTopology().dddConstants().layers(true);

  trigger_nose_layers_.resize(noseLayers_ + 1);
  trigger_nose_layers_[0] = 0;  // layer number 0 doesn't exist
  unsigned trigger_nose_layer = 1;
  for (unsigned layer = 1; layer < trigger_nose_layers_.size(); layer++) {
    trigger_nose_layers_[layer] = trigger_nose_layer;
    trigger_nose_layer++;
  }
}

unsigned HGCalTriggerGeometryV9Imp4::getTriggerCellFromCell(const unsigned cell_id) const {
  unsigned det = DetId(cell_id).det();
  unsigned trigger_cell_id = 0;
  // Scintillator
  if (det == DetId::HGCalHSc) {
    // Very rough mapping from cells to TC
    HGCScintillatorDetId cell_sc_id(cell_id);
    int ieta = ((cell_sc_id.ietaAbs() - 1) / hSc_triggercell_size_ + 1) * cell_sc_id.zside();
    int iphi = (cell_sc_id.iphi() - 1) / hSc_triggercell_size_ + 1;
    trigger_cell_id = HGCScintillatorDetId(cell_sc_id.type(), cell_sc_id.layer(), ieta, iphi);
  }
  // HFNose
  else if (det == DetId::Forward && DetId(cell_id).subdetId() == ForwardSubdetector::HFNose) {
    HFNoseDetId cell_nose_id(cell_id);
    trigger_cell_id = HFNoseTriggerDetId(HGCalTriggerSubdetector::HFNoseTrigger,
                                         cell_nose_id.zside(),
                                         cell_nose_id.type(),
                                         cell_nose_id.layer(),
                                         cell_nose_id.waferU(),
                                         cell_nose_id.waferV(),
                                         cell_nose_id.triggerCellU(),
                                         cell_nose_id.triggerCellV());
  }
  // Silicon
  else if (det == DetId::HGCalEE || det == DetId::HGCalHSi) {
    HGCSiliconDetId cell_si_id(cell_id);
    trigger_cell_id = HGCalTriggerDetId(
        det == DetId::HGCalEE ? HGCalTriggerSubdetector::HGCalEETrigger : HGCalTriggerSubdetector::HGCalHSiTrigger,
        cell_si_id.zside(),
        cell_si_id.type(),
        cell_si_id.layer(),
        cell_si_id.waferU(),
        cell_si_id.waferV(),
        cell_si_id.triggerCellU(),
        cell_si_id.triggerCellV());
  }
  return trigger_cell_id;
}

unsigned HGCalTriggerGeometryV9Imp4::getModuleFromCell(const unsigned cell_id) const {
  return getModuleFromTriggerCell(getTriggerCellFromCell(cell_id));
}

unsigned HGCalTriggerGeometryV9Imp4::getModuleFromTriggerCell(const unsigned trigger_cell_id) const {
  unsigned det = DetId(trigger_cell_id).det();
  int zside = 0;
  unsigned tc_type = 1;
  unsigned layer = 0;
  unsigned module_id = 0;

  // Scintillator
  if (det == DetId::HGCalHSc) {
    HGCScintillatorDetId trigger_cell_sc_id(trigger_cell_id);
    tc_type = trigger_cell_sc_id.type();
    layer = trigger_cell_sc_id.layer();
    zside = trigger_cell_sc_id.zside();
    int tc_eta = trigger_cell_sc_id.ietaAbs();
    int tc_phi = trigger_cell_sc_id.iphi();
    unsigned sector = tcEtaphiMappingToSector0(tc_eta, tc_phi);
    int ieta = 0;
    int iphi = 0;
    getScintillatoriEtaiPhi(ieta, iphi, tc_eta, tc_phi, layer);
    module_id =
        HGCalTriggerModuleDetId(HGCalTriggerSubdetector::HGCalHScTrigger, zside, tc_type, layer, sector, ieta, iphi);
  }
  // HFNose
  else if (det == DetId::HGCalTrigger and
           HGCalTriggerDetId(trigger_cell_id).subdet() == HGCalTriggerSubdetector::HFNoseTrigger) {
    HFNoseTriggerDetId trigger_cell_trig_id(trigger_cell_id);
    tc_type = trigger_cell_trig_id.type();
    layer = trigger_cell_trig_id.layer();
    zside = trigger_cell_trig_id.zside();
    int waferu = trigger_cell_trig_id.waferU();
    int waferv = trigger_cell_trig_id.waferV();
    unsigned sector = geom_rotation_120_.uvMappingToSector0(
        getWaferCentring(layer, HGCalTriggerSubdetector::HFNoseTrigger), waferu, waferv);
    module_id =
        HGCalTriggerModuleDetId(HGCalTriggerSubdetector::HFNoseTrigger, zside, tc_type, layer, sector, waferu, waferv);
  }
  // Silicon
  else {
    HGCalTriggerDetId trigger_cell_trig_id(trigger_cell_id);
    unsigned subdet = trigger_cell_trig_id.subdet();
    tc_type = trigger_cell_trig_id.type();
    layer = trigger_cell_trig_id.layer();
    zside = trigger_cell_trig_id.zside();
    int waferu = trigger_cell_trig_id.waferU();
    int waferv = trigger_cell_trig_id.waferV();
    unsigned sector = geom_rotation_120_.uvMappingToSector0(getWaferCentring(layer, subdet), waferu, waferv);
    module_id = HGCalTriggerModuleDetId(HGCalTriggerSubdetector(subdet), zside, tc_type, layer, sector, waferu, waferv);
  }
  return module_id;
}

HGCalTriggerGeometryBase::geom_set HGCalTriggerGeometryV9Imp4::getCellsFromTriggerCell(
    const unsigned trigger_cell_id) const {
  DetId trigger_cell_det_id(trigger_cell_id);
  unsigned det = trigger_cell_det_id.det();
  geom_set cell_det_ids;
  // Scintillator
  if (det == DetId::HGCalHSc) {
    HGCScintillatorDetId trigger_cell_sc_id(trigger_cell_id);
    int ieta0 = (trigger_cell_sc_id.ietaAbs() - 1) * hSc_triggercell_size_ + 1;
    int iphi0 = (trigger_cell_sc_id.iphi() - 1) * hSc_triggercell_size_ + 1;
    for (int ietaAbs = ieta0; ietaAbs < ieta0 + (int)hSc_triggercell_size_; ietaAbs++) {
      int ieta = ietaAbs * trigger_cell_sc_id.zside();
      for (int iphi = iphi0; iphi < iphi0 + (int)hSc_triggercell_size_; iphi++) {
        unsigned cell_id = HGCScintillatorDetId(trigger_cell_sc_id.type(), trigger_cell_sc_id.layer(), ieta, iphi);
        if (validCellId(DetId::HGCalHSc, cell_id))
          cell_det_ids.emplace(cell_id);
      }
    }
  }
  // HFNose
  else if (det == DetId::HGCalTrigger and
           HGCalTriggerDetId(trigger_cell_id).subdet() == HGCalTriggerSubdetector::HFNoseTrigger) {
    HFNoseTriggerDetId trigger_cell_nose_id(trigger_cell_id);
    int layer = trigger_cell_nose_id.layer();
    int zside = trigger_cell_nose_id.zside();
    int type = trigger_cell_nose_id.type();
    int waferu = trigger_cell_nose_id.waferU();
    int waferv = trigger_cell_nose_id.waferV();
    std::vector<int> cellus = trigger_cell_nose_id.cellU();
    std::vector<int> cellvs = trigger_cell_nose_id.cellV();
    for (unsigned ic = 0; ic < cellus.size(); ic++) {
      HFNoseDetId cell_det_id(zside, type, layer, waferu, waferv, cellus[ic], cellvs[ic]);
      cell_det_ids.emplace(cell_det_id);
    }
  }
  // Silicon
  else {
    HGCalTriggerDetId trigger_cell_trig_id(trigger_cell_id);
    unsigned subdet = trigger_cell_trig_id.subdet();
    if (subdet == HGCalTriggerSubdetector::HGCalEETrigger || subdet == HGCalTriggerSubdetector::HGCalHSiTrigger) {
      DetId::Detector cell_det = (subdet == HGCalTriggerSubdetector::HGCalEETrigger ? DetId::HGCalEE : DetId::HGCalHSi);
      int layer = trigger_cell_trig_id.layer();
      int zside = trigger_cell_trig_id.zside();
      int type = trigger_cell_trig_id.type();
      int waferu = trigger_cell_trig_id.waferU();
      int waferv = trigger_cell_trig_id.waferV();
      std::vector<int> cellus = trigger_cell_trig_id.cellU();
      std::vector<int> cellvs = trigger_cell_trig_id.cellV();
      for (unsigned ic = 0; ic < cellus.size(); ic++) {
        HGCSiliconDetId cell_det_id(cell_det, zside, type, layer, waferu, waferv, cellus[ic], cellvs[ic]);
        cell_det_ids.emplace(cell_det_id);
      }
    }
  }
  return cell_det_ids;
}

HGCalTriggerGeometryBase::geom_set HGCalTriggerGeometryV9Imp4::getCellsFromModule(const unsigned module_id) const {
  geom_set cell_det_ids;
  geom_set trigger_cells = getTriggerCellsFromModule(module_id);

  for (auto trigger_cell_id : trigger_cells) {
    geom_set cells = getCellsFromTriggerCell(trigger_cell_id);
    cell_det_ids.insert(cells.begin(), cells.end());
  }
  return cell_det_ids;
}

HGCalTriggerGeometryBase::geom_ordered_set HGCalTriggerGeometryV9Imp4::getOrderedCellsFromModule(
    const unsigned module_id) const {
  geom_ordered_set cell_det_ids;
  geom_ordered_set trigger_cells = getOrderedTriggerCellsFromModule(module_id);
  for (auto trigger_cell_id : trigger_cells) {
    geom_set cells = getCellsFromTriggerCell(trigger_cell_id);
    cell_det_ids.insert(cells.begin(), cells.end());
  }
  return cell_det_ids;
}

HGCalTriggerGeometryBase::geom_set HGCalTriggerGeometryV9Imp4::getTriggerCellsFromModule(
    const unsigned module_id) const {
  HGCalTriggerModuleDetId hgc_module_id(module_id);
  unsigned subdet = hgc_module_id.triggerSubdetId();

  geom_set trigger_cell_det_ids;
  // Scintillator
  if (subdet == HGCalTriggerSubdetector::HGCalHScTrigger) {
    int ieta0 = hgc_module_id.eta();
    int iphi0 = hgc_module_id.phi();

    unsigned layer = hgc_module_id.layer();
    etaphiMappingFromSector0(ieta0, iphi0, hgc_module_id.sector());
    int split = (layer > hSc_layer_for_split_) ? hSc_back_layers_split_ : hSc_front_layers_split_;
    if (ieta0 == 1) {
      ieta0 = ieta0 + split;
    } else {
      ieta0 = ieta0 + 1;
    }
    iphi0 = (iphi0 * hSc_tcs_per_module_phi_) + hSc_tc_layer0_min_ + 1;
    int total_tcs = hSc_num_panels_per_sector_ * hSc_tcs_per_module_phi_ * nSectors_;
    if (iphi0 > total_tcs) {
      iphi0 = iphi0 - total_tcs;
    }

    int hSc_tcs_per_module_eta = (layer > hSc_layer_for_split_) ? hSc_back_layers_split_ : hSc_front_layers_split_;

    for (int ietaAbs = ieta0; ietaAbs < ieta0 + (int)hSc_tcs_per_module_eta; ietaAbs++) {
      int ieta = ietaAbs * hgc_module_id.zside();
      for (int iphi = iphi0; iphi < iphi0 + (int)hSc_tcs_per_module_phi_; iphi++) {
        unsigned trigger_cell_id = HGCScintillatorDetId(hgc_module_id.type(), hgc_module_id.layer(), ieta, iphi);
        if (validTriggerCellFromCells(trigger_cell_id))
          trigger_cell_det_ids.emplace(trigger_cell_id);
      }
    }
  }
  // HFNose
  else if (subdet == HGCalTriggerSubdetector::HFNoseTrigger) {
    HFNoseDetId module_nose_id(module_id);
    HFNoseDetIdToModule hfn;
    std::vector<HFNoseTriggerDetId> ids = hfn.getTriggerDetIds(module_nose_id);
    for (auto const& idx : ids) {
      if (validTriggerCellFromCells(idx.rawId()))
        trigger_cell_det_ids.emplace(idx);
    }
  }
  // Silicon
  else {
    HGCSiliconDetIdToROC tc2roc;
    int moduleU = hgc_module_id.moduleU();
    int moduleV = hgc_module_id.moduleV();
    unsigned layer = hgc_module_id.layer();

    //Rotate to sector
    geom_rotation_120_.uvMappingFromSector0(getWaferCentring(layer, subdet), moduleU, moduleV, hgc_module_id.sector());

    DetId::Detector det = (subdet == HGCalTriggerSubdetector::HGCalEETrigger ? DetId::HGCalEE : DetId::HGCalHSi);

    unsigned wafer_type = detIdWaferType(det, layer, moduleU, moduleV);
    int nroc = (wafer_type == HGCSiliconDetId::HGCalFine ? 6 : 3);
    // Loop on ROCs in wafer
    for (int roc = 1; roc <= nroc; roc++) {
      // loop on TCs in ROC
      auto tc_uvs = tc2roc.getTriggerId(roc, wafer_type);
      for (const auto& tc_uv : tc_uvs) {
        HGCalTriggerDetId trigger_cell_id(
            subdet, hgc_module_id.zside(), wafer_type, layer, moduleU, moduleV, tc_uv.first, tc_uv.second);
        if (validTriggerCellFromCells(trigger_cell_id.rawId()))
          trigger_cell_det_ids.emplace(trigger_cell_id);
      }
    }
  }

  return trigger_cell_det_ids;
}

HGCalTriggerGeometryBase::geom_ordered_set HGCalTriggerGeometryV9Imp4::getOrderedTriggerCellsFromModule(
    const unsigned module_id) const {
  HGCalTriggerModuleDetId hgc_module_id(module_id);
  unsigned subdet = hgc_module_id.triggerSubdetId();

  geom_ordered_set trigger_cell_det_ids;

  // Scintillator
  if (subdet == HGCalTriggerSubdetector::HGCalHScTrigger) {
    int ieta0 = hgc_module_id.eta();
    int iphi0 = hgc_module_id.phi();

    unsigned layer = hgc_module_id.layer();
    etaphiMappingFromSector0(ieta0, iphi0, hgc_module_id.sector());
    int split = (layer > hSc_layer_for_split_) ? hSc_back_layers_split_ : hSc_front_layers_split_;
    if (ieta0 == 1) {
      ieta0 = ieta0 + split;
    } else {
      ieta0 = ieta0 + 1;
    }
    iphi0 = (iphi0 * hSc_tcs_per_module_phi_) + hSc_tc_layer0_min_ + 1;
    int total_tcs = hSc_num_panels_per_sector_ * hSc_tcs_per_module_phi_ * nSectors_;
    if (iphi0 > total_tcs) {
      iphi0 = iphi0 - total_tcs;
    }

    int hSc_tcs_per_module_eta = (layer > hSc_layer_for_split_) ? hSc_back_layers_split_ : hSc_front_layers_split_;

    for (int ietaAbs = ieta0; ietaAbs < ieta0 + (int)hSc_tcs_per_module_eta; ietaAbs++) {
      int ieta = ietaAbs * hgc_module_id.zside();
      for (int iphi = iphi0; iphi < iphi0 + (int)hSc_tcs_per_module_phi_; iphi++) {
        unsigned trigger_cell_id = HGCScintillatorDetId(hgc_module_id.type(), hgc_module_id.layer(), ieta, iphi);
        if (validTriggerCellFromCells(trigger_cell_id))
          trigger_cell_det_ids.emplace(trigger_cell_id);
      }
    }
  }

  // HFNose
  else if (subdet == HGCalTriggerSubdetector::HFNoseTrigger) {
    HFNoseDetId module_nose_id(module_id);
    HFNoseDetIdToModule hfn;
    std::vector<HFNoseTriggerDetId> ids = hfn.getTriggerDetIds(module_nose_id);
    for (auto const& idx : ids) {
      if (validTriggerCellFromCells(idx.rawId()))
        trigger_cell_det_ids.emplace(idx);
    }
  }
  // Silicon
  else {
    HGCSiliconDetIdToROC tc2roc;
    int moduleU = hgc_module_id.moduleU();
    int moduleV = hgc_module_id.moduleV();
    unsigned layer = hgc_module_id.layer();

    //Rotate to sector
    geom_rotation_120_.uvMappingFromSector0(getWaferCentring(layer, subdet), moduleU, moduleV, hgc_module_id.sector());

    DetId::Detector det = (subdet == HGCalTriggerSubdetector::HGCalEETrigger ? DetId::HGCalEE : DetId::HGCalHSi);

    unsigned wafer_type = detIdWaferType(det, layer, moduleU, moduleV);
    int nroc = (wafer_type == HGCSiliconDetId::HGCalFine ? 6 : 3);
    // Loop on ROCs in wafer
    for (int roc = 1; roc <= nroc; roc++) {
      // loop on TCs in ROC
      auto tc_uvs = tc2roc.getTriggerId(roc, wafer_type);
      for (const auto& tc_uv : tc_uvs) {
        HGCalTriggerDetId trigger_cell_id(
            subdet, hgc_module_id.zside(), wafer_type, layer, moduleU, moduleV, tc_uv.first, tc_uv.second);
        if (validTriggerCellFromCells(trigger_cell_id.rawId()))
          trigger_cell_det_ids.emplace(trigger_cell_id);
      }
    }
  }

  return trigger_cell_det_ids;
}

HGCalTriggerGeometryBase::geom_set HGCalTriggerGeometryV9Imp4::getNeighborsFromTriggerCell(
    const unsigned trigger_cell_id) const {
  throw cms::Exception("FeatureNotImplemented") << "Neighbor search is not implemented in HGCalTriggerGeometryV9Imp4";
}

unsigned HGCalTriggerGeometryV9Imp4::getLinksInModule(const unsigned module_id) const {
  HGCalTriggerModuleDetId module_det_id(module_id);
  unsigned subdet = module_det_id.triggerSubdetId();

  unsigned links = 0;
  // HF Nose
  if (subdet == HGCalTriggerSubdetector::HFNoseTrigger) {
    links = 1;
  }
  // TO ADD HFNOSE : getLinksInModule
  // Silicon and Scintillator
  else {
    int packed_module =
        packLayerSubdetWaferId(module_det_id.layer(), subdet, module_det_id.moduleU(), module_det_id.moduleV());
    links = links_per_module_.at(packed_module);
  }
  return links;
}

unsigned HGCalTriggerGeometryV9Imp4::getModuleSize(const unsigned module_id) const {
  unsigned nWafers = 1;
  return nWafers;
}

unsigned HGCalTriggerGeometryV9Imp4::getNextSector(const unsigned sector) const {
  unsigned next_sector = 0;
  if (sector < 2) {
    next_sector = sector + 1;
  }
  return next_sector;
}

unsigned HGCalTriggerGeometryV9Imp4::getPreviousSector(const unsigned sector) const {
  unsigned previous_sector = 2;
  if (sector > 0) {
    previous_sector = sector - 1;
  }
  return previous_sector;
}

HGCalTriggerGeometryBase::geom_set HGCalTriggerGeometryV9Imp4::getStage1FpgasFromStage2Fpga(
    const unsigned stage2_id) const {
  geom_set stage1_ids;
  geom_set stage1_links = getStage1LinksFromStage2Fpga(stage2_id);
  for (const auto& stage1_link : stage1_links) {
    HGCalTriggerBackendDetId id(stage1_link);
    const auto s1FPGA = getStage1FpgaFromStage1Link(stage1_link);
    stage1_ids.emplace(getStage1FpgaFromStage1Link(stage1_link));
  }
  return stage1_ids;
}

HGCalTriggerGeometryBase::geom_set HGCalTriggerGeometryV9Imp4::getStage2FpgasFromStage1Fpga(
    const unsigned stage1_id) const {
  geom_set stage2_ids;

  geom_set stage1_links = getStage1LinksFromStage1Fpga(stage1_id);
  for (const auto& stage1_link : stage1_links) {
    stage2_ids.emplace(getStage2FpgaFromStage1Link(stage1_link));
  }

  return stage2_ids;
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

unsigned HGCalTriggerGeometryV9Imp4::getStage1FpgaFromStage1Link(const unsigned link_id) const {
  HGCalTriggerBackendDetId id(link_id);
  unsigned stage1_label = stage1link_to_stage1_.at(id.label());

  return HGCalTriggerBackendDetId(
      id.zside(), HGCalTriggerBackendDetId::BackendType::Stage1FPGA, id.sector(), stage1_label);
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

HGCalTriggerGeometryBase::geom_set HGCalTriggerGeometryV9Imp4::getStage1LinksFromStage1Fpga(
    const unsigned stage1_id) const {
  geom_set stage1link_ids;
  HGCalTriggerBackendDetId id(stage1_id);
  auto stage1_itrs = stage1_to_stage1links_.equal_range(id.label());
  for (auto stage1_itr = stage1_itrs.first; stage1_itr != stage1_itrs.second; stage1_itr++) {
    stage1link_ids.emplace(HGCalTriggerBackendDetId(
        id.zside(), HGCalTriggerBackendDetId::BackendType::Stage1Link, id.sector(), stage1_itr->second));
  }
  return stage1link_ids;
}

std::vector<unsigned> HGCalTriggerGeometryV9Imp4::getLpgbtsFromStage1Fpga(const unsigned stage1_id) const {
  std::vector<unsigned> lpgbt_ids;
  HGCalTriggerBackendDetId id(stage1_id);
  const auto stage1_lpgbts = stage1_to_lpgbts_.at(id.label());
  lpgbt_ids.reserve(stage1_lpgbts.size());
  for (const auto& stage1_lpgbt : stage1_lpgbts) {
    lpgbt_ids.emplace_back(
        HGCalTriggerBackendDetId(id.zside(), HGCalTriggerBackendDetId::BackendType::LpGBT, id.sector(), stage1_lpgbt));
  }

  return lpgbt_ids;
}

unsigned HGCalTriggerGeometryV9Imp4::getStage1FpgaFromLpgbt(const unsigned lpgbt_id) const {
  HGCalTriggerBackendDetId id(lpgbt_id);
  unsigned stage1_label = lpgbt_to_stage1_.at(id.label());

  return HGCalTriggerBackendDetId(
      id.zside(), HGCalTriggerBackendDetId::BackendType::Stage1FPGA, id.sector(), stage1_label);
}

HGCalTriggerGeometryBase::geom_set HGCalTriggerGeometryV9Imp4::getModulesFromLpgbt(const unsigned lpgbt_id) const {
  geom_set modules;
  HGCalTriggerBackendDetId id(lpgbt_id);

  auto lpgbt_itrs = lpgbt_to_modules_.equal_range(id.label());
  for (auto lpgbt_itr = lpgbt_itrs.first; lpgbt_itr != lpgbt_itrs.second; lpgbt_itr++) {
    unsigned layer = 0;
    int moduleU = 0;
    int moduleV = 0;
    int subdet = 0;
    unpackLayerSubdetWaferId(lpgbt_itr->second, layer, subdet, moduleU, moduleV);
    unsigned det = 0;
    switch (subdet) {
      case HGCalTriggerSubdetector::HGCalEETrigger:
        det = DetId::HGCalEE;
        break;
      case HGCalTriggerSubdetector::HGCalHSiTrigger:
        det = DetId::HGCalHSi;
        break;
      case HGCalTriggerSubdetector::HGCalHScTrigger:
        det = DetId::HGCalHSc;
        break;
      default:
        det = DetId::HGCalEE;
        break;
    }

    int type = detIdWaferType(det, layer, moduleU, moduleV);
    modules.emplace(HGCalTriggerModuleDetId(
        HGCalTriggerSubdetector(subdet), id.zside(), type, layer, id.sector(), moduleU, moduleV));
  }
  return modules;
}

HGCalTriggerGeometryV9Imp4::geom_set HGCalTriggerGeometryV9Imp4::getLpgbtsFromModule(const unsigned module_id) const {
  geom_set lpgbt_ids;
  HGCalTriggerModuleDetId id(module_id);

  auto module_itrs = module_to_lpgbts_.equal_range(
      packLayerSubdetWaferId(id.layer(), id.triggerSubdetId(), id.moduleU(), id.moduleV()));
  for (auto module_itr = module_itrs.first; module_itr != module_itrs.second; module_itr++) {
    lpgbt_ids.emplace(HGCalTriggerBackendDetId(
        id.zside(), HGCalTriggerBackendDetId::BackendType::LpGBT, id.sector(), module_itr->second));
  }

  return lpgbt_ids;
}

unsigned HGCalTriggerGeometryV9Imp4::getStage1FpgaFromModule(const unsigned module_id) const {
  HGCalTriggerModuleDetId id(module_id);

  unsigned stage1_label =
      module_to_stage1_.at(packLayerSubdetWaferId(id.layer(), id.triggerSubdetId(), id.moduleU(), id.moduleV()));

  return HGCalTriggerBackendDetId(
      id.zside(), HGCalTriggerBackendDetId::BackendType::Stage1FPGA, id.sector(), stage1_label);
}

GlobalPoint HGCalTriggerGeometryV9Imp4::getTriggerCellPosition(const unsigned trigger_cell_det_id) const {
  unsigned det = DetId(trigger_cell_det_id).det();

  // Position: barycenter of the trigger cell.
  Basic3DVector<float> triggerCellVector(0., 0., 0.);
  const auto cell_ids = getCellsFromTriggerCell(trigger_cell_det_id);
  // Scintillator
  if (det == DetId::HGCalHSc) {
    for (const auto& cell : cell_ids) {
      triggerCellVector += hscGeometry()->getPosition(cell).basicVector();
    }
  }
  // HFNose
  else if (det == DetId::HGCalTrigger and
           HGCalTriggerDetId(trigger_cell_det_id).subdet() == HGCalTriggerSubdetector::HFNoseTrigger) {
    for (const auto& cell : cell_ids) {
      HFNoseDetId cellDetId(cell);
      triggerCellVector += noseGeometry()->getPosition(cellDetId).basicVector();
    }
  }
  // Silicon
  else {
    for (const auto& cell : cell_ids) {
      HGCSiliconDetId cellDetId(cell);
      triggerCellVector += (cellDetId.det() == DetId::HGCalEE ? eeGeometry()->getPosition(cellDetId)
                                                              : hsiGeometry()->getPosition(cellDetId))
                               .basicVector();
    }
  }
  return GlobalPoint(triggerCellVector / cell_ids.size());
}

GlobalPoint HGCalTriggerGeometryV9Imp4::getModulePosition(const unsigned module_det_id) const {
  unsigned subdet = HGCalTriggerModuleDetId(module_det_id).triggerSubdetId();
  // Position: barycenter of the module.
  Basic3DVector<float> moduleVector(0., 0., 0.);
  const auto cell_ids = getCellsFromModule(module_det_id);
  // Scintillator
  if (subdet == HGCalTriggerSubdetector::HGCalHScTrigger) {
    for (const auto& cell : cell_ids) {
      moduleVector += hscGeometry()->getPosition(cell).basicVector();
    }
  }
  // HFNose
  else if (subdet == HGCalTriggerSubdetector::HFNoseTrigger) {
    for (const auto& cell : cell_ids) {
      HFNoseDetId cellDetId(cell);
      moduleVector += noseGeometry()->getPosition(cellDetId).basicVector();
    }
  }  // Silicon
  else {
    for (const auto& cell : cell_ids) {
      HGCSiliconDetId cellDetId(cell);
      moduleVector += (cellDetId.det() == DetId::HGCalEE ? eeGeometry()->getPosition(cellDetId)
                                                         : hsiGeometry()->getPosition(cellDetId))
                          .basicVector();
    }
  }

  return GlobalPoint(moduleVector / cell_ids.size());
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
  const unsigned int stage2_id = 0;

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
      const auto s1Subsector = l1thgcmapping::get_subsector(belink);
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

unsigned HGCalTriggerGeometryV9Imp4::packLayerSubdetWaferId(unsigned layer, int subdet, int waferU, int waferV) const {
  unsigned packed_value = 0;

  packed_value |=
      ((waferU & HGCalTriggerModuleDetId::kHGCalModuleUMask) << HGCalTriggerModuleDetId::kHGCalModuleUOffset);
  packed_value |=
      ((waferV & HGCalTriggerModuleDetId::kHGCalModuleVMask) << HGCalTriggerModuleDetId::kHGCalModuleVOffset);
  packed_value |= ((subdet & HGCalTriggerModuleDetId::kHGCalTriggerSubdetMask)
                   << HGCalTriggerModuleDetId::kHGCalTriggerSubdetOffset);
  packed_value |= ((layer & HGCalTriggerModuleDetId::kHGCalLayerMask) << HGCalTriggerModuleDetId::kHGCalLayerOffset);
  return packed_value;
}

void HGCalTriggerGeometryV9Imp4::unpackLayerSubdetWaferId(
    unsigned wafer, unsigned& layer, int& subdet, int& waferU, int& waferV) const {
  waferU = (wafer >> HGCalTriggerModuleDetId::kHGCalModuleUOffset) & HGCalTriggerModuleDetId::kHGCalModuleUMask;
  waferV = (wafer >> HGCalTriggerModuleDetId::kHGCalModuleVOffset) & HGCalTriggerModuleDetId::kHGCalModuleVMask;
  subdet =
      (wafer >> HGCalTriggerModuleDetId::kHGCalTriggerSubdetOffset) & HGCalTriggerModuleDetId::kHGCalTriggerSubdetMask;
  layer = (wafer >> HGCalTriggerModuleDetId::kHGCalLayerOffset) & HGCalTriggerModuleDetId::kHGCalLayerMask;
}

void HGCalTriggerGeometryV9Imp4::etaphiMappingFromSector0(int& ieta, int& iphi, unsigned sector) const {
  if (sector == 0) {
    return;
  }
  if (sector == 2) {
    iphi = iphi + hSc_num_panels_per_sector_;
  } else if (sector == 1) {
    iphi = iphi + (2 * hSc_num_panels_per_sector_);
  }
}

HGCalGeomRotation::WaferCentring HGCalTriggerGeometryV9Imp4::getWaferCentring(unsigned layer, int subdet) const {
  if (subdet == HGCalTriggerSubdetector::HGCalEETrigger) {  // CE-E
    return HGCalGeomRotation::WaferCentring::WaferCentred;
  } else if (subdet == HGCalTriggerSubdetector::HGCalHSiTrigger) {
    if ((layer % 2) == 1) {  // CE-H Odd
      return HGCalGeomRotation::WaferCentring::CornerCentredY;
    } else {  // CE-H Even
      return HGCalGeomRotation::WaferCentring::CornerCentredMercedes;
    }
  } else if (subdet == HGCalTriggerSubdetector::HFNoseTrigger) {  //HFNose
    return HGCalGeomRotation::WaferCentring::WaferCentred;
  } else {
    edm::LogError("HGCalTriggerGeometryV9Imp4")
        << "HGCalTriggerGeometryV9Imp4: trigger sub-detector expected to be silicon";
    return HGCalGeomRotation::WaferCentring::WaferCentred;
  }
}

unsigned HGCalTriggerGeometryV9Imp4::tcEtaphiMappingToSector0(int& tc_ieta, int& tc_iphi) const {
  unsigned sector = 0;

  if (tc_iphi > hSc_tc_layer0_min_ && tc_iphi <= hSc_tc_layer0_min_ + ntc_per_wafer_) {
    sector = 0;
  } else if (tc_iphi > hSc_tc_layer0_min_ + ntc_per_wafer_ && tc_iphi <= hSc_tc_layer0_min_ + 2 * ntc_per_wafer_) {
    sector = 2;
  } else {
    sector = 1;
  }

  if (sector == 0) {
    tc_iphi = tc_iphi - hSc_tc_layer0_min_;
  } else if (sector == 2) {
    tc_iphi = tc_iphi - (hSc_tc_layer0_min_ + ntc_per_wafer_);
  } else if (sector == 1) {
    if (tc_iphi <= hSc_tc_layer0_min_) {
      tc_iphi = tc_iphi + nSectors_ * ntc_per_wafer_;
    }
    tc_iphi = tc_iphi - (nSectors_ * ntc_per_wafer_ - hSc_tc_layer0_min_);
  }

  return sector;
}

void HGCalTriggerGeometryV9Imp4::getScintillatoriEtaiPhi(
    int& ieta, int& iphi, int tc_eta, int tc_phi, unsigned layer) const {
  iphi = (tc_phi - 1) / hSc_tcs_per_module_phi_;  //Phi index 1-12

  int split = hSc_front_layers_split_;
  if (layer > hSc_layer_for_split_) {
    split = hSc_back_layers_split_;
  }
  if (tc_eta <= split) {
    ieta = 0;
  } else {
    ieta = 1;
  }
}

bool HGCalTriggerGeometryV9Imp4::validTriggerCell(const unsigned trigger_cell_id) const {
  return validTriggerCellFromCells(trigger_cell_id);
}

bool HGCalTriggerGeometryV9Imp4::disconnectedModule(const unsigned module_id) const {
  bool disconnected = false;
  HGCalTriggerModuleDetId id(module_id);
  if (module_to_stage1_.find(packLayerSubdetWaferId(id.layer(), id.triggerSubdetId(), id.moduleU(), id.moduleV())) ==
      module_to_stage1_.end()) {
    disconnected = true;
  }
  if (disconnected_layers_.find(layerWithOffset(module_id)) != disconnected_layers_.end()) {
    disconnected = true;
  }
  return disconnected;
}

unsigned HGCalTriggerGeometryV9Imp4::triggerLayer(const unsigned id) const {
  unsigned layer = layerWithOffset(id);

  if (DetId(id).det() == DetId::HGCalTrigger and
      HGCalTriggerDetId(id).subdet() == HGCalTriggerSubdetector::HFNoseTrigger) {
    if (layer >= trigger_nose_layers_.size())
      return 0;
    return trigger_nose_layers_[layer];
  }
  if (layer >= trigger_layers_.size())
    return 0;
  return trigger_layers_[layer];
}

bool HGCalTriggerGeometryV9Imp4::validCell(unsigned cell_id) const {
  bool is_valid = false;
  unsigned det = DetId(cell_id).det();
  switch (det) {
    case DetId::HGCalEE:
      is_valid = eeTopology().valid(cell_id);
      break;
    case DetId::HGCalHSi:
      is_valid = hsiTopology().valid(cell_id);
      break;
    case DetId::HGCalHSc:
      is_valid = hscTopology().valid(cell_id);
      break;
    case DetId::Forward:
      is_valid = noseTopology().valid(cell_id);
      break;
    default:
      is_valid = false;
      break;
  }
  return is_valid;
}

bool HGCalTriggerGeometryV9Imp4::validTriggerCellFromCells(const unsigned trigger_cell_id) const {
  // Check the validity of a trigger cell with the
  // validity of the cells. One valid cell in the
  // trigger cell is enough to make the trigger cell
  // valid.
  const geom_set cells = getCellsFromTriggerCell(trigger_cell_id);
  bool is_valid = false;
  for (const auto cell_id : cells) {
    unsigned det = DetId(cell_id).det();
    is_valid |= validCellId(det, cell_id);
    if (is_valid)
      break;
  }
  return is_valid;
}

bool HGCalTriggerGeometryV9Imp4::validCellId(unsigned subdet, unsigned cell_id) const {
  bool is_valid = false;
  switch (subdet) {
    case DetId::HGCalEE:
      is_valid = eeTopology().valid(cell_id);
      break;
    case DetId::HGCalHSi:
      is_valid = hsiTopology().valid(cell_id);
      break;
    case DetId::HGCalHSc:
      is_valid = hscTopology().valid(cell_id);
      break;
    case DetId::Forward:
      is_valid = noseTopology().valid(cell_id);
      break;
    default:
      is_valid = false;
      break;
  }
  return is_valid;
}

int HGCalTriggerGeometryV9Imp4::detIdWaferType(unsigned det, unsigned layer, short waferU, short waferV) const {
  int wafer_type = 0;
  switch (det) {
    case DetId::HGCalEE:
      wafer_type = eeTopology().dddConstants().getTypeHex(layer, waferU, waferV);
      break;
    case DetId::HGCalHSi:
      wafer_type = hsiTopology().dddConstants().getTypeHex(layer, waferU, waferV);
      break;
    case DetId::HGCalHSc:
      wafer_type = hscTopology().dddConstants().getTypeTrap(layer);
      break;
    default:
      break;
  };
  return wafer_type;
}

void HGCalTriggerGeometryV9Imp4::layerWithoutOffsetAndSubdetId(unsigned& layer, int& subdetId, bool isSilicon) const {
  if (!isSilicon) {
    layer = layer - heOffset_;
    subdetId = HGCalTriggerSubdetector::HGCalHScTrigger;
  } else {
    if (layer > heOffset_) {
      subdetId = HGCalTriggerSubdetector::HGCalHSiTrigger;
      layer = layer - heOffset_;
    } else {
      subdetId = HGCalTriggerSubdetector::HGCalEETrigger;
    }
  }
}

unsigned HGCalTriggerGeometryV9Imp4::layerWithOffset(unsigned id) const {
  unsigned det = DetId(id).det();
  unsigned layer = 0;

  if (det == DetId::HGCalTrigger) {
    unsigned subdet = HGCalTriggerDetId(id).subdet();
    if (subdet == HGCalTriggerSubdetector::HGCalEETrigger) {
      layer = HGCalTriggerDetId(id).layer();
    } else if (subdet == HGCalTriggerSubdetector::HGCalHSiTrigger) {
      layer = heOffset_ + HGCalTriggerDetId(id).layer();
    } else if (subdet == HGCalTriggerSubdetector::HFNoseTrigger) {
      layer = HFNoseTriggerDetId(id).layer();
    }
  } else if (det == DetId::HGCalHSc) {
    layer = heOffset_ + HGCScintillatorDetId(id).layer();
  } else if (det == DetId::Forward) {
    unsigned subdet = HGCalTriggerModuleDetId(id).triggerSubdetId();
    if (subdet == HGCalTriggerSubdetector::HGCalEETrigger) {
      layer = HGCalTriggerModuleDetId(id).layer();
    } else if (subdet == HGCalTriggerSubdetector::HGCalHSiTrigger ||
               subdet == HGCalTriggerSubdetector::HGCalHScTrigger) {
      layer = heOffset_ + HGCalDetId(id).layer();
    } else if (subdet == HGCalTriggerSubdetector::HFNoseTrigger) {
      layer = HGCalTriggerModuleDetId(id).layer();
    }
  }
  return layer;
}

DEFINE_EDM_PLUGIN(HGCalTriggerGeometryFactory, HGCalTriggerGeometryV9Imp4, "HGCalTriggerGeometryV9Imp4");
