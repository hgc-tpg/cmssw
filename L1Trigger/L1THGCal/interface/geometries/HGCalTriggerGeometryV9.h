#include "Geometry/HGCalCommonData/interface/HGCalGeomRotation.h"
#include "DataFormats/ForwardDetId/interface/ForwardSubdetector.h"
#include "DataFormats/ForwardDetId/interface/HFNoseTriggerDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCScintillatorDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCSiliconDetIdToROC.h"
#include "DataFormats/ForwardDetId/interface/HGCalDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCalTriggerDetId.h"
#include "L1Trigger/L1THGCal/interface/HGCalTriggerGeometryBase.h"
#include "DataFormats/ForwardDetId/interface/HGCalTriggerModuleDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCalTriggerBackendDetId.h"
#include "DataFormats/ForwardDetId/interface/HFNoseDetIdToModule.h"

#include <vector>

class HGCalTriggerGeometryV9 : public HGCalTriggerGeometryBase {
public:
  HGCalTriggerGeometryV9(const edm::ParameterSet& conf);

  void initialize(const HGCalGeometry*, const HGCalGeometry*, const HGCalGeometry*) final;
  void initialize(const HGCalGeometry*, const HGCalGeometry*, const HGCalGeometry*, const HGCalGeometry*) final;
  void reset() override;

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

  geom_set getStage1LinksFromStage2Fpga(const unsigned) const override = 0;
  unsigned getStage1FpgaFromStage1Link(const unsigned) const final;
  unsigned getStage2FpgaFromStage1Link(const unsigned) const override = 0;
  geom_set getStage1LinksFromStage1Fpga(const unsigned) const final;
  std::vector<unsigned> getLpgbtsFromStage1Fpga(const unsigned) const final;
  geom_set getModulesFromStage1Fpga(const unsigned) const final;
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

protected:
  // module related maps
  std::unordered_map<unsigned, unsigned> links_per_module_;

  std::unordered_multimap<unsigned, unsigned> stage2_to_stage1links_;
  std::unordered_map<unsigned, unsigned> stage1link_to_stage2_;
  std::unordered_map<unsigned, unsigned> stage1link_to_stage1_;
  std::unordered_multimap<unsigned, unsigned> stage1_to_stage1links_;
  std::unordered_map<unsigned, std::vector<unsigned>> stage1_to_lpgbts_;
  std::unordered_map<unsigned, unsigned> lpgbt_to_stage1_;
  std::unordered_multimap<unsigned, unsigned> lpgbt_to_modules_;
  std::unordered_multimap<unsigned, unsigned> module_to_lpgbts_;
  std::unordered_map<unsigned, unsigned> module_to_stage1_;

  virtual void fillMaps() = 0;

  unsigned packLayerSubdetWaferId(unsigned layer, int subdet, int waferU, int waferV) const;
  void layerWithoutOffsetAndSubdetId(unsigned& layer, int& subdetId, bool isSilicon) const;
  unsigned getNextSector(const unsigned sector) const;
  unsigned getPreviousSector(const unsigned sector) const;

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

  // rotation class
  HGCalGeomRotation geom_rotation_120_ = {HGCalGeomRotation::SectorType::Sector120Degrees};

  // Disconnected modules and layers
  std::unordered_set<unsigned> disconnected_layers_;
  std::vector<unsigned> trigger_layers_;
  std::vector<unsigned> trigger_nose_layers_;
  unsigned last_trigger_layer_ = 0;

  // layer offsets
  unsigned heOffset_ = 0;
  unsigned noseLayers_ = 0;
  unsigned totalLayers_ = 0;

  bool validCellId(unsigned det, unsigned cell_id) const;
  bool validTriggerCellFromCells(const unsigned) const;

  int detIdWaferType(unsigned det, unsigned layer, short waferU, short waferV) const;
  void unpackLayerSubdetWaferId(unsigned wafer, unsigned& layer, int& subdet, int& waferU, int& waferV) const;
  HGCalGeomRotation::WaferCentring getWaferCentring(unsigned layer, int subdet) const;
  void etaphiMappingFromSector0(int& ieta, int& iphi, unsigned sector) const;
  unsigned tcEtaphiMappingToSector0(int& tc_ieta, int& tc_iphi) const;
  void getScintillatoriEtaiPhi(int& ieta, int& iphi, int tc_eta, int tc_phi, unsigned layer) const;
  unsigned layerWithOffset(unsigned) const;
};
