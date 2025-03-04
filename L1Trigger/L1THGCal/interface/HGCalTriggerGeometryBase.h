#ifndef __L1Trigger_L1THGCal_HGCalTriggerGeometryBase_h__
#define __L1Trigger_L1THGCal_HGCalTriggerGeometryBase_h__

#include <iostream>
#include <unordered_map>
#include <unordered_set>

#include "DataFormats/GeometryVector/interface/GlobalPoint.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "Geometry/CaloTopology/interface/HGCalTopology.h"
#include "Geometry/HGCalGeometry/interface/HGCalGeometry.h"

// Pure virtual trigger geometry class
// Provides the interface to access trigger cell and module mappings
class HGCalTriggerGeometryBase {
public:
  typedef std::unordered_map<unsigned, unsigned> geom_map;
  typedef std::unordered_set<unsigned> geom_set;
  typedef std::set<unsigned> geom_ordered_set;

  HGCalTriggerGeometryBase(const edm::ParameterSet& conf);
  virtual ~HGCalTriggerGeometryBase() {}

  const std::string& name() const { return name_; }

  bool isWithNoseGeometry() const { return isNose_; }

  const HGCalGeometry* noseGeometry() const { return hgc_nose_geometry_; }
  const HGCalGeometry* eeGeometry() const { return hgc_ee_geometry_; }
  const HGCalGeometry* fhGeometry() const { return hgc_hsi_geometry_; }
  const HGCalGeometry* hsiGeometry() const { return fhGeometry(); }
  const HGCalGeometry* hscGeometry() const { return hgc_hsc_geometry_; }
  const HGCalTopology& noseTopology() const { return noseGeometry()->topology(); }
  const HGCalTopology& eeTopology() const { return eeGeometry()->topology(); }
  const HGCalTopology& fhTopology() const { return fhGeometry()->topology(); }
  const HGCalTopology& hsiTopology() const { return hsiGeometry()->topology(); }
  const HGCalTopology& hscTopology() const { return hscGeometry()->topology(); }

  void setWithNoseGeometry(const bool isNose) { isNose_ = isNose; }

  // non-const access to the geometry class
  virtual void initialize(const HGCalGeometry*, const HGCalGeometry*, const HGCalGeometry*) = 0;
  virtual void initialize(const HGCalGeometry*, const HGCalGeometry*, const HGCalGeometry*, const HGCalGeometry*) = 0;
  virtual void reset();

  // const access to the geometry class
  virtual unsigned getTriggerCellFromCell(const unsigned cell_det_id) const = 0;
  virtual unsigned getModuleFromCell(const unsigned cell_det_id) const = 0;
  virtual unsigned getModuleFromTriggerCell(const unsigned trigger_cell_det_id) const = 0;

  virtual geom_set getCellsFromTriggerCell(const unsigned cell_det_id) const = 0;
  virtual geom_set getCellsFromModule(const unsigned cell_det_id) const = 0;
  virtual geom_set getTriggerCellsFromModule(const unsigned trigger_cell_det_id) const = 0;

  virtual geom_set getStage1FpgasFromStage2Fpga(const unsigned stage2_id) const = 0;
  virtual geom_set getStage2FpgasFromStage1Fpga(const unsigned stage1_id) const = 0;

  virtual geom_set getStage1LinksFromStage2Fpga(const unsigned) const = 0;
  virtual unsigned getStage1FpgaFromStage1Link(const unsigned) const = 0;
  virtual unsigned getStage2FpgaFromStage1Link(const unsigned) const = 0;
  virtual geom_set getStage1LinksFromStage1Fpga(const unsigned) const = 0;
  virtual std::vector<unsigned> getLpgbtsFromStage1Fpga(const unsigned stage1_id) const = 0;
  virtual geom_set getModulesFromStage1Fpga(const unsigned stage1_id) const = 0;
  virtual unsigned getStage1FpgaFromLpgbt(const unsigned lpgbt_id) const = 0;
  virtual geom_set getModulesFromLpgbt(const unsigned lpgbt_id) const = 0;
  virtual geom_set getLpgbtsFromModule(const unsigned module_id) const = 0;
  virtual unsigned getStage1FpgaFromModule(const unsigned module_id) const = 0;

  virtual geom_ordered_set getOrderedCellsFromModule(const unsigned cell_det_id) const = 0;
  virtual geom_ordered_set getOrderedTriggerCellsFromModule(const unsigned trigger_cell_det_id) const = 0;

  virtual geom_set getNeighborsFromTriggerCell(const unsigned trigger_cell_det_id) const = 0;

  virtual unsigned getLinksInModule(const unsigned module_id) const = 0;
  virtual unsigned getModuleSize(const unsigned module_id) const = 0;

  virtual GlobalPoint getTriggerCellPosition(const unsigned trigger_cell_det_id) const = 0;
  virtual GlobalPoint getModulePosition(const unsigned module_det_id) const = 0;

  virtual bool validCell(const unsigned cell_id) const = 0;
  virtual bool validTriggerCell(const unsigned trigger_cell_id) const = 0;
  virtual bool disconnectedModule(const unsigned module_id) const = 0;
  virtual unsigned lastTriggerLayer() const = 0;
  virtual unsigned triggerLayer(const unsigned id) const = 0;
  virtual const std::vector<unsigned>& triggerLayers() const = 0;

protected:
  void setEEGeometry(const HGCalGeometry* geom) { hgc_ee_geometry_ = geom; }
  void setHSiGeometry(const HGCalGeometry* geom) { hgc_hsi_geometry_ = geom; }
  void setHScGeometry(const HGCalGeometry* geom) { hgc_hsc_geometry_ = geom; }
  void setNoseGeometry(const HGCalGeometry* geom) { hgc_nose_geometry_ = geom; }

private:
  const std::string name_;

  bool isNose_ = false;
  const HGCalGeometry* hgc_ee_geometry_ = nullptr;
  const HGCalGeometry* hgc_hsi_geometry_ = nullptr;
  const HGCalGeometry* hgc_hsc_geometry_ = nullptr;
  const HGCalGeometry* hgc_nose_geometry_ = nullptr;
};

#include "FWCore/PluginManager/interface/PluginFactory.h"
typedef edmplugin::PluginFactory<HGCalTriggerGeometryBase*(const edm::ParameterSet&)> HGCalTriggerGeometryFactory;

#endif
