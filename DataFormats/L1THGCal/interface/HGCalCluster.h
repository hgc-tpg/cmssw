#ifndef DataFormats_L1Trigger_HGCalCluster_h
#define DataFormats_L1Trigger_HGCalCluster_h

#include "DataFormats/Common/interface/Ptr.h"
#include "DataFormats/L1Trigger/interface/BXVector.h"
#include "DataFormats/L1THGCal/interface/HGCalTriggerCell.h"
#include "DataFormats/L1THGCal/interface/HGCalClusterT.h"

namespace l1t {

  class HGCalCluster : public HGCalClusterT<l1t::HGCalTriggerCell> {
  public:
    HGCalCluster() {}
    HGCalCluster(const LorentzVector p4, int pt = 0, int eta = 0, int phi = 0);

    HGCalCluster(const edm::Ptr<l1t::HGCalTriggerCell> &tc);

    ~HGCalCluster() override;

    void setModule(uint32_t module) { module_ = module; }
    uint32_t module() const { return module_; }

    void setColumn(int column) { column_ = column; }
    int column() const { return column_; }

    void setFrame(unsigned frame) { frame_ = frame; }
    unsigned frame() const { return frame_; }

    void setChannel(unsigned channel) { channel_ = channel; }
    unsigned channel() const { return channel_; }

  private:
    uint32_t module_;
    unsigned channel_;
    unsigned frame_;
    int column_;
  };

  typedef BXVector<HGCalCluster> HGCalClusterBxCollection;

}  // namespace l1t

#endif
