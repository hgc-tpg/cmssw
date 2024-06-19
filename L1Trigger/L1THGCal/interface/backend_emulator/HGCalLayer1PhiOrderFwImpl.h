#ifndef __L1Trigger_L1THGCal_HGCalLayer1PhiOrderFwImpl_h__
#define __L1Trigger_L1THGCal_HGCalLayer1PhiOrderFwImpl_h__

#include "L1Trigger/L1THGCal/interface/backend_emulator/HGCalTriggerCell_SA.h"
#include "L1Trigger/L1THGCal/interface/backend_emulator/HGCalLayer1PhiOrderFwConfig.h"
#include "L1Trigger/L1THGCal/interface/backend_emulator/BatcherSorter.h"

#include <vector>
#include <cstdint>        // uint32_t, unsigned
#include <unordered_map>  // std::unordered_map

namespace l1thgcfirmware {

  class HGCalLayer1PhiOrderFwImpl {
  public:
    HGCalLayer1PhiOrderFwImpl();
    ~HGCalLayer1PhiOrderFwImpl() {}

    void runAlgorithm() const;

    unsigned run(const l1thgcfirmware::HGCalTriggerCellSACollection& tcs_in,
                 const l1thgcfirmware::HGCalLayer1PhiOrderFwConfig& theConf,
                 l1thgcfirmware::HGCalTriggerCellSACollection& tcs_out) const;

  private:
    const int columnMask_ = 0x1F;  //
    const int columnOffset_ = 6;   //Frame and channel both have 3 bits
    const int channelMask_ = 7;    //3 bits for channel
    const int channelOffset_ = 3;  //Frame has 3 bits
    const int frameMask_ = 7;      //3 bits for frame
    const int frameOffset_ = 0;

    int packColChnFrame(int column, unsigned channel, unsigned frame) const;
    void unpackColChnFrame(int packedbin, int& column, unsigned& channel, unsigned& frame) const;

    std::vector<std::pair<l1thgcfirmware::HGCalTriggerCell, int>> assignTCToCol(
        const l1thgcfirmware::HGCalLayer1PhiOrderFwConfig& theConf,
        std::vector<l1thgcfirmware::HGCalTriggerCell> tcs) const;
    std::vector<std::pair<l1thgcfirmware::HGCalTriggerCell, int>> assignTCToChnAndFrame(
        const l1thgcfirmware::HGCalLayer1PhiOrderFwConfig& theConf,
        std::vector<std::pair<l1thgcfirmware::HGCalTriggerCell, int>> ord_tcs) const;
  };

}  // namespace l1thgcfirmware

#endif
