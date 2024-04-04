#ifndef __L1Trigger_L1THGCal_HGCalLayer1TruncationFwImpl_h__
#define __L1Trigger_L1THGCal_HGCalLayer1TruncationFwImpl_h__

#include "L1Trigger/L1THGCal/interface/backend_emulator/HGCalTriggerCell_SA.h"
#include "L1Trigger/L1THGCal/interface/backend_emulator/HGCalLayer1TruncationFwConfig.h"
#include "L1Trigger/L1THGCal/interface/backend_emulator/BatcherSorter.h"

#include <vector>
#include <cstdint>        // uint32_t, unsigned
#include <unordered_map>  // std::unordered_map

namespace l1thgcfirmware {

  class HGCalLayer1TruncationFwImpl {
  public:
    HGCalLayer1TruncationFwImpl();
    ~HGCalLayer1TruncationFwImpl() {}

    void runAlgorithm() const;

    unsigned run(const l1thgcfirmware::HGCalTriggerCellSACollection& tcs_in,
                 const l1thgcfirmware::HGCalLayer1TruncationFwConfig& theConf,
                 l1thgcfirmware::HGCalTriggerCellSACollection& tcs_out) const;

  private:
    static constexpr unsigned offset_roz_ = 1;
    static constexpr unsigned mask_roz_ = 0x3f;  // 6 bits, max 64 bins
    static constexpr unsigned mask_phi_ = 1;

    bool do_truncate_;
    double roz_min_ = 0.;
    double roz_max_ = 0.;
    unsigned roz_bins_ = 42;
    std::vector<unsigned> max_tcs_per_bin_;
    std::vector<double> phi_edges_;

    uint32_t packBin(unsigned roverzbin, unsigned phibin) const;
    void unpackBin(unsigned packedbin, unsigned& roverzbin, unsigned& phibin) const;
    int phiBin(unsigned roverzbin, double phi, const std::vector<double>& phiedges) const;
    double rotatedphi(double x, double y, double z, unsigned sector) const;
    double rotatedphi(double phi, unsigned sector) const;

    unsigned smallerMultOfFourGreaterThan(unsigned n) const;
  };

}  // namespace l1thgcfirmware

#endif
