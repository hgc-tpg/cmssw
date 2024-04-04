#include "L1Trigger/L1THGCal/interface/backend_emulator/HGCalLayer1TruncationFwImpl.h"
#include <cmath>
#include <iostream>

using namespace l1thgcfirmware;

HGCalLayer1TruncationFwImpl::HGCalLayer1TruncationFwImpl() {}

unsigned HGCalLayer1TruncationFwImpl::run(const l1thgcfirmware::HGCalTriggerCellSACollection& tcs_in,
                                          const l1thgcfirmware::HGCalLayer1TruncationFwConfig& theConf,
                                          l1thgcfirmware::HGCalTriggerCellSACollection& tcs_out) const {
  std::unordered_map<unsigned, std::vector<l1thgcfirmware::HGCalTriggerCell>> tcs_per_bin;

  // configuation:
  bool do_truncate = theConf.doTruncate();
  const std::vector<unsigned>& maxtcsperbin = theConf.maxTcsPerBin();

  // group TCs per (r/z, phi) bins
  for (const auto& tc : tcs_in) {
    unsigned roverzbin = tc.rOverZ();
    int phibin = tc.phi();
    if (phibin < 0)
      return 1;
    unsigned packed_bin = packBin(roverzbin, phibin);

    tcs_per_bin[packed_bin].push_back(tc);
  }

  for (auto& bin_tcs : tcs_per_bin) {
    unsigned roverzbin = 0;
    unsigned phibin = 0;
    unpackBin(bin_tcs.first, roverzbin, phibin);

    const unsigned ntcin = std::max(smallerMultOfFourGreaterThan(bin_tcs.second.size()),
                                    smallerMultOfFourGreaterThan(maxtcsperbin[roverzbin]));
    const unsigned ntcout = (do_truncate ? maxtcsperbin[roverzbin] : bin_tcs.second.size());

    l1thgcfirmware::BatcherSorter tcSorter(ntcin, ntcout);

    std::vector<unsigned> theTCsIn_mipt(ntcin);
    std::vector<unsigned> theTCsOut_mipt(ntcout);
    std::vector<unsigned> theTCsOut_addr(ntcout);

    for (unsigned i = 0; i < ntcin; ++i) {
      if (i < bin_tcs.second.size())
        theTCsIn_mipt[i] = bin_tcs.second.at(i).energy();
      else
        theTCsIn_mipt[i] = 0;
    }

    tcSorter.sorting(theTCsIn_mipt, theTCsOut_mipt, theTCsOut_addr);

    for (const unsigned& tcid : theTCsOut_addr) {
      if (tcid < bin_tcs.second.size()) {
        tcs_out.push_back(bin_tcs.second.at(tcid));
      } else {  // Create and push a dummy TC with 0 energy and identifiable address
        l1thgcfirmware::HGCalTriggerCell dummytc(false, false, roverzbin, phibin, 0, 0);
        dummytc.setIndex(255);  // set moduleID to dummy value from FW
        dummytc.setCmsswIndex(
            std::make_pair(std::numeric_limits<unsigned int>::max(),
                           63));  // set TC ID to dummy value to be read in CMSSW (1) and in emulator (2)
        tcs_out.push_back(dummytc);
      }
    }
  }

  return 0;
}

uint32_t HGCalLayer1TruncationFwImpl::packBin(unsigned roverzbin, unsigned phibin) const {
  unsigned packed_bin = 0;
  packed_bin |= ((roverzbin & mask_roz_) << offset_roz_);
  packed_bin |= (phibin & mask_phi_);
  return packed_bin;
}

void HGCalLayer1TruncationFwImpl::unpackBin(unsigned packedbin, unsigned& roverzbin, unsigned& phibin) const {
  roverzbin = ((packedbin >> offset_roz_) & mask_roz_);
  phibin = (packedbin & mask_phi_);
}

int HGCalLayer1TruncationFwImpl::phiBin(unsigned roverzbin, double phi, const std::vector<double>& phiedges) const {
  unsigned phi_bin = 0;
  if (roverzbin >= phiedges.size())
    return -1;
  double phi_edge = phiedges[roverzbin];
  if (phi > phi_edge)
    phi_bin = 1;
  return phi_bin;
}

double HGCalLayer1TruncationFwImpl::rotatedphi(double phi, unsigned sector) const {
  if (sector == 1) {
    if (phi < M_PI and phi > 0)
      phi = phi - (2. * M_PI / 3.);
    else
      phi = phi + (4. * M_PI / 3.);
  } else if (sector == 2) {
    phi = phi + (2. * M_PI / 3.);
  }
  return phi;
}

unsigned HGCalLayer1TruncationFwImpl::smallerMultOfFourGreaterThan(unsigned N) const {
  unsigned remnant = (N + 4) % 4;
  if (remnant == 0)
    return N;
  else
    return (N + 4 - remnant);
}
