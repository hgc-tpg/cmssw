#include "L1Trigger/L1THGCal/interface/backend_emulator/HGCalLayer1PhiOrderFwImpl.h"
#include <cmath>
#include <algorithm>

using namespace l1thgcfirmware;

HGCalLayer1PhiOrderFwImpl::HGCalLayer1PhiOrderFwImpl() {}

struct {
  bool operator()(l1thgcfirmware::HGCalTriggerCell a, l1thgcfirmware::HGCalTriggerCell b) const {
    return a.phi() < b.phi();
  }
} sortByPhi;

unsigned HGCalLayer1PhiOrderFwImpl::run(const l1thgcfirmware::HGCalTriggerCellSACollection& tcs_in,
                                        const l1thgcfirmware::HGCalLayer1PhiOrderFwConfig& theConf,
                                        l1thgcfirmware::HGCalTriggerCellSACollection& tcs_out) const {
  std::unordered_map<unsigned, std::vector<l1thgcfirmware::HGCalTriggerCell>> tcs_per_bin;

  // group TCs per unique module
  for (const auto& tc : tcs_in) {
    unsigned module_id = tc.moduleId();
    tcs_per_bin[module_id].push_back(tc);
  }

  for (auto& bin_tcs : tcs_per_bin) {
    std::vector<l1thgcfirmware::HGCalTriggerCell> sorted_tcs;

    sorted_tcs = bin_tcs.second;

    std::vector<std::pair<l1thgcfirmware::HGCalTriggerCell, int>> tcs_per_col = assignTCToCol(theConf, sorted_tcs);
    std::vector<std::pair<l1thgcfirmware::HGCalTriggerCell, int>> tcs_with_ccf =
        assignTCToChnAndFrame(theConf, tcs_per_col);

    for (auto& tcobj : tcs_with_ccf) {
      tcs_out.push_back(tcobj.first);
      int col = 0;
      unsigned ch = 0;
      unsigned fr = 0;
      unpackColChnFrame(tcobj.second, col, ch, fr);
      tcs_out.back().setColumn(col);
      tcs_out.back().setChannel(ch);
      tcs_out.back().setFrame(fr);
    }
  }

  return 0;
}

std::vector<std::pair<l1thgcfirmware::HGCalTriggerCell, int>> HGCalLayer1PhiOrderFwImpl::assignTCToCol(
    const l1thgcfirmware::HGCalLayer1PhiOrderFwConfig& theConf,
    std::vector<l1thgcfirmware::HGCalTriggerCell> tcs) const {
  std::vector<std::pair<l1thgcfirmware::HGCalTriggerCell, int>> theOrderedTCs;
  std::sort(tcs.begin(), tcs.end(), sortByPhi);
  unsigned theColumnIndex =
      0;  //start filling the first associated column. This assumes the columns are already ordered correctly! NB this is an index, not the column ID
  unsigned nTCinColumn = 0;  //Number of TCs already in column
  for (auto& tc : tcs) {
    uint32_t theModuleId = tc.moduleId();
    while (!(nTCinColumn < theConf.getColBudgetAtIndex(theModuleId, theColumnIndex))) {
      theColumnIndex += 1;
      nTCinColumn = 0;
    }
    theOrderedTCs.push_back(std::make_pair(tc, theConf.getColFromBudgetMapAtIndex(theModuleId, theColumnIndex)));
    nTCinColumn += 1;
  }
  return theOrderedTCs;
}

std::vector<std::pair<l1thgcfirmware::HGCalTriggerCell, int>> HGCalLayer1PhiOrderFwImpl::assignTCToChnAndFrame(
    const l1thgcfirmware::HGCalLayer1PhiOrderFwConfig& theConf,
    std::vector<std::pair<l1thgcfirmware::HGCalTriggerCell, int>> ord_tcs) const {
  std::vector<std::pair<l1thgcfirmware::HGCalTriggerCell, int>> theTCsWithChnFrame;
  //std::unordered_map<int, std::pair<int,int>> chnAndFrameCounterForCol; //Need to track channel and frame counter for particular column.
  //need map from column to index for that specific column. Index should be the index in a vector that already contains the possible combinations of chn,frame (could add the slot)
  std::unordered_map<unsigned, std::unordered_map<int, unsigned>>
      chnAndFrameCounterForCol;  //This will track the index 'per column', index is index in a vector of tuples;
  for (auto& tc : ord_tcs) {
    int theCol = tc.second;
    uint32_t theModuleId = tc.first.moduleId();
    if (chnAndFrameCounterForCol[theModuleId].count(theCol) == 0)
      chnAndFrameCounterForCol[theModuleId][theCol] = 0;
    unsigned theChnFrameIndex = chnAndFrameCounterForCol[theModuleId][theCol];
    chnAndFrameCounterForCol[theModuleId][theCol] += 1;
    int thePackedCCF = packColChnFrame(theCol,
                                       theConf.getChannelAtIndex(theModuleId, theCol, theChnFrameIndex),
                                       theConf.getFrameAtIndex(theModuleId, theCol, theChnFrameIndex));
    theTCsWithChnFrame.push_back(std::make_pair(tc.first, thePackedCCF));
  }
  return theTCsWithChnFrame;
}

int HGCalLayer1PhiOrderFwImpl::packColChnFrame(int column, unsigned channel, unsigned frame) const {
  int packed_bin = 0;
  packed_bin |= (frame & frameMask_) << frameOffset_;
  packed_bin |= (channel & channelMask_) << channelOffset_;
  packed_bin |= (std::abs(column) & columnMask_) << columnOffset_;
  if (column < 0)
    packed_bin = -packed_bin;
  return packed_bin;
}

void HGCalLayer1PhiOrderFwImpl::unpackColChnFrame(int packedBin, int& column, unsigned& channel, unsigned& frame) const {
  int colsign = 1;
  if (packedBin < 0)
    colsign = -1;
  frame = (std::abs(packedBin) >> frameOffset_) & frameMask_;
  channel = (std::abs(packedBin) >> channelOffset_) & channelMask_;
  column = (std::abs(packedBin) >> columnOffset_) & columnMask_;
  column = colsign * column;
}
