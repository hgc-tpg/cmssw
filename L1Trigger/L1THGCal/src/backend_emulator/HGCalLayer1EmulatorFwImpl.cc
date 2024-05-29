#include "L1Trigger/L1THGCal/interface/backend_emulator/HGCalLayer1EmulatorFwImpl.h"
#include <cmath>
#include <algorithm>

using namespace l1thgcfirmware;

HGCalLayer1EmulatorFwImpl::HGCalLayer1EmulatorFwImpl() {}

struct {
   bool operator()(l1thgcfirmware::HGCalTriggerCell a, l1thgcfirmware::HGCalTriggerCell b) const { return a.phi() < b.phi(); }
} sortByPhi;

unsigned HGCalLayer1EmulatorFwImpl::run(const l1thgcfirmware::HGCalTriggerCellSACollection& tcs_in,
                                          const l1thgcfirmware::HGCalLayer1TruncationFwConfig& theConf,
                                          l1thgcfirmware::HGCalTriggerCellSACollection& tcs_out) const {
  std::unordered_map<unsigned, std::vector<l1thgcfirmware::HGCalTriggerCell>> tcs_per_bin;

  // configuation:
  //bool do_truncate = theConf.doTruncate();
  //const std::vector<unsigned>& maxtcsperbin = theConf.maxTcsPerBin();

  //const std::vector<std::pair<int,int>> cols_budget = theConf.maxTcsPerColumn();
  //const std::unordered_map<int,std::vector<std::pair<int,int>>> chns_frs_percol = theConf.channelsAndFramesPerColumn();

  // group TCs per unique module
  for (const auto& tc : tcs_in) {
    unsigned module_id = tc.moduleId(); 
    //unsigned roverzbin = tc.rOverZ();
    //int phibin = tc.phi();
    //if (phibin < 0)
    //  return 1;
    //unsigned packed_bin = packBin(roverzbin, phibin);

    tcs_per_bin[module_id].push_back(tc);
  }

  for (auto& bin_tcs : tcs_per_bin) {

    unsigned module_id = bin_tcs.first;

    std::vector<l1thgcfirmware::HGCalTriggerCell> sorted_tcs;
    
    sorted_tcs = bin_tcs.second;

    std::vector<std::pair<l1thgcfirmware::HGCalTriggerCell,int>> tcs_per_col = assignTCToCol(theConf,sorted_tcs);
    //std::vector<std::pair<l1thgcfirmware::HGCalTriggerCell,int>> tcs_per_col = assignTCToCol(cols_budget,sorted_tcs);
    std::vector<std::pair<l1thgcfirmware::HGCalTriggerCell,int>> tcs_with_ccf = assignTCToChnAndFrame(theConf,tcs_per_col);

    for (auto & tcobj : tcs_with_ccf){
        tcs_out.push_back(tcobj.first);
        int col;
        int ch;
        int fr;
        unpackColChnFrame(tcobj.second,col,ch,fr);
        tcs_out.back().setColumn(col);
        tcs_out.back().setChannel(ch);
        tcs_out.back().setFrame(fr);
    }

  }

  return 0;
}

/*std::vector<std::pair<l1thgcfirmware::HGCalTriggerCell,int>> HGCalLayer1EmulatorFwImpl::assignTCToCol(std::vector<std::pair<int,int>> theCols, std::vector<l1thgcfirmware::HGCalTriggerCell> tcs) const {
  std::vector<std::pair<l1thgcfirmware::HGCalTriggerCell,int>> theOrderedTCs;
  std::sort(tcs.begin(), tcs.end(),sortByPhi);
  int theColumnIndex=0;//start filling the first associated column. This assumes the columns are already ordered correctly!
  int nTCinColumn=0;//Number of TCs already in column
  for(auto & tc: tcs){
    while(!(nTCinColumn < theCols.at(theColumnIndex).second)){
      theColumnIndex+=1;
      nTCinColumn=0;
    }
    theOrderedTCs.push_back(std::make_pair(tc,theCols.at(theColumnIndex).first));
    nTCinColumn+=1;
  }
  return theOrderedTCs;
}*/

std::vector<std::pair<l1thgcfirmware::HGCalTriggerCell,int>> HGCalLayer1EmulatorFwImpl::assignTCToCol(const l1thgcfirmware::HGCalLayer1TruncationFwConfig& theConf, std::vector<l1thgcfirmware::HGCalTriggerCell> tcs) const {
  std::vector<std::pair<l1thgcfirmware::HGCalTriggerCell,int>> theOrderedTCs;
  std::sort(tcs.begin(), tcs.end(),sortByPhi);
  int theColumnIndex=0;//start filling the first associated column. This assumes the columns are already ordered correctly!
  int nTCinColumn=0;//Number of TCs already in column
  for(auto & tc: tcs){
    while(!(nTCinColumn < theConf.getColBudgetAtIndex(theColumnIndex))){
   // while(!(nTCinColumn < theCols.at(theColumnIndex).second)){
      theColumnIndex+=1;
      nTCinColumn=0;
    }
    theOrderedTCs.push_back(std::make_pair(tc,theConf.getColFromBudgetMapAtIndex(theColumnIndex)));
    //theOrderedTCs.push_back(std::make_pair(tc,theCols.at(theColumnIndex).first));
    nTCinColumn+=1;
  }
  return theOrderedTCs;
}


std::vector<std::pair<l1thgcfirmware::HGCalTriggerCell,int>> HGCalLayer1EmulatorFwImpl::assignTCToChnAndFrame(const l1thgcfirmware::HGCalLayer1TruncationFwConfig& theConf, std::vector<std::pair<l1thgcfirmware::HGCalTriggerCell,int>> ord_tcs) const {
  std::vector<std::pair<l1thgcfirmware::HGCalTriggerCell,int>> theTCsWithChnFrame;
  //std::unordered_map<int, std::pair<int,int>> chnAndFrameCounterForCol; //Need to track channel and frame counter for particular column.
  //need map from column to index for that specific column. Index should be the index in a vector that already contains the possible combinations of chn,frame (posibly slot)
  std::unordered_map<int, int> chnAndFrameCounterForCol;//This will track the index 'per column', index is index in a vector of tuples;
  for (auto & tc: ord_tcs){
    int theCol = tc.second;
    if (chnAndFrameCounterForCol.count(theCol) == 0) chnAndFrameCounterForCol[theCol] = 0;
    int theChnFrameIndex = chnAndFrameCounterForCol[theCol];
    chnAndFrameCounterForCol[theCol]+=1;
    int thePackedCCF = packColChnFrame(theCol,theConf.getChannelAtIndex(theCol,theChnFrameIndex),theConf.getFrameAtIndex(theCol,theChnFrameIndex));
    theTCsWithChnFrame.push_back(std::make_pair(tc.first,thePackedCCF));
  }
  return theTCsWithChnFrame;
}


/*std::vector<std::pair<l1thgcfirmware::HGCalTriggerCell,int>> HGCalLayer1EmulatorFwImpl::assignTCToChnAndFrame(std::unordered_map<int,std::vector<std::pair<int,int>>> chnsAndFrames, std::vector<std::pair<l1thgcfirmware::HGCalTriggerCell,int>> ord_tcs) const {
  std::vector<std::pair<l1thgcfirmware::HGCalTriggerCell,int>> theTCsWithChnFrame;
  //std::unordered_map<int, std::pair<int,int>> chnAndFrameCounterForCol; //Need to track channel and frame counter for particular column.
  //need map from column to index for that specific column. Index should be the index in a vector that already contains the possible combinations of chn,frame (posibly slot)
  std::unordered_map<int, int> chnAndFrameCounterForCol;//This will track the index 'per column', index is index in a vector of tuples;
  for (auto & tc: ord_tcs){
    int theCol = tc.second;
    if (chnAndFrameCounterForCol.count(theCol) == 0) chnAndFrameCounterForCol[theCol] = 0;
    int theChnFrameIndex = chnAndFrameCounterForCol[theCol];
    chnAndFrameCounterForCol[theCol]+=1;
    int thePackedCCF = packColChnFrame(theCol,chnsAndFrames[theCol].at(theChnFrameIndex).first,chnsAndFrames[theCol].at(theChnFrameIndex).second);
    theTCsWithChnFrame.push_back(std::make_pair(tc.first,thePackedCCF));
  }
  return theTCsWithChnFrame;
}*/


int HGCalLayer1EmulatorFwImpl::packColChnFrame(int column, int channel, int frame) const {//temporary very dumb 3 int -> 1 int conversion, can make this much smarter
  int packed_bin = 0;
  packed_bin = column*100000+channel*1000+frame;
  return packed_bin;
}

void HGCalLayer1EmulatorFwImpl::unpackColChnFrame(int packedbin, int& column, int& channel, int& frame) const {
  column = packedbin/100000;
  channel = (packedbin-column*100000)/1000;
  frame = (packedbin-column*100000-channel*1000);
}

uint32_t HGCalLayer1EmulatorFwImpl::packBin(unsigned roverzbin, unsigned phibin) const {
  unsigned packed_bin = 0;
  packed_bin |= ((roverzbin & mask_roz_) << offset_roz_);
  packed_bin |= (phibin & mask_phi_);
  return packed_bin;
}

void HGCalLayer1EmulatorFwImpl::unpackBin(unsigned packedbin, unsigned& roverzbin, unsigned& phibin) const {
  roverzbin = ((packedbin >> offset_roz_) & mask_roz_);
  phibin = (packedbin & mask_phi_);
}

int HGCalLayer1EmulatorFwImpl::phiBin(unsigned roverzbin, double phi, const std::vector<double>& phiedges) const {
  unsigned phi_bin = 0;
  if (roverzbin >= phiedges.size())
    return -1;
  double phi_edge = phiedges[roverzbin];
  if (phi > phi_edge)
    phi_bin = 1;
  return phi_bin;
}

double HGCalLayer1EmulatorFwImpl::rotatedphi(double phi, unsigned sector) const {
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

double HGCalLayer1EmulatorFwImpl::rotatedphi(double x, double y, double z, unsigned sector) const {
  if (z > 0)
    x = -x;
  double phi = std::atan2(y, x);
  return this->rotatedphi(phi, sector);
}

unsigned HGCalLayer1EmulatorFwImpl::rozBin(double roverz, double rozmin, double rozmax, unsigned rozbins) const {
  constexpr double margin = 1.001;
  double roz_bin_size = (rozbins > 0 ? (rozmax - rozmin) * margin / double(rozbins) : 0.);
  unsigned roverzbin = 0;
  if (roz_bin_size > 0.) {
    roverz -= rozmin;
    roverz = std::clamp(roverz, 0., rozmax - rozmin);
    roverzbin = unsigned(roverz / roz_bin_size);
  }

  return roverzbin;
}

unsigned HGCalLayer1EmulatorFwImpl::smallerMultOfFourGreaterThan(unsigned N) const {
  unsigned remnant = (N + 4) % 4;
  if (remnant == 0)
    return N;
  else
    return (N + 4 - remnant);
}
