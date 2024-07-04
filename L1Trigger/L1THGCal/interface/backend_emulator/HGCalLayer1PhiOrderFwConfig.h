#ifndef __L1Trigger_L1THGCal_HGCalLayer1PhiOrderFwConfig_h__
#define __L1Trigger_L1THGCal_HGCalLayer1PhiOrderFwConfig_h__

#include <vector>
#include <cstdint>  // uint32_t
#include <unordered_map>

namespace l1thgcfirmware {

  class HGCalLayer1PhiOrderFwConfig {
  public:
    HGCalLayer1PhiOrderFwConfig() {}

    void setSector120(const unsigned sector) { sector120_ = sector; }
    void setFPGAID(const uint32_t fpga_id) { fpga_id_ = fpga_id; }

    void
    configureMappingInfo() {  //creates a map between column number and a vector of pairs containing the available (channel,frame) combinations. Each frame accepts 4 TCs, so each (channel,frame) combination appears 4 times in the vector. We could optionally add a unique "slot" index to this vector, but it is not strictly necessary. The function also creates a vector of pairs between column number and number of TCs the column can accept.
      for (unsigned i = 0; i < 20; ++i) {
        //std::pair<unsigned, int> mod_col_pair = std::make_pair(dummyModId_,i);//
        for (unsigned j = 0; j < 2; ++j) {
          for (unsigned k = 0; k < 4; ++k) {  //number of slots per frame
            chn_frame_slots_per_mod_and_col_[dummyModId_][i].push_back(std::make_pair(i + 3, j));
          }
          for (unsigned k = 0; k < 4; ++k) {  //number of slots per frame
            chn_frame_slots_per_mod_and_col_[dummyModId_][i].push_back(std::make_pair(i + 4, j));
          }
        }
        max_tcs_per_module_and_column_[dummyModId_].push_back(
            std::make_pair(i, chn_frame_slots_per_mod_and_col_[dummyModId_][i].size()));
      }
    }

    unsigned phiSector() const { return sector120_; }
    uint32_t fpgaID() const { return fpga_id_; }
    unsigned getColBudgetAtIndex(unsigned moduleId, unsigned theColumnIndex) const {
      return max_tcs_per_module_and_column_.at(dummyModId_).at(theColumnIndex).second;
    }  //Get TC budget for column at index theColumnIndex in the vector
    int getColFromBudgetMapAtIndex(unsigned moduleId, unsigned theColumnIndex) const {
      return max_tcs_per_module_and_column_.at(dummyModId_).at(theColumnIndex).first;
    }  //Get column number at index theColumnIndex in the vector
    unsigned getChannelAtIndex(unsigned moduleId, int theColumn, unsigned theChnFrameIndex) const {
      return chn_frame_slots_per_mod_and_col_.at(dummyModId_).at(theColumn).at(theChnFrameIndex).first;
    }  //Extract channel number for colnr theColumn, at given channel+frame index in the vector
    unsigned getFrameAtIndex(unsigned moduleId, int theColumn, unsigned theChnFrameIndex) const {
      return chn_frame_slots_per_mod_and_col_.at(dummyModId_).at(theColumn).at(theChnFrameIndex).second;
    }  //Extract frame number for colnr theColumn, at given channel+frame index in the vector

  private:
    unsigned sector120_;
    uint32_t fpga_id_;
    std::unordered_map<unsigned, std::vector<std::pair<int, unsigned>>> max_tcs_per_module_and_column_;
    std::unordered_map<unsigned, std::unordered_map<int, std::vector<std::pair<unsigned, unsigned>>>>
        chn_frame_slots_per_mod_and_col_;
    const uint32_t dummyModId_ = 1879048191;  // Just to avoid filling maps for random module ID values. Temporary!!
  };

}  // namespace l1thgcfirmware

#endif
