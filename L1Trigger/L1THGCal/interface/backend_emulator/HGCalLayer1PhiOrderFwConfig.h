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

    void configureMappingInfo(){
      for(unsigned i=0; i<20; ++i){
        for(unsigned j=0;j<2;++j){
          for(unsigned k=0; k<4; ++k){//number of slots per frame
              chn_frame_slots_per_col_[i].push_back(std::make_pair(i+3,j));
          }
          for(unsigned k=0; k<4; ++k){//number of slots per frame
              chn_frame_slots_per_col_[i].push_back(std::make_pair(i+4,j));
          }
        }  
        max_tcs_per_column_.push_back(std::make_pair(i,chn_frame_slots_per_col_[i].size()));
      }
    }

    unsigned phiSector() const { return sector120_; }
    uint32_t fpgaID() const { return fpga_id_; }
    unsigned getColBudgetAtIndex(unsigned theColumnIndex ) const { return max_tcs_per_column_.at(theColumnIndex).second; }
    int getColFromBudgetMapAtIndex(unsigned theColumnIndex ) const  {return max_tcs_per_column_.at(theColumnIndex).first; }
    unsigned getChannelAtIndex(unsigned theColumnIndex, unsigned theChnFrameIndex ) const { return chn_frame_slots_per_col_.at(theColumnIndex).at(theChnFrameIndex).first; }
    unsigned getFrameAtIndex(unsigned theColumnIndex, unsigned theChnFrameIndex ) const {return chn_frame_slots_per_col_.at(theColumnIndex).at(theChnFrameIndex).second; }


  private:
    unsigned sector120_;
    uint32_t fpga_id_;
    std::vector<std::pair<int,unsigned>> max_tcs_per_column_;
    std::unordered_map<int,std::vector<std::pair<unsigned,unsigned>>> chn_frame_slots_per_col_; 
  };

}  // namespace l1thgcfirmware

#endif
