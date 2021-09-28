#ifndef __L1Trigger_L1THGCal_DistServer_h__
#define __L1Trigger_L1THGCal_DistServer_h__

#include "L1Trigger/L1THGCal/interface/backend_emulator/HGCalTriggerCell_SA.h"

#include <vector>

namespace l1thgcfirmware {

    class DistServer {
    public:
        DistServer( unsigned int nInputs, unsigned int nOutputs, unsigned int nInterleaving );
        ~DistServer() {}

        HGCalTriggerCellSAPtrCollection clock(HGCalTriggerCellSAPtrCollection& inputs);

        unsigned int nInputs() const { return nInputs_; }
        unsigned int nOutputs() const { return nOutputs_; }
        unsigned int nInterleaving() const { return nInterleaving_; }
        std::vector< std::vector< unsigned int> >& addr() { return addr_; }
        l1thgcfirmware::HGCalTriggerCellSAPtrCollections& inputs() { return inputs_; }

    private:
        unsigned int nInputs_;
        unsigned int nOutputs_;
        unsigned int nInterleaving_;

        l1thgcfirmware::HGCalTriggerCellSAPtrCollections inputs_;
        std::vector< std::vector< unsigned int> > addr_;        
    };
}

#endif
