#include "L1Trigger/L1THGCal/interface/HGCalAlgoWrapperBase.h"

#include "DataFormats/L1THGCal/interface/HGCalTriggerCell.h"
#include "L1Trigger/L1THGCal/interface/backend_emulator/HGCalTriggerCell_SA.h"
#include "L1Trigger/L1THGCal/interface/backend_emulator/HGCalLayer1PhiOrderFwImpl.h"
#include "L1Trigger/L1THGCal/interface/backend_emulator/HGCalLayer1PhiOrderFwConfig.h"

#include "L1Trigger/L1THGCal/interface/HGCalTriggerTools.h"

class HGCalLayer1PhiOrderWrapper : public HGCalLayer1PhiOrderWrapperBase {
public:
  HGCalLayer1PhiOrderWrapper(const edm::ParameterSet& conf);
  ~HGCalLayer1PhiOrderWrapper() override = default;

  void configure(
      const std::tuple<const HGCalTriggerGeometryBase* const, const unsigned&, const uint32_t&>& configuration) override;

  void process(const std::vector<edm::Ptr<l1t::HGCalTriggerCell>>& fpga_tcs,
               l1t::HGCalClusterBxCollection& clusters) const override;

private:
  void convertCMSSWInputs(const std::vector<edm::Ptr<l1t::HGCalTriggerCell>>& fpga_tcs,
                          l1thgcfirmware::HGCalTriggerCellSACollection& fpga_tcs_SA) const;

  void convertAlgorithmOutputs(const l1thgcfirmware::HGCalTriggerCellSACollection& fpga_tcs_out,
                               const std::vector<edm::Ptr<l1t::HGCalTriggerCell>>& fpga_tcs_original,
                               l1t::HGCalClusterBxCollection& clusters) const;

  void setGeometry(const HGCalTriggerGeometryBase* const geom) { triggerTools_.setGeometry(geom); }

  double rotatedphi(double phi, unsigned sector) const;

  HGCalTriggerTools triggerTools_;
  l1thgcfirmware::HGCalLayer1PhiOrderFwImpl theAlgo_;
  l1thgcfirmware::HGCalLayer1PhiOrderFwConfig theConfiguration_;

  // Scale factor for quantities sent to emulator to keep floating point precision. Value is arbitrary and could be set to relevant value.
  const unsigned int FWfactor_ = 10000;
};

HGCalLayer1PhiOrderWrapper::HGCalLayer1PhiOrderWrapper(const edm::ParameterSet& conf)
    : HGCalLayer1PhiOrderWrapperBase(conf),
      theAlgo_(),
      theConfiguration_() {}

void HGCalLayer1PhiOrderWrapper::convertCMSSWInputs(const std::vector<edm::Ptr<l1t::HGCalTriggerCell>>& fpga_tcs,
                                                      l1thgcfirmware::HGCalTriggerCellSACollection& fpga_tcs_SA) const {
  fpga_tcs_SA.clear();
  fpga_tcs_SA.reserve(fpga_tcs.size());
  unsigned int itc = 0;
  for (auto& tc : fpga_tcs) {
    double phi = rotatedphi(tc->phi(), theConfiguration_.phiSector());
    phi += (phi < 0) ? M_PI : 0;
    unsigned int digi_phi = phi * FWfactor_;
    unsigned int digi_energy = (tc->mipPt()) * FWfactor_;
    fpga_tcs_SA.emplace_back(true, true, 1, digi_phi, triggerTools_.layerWithOffset(tc->detId()), digi_energy); //Assign dummy RoverZbin to not change the dataformat/creators. 
    fpga_tcs_SA.back().setModuleId(triggerTools_.getModuleId(tc->detId()));
    fpga_tcs_SA.back().setCmsswIndex(std::make_pair(itc, 0));
    ++itc;
  }
}

void HGCalLayer1PhiOrderWrapper::convertAlgorithmOutputs(
    const l1thgcfirmware::HGCalTriggerCellSACollection& fpga_tcs_out,
    const std::vector<edm::Ptr<l1t::HGCalTriggerCell>>& fpga_tcs_original,
    l1t::HGCalClusterBxCollection& clusters) const {

  std::vector<l1t::HGCalCluster> clustersTmp;

  for (std::vector<edm::Ptr<l1t::HGCalTriggerCell>>::const_iterator otc = fpga_tcs_original.begin(); otc!=fpga_tcs_original.end(); ++otc) {
    clustersTmp.emplace_back(*otc);
  }


  clusters.resize(0,clustersTmp.size());
  for (auto& tc : fpga_tcs_out) {
    unsigned tc_cmssw_id = tc.cmsswIndex().first;
    if (tc_cmssw_id < fpga_tcs_original.size()){
      clustersTmp.at(tc_cmssw_id).setColumn(tc.column());
      clustersTmp.at(tc_cmssw_id).setFrame(tc.frame());
      clustersTmp.at(tc_cmssw_id).setChannel(tc.channel());
      clustersTmp.at(tc_cmssw_id).setModule(tc.moduleId());
      clusters.set(0,tc_cmssw_id,clustersTmp.at(tc_cmssw_id));
    }
  }
}

void HGCalLayer1PhiOrderWrapper::process(const std::vector<edm::Ptr<l1t::HGCalTriggerCell>>& fpga_tcs,
                                           l1t::HGCalClusterBxCollection& clusters) const {
  l1thgcfirmware::HGCalTriggerCellSACollection fpga_tcs_SA;
  convertCMSSWInputs(fpga_tcs, fpga_tcs_SA);

  l1thgcfirmware::HGCalTriggerCellSACollection tcs_out_SA;
  unsigned error_code = theAlgo_.run(fpga_tcs_SA, theConfiguration_, tcs_out_SA);

  if (error_code == 1)
    throw cms::Exception("HGCalLayer1PhiOrderFwImpl::Error") << " This shouldn't happen ";

  convertAlgorithmOutputs(tcs_out_SA, fpga_tcs, clusters);
}

void HGCalLayer1PhiOrderWrapper::configure(
    const std::tuple<const HGCalTriggerGeometryBase* const, const unsigned&, const uint32_t&>& configuration) {
  setGeometry(std::get<0>(configuration));

  theConfiguration_.setSector120(std::get<1>(configuration));
  theConfiguration_.setFPGAID(std::get<2>(configuration));
  theConfiguration_.configureMappingInfo();
};

double HGCalLayer1PhiOrderWrapper::rotatedphi(double phi, unsigned sector) const {
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

DEFINE_EDM_PLUGIN(HGCalLayer1PhiOrderWrapperBaseFactory,
                  HGCalLayer1PhiOrderWrapper,
                  "HGCalLayer1PhiOrderWrapper");
