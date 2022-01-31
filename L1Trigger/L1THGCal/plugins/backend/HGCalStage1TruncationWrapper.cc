#include "L1Trigger/L1THGCal/interface/HGCalAlgoWrapperBase.h"

#include "DataFormats/L1THGCal/interface/HGCalTriggerCell.h"
#include "L1Trigger/L1THGCal/interface/backend_emulator/HGCalTriggerCell_SA.h"
#include "L1Trigger/L1THGCal/interface/backend_emulator/HGCalStage1TruncationImpl_SA.h"
#include "L1Trigger/L1THGCal/interface/backend_emulator/HGCalStage1TruncationConfig_SA.h"

#include "L1Trigger/L1THGCal/interface/HGCalTriggerTools.h"

class HGCalStage1TruncationWrapper : public HGCalStage1TruncationWrapperBase {
public:
  HGCalStage1TruncationWrapper(const edm::ParameterSet& conf);
  ~HGCalStage1TruncationWrapper() override = default;

  void configure(
      const std::tuple<const HGCalTriggerGeometryBase* const, const unsigned&, const uint32_t&>& configuration) override;

  void process(const std::vector<edm::Ptr<l1t::HGCalTriggerCell>>& fpga_tcs,
               std::vector<edm::Ptr<l1t::HGCalTriggerCell>>& tcs_out) const override;

private:
  void convertCMSSWInputs(const std::vector<edm::Ptr<l1t::HGCalTriggerCell>>& fpga_tcs,
                          l1thgcfirmware::HGCalTriggerCellSACollection& fpga_tcs_SA) const;

  void convertAlgorithmOutputs(const l1thgcfirmware::HGCalTriggerCellSACollection& fpga_tcs_out,
                               const std::vector<edm::Ptr<l1t::HGCalTriggerCell>>& fpga_tcs_original,
                               std::vector<edm::Ptr<l1t::HGCalTriggerCell>>& fpga_tcs_trunc) const;

  void setGeometry(const HGCalTriggerGeometryBase* const geom) { triggerTools_.setGeometry(geom); }

  double rotatedphi(double phi, unsigned sector) const;

  HGCalTriggerTools triggerTools_;
  l1thgcfirmware::HGCalStage1TruncationImplSA theAlgo_;
  l1thgcfirmware::Stage1TruncationConfig theConfiguration_;
};

HGCalStage1TruncationWrapper::HGCalStage1TruncationWrapper(const edm::ParameterSet& conf)
    : HGCalStage1TruncationWrapperBase(conf),
      theAlgo_(),
      theConfiguration_(conf.getParameter<bool>("doTruncation"),
                        conf.getParameter<double>("rozMin"),
                        conf.getParameter<double>("rozMax"),
                        conf.getParameter<unsigned>("rozBins"),
                        conf.getParameter<std::vector<unsigned>>("maxTcsPerBin"),
                        conf.getParameter<std::vector<double>>("phiSectorEdges")) {}

void HGCalStage1TruncationWrapper::convertCMSSWInputs(const std::vector<edm::Ptr<l1t::HGCalTriggerCell>>& fpga_tcs,
                                                      l1thgcfirmware::HGCalTriggerCellSACollection& fpga_tcs_SA) const {
  fpga_tcs_SA.clear();
  fpga_tcs_SA.reserve(fpga_tcs.size());
  unsigned int itc = 0;
  for (auto& tc : fpga_tcs) {
    //    fpga_tcs_SA.emplace_back(tc->position().x(),
    //                             tc->position().y(),
    //                             tc->position().z(),
    //                             triggerTools_.zside(tc->detId()),
    //                             triggerTools_.layerWithOffset(tc->detId()),
    //                             tc->eta(),
    //                             tc->phi(),
    //                             tc->pt(),
    //                             tc->mipPt(),
    //                             itc);
    const GlobalPoint& position = tc->position();
    double x = position.x();
    double y = position.y();
    double z = position.z();
    unsigned int digi_rOverZ = (std::sqrt(x * x + y * y) / std::abs(z)) * 10000;  // Magic numbers
    double phi = rotatedphi(tc->phi(), theConfiguration_.phiSector());
    phi += (phi < 0) ? M_PI : 0;
    unsigned int digi_phi = phi * 10000 /* * 1944 / M_PI*/;  // Magic numbers
    unsigned int digi_energy = (tc->mipPt()) * 10000;        // Magic numbers
    fpga_tcs_SA.emplace_back(
        true, true, digi_rOverZ, digi_phi, triggerTools_.layerWithOffset(tc->detId()), digi_energy, itc);
    ++itc;
  }
}

void HGCalStage1TruncationWrapper::convertAlgorithmOutputs(
    const l1thgcfirmware::HGCalTriggerCellSACollection& fpga_tcs_out,
    const std::vector<edm::Ptr<l1t::HGCalTriggerCell>>& fpga_tcs_original,
    std::vector<edm::Ptr<l1t::HGCalTriggerCell>>& fpga_tcs_trunc) const {
  for (auto& tc : fpga_tcs_out) {
    unsigned tc_cmssw_id = tc.index_cmssw();
    fpga_tcs_trunc.push_back(fpga_tcs_original[tc_cmssw_id]);
  }
}

void HGCalStage1TruncationWrapper::process(const std::vector<edm::Ptr<l1t::HGCalTriggerCell>>& fpga_tcs,
                                           std::vector<edm::Ptr<l1t::HGCalTriggerCell>>& tcs_out) const {
  l1thgcfirmware::HGCalTriggerCellSACollection fpga_tcs_SA;
  convertCMSSWInputs(fpga_tcs, fpga_tcs_SA);

  l1thgcfirmware::HGCalTriggerCellSACollection tcs_out_SA;
  unsigned error_code = theAlgo_.run(fpga_tcs_SA, theConfiguration_, tcs_out_SA);

  if (error_code == 1)
    throw cms::Exception("HGCalStage1TruncationImpl::OutOfRange") << "roverzbin index out of range";

  convertAlgorithmOutputs(tcs_out_SA, fpga_tcs, tcs_out);
}

void HGCalStage1TruncationWrapper::configure(
    const std::tuple<const HGCalTriggerGeometryBase* const, const unsigned&, const uint32_t&>& configuration) {
  setGeometry(std::get<0>(configuration));

  theConfiguration_.setSector120(std::get<1>(configuration));
  theConfiguration_.setFPGAID(std::get<2>(configuration));
};

double HGCalStage1TruncationWrapper::rotatedphi(double phi, unsigned sector) const {
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

DEFINE_EDM_PLUGIN(HGCalStage1TruncationWrapperBaseFactory,
                  HGCalStage1TruncationWrapper,
                  "HGCalStage1TruncationWrapper");
