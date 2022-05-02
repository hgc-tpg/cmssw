#include "L1Trigger/L1THGCal/interface/HGCalTriggerCellCalibration.h"

#include <cmath>

HGCalTriggerCellCalibration::HGCalTriggerCellCalibration(const edm::ParameterSet& conf)
    : lsb_(conf.getParameter<double>("lsb")),
      fCperMIP_(conf.getParameter<std::vector<double>>("fCperMIP")),
      chargeCollectionEfficiency_(conf.getParameter<edm::ParameterSet>("chargeCollectionEfficiency")
                                      .getParameter<std::vector<double>>("values")),
      thicknessCorrection_(conf.getParameter<std::vector<double>>("thicknessCorrection")),
      dEdX_weights_(conf.getParameter<std::vector<double>>("dEdXweights")),
      new_digi_(conf.getParameter<bool>("newDigi")) {
  for (const auto& fCperMIP : fCperMIP_) {
    if (fCperMIP <= 0) {
      edm::LogWarning("DivisionByZero") << "WARNING: zero or negative MIP->fC correction factor. It won't be "
                                           "applied to correct trigger cell energies.";
    }
  }
  for (const auto& cce : chargeCollectionEfficiency_) {
    if (cce <= 0) {
      edm::LogWarning("DivisionByZero") << "WARNING: zero or negative cell-thickness correction factor. It won't be "
                                           "applied to correct trigger cell energies.";
    }
  }
  for (const auto& thickCorr : thicknessCorrection_) {
    if (thickCorr <= 0) {
      edm::LogWarning("DivisionByZero") << "WARNING: zero or negative cell-thickness correction factor. It won't be "
                                           "applied to correct trigger cell energies.";
    }
  }

  if (new_digi_) {
    noise_map_.setDoseMap(conf.getParameter<edm::ParameterSet>("noise").getParameter<std::string>("doseMap"),
        conf.getParameter<edm::ParameterSet>("noise").getParameter<uint32_t>("scaleByDoseAlgo"));
    noise_map_.setFluenceScaleFactor(conf.getParameter<edm::ParameterSet>("noise").getParameter<double>("scaleByDoseFactor"));
    noise_map_.setIleakParam(
        conf.getParameter<edm::ParameterSet>("ileakParam").getParameter<std::vector<double>>("ileakParam"));
    noise_map_.setCceParam(
        conf.getParameter<edm::ParameterSet>("cceParams").getParameter<std::vector<double>>("cceParamFine"),
        conf.getParameter<edm::ParameterSet>("cceParams").getParameter<std::vector<double>>("cceParamThin"),
        conf.getParameter<edm::ParameterSet>("cceParams").getParameter<std::vector<double>>("cceParamThick"));
  }
}

void HGCalTriggerCellCalibration::setGeometry(const HGCalTriggerGeometryBase* const geom, DetId::Detector det) {
  triggerTools_.setGeometry(geom);
  if (new_digi_) {
    //assign the geometry and tell the tool that the gain is automatically set to have the MIP close to 10ADC counts
    switch (det) {
      case DetId::HGCalEE:
        noise_map_.setGeometry(
            triggerTools_.getTriggerGeometry()->eeGeometry(), HGCalSiNoiseMap<HGCSiliconDetId>::AUTO, 10);
        break;
      case DetId::HGCalHSi:
        noise_map_.setGeometry(
            triggerTools_.getTriggerGeometry()->hsiGeometry(), HGCalSiNoiseMap<HGCSiliconDetId>::AUTO, 10);
        break;
      default:
        throw cms::Exception("SetupError") << "Non supported detector type " << det << " for HGCalSiNoiseMap setup";
    }
  }
}

void HGCalTriggerCellCalibration::calibrateInMipT(l1t::HGCalTriggerCell& trgCell) const {
  DetId trgdetid(trgCell.detId());
  bool isSilicon = triggerTools_.isSilicon(trgdetid);
  constexpr int kScintillatorIndex = 0;
  unsigned thickness = isSilicon ? triggerTools_.thicknessIndex(trgdetid) : kScintillatorIndex;
  if (thickness >= fCperMIP_.size()) {
    throw cms::Exception("OutOfBound") << "Trying to access thickness index " << thickness
                                       << " in fCperMIP, which is of size " << fCperMIP_.size();
  }
  if (thickness >= chargeCollectionEfficiency_.size()) {
    throw cms::Exception("OutOfBound") << "Trying to access thickness index " << thickness
                                       << " in chargeCollectionEfficiency, which is of size "
                                       << chargeCollectionEfficiency_.size();
  }

  /* get the hardware pT in ADC counts: */
  int hwPt = trgCell.hwPt();

  // Convert ADC to charge in fC (in EE+FH) or in MIPs (in BH)
  double amplitude = hwPt * lsb_;
  double trgCellMipP = amplitude;

  if (new_digi_) {
    // double mipfC = 0.;
    double cce = 0.;
    auto cells = triggerTools_.getTriggerGeometry()->getCellsFromTriggerCell(trgdetid);
    for (const auto& cellid : cells) {
      HGCalSiNoiseMap<HGCSiliconDetId>::SiCellOpCharacteristics siop = noise_map_.getSiCellOpCharacteristics(cellid);
      // mipfC += double(siop.mipfC);
      cce += siop.core.cce;
    }
    // mipfC /= cells.size();
    cce /= cells.size();
    trgCellMipP /= cce;
    if (fCperMIP_[thickness] > 0) {
      trgCellMipP /= fCperMIP_[thickness];
    }
    // trgCellMipP /= mipfC;
  } else {
    if (chargeCollectionEfficiency_[thickness] > 0) {
      trgCellMipP /= chargeCollectionEfficiency_[thickness];
    }
    if (fCperMIP_[thickness] > 0) {
      trgCellMipP /= fCperMIP_[thickness];
    }
  }

  /* compute the transverse-mip */
  double trgCellMipPt = trgCellMipP / std::cosh(trgCell.eta());

  /* setting pT [mip] */
  trgCell.setMipPt(trgCellMipPt);
}

void HGCalTriggerCellCalibration::calibrateMipTinGeV(l1t::HGCalTriggerCell& trgCell) const {
  constexpr double MevToGeV(0.001);
  double trgCellEt(0.);

  DetId trgdetid(trgCell.detId());
  unsigned trgCellLayer = triggerTools_.layerWithOffset(trgdetid);
  bool isSilicon = triggerTools_.isSilicon(trgdetid);
  constexpr int kScintillatorIndex = 0;
  unsigned thickness = isSilicon ? triggerTools_.thicknessIndex(trgdetid) : kScintillatorIndex;
  if (thickness >= thicknessCorrection_.size()) {
    throw cms::Exception("OutOfBound") << "Trying to access thickness index " << thickness
                                       << " in thicknessCorrection, which is of size " << thicknessCorrection_.size();
  }

  /* weight the amplitude by the absorber coefficient in MeV/mip + bring it in
   * GeV */
  trgCellEt = trgCell.mipPt() * MevToGeV;
  trgCellEt *= dEdX_weights_.at(trgCellLayer);
  // [> correct for the cell-thickness <]
  if (thicknessCorrection_[thickness] > 0) {
    trgCellEt /= thicknessCorrection_[thickness];
  }

  math::PtEtaPhiMLorentzVector calibP4(trgCellEt, trgCell.eta(), trgCell.phi(), 0.);
  trgCell.setP4(calibP4);
}

void HGCalTriggerCellCalibration::calibrateInGeV(l1t::HGCalTriggerCell& trgCell) const {
  /* calibrate from ADC count to transverse mip */
  calibrateInMipT(trgCell);

  /* calibrate from mip count to GeV */
  calibrateMipTinGeV(trgCell);
}
