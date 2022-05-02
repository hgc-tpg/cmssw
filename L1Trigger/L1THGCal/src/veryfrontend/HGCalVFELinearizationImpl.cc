#include "L1Trigger/L1THGCal/interface/veryfrontend/HGCalVFELinearizationImpl.h"
#include "SimCalorimetry/HGCalSimAlgos/interface/HGCalSiNoiseMap.h"

#include <cmath>

HGCalVFELinearizationImpl::HGCalVFELinearizationImpl(const edm::ParameterSet& conf)
    : linLSB_(conf.getParameter<double>("linLSB")),
      adcsaturation_(conf.getParameter<double>("adcsaturation")),
      tdcnBits_(conf.getParameter<uint32_t>("tdcnBits")),
      tdcOnset_(conf.getParameter<double>("tdcOnset")),
      adcnBits_(conf.getParameter<uint32_t>("adcnBits")),
      tdcsaturation_(conf.getParameter<double>("tdcsaturation")),
      linnBits_(conf.getParameter<uint32_t>("linnBits")),
      oot_coefficients_(conf.getParameter<std::vector<double>>("oot_coefficients")),
      new_digi_(conf.getParameter<bool>("newDigi")) {
  constexpr int kOot_order = 2;
  if (oot_coefficients_.size() != kOot_order) {
    throw cms::Exception("BadConfiguration") << "OOT subtraction needs " << kOot_order << " coefficients";
  }
  adcLSB_ = std::ldexp(adcsaturation_, -adcnBits_);
  tdcLSB_ = std::ldexp(tdcsaturation_, -tdcnBits_);
  linMax_ = (0x1 << linnBits_) - 1;

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

void HGCalVFELinearizationImpl::setGeometry(const HGCalTriggerGeometryBase* const geom, DetId::Detector det) {
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

void HGCalVFELinearizationImpl::linearize(const std::vector<HGCalDataFrame>& dataframes,
                                          std::vector<std::pair<DetId, uint32_t>>& linearized_dataframes) {
  constexpr int kIntimeSample = 2;
  constexpr int kOuttime1Sample = 1;  // in time - 1;
  constexpr int kOuttime2Sample = 0;  // in time - 2;

  for (const auto& frame : dataframes) {  //loop on DIGI
    bool isTDC(frame[kIntimeSample].mode());
    double rawData(double(frame[kIntimeSample].data()));
    bool isBusy(isTDC && rawData == 0);

    double adcLSB = adcLSB_;
    double noise = 0.;
    if (new_digi_) {
      HGCalSiNoiseMap<HGCSiliconDetId>::SiCellOpCharacteristics siop =
          noise_map_.getSiCellOpCharacteristics(frame.id());
      HGCalSiNoiseMap<HGCSiliconDetId>::GainRange_t gain((HGCalSiNoiseMap<HGCSiliconDetId>::GainRange_t)siop.core.gain);
      adcLSB = noise_map_.getLSBPerGain()[gain];
      noise = siop.core.noise;
    }

    double amplitude = 0.;
    if (isBusy) {
      continue;
    }
    if (isTDC) {
      amplitude = (std::floor(tdcOnset_ / adcLSB) + 1.0) * adcLSB + rawData * tdcLSB_;
    } else {  //ADC mode
      double data = rawData;
      // applies OOT PU subtraction only in the ADC mode
      if (!frame[kOuttime1Sample].mode()) {
        data += oot_coefficients_[kOuttime1Sample] * frame[kOuttime1Sample].data();
        if (!frame[kOuttime2Sample].mode()) {
          data += oot_coefficients_[kOuttime2Sample] * frame[kOuttime2Sample].data();
        }
      }
      amplitude = std::max(0., data) * adcLSB;
    }
    if (amplitude < 3. * noise) {
      continue;
    }
    uint32_t amplitude_int = uint32_t(std::floor(amplitude / linLSB_ + 0.5));
    if (amplitude_int == 0)
      continue;
    if (amplitude_int > linMax_)
      amplitude_int = linMax_;

    linearized_dataframes.emplace_back(frame.id(), amplitude_int);
  }
}
