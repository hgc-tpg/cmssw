#include "L1Trigger/L1THGCal/interface/HGCalProcessorBase.h"

#include "DataFormats/L1THGCal/interface/HGCalCluster.h"
#include "DataFormats/L1THGCal/interface/HGCalMulticluster.h"
#include "L1Trigger/L1THGCal/interface/HGCalTriggerTools.h"
#include "L1Trigger/L1THGCal/interface/backend/HGCalShowerShape.h"
#include "L1Trigger/L1THGCal/interface/backend_semiemulator/Stage2.hh"

#include "FWCore/ParameterSet/interface/FileInPath.h"

#include <array>
#include <cmath>
#include <cstdint>
#include <unordered_map>
#include <utility>
#include <vector>

class HGCalBackendLayer2ProcessorSemiEmulator : public HGCalBackendLayer2ProcessorBase {
public:
  explicit HGCalBackendLayer2ProcessorSemiEmulator(const edm::ParameterSet& config)
      : HGCalBackendLayer2ProcessorBase(config),
        minimumClusterPt_(config.getParameter<double>("minimumClusterPt")),
        algorithm_(config.getParameter<double>("triangleSideLength") * std::sqrt(3.)) {
    const edm::FileInPath meanEtaLut(config.getParameter<edm::FileInPath>("meanEtaLUT"));
    const edm::FileInPath sigmaEtaLut(config.getParameter<edm::FileInPath>("sigmaEtaLUT"));
    clusterPropertyLut_.readMuEtaLUT(meanEtaLut.fullPath().c_str());
    clusterPropertyLut_.readSigmaEtaLUT(sigmaEtaLut.fullPath().c_str());
    algorithm_.setClusPropLUT(&clusterPropertyLut_);
  }

  void setGeometry(const HGCalTriggerGeometryBase* const geometry) override {
    HGCalBackendLayer2ProcessorBase::setGeometry(geometry);
    triggerTools_.setGeometry(geometry);
    showerShape_.setGeometry(geometry);
  }

  void run(const edm::Handle<l1t::HGCalClusterBxCollection>& input,
           std::pair<l1t::HGCalMulticlusterBxCollection, l1t::HGCalClusterBxCollection>& output) override {
    std::array<std::vector<edm::Ptr<l1t::HGCalCluster>>, 6> clustersPerSector;

    for (unsigned index = 0; index < input->size(); ++index) {
      edm::Ptr<l1t::HGCalCluster> cluster(input, index);
      const double absZ = std::abs(cluster->position().z());
      if (absZ == 0.)
        continue;

      const unsigned endcapOffset = cluster->position().z() > 0. ? 3 : 0;
      for (unsigned sector = 0; sector < 3; ++sector) {
        const unsigned sectorAndEndcap = sector + endcapOffset;
        TPGTCFloats rotated;
        rotated.setROverZPhiF(
            cluster->position().x() / absZ, cluster->position().y() / absZ, sectorAndEndcap);
        if (rotated.getXOverZF() >= 0.)
          clustersPerSector[sectorAndEndcap].push_back(cluster);
      }
    }

    for (unsigned sector = 0; sector < clustersPerSector.size(); ++sector)
      runSector(clustersPerSector[sector], sector, output.first);
  }

private:
  void runSector(const std::vector<edm::Ptr<l1t::HGCalCluster>>& inputs,
                 const unsigned sector,
                 l1t::HGCalMulticlusterBxCollection& output) {
    if (inputs.empty())
      return;

    std::vector<TPGTCBits> encodedInputs;
    encodedInputs.reserve(inputs.size());
    for (unsigned index = 0; index < inputs.size(); ++index) {
      const auto& input = inputs[index];
      const double absZ = std::abs(input->position().z());
      TPGTCFloats encoded;
      encoded.setZero();
      encoded.setROverZPhiF(input->position().x() / absZ, input->position().y() / absZ, sector);
      encoded.setEnergyGeV(input->pt());
      encoded.setLayer(triggerTools_.layerWithOffset(input->detId()));
      encoded.setCMSSWIndex(index);
      encodedInputs.push_back(encoded);
    }

    std::vector<TPGCluster> clusters;
    algorithm_.run(encodedInputs, clusters);
    for (const auto& cluster : clusters) {
      const double pt = cluster.getEnergyGeV();
      if (pt < minimumClusterPt_)
        continue;

      const double eta = cluster.getGlobalEtaRad(sector);
      const double phi = cluster.getGlobalPhiRad(sector);
      const math::PtEtaPhiMLorentzVector p4(pt, eta, phi, 0.);
      l1t::HGCalMulticluster outputCluster;
      outputCluster.setP4(p4);

      for (const int inputIndex : cluster.getCMSSWIndices()) {
        if (inputIndex >= 0 && static_cast<unsigned>(inputIndex) < inputs.size())
          outputCluster.addConstituent(inputs[inputIndex], false, 0.);
      }
      showerShape_.fillShapes(outputCluster, *geometry());

      const auto& hardware = cluster.getClData();
      outputCluster.setHwPt(hardware.e.to_uint());
      outputCluster.setHwEta(hardware.w_eta.to_uint());
      outputCluster.setHwPhi(hardware.w_phi.to_int());
      outputCluster.setHwQual(hardware.qualFlags.to_uint());
      outputCluster.setFirstLayer(hardware.firstLayer.to_uint());
      outputCluster.setShowerLength(hardware.showerLength.to_uint());
      outputCluster.setCoreShowerLength(hardware.coreShowerLength.to_uint());
      outputCluster.setZBarycenter(l1thgcfirmware::Scales::floatZ(hardware.w_z));
      outputCluster.setSigmaZZ(l1thgcfirmware::Scales::floatSigmaZ(hardware.sigma_z));
      outputCluster.setSigmaPhiPhiTot(l1thgcfirmware::Scales::floatSigmaPhi(hardware.sigma_phi));
      outputCluster.setSigmaEtaEtaTot(l1thgcfirmware::Scales::floatSigmaEta(hardware.sigma_eta));
      outputCluster.setSigmaRRTot(l1thgcfirmware::Scales::floatSigmaRozRoz(hardware.sigma_roz));

      if (hardware.e != 0) {
        const double electromagneticFraction = l1thgcfirmware::Scales::floatFrac(hardware.fractionInCE_E);
        outputCluster.saveEnergyInterpretation(l1t::HGCalMulticluster::EnergyInterpretation::EM,
                                               electromagneticFraction * outputCluster.energy());
      }
      output.push_back(0, outputCluster);
    }
  }

  double minimumClusterPt_;
  HGCalTriggerTools triggerTools_;
  HGCalShowerShape showerShape_;
  TPGStage2Configuration::ClusPropLUT clusterPropertyLut_;
  TPGStage2Emulation::Stage2 algorithm_;
};

DEFINE_EDM_PLUGIN(HGCalBackendLayer2Factory,
                  HGCalBackendLayer2ProcessorSemiEmulator,
                  "HGCalBackendLayer2ProcessorSemiEmulator");
