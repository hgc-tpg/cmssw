#include "L1Trigger/L1THGCal/interface/HGCalProcessorBase.h"

#include "DataFormats/ForwardDetId/interface/HGCalTriggerBackendDetId.h"
#include "DataFormats/L1THGCal/interface/HGCalCluster.h"
#include "DataFormats/L1THGCal/interface/HGCalMulticluster.h"
#include "L1Trigger/L1THGCal/interface/HGCalTriggerGeometryBase.h"
#include "L1Trigger/L1THGCal/interface/backend/HGCalStage2ClusterDistribution.h"

#include "L1Trigger/L1THGCal/interface/HGC/HGC.hpp"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <unordered_map>
#include <utility>
#include <vector>

namespace {
  constexpr unsigned kFrames = 162;
  constexpr unsigned kDecoderLanes = 154;
  constexpr double kPi = 3.14159265358979323846;
  constexpr unsigned kHistogramRows = 40;
  constexpr int kMinHistogramColumn = -21;
  constexpr int kMaxHistogramColumn = 22;

  uint64_t mix(uint64_t value) {
    value += 0x9e3779b97f4a7c15ULL;
    value = (value ^ (value >> 30)) * 0xbf58476d1ce4e5b9ULL;
    value = (value ^ (value >> 27)) * 0x94d049bb133111ebULL;
    return value ^ (value >> 31);
  }

  double wrapPhi(double phi) {
    while (phi > kPi)
      phi -= 2. * kPi;
    while (phi <= -kPi)
      phi += 2. * kPi;
    return phi;
  }

  double localPhi(const GlobalPoint& position, const int sector) {
    const double x = position.z() > 0. ? -position.x() : position.x();
    double phi = std::atan2(position.y(), x);
    if (sector == 1)
      phi += (phi > 0. ? -2. : 4.) * kPi / 3.;
    else if (sector == 2)
      phi += 2. * kPi / 3.;
    return phi;
  }

  double globalPhi(const double sectorZeroPhi, const short zside, const int sector) {
    double phi = sectorZeroPhi;
    if (sector == 1)
      phi += 2. * kPi / 3.;
    else if (sector == 2)
      phi += 4. * kPi / 3.;
    if (zside > 0)
      phi = kPi - phi;
    return wrapPhi(phi);
  }

  unsigned histogramRow(const GlobalPoint& position) {
    const double roverz = std::hypot(position.x(), position.y()) / std::abs(position.z());
    const double fraction = (roverz - HGC::c_roz_global_min) / (HGC::c_roz_global_max - HGC::c_roz_global_min);
    return std::clamp(static_cast<int>(std::floor(fraction * kHistogramRows)),
                      0,
                      static_cast<int>(kHistogramRows - 1));
  }

  int histogramColumn(const double phi) {
    return std::clamp(static_cast<int>(std::lround((phi - kPi / 2.) / HGC::c_HistoPhiBin)),
                      kMinHistogramColumn,
                      kMaxHistogramColumn);
  }

  unsigned histogramCell(const unsigned row, const int column) {
    int phase = (column - 3) % 4;
    if (phase < 0)
      phase += 4;
    return 4 * row + phase;
  }

  int outputColumn(const unsigned lane, const unsigned cell) {
    return DecoderColumns[lane] + ((447 + static_cast<int>(cell) - DecoderColumns[lane]) % 4);
  }

  LinkTriggerCell encodeTriggerCell(const l1t::HGCalCluster& cluster, const unsigned tcId, const double energyScale) {
    const auto energy = static_cast<uint64_t>(std::max(0., std::round(cluster.pt() * energyScale)));
    uint8_t exponent = 0;
    uint8_t mantissa = 0;
    if (energy < 16) {
      mantissa = energy;
    } else {
      uint64_t shiftedEnergy = energy;
      while (shiftedEnergy > 31 && exponent < 31) {
        shiftedEnergy >>= 1;
        ++exponent;
      }
      mantissa = shiftedEnergy & 0xf;
    }
    return LinkTriggerCell(mantissa, exponent, tcId, false, true);
  }

  TcDecoder makeDecoder(const l1t::HGCalCluster& cluster,
                        const HGCalTriggerGeometryBase& geometry,
                        const unsigned cell,
                        const double phi) {
    const auto& position = cluster.position();
    const double absZ = std::abs(position.z());
    const double roverz = absZ > 0. ? std::hypot(position.x(), position.y()) / absZ : 0.;
    const unsigned layer = geometry.triggerLayer(cluster.detId());

    TcDecoder decoded;
    // The current hardware decoder stores values after the fixed-point shifts
    // applied in Cluster_Decoders. These are temporary V19-derived values until
    // the generated FE/BE decoder ROM replaces the virtual transport map.
    decoded.RoverZ = std::min(8191u, static_cast<unsigned>(std::round(roverz / HGC::c_LSB_roz_TC)));
    const double boundedPhi = std::clamp(phi, 0., kPi);
    decoded.Phi = std::min(4095u,
                           static_cast<unsigned>(std::round(boundedPhi / HGC::c_LSB_phi_TC)));
    decoded.TriggerLayer = layer;
    decoded.LayerDepth = std::min(4095u, static_cast<unsigned>(std::round(absZ / 0.5)));
    decoded.Cell = cell;
    decoded.LayerWeight = layer > 0 && layer < 15 ? 262144u : 0u;
    decoded.DataValid = true;
    return decoded;
  }

  l1t::HGCalMulticluster makeMulticluster(const ClusterProperty& property,
                                          const HGCalTriggerBackendDetId& fpga) {
    const int localPhi = property.ET_Weighted_Phi & 0x100 ? property.ET_Weighted_Phi - 0x200
                                                           : property.ET_Weighted_Phi;
    const double eta = fpga.zside() * property.ET_Weighted_Eta * kPi / 720.;
    const double phi = globalPhi(localPhi * kPi / 720. + kPi / 2., fpga.zside(), fpga.sector());
    const double pt = property.ET * 0.25;
    const math::PtEtaPhiMLorentzVector polarP4(pt, eta, phi, 0.);
    const l1t::HGCalMulticluster::LorentzVector p4(
        polarP4.Px(), polarP4.Py(), polarP4.Pz(), polarP4.E());
    l1t::HGCalMulticluster output(p4, property.ET, property.ET_Weighted_Eta, localPhi);
    output.setHwQual((property.Saturated_Trigger_Cell << 0) |
                     (property.Quality_of_Sigmas_and_Means << 1) |
                     (property.Saturated_Phi << 2));
    output.setFirstLayer(property.First_Layer);
    output.setMaxLayer(property.Last_Layer);
    output.setShowerLength(property.Shower_Length);
    output.setCoreShowerLength(property.CoreShowerLen);
    output.setSigmaZZ(property.Sigma_ZZ);
    output.setSigmaPhiPhiTot(property.Sigma_PhiPhi);
    output.setSigmaEtaEtaTot(property.Sigma_EtaEta);
    output.setSigmaRRTot(property.Sigma_RoZRoZ);
    return output;
  }
}  // namespace

class HGCalBackendLayer2ProcessorStage2Emulator : public HGCalBackendLayer2ProcessorBase {
public:
  explicit HGCalBackendLayer2ProcessorStage2Emulator(const edm::ParameterSet& config)
      : HGCalBackendLayer2ProcessorBase(config),
        distributor_(config.getParameterSet("DistributionParameters")),
        energyScale_(config.getParameter<double>("virtualTransportEnergyScale")),
        debug_(config.getUntrackedParameter<bool>("debug", false)) {}

  void run(const edm::Handle<l1t::HGCalClusterBxCollection>& input,
           std::pair<l1t::HGCalMulticlusterBxCollection, l1t::HGCalClusterBxCollection>& output) override {
    std::unordered_map<uint32_t, std::vector<edm::Ptr<l1t::HGCalCluster>>> clustersPerFpga;
    for (unsigned index = 0; index < input->size(); ++index) {
      edm::Ptr<l1t::HGCalCluster> cluster(input, index);
      const unsigned module = geometry()->getModuleFromTriggerCell(cluster->detId());
      const unsigned stage1Fpga = geometry()->getStage1FpgaFromModule(module);
      const auto candidateFpgas = geometry()->getStage2FpgasFromStage1Fpga(stage1Fpga);
      const auto stage2Fpgas = distributor_.getStage2FPGAs(stage1Fpga, candidateFpgas, cluster);
      for (const auto fpga : stage2Fpgas)
        clustersPerFpga[fpga].push_back(cluster);
    }

    for (const auto& [fpga, clusters] : clustersPerFpga) {
      runSector(clusters, HGCalTriggerBackendDetId(fpga), output.first);
    }
  }

private:
  void runSector(const std::vector<edm::Ptr<l1t::HGCalCluster>>& clusters,
                 const HGCalTriggerBackendDetId& fpga,
                 l1t::HGCalMulticlusterBxCollection& output) const {
    std::array<std::array<LinkTriggerCell, kDecoderLanes>, kFrames> virtualInputs{};
    std::array<std::array<bool, kDecoderLanes>, kFrames> occupiedLanes{};
    std::array<std::array<Cluster, kDecoderLanes>, kFrames> protoClusters{};
    unsigned assignedCells = 0;
    unsigned unassignedCells = 0;
    double inputPt = 0.;

    for (const auto& cluster : clusters) {
      const double phi = localPhi(cluster->position(), fpga.sector());
      if (phi < 0. || phi > kPi)
        continue;
      const int column = histogramColumn(phi);
      const unsigned cell = histogramCell(histogramRow(cluster->position()), column);
      const uint64_t hash = mix(cluster->detId());
      bool assigned = false;
      for (unsigned laneOffset = 0; laneOffset < kDecoderLanes && !assigned; ++laneOffset) {
        const unsigned lane = (hash + laneOffset) % kDecoderLanes;
        if (outputColumn(lane, cell) != column)
          continue;
        for (unsigned frameOffset = 0; frameOffset < kFrames; ++frameOffset) {
          const unsigned frame = ((hash >> 16) + frameOffset) % kFrames;
          if (occupiedLanes[frame][lane])
            continue;
          virtualInputs[frame][lane] = encodeTriggerCell(*cluster, hash & 0x3f, energyScale_);
          const auto decoder = makeDecoder(*cluster, *geometry(), cell, phi);
          HGC::Cluster_Decode(virtualInputs[frame][lane], decoder, protoClusters[frame][lane]);
          occupiedLanes[frame][lane] = true;
          assigned = true;
          ++assignedCells;
          inputPt += cluster->pt();
          break;
        }
      }
      // A future real transport map replaces this bounded virtual allocator.
      // Until then, unassigned cells are the explicitly modelled transport overflow.
      if (!assigned)
        ++unassignedCells;
    }

    std::array<ClusterProperty, kFrames> properties{};
    HGC::Clusters_Step2(protoClusters, properties);
    unsigned outputClusters = 0;
    double outputPt = 0.;
    for (const auto& property : properties) {
      // Each Stage-2 sector processes an overlapping ~180 degree region, but
      // owns only its central 120 degrees.  Keep the overlap for clustering
      // and publish the cluster only from its nominal sector.
      if (property.DataValid && property.Nominal_Phi && property.ET > 0) {
        output.push_back(0, makeMulticluster(property, fpga));
        ++outputClusters;
        outputPt += property.ET * 0.25;
      }
    }
    if (debug_)
      edm::LogVerbatim("HGCalStage2Emulator")
          << "zside=" << fpga.zside() << " sector=" << fpga.sector() << " input=" << clusters.size()
          << " accepted=" << assignedCells << " overflow=" << unassignedCells << " inputPt=" << inputPt
          << " output=" << outputClusters << " outputPt=" << outputPt;
  }

  HGCalStage2ClusterDistribution distributor_;
  double energyScale_;
  bool debug_;
};

DEFINE_EDM_PLUGIN(HGCalBackendLayer2Factory,
                  HGCalBackendLayer2ProcessorStage2Emulator,
                  "HGCalBackendLayer2ProcessorStage2Emulator");
