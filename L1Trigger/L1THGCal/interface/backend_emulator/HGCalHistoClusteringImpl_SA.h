#ifndef __L1Trigger_L1THGCal_HGCalHistoClusteringImplSA_h__
#define __L1Trigger_L1THGCal_HGCalHistoClusteringImplSA_h__

#include "L1Trigger/L1THGCal/interface/backend_emulator/HGCalTriggerCell_SA.h"
#include "L1Trigger/L1THGCal/interface/backend_emulator/HGCalHistogramCell_SA.h"
#include "L1Trigger/L1THGCal/interface/backend_emulator/CentroidHelper.h"
#include "L1Trigger/L1THGCal/interface/backend_emulator/HGCalCluster_SA.h"
#include "L1Trigger/L1THGCal/interface/backend_emulator/HGCalHistoClusteringConfig_SA.h"

#include <vector>

namespace l1thgcfirmware {

class HGCalHistoClusteringImplSA {
public:
  HGCalHistoClusteringImplSA(l1thgcfirmware::ClusterAlgoConfig& config);
  ~HGCalHistoClusteringImplSA() {}

  void runAlgorithm(HGCalTriggerCellSAPtrCollections& inputs, HGCalTriggerCellSAPtrCollection& clusteredTCs, HGCalTriggerCellSAPtrCollection& unclusteredTCs, CentroidHelperPtrCollection& prioritizedMaxima, CentroidHelperPtrCollection& readoutFlags,
  HGCalClusterSAPtrCollection& clusterSums ) const;

private:
  // TC input step
  l1thgcfirmware::HGCalTriggerCellSAPtrCollection triggerCellInput( l1thgcfirmware::HGCalTriggerCellSAPtrCollections& inputs ) const;

  // TC distribution steps
  void triggerCellDistribution0( l1thgcfirmware::HGCalTriggerCellSAPtrCollection& triggerCellsIn ) const;
  l1thgcfirmware::HGCalTriggerCellSAPtrCollections triggerCellDistribution1( l1thgcfirmware::HGCalTriggerCellSAPtrCollection& triggerCellsIn ) const;
  l1thgcfirmware::HGCalTriggerCellSAPtrCollections triggerCellDistribution2( l1thgcfirmware::HGCalTriggerCellSAPtrCollection& triggerCellsIn, l1thgcfirmware::HGCalTriggerCellSAPtrCollections& inTriggerCellDistributionGrid ) const;
  l1thgcfirmware::HGCalTriggerCellSAPtrCollections triggerCellDistribution3( l1thgcfirmware::HGCalTriggerCellSAPtrCollection& triggerCellsIn, l1thgcfirmware::HGCalTriggerCellSAPtrCollections& inTriggerCellDistributionGrid ) const;
  void triggerCellDistribution4( l1thgcfirmware::HGCalTriggerCellSAPtrCollection& triggerCellsIn ) const;
  void triggerCellDistribution5( l1thgcfirmware::HGCalTriggerCellSAPtrCollection& triggerCellsIn, l1thgcfirmware::HGCalTriggerCellSAPtrCollections& inTriggerCellDistributionGrid ) const;

  // Histogram steps
  l1thgcfirmware::HGCalHistogramCellSAPtrCollection triggerCellToHistogramCell( l1thgcfirmware::HGCalTriggerCellSAPtrCollection& triggerCellsIn ) const;
  l1thgcfirmware::HGCalHistogramCellSAPtrCollection makeHistogram( l1thgcfirmware::HGCalHistogramCellSAPtrCollection histogramCells ) const;

  // Smearing steps
  void smearing1D( l1thgcfirmware::HGCalHistogramCellSAPtrCollection& histogram ) const;
  void areaNormalization( l1thgcfirmware::HGCalHistogramCellSAPtrCollection& histogram ) const;
  void smearing2D( l1thgcfirmware::HGCalHistogramCellSAPtrCollection& histogram ) const;

  // Maxima finding
  void thresholdMaximaFinder( l1thgcfirmware::HGCalHistogramCellSAPtrCollection& histogram ) const;
  void calculateAveragePosition( l1thgcfirmware::HGCalHistogramCellSAPtrCollection& histogram ) const;

  // Clustering
  // l1thgcfirmware::DOICollection domainOfInfluence( l1thgcfirmware::HGCalHistogramCellSAPtrCollection& histogram ) const;
  void clusterizer( l1thgcfirmware::HGCalTriggerCellSAPtrCollection& triggerCellsIn, l1thgcfirmware::HGCalHistogramCellSAPtrCollection& histogram, l1thgcfirmware::HGCalTriggerCellSAPtrCollection& clusteredTriggerCells, l1thgcfirmware::HGCalTriggerCellSAPtrCollection& unclusteredTriggerCells, l1thgcfirmware::CentroidHelperPtrCollection& prioritizedMaxima, l1thgcfirmware::CentroidHelperPtrCollection& readoutFlags ) const;

  // Cluster properties
  l1thgcfirmware::HGCalClusterSAPtrCollection triggerCellToCluster( l1thgcfirmware::HGCalTriggerCellSAPtrCollection& clusteredTriggerCells ) const;
  void clusterSum( l1thgcfirmware::HGCalClusterSAPtrCollection& protoClusters, l1thgcfirmware::CentroidHelperPtrCollection& readoutFlags, l1thgcfirmware::HGCalClusterSAPtrCollection& clusterAccumulation, l1thgcfirmware::HGCalClusterSAPtrCollection& clusterSums ) const;
  std::pair< unsigned int, unsigned int > sigma_Energy(unsigned int N_TC_W, unsigned long int Sum_W2, unsigned int Sum_W) const;
  std::pair< unsigned int, unsigned int > mean_coordinate(unsigned int Sum_Wc, unsigned int Sum_W) const;
  std::pair< unsigned int, unsigned int > sigma_Coordinate(unsigned int Sum_W, unsigned long int Sum_Wc2, unsigned int Sum_Wc) const;
  std::pair< unsigned int, unsigned int > energy_ratio(unsigned int E_N, unsigned int E_D) const;
  std::vector<int> showerLengthProperties(unsigned long int layerBits) const;
  void clusterProperties(l1thgcfirmware::HGCalClusterSAPtrCollection& clusterSums) const;


  // Useful functions
  void initializeTriggerCellDistGrid( l1thgcfirmware::HGCalTriggerCellSAPtrCollections& grid, unsigned int nX, unsigned int nY ) const;

  void runDistServers( const l1thgcfirmware::HGCalTriggerCellSAPtrCollections& gridIn,
                       l1thgcfirmware::HGCalTriggerCellSAPtrCollections& gridOut,
                       l1thgcfirmware::HGCalTriggerCellSAPtrCollection& tcsOut,
                       unsigned int latency,
                       unsigned int nDistServers,
                       unsigned int nInputs,
                       unsigned int nOutputs,
                       unsigned int nInterleave,
                      bool setOutputGrid ) const;

  l1thgcfirmware::ClusterAlgoConfig& config_;
};

}

#endif
