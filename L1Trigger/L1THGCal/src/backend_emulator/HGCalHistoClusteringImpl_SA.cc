#include "L1Trigger/L1THGCal/interface/backend_emulator/HGCalHistoClusteringImpl_SA.h"

using namespace std;
using namespace l1thgcfirmware;
HGCalHistoClusteringImplSA::HGCalHistoClusteringImplSA( ClusterAlgoConfig& config ) : config_(config), tcDistribution_(config), seeding_(config), clustering_(config), clusterProperties_(config) {}

void HGCalHistoClusteringImplSA::runAlgorithm(const HGCalTriggerCellSAPtrCollections& inputs, HGCalTriggerCellSAPtrCollection& clusteredTCs, HGCalTriggerCellSAPtrCollection& unclusteredTCs, HGCalClusterSAPtrCollection& clusterSums ) const {
  // config_.printConfiguration();

  // TC distribution
  HGCalTriggerCellSAPtrCollection distributedTCs;
  tcDistribution_.runTriggerCellDistribution(inputs,distributedTCs);
 
  // Histogramming and seeding
  HGCalHistogramCellSAPtrCollection histogram;
  seeding_.runSeeding(distributedTCs, histogram);

  // Clustering
  CentroidHelperPtrCollection readoutFlags;
  clustering_.runClustering( distributedTCs, histogram, clusteredTCs, readoutFlags );

  // Cluster properties
  clusterProperties_.runClusterProperties( clusteredTCs, readoutFlags, clusterSums);

}


