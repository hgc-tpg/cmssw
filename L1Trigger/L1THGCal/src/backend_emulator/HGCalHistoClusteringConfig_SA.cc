#include "L1Trigger/L1THGCal/interface/backend_emulator/HGCalHistoClusteringConfig_SA.h"
#include <math.h>
#include <iostream>
using namespace std;
using namespace l1thgcfirmware;
ClusterAlgoConfig::ClusterAlgoConfig() :
  clusterizerOffset_(0),
  cClocks_(0),
  cInputs_(0),
  cInputs2_(0),
  cInt_(0),
  cColumns_(0),
  cRows_(0),
  rOverZHistOffset_(0),
  rOverZBinSize_(0),
  kernelWidths_(),
  areaNormalizations_(),
  thresholdMaximaConstants_(),
  nBinsCosLUT_(0),
  cosLUT_(),
  clusterizerMagicTime_(0),
  stepLatency_({
    { Input  ,  0 },
    { Dist0  ,  0 },
    { Dist1  ,  0 },
    { Dist2  ,  0 },
    { Dist3  ,  0 },
    { Dist4  ,  0 },
    { Dist5  ,  0 },
    { TcToHc , 0 },
    { Hist   , 0 },
    { Smearing1D , 0 },
    { NormArea   , 0 },
    { Smearing2D , 0 },
    { Maxima1D   , 0 },
    { Maxima2D   , 0 },
    { CalcAverage , 0 },
    { Clusterizer , 0 },
    { TriggerCellToCluster , 0 },
    { ClusterSum , 0 }
  }),
  depths_(),
  triggerLayers_(),
  layerWeights_E_(),
  layerWeights_E_EM_(),
  layerWeights_E_EM_core_(),
  layerWeights_E_H_early_(),
  correction_(), // 0b011111111111111111
  saturation_() // (2 ** 19) - 1
{
}

void ClusterAlgoConfig::setStepLatencies( const std::vector<unsigned int> latencies ) {
  // Add check that stepLatency is at least same size as latencies
  // But not as cms.exception
  for ( unsigned int iStep = 0; iStep < latencies.size(); ++iStep ) {
    stepLatency_.at(Step(iStep)) = latencies.at(iStep);
  }
}

unsigned int ClusterAlgoConfig::getLatencyUpToAndIncluding( const Step step ) {
  unsigned int latency = 0;
  for ( int iStep = 0; iStep <= step; ++iStep ) latency += getStepLatency(Step(iStep));
  return latency;
}

void ClusterAlgoConfig::initializeLUTs() {
  initializeSmearingKernelConstants( cRows_, rOverZHistOffset_, rOverZBinSize_ );
  initializeThresholdMaximaConstants( cRows_, thresholdMaximaParam_a_, thresholdMaximaParam_b_, thresholdMaximaParam_c_  );
  initializeCosLUT();
}

void ClusterAlgoConfig::initializeSmearingKernelConstants( unsigned int bins, unsigned int offset, unsigned int height ) {
  const unsigned int lWidth0 = offset + (0.5*height);
  const unsigned int lTarget = int( (maxBinsSmearing1D_+0.5)*lWidth0 - 0.5 );
  for ( unsigned int iBin = 0; iBin < bins; ++iBin ) {
    unsigned int lCentre = lWidth0 + ( height * iBin );
    const unsigned int lBins = int( round(1.0 * lTarget / lCentre) );

    kernelWidths_.push_back( lBins );

    lCentre *= lBins;

    const unsigned int lRatio = int( round(1.0*lTarget/lCentre * pow(2,nBitsAreaNormLUT_) ) );

    areaNormalizations_.push_back( lRatio );
  }
}

void ClusterAlgoConfig::initializeThresholdMaximaConstants( unsigned int bins, unsigned int a, unsigned int b, int c ) {
  for ( unsigned int iBin = 0; iBin < bins; ++iBin ) {
    int threshold = a + b*iBin + c*iBin*iBin;
    thresholdMaximaConstants_.push_back( threshold );
  }
}

void ClusterAlgoConfig::initializeCosLUT() {
  for ( unsigned int iBin = 0; iBin < nBinsCosLUT_+1; ++iBin ) { 
    unsigned int cosBin = round( pow(2,nBitsCosLUT_) * ( 1 - cos(iBin*M_PI/phiNValues_) ) );
    cosLUT_.push_back( cosBin );
  }
}

void ClusterAlgoConfig::printConfiguration() {
  cout << "Running the algorithm" << endl;
  cout << "Config params" << endl;
  cout << "Clusterizer offset : " << clusterizerOffset() << endl;
  cout << "Latencies : ";
  for ( const auto& latency : stepLatency_ ) cout << latency.first << " " << latency.second << " ";
  cout << endl;
  cout << "cClocks : " << cClocks() << endl;
  cout << "cInputs : " << cInputs() << endl;
  cout << "cInputs2 : " << cInputs2() << endl;
  cout << "cColumns : " << cColumns() << endl;
  cout << "cRows : " << cRows() << endl;
  cout << "rOverZHistOffset : " << rOverZHistOffset() << endl;
  cout << "rOverZBinSize : " << rOverZBinSize() << endl;
  cout << "nBinsCosLUT : " << nBinsCosLUT() << endl;
  cout << "clusterizerMagicTime : " << clusterizerMagicTime() << endl;
  cout << "depths : ";
  for ( const auto& depth : depths() ) cout << depth << " ";
  cout << endl;
  cout << "triggerLayers : ";
  for ( const auto& triggerLayer : triggerLayers() ) cout << triggerLayer << " ";
  cout << endl;
  cout << "layerWeights_E : ";
  for ( const auto& weight : layerWeights_E() ) cout << weight << " ";
  cout << endl;
  cout << "layerWeights_E_EM : ";
  for ( const auto& weight : layerWeights_E_EM() ) cout << weight << " ";
  cout << endl;
  cout << "layerWeights_E_EM_core : ";
  for ( const auto& weight : layerWeights_E_EM_core() ) cout << weight << " ";
  cout << endl;
  cout << "layerWeights_E_H_early : ";
  for ( const auto& weight : layerWeights_E_H_early() ) cout << weight << " ";
  cout << endl;
  cout << "correction : " << correction() << endl;
  cout << "saturation : " << saturation() << endl;
}


