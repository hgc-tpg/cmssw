#ifndef __L1Trigger_L1THGCal_HGCalHistoCluteringConfig_SA_h__
#define __L1Trigger_L1THGCal_HGCalHistoCluteringConfig_SA_h__

#include <map>
#include <vector>

namespace l1thgcfirmware {

  enum Step { Uninitialized = -1, 
              Input = 0,
              Dist0 = 1,
              Dist1 = 2,
              Dist2 = 3,
              Dist3 = 4,
              Dist4 = 5,
              Dist5 = 6,
              TcToHc = 7,
              Hist = 8,
              Smearing1D = 9,
              NormArea = 10,
              Smearing2D = 11,
              Maxima1D = 12, // Not actually used currently
              Maxima2D = 13,
              CalcAverage = 14,
              Clusterizer = 15,
              TriggerCellToCluster = 16,
              ClusterSum = 17
            };

  class ClusterAlgoConfig {
  public:
    ClusterAlgoConfig();

    void setParameters() {}

    void setSector( const unsigned int sector ) { sector_ = sector; }
    unsigned int sector() const { return sector_; }

    void setZSide( const int zSide ) { zSide_ = zSide; }
    int zSide() const { return zSide_; }

    void setStepLatencies( const std::vector<unsigned int> latencies );
    unsigned int getStepLatency( const Step step ) const { return stepLatency_.at(step); }
    unsigned int getLatencyUpToAndIncluding( const Step step );

    void setClusterizerOffset( const unsigned clusterizerOffset ) { clusterizerOffset_ = clusterizerOffset ;}
    unsigned int clusterizerOffset() const { return clusterizerOffset_; }

    void setCClocks( const unsigned cClocks ) { cClocks_ = cClocks;}
    unsigned int cClocks() const { return cClocks_; }

    void setCInputs( const unsigned cInputs ) { cInputs_ = cInputs;}
    unsigned int cInputs() const { return cInputs_; }

    void setCInputs2( const unsigned cInputs2 ) { cInputs2_ = cInputs2;}
    unsigned int cInputs2() const { return cInputs2_; }

    void setCInt( const unsigned cInt ) { cInt_ = cInt;}
    unsigned int cInt() const { return cInt_; }

    void setCColumns( const unsigned cColumns ) { cColumns_ = cColumns;}
    unsigned int cColumns() const { return cColumns_; }

    void setCRows( const unsigned cRows ) { cRows_ = cRows;}
    unsigned int cRows() const { return cRows_; }

    void setROverZHistOffset( const unsigned rOverZHistOffset ) { rOverZHistOffset_ = rOverZHistOffset;}
    unsigned int rOverZHistOffset() const { return rOverZHistOffset_; }

    void setROverZBinSize( const unsigned rOverZBinSize ) { rOverZBinSize_ = rOverZBinSize;}
    unsigned int rOverZBinSize() const { return rOverZBinSize_; }

    unsigned int kernelWidth( unsigned int iBin ) const { return kernelWidths_.at(iBin); }
    unsigned int areaNormalization( unsigned int iBin ) const { return areaNormalizations_.at(iBin); }
    
    void setROverZRange( const float rOverZRange ) { rOverZRange_ = rOverZRange_; }
    unsigned int rOverZRange() const { return rOverZRange_; }
    void setROverZNValues( const float rOverZNValues ) { rOverZNValues_ = rOverZNValues_; }
    unsigned int rOverZNValues() const { return rOverZNValues_; }
    void setPhiRange( const float phiRange ) { phiRange_ = phiRange_; }
    unsigned int phiRange() const { return phiRange_; }
    void setPhiNValues( const float phiNValues ) { phiNValues_ = phiNValues_; }
    unsigned int phiNValues() const { return phiNValues_; }

    void setThresholdParams( const unsigned int a, const unsigned int b, const int c ) { 
      thresholdMaximaParam_a_ = a;
      thresholdMaximaParam_b_ = b;
      thresholdMaximaParam_c_ = c;
    }
    unsigned int thresholdMaxima( unsigned int iBin ) const { return thresholdMaximaConstants_.at(iBin); }

    void setNBinsCosLUT( const unsigned int nBinsCosLUT ) { nBinsCosLUT_ = nBinsCosLUT; }
    unsigned int nBinsCosLUT() const { return nBinsCosLUT_; }

    unsigned int cosLUT( unsigned int iBin ) const { return cosLUT_.at(iBin); }

    void setClusterizerMagicTime( const unsigned clusterizerMagicTime ) { clusterizerMagicTime_ = clusterizerMagicTime;}
    unsigned int clusterizerMagicTime() const { return clusterizerMagicTime_; }

    void setDepths( const std::vector<unsigned int> depths ) {
      depths_.clear();
      for ( const auto& depth : depths ) depths_.push_back( depth );
    }
    std::vector<unsigned int> depths() const { return depths_; }
    unsigned int depth( unsigned int iLayer ) const { return depths_.at(iLayer); }

    void setTriggerLayers( const std::vector<unsigned int> triggerLayers ) {
      triggerLayers_.clear();
      for ( const auto& triggerLayer : triggerLayers ) triggerLayers_.push_back( triggerLayer );
    }
    std::vector<unsigned int> triggerLayers() const { return triggerLayers_; }
    unsigned int triggerLayer( unsigned int iLayer ) const { return triggerLayers_.at(iLayer); }
    
    void setLayerWeights_E( const std::vector<unsigned int> layerWeights_E ) {
      layerWeights_E_.clear();
      for ( const auto& weight : layerWeights_E ) layerWeights_E_.push_back( weight );
    }
    std::vector<unsigned int> layerWeights_E() const { return layerWeights_E_; }
    unsigned int layerWeight_E( unsigned int iTriggerLayer ) const { return layerWeights_E_.at(iTriggerLayer); }

    void setLayerWeights_E_EM( const std::vector<unsigned int> layerWeights_E_EM ) {
      layerWeights_E_EM_.clear();
      for ( const auto& weight : layerWeights_E_EM ) layerWeights_E_EM_.push_back( weight );
    }
    std::vector<unsigned int> layerWeights_E_EM() const { return layerWeights_E_EM_; }
    unsigned int layerWeight_E_EM( unsigned int iTriggerLayer ) const { return layerWeights_E_EM_.at(iTriggerLayer); }
    
    void setLayerWeights_E_EM_core( const std::vector<unsigned int> layerWeights_E_EM_core ) {
      layerWeights_E_EM_core_.clear();
      for ( const auto& weight : layerWeights_E_EM_core ) layerWeights_E_EM_core_.push_back( weight );
    }
    std::vector<unsigned int> layerWeights_E_EM_core() const { return layerWeights_E_EM_core_; }
    unsigned int layerWeight_E_EM_core( unsigned int iTriggerLayer ) const { return layerWeights_E_EM_core_.at(iTriggerLayer); }


    void setLayerWeights_E_H_early( const std::vector<unsigned int> layerWeights_E_H_early ) {
      layerWeights_E_H_early_.clear();
      for ( const auto& weight : layerWeights_E_H_early ) layerWeights_E_H_early_.push_back( weight );
    }
    std::vector<unsigned int> layerWeights_E_H_early() const { return layerWeights_E_H_early_; }
    unsigned int layerWeight_E_H_early( unsigned int iTriggerLayer ) const { return layerWeights_E_H_early_.at(iTriggerLayer); }

    void setCorrection( const unsigned correction ) { correction_ = correction;}
    unsigned int correction() const { return correction_; }

    void setSaturation( const unsigned saturation ) { saturation_ = saturation;}
    unsigned int saturation() const { return saturation_; }

    void initializeLUTs();

    void printConfiguration(); // For debugging

  private:
    void initializeSmearingKernelConstants( unsigned int bins, unsigned int offset, unsigned int height );
    void initializeThresholdMaximaConstants( unsigned int bins, unsigned int a, unsigned int b, int c );
    void initializeCosLUT();

    unsigned int histogramOffset_;
    unsigned int clusterizerOffset_;
    unsigned int cClocks_;
    unsigned int cInputs_;
    unsigned int cInputs2_; // Better name for variable?
    unsigned int cInt_;
    unsigned int cColumns_;
    unsigned int cRows_;
    unsigned int rOverZHistOffset_;
    unsigned int rOverZBinSize_;

    float rOverZRange_;
    float rOverZNValues_;
    float phiRange_;
    float phiNValues_;

    std::vector<unsigned int> kernelWidths_;
    std::vector<unsigned int> areaNormalizations_;

    unsigned int thresholdMaximaParam_a_;
    unsigned int thresholdMaximaParam_b_;
    int thresholdMaximaParam_c_;
    std::vector<int> thresholdMaximaConstants_;

    unsigned int nBinsCosLUT_;
    std::vector<unsigned int> cosLUT_;

    unsigned int clusterizerMagicTime_;

    std::map<Step,unsigned int> stepLatency_;

    // Parameters for triggerCellToCluster
    std::vector<unsigned int> depths_;
    std::vector<unsigned int> triggerLayers_;
    std::vector<unsigned int> layerWeights_E_;
    std::vector<unsigned int> layerWeights_E_EM_;
    std::vector<unsigned int> layerWeights_E_EM_core_;
    std::vector<unsigned int> layerWeights_E_H_early_;
    unsigned int correction_;
    unsigned int saturation_;

    unsigned int sector_;
    int zSide_;

  };

}  // namespace l1thgcfirmware

#endif
