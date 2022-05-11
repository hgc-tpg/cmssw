#include "L1Trigger/L1THGCal/interface/backend_emulator/HGCalTCDistribution_SA.h"

using namespace std;
using namespace l1thgcfirmware;

HGCalTCDistribution::HGCalTCDistribution(ClusterAlgoConfig& config) : config_(config) {}

void HGCalTCDistribution::runTriggerCellDistribution(const HGCalTriggerCellSAPtrCollections& triggerCellsIn, HGCalTriggerCellSAPtrCollection& triggerCellsOut ) const {
 
    HGCalTriggerCellSAShrPtrCollection triggerCellsOutTemp;
    triggerCellInput( triggerCellsIn, triggerCellsOutTemp );
    triggerCellDistribution0( triggerCellsOutTemp );
    HGCalTriggerCellSAShrPtrCollections tcDistGrid1;
    triggerCellDistribution1( triggerCellsOutTemp, tcDistGrid1 );
    HGCalTriggerCellSAShrPtrCollections tcDistGrid2;
    triggerCellDistribution2( tcDistGrid1, triggerCellsOutTemp, tcDistGrid2 );
    HGCalTriggerCellSAShrPtrCollections tcDistGrid3;
    triggerCellDistribution3( tcDistGrid2, triggerCellsOutTemp, tcDistGrid3 );
    triggerCellDistribution4( triggerCellsOutTemp );
    triggerCellDistribution5( tcDistGrid3, triggerCellsOut );

}

void HGCalTCDistribution::triggerCellInput( const HGCalTriggerCellSAPtrCollections& inputTCs, HGCalTriggerCellSAShrPtrCollection& outputTCs ) const {

  outputTCs.clear();
  for (unsigned int iFrame = 0; iFrame < inputTCs.size(); ++iFrame ) {
    for (unsigned int iInput = 0; iInput < inputTCs[iFrame].size(); ++iInput ) {
      auto& tc = inputTCs[iFrame][iInput];

      tc->setIndex(iInput);
      tc->setClock(iFrame+1);

      if ( tc->dataValid() ) {
        outputTCs.push_back(make_shared<HGCalTriggerCell>(*tc));
      }
    }
  }
}

void HGCalTCDistribution::triggerCellDistribution0( const HGCalTriggerCellSAShrPtrCollection& triggerCellsIn ) const {
  for ( auto& tc : triggerCellsIn ) {
    unsigned int newIndex = tc->index() + int( tc->index() / 24 ); // Magic numbers
    tc->setIndex( newIndex );
  }
}

void HGCalTCDistribution::triggerCellDistribution1( const HGCalTriggerCellSAShrPtrCollection& triggerCellsIn, HGCalTriggerCellSAShrPtrCollections& outTriggerCellDistributionGrid ) const {

  initializeTriggerCellDistGrid( outTriggerCellDistributionGrid, config_.cClocks(), config_.cInputs2() );

  const unsigned int stepLatency = config_.getStepLatency( Dist1 );
  for ( auto& tc : triggerCellsIn ) {
    tc->addLatency( stepLatency );
    unsigned int sector = int( tc->index() / 25 ); // Magic numbers
    outTriggerCellDistributionGrid[tc->clock()-2][tc->index()] = tc; // Magic numbers
    for ( int iSortKey = 5; iSortKey >= 0; --iSortKey ) { // Magic numbers
      // Split each 60 degree sector into 6 phi region
      // Sort key is index of small phi region
      if ( int( tc->phi() % 1944 ) > int( 108 * iSortKey + 648 * sector ) ) { // Magic numbers
        tc->setSortKey( iSortKey );
        break;
      }
    }
  }
}

void HGCalTCDistribution::triggerCellDistribution2( const HGCalTriggerCellSAShrPtrCollections& inTriggerCellDistributionGrid, HGCalTriggerCellSAShrPtrCollection& triggerCellsOut, HGCalTriggerCellSAShrPtrCollections& outTriggerCellDistributionGrid ) const {

  const unsigned int latency = config_.getLatencyUpToAndIncluding( Dist2 );

  triggerCellsOut.clear();
  initializeTriggerCellDistGrid( outTriggerCellDistributionGrid, config_.cClocks(), config_.cInt() );

  runDistServers(inTriggerCellDistributionGrid,
                  outTriggerCellDistributionGrid,
                  triggerCellsOut,
                  latency,
                  15, 5, 6, 4, true); // Magic numbers

}

void HGCalTCDistribution::triggerCellDistribution3( const HGCalTriggerCellSAShrPtrCollections& inTriggerCellDistributionGrid, HGCalTriggerCellSAShrPtrCollection& triggerCellsOut, HGCalTriggerCellSAShrPtrCollections& outTriggerCellDistributionGrid ) const {
  triggerCellsOut.clear();
  initializeTriggerCellDistGrid( outTriggerCellDistributionGrid, config_.cClocks(), config_.cInt() );
  for ( unsigned int iClock = 0; iClock < config_.cClocks(); ++iClock ) {
    for ( unsigned int i = 0; i < 3; ++i ) { // Magic numbers
      for ( unsigned int k = 0; k < 6; ++k ) { // Magic numbers
        for ( unsigned int j = 0; j < 5; ++j ) { // Magic numbers
          auto& tc = inTriggerCellDistributionGrid[iClock][ 30*i + 6*j + k ]; // Magic numbers
          if ( tc->dataValid() ) {
              tc->setIndex( (30*i) + (5*k) + j ); // Magic numbers
              triggerCellsOut.push_back( tc );
              outTriggerCellDistributionGrid[iClock-2][tc->index()] = tc;
          }
        }
      }
    }
  }
}

void HGCalTCDistribution::triggerCellDistribution4( const HGCalTriggerCellSAShrPtrCollection& triggerCellsIn ) const {
  const unsigned int stepLatency = config_.getStepLatency( Dist4 );
  for ( auto& lCell : triggerCellsIn ) {
    lCell->addLatency( stepLatency );
    unsigned int sector = int( lCell->index() / 5 ); // Magic numbers
    for ( int iSortKey = 5; iSortKey >= 0; --iSortKey ) {  // Magic numbers
      if ( int(lCell->phi() % 1944) > int( ( 18 * iSortKey ) + ( 108 * sector ) ) ) { // Magic numbers
        lCell->setSortKey(iSortKey);
        break;
      }
    }
  }
}

void HGCalTCDistribution::triggerCellDistribution5( const HGCalTriggerCellSAShrPtrCollections& inTriggerCellDistributionGrid, HGCalTriggerCellSAPtrCollection& triggerCellsOut ) const {
  const unsigned int latency = config_.getLatencyUpToAndIncluding( Dist5 );

  // Dummy distribution grid?  After writing the runDistServers function to avoid duplicating identical code from triggerCellDistribution2, I realised the second possible use case of runDistServers isn't completely identical i.e. the triggerCellDistributionGrid isn't used/set
  HGCalTriggerCellSAShrPtrCollections outTriggerCellDistributionGrid;
  HGCalTriggerCellSAShrPtrCollection triggerCellsOutTemp;
  runDistServers(inTriggerCellDistributionGrid,
                  outTriggerCellDistributionGrid,
                  triggerCellsOutTemp,
                  latency,
                  18, 5, 6, 4, false); // Magic numbers

  // Temp fix whilst moving from shared to unique ptr
  triggerCellsOut.clear();
  for ( const auto& tc : triggerCellsOutTemp ) {
    triggerCellsOut.push_back( make_unique<HGCalTriggerCell>(*tc));
  }

}

void HGCalTCDistribution::initializeTriggerCellDistGrid( HGCalTriggerCellSAShrPtrCollections& grid, unsigned int nX, unsigned int nY ) const {
  for (unsigned int iX = 0; iX < nX; ++iX ) {
    HGCalTriggerCellSAShrPtrCollection temp;
    for (unsigned int iY = 0; iY < nY; ++iY ) {
      temp.emplace_back(make_shared<HGCalTriggerCell>());
    }
    grid.push_back(temp);
  }
}

void HGCalTCDistribution::runDistServers( const HGCalTriggerCellSAShrPtrCollections& gridIn,
                      HGCalTriggerCellSAShrPtrCollections& gridOut,
                      HGCalTriggerCellSAShrPtrCollection& tcsOut,
                      unsigned int latency,
                      unsigned int nDistServers,
                      unsigned int nInputs,
                      unsigned int nOutputs,
                      unsigned int nInterleave,
                      bool setOutputGrid ) const {
  vector< DistServer > distServers(nDistServers, DistServer(nInputs, nOutputs, nInterleave));

  for ( unsigned int iClock = 0; iClock < config_.cClocks(); ++iClock ) {
    for ( unsigned int iDistServer = 0; iDistServer < nDistServers; ++iDistServer ) {
      auto first = gridIn[iClock].cbegin() + nInputs*iDistServer;
      auto last = gridIn[iClock].cbegin() + nInputs*(iDistServer+1);
      HGCalTriggerCellSAShrPtrCollection inCells(first, last);
      HGCalTriggerCellSAShrPtrCollection lCells = distServers[iDistServer].clock(inCells);

      for ( unsigned int iOutput = 0; iOutput<lCells.size(); ++iOutput ) {
        auto& tc = lCells[iOutput];
        if ( tc->dataValid() ){
          tc->setIndex(nOutputs*iDistServer+iOutput);
          tc->setClock(iClock+latency);

          tcsOut.push_back(tc);
          if ( setOutputGrid ) {
            gridOut[iClock][tc->index()] = tc;
          }
        }
      }
    }
  }
}
