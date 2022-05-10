#include "L1Trigger/L1THGCal/interface/backend_emulator/HGCalTCDistribution_SA.h"

using namespace std;
using namespace l1thgcfirmware;

HGCalTCDistribution::HGCalTCDistribution(ClusterAlgoConfig& config) : config_(config) {}

void HGCalTCDistribution::runTriggerCellDistribution(const HGCalTriggerCellSAPtrCollections& triggerCellsIn, HGCalTriggerCellSAPtrCollection& triggerCellsOut ) const {
 
    triggerCellsOut = triggerCellInput( triggerCellsIn );
    triggerCellDistribution0( triggerCellsOut );
    HGCalTriggerCellSAPtrCollections tcDistGrid1 = triggerCellDistribution1( triggerCellsOut );
    HGCalTriggerCellSAPtrCollections tcDistGrid2 = triggerCellDistribution2( triggerCellsOut, tcDistGrid1 );
    HGCalTriggerCellSAPtrCollections tcDistGrid3 = triggerCellDistribution3( triggerCellsOut, tcDistGrid2 );
    triggerCellDistribution4( triggerCellsOut );
    triggerCellDistribution5( triggerCellsOut, tcDistGrid3 );

}

HGCalTriggerCellSAPtrCollection HGCalTCDistribution::triggerCellInput( const HGCalTriggerCellSAPtrCollections& inputs ) const {

  HGCalTriggerCellSAPtrCollection triggerCellsIn;
  for (unsigned int iFrame = 0; iFrame < inputs.size(); ++iFrame ) {
    for (unsigned int iInput = 0; iInput < inputs[iFrame].size(); ++iInput ) {
      auto& tc = inputs[iFrame][iInput];
      tc->setIndex(iInput);
      tc->setClock(iFrame+1);
      if ( tc->dataValid() ) {
        triggerCellsIn.push_back(tc);
      }
    }
  }
  return triggerCellsIn;
}

void HGCalTCDistribution::triggerCellDistribution0( HGCalTriggerCellSAPtrCollection& triggerCellsIn ) const {
  for ( auto& tc : triggerCellsIn ) {
    unsigned int newIndex = tc->index() + int( tc->index() / 24 ); // Magic numbers
    tc->setIndex( newIndex );
  }
}

HGCalTriggerCellSAPtrCollections HGCalTCDistribution::triggerCellDistribution1( HGCalTriggerCellSAPtrCollection& triggerCellsIn ) const {

  HGCalTriggerCellSAPtrCollections triggerCellDistributionGrid;
  initializeTriggerCellDistGrid( triggerCellDistributionGrid, config_.cClocks(), config_.cInputs2() );

  const unsigned int stepLatency = config_.getStepLatency( Dist1 );
  for ( auto& tc : triggerCellsIn ) {
    tc->addLatency( stepLatency );
    unsigned int sector = int( tc->index() / 25 ); // Magic numbers
    triggerCellDistributionGrid[tc->clock()-2][tc->index()] = tc; // Magic numbers
    for ( int iSortKey = 5; iSortKey >= 0; --iSortKey ) { // Magic numbers
      // Split each 60 degree sector into 6 phi region
      // Sort key is index of small phi region
      if ( int( tc->phi() % 1944 ) > int( 108 * iSortKey + 648 * sector ) ) { // Magic numbers
        tc->setSortKey( iSortKey );
        break;
      }
    }
  }

  return triggerCellDistributionGrid;
}

HGCalTriggerCellSAPtrCollections HGCalTCDistribution::triggerCellDistribution2( HGCalTriggerCellSAPtrCollection& triggerCellsIn, HGCalTriggerCellSAPtrCollections& inTriggerCellDistributionGrid ) const {

  const unsigned int latency = config_.getLatencyUpToAndIncluding( Dist2 );

  HGCalTriggerCellSAPtrCollection triggerCells;
  HGCalTriggerCellSAPtrCollections triggerCellDistributionGrid;
  initializeTriggerCellDistGrid( triggerCellDistributionGrid, config_.cClocks(), config_.cInt() );

  runDistServers(inTriggerCellDistributionGrid,
                  triggerCellDistributionGrid,
                  triggerCells,
                  latency,
                  15, 5, 6, 4, true); // Magic numbers

  triggerCellsIn = triggerCells;
  return triggerCellDistributionGrid;
}

HGCalTriggerCellSAPtrCollections HGCalTCDistribution::triggerCellDistribution3( HGCalTriggerCellSAPtrCollection& triggerCellsIn, HGCalTriggerCellSAPtrCollections& inTriggerCellDistributionGrid ) const {
  HGCalTriggerCellSAPtrCollection triggerCells;
  HGCalTriggerCellSAPtrCollections triggerCellDistributionGrid;
  initializeTriggerCellDistGrid( triggerCellDistributionGrid, config_.cClocks(), config_.cInt() );
  for ( unsigned int iClock = 0; iClock < config_.cClocks(); ++iClock ) {
    for ( unsigned int i = 0; i < 3; ++i ) { // Magic numbers
      for ( unsigned int k = 0; k < 6; ++k ) { // Magic numbers
        for ( unsigned int j = 0; j < 5; ++j ) { // Magic numbers
          auto& tc = inTriggerCellDistributionGrid[iClock][ 30*i + 6*j + k ]; // Magic numbers
          if ( tc->dataValid() ) {
              tc->setIndex( (30*i) + (5*k) + j ); // Magic numbers
              triggerCells.push_back( tc );
              triggerCellDistributionGrid[iClock-2][tc->index()] = tc;
          }
        }
      }
    }
  }
  triggerCellsIn = triggerCells;

  return triggerCellDistributionGrid;
}

void HGCalTCDistribution::triggerCellDistribution4( HGCalTriggerCellSAPtrCollection& triggerCellsIn ) const {
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

void HGCalTCDistribution::triggerCellDistribution5( HGCalTriggerCellSAPtrCollection& triggerCellsIn, HGCalTriggerCellSAPtrCollections& inTriggerCellDistributionGrid ) const {
  const unsigned int latency = config_.getLatencyUpToAndIncluding( Dist5 );

  HGCalTriggerCellSAPtrCollection triggerCells;
  // Dummy distribution grid?  After writing the runDistServers function to avoid duplicating identical code from triggerCellDistribution2, I realised the second possible use case of runDistServers isn't completely identical i.e. the triggerCellDistributionGrid isn't used/set
  HGCalTriggerCellSAPtrCollections triggerCellDistributionGrid;
  initializeTriggerCellDistGrid( triggerCellDistributionGrid, 1, 1 );

  runDistServers(inTriggerCellDistributionGrid,
                  triggerCellDistributionGrid,
                  triggerCells,
                  latency,
                  18, 5, 6, 4, false); // Magic numbers

  triggerCellsIn = triggerCells;
}

void HGCalTCDistribution::initializeTriggerCellDistGrid( HGCalTriggerCellSAPtrCollections& grid, unsigned int nX, unsigned int nY ) const {
  for (unsigned int iX = 0; iX < nX; ++iX ) {
    HGCalTriggerCellSAPtrCollection temp;
    for (unsigned int iY = 0; iY < nY; ++iY ) {
      temp.emplace_back(make_shared<HGCalTriggerCell>());
    }
    grid.push_back(temp);
  }
}

void HGCalTCDistribution::runDistServers( const HGCalTriggerCellSAPtrCollections& gridIn,
                      HGCalTriggerCellSAPtrCollections& gridOut,
                      HGCalTriggerCellSAPtrCollection& tcsOut,
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
      HGCalTriggerCellSAPtrCollection inCells(first, last);
      HGCalTriggerCellSAPtrCollection lCells = distServers[iDistServer].clock(inCells);

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
