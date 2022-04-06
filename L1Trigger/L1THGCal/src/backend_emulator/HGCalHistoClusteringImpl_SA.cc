#include "L1Trigger/L1THGCal/interface/backend_emulator/HGCalHistoClusteringImpl_SA.h"
#include "L1Trigger/L1THGCal/interface/backend_emulator/DistServer.h"
#include "L1Trigger/L1THGCal/interface/backend_emulator/CentroidHelper.h"

#include <math.h>
#include <map>
#include <vector>
#include <algorithm>
#include <cmath>
#include <iostream>

using namespace std;
using namespace l1thgcfirmware;
HGCalHistoClusteringImplSA::HGCalHistoClusteringImplSA( ClusterAlgoConfig& config ) : config_(config) {}

void HGCalHistoClusteringImplSA::runAlgorithm(HGCalTriggerCellSAPtrCollections& inputs, HGCalTriggerCellSAPtrCollection& clusteredTCs, HGCalTriggerCellSAPtrCollection& unclusteredTCs, CentroidHelperPtrCollection& prioritizedMaxima, CentroidHelperPtrCollection& readoutFlags, HGCalClusterSAPtrCollection& clusterSums ) const {
  // config_.printConfiguration();

  HGCalTriggerCellSAPtrCollection triggerCellsIn = triggerCellInput( inputs );

  // TC Distribution
  triggerCellDistribution0( triggerCellsIn );
  HGCalTriggerCellSAPtrCollections tcDistGrid1 = triggerCellDistribution1( triggerCellsIn );
  HGCalTriggerCellSAPtrCollections tcDistGrid2 = triggerCellDistribution2( triggerCellsIn, tcDistGrid1 );
  HGCalTriggerCellSAPtrCollections tcDistGrid3 = triggerCellDistribution3( triggerCellsIn, tcDistGrid2 );
  triggerCellDistribution4( triggerCellsIn );
  triggerCellDistribution5( triggerCellsIn, tcDistGrid3 );

  // Histogram
  HGCalHistogramCellSAPtrCollection histoCells = triggerCellToHistogramCell( triggerCellsIn );
  HGCalHistogramCellSAPtrCollection histogram = makeHistogram( histoCells );

  // Smearing
  smearing1D( histogram );
  areaNormalization( histogram );
  smearing2D( histogram );

  //Maxima finding
  thresholdMaximaFinder( histogram );
  localMaximaFinder( histogram );
  calculateAveragePosition( histogram );

  // Clustering
  clusterizer( triggerCellsIn, histogram, clusteredTCs, unclusteredTCs, prioritizedMaxima, readoutFlags );

  // Cluster properties
  HGCalClusterSAPtrCollection protoClusters = triggerCellToCluster( clusteredTCs );
  HGCalClusterSAPtrCollection clusterAccumulation;
  clusterSum( protoClusters, readoutFlags, clusterAccumulation, clusterSums );
  clusterProperties(clusterSums);

}

HGCalTriggerCellSAPtrCollection HGCalHistoClusteringImplSA::triggerCellInput( HGCalTriggerCellSAPtrCollections& inputs ) const {

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

void HGCalHistoClusteringImplSA::triggerCellDistribution0( HGCalTriggerCellSAPtrCollection& triggerCellsIn ) const {
  for ( auto& tc : triggerCellsIn ) {
    unsigned int newIndex = tc->index() + int( tc->index() / 24 ); // Magic numbers
    tc->setIndex( newIndex );
  }
}

HGCalTriggerCellSAPtrCollections HGCalHistoClusteringImplSA::triggerCellDistribution1( HGCalTriggerCellSAPtrCollection& triggerCellsIn ) const {

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

HGCalTriggerCellSAPtrCollections HGCalHistoClusteringImplSA::triggerCellDistribution2( HGCalTriggerCellSAPtrCollection& triggerCellsIn, HGCalTriggerCellSAPtrCollections& inTriggerCellDistributionGrid ) const {

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

HGCalTriggerCellSAPtrCollections HGCalHistoClusteringImplSA::triggerCellDistribution3( HGCalTriggerCellSAPtrCollection& triggerCellsIn, HGCalTriggerCellSAPtrCollections& inTriggerCellDistributionGrid ) const {
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

void HGCalHistoClusteringImplSA::triggerCellDistribution4( HGCalTriggerCellSAPtrCollection& triggerCellsIn ) const {
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

void HGCalHistoClusteringImplSA::triggerCellDistribution5( HGCalTriggerCellSAPtrCollection& triggerCellsIn, HGCalTriggerCellSAPtrCollections& inTriggerCellDistributionGrid ) const {
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

HGCalHistogramCellSAPtrCollection HGCalHistoClusteringImplSA::triggerCellToHistogramCell( HGCalTriggerCellSAPtrCollection& triggerCellsIn ) const {

  const unsigned int latency = config_.getStepLatency( TcToHc );

  HGCalHistogramCellSAPtrCollection histoCells;
  for ( auto& tc : triggerCellsIn ) {
    auto hc = make_shared<HGCalHistogramCell>( tc->clock() + latency,
                                               tc->index(),
                                               tc->energy(),
                                               tc->phi(),
                                               tc->rOverZ(),
                                               1,
                                               int( ( tc->rOverZ() - config_.rOverZHistOffset() )/ config_.rOverZBinSize() ) // Magic numbers 
                                              );
    tc->setClock( hc->clock() );
    tc->setSortKey( hc->sortKey() );
    histoCells.push_back( hc );
  }

  return histoCells;
}

HGCalHistogramCellSAPtrCollection HGCalHistoClusteringImplSA::makeHistogram( HGCalHistogramCellSAPtrCollection histogramCells ) const {

  HGCalHistogramCellSAPtrCollection histogram;
  const unsigned int latency = config_.getLatencyUpToAndIncluding( Hist );
  for ( unsigned int iRow = 0; iRow < config_.cRows(); ++iRow ) {
    for ( unsigned int iColumn = 0; iColumn < config_.cColumns(); ++iColumn ) {

      auto hc = make_shared<HGCalHistogramCell>( latency, iColumn, iRow );
      histogram.push_back( hc );
    }
  }

  for ( const auto& hc : histogramCells ) {
    const unsigned int binIndex = config_.cColumns() * hc->sortKey() + hc->index() ;
    *histogram.at( binIndex ) += *hc;
  }

  return histogram;

}

void HGCalHistoClusteringImplSA::smearing1D( HGCalHistogramCellSAPtrCollection& histogram ) const {

  HGCalHistogramCellSACollection lHistogram;
  for ( unsigned int iBin = 0; iBin < histogram.size(); ++iBin ) {
    lHistogram.emplace_back( *histogram.at(iBin) );
  }

  const unsigned int stepLatency = config_.getStepLatency( Step::Smearing1D );
  for ( unsigned int iBin = 0; iBin < lHistogram.size(); ++iBin ) {
    auto& hc = histogram.at(iBin);
    hc->addLatency( stepLatency );

    const unsigned int col = hc->index();
    const unsigned int row = hc->sortKey();
    const unsigned int binIndex = config_.cColumns() * row + col;
    unsigned int scale = 1;
    int width = config_.kernelWidth( row );
    unsigned int offset = 1;
    while ( width > 0 ) {
      shared_ptr<HGCalHistogramCell> l1 = make_shared<HGCalHistogramCell>(HGCalHistogramCell());
      shared_ptr<HGCalHistogramCell> l2 = make_shared<HGCalHistogramCell>(HGCalHistogramCell());
      shared_ptr<HGCalHistogramCell> r1 = make_shared<HGCalHistogramCell>(HGCalHistogramCell());
      shared_ptr<HGCalHistogramCell> r2 = make_shared<HGCalHistogramCell>(HGCalHistogramCell());


      if ( width >= 2 ) {
        if ( int(col - offset - 1)  >= 0 ) {
          l2 = make_shared<HGCalHistogramCell>(lHistogram[binIndex - offset - 1]/4); // Magic numbers (?)
        }
        if ( int(col + offset + 1) <= int(config_.cColumns()-1) ) {
          r2 = make_shared<HGCalHistogramCell>(lHistogram[binIndex + offset + 1]/4); // Magic numbers (?)
        }
      }
      
      if ( int(col - offset)  >= 0 ) {
        l1 = make_shared<HGCalHistogramCell>(lHistogram[binIndex - offset]/2); // Magic numbers (?)
      }
      if ( int(col + offset)  <= int(config_.cColumns()-1) ) {
        r1 = make_shared<HGCalHistogramCell>(lHistogram[binIndex + offset]/2); // Magic numbers (?)
      }
      *hc += ( ( *l2 + *l1 ) / scale + ( *r2 + *r1 ) / scale );
      scale *= 4; // Magic numbers
      width -= 2; // Magic numbers
      offset += 2; // Magic numbers
    }
  }
}

void HGCalHistoClusteringImplSA::areaNormalization( HGCalHistogramCellSAPtrCollection& histogram ) const {
  const unsigned int stepLatency = config_.getStepLatency( NormArea );
  for ( unsigned int iBin = 0; iBin < histogram.size(); ++iBin ) {
    HGCalHistogramCell& hc = *histogram.at(iBin);
    hc.addLatency( stepLatency );
    hc *= config_.areaNormalization(hc.sortKey());
  }
}

void HGCalHistoClusteringImplSA::smearing2D( HGCalHistogramCellSAPtrCollection& histogram ) const {
  HGCalHistogramCellSACollection lHistogram;
  for ( unsigned int iBin = 0; iBin < histogram.size(); ++iBin ) {
    lHistogram.emplace_back( *histogram.at(iBin) );
  }
  const unsigned int stepLatency = config_.getStepLatency( Step::Smearing2D );
  for ( unsigned int iBin = 0; iBin < lHistogram.size(); ++iBin ) {
    auto& hc = histogram.at(iBin);
    hc->addLatency( stepLatency );

    const unsigned int col = hc->index();
    const int row = hc->sortKey();
    const unsigned int binIndex = config_.cColumns() * row + col;
    if ( row - 1 >= 0 ) {
      *hc += (lHistogram[binIndex - config_.cColumns()] / 2 ); // Magic numbers (?)
    }
    if ( row + 1 <= int(config_.cRows()-1) ) {
      *hc += (lHistogram[binIndex + config_.cColumns()] / 2 ); // Magic numbers (?)
    }
  }
}

void HGCalHistoClusteringImplSA::thresholdMaximaFinder( HGCalHistogramCellSAPtrCollection& histogram ) const {
  const unsigned int stepLatency = config_.getStepLatency( Maxima2D );

  // std::cout << "Histogram for maxima finding" << std::endl;
  // printHistogram( histogram );

  for ( auto& hc : histogram ) {
    hc->addLatency( stepLatency );
    if ( hc->S() <= config_.thresholdMaxima( hc->sortKey() ) ) {
      hc->setS(0);
      hc->setX(0);
      hc->setY(0);
      hc->setN(0);
    }
  }
  // std::cout << "Threshold Maxima" << std::endl;
  // printHistogram( histogram );
}

// Temporary simulation of local maxima finder
// Not an emulation of any firmware
void HGCalHistoClusteringImplSA::localMaximaFinder( l1thgcfirmware::HGCalHistogramCellSAPtrCollection& histogram ) const {
  const std::vector<unsigned> maximaWidths{ 6, 5, 5, 5, 4, 4, 4, 4, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1 };


  for ( auto& hc : histogram ) {
    if ( hc->S() > 0 ) {
      const int colRef = hc->index();
      const int rowRef = hc->sortKey();
      bool isMaxima = true;
      // std::cout << "Got a histo cell : " << hc->index() << " " << hc->sortKey() << " " << hc->S() << std::endl;
      // std::cout << "Phi range : " << maximaWidths.at(hc->sortKey()) << std::endl;
      const int phiRange = maximaWidths.at(hc->sortKey());
      for ( int colOffset = -1 * phiRange; colOffset <= phiRange; ++colOffset ) {
        const int col = colRef + colOffset;
        if ( col < 0 || col >= (int) config_.cColumns() ) continue;
        for ( int rowOffset = -1; rowOffset <= 1; ++ rowOffset ) {
          const int row = rowRef + rowOffset;
          if ( row < 0 || row >= (int) config_.cRows() ) continue;
          const unsigned int binIndex = config_.cColumns() * row + col;
          const auto& bin = histogram.at(binIndex);
          // std::cout << "Comparing to : " << col << " " << row << " " << bin->S() << std::endl;

          if ( colOffset == 0 && rowOffset == 0 ) continue;
          else if ( ( col < colRef ) ||
               ( col == colRef && rowOffset == -1 ) ||
               ( col == colRef+1 && rowOffset == -1 ) ) {
            if ( !(hc->S() >= bin->S() ) ) isMaxima = false;
          }
          else {
            if ( !(hc->S() > bin->S() ) ) isMaxima = false;
          }
        }
      }
      if ( !isMaxima ) {
        hc->setS(0);
        hc->setX(0);
        hc->setY(0);
        hc->setN(0);
      }
    }
  }
  // std::cout << "Local Maxima" << std::endl;
  // printHistogram( histogram );
}

void HGCalHistoClusteringImplSA::calculateAveragePosition( HGCalHistogramCellSAPtrCollection& histogram ) const {
  const unsigned int stepLatency = config_.getStepLatency( CalcAverage );
  for ( auto& hc : histogram ) {
    hc->addLatency( stepLatency );
    if ( hc->N() > 0 ) {
      unsigned int inv_N = int( round(1.0 * 0x1FFFF / hc->N() ) ); // Magic numbers
      hc->setX( ( hc->X() * inv_N ) >> 17 ); //Magic numbers
      hc->setY( ( hc->Y() * inv_N ) >> 17 ); //Magic numbers
    }
  }
}

void HGCalHistoClusteringImplSA::clusterizer( HGCalTriggerCellSAPtrCollection& triggerCellsIn, HGCalHistogramCellSAPtrCollection& histogram, HGCalTriggerCellSAPtrCollection& clusteredTriggerCellsOut, HGCalTriggerCellSAPtrCollection& unclusteredTriggerCellsOut, CentroidHelperPtrCollection& prioritizedMaxima, CentroidHelperPtrCollection& readoutFlagsOut ) const {
  
  unsigned int seedCounter = 0;
  vector< CentroidHelperPtrCollection > fifos( 18, CentroidHelperPtrCollection() ); // Magic numbers
  vector<unsigned int> clock( config_.cColumns(), config_.clusterizerMagicTime() );
  CentroidHelperPtrCollection latched( 18+1, make_shared<CentroidHelper>() ); // Magic numbers

  HGCalTriggerCellSAPtrCollections clusteredTriggerCells( config_.cColumns(), HGCalTriggerCellSAPtrCollection() );
  HGCalTriggerCellSAPtrCollections unclusteredTriggerCells( config_.cColumns(), HGCalTriggerCellSAPtrCollection() );
  CentroidHelperPtrCollections readoutFlags( config_.cColumns(), CentroidHelperPtrCollection() );

  HGCalTriggerCellSAPtrCollectionss triggerCellBuffers( config_.cColumns(), HGCalTriggerCellSAPtrCollections( config_.cRows(), HGCalTriggerCellSAPtrCollection() ) );
  for (const auto& tc : triggerCellsIn ) {
    triggerCellBuffers.at( tc->index() ).at( tc->sortKey() ).push_back( tc );
  }

  for ( unsigned int iRow = 0; iRow < config_.cRows(); ++iRow ) {
    for ( unsigned int j = 0; j < 4; ++j ) { // Magic numbers
      for ( unsigned int k = 0; k < 18; ++k ) { // Magic numbers
        unsigned int col = 18 + (4*k) + j; // Magic numbers
        const auto& cell = histogram.at( config_.cColumns() * iRow + col );
        if ( cell->S() > 0 ) {
          auto ch = make_shared<CentroidHelper>( cell->clock() + 1 + j,
                                                4*k + j, // Magic numbers
                                                cell->index(),
                                                cell->sortKey(),
                                                cell->S(),
                                                cell->X(),
                                                cell->Y(),
                                                true
                                                );
          fifos[k].push_back( ch );
          ++seedCounter;
        }
      }
    }
  }

  while ( seedCounter > 0 ) {
    for ( unsigned int i = 0; i < 18; ++i ) { // Magic numbers
      if ( !latched[i]->dataValid() ) {
        if ( fifos[i].size() > 0 ) {
          latched[i] = fifos[i][0];
          fifos[i].erase(fifos.at(i).begin());
        }
      }
    }

    CentroidHelperPtrCollection accepted( 20, make_shared<CentroidHelper>() ); // Magic numbers
    CentroidHelperPtrCollection lastLatched( latched );

    for ( unsigned int i = 0; i < 18; ++i ) { // Magic numbers
      // Different implementation to python emulator
      // For i=0, i-1=-1, which would give the last element of lastLatched in python, but is out of bounds in C++
      // Similar for i=17
      // Need to find out intended behaviour
      bool deltaMinus = (i>0) ? ( lastLatched[i]->column() - lastLatched[i-1]->column() ) > 6 : false ; // Magic numbers
      bool deltaPlus = (i<17) ? ( lastLatched[i+1]->column() - lastLatched[i]->column() ) > 6 : false ; // Magic numbers

      bool compareEMinus = (i>0) ? ( lastLatched[i]->energy() > lastLatched[i-1]->energy() ) : false; // Magic numbers
      bool compareEPlus = (i<17) ? ( lastLatched[i]->energy() >= lastLatched[i+1]->energy() ) : false; // Magic numbers

      if ( lastLatched[i]->dataValid() ) {
        // Similar out of bounds issue here
        // if ( ( !lastLatched[i+1]->dataValid() || compareEPlus || deltaPlus ) && ( !lastLatched[i-1]->dataValid() || compareEMinus || deltaMinus ) ) {

        bool accept = true;
        if ( lastLatched.size() > i+1 ) {
          if ( !lastLatched[i+1]->dataValid() || compareEPlus || deltaPlus ) {
            accept = true;
          }
          else {
            accept = false;
          }
        }

        if ( i > 0 ) {
          if ( !lastLatched[i-1]->dataValid() || compareEMinus || deltaMinus  ) {
            accept = true;
          }
          else {
            accept = false;
          }
        }

        if ( accept ) {
          accepted[i] = latched[i];
          latched[i] = make_shared<CentroidHelper>();
          --seedCounter;
        }
      }
    }

    CentroidHelperPtrCollection output( config_.cColumns(), make_shared<CentroidHelper>() );
    for ( const auto& a : accepted ) {
      if ( a->dataValid() ) {
        for ( unsigned int iCol = a->column() - 3; iCol < a->column() + 4; ++iCol ) { // Magic numbers
          clock[iCol] = clock[a->column()];
          output[iCol] = make_shared<CentroidHelper>(*a);
          output[iCol]->setIndex( iCol );
          output[iCol]->setClock( clock[iCol] );
          prioritizedMaxima.push_back( output[iCol] );
        }
      }
    }

    for ( const auto& a : accepted ) {
      if ( a->dataValid() ) {
        unsigned int dR2Cut = 20000; // Magic numbers
        unsigned int T=0; // Magic numbers

        for ( unsigned int iCol = a->column() - 3; iCol < a->column() + 4; ++iCol ) { // Magic numbers
          clock[ iCol ] += 8; // Magic numbers
          for ( int k = -2; k < 3; ++k ) { // Magic numbers
            int row = a->row() + k;
            if ( row < 0 ) continue;
            if ( row >= int(config_.cRows()) ) continue; // Not in python emulator, but required to avoid out of bounds access
            if ( triggerCellBuffers[iCol][row].size() == 0 ) {
              clock[iCol] += 1 ;
              continue;
            }

            for ( auto& tc : triggerCellBuffers[iCol][row] ) {
              clock[iCol] += 1 ;

              unsigned int r1 = tc->rOverZ();
              unsigned int r2 = a->Y();
              int dR = r1 - r2;
              int dPhi = tc->phi() - a->X();
              unsigned int dR2 = dR * dR;
              unsigned int cosTerm = ( abs(dPhi) > config_.nBinsCosLUT() ) ? 2047 : config_.cosLUT( abs(dPhi) ); // Magic numbers
              dR2 += int( r1 * r2 / pow(2,7) ) * cosTerm / pow(2,10); // Magic numbers
              tc->setClock( clock[iCol] + 1 );
              if ( clock[iCol] > T ) T = clock[iCol];

              if ( dR2 < dR2Cut ) {
                // std::cout << "Clustered a TC to a seed : " << tc->energy() << " " << a->row() << " " << a->column() << " " << a->energy() << std::endl;
                clusteredTriggerCells[iCol].push_back(tc);
              }
              else {
                unclusteredTriggerCells[iCol].push_back(tc);
              }
            }
          }

          for ( const auto& tc : clusteredTriggerCells[iCol] ) {
            auto tcMatch = std::find_if(triggerCellBuffers[iCol][tc->sortKey()].begin(), triggerCellBuffers[iCol][tc->sortKey()].end(), [&](const HGCalTriggerCellSAPtr tcToMatch) {
              bool isMatch = tc->index() == tcToMatch->index() &&
                             tc->rOverZ() == tcToMatch->rOverZ() &&
                             tc->layer() == tcToMatch->layer() &&
                             tc->energy() == tcToMatch->energy() &&
                             tc->phi() == tcToMatch->phi() &&
                             tc->sortKey() == tcToMatch->sortKey() &&
                             tc->deltaR2() == tcToMatch->deltaR2() &&
                             tc->dX() == tcToMatch->dX() &&
                             tc->Y() == tcToMatch->Y() &&
                             tc->frameValid() == tcToMatch->frameValid() &&
                             tc->dataValid() == tcToMatch->dataValid() &&
                             tc->clock() == tcToMatch->clock();
              return isMatch;
            });

            if ( tcMatch != triggerCellBuffers[iCol][tc->sortKey()].end() ) {
              triggerCellBuffers[iCol][tc->sortKey()].erase(tcMatch);
            }
          }
        }
        for ( unsigned int iCol = a->column() - 3; iCol < a->column() + 4; ++iCol ) { // Magic numbers
          clock[iCol] = T+1;

          CentroidHelperPtr readoutFlag = make_shared<CentroidHelper>(T-2, iCol, true);
          if ( readoutFlag->clock() == 448 ) { // Magic numbers
            readoutFlag->setClock( readoutFlag->clock() + 1 );
          }

          readoutFlags[iCol].push_back( readoutFlag );
        }
      }
    }
  }

  for ( unsigned int i = 0; i <1000; ++i ) { // Magic numbers
    for ( unsigned int iCol = 0; iCol < config_.cColumns(); ++iCol ) {
      for ( const auto& clustered : clusteredTriggerCells[iCol] ) {
        if ( clustered->clock() == config_.clusterizerMagicTime() + i ) {
          clusteredTriggerCellsOut.push_back( clustered );
        }
      }

      for ( const auto& unclustered : unclusteredTriggerCells[iCol] ) {
        if ( unclustered->clock() == config_.clusterizerMagicTime() + i ) {
          unclusteredTriggerCellsOut.push_back( unclustered );
        }
      }

      for ( const auto& readoutFlag : readoutFlags[iCol] ) {
        if ( readoutFlag->clock() == config_.clusterizerMagicTime() + i ) {
          readoutFlagsOut.push_back( readoutFlag );
        }
      }
    }
  }

  // std::cout << "Output from Clusterizer" << std::endl;
  // std::cout << "Number of clustered TCs : " << clusteredTriggerCellsOut.size() << std::endl;
  // for ( const auto& tc : clusteredTriggerCellsOut ) {
  //   std::cout << tc->clock() << " " << tc->index() << " " << tc->rOverZ() << " " << tc->layer() << " " << tc->energy() << " " << tc->phi() << " " << tc->sortKey() << " " << tc->deltaR2() << " " << tc->dX() << " " << tc->Y() << " " << tc->dataValid() << std::endl;
  // }
  // std::cout << "Number of unclustered TCs : " << unclusteredTriggerCellsOut.size() << std::endl;
  // std::cout << "Number of readoutFlags : " << readoutFlagsOut.size() << std::endl;
  // for ( const auto& f : readoutFlagsOut ) {
  //   std::cout << f->clock() << " " << f->index() << " " << f->column() << " " << f->row() << " " << f->energy() << " " << f->X() << " " << f->Y() << " " << f->dataValid() << std::endl;
  // }
}

HGCalClusterSAPtrCollection HGCalHistoClusteringImplSA::triggerCellToCluster( HGCalTriggerCellSAPtrCollection& clusteredTriggerCells ) const {

  const unsigned int stepLatency = config_.getStepLatency( TriggerCellToCluster );

  HGCalClusterSAPtrCollection protoClusters;

  for ( const auto& tc : clusteredTriggerCells ) {

    auto cluster = make_shared<HGCalCluster>( tc->clock() + stepLatency,
                                              tc->index(),
                                              true, true
                                            );

    // Cluster from single TC
    // Does this ever happen?
    if ( tc->deltaR2() >= 25000 ) { // Magic numbers
      protoClusters.push_back( cluster );
      continue;
    }

    unsigned long int s_TC_W = ( int( tc->energy() / 4 ) == 0 ) ? 1 : tc->energy() / 4;
    unsigned long int s_TC_Z = config_.depth( tc->layer() );

    unsigned int triggerLayer = config_.triggerLayer( tc->layer() );
    unsigned int s_E_EM = ( (  ( (unsigned long int) tc->energy() * config_.layerWeight_E_EM( triggerLayer ) ) + config_.correction() ) >> 18 );
    if ( s_E_EM > config_.saturation() ) s_E_EM = config_.saturation();



    unsigned int s_E_EM_core = ( ( (unsigned long int) tc->energy() * config_.layerWeight_E_EM_core( triggerLayer ) + config_.correction() ) >> 18 );
    if ( s_E_EM_core > config_.saturation() ) s_E_EM_core = config_.saturation();

    // Alternative constructor perhaps?
    cluster->set_n_tc( 1 ); // Magic numbers
    cluster->set_n_tc_w( 1 ); // Magic numbers
    
    cluster->set_e( ( config_.layerWeight_E( triggerLayer ) == 1 ) ? tc->energy() : 0  );
    cluster->set_e_h_early( ( config_.layerWeight_E_H_early( triggerLayer ) == 1 ) ? tc->energy() : 0  );

    cluster->set_e_em( s_E_EM );
    cluster->set_e_em_core( s_E_EM_core );

    cluster->set_w( s_TC_W );
    cluster->set_w2( s_TC_W * s_TC_W );

    cluster->set_wz( s_TC_W * s_TC_Z );
    cluster->set_weta( 0 );
    cluster->set_wphi( s_TC_W * tc->phi() );
    cluster->set_wroz( s_TC_W * tc->rOverZ() );

    cluster->set_wz2( s_TC_W * s_TC_Z * s_TC_Z );
    cluster->set_weta2( 0 );
    cluster->set_wphi2( s_TC_W * tc->phi() * tc->phi() );
    cluster->set_wroz2( s_TC_W * tc->rOverZ() * tc->rOverZ() );

    cluster->set_layerbits( cluster->layerbits() | ( ( (unsigned long int) 1) << ( 36 - triggerLayer ) ) ); // Magic numbers
    cluster->set_sat_tc( cluster->e() == config_.saturation() || cluster->e_em() == config_.saturation() );
    cluster->set_shapeq(1);

    protoClusters.push_back( cluster );
  }

  // std::cout << "Output from triggerCellToCluster" << std::endl;
  // std::cout << "Protoclusters : " << protoClusters.size() << std::endl;
  // for ( const auto& pclus : protoClusters ) {

  //     std::cout << pclus->clock() << " " << pclus->index() << " " << pclus->n_tc() << " " << pclus->e() << " " << pclus->e_em() << " " << pclus->e_em_core() << " " << pclus->e_h_early() << " " << pclus->w() << " " << pclus->n_tc_w() << " " << pclus->weta2() << " " << pclus->wphi2() << std::endl;
  // }
  return protoClusters;
}

void HGCalHistoClusteringImplSA::clusterSum( HGCalClusterSAPtrCollection& protoClusters, CentroidHelperPtrCollection& readoutFlags, HGCalClusterSAPtrCollection& clusterAccumulation, HGCalClusterSAPtrCollection& clusterSums ) const {

  HGCalClusterSAPtrCollections protoClustersPerColumn( config_.cColumns(), HGCalClusterSAPtrCollection() );
  vector<unsigned int> clock( config_.cColumns(), 0 );
  for ( const auto& protoCluster : protoClusters ) {
    protoClustersPerColumn.at( protoCluster->index() ).push_back( protoCluster );
  }

  map<unsigned int, HGCalClusterSAPtr> sums;

  for ( const auto& flag : readoutFlags ) {
    auto accumulator = make_shared<HGCalCluster>( 0,
                                                  0,
                                                  true, true
                                                );
    flag->setClock( flag->clock() + 23 ); // Magic numbers

    for ( const auto& protoCluster : protoClustersPerColumn.at( flag->index() ) ) {
      if ( protoCluster->clock() <= clock.at( flag->index() ) ) continue;
      if ( protoCluster->clock() > flag->clock() ) continue;
      *accumulator += *protoCluster;
    }

    clock.at( flag->index() ) = flag->clock();
    accumulator->setClock( flag->clock() );
    accumulator->setIndex( flag->index() );
    accumulator->setDataValid( true );
    clusterAccumulation.push_back( accumulator );

    if ( sums.find( flag->clock() ) == sums.end() ) {
      auto sum = make_shared<HGCalCluster>( flag->clock() + 7, // Magic numbers
                                            0,
                                            true, true
                                          );
      sums[flag->clock()] = sum;
    }

    *(sums.at( flag->clock() )) += *accumulator;
  }

  for (const auto& sum: sums) {
    clusterSums.push_back( sum.second );
  }

  // std::cout << "Output from ClusterSum" << std::endl;
  // unsigned int nTCs = 0;
  // for ( const auto& c : clusterAccumulation ) {
  //   nTCs += c->n_tc();
  // }
  // std::cout << nTCs << std::endl;
  // nTCs = 0;
  //unsigned int ID = 0;
  // Print the info as the Python emu format
  //std::cout << "ID,event_ID,cluster_ID,N_TC,E,E_EM,E_EM_core,E_H_early,W,N_TC_W,W2,Wz,Weta,Wphi,Wroz,Wz2,Weta2,Wphi2,Wroz2,LayerBits,Sat_TC,ShapeQ,SortKey" << std::endl;
  //for ( const auto& c : clusterSums ) {
     //std::cout << ID << ",0," << c->index() << "," << c->n_tc() << "," << c->e() << "," << c->e_em() << "," << c->e_em_core() << "," << c->e_h_early() << "," << c->w() << "," << c->n_tc_w() << "," << c->w2() << "," << c->wz() << "," << c->weta() << "," << c->wphi() << "," << c->wroz() << "," << c->wz2() << "," << c->weta2() << "," << c->wphi2() << "," << c->wroz2() << "," << c->layerbits() << "," << c->sat_tc() << "," << c->shapeq() << ",0" << std::endl;
     //ID++;
     //nTCs += c->n_tc();
  //}
  // std::cout << nTCs << std::endl;


}

std::pair< unsigned int, unsigned int > HGCalHistoClusteringImplSA::sigma_Energy(unsigned int N_TC_W, unsigned long int Sum_W2, unsigned int Sum_W) const {
  unsigned long int N = N_TC_W*Sum_W2 - pow(Sum_W,2);
  unsigned long int D = pow(N_TC_W,2);
  double intpart;
  double frac =  modf(sqrt(N/D),&intpart)*pow(2,1);
  return { (unsigned int)intpart, (unsigned int)frac };
}

std::pair< unsigned int, unsigned int > HGCalHistoClusteringImplSA::mean_coordinate(unsigned int Sum_Wc, unsigned int Sum_W) const {
  double intpart;
  double frac =  modf((double)Sum_Wc/Sum_W,&intpart)*pow(2,2);
  return { (unsigned int)intpart, (unsigned int)frac };
}

std::pair< unsigned int, unsigned int > HGCalHistoClusteringImplSA::sigma_Coordinate(unsigned int Sum_W, unsigned long int Sum_Wc2, unsigned int Sum_Wc) const {
  unsigned long int N = Sum_W*Sum_Wc2 - pow(Sum_Wc,2);
  unsigned long int D = pow(Sum_W,2);
  double intpart;
  double frac =  modf((double)sqrt(N/D),&intpart)*pow(2,1);
  return { (unsigned int)intpart, (unsigned int)frac };
}

std::pair< unsigned int, unsigned int > HGCalHistoClusteringImplSA::energy_ratio(unsigned int E_N, unsigned int E_D) const {
  if ( E_D == 0 ) {
    return { 0 , 0 };
  } else {
    double intpart;
    double frac =  modf((double)E_N/E_D,&intpart)*pow(2,8);
    return { (unsigned int)intpart, (unsigned int)frac };
  }
}

std::vector<int> HGCalHistoClusteringImplSA::showerLengthProperties(unsigned long int layerBits) const {
  int counter = 0;
  int firstLayer = 0;
  bool firstLayerFound = false;
  int lastLayer = 0;
  std::vector<int> layerBits_array;

  for (int idx = 0; idx<36; idx++) {
    if ( (layerBits&(1L<<(35-idx)) ) >= 1L ) {
      if ( !firstLayerFound ) {
        firstLayer = idx+1;
        firstLayerFound=true;
      }
      lastLayer = idx+1;
      counter += 1;
    } else {
      layerBits_array.push_back(counter);
      counter = 0;
    }
  }
  int showerLen = lastLayer - firstLayer + 1;
  int coreShowerLen = 36;
  if ( layerBits_array.size()>0 ) {
    coreShowerLen = *std::max_element(layerBits_array.begin(), layerBits_array.end());
  }

  std::vector<int> output = {firstLayer,lastLayer, showerLen, coreShowerLen};

  return output;
}

void HGCalHistoClusteringImplSA::clusterProperties(  HGCalClusterSAPtrCollection& clusterSums ) const {
   unsigned int nTCs = 0;
  //  std::cout << "Cluster Prop" << std::endl;
   for ( const auto& c : clusterSums ) {
     if ( c->n_tc_w() == 0 ) continue;
      std::pair< unsigned int, unsigned int > sigmaEnergy = sigma_Energy( c->n_tc_w(), c->w2(), c->w() );
      c->set_Sigma_E_Quotient( sigmaEnergy.first );
      c->set_Sigma_E_Fraction( sigmaEnergy.second );
      std::pair< unsigned int, unsigned int > Mean_z = mean_coordinate( c->wz(), c->w() );
      c->set_Mean_z_Quotient( Mean_z.first );
      c->set_Mean_z_Fraction( Mean_z.second );
      std::pair< unsigned int, unsigned int > Mean_phi = mean_coordinate( c->wphi(), c->w() );      
      c->set_Mean_phi_Quotient( Mean_phi.first );
      c->set_Mean_phi_Fraction( Mean_phi.second );
      std::pair< unsigned int, unsigned int > Mean_eta = mean_coordinate( c->weta(), c->w() );
      c->set_Mean_eta_Quotient( Mean_eta.first );
      c->set_Mean_eta_Fraction( Mean_eta.second );
      std::pair< unsigned int, unsigned int > Mean_roz = mean_coordinate( c->wroz(), c->w() );
      c->set_Mean_roz_Quotient( Mean_roz.first );
      c->set_Mean_roz_Fraction( Mean_roz.second );
      std::pair< unsigned int, unsigned int > Sigma_z = sigma_Coordinate( c->w(), c->wz2(), c->wz() );
      c->set_Sigma_z_Quotient( Sigma_z.first );
      c->set_Sigma_z_Fraction( Sigma_z.second );
      std::pair< unsigned int, unsigned int > Sigma_phi = sigma_Coordinate( c->w(), c->wphi2(), c->wphi() );
      c->set_Sigma_phi_Quotient( Sigma_phi.first );
      c->set_Sigma_phi_Fraction( Sigma_phi.second );
      std::pair< unsigned int, unsigned int > Sigma_eta = sigma_Coordinate( c->w(), c->weta2(), c->weta() );
      c->set_Sigma_eta_Quotient( Sigma_eta.first );
      c->set_Sigma_eta_Fraction( Sigma_eta.second );
      std::pair< unsigned int, unsigned int > Sigma_roz = sigma_Coordinate( c->w(), c->wroz2(), c->wroz() );
      c->set_Sigma_roz_Quotient( Sigma_roz.first );
      c->set_Sigma_roz_Fraction( Sigma_roz.second );
      std::vector<int> layeroutput = showerLengthProperties( c->layerbits() ); 
      c->set_FirstLayer( layeroutput[0] );
      c->set_LastLayer( layeroutput[1] );
      c->set_ShowerLen( layeroutput[2] );
      c->set_CoreShowerLen( layeroutput[3] );
      std::pair< unsigned int, unsigned int > E_EM_over_E = energy_ratio( c->e_em() , c->e() ); 
      c->set_E_EM_over_E_Quotient( E_EM_over_E.first );
      c->set_E_EM_over_E_Fraction( E_EM_over_E.second );
      std::pair< unsigned int, unsigned int > E_EM_core_over_E_EM = energy_ratio( c->e_em_core() , c->e() );
      c->set_E_EM_core_over_E_EM_Quotient( E_EM_core_over_E_EM.first );
      c->set_E_EM_core_over_E_EM_Fraction( E_EM_core_over_E_EM.second );
      std::pair< unsigned int, unsigned int > E_H_early_over_E =  energy_ratio( c->e_h_early() , c->e() );
      c->set_E_H_early_over_E_Quotient( E_H_early_over_E.first );
      c->set_E_H_early_over_E_Fraction( E_H_early_over_E.second );


      // std::cout << c->clock() << " " << c->index() << " " << c->n_tc() << " " << c->e() << " " << c->e_em() << " " << c->e_em_core() << " " << c->e_h_early() << " N "<< c->Sigma_E_Quotient()<< " " <<c->Sigma_E_Fraction() << " " << c->Mean_z_Quotient() << " " << c->Mean_z_Fraction() << " " << c->Mean_phi_Quotient() << " "<< c->Mean_phi_Fraction() << " " << c->Mean_eta_Quotient() << " " << c->Mean_eta_Fraction() << " " << c->Mean_roz_Quotient() << " " << c->Mean_roz_Fraction() <<  " "<< c->Sigma_z_Quotient() << " "<< c->Sigma_z_Fraction() << " "<< c->Sigma_phi_Quotient() << " " << c-> Sigma_phi_Fraction() << " " << c->Sigma_eta_Quotient() << " " << c->Sigma_eta_Fraction() << " " << c->Sigma_roz_Quotient() << " "<< c->Sigma_roz_Fraction() << " "<< c->FirstLayer() <<" "<< c->LastLayer() << " "<< c->ShowerLen() << " " << c->CoreShowerLen() << " "<< c->E_EM_over_E_Quotient() << " " << c->E_EM_over_E_Fraction() << " " << c->E_EM_core_over_E_EM_Quotient() << " "<< c->E_EM_core_over_E_EM_Fraction() << " " << c->E_H_early_over_E_Quotient() << " " << c->E_H_early_over_E_Fraction() << std::endl;
      nTCs += c->n_tc();
   }
  //  std::cout << nTCs << std::endl;
}

void HGCalHistoClusteringImplSA::initializeTriggerCellDistGrid( HGCalTriggerCellSAPtrCollections& grid, unsigned int nX, unsigned int nY ) const {
  for (unsigned int iX = 0; iX < nX; ++iX ) {
    HGCalTriggerCellSAPtrCollection temp;
    for (unsigned int iY = 0; iY < nY; ++iY ) {
      temp.emplace_back(make_shared<HGCalTriggerCell>());
    }
    grid.push_back(temp);
  }
}

void HGCalHistoClusteringImplSA::runDistServers( const HGCalTriggerCellSAPtrCollections& gridIn,
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

void HGCalHistoClusteringImplSA::printHistogram( HGCalHistogramCellSAPtrCollection& histogram ) const {
  for ( unsigned int iRow = 0; iRow < config_.cRows();  ++iRow ) {
    for ( unsigned int iCol = 0; iCol < config_.cColumns();  ++iCol ) {
      unsigned binIndex = config_.cColumns() * iRow + iCol;
      std::cout << histogram.at(binIndex)->S() << " ";
    }
    std::cout << std::endl;
  }
}


