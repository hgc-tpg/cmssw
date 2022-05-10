#include "L1Trigger/L1THGCal/interface/backend_emulator/HGCalHistoClustering_SA.h"

#include <cmath>
#include <algorithm>

using namespace std;
using namespace l1thgcfirmware;

HGCalHistoClustering::HGCalHistoClustering(ClusterAlgoConfig& config) : config_(config) {}

void HGCalHistoClustering::runClustering(const HGCalTriggerCellSAPtrCollection& triggerCellsIn, const HGCalHistogramCellSAPtrCollection& histogramIn, HGCalTriggerCellSAPtrCollection& clusteredTriggerCellsOut, CentroidHelperPtrCollection& readoutFlagsOut  ) const {

  HGCalTriggerCellSAPtrCollection unclusteredTriggerCells;
  CentroidHelperPtrCollection prioritizedMaxima;
  clusterizer(triggerCellsIn, histogramIn, clusteredTriggerCellsOut, unclusteredTriggerCells, prioritizedMaxima, readoutFlagsOut );

}


void HGCalHistoClustering::clusterizer( const HGCalTriggerCellSAPtrCollection& triggerCellsIn, const HGCalHistogramCellSAPtrCollection& histogram, HGCalTriggerCellSAPtrCollection& clusteredTriggerCellsOut, HGCalTriggerCellSAPtrCollection& unclusteredTriggerCellsOut, CentroidHelperPtrCollection& prioritizedMaxima, CentroidHelperPtrCollection& readoutFlagsOut ) const {
  
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
      bool deltaMinus = (i>0) ? ( lastLatched[i]->column() - lastLatched[i-1]->column() ) > 6 : true ; // Magic numbers
      bool deltaPlus = (i<17) ? ( lastLatched[i+1]->column() - lastLatched[i]->column() ) > 6 : true ; // Magic numbers

      bool compareEMinus = (i>0) ? ( lastLatched[i]->energy() > lastLatched[i-1]->energy() ) : true; // Magic numbers
      bool compareEPlus = (i<17) ? ( lastLatched[i]->energy() >= lastLatched[i+1]->energy() ) : true; // Magic numbers

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