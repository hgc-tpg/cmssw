#include "L1Trigger/L1THGCal/interface/backend_emulator/HGCalHistoSeeding_SA.h"

#include <cmath>

using namespace std;
using namespace l1thgcfirmware;

HGCalHistoSeeding::HGCalHistoSeeding(ClusterAlgoConfig& config) : config_(config) {}

void HGCalHistoSeeding::runSeeding(const HGCalTriggerCellSAPtrCollection& triggerCellsIn, HGCalHistogramCellSAPtrCollection& histogramOut ) const {

  HGCalHistogramCellSAPtrCollection histoCells;
  triggerCellToHistogramCell( triggerCellsIn, histoCells );
  makeHistogram( histoCells, histogramOut );

  // Smearing
  smearHistogram1D( histogramOut );
  normalizeArea( histogramOut );
  smearHistogram2D( histogramOut );

  //Maxima finding
  thresholdMaximaFinder( histogramOut );
  localMaximaFinder( histogramOut );
  calculateAveragePosition( histogramOut );
}

// Convert each trigger cell into a histogram cell
// Phi bin corresponds to the stream the TC is in, r/z bin is calculated from the TC r/z 
void HGCalHistoSeeding::triggerCellToHistogramCell( const HGCalTriggerCellSAPtrCollection& triggerCellsIn, HGCalHistogramCellSAPtrCollection& histogramOut ) const {

  const unsigned int latency = config_.getStepLatency( TcToHc );

  histogramOut.clear();
  for ( auto& tc : triggerCellsIn ) {
    auto hc = make_unique<HGCalHistogramCell>( tc->clock() + latency,
                                               tc->index(),
                                               tc->energy(),
                                               tc->phi(),
                                               tc->rOverZ(),
                                               1,
                                               int( ( tc->rOverZ() - config_.rOverZHistOffset() )/ config_.rOverZBinSize() )
                                              );
    tc->setClock( hc->clock() );
    tc->setSortKey( hc->sortKey() );
    histogramOut.push_back( std::move(hc) );
  }

}

// Sum all input istogram cells (each representing one TC) to produce one histogram
void HGCalHistoSeeding::makeHistogram( const HGCalHistogramCellSAPtrCollection& histogramCells, HGCalHistogramCellSAPtrCollection& histogramOut ) const {

  const unsigned int latency = config_.getLatencyUpToAndIncluding( Hist );

  histogramOut.clear();
  for ( unsigned int iRow = 0; iRow < config_.cRows(); ++iRow ) {
    for ( unsigned int iColumn = 0; iColumn < config_.cColumns(); ++iColumn ) {
      histogramOut.push_back( make_unique<HGCalHistogramCell>( latency, iColumn, iRow ) );
    }
  }

  for ( const auto& hc : histogramCells ) {
    const unsigned int binIndex = config_.cColumns() * hc->sortKey() + hc->index() ;
    *histogramOut.at( binIndex ) += *hc;
  }
}

// Smearing at constant r/z using exponentially falling 2^-n kernel
// Maximum extent of smearing (how many neighbouring bins included) can vary for each r/z bin, and is specified in the configuration
void HGCalHistoSeeding::smearHistogram1D( HGCalHistogramCellSAPtrCollection& histogram ) const {

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
      unique_ptr<HGCalHistogramCell> l1 = make_unique<HGCalHistogramCell>(HGCalHistogramCell());
      unique_ptr<HGCalHistogramCell> l2 = make_unique<HGCalHistogramCell>(HGCalHistogramCell());
      unique_ptr<HGCalHistogramCell> r1 = make_unique<HGCalHistogramCell>(HGCalHistogramCell());
      unique_ptr<HGCalHistogramCell> r2 = make_unique<HGCalHistogramCell>(HGCalHistogramCell());


      if ( width >= 2 ) {
        if ( int(col - offset - 1)  >= 0 ) {
          l2 = make_unique<HGCalHistogramCell>(lHistogram[binIndex - offset - 1]/4);
        }
        if ( int(col + offset + 1) <= int(config_.cColumns()-1) ) {
          r2 = make_unique<HGCalHistogramCell>(lHistogram[binIndex + offset + 1]/4);
        }
      }
      
      if ( int(col - offset)  >= 0 ) {
        l1 = make_unique<HGCalHistogramCell>(lHistogram[binIndex - offset]/2);
      }
      if ( int(col + offset)  <= int(config_.cColumns()-1) ) {
        r1 = make_unique<HGCalHistogramCell>(lHistogram[binIndex + offset]/2);
      }
      *hc += ( ( *l2 + *l1 ) / scale + ( *r2 + *r1 ) / scale );
      scale *= 4;
      width -= 2;
      offset += 2;
    }
  }
}

void HGCalHistoSeeding::normalizeArea( HGCalHistogramCellSAPtrCollection& histogram ) const {
  const unsigned int stepLatency = config_.getStepLatency( NormArea );
  for ( unsigned int iBin = 0; iBin < histogram.size(); ++iBin ) {
    HGCalHistogramCell& hc = *histogram.at(iBin);
    hc.addLatency( stepLatency );
    hc.setS( int( 1.0*hc.S() * config_.areaNormalization(hc.sortKey()) / pow(2,config_.nBitsAreaNormLUT()+1) ) );
  }
}

// Smearing at constant phi with exponentially falling 2^-n kernel
// Limited to +/- 1 bin in phi (only consider nearest neighbour), so adding half of the row above and below
void HGCalHistoSeeding::smearHistogram2D( HGCalHistogramCellSAPtrCollection& histogram ) const {
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
      *hc += (lHistogram[binIndex - config_.cColumns()] / 2 );
    }
    if ( row + 1 <= int(config_.cRows()-1) ) {
      *hc += (lHistogram[binIndex + config_.cColumns()] / 2 );
    }
  }
}

// Simple threshold maxima finder
// Keep only bins passing threshold
// Threshold can vary vs r/z
void HGCalHistoSeeding::thresholdMaximaFinder( HGCalHistogramCellSAPtrCollection& histogram ) const {
  const unsigned int stepLatency = config_.getStepLatency( Maxima2D );
  for ( auto& hc : histogram ) {
    hc->addLatency( stepLatency );
    if ( hc->S() <= config_.thresholdMaxima( hc->sortKey() ) ) {
      hc->setS(0);
      hc->setX(0);
      hc->setY(0);
      hc->setN(0);
    }
  }
}

// Temporary simulation of local maxima finder
// Not an emulation of any firmware
void HGCalHistoSeeding::localMaximaFinder( HGCalHistogramCellSAPtrCollection& histogram ) const {
  // const std::vector<unsigned> maximaWidths{ 6, 5, 5, 5, 4, 4, 4, 4, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 }; // Padded this with 4*1 at the end.  Firmware this was based on has 40 entries, but there are 44 bins.
  const std::vector<unsigned> maximaWidths(config_.cRows(), 1);

  for ( auto& hc : histogram ) {
    if ( hc->S() > 0 ) {
      const int colRef = hc->index();
      const int rowRef = hc->sortKey();
      bool isMaxima = true;
      const int phiRange = maximaWidths.at(hc->sortKey());
      for ( int colOffset = -1 * phiRange; colOffset <= phiRange; ++colOffset ) {
        const int col = colRef + colOffset;
        if ( col < 0 || col >= (int) config_.cColumns() ) continue;
        for ( int rowOffset = -1; rowOffset <= 1; ++ rowOffset ) {
          const int row = rowRef + rowOffset;
          if ( row < 0 || row >= (int) config_.cRows() ) continue;
          const unsigned int binIndex = config_.cColumns() * row + col;
          const auto& bin = histogram.at(binIndex);

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
}


// Calculate average of phi and r/z (now labelled X and Y) of TCs within a bin
// As performed in firmware, where the sum of TC phi or r/z is multiplied by 2^17/N_TCs (factors stored in a LUT)
// Then result is shifted by 17 bits 
void HGCalHistoSeeding::calculateAveragePosition( HGCalHistogramCellSAPtrCollection& histogram ) const {
  const unsigned int stepLatency = config_.getStepLatency( CalcAverage );
  for ( auto& hc : histogram ) {
    hc->addLatency( stepLatency );
    if ( hc->N() > 0 ) {
      unsigned int inv_N = int( round(1.0 * 0x1FFFF / hc->N() ) );
      hc->setX( ( hc->X() * inv_N ) >> 17 );
      hc->setY( ( hc->Y() * inv_N ) >> 17 );
    }
  }
}

// Useful debugging function for printing histogram contents
void HGCalHistoSeeding::printHistogram( const HGCalHistogramCellSAPtrCollection& histogram ) const {
  for ( unsigned int iRow = 0; iRow < config_.cRows();  ++iRow ) {
    for ( unsigned int iCol = 0; iCol < config_.cColumns();  ++iCol ) {
      unsigned binIndex = config_.cColumns() * iRow + iCol;
      std::cout << histogram.at(binIndex)->S() << " ";
    }
    std::cout << std::endl;
  }
}