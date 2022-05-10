#include "L1Trigger/L1THGCal/interface/backend_emulator/HGCalHistoSeeding_SA.h"

#include <cmath>

using namespace std;
using namespace l1thgcfirmware;

HGCalHistoSeeding::HGCalHistoSeeding(ClusterAlgoConfig& config) : config_(config) {}

void HGCalHistoSeeding::runSeeding(const HGCalTriggerCellSAPtrCollection& triggerCellsIn, HGCalHistogramCellSAPtrCollection& histogramOut ) const {

  HGCalHistogramCellSAPtrCollection histoCells = triggerCellToHistogramCell( triggerCellsIn );
  histogramOut = makeHistogram( histoCells );

  // Smearing
  smearing1D( histogramOut );
  areaNormalization( histogramOut );
  smearing2D( histogramOut );

  //Maxima finding
  thresholdMaximaFinder( histogramOut );
  localMaximaFinder( histogramOut );
  calculateAveragePosition( histogramOut );
}

HGCalHistogramCellSAPtrCollection HGCalHistoSeeding::triggerCellToHistogramCell( const HGCalTriggerCellSAPtrCollection& triggerCellsIn ) const {

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

HGCalHistogramCellSAPtrCollection HGCalHistoSeeding::makeHistogram( HGCalHistogramCellSAPtrCollection histogramCells ) const {

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

void HGCalHistoSeeding::smearing1D( HGCalHistogramCellSAPtrCollection& histogram ) const {

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

void HGCalHistoSeeding::areaNormalization( HGCalHistogramCellSAPtrCollection& histogram ) const {
  const unsigned int stepLatency = config_.getStepLatency( NormArea );
  for ( unsigned int iBin = 0; iBin < histogram.size(); ++iBin ) {
    HGCalHistogramCell& hc = *histogram.at(iBin);
    hc.addLatency( stepLatency );
    hc *= config_.areaNormalization(hc.sortKey());
  }
}

void HGCalHistoSeeding::smearing2D( HGCalHistogramCellSAPtrCollection& histogram ) const {
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

void HGCalHistoSeeding::thresholdMaximaFinder( HGCalHistogramCellSAPtrCollection& histogram ) const {
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
void HGCalHistoSeeding::localMaximaFinder( HGCalHistogramCellSAPtrCollection& histogram ) const {
  // const std::vector<unsigned> maximaWidths{ 6, 5, 5, 5, 4, 4, 4, 4, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 }; // Padded this with 4*1 at the end.  Firmware this was based on has 40 entries, but there are 44 bins.
  const std::vector<unsigned> maximaWidths(config_.cRows(), 1);

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

void HGCalHistoSeeding::calculateAveragePosition( HGCalHistogramCellSAPtrCollection& histogram ) const {
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

void HGCalHistoSeeding::printHistogram( HGCalHistogramCellSAPtrCollection& histogram ) const {
  for ( unsigned int iRow = 0; iRow < config_.cRows();  ++iRow ) {
    for ( unsigned int iCol = 0; iCol < config_.cColumns();  ++iCol ) {
      unsigned binIndex = config_.cColumns() * iRow + iCol;
      std::cout << histogram.at(binIndex)->S() << " ";
    }
    std::cout << std::endl;
  }
}