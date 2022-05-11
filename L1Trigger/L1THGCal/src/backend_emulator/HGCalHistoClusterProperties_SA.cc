#include "L1Trigger/L1THGCal/interface/backend_emulator/HGCalHistoClusterProperties_SA.h"

#include <cmath>
#include <algorithm>

using namespace std;
using namespace l1thgcfirmware;

HGCalHistoClusterProperties::HGCalHistoClusterProperties(ClusterAlgoConfig& config) : config_(config) {}

void HGCalHistoClusterProperties::runClusterProperties(const l1thgcfirmware::HGCalClusterSAPtrCollection& protoClustersIn, const CentroidHelperPtrCollection& readoutFlags, HGCalClusterSAPtrCollection& clustersOut ) const {
  // Cluster properties
  HGCalClusterSAPtrCollection clusterAccumulation;
  clusterSum( protoClustersIn, readoutFlags, clusterAccumulation, clustersOut );
  clusterProperties(clustersOut);
}

void HGCalHistoClusterProperties::clusterSum( const HGCalClusterSAPtrCollection& protoClusters, const CentroidHelperPtrCollection& readoutFlags, HGCalClusterSAPtrCollection& clusterAccumulation, HGCalClusterSAPtrCollection& clusterSums ) const {

  HGCalClusterSAPtrCollections protoClustersPerColumn( config_.cColumns() );
  vector<unsigned int> clock( config_.cColumns(), 0 );
  for ( const auto& protoCluster : protoClusters ) {
    auto index = protoCluster->index();
    // Do we need to make a copy of protoCluster here?
    protoClustersPerColumn.at( index ).push_back( make_unique<HGCalCluster>( *protoCluster ) );
  }

  map<unsigned int, HGCalClusterSAPtr> sums;

  for ( const auto& flag : readoutFlags ) {
    auto accumulator = make_unique<HGCalCluster>( 0,
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

    if ( sums.find( flag->clock() ) == sums.end() ) {
      auto sum = make_unique<HGCalCluster>( flag->clock() + 7, // Magic numbers
                                            0,
                                            true, true
                                          );
      sums[flag->clock()] = move(sum);
    }

    *(sums.at( flag->clock() )) += *accumulator;

    clusterAccumulation.push_back( move(accumulator) );

  }

  for (auto& sum: sums) {
    clusterSums.push_back( move( sum.second ) );
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

void HGCalHistoClusterProperties::clusterProperties(  HGCalClusterSAPtrCollection& clusterSums ) const {
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

std::pair< unsigned int, unsigned int > HGCalHistoClusterProperties::sigma_Energy(unsigned int N_TC_W, unsigned long int Sum_W2, unsigned int Sum_W) const {
  unsigned long int N = N_TC_W*Sum_W2 - pow(Sum_W,2);
  unsigned long int D = pow(N_TC_W,2);
  double intpart;
  double frac =  modf(sqrt(N/D),&intpart)*pow(2,1);
  return { (unsigned int)intpart, (unsigned int)frac };
}

std::pair< unsigned int, unsigned int > HGCalHistoClusterProperties::mean_coordinate(unsigned int Sum_Wc, unsigned int Sum_W) const {
  double intpart;
  double frac =  modf((double)Sum_Wc/Sum_W,&intpart)*pow(2,2);
  return { (unsigned int)intpart, (unsigned int)frac };
}

std::pair< unsigned int, unsigned int > HGCalHistoClusterProperties::sigma_Coordinate(unsigned int Sum_W, unsigned long int Sum_Wc2, unsigned int Sum_Wc) const {
  unsigned long int N = Sum_W*Sum_Wc2 - pow(Sum_Wc,2);
  unsigned long int D = pow(Sum_W,2);
  double intpart;
  double frac =  modf((double)sqrt(N/D),&intpart)*pow(2,1);
  return { (unsigned int)intpart, (unsigned int)frac };
}

std::pair< unsigned int, unsigned int > HGCalHistoClusterProperties::energy_ratio(unsigned int E_N, unsigned int E_D) const {
  if ( E_D == 0 ) {
    return { 0 , 0 };
  } else {
    double intpart;
    double frac =  modf((double)E_N/E_D,&intpart)*pow(2,8);
    return { (unsigned int)intpart, (unsigned int)frac };
  }
}

std::vector<int> HGCalHistoClusterProperties::showerLengthProperties(unsigned long int layerBits) const {
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

