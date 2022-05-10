#include "L1Trigger/L1THGCal/interface/backend_emulator/HGCalHistoClusterProperties_SA.h"

#include <cmath>
#include <algorithm>

using namespace std;
using namespace l1thgcfirmware;

HGCalHistoClusterProperties::HGCalHistoClusterProperties(ClusterAlgoConfig& config) : config_(config) {}

void HGCalHistoClusterProperties::runClusterProperties(const HGCalTriggerCellSAPtrCollection& triggerCellsIn, const CentroidHelperPtrCollection& readoutFlags, HGCalClusterSAPtrCollection& clustersOut ) const {
  // Cluster properties
  HGCalClusterSAPtrCollection protoClusters = triggerCellToCluster( triggerCellsIn );
  HGCalClusterSAPtrCollection clusterAccumulation;
  clusterSum( protoClusters, readoutFlags, clusterAccumulation, clustersOut );
  clusterProperties(clustersOut);
}

HGCalClusterSAPtrCollection HGCalHistoClusterProperties::triggerCellToCluster( const HGCalTriggerCellSAPtrCollection& clusteredTriggerCells ) const {

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

    cluster->add_constituent( tc );

    protoClusters.push_back( cluster );
  }

  // std::cout << "Output from triggerCellToCluster" << std::endl;
  // std::cout << "Protoclusters : " << protoClusters.size() << std::endl;
  // for ( const auto& pclus : protoClusters ) {

  //     std::cout << pclus->clock() << " " << pclus->index() << " " << pclus->n_tc() << " " << pclus->e() << " " << pclus->e_em() << " " << pclus->e_em_core() << " " << pclus->e_h_early() << " " << pclus->w() << " " << pclus->n_tc_w() << " " << pclus->weta2() << " " << pclus->wphi2() << std::endl;
  // }
  return protoClusters;
}

void HGCalHistoClusterProperties::clusterSum( const HGCalClusterSAPtrCollection& protoClusters, const CentroidHelperPtrCollection& readoutFlags, HGCalClusterSAPtrCollection& clusterAccumulation, HGCalClusterSAPtrCollection& clusterSums ) const {

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