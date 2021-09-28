#ifndef L1Trigger_L1THGCal_HGCalCluster_SA_h
#define L1Trigger_L1THGCal_HGCalCluster_SA_h

#include <vector>
#include <memory>
#include <iostream>
namespace l1thgcfirmware {

  class HGCalCluster {
  public:
    HGCalCluster( unsigned int clock,
                  unsigned int index,
                  bool frameValid,
                  bool dataValid)
                : clock_( clock ),
                  index_(index),
                  layerBits_(0),
                  frameValid_(frameValid),
                  dataValid_(dataValid),
                  n_tc_(0),
                  e_(0),
                  e_em_(0),
                  e_em_core_(0),
                  e_h_early_(0),
                  w_(0),
                  n_tc_w_(0),
                  w2_(0),
                  wz_(0),
                  weta_(0),
                  wphi_(0),
                  wroz_(0),
                  wz2_(0),
                  weta2_(0),
                  wphi2_(0),
                  wroz2_(0),
                  layerbits_(0),
                  sat_tc_(false),
                  shapeq_(1),
                  Sigma_E_Quotient_(0),
                  Sigma_E_Fraction_(0),
                  Mean_z_Quotient_(0),
                  Mean_z_Fraction_(0),
                  Mean_phi_Quotient_(0),
                  Mean_phi_Fraction_(0),
                  Mean_eta_Quotient_(0),
                  Mean_eta_Fraction_(0),
                  Mean_roz_Quotient_(0),
                  Mean_roz_Fraction_(0),
                  Sigma_z_Quotient_(0),
                  Sigma_z_Fraction_(0),
                  Sigma_phi_Quotient_(0),
                  Sigma_phi_Fraction_(0),
                  Sigma_eta_Quotient_(0),
                  Sigma_eta_Fraction_(0),
                  Sigma_roz_Quotient_(0),
                  Sigma_roz_Fraction_(0),
                  FirstLayer_(0),
                  LastLayer_(0),
                  ShowerLen_(0),
                  CoreShowerLen_(0),
                  E_EM_over_E_Quotient_(0),
                  E_EM_over_E_Fraction_(0),
                  E_EM_core_over_E_EM_Quotient_(0),
                  E_EM_core_over_E_EM_Fraction_(0),
                  E_H_early_over_E_Quotient_(0),
                  E_H_early_over_E_Fraction_(0)
                  {}

    ~HGCalCluster(){}

    std::pair< unsigned int, unsigned int > Sigma_Energy(unsigned int N_TC_W, unsigned long int Sum_W2, unsigned int Sum_W);
    std::pair< unsigned int, unsigned int > Mean_coordinate(unsigned int Sum_Wc, unsigned int Sum_W);
    std::pair< unsigned int, unsigned int > Sigma_Coordinate(unsigned int Sum_W, unsigned long int Sum_Wc2, unsigned int Sum_Wc);
    std::pair< unsigned int, unsigned int > Energy_ratio(unsigned int E_N, unsigned int E_D);
    std::vector<int> ShowerLengthProperties(unsigned long int layerBits);

    // Setters
    void setClock( const unsigned int clock ) { clock_ = clock; }
    void setIndex( const unsigned int index ) { index_ = index; }
    void setDataValid( const bool dataValid ) { dataValid_ = dataValid; }

    void set_n_tc( unsigned int n_tc ) { n_tc_ = n_tc; }
    void set_e( unsigned int e ) { e_ = e; }
    void set_e_em( unsigned int e_em ) { e_em_ = e_em; }
    void set_e_em_core( unsigned int e_em_core ) { e_em_core_ = e_em_core; }
    void set_e_h_early( unsigned int e_h_early ) { e_h_early_ = e_h_early; }
    void set_w( unsigned int w ) { w_ = w; }
    void set_n_tc_w( unsigned int n_tc_w ) { n_tc_w_ = n_tc_w; }
    void set_w2( unsigned long int w2 ) { w2_ = w2; }
    void set_wz( unsigned int wz ) { wz_ = wz; }
    void set_weta( unsigned int weta ) { weta_ = weta; }
    void set_wphi( unsigned int wphi ) { wphi_ = wphi; }
    void set_wroz( unsigned int wroz ) { wroz_ = wroz; }
    void set_wz2( unsigned long int wz2 ) { wz2_ = wz2; }
    void set_weta2( unsigned long int weta2 ) { weta2_ = weta2; }
    void set_wphi2( unsigned long int wphi2 ) { wphi2_ = wphi2; }
    void set_wroz2( unsigned long int wroz2 ) { wroz2_ = wroz2; }
    void set_layerbits( unsigned long int layerbits ) { layerbits_ = layerbits; }
    void set_sat_tc( bool sat_tc ) { sat_tc_ = sat_tc; }
    void set_shapeq( unsigned int shapeq ) { shapeq_ = shapeq; }

    void set_Sigma_E_Quotient( unsigned long int Sigma_E_Quotient ) { Sigma_E_Quotient_ = Sigma_E_Quotient; } 
    void set_Sigma_E_Fraction( unsigned long int Sigma_E_Fraction ) { Sigma_E_Fraction_ = Sigma_E_Fraction; }
    void set_Mean_z_Quotient( unsigned long int Mean_z_Quotient ) { Mean_z_Quotient_ = Mean_z_Quotient; }
    void set_Mean_z_Fraction( unsigned long int Mean_z_Fraction ) { Mean_z_Fraction_ = Mean_z_Fraction; }
    void set_Mean_phi_Quotient( unsigned long int Mean_phi_Quotient ) { Mean_phi_Quotient_ = Mean_phi_Quotient; }
    void set_Mean_phi_Fraction( unsigned long int Mean_phi_Fraction ) { Mean_phi_Fraction_ = Mean_phi_Fraction; }
    void set_Mean_eta_Quotient( unsigned long int Mean_eta_Quotient ) { Mean_eta_Quotient_ = Mean_eta_Quotient; }
    void set_Mean_eta_Fraction( unsigned long int Mean_eta_Fraction ) { Mean_eta_Fraction_ = Mean_eta_Fraction; }
    void set_Mean_roz_Quotient( unsigned long int Mean_roz_Quotient ) { Mean_roz_Quotient_ = Mean_roz_Quotient; }
    void set_Mean_roz_Fraction( unsigned long int Mean_roz_Fraction ) { Mean_roz_Fraction_ = Mean_roz_Fraction; }
    void set_Sigma_z_Quotient( unsigned long int Sigma_z_Quotient ) { Sigma_z_Quotient_ = Sigma_z_Quotient; }
    void set_Sigma_z_Fraction( unsigned long int Sigma_z_Fraction ) { Sigma_z_Fraction_ = Sigma_z_Fraction; }
    void set_Sigma_phi_Quotient( unsigned long int Sigma_phi_Quotient ) { Sigma_phi_Quotient_ = Sigma_phi_Quotient; }
    void set_Sigma_phi_Fraction( unsigned long int Sigma_phi_Fraction ) { Sigma_phi_Fraction_ = Sigma_phi_Fraction; }
    void set_Sigma_eta_Quotient( unsigned long int Sigma_eta_Quotient ) { Sigma_eta_Quotient_ = Sigma_eta_Quotient; }
    void set_Sigma_eta_Fraction( unsigned long int Sigma_eta_Fraction ) { Sigma_eta_Fraction_ = Sigma_eta_Fraction; }
    void set_Sigma_roz_Quotient( unsigned long int Sigma_roz_Quotient ) { Sigma_roz_Quotient_ = Sigma_roz_Quotient; }
    void set_Sigma_roz_Fraction( unsigned long int Sigma_roz_Fraction ) { Sigma_roz_Fraction_ = Sigma_roz_Fraction; }
    void set_FirstLayer( unsigned long int FirstLayer ) { FirstLayer_ = FirstLayer; }
    void set_LastLayer( unsigned long int LastLayer ) { LastLayer_ = LastLayer; }
    void set_ShowerLen( unsigned long int ShowerLen ) { ShowerLen_ = ShowerLen; }
    void set_CoreShowerLen( unsigned long int CoreShowerLen ) { CoreShowerLen_ = CoreShowerLen; }
    void set_E_EM_over_E_Quotient( unsigned long int E_EM_over_E_Quotient ) { E_EM_over_E_Quotient_ = E_EM_over_E_Quotient; }
    void set_E_EM_over_E_Fraction( unsigned long int E_EM_over_E_Fraction ) { E_EM_over_E_Fraction_ = E_EM_over_E_Fraction; }
    void set_E_EM_core_over_E_EM_Quotient( unsigned long int E_EM_core_over_E_EM_Quotient ) { E_EM_core_over_E_EM_Quotient_ = E_EM_core_over_E_EM_Quotient; }
    void set_E_EM_core_over_E_EM_Fraction( unsigned long int E_EM_core_over_E_EM_Fraction ) { E_EM_core_over_E_EM_Fraction_ = E_EM_core_over_E_EM_Fraction; }
    void set_E_H_early_over_E_Quotient( unsigned long int E_H_early_over_E_Quotient ) { E_H_early_over_E_Quotient_ = E_H_early_over_E_Quotient; }
    void set_E_H_early_over_E_Fraction( unsigned long int E_H_early_over_E_Fraction ) { E_H_early_over_E_Fraction_ = E_H_early_over_E_Fraction; }



    // Getters
    unsigned int clock() const { return clock_; }
    unsigned int index() const { return index_; }
    bool frameValid() const { return frameValid_; }
    bool dataValid() const { return dataValid_; }

    unsigned int n_tc() const { return n_tc_; }
    unsigned int e() const { return e_; }
    unsigned int e_em() const { return e_em_; }
    unsigned int e_em_core() const { return e_em_core_; }
    unsigned int e_h_early() const { return e_h_early_; }
    unsigned int w() const { return w_; }
    unsigned int n_tc_w() const { return n_tc_w_; }
    unsigned long int w2() const { return w2_; }
    unsigned int wz() const { return wz_; }
    unsigned int weta() const { return weta_; }
    unsigned int wphi() const { return wphi_; }
    unsigned int wroz() const { return wroz_; }
    unsigned long int wz2() const { return wz2_; }
    unsigned long int weta2() const { return weta2_; }
    unsigned long int wphi2() const { return wphi2_; }
    unsigned long int wroz2() const { return wroz2_; }
    unsigned long int layerbits() const { return layerbits_; }
    bool sat_tc() const { return sat_tc_; }
    unsigned int shapeq() const { return shapeq_; }

    unsigned long int Sigma_E_Quotient()  const { return Sigma_E_Quotient_; }
    unsigned long int Sigma_E_Fraction() const { return Sigma_E_Fraction_; }
    unsigned long int Mean_z_Quotient() const { return Mean_z_Quotient_; }
    unsigned long int Mean_z_Fraction() const { return Mean_z_Fraction_; }
    unsigned long int Mean_phi_Quotient() const { return Mean_phi_Quotient_; }
    unsigned long int Mean_phi_Fraction() const { return Mean_phi_Fraction_; }
    unsigned long int Mean_eta_Quotient() const { return Mean_eta_Quotient_; }
    unsigned long int Mean_eta_Fraction() const { return Mean_eta_Fraction_; }
    unsigned long int Mean_roz_Quotient() const { return Mean_roz_Quotient_; }
    unsigned long int Mean_roz_Fraction() const { return Mean_roz_Fraction_; }
    unsigned long int Sigma_z_Quotient() const { return Sigma_z_Quotient_; }
    unsigned long int Sigma_z_Fraction() const { return Sigma_z_Fraction_; }
    unsigned long int Sigma_phi_Quotient() const { return Sigma_phi_Quotient_; }
    unsigned long int Sigma_phi_Fraction() const { return Sigma_phi_Fraction_; }
    unsigned long int Sigma_eta_Quotient() const { return Sigma_eta_Quotient_; }
    unsigned long int Sigma_eta_Fraction() const { return Sigma_eta_Fraction_; }
    unsigned long int Sigma_roz_Quotient() const { return Sigma_roz_Quotient_; }
    unsigned long int Sigma_roz_Fraction() const { return Sigma_roz_Fraction_; }
    unsigned long int FirstLayer() const { return FirstLayer_; }
    unsigned long int LastLayer() const { return LastLayer_; }
    unsigned long int ShowerLen() const { return ShowerLen_; }
    unsigned long int CoreShowerLen() const { return CoreShowerLen_; }
    unsigned long int E_EM_over_E_Quotient() const { return E_EM_over_E_Quotient_; }
    unsigned long int E_EM_over_E_Fraction() const { return E_EM_over_E_Fraction_; }
    unsigned long int E_EM_core_over_E_EM_Quotient() const { return E_EM_core_over_E_EM_Quotient_; }
    unsigned long int E_EM_core_over_E_EM_Fraction() const { return E_EM_core_over_E_EM_Fraction_; }
    unsigned long int E_H_early_over_E_Quotient() const { return E_H_early_over_E_Quotient_; }
    unsigned long int E_H_early_over_E_Fraction() const { return E_H_early_over_E_Fraction_; }

    // Operators
    const HGCalCluster& operator+=(const HGCalCluster& hc);

  private:
    unsigned int clock_;
    unsigned int index_;
    unsigned int layerBits_;
    bool frameValid_;
    bool dataValid_;

    unsigned int n_tc_;
    unsigned int e_;
    unsigned int e_em_;
    unsigned int e_em_core_;
    unsigned int e_h_early_;
    unsigned int w_;
    unsigned int n_tc_w_;
    unsigned long int w2_;
    unsigned int wz_;
    unsigned int weta_;
    unsigned int wphi_;
    unsigned int wroz_;
    unsigned long int wz2_;
    unsigned long int weta2_;
    unsigned long int wphi2_;
    unsigned long int wroz2_;
    unsigned long int layerbits_;
    bool sat_tc_;
    unsigned int shapeq_;
    unsigned long int Sigma_E_Quotient_;
    unsigned long int Sigma_E_Fraction_;
    unsigned long int Mean_z_Quotient_;
    unsigned long int Mean_z_Fraction_;
    unsigned long int Mean_phi_Quotient_;
    unsigned long int Mean_phi_Fraction_;
    unsigned long int Mean_eta_Quotient_;
    unsigned long int Mean_eta_Fraction_;
    unsigned long int Mean_roz_Quotient_;
    unsigned long int Mean_roz_Fraction_;
    unsigned long int Sigma_z_Quotient_;
    unsigned long int Sigma_z_Fraction_;
    unsigned long int Sigma_phi_Quotient_;
    unsigned long int Sigma_phi_Fraction_;
    unsigned long int Sigma_eta_Quotient_;
    unsigned long int Sigma_eta_Fraction_;
    unsigned long int Sigma_roz_Quotient_;
    unsigned long int Sigma_roz_Fraction_;
    unsigned long int FirstLayer_;
    unsigned long int LastLayer_;
    unsigned long int ShowerLen_;
    unsigned long int CoreShowerLen_;
    unsigned long int E_EM_over_E_Quotient_;
    unsigned long int E_EM_over_E_Fraction_;
    unsigned long int E_EM_core_over_E_EM_Quotient_;
    unsigned long int E_EM_core_over_E_EM_Fraction_;
    unsigned long int E_H_early_over_E_Quotient_;
    unsigned long int E_H_early_over_E_Fraction_;

 
  };

  typedef std::vector<HGCalCluster> HGCalClusterSACollection;
  typedef std::shared_ptr<HGCalCluster> HGCalClusterSAPtr;
  typedef std::vector<HGCalClusterSAPtr > HGCalClusterSAPtrCollection;
  typedef std::vector< HGCalClusterSAPtrCollection > HGCalClusterSAPtrCollections;

}  // namespace l1thgcfirmware

#endif
