/**********************************************************************
 Created on : 14/04/2025
 Purpose    : C++ implementation of Stage 2 cluster properties calculation following Milos's thesis
              https://cds.cern.ch/record/2918629?ln=en
 Author     : Indranil Das, Research Associate, Imperial
 Email      : indranil.das@cern.ch | indra.ehep@gmail.com
**********************************************************************/
#ifndef TPGClusterProperties_h
#define TPGClusterProperties_h

#include "TPGLSBScales.hh"
#include "TPGBEDataformat.hh"
#ifdef CMSSW_GIT_HASH
#include "L1Trigger/L1THGCal/interface/backend_semiemulator/HGCalCluster_HW.h"
#else
#include "CMSSWCode/DataFormats/L1THGCal/interface/HGCalCluster_HW.h"
#endif

//Class to calculate the cluster properties
class TPGClusterProperties{
public:
  TPGClusterProperties() {}

  std::vector<int> showerLengthProperties(unsigned long int layerBits) const {
    //Collected from HGCalHistoClusterProperties::showerLengthProperties of L1Trigger/L1THGCal/src/backend_emulator/HGCalHistoClusterProperties_SA.cc
    // and updated for coreShowerLen as per PropertyCalculator.py of Milos

    int counter = 0;
    int firstLayer = 0;
    bool firstLayerFound = false;
    int lastLayer = 0;
    std::vector<int> layerBits_array;

    std::bitset<34> layerBitsBitset(layerBits);
    for (size_t i = 0; i < layerBitsBitset.size(); ++i) {
      bool bit = layerBitsBitset[34-1-i];
      if ( bit ) {
	if ( !firstLayerFound ) {
          firstLayer = i + 1;
          firstLayerFound = true;
	}
	lastLayer = i+1;
	counter += 1;
      } else {
	layerBits_array.push_back(counter);
	counter = 0;
      }
    }
    layerBits_array.push_back(counter);

    int showerLen = lastLayer - firstLayer + 1;
    //int coreShowerLen = config_.nTriggerLayers();
    int coreShowerLen = 0;
    if (!layerBits_array.empty()) {
      coreShowerLen = *std::max_element(layerBits_array.begin(), layerBits_array.end());
    }
    return {firstLayer, lastLayer, showerLen, coreShowerLen};
  }

  uint32_t convertRozToEta(uint32_t wroz, uint32_t w) {
    //Collected from HGCalHistoClusterProperties::convertRozToEta of L1Trigger/L1THGCal/src/backend_emulator/HGCalHistoClusterProperties_SA.cc
    //               and PropertyCalculator.py of Milos

    if(w==0) return 0;
    //========================== Constants =======================
    float mu_roz_lower_single = lsbScales->roz_min_L1T() / lsbScales->LSB_roz();
    float mu_roz_upper_single = lsbScales->roz_max_L1T() / lsbScales->LSB_roz();
    //=================================================================

    float muroz = float(wroz)/float(w);
    muroz = (muroz<mu_roz_lower_single)?mu_roz_lower_single:muroz;
    muroz = (muroz>mu_roz_upper_single)?mu_roz_upper_single:muroz;
    float mu_roz_local_single = muroz - mu_roz_lower_single;
    ap_ufixed<32,20, AP_RND, AP_SAT> mu_roz_local_scaled = mu_roz_local_single * lsbScales->c_roz_scaler_0();
    muroz = mu_roz_local_scaled.to_float();
    uint32_t roz = uint32_t(std::round(muroz));
    if ( roz > 1023 ) roz = 1023;
    uint32_t eta = clusPropLUT->getMuEta(uint32_t(roz));
    return eta;


  }

  ap_int<9> calc_phi(uint32_t wphi, uint32_t w, bool& saturatedPhi, bool& nominalPhi){

    if(w==0) return 0;

    float muphi = float(wphi)/float(w);
    ap_ufixed<32,20, AP_RND, AP_SAT> mu_phi_scaled = muphi * lsbScales->c_phi_scaler();
    muphi = mu_phi_scaled.to_float();
    int mu_phi_int = std::round(muphi);
    int wphi_et_tmp = mu_phi_int - lsbScales->c_Phi_Offset() ;
    saturatedPhi = ((wphi_et_tmp < lsbScales->maxPhiM()) or (wphi_et_tmp > lsbScales->maxPhiP()))?true:false ;
    wphi_et_tmp = (wphi_et_tmp < lsbScales->maxPhiM())?lsbScales->maxPhiM():wphi_et_tmp ;
    wphi_et_tmp = (wphi_et_tmp > lsbScales->maxPhiP())?lsbScales->maxPhiP():wphi_et_tmp ;
    nominalPhi = ((wphi_et_tmp > lsbScales->nomPhiL()) and (wphi_et_tmp < lsbScales->nomPhiH()))?true:false ;
    ap_int<9> wphi_et = wphi_et_tmp ;

    // std::cout << "TPGClusterProperties::calc_phi muphi: " << muphi
    // 	      << ", c_phi_scaler: " << lsbScales->c_phi_scaler()
    // 	      << ", c_Phi_Offset: " << lsbScales->c_Phi_Offset()
    // 	      << ", mu_phi_scaled: " << mu_phi_scaled
    // 	      << ", mu_phi_int: " << mu_phi_int
    // 	      << ", wphi_et: " << wphi_et
    // 	      << std::endl;

    return wphi_et;
  }

  ap_uint<5> convertSigmaRozRozToSigmaEtaEta( uint64_t wroz2, uint32_t wroz, uint32_t w) {
    //Collected from HGCalHistoClusterProperties::convertSigmaRozRozToSigmaEtaEta of L1Trigger/L1THGCal/src/backend_emulator/HGCalHistoClusterProperties_SA.cc
    // Sigma eta eta calculation

    if(w==0) return 0;

    double muroz = float(wroz)/float(w);
    double mu_roz_ds_sat = (muroz>lsbScales->mu_roz_ds_upper())?lsbScales->mu_roz_ds_upper():muroz;
    mu_roz_ds_sat = (mu_roz_ds_sat<lsbScales->mu_roz_ds_lower())?lsbScales->mu_roz_ds_lower():mu_roz_ds_sat;
    double mu_roz_ds_local_single = mu_roz_ds_sat - lsbScales->mu_roz_ds_lower();
    double mu_roz_ds_local_scaled = mu_roz_ds_local_single * lsbScales->c_roz_scaler_1();
    ap_ufixed<32,20, AP_RND, AP_SAT> mu_roz_ds_local_fxp = mu_roz_ds_local_scaled;
    float mu_roz_ds_float = mu_roz_ds_local_fxp.to_float() ;
    int mu_roz_ds_int = std::round(mu_roz_ds_float);
    int mu_roz_addr = (mu_roz_ds_int>63)?63:mu_roz_ds_int;

    // # sigma_rozroz_single already calculated previously
    uint32_t sigma_rozroz_ds_int = sigma_coordinate(w, wroz2, wroz, lsbScales->c_sigma_roz_scaler_0());
    uint32_t sig_roz_addr = (sigma_rozroz_ds_int>63)?63:sigma_rozroz_ds_int;
    uint32_t sigma_eta_LUT_addr = mu_roz_addr * 64 + sig_roz_addr ; //# 64 because sig_roz is 6 bits wide
    uint32_t sigma_etaeta = clusPropLUT->getSigmaEtaEta(sigma_eta_LUT_addr);
    ap_uint<5> ret = sigma_etaeta;

    // std::cout << "TPGClusterProperties::convertSigmaRozRozToSigmaEtaEta : muroz: " << muroz
    // 		<<", mu_roz_ds_sat: " << mu_roz_ds_sat
    // 		<<", mu_roz_ds_local_scaled: " << mu_roz_ds_local_scaled
    // 		<<", mu_roz_ds_local_fxp: " << mu_roz_ds_local_fxp
    // 		<<", mu_roz_addr: " << mu_roz_addr
    // 		<<std::endl;
    // std::cout << "TPGClusterProperties::convertSigmaRozRozToSigmaEtaEta sigma_rozroz_ds_int: " << sigma_rozroz_ds_int
    // 		<<", sig_roz_addr: " << sig_roz_addr
    // 		<<", sigma_eta_LUT_addr: " << sigma_eta_LUT_addr
    // 		<<", sigma_etaeta: " << sigma_etaeta
    // 		<<", ret: " << ret
    // 		<<std::endl;

    return ret;
  }

  uint32_t sigma_coordinate(uint32_t w, uint64_t wc2, uint32_t wc, double scale ) const {
    //Collected from HGCalHistoClusterProperties::sigma_coordinate of L1Trigger/L1THGCal/src/backend_emulator/HGCalHistoClusterProperties_SA.cc

    if ( w == 0 ) return 0;
    float wt_obs = (float(w)*float(wc2) - float(wc) * float(wc))/(float(w) * float(w)) ;
    if (wt_obs<0.) return 0;
    float wt_obs_sqr = sqrt( wt_obs );
    float wt_obs_scaled = wt_obs_sqr * scale;
    ap_ufixed<32,20, AP_RND, AP_SAT> wt_obs_scaled_fxp = wt_obs_scaled;
    float wt_obs_fp1 = wt_obs_scaled_fxp.to_float() ;
    uint32_t sigma = std::round( wt_obs_fp1 );

    // std::cout << "TPGClusterProperties::sigma_coordinate  w: " << w << ", wc2: " << wc2 << ", wc: " << wc << ", scale: " << scale << std::endl;
    // std::cout << "TPGClusterProperties::sigma_coordinate wt_obs: " << wt_obs << ", wt_obs_sqr : " << wt_obs_sqr << ", wt_obs_scaled: " << wt_obs_scaled << ", wt_obs_fp1: " << wt_obs_fp1 << ", sigma : " << sigma  << std::endl;

    return sigma;

  }

  void ClusterProperties(const TPGBEDataformat::TcAccumulatorFW& accuInput, l1thgcfirmware::HGCalCluster_HW& l1TOutput, bool isPrint = false){
    //Collected from HGCalHistoClusterProperties::calcProperties of L1Trigger/L1THGCal/src/backend_emulator/HGCalHistoClusterProperties_SA.cc

    l1TOutput.clear();

    if(isPrint) std::cout<<"Calculating cluster properties" << std::endl;
    //// ================== First Word ===========================
    l1TOutput.e = accuInput.totE() * lsbScales->c_ET_scaler();
    l1TOutput.e_em = accuInput.ceeE() * lsbScales->c_ET_scaler();
    if(isPrint) std::cout<<"Set Energies" << std::endl;
    // std::cout << "Fractions : " << accuInput.totE() << " " << accuInput.totE() * lsbScales->c_ET_scaler() << " " << accuInput.ceeE() << " " << accuInput.ceeECore() << " " << accuInput.totE() << " " << lsbScales->c_frac_scaler() << " " << accuInput.ceeE()/accuInput.totE() << " " << lsbScales->c_frac_scaler() * accuInput.ceeE()/accuInput.totE() << std::endl;
    l1TOutput.fractionInCE_E = (accuInput.totE()==0)?0:(lsbScales->c_frac_scaler() * accuInput.ceeE()/accuInput.totE()) ;
    l1TOutput.fractionInCE_E = ( accuInput.ceeE() == accuInput.totE() ) ? l1thgcfirmware::eFraction_t(lsbScales->c_frac_scaler()-1) : l1TOutput.fractionInCE_E;
    l1TOutput.fractionInCoreCE_E = (accuInput.ceeE()==0)?0:(lsbScales->c_frac_scaler() * accuInput.ceeECore() / accuInput.ceeE()) ;
    l1TOutput.fractionInCoreCE_E = ( accuInput.ceeECore() == accuInput.ceeE() ) ? l1thgcfirmware::eFraction_t(lsbScales->c_frac_scaler()-1) : l1TOutput.fractionInCoreCE_E;
    l1TOutput.fractionInEarlyCE_E = (accuInput.totE()==0)?0:(lsbScales->c_frac_scaler() * accuInput.ceHEarly() /  accuInput.totE()) ;
    l1TOutput.fractionInEarlyCE_E = ( accuInput.ceHEarly() == accuInput.totE() ) ? l1thgcfirmware::eFraction_t(lsbScales->c_frac_scaler()-1) : l1TOutput.fractionInEarlyCE_E;
    // std::cout << "Fractions set: " << l1TOutput.fractionInCE_E << " " << l1TOutput.fractionInCoreCE_E << " " << l1TOutput.fractionInEarlyCE_E << std::endl;
    if(isPrint) std::cout<<"Set Fractions" << std::endl;
    l1TOutput.setGCTBits();
    if(isPrint) std::cout<<"Set GCT bits" << std::endl;
    std::vector<int> layeroutput = showerLengthProperties(accuInput.layerBits());
    if(isPrint) std::cout<<"Set Layer outputs" << std::endl;
    l1TOutput.firstLayer = layeroutput[0];
    l1TOutput.lastLayer = layeroutput[1];
    l1TOutput.showerLength = layeroutput[2];
    l1TOutput.coreShowerLength = layeroutput[3];
    l1TOutput.nTC = accuInput.numberOfTcs();
    //// ================== First Word ===========================
    if(isPrint) std::cout<<"Completed first word" << std::endl;

    //// ================== Second Word ===========================

    l1TOutput.w_eta = convertRozToEta( accuInput.sumWRoZ(), accuInput.sumW() );
    if (isPrint) std::cout<<"Completed w_eta" << std::endl;
    bool saturatedPhi = false, nominalPhi = false;
    l1TOutput.w_phi = calc_phi(accuInput.sumWPhi(), accuInput.sumW(), saturatedPhi, nominalPhi);
    if (isPrint) std::cout<<"Completed w_phi" << std::endl;
    float zratio = float(accuInput.sumWZ()) / float(accuInput.sumW()) ;
    ap_ufixed<32,20, AP_RND> wt_muz_fxp = zratio;
    uint32_t muz = std::round(wt_muz_fxp.to_float());
    l1TOutput.w_z = (muz>0xfff)?0:muz;
    if (isPrint) std::cout<<"Completed w_z" << std::endl;
    //l1TOutput.setQualityFlags(l1thgcfirmware::Scales::HGCaltoL1_et(accuInput.ceeECore()), l1thgcfirmware::Scales::HGCaltoL1_et(accuInput.ceHEarly()), accuInput.issatTC(), accuInput.shapeQ(), saturatedPhi, nominalPhi);
    l1TOutput.setQualityFlags(accuInput.ceeECore(), accuInput.ceHEarly(), accuInput.issatTC(), accuInput.shapeQ(), saturatedPhi, nominalPhi);
    //// ================== Second Word ===========================
    if(isPrint) std::cout<<"Completed second word" << std::endl;

    //// ================== Third Word ===========================
    uint32_t sigmaE = sigma_coordinate( accuInput.numberOfTcsW(), accuInput.sumW2(), accuInput.sumW(), lsbScales->c_sigma_E_scaler());
    l1TOutput.sigma_E = (sigmaE>0x7f)?0x7f:sigmaE;
    if(isPrint) std::cout<<"Completed sigma_E" << std::endl;

    uint32_t sigmaz = sigma_coordinate(accuInput.sumW(), accuInput.sumWZ2(), accuInput.sumWZ(), lsbScales->c_sigma_z_scaler());
    l1TOutput.sigma_z = (sigmaz>0x7f)?0x7f:sigmaz;
    if(isPrint) std::cout<<"Completed sigma_Z" << std::endl;

    uint32_t sigmaPhi = sigma_coordinate(accuInput.sumW(), accuInput.sumWPhi2(), accuInput.sumWPhi(), lsbScales->c_sigma_phi_scaler());
    l1TOutput.sigma_phi = (sigmaPhi>0x7f)?0x7f:sigmaPhi;
    if(isPrint) std::cout<<"Completed sigma_phi" << std::endl;

    unsigned int sigma_roz = sigma_coordinate(accuInput.sumW(), accuInput.sumWRoZ2(),  accuInput.sumWRoZ(), lsbScales->c_sigma_roz_scaler_1());
    l1TOutput.sigma_roz = (sigma_roz<127)?sigma_roz:127;
    if(isPrint) std::cout<<"Completed sigma_roz" << std::endl;

    uint32_t sigmaEta = convertSigmaRozRozToSigmaEtaEta( accuInput.sumWRoZ2(), accuInput.sumWRoZ(), accuInput.sumW());
    l1TOutput.sigma_eta = (sigmaEta>0x1f)?0x1f:sigmaEta;
    //// ================== Third Word ===========================
    if(isPrint) std::cout<<"Completed sigma_eta and third word" << std::endl;
  }

  void setClusPropLUT(const TPGStage2Configuration::ClusPropLUT *cplut) { clusPropLUT = cplut;}
  void setLSBScales(const TPGLSBScales::TPGStage2ClusterLSB *lsbSc) { lsbScales = lsbSc;}

private:
  const TPGStage2Configuration::ClusPropLUT *clusPropLUT;
  const TPGLSBScales::TPGStage2ClusterLSB *lsbScales;
};
#endif
