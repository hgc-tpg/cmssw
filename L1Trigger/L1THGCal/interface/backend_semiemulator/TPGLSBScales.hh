/**********************************************************************
 Created on : 17/04/2025
 Purpose    : Classes for LSB, scaling values
 Author     : Indranil Das, Research Associate, Imperial
 Email      : indranil.das@cern.ch | indra.ehep@gmail.com
**********************************************************************/
#ifndef TPGLSBScales_h
#define TPGLSBScales_h

namespace TPGLSBScales{

  // class TPGFEHGCROCLSB{
  // public:
  //   TPGFEEHGCROCLSB() {}
  // private:
  // };

  // class TPGFEECONTLSB{
  // public:
  //   TPGFEECONTLSB() {}
  // private:
  // };

  // class TPGStage1TCLSB{
  // public:
  //   TPGStage1TCLSB() {}
  // private:
  // };

  // class TPGStage1TowerLSB{
  // public:
  //   TPGStage1TowerLSB() {}
  // private:
  // };

  //The source of following LSB, scaling and constants is https://gitlab.cern.ch/hgcal-tpg/Stage2/-/blob/master/ClusterProperties-Python/Emulator/MyConstants.py
  class TPGStage2ClusterLSB{
  public:
    TPGStage2ClusterLSB() {

      abs_z_global_min_    = 3210.5 ;  //in mm
      abs_z_global_max_    = 5140.3 ;  //in mm
      r_global_min_        = 313.6 ;   //in mm
      r_global_max_        = 2624.6 ;  //in mm
      roz_global_min_      = 364.3/4429 ;
      roz_global_max_      = 2624.6/4565.6 ; //0.574864202
      /*
	# Layer indices [not defined, use local function like TPGBEDataformat::TcAccumulator_FW::getTriggerLayer()]
	double triggerLayersGlobal[] = {1,3,5,7,9,11,13,15,17,19,21,23,25,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47};
	double triggerLayers[] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34};
      */
      double zLayer[]={ //in cm
		       322.155,     302, 325.212,     304, 328.269,     306, 331.326,     308, 334.383,     310,
		       337.440,     312, 340.497,     314, 343.554,     316, 346.611,     318, 349.993,     320,
		       353.375,     322, 356.757,     324, 360.139,     326, 367.976,     374, 380.586, 386.891,
		       393.196, 399.501, 405.806, 412.111, 418.416, 424.721, 431.026, 439.251, 447.476, 455.701,
		       463.926, 472.151, 480.376, 488.406, 496.826, 505.051, 513.276
       };
      for(int il=0;il<47;il++) _zLayer[il] = zLayer[il];

      maxTCET_	       = 512.;	//in GeV
      maxTCETbits_     = 18; //19 used in FW, 18 in current setup
      maxTCRoZbits_    = 12; //13 used in FW, 12 in current setup
      //maxTCRoZbits_    = 10;
      maxTCHistoColNo_ = 108;

      maxET_	     = 4097.5;	//in GeV
      maxETbits_     = 14;
      maxETFracbits_ = 8;
      maxZbits_	     = 12;
      maxEtabits_    = 10;
      maxPhibits_    = 9;	//signed

      maxSigmaET_  = 13.91526;  //in GeV
      maxSigmaZ_   = 778.098493;	//in mm
      maxSigmaPhi_ = 0.12822;
      maxSigmaRoZ_ = 0.024584;
      maxSigmaEta_ = 0.148922;

      maxSigmaETbits_  = 7;
      maxSigmaZbits_   = 7;
      maxSigmaPhibits_ = 7;
      maxSigmaRoZbits_ = 7;
      maxSigmaEtabits_ = 5;

      /*
      # Constants for GCT
      c_GCT_0 = 128 # placeholder (0.5)
      c_GCT_1 = 128 # placeholder (0.5)
      c_GCT_2 = 128 # placeholder (0.5)
      c_GCT_3 = 64 # placeholder (16 GeV)
      */

      c_Phi_Offset_ = 360 ;

      LSB_mean_roz_LUT_	 = 0.006374238205464194 ;
      LSB_sigma_roz_LUT_ = 0.000441991619167075 ;
      mu_roz_ds_lower_	 = 809.9324324324323924884083680808544158935546875 ;
      mu_roz_ds_upper_	 = 4996.798250844762151245959103107452392578125 ;

    }

    static float round(float frac, int precision){
      return float(int(frac*10*precision + 0.5))/float(10*precision);
    }

    void setMaxTCETbits(uint16_t maxTCETbits) { maxTCETbits_ = maxTCETbits ; }
    void setMaxTCRoZbits(uint16_t maxTCRoZbits) { maxTCRoZbits_ = maxTCRoZbits ; }

    float getZcm(uint16_t l)            const   { return (l>0 and l<48)?_zLayer[l-1]:0.0; }
    float roz_global_min()		const	{ return roz_global_min_; }
    float roz_global_max()		const	{ return roz_global_max_; }

    float eta_global_min()		const	{ return asinh(1/roz_global_max_); }
    float eta_global_max()		const	{ return asinh(1/roz_global_min_); }

    /*TC LSBs*/
    float LSB_E_TC()			const   { return maxTCET_/float(1<<maxTCETbits_); }
    float LSB_z_TC()			const	{ return LSB_z(); }
    float LSB_roz_TC()			const	{ return roz_global_max_/float(1<<maxTCRoZbits_); }
    float LSB_phi_TC()			const	{ return acos(-1)/float(maxTCHistoColNo_*(1<<5)); }
    //float LSB_phi_TC()			const	{ return 2*acos(-1)/float(1<<12); }

    /*Cluster LSBs*/
    float LSB_E()			const	{ return round(maxET_/float(1<<maxETbits_),2); }
    float LSB_Frac()			const	{ return 1/float(1<<maxETFracbits_); }
    float LSB_eta()			const	{ return acos(-1)/float(2*c_Phi_Offset_); }				//corresponding unsigned range 320 to 687 for (1.4<=eta<=3)
    float LSB_phi()			const	{ return acos(-1)/float(2*c_Phi_Offset_); }
    int maxPhiP()                       const   { return int(1<<(maxPhibits_-1))-1 ; }
    int maxPhiM()                       const   { return -1 * int(1<<(maxPhibits_-1)) ; }
    int nomPhiH()                       const   { return (2*c_Phi_Offset_/180)*60; } //max bins per 180 degree times 120/2
    int nomPhiL()                       const   { return -1 * ((2*c_Phi_Offset_/180)*60+1); } //max bins per 180 degree times 120/2

    float LSB_z()			const	{ return round((abs_z_global_max_-abs_z_global_min_)/float(1<<maxZbits_),1); }
    float LSB_roz()			const	{ return LSB_roz_TC(); }

    float c_Cluster_E_Sat()		const	{ return ((1<<maxETbits_)-1)*LSB_E()/LSB_E_TC(); }

    float LSB_sigma_E()			const	{ return maxSigmaET_/float(1<<maxSigmaETbits_); }
    float LSB_sigma_z()			const	{ return maxSigmaZ_/float(1<<maxSigmaZbits_); }
    float LSB_sigma_phi()		const	{ return maxSigmaPhi_/float(1<<maxSigmaPhibits_); }
    float LSB_sigma_roz()		const	{ return maxSigmaRoZ_/float(1<<maxSigmaRoZbits_); }
    float LSB_sigma_eta()		const	{ return maxSigmaEta_/float(1<<maxSigmaEtabits_); }

    uint16_t c_Phi_Offset()		const	{ return c_Phi_Offset_; } ;

    float eta_min_L1T()			const	{ return 320*LSB_phi(); }					//~1.4
    float eta_max_L1T()			const	{ return 687*LSB_phi(); }					//~3.0

    float roz_min_L1T()			const	{ return  1/sinh(eta_max_L1T()); }
    float roz_max_L1T()			const	{ return  1/sinh(eta_min_L1T()); }
    float LSB_roz_LUT()			const	{ return (roz_max_L1T()-roz_min_L1T())/float(1<<maxEtabits_); }

    float LSB_mean_roz_LUT()		const	{ return LSB_mean_roz_LUT_ ; }
    float LSB_sigma_roz_LUT()		const	{ return LSB_sigma_roz_LUT_; }
    float mu_roz_ds_lower()		const	{ return mu_roz_ds_lower_; }
    float mu_roz_ds_upper()		const	{ return mu_roz_ds_upper_; }

    float c_ET_scaler()			const	{ return LSB_E_TC()/LSB_E(); }				//L1T LSB Rounding
    float c_frac_scaler()		const	{ return 1/LSB_Frac(); }					//L1T LSB Rounding
    float c_roz_scaler_0()		const	{ return LSB_roz_TC()/LSB_roz_LUT(); }			// mu(eta) LUT
    float c_phi_scaler()		const	{ return LSB_phi_TC()/LSB_phi(); }			// L1T LSB Rounding
    float c_sigma_E_scaler()		const	{ return LSB_E_TC()/LSB_sigma_E(); }			// L1T LSB Rounding
    float c_sigma_z_scaler()		const	{ return LSB_z_TC()/LSB_sigma_z(); }			// L1T LSB Rounding
    float c_sigma_phi_scaler()		const	{ return LSB_phi_TC()/LSB_sigma_phi(); }			// L1T LSB Rounding
    float c_roz_scaler_1()		const	{ return LSB_roz_TC()/LSB_mean_roz_LUT(); }		// sigma(eta) LUT
    float c_sigma_roz_scaler_0()	const	{ return LSB_roz_TC()/LSB_sigma_roz_LUT(); }		// sigma(eta) LUT
    float c_sigma_roz_scaler_1()	const	{ return LSB_roz_TC()/LSB_sigma_roz(); }			// L1T LSB Rounding

    void print(){
      std::cout<< "roz_global_min :" << roz_global_min()
	       << ", roz_global_max :" << roz_global_max()
	       << ", eta_global_min :" << eta_global_min()
	       << ", eta_global_max :" << eta_global_max()
	       << std::endl
	       << "LSB_E_TC :" << LSB_E_TC()
	       << ", LSB_z_TC :" << LSB_z_TC()
	       << ", LSB_roz_TC :" << LSB_roz_TC()
	       << ", LSB_phi_TC :" << LSB_phi_TC()
	       << std::endl
	       << "LSB_E :" << LSB_E()
	       << ", LSB_Frac :" << LSB_Frac()
	       << ", LSB_eta :" << LSB_eta()
	       << ", LSB_phi :" << LSB_phi()
	       << ", LSB_z :" << LSB_z()
	       << ", LSB_roz :" << LSB_roz()
	       << std::endl
	       << "maxPhiP : " << maxPhiP()
	       << ", maxPhiM : " << maxPhiM()
	       << ", nomPhiL : " << nomPhiL()
	       << ", nomPhiH : " << nomPhiH()
	       << std::endl
	       << "c_Cluster_E_Sat :" << c_Cluster_E_Sat()
	       << ", c_Phi_Offset :" << c_Phi_Offset()
	       << std::endl
	       << "LSB_sigma_E :" << LSB_sigma_E()
	       << ", LSB_sigma_z :" << LSB_sigma_z()
	       << ", LSB_sigma_phi :" << LSB_sigma_phi()
	       << ", LSB_sigma_roz :" << LSB_sigma_roz()
	       << ", LSB_sigma_eta :" << LSB_sigma_eta()
	       << std::endl
	       << "eta_min_L1T :" << eta_min_L1T()
	       << ", eta_max_L1T :" << eta_max_L1T()
	       << ", roz_min_L1T :" << roz_min_L1T()
	       << ", roz_max_L1T :" << roz_max_L1T()
	       << ", LSB_roz_LUT :" << LSB_roz_LUT()
	       << std::endl
	       << "LSB_mean_roz_LUT :" << LSB_mean_roz_LUT()
	       << ", LSB_sigma_roz_LUT :" << LSB_sigma_roz_LUT()
	       << ", mu_roz_ds_lower :" << mu_roz_ds_lower()
	       << ", mu_roz_ds_upper :" << mu_roz_ds_upper()
	       << std::endl
	       << "c_ET_scaler :" << c_ET_scaler()
	       << ", c_frac_scaler :" << c_frac_scaler()
	       << ", c_roz_scaler_0 :" << c_roz_scaler_0()
	       << ", c_phi_scaler :" << c_phi_scaler()
	       << ", c_sigma_E_scaler :" << c_sigma_E_scaler()
	       << ", c_sigma_z_scaler :" << c_sigma_z_scaler()
	       << ", c_sigma_phi_scaler :" << c_sigma_phi_scaler()
	       << ", c_roz_scaler_1 :" << c_roz_scaler_1()
	       << ", c_sigma_roz_scaler_0 :" << c_sigma_roz_scaler_0()
	       << ", c_sigma_roz_scaler_1 :" << c_sigma_roz_scaler_1()
	       << std::endl;
    }

  private:
    //# FIRST PRINCIPLE GEOMETRY CONSTRAINTS. Coming from the HGCAL mechanical envelope.
    float abs_z_global_min_, abs_z_global_max_; //# measured from IP, CE-E first layer
    float r_global_min_, r_global_max_;         //# measured from beam axis
    float roz_global_min_, roz_global_max_;     //# coming from sensitive layer ranges (just outside of them)
    double _zLayer[47];

    //TC parameters
    float maxTCET_;           //in GeV
    uint16_t maxTCETbits_;
    uint16_t maxTCRoZbits_;
    uint16_t maxTCHistoColNo_;

    //Cluster parameters
    float maxET_ ;            //in GeV
    uint16_t maxETbits_;      //max allocated bits for ET
    uint16_t maxETFracbits_;  //max allocated bits for ET fractions
    uint16_t maxZbits_;
    uint16_t maxEtabits_;
    uint16_t maxPhibits_;

    //Sigma paramaters
    float  maxSigmaET_;
    float  maxSigmaZ_;
    float  maxSigmaPhi_;
    float  maxSigmaRoZ_;
    float  maxSigmaEta_;

    uint16_t maxSigmaETbits_;
    uint16_t maxSigmaZbits_;
    uint16_t maxSigmaPhibits_;
    uint16_t maxSigmaRoZbits_;
    uint16_t maxSigmaEtabits_;

    //Phi Offset
    uint16_t c_Phi_Offset_; //90 degrees equivalent at LSB = pi/720

    // # sig_eta Look-up control (placeholders)
    float LSB_mean_roz_LUT_ ;   // what was used during LUT formation, 6b mean addr
    float LSB_sigma_roz_LUT_ ;  //what was used during LUT formation, 6b sigma addr
    // #mu(roz) after dividing integer sums in EMulator (W16) appears in this float range:
    float mu_roz_ds_lower_ ;
    float mu_roz_ds_upper_ ;

  };

}
#endif
