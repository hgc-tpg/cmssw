#include "DataFormats/L1THGCal/interface/HGCalMulticluster.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "L1Trigger/L1THGCal/interface/HGCalTriggerGeometryBase.h"
#include "L1Trigger/L1THGCalUtilities/interface/HGCalTriggerNtupleBase.h"

class HGCalTriggerNtupleHGCMulticlustersSA : public HGCalTriggerNtupleBase {
public:
  HGCalTriggerNtupleHGCMulticlustersSA(const edm::ParameterSet& conf);
  ~HGCalTriggerNtupleHGCMulticlustersSA() override{};
  void initialize(TTree&, const edm::ParameterSet&, edm::ConsumesCollector&&) final;
  void fill(const edm::Event& e, const HGCalTriggerNtupleEventSetup& es) final;

private:
  void clear() final;

  edm::EDGetToken multiclusters_token_;

  int cl3d_n_;
  std::vector<uint32_t> cl3d_id_;
  std::vector<float> cl3d_pt_;
  std::vector<float> cl3d_eta_;
  std::vector<float> cl3d_phi_;
  std::vector<float> cl3d_energy_;
  // cluster shower shapes
  std::vector<int> cl3d_showerlength_;
  std::vector<int> cl3d_coreshowerlength_;
  std::vector<int> cl3d_firstlayer_;
  std::vector<unsigned long int> cl3d_Sigma_E_Quotient_;
  std::vector<unsigned long int> cl3d_Sigma_E_Fraction_;
  std::vector<unsigned long int> cl3d_Mean_z_Quotient_;
  std::vector<unsigned long int> cl3d_Mean_z_Fraction_;
  std::vector<unsigned long int> cl3d_Mean_phi_Quotient_;
  std::vector<unsigned long int> cl3d_Mean_phi_Fraction_;
  std::vector<unsigned long int> cl3d_Mean_eta_Quotient_;
  std::vector<unsigned long int> cl3d_Mean_eta_Fraction_;
  std::vector<unsigned long int> cl3d_Mean_roz_Quotient_;
  std::vector<unsigned long int> cl3d_Mean_roz_Fraction_;
  std::vector<unsigned long int> cl3d_Sigma_z_Quotient_;
  std::vector<unsigned long int> cl3d_Sigma_z_Fraction_;
  std::vector<unsigned long int> cl3d_Sigma_phi_Quotient_;
  std::vector<unsigned long int> cl3d_Sigma_phi_Fraction_;
  std::vector<unsigned long int> cl3d_Sigma_eta_Quotient_;
  std::vector<unsigned long int> cl3d_Sigma_eta_Fraction_;
  std::vector<unsigned long int> cl3d_Sigma_roz_Quotient_;
  std::vector<unsigned long int> cl3d_Sigma_roz_Fraction_;
  std::vector<unsigned long int> cl3d_E_EM_over_E_Quotient_;
  std::vector<unsigned long int> cl3d_E_EM_over_E_Fraction_;
  std::vector<unsigned long int> cl3d_E_EM_core_over_E_EM_Quotient_;
  std::vector<unsigned long int> cl3d_E_EM_core_over_E_EM_Fraction_;
  std::vector<unsigned long int> cl3d_E_H_early_over_E_Quotient_;
  std::vector<unsigned long int> cl3d_E_H_early_over_E_Fraction_;
};

DEFINE_EDM_PLUGIN(HGCalTriggerNtupleFactory,
                  HGCalTriggerNtupleHGCMulticlustersSA,
                  "HGCalTriggerNtupleHGCMulticlustersSA");

HGCalTriggerNtupleHGCMulticlustersSA::HGCalTriggerNtupleHGCMulticlustersSA(const edm::ParameterSet& conf)
    : HGCalTriggerNtupleBase(conf) {
  accessEventSetup_ = false;
}

void HGCalTriggerNtupleHGCMulticlustersSA::initialize(TTree& tree,
                                                      const edm::ParameterSet& conf,
                                                      edm::ConsumesCollector&& collector) {
  multiclusters_token_ =
      collector.consumes<l1t::HGCalMulticlusterBxCollection>(conf.getParameter<edm::InputTag>("Multiclusters"));

  std::string prefix(conf.getUntrackedParameter<std::string>("Prefix", "cl3dSA"));

  std::string bname;
  auto withPrefix([&prefix, &bname](char const* vname) -> char const* {
    bname = prefix + "_" + vname;
    return bname.c_str();
  });

  tree.Branch(withPrefix("n"), &cl3d_n_, (prefix + "_n/I").c_str());
  tree.Branch(withPrefix("id"), &cl3d_id_);
  tree.Branch(withPrefix("pt"), &cl3d_pt_);
  tree.Branch(withPrefix("eta"), &cl3d_eta_);
  tree.Branch(withPrefix("phi"), &cl3d_phi_);
  tree.Branch(withPrefix("energy"), &cl3d_energy_);
  tree.Branch(withPrefix("showerlength"), &cl3d_showerlength_);
  tree.Branch(withPrefix("coreshowerlength"), &cl3d_coreshowerlength_);
  tree.Branch(withPrefix("firstlayer"), &cl3d_firstlayer_);
  tree.Branch(withPrefix("meanZQuotient"), &cl3d_Mean_z_Quotient_);
  tree.Branch(withPrefix("meanZFraction"), &cl3d_Mean_z_Fraction_);
  tree.Branch(withPrefix("meanPhiQuotient"), &cl3d_Mean_phi_Quotient_);
  tree.Branch(withPrefix("meanPhiFraction"), &cl3d_Mean_phi_Fraction_);
  tree.Branch(withPrefix("meanEtaQuotient"), &cl3d_Mean_eta_Quotient_);
  tree.Branch(withPrefix("meanEtaFraction"), &cl3d_Mean_eta_Fraction_);
  tree.Branch(withPrefix("meanRoZQuotient"), &cl3d_Mean_roz_Quotient_);
  tree.Branch(withPrefix("meanRoZFraction"), &cl3d_Mean_roz_Fraction_);
  tree.Branch(withPrefix("sigmaEQuotient"), &cl3d_Sigma_E_Quotient_);
  tree.Branch(withPrefix("sigmaEFraction"), &cl3d_Sigma_E_Fraction_);
  tree.Branch(withPrefix("sigmaZQuotient"), &cl3d_Sigma_z_Quotient_);
  tree.Branch(withPrefix("sigmaZFraction"), &cl3d_Sigma_z_Fraction_);
  tree.Branch(withPrefix("sigmaPhiQuotient"), &cl3d_Sigma_phi_Quotient_);
  tree.Branch(withPrefix("sigmaPhiFraction"), &cl3d_Sigma_phi_Fraction_);
  tree.Branch(withPrefix("sigmaEtaQuotient"), &cl3d_Sigma_eta_Quotient_);
  tree.Branch(withPrefix("sigmaEtaFraction"), &cl3d_Sigma_eta_Fraction_);
  tree.Branch(withPrefix("sigmaRoZQuotient"), &cl3d_Sigma_roz_Quotient_);
  tree.Branch(withPrefix("sigmaRoZFraction"), &cl3d_Sigma_roz_Fraction_);
  tree.Branch(withPrefix("energyEMOverEnergyQuotient"), &cl3d_E_EM_over_E_Quotient_);
  tree.Branch(withPrefix("energyEMOverEnergyFraction"), &cl3d_E_EM_over_E_Fraction_);
  tree.Branch(withPrefix("energyEMCoreOverEnergyEMQuotient"), &cl3d_E_EM_core_over_E_EM_Quotient_);
  tree.Branch(withPrefix("energyEMCoreOverEnergyEMFraction"), &cl3d_E_EM_core_over_E_EM_Fraction_);
  tree.Branch(withPrefix("energyHEarlyOverEnergyQuotient"), &cl3d_E_H_early_over_E_Quotient_);
  tree.Branch(withPrefix("energyHEarlyOverEnergyFraction"), &cl3d_E_H_early_over_E_Fraction_);
}

void HGCalTriggerNtupleHGCMulticlustersSA::fill(const edm::Event& e, const HGCalTriggerNtupleEventSetup& es) {
  // retrieve clusters 3D
  edm::Handle<l1t::HGCalMulticlusterBxCollection> multiclusters_h;
  e.getByToken(multiclusters_token_, multiclusters_h);
  const l1t::HGCalMulticlusterBxCollection& multiclusters = *multiclusters_h;

  clear();
  for (auto cl3d_itr = multiclusters.begin(0); cl3d_itr != multiclusters.end(0); cl3d_itr++) {
    cl3d_n_++;
    cl3d_id_.emplace_back(cl3d_itr->detId());
    // physical values
    cl3d_pt_.emplace_back(cl3d_itr->pt());
    cl3d_eta_.emplace_back(cl3d_itr->eta());
    cl3d_phi_.emplace_back(cl3d_itr->phi());
    cl3d_energy_.emplace_back(cl3d_itr->energy());
    cl3d_showerlength_.emplace_back(cl3d_itr->showerLength());
    cl3d_coreshowerlength_.emplace_back(cl3d_itr->coreShowerLength());
    cl3d_firstlayer_.emplace_back(cl3d_itr->firstLayer());
    cl3d_Sigma_E_Quotient_.emplace_back(cl3d_itr->Sigma_E_Quotient());
    cl3d_Sigma_E_Fraction_.emplace_back(cl3d_itr->Sigma_E_Fraction());
    cl3d_Mean_z_Quotient_.emplace_back(cl3d_itr->Mean_z_Quotient());
    cl3d_Mean_z_Fraction_.emplace_back(cl3d_itr->Mean_z_Fraction());
    cl3d_Mean_phi_Quotient_.emplace_back(cl3d_itr->Mean_phi_Quotient());
    cl3d_Mean_phi_Fraction_.emplace_back(cl3d_itr->Mean_phi_Fraction());
    cl3d_Mean_eta_Quotient_.emplace_back(cl3d_itr->Mean_eta_Quotient());
    cl3d_Mean_eta_Fraction_.emplace_back(cl3d_itr->Mean_eta_Fraction());
    cl3d_Mean_roz_Quotient_.emplace_back(cl3d_itr->Mean_roz_Quotient());
    cl3d_Mean_roz_Fraction_.emplace_back(cl3d_itr->Mean_roz_Fraction());
    cl3d_Sigma_z_Quotient_.emplace_back(cl3d_itr->Sigma_z_Quotient());
    cl3d_Sigma_z_Fraction_.emplace_back(cl3d_itr->Sigma_z_Fraction());
    cl3d_Sigma_phi_Quotient_.emplace_back(cl3d_itr->Sigma_phi_Quotient());
    cl3d_Sigma_phi_Fraction_.emplace_back(cl3d_itr->Sigma_phi_Fraction());
    cl3d_Sigma_eta_Quotient_.emplace_back(cl3d_itr->Sigma_eta_Quotient());
    cl3d_Sigma_eta_Fraction_.emplace_back(cl3d_itr->Sigma_eta_Fraction());
    cl3d_Sigma_roz_Quotient_.emplace_back(cl3d_itr->Sigma_roz_Quotient());
    cl3d_Sigma_roz_Fraction_.emplace_back(cl3d_itr->Sigma_roz_Fraction());
    cl3d_E_EM_over_E_Quotient_.emplace_back(cl3d_itr->E_EM_over_E_Quotient());
    cl3d_E_EM_over_E_Fraction_.emplace_back(cl3d_itr->E_EM_over_E_Fraction());
    cl3d_E_EM_core_over_E_EM_Quotient_.emplace_back(cl3d_itr->E_EM_core_over_E_EM_Quotient());
    cl3d_E_EM_core_over_E_EM_Fraction_.emplace_back(cl3d_itr->E_EM_core_over_E_EM_Fraction());
    cl3d_E_H_early_over_E_Quotient_.emplace_back(cl3d_itr->E_H_early_over_E_Quotient());
    cl3d_E_H_early_over_E_Fraction_.emplace_back(cl3d_itr->E_H_early_over_E_Fraction());
  }
}

void HGCalTriggerNtupleHGCMulticlustersSA::clear() {
  cl3d_n_ = 0;
  cl3d_id_.clear();
  cl3d_pt_.clear();
  cl3d_eta_.clear();
  cl3d_phi_.clear();
  cl3d_energy_.clear();
  cl3d_showerlength_.clear();
  cl3d_coreshowerlength_.clear();
  cl3d_firstlayer_.clear();
  cl3d_Sigma_E_Quotient_.clear();
  cl3d_Sigma_E_Fraction_.clear();
  cl3d_Mean_z_Quotient_.clear();
  cl3d_Mean_z_Fraction_.clear();
  cl3d_Mean_phi_Quotient_.clear();
  cl3d_Mean_phi_Fraction_.clear();
  cl3d_Mean_eta_Quotient_.clear();
  cl3d_Mean_eta_Fraction_.clear();
  cl3d_Mean_roz_Quotient_.clear();
  cl3d_Mean_roz_Fraction_.clear();
  cl3d_Sigma_z_Quotient_.clear();
  cl3d_Sigma_z_Fraction_.clear();
  cl3d_Sigma_phi_Quotient_.clear();
  cl3d_Sigma_phi_Fraction_.clear();
  cl3d_Sigma_eta_Quotient_.clear();
  cl3d_Sigma_eta_Fraction_.clear();
  cl3d_Sigma_roz_Quotient_.clear();
  cl3d_Sigma_roz_Fraction_.clear();
  cl3d_E_EM_over_E_Quotient_.clear();
  cl3d_E_EM_over_E_Fraction_.clear();
  cl3d_E_EM_core_over_E_EM_Quotient_.clear();
  cl3d_E_EM_core_over_E_EM_Fraction_.clear();
  cl3d_E_H_early_over_E_Quotient_.clear();
  cl3d_E_H_early_over_E_Fraction_.clear();
}
