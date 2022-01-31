#ifndef __L1Trigger_L1THGCal_HGCalStage1SortingAlg_SA_h__
#define __L1Trigger_L1THGCal_HGCalStage1SortingAlg_SA_h__

#include <vector>
#include <memory>
#include <cstdint>  // unsigned

typedef std::vector<unsigned> datain_t;       // elements to sort
typedef std::vector<unsigned> datasorter_t;   // elements to process in sorter
typedef std::vector<unsigned> datamergera_t;  // elements to process in merger a
typedef std::vector<unsigned> datamergerb_t;  // elements to process in merger b
typedef std::vector<unsigned> dataout_t;      // elements selected

typedef std::vector<unsigned> adressin_t;       //elements addresses to sort
typedef std::vector<unsigned> adresssorter_t;   //elements addresses to process in sorter
typedef std::vector<unsigned> adressmergera_t;  //elements addresses to process in merger a
typedef std::vector<unsigned> adressmergerb_t;  //elements addresses to process in merger b
typedef std::vector<unsigned> adressout_t;      //elements addresses selected

typedef std::vector<unsigned> data_to_sort_t;   //Data Pointer for sorter function
typedef std::vector<unsigned> data_to_merge_t;  //Data Pointer for merger function
typedef std::vector<unsigned> adress_t;         //Address Pointer for sorter & merger functions

namespace l1thgcfirmware {

  class HGCalStage1SortingAlg_SA {
  public:
    HGCalStage1SortingAlg_SA();

    HGCalStage1SortingAlg_SA(const unsigned NTCin, const unsigned NTCout) : N(NTCin), M(NTCout) {
      unsigned ns = NTCin / 4;
      unsigned nma = (NTCout < (NTCin / 4)) ? 2 * NTCout : NTCin / 2;
      unsigned nmb = (NTCout < nma) ? 2 * NTCout : 2 * nma;
      setNS(ns);
      setNMA(nma);
      setNMB(nmb);
    }

    ~HGCalStage1SortingAlg_SA() {}

    //Generates the sorting network
    void sorting(datain_t& arr_input, dataout_t& arr_output, adressout_t& arr_adresses) const;

    void setNS(const unsigned NSval) { NS = NSval; }
    void setNMA(const unsigned NMAval) { NMA = NMAval; }
    void setNMB(const unsigned NMBval) { NMB = NMBval; }

  private:
    //Sorter algorithm
    void sorter(data_to_sort_t& list_sorter, adress_t& list_adresses) const;
    //Merger algorithm A
    void mergerA(data_to_merge_t& list_merger, adress_t& list_adresses) const;
    //Merger algorithm B
    void mergerB(data_to_merge_t& list_merger, adress_t& list_adresses) const;
    //Calculate base 2 logarithm , and round to next integer
    unsigned loop_bound_for_mergera(void) const;
    unsigned loop_bound_for_mergerb(void) const;
    unsigned loop_bound_for_sorter(void) const;
    unsigned log2_rounded(unsigned i) const;

    const unsigned N;  // number of total elements to sort
    const unsigned M;  // number of elements to select
    unsigned NS;       // number of elements to sort sequentially (block sorter input size)
    unsigned NMA;      // number of elements to merge in block merger a (block merger a input size)
    unsigned NMB;      // number of elements to merge in block merger b (block merger b input size)
  };
}  // namespace l1thgcfirmware

#endif
