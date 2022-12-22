/*************************************************************************************************/
// C++ adaptation (w/ minimal changes) of C code from:
// ******************************
// Author L. Pacheco Rodriguez, laura.pacheco-rodriguez@llr.in2p3.fr
// Collaborators: T. Romanteau, J.B. Sauvan, F. Thiant
// Sept-2018
// Laboratoire Leprince-Ringuet - CNRS
// Universite Paris-Saclay
//************************************************************************************************/

#include "L1Trigger/L1THGCal/interface/backend_emulator/HGCalStage1SortingAlg_SA.h"

#include <iostream>

using namespace l1thgcfirmware;

//****************************************TOP LEVEL FUNCTION****************************************
//-------------------------------------------------------------------------------------------------/
// It performs Sorting and selection: sorts N input energies, and selects the M highest
//
// @param arr_input - list of elements to sort (input)
// @param arr_output - list of elements sorted (output)
// @param arr_adresses - list of adresses of elements selected (output). Address goes from 0 to N
//
//-------------------------------------------------------------------------------------------------/
void HGCalStage1SortingAlg_SA::sorting(const datain_t& arr_input, dataout_t& arr_output, adressout_t& arr_adresses) const {
  //Variables declaration
  datasorter_t arr_tmp(NS);
  datamerger_t arr_tmp_so(NMA), arr_tmp_so1(NMA), arr_tmp_so2(NMA);
  datamerger_t arr_tmp_mo(NMB);
  adresssorter_t arr_tmp_adresses(NS);
  adressmerger_t arr_tmp_so_adresses(NMA), arr_tmp_so1_adresses(NMA), arr_tmp_so2_adresses(NMA);
  adressmerger_t arr_tmp_mo_adresses(NMB);

  unsigned cs;  //counter for sorter blocks
  unsigned cm;  //counter for merger blocks
  unsigned i;   //Index for loops

  //Generate four "sorter blocks" of NS inputs each-----------------------------------------------
  for (cs = 0; cs < 4; cs++) {
    //Replica of entry array. Important for a correct HW synthesis
    //It builds a list of NS elements
    for (i = 0; i < NS; i++) {
      arr_tmp[i] = arr_input[i + cs * NS];
      arr_tmp_adresses[i] = i + cs * NS;
    }

    //Call to function that sorts NS inputs
    sorter(arr_tmp, arr_tmp_adresses);

    //Arrange 4 sorted lists in 2 arrays of NMA elements each
    if (cs == 0 || cs == 1) {
      for (i = 0; i < NMA / 2; i++) {
        arr_tmp_so1[i * 2 + cs] = arr_tmp[i];
        arr_tmp_so1_adresses[i * 2 + cs] = arr_tmp_adresses[i];
      }
    } else if (cs == 2 || cs == 3) {
      for (i = 0; i < NMA / 2; i++) {
        arr_tmp_so2[(i - 1) * 2 + cs] = arr_tmp[i];
        arr_tmp_so2_adresses[(i - 1) * 2 + cs] = arr_tmp_adresses[i];
      }
    }

  }  //End of block sorter---------------------------------------------------------------------------

  //Generate two "merger blocks" of NMA inputs each-------------------------------------------------
  for (cm = 0; cm < 2; cm++) {
    //It builds a list of NMA inputs
    if (cm == 0) {
      for (i = 0; i < NMA; i++) {
        arr_tmp_so[i] = arr_tmp_so1[i];
        arr_tmp_so_adresses[i] = arr_tmp_so1_adresses[i];
      }
    } else if (cm == 1) {
      for (i = 0; i < NMA; i++) {
        arr_tmp_so[i] = arr_tmp_so2[i];
        arr_tmp_so_adresses[i] = arr_tmp_so2_adresses[i];
      }
    }

    //Call to function that merges NMA inputs
    merger(NMA, arr_tmp_so, arr_tmp_so_adresses);

    //Arrange 2 sorted lists in 1 array of NMB elements
    if (cm == 0) {
      for (i = 0; i < NMB / 2; i++) {
        arr_tmp_mo[i * 2] = arr_tmp_so[i];
        arr_tmp_mo_adresses[i * 2] = arr_tmp_so_adresses[i];
      }
    } else if (cm == 1) {
      for (i = 0; i < NMB / 2; i++) {
        arr_tmp_mo[i * 2 + 1] = arr_tmp_so[i];
        arr_tmp_mo_adresses[i * 2 + 1] = arr_tmp_so_adresses[i];
      }
    }

  }  //End of block for merging-------------------------------------------------------------------

  //Generates 1 "merger block" of NMB inputs ----------------------------------------------------
  //Call to function to merge
  merger(NMB, arr_tmp_mo, arr_tmp_mo_adresses);

  //Replica of output array. Important for a correct HW synthesis
  for (i = 0; i < M; i++) {
    arr_output[i] = arr_tmp_mo[i];
    arr_adresses[i] = arr_tmp_mo_adresses[i];
  }
}

//*********************************************FUNCTIONS*******************************************
//-------------------------------------------------------------------------------------------------/
// Sorter algorithm of NS inputs
// Based on odd-even Batcher algorithm
// For more details, see:
// The Art Of Computer Programming, D. Knuth - Sorting and Searching (2nd edition Volume 3), 1998, page 111
//-------------------------------------------------------------------------------------------------/

void HGCalStage1SortingAlg_SA::sorter(data_to_sort_t& list_sorter, adress_t& list_adresses) const {
  unsigned tmp_data, tmp_adress;
  unsigned i, p, q, r, d, constant_bound;
  datasorter_t list_sorter_tmp(NS);
  adresssorter_t list_adresses_tmp(NS);

  //Replica of entry array. Important for a correct HW synthesis
  for (i = 0; i < NS; i++) {
    list_sorter_tmp[i] = list_sorter[i];
    list_adresses_tmp[i] = list_adresses[i];
  }

  //Calculate bound for loop, which is 2^[log2(NS)]
  constant_bound = log2_rounded(NS);

  //Performed for p=2^(n-1), 2^(n-2), ... , 1
  for (p = constant_bound, q = constant_bound, r = 0, d = p; p >= 1; p = p / 2, q = constant_bound, r = 0, d = p) {
    //Performed for q=2^(n-1), 2^(n-2), ... , p
    for (; q >= p; d = q - p, r = p, q = q / 2) {
      //counter_levels++;

      //Parcour the list_sorter of elements to do exchanges
      for (i = 0; i < (NS - d); i++) {
        if ((i & p) == r) {
          //counter_comparators++;

          //Compare & Exchange
          if (list_sorter_tmp[i] < list_sorter_tmp[i + d]) {
            tmp_data = list_sorter_tmp[i];
            list_sorter_tmp[i] = list_sorter_tmp[i + d];
            list_sorter_tmp[i + d] = tmp_data;

            tmp_adress = list_adresses_tmp[i];
            list_adresses_tmp[i] = list_adresses_tmp[i + d];
            list_adresses_tmp[i + d] = tmp_adress;
          }
        }

      }  //Label2
    }    //Label1
  }      //Label0

  //Replica of output array. Important for correct HW synthesis
  for (i = 0; i < NS; i++) {
    list_sorter[i] = list_sorter_tmp[i];
    list_adresses[i] = list_adresses_tmp[i];
  }
}

//-------------------------------------------------------------------------------------------------/
// Merger algorithm of NM inputs
// Based on odd-even Batcher algorithm
// For more details, see:
// The Art Of Computer Programming, D. Knuth - Sorting and Searching (2nd edition Volume 3), 1998, page 111
//
//-------------------------------------------------------------------------------------------------/
void HGCalStage1SortingAlg_SA::merger(const unsigned& NM, data_to_merge_t& list_merger, adress_t& list_adresses) const {
  unsigned tmp, tmp_adress;
  unsigned i, p, q, r, d;
  datamerger_t list_merger_tmp(NM);
  adressmerger_t list_adresses_tmp(NM);

  //Replica of entry array. Important for correct HW synthesis
  for (i = 0; i < NM; i++) {
    list_merger_tmp[i] = list_merger[i];
    list_adresses_tmp[i] = list_adresses[i];
  }

  //Variables initialization for nested loops
  p = 1;
  q = log2_rounded(NM);
  r = 0;
  d = p;

  //Performed for q=2^(n-1), 2^(n-2), ... , p
  for (; q >= p; d = q - p, r = p, q = q / 2) {
    //counter_levels++;

    //Parcour the list_sorter of elements to do exchanges
    for (i = 0; i < (NM - d); i++) {
      if ((i & p) == r) {
        //counter_comparators++;

        //Compare & Exchange
        if (list_merger_tmp[i] < list_merger_tmp[i + d]) {
          tmp = list_merger_tmp[i];
          list_merger_tmp[i] = list_merger_tmp[i + d];
          list_merger_tmp[i + d] = tmp;

          tmp_adress = list_adresses_tmp[i];
          list_adresses_tmp[i] = list_adresses_tmp[i + d];
          list_adresses_tmp[i + d] = tmp_adress;
        }
      }

    }  //Label2
  }    //Label1

  //Replica of output array. Important for correct HW synthesis
  for (i = 0; i < NM; i++) {
    list_merger[i] = list_merger_tmp[i];
    list_adresses[i] = list_adresses_tmp[i];
  }
}

//-------------------------------------------------------------------------------------------------/
//  Calculate log base 2 of i and round the result to next integer
//  Equivalent to res=ceil(log(i)/log(2))
//
//-------------------------------------------------------------------------------------------------/
unsigned HGCalStage1SortingAlg_SA::log2_rounded(unsigned i) const {
  unsigned res = 0;
  unsigned tmp = i;

  //Shift to right untill i == 1 for searching position of MSB
  while (i >>= 1) {
    ++res;
  }

  //If 2^res < i, round
  if ((unsigned)1 << res < tmp) {
    res++;
  }

  return res;
}
