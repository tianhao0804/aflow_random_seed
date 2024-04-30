// ***************************************************************************
// * *
// * AFlow COREY OSES Johns Hopkins University 2012-2024 *
// * AFlow Yuxiang Liu Johns Hopkins University 2024 *
// * *
// ***************************************************************************
// aflow_pocc_sampling.h and aflow_random_pocc_sampling.cpp*

#ifndef _AFLOW_POCC_SAMPLING_H_
#define _AFLOW_POCC_SAMPLING_H_
namespace pocc {
    vector<unsigned long long int> FirstRandomSampling(unsigned long long int hnf_count, unsigned long long int types_config_permutations_count);
    std::list<pocc::POccSuperCellSet> SecondRandomSamplingWithRate(std::list<pocc::POccSuperCellSet> l_supercell_sets, unsigned long long int hnf_count);
    std::list<pocc::POccSuperCellSet> SecondRandomSampling(std::list<pocc::POccSuperCellSet> l_supercell_sets);
    double setPOccSampleRate(const string& pocc_sample_rate_string, int sample_round);   
}

#endif
