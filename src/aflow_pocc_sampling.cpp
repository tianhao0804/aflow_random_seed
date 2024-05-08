// ***************************************************************************
// *                                                                         *
// *       AFlow COREY OSES Johns Hopkins University 2012-2024               *
// *       AFlow Yuxiang Liu Johns Hopkins University 2024                   *
// *                                                                         *                         
// ***************************************************************************
//aflow_pocc_sampling.cpp
//Reading Pocc sampling rate from --pocc_sample_rate and --pocc_sample_number
//Two random seed sampling functions are provided for configuration selections.
//First random seed function is for unique configuration calculations.
//Second random seed function is for DFT calculations.
//2024- yliu597@jh.edu


#ifndef _AFLOW_POCC_SAMPLING_CPP_
#define _AFLOW_POCC_SAMPLING_CPP_

#include "aflow.h"
#include "aflow_pocc.h"
#include <unordered_set>

#define _DEBUG_POCC_SAMPLE_ false 

namespace pocc {
    bool XHOST_POCC_SAMPLE_RATE=XHOST.vflag_pflow.flag("POCC_SAMPLE_RATE");
    bool XHOST_POCC_SAMPLE_NUMBER=XHOST.vflag_pflow.flag("POCC_SAMPLE_NUMBER");

    double POccCalculator::setPOccSampleRate(const string& pocc_sample_rate_string, int sample_round){//give the warning message if --pocc_sample_rate format is wrong
        double pocc_sample_rate=0.0;
        stringstream message;
        //if(pocc_sample_rate_string.empty()){
        //   message << "No POCC random seed sampling rate is found" << endl;
        //   throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_INPUT_ILLEGAL_);
        //}   
        if(aurostd::substring2bool(pocc_sample_rate_string,",")){
           message << "Cannot handle more than one pocc_sample_rate specification";
           throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_INPUT_ILLEGAL_);
        }   
 
        vector<string> tokens;
        aurostd::string2tokens(pocc_sample_rate_string,tokens,":");
        if(tokens.size()<2){
           message << "No second random seed sampling rate is found" << endl;
           throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_INPUT_ILLEGAL_);
        }   
        else{
           if(!aurostd::isfloat(tokens[sample_round])){
             message << "Input is not a float [input=" << tokens[sample_round] << "]";
             throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_INPUT_ILLEGAL_);
           }   
           pocc_sample_rate=aurostd::string2utype<double>(tokens[sample_round]);
           if(aurostd::isequal(pocc_sample_rate,0.0,_AFLOW_POCC_ZERO_TOL_)){
             message << "Sample rate is too small (0) [input=" << tokens[sample_round] << "]";
             throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_INPUT_ILLEGAL_);
           }   
           if(std::signbit(pocc_sample_rate)){
              message << "Sample rate cannot be negative [input=" << tokens[sample_round] << "]";
              throw aurostd::xerror(__AFLOW_FILE__,__AFLOW_FUNC__,message,_INPUT_ILLEGAL_);
           }   
        }  
        return pocc_sample_rate;
    }

    vector<unsigned long long int> POccCalculator::FirstRandomSampling(unsigned long long int hnf_count, unsigned long long int types_config_permutations_count){
         size_t sample_config_count=0;
         bool LDEBUG = (FALSE || _DEBUG_POCC_SAMPLE_ || XHOST.DEBUG);
         if(LDEBUG){
             cerr << __AFLOW_FUNC__ << " This is a test message" << endl;
             cerr << __AFLOW_FUNC__ << " This is first random seed sampling for unique configuration calculations" << endl;
         }
         int sample_seed = DEFAULT_POCC_SAMPLE_SEED;
         srand(sample_seed); //initial random seed sampling
         if(XHOST_POCC_SAMPLE_RATE){
             double sample_rate = setPOccSampleRate(XHOST.vflag_pflow.getattachedscheme("POCC_SAMPLE_RATE"),0);
             sample_config_count = static_cast<size_t>(types_config_permutations_count * sample_rate);
         } // calculate total amount of selected decoration permutation configurations for each derivative supperlattice
         if(XHOST_POCC_SAMPLE_NUMBER){sample_config_count = static_cast<size_t>(aurostd::string2utype<unsigned long long int>(XHOST.vflag_pflow.getattachedscheme("POCC_SAMPLE_NUMBER"),0));} // calculate total amount of selected decoration permutation configurations for each derivative supperlattice
         if(LDEBUG){cerr << __AFLOW_FUNC__ << " aflow_pocc_random_seed: " << "sample_config_count" << sample_config_count << endl;}
         std::unordered_set<unsigned long long int> selected_site_config_indices_set; // unordered_set could auotmatically select unique index avoiding repetion and is cheaper than find
         while(selected_site_config_indices_set.size() < sample_config_count) {
             unsigned long long int index = rand() % types_config_permutations_count; //select permuation decoration configuration index
             selected_site_config_indices_set.insert(index);
         }
         std::vector<unsigned long long int> selected_site_config_indices(selected_site_config_indices_set.begin(), selected_site_config_indices_set.end());
         sort(selected_site_config_indices.begin(), selected_site_config_indices.end());//set selected_indices to save sorted indices
         if(LDEBUG)
         {
            cerr << __AFLOW_FUNC__ << "selected_site_config_indices.size: " << selected_site_config_indices.size() << endl;
            cerr << __AFLOW_FUNC__ << "selected_site_config_indices.size: " << selected_site_config_indices.size() << endl;
         }
         if(LDEBUG){
            for(size_t i=0; i< selected_site_config_indices.size(); i++)
            {
                cerr << __AFLOW_FUNC__ << " This is a test message" << endl;
                cerr << __AFLOW_FUNC__ << " Selected permutation configurations in each derivative superlattice: C" << selected_site_config_indices[i] << endl;
                cerr << __AFLOW_FUNC__ << " sample_config_count" <<  sample_config_count << endl;
                cerr << __AFLOW_FUNC__ << "selected_site_config_indices.size: " << selected_site_config_indices.size() << endl;
            }
         }
         std::vector<unsigned long long int> selected_indices_set;
         for(size_t i = 0; i < hnf_count; i++){
             for(size_t j : selected_site_config_indices){
                 selected_indices_set.push_back(i * types_config_permutations_count + j);
                 if(LDEBUG)
                 {
                    cerr << __AFLOW_FUNC__ << "selected_indices_set.size: " << selected_indices_set.size() << endl;
                    cerr << __AFLOW_FUNC__ << "i_j = " << i << j  << endl;
                 }
             }
         }
         if(LDEBUG){cerr << __AFLOW_FUNC__ << "size of selected_indices_set = " << selected_indices_set.size() << endl;}
         return selected_indices_set;
    }
    
    std::list<pocc::POccSuperCellSet> POccCalculator::SecondRandomSamplingWithRate(std::list<pocc::POccSuperCellSet> l_supercell_sets, unsigned long long int hnf_count){
        bool LDEBUG = (FALSE || _DEBUG_POCC_SAMPLE_ || XHOST.DEBUG);
        double sample_rate = setPOccSampleRate(XHOST.vflag_pflow.getattachedscheme("POCC_SAMPLE_RATE"),1);
        vector<unsigned long long int> ihnf_unique_supercell_num(hnf_count,0);
        l_supercell_sets.sort();
        for(unsigned long long int i=0;i<l_supercell_sets.size();i++) {
            std::list<pocc::POccSuperCellSet>::iterator it = l_supercell_sets.begin();
            std::advance(it,i);
            pocc::POccSuperCell pscs=(*it).getSuperCell();
            ihnf_unique_supercell_num[pscs.m_hnf_index]++; // calculate unique decoration configurations for each superlattices and saved in ihnf_unique_supercell_num
            if(LDEBUG){
               cerr << __AFLOW_FUNC__ << " hnf_index: " << pscs.m_hnf_index << endl;
               cerr << __AFLOW_FUNC__ << " hnf_index: " <<  pscs.m_hnf_index << " hnf_unique_num: " << ihnf_unique_supercell_num[pscs.m_hnf_index] << endl;

               cerr << __AFLOW_FUNC__ << " H" << pscs.m_hnf_index << "C" << pscs.m_site_config_index << endl;
            }
        }
        if(LDEBUG){
            cerr << __AFLOW_FUNC__ << " This is a test message" << endl;
            cerr << __AFLOW_FUNC__ << " This is second random seed sampling for DFT calculations" << endl;
            cerr << __AFLOW_FUNC__ << " Random seed sampling rate is: " << 100*sample_rate << "%" << endl;
        } 
        std::list<pocc::POccSuperCellSet> selected_random_sampling_supercell_sets;
        for(unsigned long long int i=0; i<hnf_count; i++){ // loop for each superlattice to generate the sample index for each superlattice
            vector<int> sample_index(ihnf_unique_supercell_num[i]); //sample_index save the unique configuration indices for each superlattice
            if(i==0){
               iota(sample_index.begin(),sample_index.end(),0); //iota create integer vector from 0 for 1st superlattice
            }else{
               unsigned long long int summation=0;
               for(unsigned long long int j=0; j<i; j++){ summation+=ihnf_unique_supercell_num[j]; }
               if(LDEBUG) {cerr << __AFLOW_FUNC__ << "summation: " << summation << endl;}
               iota(sample_index.begin(),sample_index.end(),summation); //iota create integer vector from summation for the other supperlattice exclude from 1st
            }
            int sample_seed = DEFAULT_POCC_SAMPLE_SEED;
            srand(sample_seed); //initial random seed sampling
            size_t sample_config_count = static_cast<size_t>(ihnf_unique_supercell_num[i]*sample_rate); //calculate total amount of selected configurations for each superlattice
            if(LDEBUG)
            {
               cerr << __AFLOW_FUNC__ << " ihnf_nique_supercell_num[" << i << "] = " << ihnf_unique_supercell_num[i]<< endl;
               cerr << __AFLOW_FUNC__ << " sample_config_count[" << i << "] = "  << sample_config_count << endl;
            }
            std::unordered_set<unsigned long long int> selected_config_indices_set; // unordered_set could auotmatically select unique index avoiding repetion and is cheaper than find
            while(selected_config_indices_set.size() < sample_config_count) {
                unsigned long long int index = rand() % ihnf_unique_supercell_num[i] ; //select permuation decoration configuration index
                selected_config_indices_set.insert(index);
            }
            std::vector<unsigned long long int> selected_config_indices(selected_config_indices_set.begin(), selected_config_indices_set.end());
            sort(selected_config_indices.begin(), selected_config_indices.end());
            for(size_t i=0; i< sample_config_count; i++){
                std::list<pocc::POccSuperCellSet>::iterator it = l_supercell_sets.begin();
                std::advance(it,sample_index[selected_config_indices[i]]);
                selected_random_sampling_supercell_sets.push_back(*it);//push back the selected supercell
                if(LDEBUG){
                   pocc::POccSuperCell pscs=(*it).getSuperCell();
                   cerr << __AFLOW_FUNC__ << " This is a test message" << endl;
                   cerr << __AFLOW_FUNC__ << " Second sampling is: " << "H" << pscs.m_hnf_index << "C" << pscs.m_site_config_index << endl;
                }
            }
         }
         return selected_random_sampling_supercell_sets;
    }

    std::list<pocc::POccSuperCellSet> POccCalculator::SecondRandomSampling(std::list<pocc::POccSuperCellSet> l_supercell_sets){
        bool LDEBUG = (FALSE || _DEBUG_POCC_SAMPLE_ || XHOST.DEBUG);
        std::unordered_set<unsigned long long int> selected_config_indices_set; // unordered_set could auotmatically select unique index avoiding repetion and is cheaper than find
        unsigned long long int sample_number;
        if(XHOST_POCC_SAMPLE_RATE){
           double sample_rate = setPOccSampleRate(XHOST.vflag_pflow.getattachedscheme("POCC_SAMPLE_RATE"),1);
           sample_number = static_cast<size_t>(l_supercell_sets.size() * sample_rate);
        } 
        if(XHOST_POCC_SAMPLE_NUMBER){sample_number = static_cast<size_t>(aurostd::string2utype<unsigned long long int>(XHOST.vflag_pflow.getattachedscheme("POCC_SAMPLE_NUMBER"),0));} // calculate total amount   of selected decoration permutation configurations for each derivative supperlattice

        int sample_seed = DEFAULT_POCC_SAMPLE_SEED;
        srand(sample_seed); //initial random seed sampling
        l_supercell_sets.sort();
        while(selected_config_indices_set.size() < sample_number) {
            unsigned long long int index = rand() % l_supercell_sets.size(); //select random number from whole l_supercell_sets
            selected_config_indices_set.insert(index);
        }
        std::vector<unsigned long long int> selected_config_indices(selected_config_indices_set.begin(), selected_config_indices_set.end());
        sort(selected_config_indices.begin(), selected_config_indices.end());//set selected_indices to save sorted indices
        if(LDEBUG){cerr << __AFLOW_FUNC__ << " selected_config_indices.size: " << selected_config_indices.size() << endl;}
        std::list<pocc::POccSuperCellSet> selected_random_sampling_supercell_sets;
        for(size_t i = 0; i < sample_number; i++){
            std::list<pocc::POccSuperCellSet>::iterator it = l_supercell_sets.begin();
            std::advance(it,selected_config_indices[i]);
            selected_random_sampling_supercell_sets.push_back(*it);//push back the selected supercell
        }

         return selected_random_sampling_supercell_sets;
    }

}
#endif
