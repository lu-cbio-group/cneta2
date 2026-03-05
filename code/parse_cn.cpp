#include "parse_cn.hpp"


// collection of functions for reading/writing copy number data


/** 
 * @brief Convert a state (0-based index of copy number combinations used for rate matrix) to total copy number.    
 * @param state The state index representing a specific copy number combination.
 * @param cn_max The maximum copy number allowed by the program.
 * @return The total copy number corresponding to the given state.
 */
int state_to_total_cn(int state, int cn_max){
    int cn = 0;
    vector<int> sums;
    // cout << "Sums of state: ";

    for(int i = 1; i <= (cn_max + 1); i++){
        int s = i * (i + 1) / 2;
        sums.push_back(s);
        // cout << "\t" << s;
    }
    // cout << endl;

    if(state < sums[0]) return 0;

    int i = 1;
    do{
        // the start and end index of copy number i
        if(state >= sums[i - 1] && state < sums[i]) return i;
        i++;
    }while(i <= cn_max);

    return 0;
}


void state_to_allele_cn(int state, int cn_max, int& cnA, int& cnB){
    int cn = 0;
    vector<int> sums;
    // cout << "Sums of state: ";
    for(int i = 1; i <= (cn_max + 1); i++){
        int s = i * (i + 1) / 2;
        sums.push_back(s);
        // cout << "\t" << s;
    }
    // cout << endl;

    if(state < sums[0]) return;

    int i = 1;
    do{
        if(state >= sums[i - 1] && state < sums[i]){
             // total copy number is i
             int diff = state - sums[i-1];
             cnA = diff;
             cnB = i - cnA;
        }
        i++;
    }while(i <= cn_max);
}


int allele_cn_to_state(int cnA, int cnB){
    int tcn = cnA + cnB;
    int s = 0;

    int nprev = 0;
    // There are i+1 combinations for a total copy number of i
    for(int i = 0; i < tcn; i++){
        nprev += i + 1;
    }
    // cout << nprev << " cases before " << cnA << "," << cnB << endl;
    s = nprev;

    for(int j = 0; j < cnA; j++){
        // cout << j << endl;
        s += 1;
    }
    // cout << "State is " << s << endl;
    return s;
}



/** 
 * @brief Convert a change state (0-based index of copy number combinations used for rate matrix) to total copy number.    
 * @param state The state index representing a specific copy number combination.
 * @param cn_max The maximum copy number allowed by the program.
 * @return The total copy number corresponding to the given state.
 * When max_haplotype_change = 1, 9 states
 * 0	1	2	3	4	5	6	7	8
// 0	1	1	2	2	2	3	3	4
// 0/0	0/1	1/0	0/2	2/0	 1/1	 1/2	 2/1	 2/2
// -2	-1	-1	0	0	0	+1	+1	+2
// -1/-1	-1/0	0/-1	-1/+1	+1/-1	0/0	0/+1	+1/0	+1/+1
 * When max_haplotype_change = 2, 16 states
// 0	1	2	3	4	5	6	7	8	9	10	11	12	13	14	15
// 0	1	1	2	2	2	3	3	3	3	4	4	4	5	5	6
// 0/0	0/1	1/0	0/2	2/0	 1/1	 1/2	 2/1	0/3	3/0	  2/2	 1/3	 3/1	2/3	3/2	3/3
// -1/-1	-1/0	0/-1	-1/+1	+1/-1	0/0	0/+1	+1/0	-1/+2	+1/+1	+1/+1	0/+2	+2/0	+1/+2	+2/+1	+2/+2
 */
int change_state_to_total_cn(int state, int max_haplotype_change){
    int cn = 0;
    vector<int> sums;
    // cout << "Sums of state: ";
    int cn_max = 4;
    vector<int> n_states = make_peak_vector(3);   // the number of states corresponding to each copy number: 1,2,3,4,3,2,1 when +2 allowed
    if(max_haplotype_change == 2){
        cn_max = 6;
        n_states = make_peak_vector(4); 
    }

    int s = 0;
    for(int i = 0; i <= cn_max; i++){
        s = s + n_states[i];
        sums.push_back(s);
        // cout << "\t" << s;
    }
    // cout << endl;

    if(state < sums[0]) return 0;

    int i = 1;
    do{
        if(state >= sums[i - 1] && state < sums[i]) return i;
        i++;
    }while(i <= cn_max);

    return 0;
}


// only cn_max 4, 6 supported for now, can be extended to higher copy number if needed
void change_state_to_allele_cn(int state, int max_haplotype_change, int& cnA, int& cnB){
    int cn = 0;
    vector<int> sums;
    // cout << "Sums of state: ";
    int cn_max = 4;
    vector<int> n_states = make_peak_vector(3);   // the number of states corresponding to each copy number: 1,2,3,4,3,2,1 when +2 allowed
    if(max_haplotype_change == 2){
        cn_max = 6;
        n_states = make_peak_vector(4); 
    }

    int s = 0;
    for(int i = 0; i <= cn_max; i++){
        s = s + n_states[i];
        sums.push_back(s);
        // cout << "\t" << s;
    }
    // cout << endl;

    if(state < sums[0]) return;

    int i = 1;
    do{
        // order of first element changed due to diagonal element grouping (not strictly increasing from 0)
        // 2: 0, 2, 1; 3: 1, 2(, 0); 4: 2(, 1, 3); 5: 2, 3; 6: 3
        if(state >= sums[i - 1] && state < sums[i]){
            int diff = state - sums[i-1];
            if(i <= 1){              
                cnA = diff;
                cnB = i - cnA;
            }else if(i == 2){
                int shift = 0;
                switch (diff){
                case 1:
                    shift = 1;
                    break;
                case 2:
                    shift = -1;
                    break;
                default:
                    break;
                }
                cnA = diff + shift;
                cnB = i - cnA;
            }else if(i == 3){
                if(cn_max > 4 && diff == 2){
                    cnA = 0;
                    cnB = 3;
                }else{
                    cnA = diff + 1;
                    cnB = i - cnA;
                }
            }else if(i == 4){
                if(cn_max == 4 || diff == 0){
                    cnA = 2;
                    cnB = 2;
                }else{ 
                    assert(diff > 0);
                    int shift = 0;
                    if(diff == 2){
                        shift = 1;
                    }
                    cnA = diff + shift;
                    cnB = i - cnA;
                }
            }
            else if(i == 5){
                cnA = diff + 2;
                cnB = i - cnA;
            }else if(i == 6){
                cnA = 3;
                cnB = 3;
            }else{ }
        }
        i++;
    }while(i <= cn_max);
}


/** 
 * @brief Decompose relative total copy number into multi-level changes (WGD, chromosome, segment) given original data and computed changes.
 * @param s_info Sample information where each innermost vector contains copy number for a segment in the format [chr, sid, cn] for each sample.
 * @param cn_change_info A vector of CN_CHANGE structs to store the decomposed copy number changes for each sample, same dimensions as s_info.
 * @param input_data INPUT_DATA struct containing sample-level information needed for decomposition.
 * @param cn_max Maximum copy number allowed by the program.
 * @param is_total Indicator whether the copy numbers are total (1) or haplotype-specific (0). Only for total copy number decomposition for now.
 * @param debug Debug flag for verbose output.
 */
void rcn_to_decomposition(const vector<vector<vector<int>>>& s_info, vector<vector<CN_CHANGE>>& s_info_change, vector<vector<int>>& sample_change_site, vector<int>& site_max_change, vector<int>& sample_max_cn, int debug){
    if(debug) cout << "\tChanging relativecopy number to multiple level changes" << endl;

    if(debug > 1){
        for(size_t i = 0; i < s_info.size(); i++){
            cout << "\nSample " << (i+1) << " original copy number:" << endl;
            for(size_t j = 0; j < s_info[i].size(); j++){
                cout << "\tSegment " << (j+1) << ": ";
                for(size_t k = 0; k < s_info[i][j].size(); k++){
                    cout << "\t" << s_info[i][j][k];
                }
                cout << endl;
            }
        }
    }

    s_info_change.resize(s_info.size());    // for samples
    for (size_t i = 0; i < s_info.size(); i++) {     // for all segments
        s_info_change[i].resize(s_info[i].size());
    }
    site_max_change.clear();
    site_max_change.reserve(s_info.size());

    for(size_t i = 0; i < s_info.size(); i++){
        vector<vector<int>> s_cn = s_info[i];
        vector<int> cns;    // store all copy numbers to get the maximum value
        vector<int> site_change;     // copy number changes across all chromosomes in sample i       
        int max_dev = 0;        // max deviation from avg_cn
        int max_chr = -1, max_seg = -1, max_cn_val = -1;

        // get multi-level CN changes for each sample    
        for(size_t j = 0; j < s_cn.size(); j++){
            int chr = s_cn[j][0];
            int cn  = s_cn[j][2];

            CN_CHANGE cc;
            cc.cn_state = cn;
            // cc.cnA = -1;
            // cc.cnB = -1;
            cc.num_wgd = 0;
            cc.cn_change_chr = 0;
            int cn_change = cn - NORM_PLOIDY;
            cc.cn_change_site = cn_change;
            s_info_change[i][j] = cc;

            cns.push_back(cn);
            site_change.push_back(cn_change);

            int dev = abs(cn_change);
            if(dev > max_dev){
                max_dev = dev;
                max_chr = s_cn[j][0];
                max_seg = s_cn[j][1];
                max_cn_val = cn;
            }           
        }

        sample_change_site.push_back(site_change);
        int mcn = *max_element(cns.begin(), cns.end());     // Getting the largest copy number for each sample
        sample_max_cn.push_back(mcn);
        site_max_change.push_back(max_dev);    
              
        if(debug){
            cout << "Sample " << (i+1)
                 << ": max_cn=" << mcn
                 << " |max_cn-avg_cn|=" << max_dev
                 << " -> site_change=" << max_dev
                 << " (at chr " << max_chr << ", seg " << max_seg << ", cn " << max_cn_val << ")"
                 << endl;
        }        
    }

    if(debug > 1){
        for(size_t i = 0; i < s_info_change.size(); i++){
            cout << "\nSample " << (i+1) << " decomposed copy number changes:" << endl;
            for(size_t j = 0; j < s_info_change[i].size(); j++){
                cout << "\tSegment " << (j+1) << ": " << s_info_change[i][j] << endl;
            }
        }
    }
}



/** 
 * @brief Decompose total copy number into multi-level changes (WGD, chromosome, segment) given original data and computed changes.
 * @param s_info Sample information where each innermost vector contains copy number for a segment in the format [chr, sid, cn] for each sample.
 * @param cn_change_info A vector of CN_CHANGE structs to store the decomposed copy number changes for each sample, same dimensions as s_info.
 * @param input_data INPUT_DATA struct containing sample-level information needed for decomposition.
 * @param cn_max Maximum copy number allowed by the program.
 * @param is_total Indicator whether the copy numbers are total (1) or haplotype-specific (0). Only for total copy number decomposition for now.
 * @param debug Debug flag for verbose output.
 */
void cn_to_decomposition(const vector<vector<vector<int>>>& s_info, vector<vector<CN_CHANGE>>& s_info_change, const INPUT_DATA& input_data, int debug){
    if(debug) cout << "\tChanging copy number to multiple level changes" << endl;
    vector<int> sample_num_wgd = input_data.sample_num_wgd;
    vector<vector<int>> sample_change_chr = input_data.sample_change_chr;   
    vector<vector<int>> sample_change_site = input_data.sample_change_site;   

    if(debug > 1){
        for(size_t i = 0; i < s_info.size(); i++){
            cout << "\nSample " << (i+1) << " original copy number:" << endl;
            for(size_t j = 0; j < s_info[i].size(); j++){
                cout << "\tSegment " << (j+1) << ": ";
                for(size_t k = 0; k < s_info[i][j].size(); k++){
                    cout << "\t" << s_info[i][j][k];
                }
                cout << endl;
            }
        }
    }

    s_info_change.resize(s_info.size());    // for samples
    for (size_t i = 0; i < s_info.size(); i++) {     // for all segments
        s_info_change[i].resize(s_info[i].size());
    }

    for(size_t i = 0; i < s_info.size(); i++){
        vector<vector<int>> s_cn = s_info[i];   // vector<int>: each segment in the sample, with format [chr, sid, cn]

        // get multi-level CN changes for each sample    
        for(size_t j = 0; j < s_cn.size(); j++){
            int chr = s_cn[j][0];
            int cn  = s_cn[j][2];

            CN_CHANGE cc;
            cc.cn_state = cn;
            // cc.cnA = -1;
            // cc.cnB = -1;
            cc.num_wgd = sample_num_wgd[i];
            cc.cn_change_chr = sample_change_chr[i][chr - 1];   // chr is 1-based in input, but sample_change_chr is 0-based
            int site_change = sample_change_site[i][j];
            if(cc.cn_change_chr != 0){   // if there is chromosome change, adjust site change accordingly to avoid double counting
                if(cc.cn_change_chr % 1000 != 0){
                    site_change = site_change - cc.cn_change_chr;
                }else{
                    site_change = site_change - cc.cn_change_chr / 1000;
                }
            }
            cc.cn_change_site = site_change;

            s_info_change[i][j] = cc;
        }
    }

    if(debug){
        for(size_t i = 0; i < s_info_change.size(); i++){
            cout << "\nSample " << (i+1) << " decomposed copy number changes:" << endl;
            for(size_t j = 0; j < s_info_change[i].size(); j++){
                cout << "\tSegment " << (j+1) << ": " << s_info_change[i][j] << endl;
            }
        }
    }
}


// segment level using absolute copy numbers with larger level changes excluded
/** 
 * @brief Decompose haplotype-specific copy number into multi-level changes (WGD, chromosome, segment) given original data and computed changes.
 * @param s_info Sample information where each innermost vector contains copy number for a segment in the format [chr, sid, cn] for each sample.
 * @param cn_change_info A vector of CN_CHANGE structs to store the decomposed copy number changes for each sample, same dimensions as s_info.
 * @param input_data INPUT_DATA struct containing sample-level information needed for decomposition.
 * @param cn_max Maximum copy number allowed by the program.
 * @param is_total Indicator whether the copy numbers are total (1) or haplotype-specific (0). Only for total copy number decomposition for now.
 * @param debug Debug flag for verbose output.
 */
void hcn_to_decomposition(const vector<vector<vector<int>>>& s_info,  vector<vector<CN_CHANGE>>& s_info_change, const INPUT_DATA& input_data, int debug){
    if(debug) cout << "\tChanging copy number to multiple level changes" << endl;
    vector<int> sample_num_wgd = input_data.sample_num_wgd;
    vector<vector<int>> sample_change_chr = input_data.sample_change_chr;   
    vector<vector<int>> sample_change_site = input_data.sample_change_site;   

    if(debug > 1){
        for(size_t i = 0; i < s_info.size(); i++){
            cout << "\nSample " << (i+1) << " original copy number:" << endl;
            for(size_t j = 0; j < s_info[i].size(); j++){
                cout << "\tSegment " << (j+1) << ": ";
                for(size_t k = 0; k < s_info[i][j].size(); k++){
                    cout << "\t" << s_info[i][j][k];
                }
                cout << endl;
            }
        }
    }

    s_info_change.resize(s_info.size());    // for samples
    for (size_t i = 0; i < s_info.size(); i++) {     // for all segments
        s_info_change[i].resize(s_info[i].size());
    }

    for(size_t i = 0; i < s_info.size(); i++){
        vector<vector<int>> s_cn = s_info[i];

        // get multi-level CN changes for each sample    
        for(size_t j = 0; j < s_cn.size(); j++){
            int chr = s_cn[j][0];
            int cn  = s_cn[j][2];

            CN_CHANGE cc;
            cc.cn_state = cn;
            // cc.cnA = -1;
            // cc.cnB = -1;
            cc.num_wgd = sample_num_wgd[i];
            cc.cn_change_chr = sample_change_chr[i][chr];
            // state for copy number values 
            cc.cn_change_site = sample_change_site[i][j];

            s_info_change[i][j] = cc;
        }
    }

    if(debug > 1){
        for(size_t i = 0; i < s_info_change.size(); i++){
            cout << "\nSample " << (i+1) << " decomposed copy number changes:" << endl;
            for(size_t j = 0; j < s_info_change[i].size(); j++){
                cout << "\tSegment " << (j+1) << ": " << s_info_change[i][j] << endl;
            }
        }
    }    
}


/********************* Reading original input files  ***************************/
/** 
 * @brief Read time information from a file
 * @param filename The name of the input file containing time data.
 * @param Ns The number of samples to read from the file.   
 * @param age A reference to an integer that will store the minimum age read from the file.
 * @param debug An integer flag for enabling debug output.  
 */
vector<double> read_time_info(const string& filename, const int Ns, int& age, int debug){
  if(debug) cout << "\tread_time_info" << endl;

  vector<double> t_info;
  vector<int> ages;
  ifstream infile(filename.c_str());
  if(infile.is_open()){
    string line;
    while(!getline(infile, line).eof()){
      if(line.empty()) continue;

      vector<string> split;
      string buf;
      stringstream ss(line);
      while(ss >> buf) split.push_back(buf);

      double dt = atof(split[1].c_str());
      //cout << "read dt: " << dt << endl;
      t_info.push_back(dt);

      if(split.size() > 2){
          int a = atoi(split[2].c_str());
          ages.push_back(a);
      }
    }
    if(ages.size() > 0){
        age = *min_element(ages.begin(), ages.end());
    }
  }else{
    cerr << "Error: open of time data unsuccessful: " << filename << endl;
    exit(EXIT_FAILURE);
  }

  if(t_info.size() != Ns){
    cerr << "Error: timing information does not contain " << Ns << " entries: " << filename << endl;
    exit(EXIT_FAILURE);
  }

  return t_info;
}



/** 
 *  @brief Read the input copy number file
 *  @param s_info A reference to a vector of vectors of vectors that will store the copy number information for each sample in the format s_info[sample - 1] = {vcn = {chr, sid, cn}}.
 *  @param filename The name of the input file containing copy number data.
 *  @param Ns The number of samples to read from the file.
 *  @param num_total_bins A reference to an integer that will store the total number of bins read from the file.
 *  @param cn_max The maximum allowed copy number.
 *  @param is_total An integer flag indicating whether the copy numbers are total (1) or haplotype-specific (0).
 *  @param is_rcn An integer flag indicating whether the copy numbers are relative copy numbers (1) or absolute copy numbers (0).
 *  @param debug An integer flag for enabling debug output.
 */ 
void read_cn(vector<vector<vector<int>>>& s_info, const string& filename, int Ns, int& num_total_bins, int cn_max, int is_total, int is_rcn, int debug){
    num_total_bins = 0;
    // data indexed by [sample][data][ chr, bid, cn ]
    for(int i = 0; i < Ns; ++i) s_info.push_back(vector<vector<int>>());

    igzstream infile(filename.c_str());
    int counter = 0;
    string line;
    int prev_sample = 1;

    while(!getline(infile, line).eof()){
      if(line.empty()) continue;

      vector<string> split;
      string buf;
      stringstream ss(line);
      while(ss >> buf) split.push_back(buf);

      const char* sstr = split[0].c_str();
      if(!isdigit(*sstr)){
          cout << "Sample ID must be an integer, ordered from 1 to the specified number of patient samples!" << endl;
          exit(EXIT_FAILURE);
      }

      int sample = atoi(sstr);
      // Read next sample
      if(sample > Ns){
          cout << "Skipping sample with ID larger than " << Ns << endl;
          break;
      }
      if(prev_sample != sample){
          num_total_bins = counter;
          counter = 0;
      }

      int chr = atoi(split[1].c_str());  // chr
      int sid = atoi(split[2].c_str());  // site ID
      int cn = -1;

      if(is_rcn){
          cn = atoi(split[3].c_str());  // copy number
          cn = cn + NORM_PLOIDY;
      }else{
        if(is_total){
            if(split.size() != 4){
                cout << "There should be 4 columns in the file for total copy numbers" << endl;
                exit(EXIT_FAILURE);
            }
            cn = atoi(split[3].c_str());  // copy number
            if(cn < 0){
                cout << "Negative copy numbers in line " << line << "!" << endl;
                exit(EXIT_FAILURE);
            }
            if(cn > cn_max){
                if(debug) cout << "copy number " << cn << " is decreased to " << cn_max << endl;
                cn = cn_max;
            }
        }else{
            if(split.size() != 5){
                cout << "Current file has " << split.size() << " columns " << endl;
                cout << "There should be 5 columns in the file for haplotype-specific copy numbers" << endl;
                exit(EXIT_FAILURE);
            }
            int cn1 = atoi(split[3].c_str());  // copy number
            int cn2 = atoi(split[4].c_str());  // copy number
            if(cn1 < 0 || cn2 < 0){
                cout << "Negative copy numbers in line " << line << "!" << endl;
                exit(EXIT_FAILURE);
            }
            int tcn = cn1 + cn2;
            if(tcn > cn_max){   // a bit hard to decease haplotype-specific CNs
                cout << "Total copy number " << tcn << " is larger than " << cn_max << "! Please adjust input accordingly!";
              //   int larger_cn = (cn1 > cn2) ? cn1 : cn2;
              //   int diff = tcn - cn_max;
                exit(EXIT_FAILURE);
            }
            // It may include states not considered by the site change matrix due to constraints on the allowed number of changes at site level alone
            cn = allele_cn_to_state(cn1, cn2);  // convert haplotype-specific copy numbers into state directly while reading
        }
      }

      vector<int> vcn{chr, sid, cn};
      s_info[sample - 1].push_back(vcn);

      counter++;

      prev_sample = sample;
    }

    if(debug) cout << "\tSuccessfully read input file with " << num_total_bins << " sites" << endl;

}



/************************* Analysing copy number changes ******************/
/** 
 *  @brief Calculate average copy number or ploidy for a single sample.
 *  @param s_cn copy number information for a sample in the format [chr, sid, cn].
 *  @param chr_cn A map to store copy numbers for each chromosome.  
 *  @param cn_max Maximum copy number.
 *  @param is_total Indicator whether the copy numbers are total (1) or haplotype-specific (0).
 *  @return The average copy number or ploidy for the sample.
 */
double compute_sample_avg_cn(const vector<vector<int>>& s_cn, map<int, vector<int>>& chr_cn, int cn_max, int is_total){
    if (s_cn.empty()) return 0.0;

    double sum_cn = 0.0;
    int num_seg = 0;

    // iterate through all records for a sample
    for(size_t j = 0; j < s_cn.size(); j++){
        int chr = s_cn[j][0];
        int cn = s_cn[j][2];
        if(!is_total) cn = state_to_total_cn(cn, cn_max);   
        chr_cn[chr].push_back(cn);
        sum_cn += cn;
        num_seg++;
    }

    double avg_cn = sum_cn / (double)num_seg;

    return avg_cn;
}


void get_sample_ploidy(const vector<vector<vector<int>>>& s_info, vector<double>& sample_avg_cn, vector<map<int, vector<int>>>& sample_chr_cn, int cn_max, int is_total, int debug){
    cout << "\nGetting the average copy number (ploidy) in each sample" << endl;

    for(size_t i = 0; i < s_info.size(); i++){
        vector<vector<int>> s_cn = s_info[i];
        map<int, vector<int>> chr_cn;
        double avg_cn = compute_sample_avg_cn(s_cn, chr_cn, cn_max, is_total);
        sample_avg_cn.push_back(avg_cn);
        sample_chr_cn.push_back(chr_cn);

        cout << "Sample " << i + 1 << " has average copy number " << avg_cn << endl;
    }
}

/** 
 *  @brief Estimate the potential number of WGDs for each sample        
 *  @param s_info A vector of vectors of vectors where each innermost vector contains copy number information for a segment in the format [chr, sid, cn] for each sample.
 *  @param sample_num_wgd A vector to store the estimated number of WGDs for each sample in the input.
 *  @param debug Debug flag for verbose output.
 */
void get_num_wgd(const vector<vector<vector<int>>>& s_info, const vector<double>& sample_avg_cn, vector<int>& sample_num_wgd, int debug){
    cout << "\nGetting the potential number of WGDs for each sample" << endl;
   
    for(size_t i = 0; i < s_info.size(); i++){
        int nwgd = 0;
        // int sum_cn = 0;
        // vector<vector<int>> s_cn = s_info[i];
        // count the presence of each copy number
        // map<int, int> cn_count;
        // for(int k = 0; k < cn_max; k++){
        //     cn_count[k] = 0;
        // }
        // for (int j = 0; j < s_cn.size(); j++){
        //     vector<int> cp = s_cn[j];
        //     int cn = cp[2];
        //     if(!is_total) cn = state_to_total_cn(cn, cn_max);
        //     // cn_count[cn] += 1;
        //     sum_cn += cn;
        // }
        // // Find the most frequent or mode copy number
        // int most_freq_cn = 0;
        // int max_count = 0;
        // for(auto cnc : cn_count){
        //     if(debug) cout << cnc.first << "\t" << cnc.second << endl;
        //     if(cnc.second > max_count){
        //         max_count = cnc.second;
        //         most_freq_cn = cnc.first;
        //     }
        // }
        //
        // int mode_logcn = ceil(log2(most_freq_cn));
        // if(mode_logcn > 1) nwgd = mode_logcn - 1;

        // use ploidy, to be applicable to real data, assume at most 1 WGD events
        // Estimate number of WGD based on ploidy
        // float ploidy = (float) sum_cn / s_cn.size();
        double ploidy = sample_avg_cn[i];
        if(ploidy > 2 * WGD_CUTOFF){
            nwgd = 2;
        }else if(ploidy > WGD_CUTOFF){
            nwgd = 1;
        }else{
            nwgd = 0;
        }

        sample_num_wgd.push_back(nwgd);
        // if(debug)
        cout << "Sample " << i + 1 << " probably has " << nwgd << " WGD event, with ploidy " << ploidy << endl;
    }
}



/** 
 *  @brief Estimate the potential number of chromosome changes for each sample        
 *  @param s_info Sample information where each innermost vector contains copy number for a segment in the format [chr, sid, cn] for each sample.
 *  @param sample_avg_cn A vector to store the average copy number for each sample in the input.
 *  @param sample_chr_cn A vector of maps to store copy numbers for each chromosome in each sample.
 *  @param sample_change_chr A vector to store the estimated number of chromosome changes for each sample in the input, indexed by [sample][chromosome].
 *  @param chr_max_change A vector to store the maximum absolute chromosome change for each sample in the input.
 *  @param cn_max Maximum copy number.
 *  @param is_total Indicator whether the copy numbers are total (1) or haplotype-specific (0).
 *  @param debug Debug flag for verbose output.
 */
void get_chr_change(const vector<vector<vector<int>>>& s_info, const vector<double>& sample_avg_cn, const vector<map<int, vector<int>>>& sample_chr_cn, vector<vector<int>>& sample_change_chr, vector<int>& chr_max_change, int cn_max, int is_total, int debug){
    cout << "\nGetting the potential number of chromosome changes for each sample" << endl;
    
    // vector<int> chr_gain_max, chr_loss_min;
    for(size_t i = 0; i < s_info.size(); i++){
        // chr, seg, CN
        vector<vector<int>> s_cn = s_info[i];
        map<int, vector<int>> chr_cn = sample_chr_cn[i];
        double avg_cn = sample_avg_cn[i];
        vector<int> chr_change;     // copy number changes across all chromosomes in sample i
          
        // avg_cn = compute_sample_avg_cn(s_cn, chr_cn, cn_max, is_total); //replace above block

        for(auto c : chr_cn){
            vector<int> cp = c.second;
            int chr_sum_cn = accumulate(cp.begin(), cp.end(), 0);
            double avg_chr_cn = (double) chr_sum_cn / cp.size();
            
            // If many chr gains, avg_cn will be large, and the difference will be small
            // this number can be very large if avg_cn is large caused by WGD
            double num_change = avg_chr_cn - avg_cn;
            // int round_num_change = (int) (num_change + 0.5 - (num_change < 0));
            int round_num_change = (int) lround(num_change);
            // use multiple of 2 to roughly determine events before or after WGD
            if(avg_cn > WGD_CUTOFF){    // one WGD 
                int reminder = round_num_change % NORM_PLOIDY;
                // likely occurred before WGD, with copy number change in multiple of 2 or ploidy
                // Another possibility is that the change occurred after WGD but with very large copy number change
                // encode ambiguity by allowing both possibilities 
                if(reminder == 0){   
                    round_num_change = round_num_change / NORM_PLOIDY;  // assume before WGD
                    round_num_change = round_num_change * 1000;  // encode ambiguous orders regarding WGD
                }
            }
            // if(round_num_change < -2) round_num_change = -2;
            // // TODO: change to dynamic threshold later
            // if(round_num_change > 2) round_num_change = 2;
            cout << "Number of segments in chromosome " << c.first << " is " << cp.size() << "; avg cn: " << avg_chr_cn << "; exact num changes: " << num_change  << "; num changes: " << round_num_change << endl;

            chr_change.push_back(round_num_change);
        }

        if(debug){
            cout << "Sample " << i+1 << endl;
            // cout << "Number of segments " << num_seg << endl;
            int num_gain = count_if(chr_change.begin(), chr_change.end(), [](int c){return c >= 1;});
            int num_loss = count_if(chr_change.begin(), chr_change.end(), [](int c){return c <= -1;});
            // cout << "There are probably " << num_gain << " chromosome gain events and " << num_loss << " chromosome loss events" << endl;
            cout << "Chromosomes with gain-like shift: " << num_gain
                 << ", loss-like shift: " << num_loss << endl;
            cout << "Average copy number in the genome " << avg_cn << endl;
        }

        // calculate max gain/loss and absolute change
        int gain_cn = 0;
        int loss_cn = 0;
        int max_abs_change = 0;
        if (!chr_change.empty()) {
            gain_cn = max(0, *max_element(chr_change.begin(), chr_change.end()));
            loss_cn = min(0, *min_element(chr_change.begin(), chr_change.end()));
            max_abs_change = max(abs(gain_cn), abs(loss_cn));
        }
        cout << "   Maximum gain on a chromosome in sample " << i+1 << " is " << gain_cn << endl;
        cout << "   Maximum loss on a chromosome in sample " << i+1 << " is " << loss_cn << endl;
        cout << "   Maximum absolute change on a chromosome in sample " << i+1 << " is " << max_abs_change << endl;

        // chr_gain_max.push_back(gain_cn);
        // chr_loss_min.push_back(loss_cn);
        chr_max_change.push_back(max_abs_change);
        sample_change_chr.push_back(chr_change);
    }
}


// get maximum copy number for each sample, not used for now, integrated into get_site_change
void get_sample_mcn(const vector<vector<vector<int>>>& s_info, vector<int>& sample_max_cn, int cn_max, int is_total, int debug){
    for(size_t i = 0; i < s_info.size(); i++){
        vector<vector<int>> s_cn = s_info[i];
        vector<int> cns;
        // iterate through all records for a sample
        for(int j = 0; j < s_cn.size(); j++){
            int cn = s_cn[j][2];
            // if(debug) cout << s_cn[j][0] << "\t" << s_cn[j][1] << "\t" << s_cn[j][2] << "\n";
            if(!is_total) cn = state_to_total_cn(cn, cn_max);
            cns.push_back(cn);
        }
        int mcn = *max_element(cns.begin(), cns.end());
        sample_max_cn.push_back(mcn);
        if(debug) cout << "Largest copy number in sample " << i + 1 << " is " << mcn << endl;
    }
}


/** 
 *  @brief Compute per-sample site-level max deviation from avg_cn (baseline), then take ceiling as an integer bound.
 *  @param s_info Sample information where each innermost vector contains copy number for a segment in the format [chr, sid, cn] for each sample.
 *  @param site_max_change A vector to store the estimated number of site changes for each sample in the input.
 *  @param sample_max_cn A vector to store the maximum copy number for each sample in the input.
 *  @param cn_max Maximum copy number.
 *  @param is_total Indicator whether the copy numbers are total (1) or haplotype-specific (0).
 *  @param debug Debug flag for verbose output.
 */         
void get_site_change(const vector<vector<vector<int>>>& s_info, const vector<double>& sample_avg_cn, vector<vector<int>>& sample_change_site, vector<int>& site_max_change, vector<int>& sample_max_cn, int cn_max, int is_total, int debug){
    cout << "\nGetting the potential number of site changes for each sample" << endl;
  
    site_max_change.clear();
    site_max_change.reserve(s_info.size());

    for(size_t i = 0; i < s_info.size(); i++){
        vector<vector<int>> s_cn = s_info[i];
        if(s_cn.empty()){
            cout << "Sample " << (i+1) << " has 0 records, site_change=0" << endl;
            exit(EXIT_FAILURE);
        }

        vector<int> cns;    // store all copy numbers to get the maximum value
        double avg_cn = sample_avg_cn[i];
        vector<int> site_change;     // copy number changes across all chromosomes in sample i       
        int max_dev = 0;        // max deviation from avg_cn
        int max_chr = -1, max_seg = -1, max_cn_val = -1;

        for(size_t j = 0; j < s_cn.size(); j++){
            int cn = s_cn[j][2];
            if(!is_total) cn = state_to_total_cn(cn, cn_max);
            cns.push_back(cn);

            int cn_change = cn - lround(avg_cn);
            // use multiple of 2 to roughly determine events before or after WGD
            if(avg_cn > WGD_CUTOFF){    // one WGD 
                int reminder = cn % NORM_PLOIDY;
                // TODO: set ambuiguity flag
                if(reminder == 0){   // likely occurred before WGD, with copy number change in multiple of 2 or ploidy
                    cn_change = cn / NORM_PLOIDY - NORM_PLOIDY;
                }
            }
            
            // TODO: set those outside limits of rate matrix
            if(cn_change < -2) cn_change = -2;
            if(cn_change > 4) cn_change = 4;       

            // used to set change_site for CN_CHANGE variable for this site on this chr in the sample 
            site_change.push_back(cn_change);

            int dev = abs(cn_change);
            if(dev > max_dev){
                max_dev = dev;
                max_chr = s_cn[j][0];
                max_seg = s_cn[j][1];
                max_cn_val = cn;
            }
        }   // finish traversing all sites in sample i
      
        sample_change_site.push_back(site_change);
        int mcn = *max_element(cns.begin(), cns.end());     // Getting the largest copy number for each sample
        sample_max_cn.push_back(mcn);
        site_max_change.push_back(max_dev);
        
        if(debug){
            cout << "Sample " << (i+1)
                 << ": max_cn=" << mcn
                 << " avg_cn=" << avg_cn
                 << " |max_cn-avg_cn|=" << max_dev
                 << " -> site_change=" << max_dev
                 << " (at chr " << max_chr << ", seg " << max_seg << ", cn " << max_cn_val << ")"
                 << endl;
        }
    }
}


/************************* Parsing segments ******************/
/** 
 * @brief Identify variable bins across all samples for segment merging.
 * @param s_info Sample information where each innermost vector contains copy number for a segment in the format [chr, sid, cn] for each sample.
 * @param Ns The number of samples.
 * @param num_total_bins The total number of bins.      
 * @param num_invar_bins A reference to an integer that will store the number of invariable bins found.
 * @param var_bins A reference to a vector that will store flags indicating whether each bin is variable (1) or invariable (0).
 * @param is_total An integer flag indicating whether the copy numbers are total (1) or haplotype-specific (0).
 * @param debug Debug flag for verbose output.
 */
void get_var_bins(const vector<vector<vector<int>>>& s_info, int Ns, int num_total_bins, int& num_invar_bins, vector<int>& var_bins, int is_total, int debug){
    for(int k = 0; k < num_total_bins; ++k){
        // using sum of CNs across samples to detect variant bins does not work for special cases
        vector<int> cns;
        for(int i = 0; i < Ns; ++i){
            int cn = s_info[i][k][2];
            cns.push_back(cn);
        }

        bool is_invar = true;

        if(is_total){
            is_invar = all_of(cns.begin(), cns.end(), [&] (int i) {return i == NORM_PLOIDY;});
        }else{ // For haplotype-specific CN, normal state is 1/1, with ID 4
            is_invar = all_of(cns.begin(), cns.end(), [&] (int i) {return i == NORM_ALLElE_STATE;});
        }

        if(is_invar){
          num_invar_bins += 1;
        }else{
          var_bins[k] = 1;
        }
    }

    if(debug > 1){
        cout << "\n\tVariable bins found:" << endl;
        for(int k = 0; k < num_total_bins; ++k){
            if(var_bins[k]){
              cout << s_info[0][k][0] << "\t" << s_info[0][k][1];
              for(int i = 0; i < Ns; ++i) cout << "\t" << s_info[i][k][2];
              cout << endl;
            }
        }
    }

    int nvar = accumulate(var_bins.begin(), var_bins.end(), 0);
    cout << "\nSummary of variable/invariable bins:" << endl;
    cout << "\tTotal number of bins:\t" << num_total_bins << endl;
    cout << "\tNumber of variable bins:\t" << nvar << endl;
    cout << "\tNumber of invariable bins:\t" << num_invar_bins << endl;
}



/** 
 * @brief Merge consecutive variable bins with the same copy number profiles across all samples into segments.
 * @param s_info Sample information where each innermost vector contains copy number for a segment in the format [chr, sid, cn] for each sample.
 * @param Ns The number of samples.
 * @param num_total_bins The total number of bins.
 * @param num_invar_bins A reference to an integer that will store the number of invariable bins found.
 * @param is_total An integer flag indicating whether the copy numbers are total (1) or haplotype-specific (0).
 * @param debug Debug flag for verbose output.
 * @return A vector of vectors representing the merged segments, where each seg is represented by: {chr, id_start, id_end, seg_start, seg_end}
 */
vector<vector<int>> get_invar_segs(const vector<vector<vector<int>>>& s_info, int Ns, int num_total_bins, int& num_invar_bins, int is_total, int debug){
    num_invar_bins = 0;
    vector<int> var_bins(num_total_bins, 0);
    if(debug) cout << "\tGetting all the variable bins" << endl;
    get_var_bins(s_info, Ns, num_total_bins, num_invar_bins, var_bins, is_total, debug);

    vector<vector<int>> segs;
    for(int k = 0; k < num_total_bins;){  
        if(var_bins[k]){ // only starting from variable bins
          int chr = s_info[0][k][0];
          int seg_start = s_info[0][k][1];
          int id_start = k;

          // hold all the sites in a bin
          vector<int> prev_bin;
          for(int j = 0; j < Ns; ++j) prev_bin.push_back(s_info[j][k][2]);
        //   if(debug) cout << "seg_start: " << chr << "\t" << seg_start << ", cn =  " << s_info[0][k][2] << endl;

          // Check the subsequent bins
          bool const_cn = true;
          k++;
          while(k < num_total_bins && var_bins[k] && const_cn){ // break if next bin is invarible or has different CN
              vector<int> curr_bin;
              for(int j = 0; j < Ns; ++j) curr_bin.push_back(s_info[j][k][2]);
              if(is_equal_vector(prev_bin, curr_bin) && s_info[0][k][0] == chr){
            	  const_cn = true;
            	  ++k;
              }else{
            	  const_cn = false;
            	//   if(debug) cout << "\tsplitting segment" << endl;
              }
          }
          int seg_end = s_info[0][k-1][1];
          int id_end = k - 1;
          // if(debug) cout << "seg_end:\t" << seg_end << "\t" << k << endl;

          // id_*: start from 0 until #segments; seg_*: original segment ID
          vector<int> seg{chr, id_start, id_end, seg_start, seg_end};
          segs.push_back(seg);

          // rewind k by one to get the split segment start correct
          if(!const_cn) k--;
        }
        ++k;
    }
    cout << "\tFound segments (after merging consecutive bins):\t" << segs.size() << endl;

    return segs;
}


/** 
 * @brief Compute average copy number for a given segment across all samples at site i
 * @param i Index of the segment across all the sites
 * @param segs Vector of segments represented by: {chr, id_start, id_end, seg_start, seg_end}, indicating segment locations
 * @param s_info Sample information where each innermost vector contains copy number for a segment in the format [chr, sid, cn] for each sample
 * @param Ns Number of samples.
 * @param cn_max Maximum copy number state, different from maximum copy number value in the input when the input is haplotype-specific copy number.
 * @param debug Debug flag for verbose output
 * @return A vector of average copy numbers for the segment across all samples
 */
void compute_segment_cn_state(vector<int>& seg_cn, const vector<int>& segs_curr, const vector<vector<vector<int>>>& s_info, int Ns, int cn_max, int debug){
    if(debug > 1){
        cout << "\nComputing segment copy number for segment on chr " << segs_curr[0] << ":" << segs_curr[1] + 1 << "\t" << segs_curr[2] + 1 << "\t" << segs_curr[3] << "\t" << segs_curr[4] << endl;
    }
    assert(s_info.size() == Ns);

    for(int j = 0; j < Ns; ++j){
        int sum_cn = 0;
        float avg_seg_cn = 0.0;
        int bin_cn = s_info[j][segs_curr[1]][2];  // copy number of the/total first bin in the segment

        // sum copy numbers across all bins in the segment to compute average, actually all bins should have the same CN
        for(int k = segs_curr[1]; k < (segs_curr[2] + 1); ++k){
            sum_cn += s_info[j][k][2];
        }
        avg_seg_cn = sum_cn / (segs_curr[2] - segs_curr[1] + 1);

        // The average should be the same as the value of each bin
        assert(avg_seg_cn == bin_cn);

        // check all CNs across the segment are integers
        if(ceil(avg_seg_cn) != floor(avg_seg_cn)){
            cout << "Fractional copy number << " << avg_seg_cn << " at " << segs_curr[1] << ", " << segs_curr[2] << endl;
            exit(EXIT_FAILURE);
        }

        if(avg_seg_cn > cn_max){
            cout << "INVALID (larger than maximum allowed) copy number << " << avg_seg_cn << " at " << segs_curr[1] << ", " << segs_curr[2] << endl;
            exit(EXIT_FAILURE);
        }

        if(debug > 1){
            cout << "\tSample " << j + 1 << " with CN : " << bin_cn << endl;

        }        
        seg_cn.push_back(bin_cn);
    }

    if(debug > 1 ){
        cout << "\nFinal segment copy number for this segment: ";
        for(int j = 0; j < Ns; ++j) cout << "\t" << seg_cn[j];
        cout << endl;
        cout << "----------------------------------------" << endl;
    }

}


/** 
 * @brief Compute average copy number for a given segment across all samples at site i
 * @param seg_cn  A vector of average copy numbers for the segment across all samples
 * @param segs_curr Vector of segments represented by: {chr, id_start, id_end, seg_start, seg_end}, indicating segment locations
 * @param s_info Sample information where each innermost vector contains copy number for a segment in the format [chr, sid, cn] for each sample
 * @param Ns Number of samples.
 * @param debug Debug flag for verbose output
 */
void compute_segment_cn_change(vector<CN_CHANGE>& seg_cn, const vector<int>& segs_curr, const vector<vector<CN_CHANGE>>& s_info_change, int Ns, int Nsite, int debug){
    if(debug > 1){
        cout << "\nComputing segment copy number for segment on chr " << segs_curr[0] << ":" << segs_curr[1] + 1 << "\t" << segs_curr[2] + 1 << "\t" << segs_curr[3] << "\t" << segs_curr[4] << endl;
    }
    assert(s_info_change.size() == Ns);

    // TODO: need to check consistency of CN_CHANGE across all bins in the segment
    for(int j = 0; j < Ns; ++j){
        assert(s_info_change[j].size() == Nsite);  // all sites for a sample
        CN_CHANGE bin_cn = s_info_change[j][segs_curr[1]];  // copy number of the/total first bin in the segment

        if(debug > 1){
            // print_vector<CN_CHANGE>(s_info[j][segs_curr[1]]);
            cout << "\tSample " << j + 1 << " with CN : " << bin_cn << endl;

        }        
        seg_cn.push_back(bin_cn);
    }

    if(debug > 1 ){
        cout << "\nFinal segment copy number for this segment: ";
        for(int j = 0; j < Ns; ++j) cout << "\t" << seg_cn[j];
        cout << endl;
        cout << "----------------------------------------" << endl;
    }

}


// each seg is represented by: {chr, id_start, id_end, seg_start, seg_end}
vector<vector<int>> get_all_segs(const vector<vector<vector<int>>>& s_info, int Ns, int num_total_bins, int& num_invar_bins, int incl_all, int is_total, int debug){
    num_invar_bins = 0;
    vector<int> var_bins(num_total_bins, 0);
    get_var_bins(s_info, Ns, num_total_bins, num_invar_bins, var_bins, is_total, debug);

    vector<vector<int>> segs;
    if(incl_all){   // bins with normal CNs across all samples are also included in the CN matrix
        for(int k = 0; k < num_total_bins; k++){
            int chr = s_info[0][k][0];
            int seg_start = s_info[0][k][1];   // original segment ID read from input, so starts from 1
            int id_start = k;  // new segment ID for tracking bin merging, starts from 0
            int seg_end = s_info[0][k][1];
            int id_end = k;
            vector<int> seg{chr, id_start, id_end, seg_start, seg_end};
            segs.push_back(seg);
        }
    }else{  // normal sites will be accounted for based on the count
        for(int k = 0; k < num_total_bins; k++){
            if(var_bins[k]){
                int chr = s_info[0][k][0];
                int seg_start = s_info[0][k][1];
                int id_start = k;
                int seg_end = s_info[0][k][1];
                int id_end = k;
                vector<int> seg{chr, id_start, id_end, seg_start, seg_end};
                segs.push_back(seg);
            }
        }
    }

    cout << "\tTotal number of segments:\t\t" << segs.size() << endl;

    return segs;
}


// cn_max: maximum total copy number or maximum state ID for haplotype-specific copy number, used to validate input CNs
vector<vector<int>> group_segs(const vector<vector<int>>& segs, const vector<vector<vector<int>>>& s_info, int Ns, int cn_max, int debug){
    vector<vector<int>> ret;

    if(debug){
        cout << "\nGrouping segments to get per-sample copy numbers" << endl;
        print_nested_vector<int>(segs);
        print_nested_vector2<int>(s_info);
    }

    for(size_t i = 0; i < segs.size(); ++i){
        vector<int> seg_cn; // for CNs of all samples in this segment
        vector<int> segs_curr = segs[i];
        compute_segment_cn_state(seg_cn, segs_curr, s_info, Ns, cn_max, debug);

        // chr, start, end
        vector<int> vals{segs[i][0], segs[i][1], segs[i][2]};
        for(int j = 0; j < Ns; ++j){
        	int cn = (int) seg_cn[j];
        	vals.push_back(cn);
        }

        ret.push_back(vals);
    }

    cout << "\nUsing segments:\t\t" << ret.size() << endl;
    if(debug > 2){
        for(int j = 0; j < ret.size(); ++j){
             for(int k = 0; k < Ns; ++k){
                 cout << "\t" << ret[j][k+3];
             }
             cout << endl;
        }
    }

    return ret;
}



/** 
 * @brief Merge consecutive bins with the same copy numbers across all samples      
 * @param segs Input segments represented by: {chr, id_start, id_end, seg_start, seg_end}
 * @param s_info Sample information where each innermost vector contains copy number state (int) or change (CN_CHANGE) for a segment in the format [chr, sid, cn] for each sample.
 * @param Ns Number of samples.
 * @param cn_max Maximum copy number.
 * @param seg_file Output file name to write the segments with copy numbers for each sample.        
 * @param debug Debug flag for verbose output.
 * @return A map where each key is a chromosome number and the value is a vector of segments with their copy numbers for each sample.
 */
void group_segs_by_chr_state(const vector<vector<int>>& segs, const vector<vector<vector<int>>>& s_info, map<int, vector<vector<int>>>& data, int Ns, int cn_max, const string& seg_file, int debug){
    int Nchar = 0;

    ofstream fcn;
    if(seg_file != "")  fcn.open(seg_file);

    for(size_t i = 0; i < segs.size(); ++i){
        vector<int> seg_cn;
        vector<int> segs_curr = segs[i];
        compute_segment_cn_state(seg_cn, segs_curr, s_info, Ns, cn_max, debug);

        // bin ID plus 1 to match original input format
        if(seg_file != "") fcn << segs_curr[0] << "\t" << segs_curr[1] + 1 << "\t" << segs_curr[2] + 1 << "\t" << segs_curr[3] << "\t" << segs_curr[4];

        vector<int> vals{segs_curr[0], segs_curr[1], segs_curr[2]};   // chr, start, end

        for(int j = 0; j < Ns; ++j){
            // for copy number values 
            int cn = seg_cn[j];
            vals.push_back(cn);

            if(seg_file != "") fcn << "\t" << cn;
        }
        if(seg_file != "") fcn << endl;

        data[segs[i][0]].push_back(vals);
        Nchar += 1;
    }

    if(seg_file != "") fcn.close();
    // output segments in a file for reference by site, which can be converted to the same format as original input
    cout << "\testmu\t\t" << Nchar << endl;
    if(debug > 2){
        // row: sites, column: samples, value: CN
        for(auto it : data){
            vector<vector<int>> sites = it.second;
            for(int j = 0; j < sites.size(); ++j){
                 cout << sites[j][0] << "\t" << sites[j][1] + 1 << "\t" << sites[j][2] + 1;
                 for(int k = 0; k < Ns; ++k){
                     cout << "\t" << sites[j][k+3];
                 }
                 cout << endl;
            }
        }
    }

}



/** 
 * @brief Merge consecutive bins with the same copy numbers across all samples      
 * @param segs Input segments represented by: {chr, id_start, id_end, seg_start, seg_end}
 * @param s_info Sample information where each innermost vector contains copy number state (int) or change (CN_CHANGE) for a segment in the format [chr, sid, cn] for each sample.
 * @param Ns Number of samples.
 * @param cn_max Maximum copy number.
 * @param seg_file Output file name to write the segments with copy numbers for each sample.        
 * @param debug Debug flag for verbose output.
 * @return A map where each key is a chromosome number and the value is a vector of segments with their copy numbers for each sample.
 */
void group_segs_by_chr_change(const vector<vector<int>>& segs, const vector<vector<CN_CHANGE>>& s_info_change, map<int, vector<vector<int>>>& data_change, int Ns, int Nsite, const string& seg_file, int debug){
    if(debug) cout << "\nGrouping segments by chromosome with CN_CHANGE type" << endl;
    int Nchar = 0;

    ofstream fcn;
    if(seg_file != ""){
        fcn.open(seg_file);
        if(!fcn.is_open()){
            cerr << "Failed to open file: " << seg_file << endl;
            exit(EXIT_FAILURE);
        }
    }  

    for(size_t i = 0; i < segs.size(); ++i){
        vector<CN_CHANGE> seg_cn;
        vector<int> segs_curr = segs[i];
        compute_segment_cn_change(seg_cn, segs_curr, s_info_change, Ns, Nsite, debug);
        assert(seg_cn.size() == Ns);

        if(seg_file != "") fcn << segs_curr[0] << "\t" << segs_curr[1] + 1 << "\t" << segs_curr[2] + 1 << "\t" << segs_curr[3] << "\t" << segs_curr[4];

        vector<int> vals{segs_curr[0], segs_curr[1], segs_curr[2]};   // chr, start, end
        // if(debug > 1) print_vector<int>(vals);

        for(int j = 0; j < Ns; ++j){
            // for copy number changes
            CN_CHANGE cn = seg_cn[j];
            vector<int> cn_change = cn.to_vector();
            vals.insert(vals.end(), cn_change.begin(), cn_change.end());
            
            if(debug > 1){
                cout << "Segment " << i + 1 << ", Sample " << j + 1 << ": " << cn << endl;
                print_vector<int>(cn_change);       
            }            
            
            if(seg_file != "") fcn << "\t" << cn;
        }
        if(seg_file != "") fcn << endl;

         if(debug > 1) print_vector<int>(vals);
        assert(vals.size() == 3 + 4 * Ns);
        
        data_change[segs_curr[0]].push_back(vals);
        Nchar += 1;
    }

    if(seg_file != "") fcn.close();
    // output segments in a file for reference by site, which can be converted to the same format as original input
    cout << "\nUsing segments:\t\t" << Nchar << endl;
    // if(debug){
    //     // row: sites, column: samples, value: CN
    //     for(auto it : res){
    //         vector<vector<int>> sites = it.second;
    //         for(int j = 0; j < sites.size(); ++j){
    //              cout << sites[j][0] << "\t" << sites[j][1] + 1 << "\t" << sites[j][2] + 1;
    //              for(int k = 0; k < Ns; ++k){
    //                  cout << "\t" << sites[j][k+3];
    //              }
    //              cout << endl;
    //         }
    //     }

}


/** 
 * @brief Read the input copy number file and merge bins into segments if specified
 * @param s_info Output parameter to store copy number information for each sample.
 * @param segs Output parameter to store segment information.
 * @param filename Input file name containing copy number data.
 * @param input_property Structure containing properties of the input data.
 * @param input_data Structure to store processed input data.
 * @param debug Debug flag for verbose output.
 * @return The maximum copy number state used in the analysis of bounded models. Different from cn_max in the input_property when haplotype-specific CNs are used.
 */
int get_segs_cn(vector<vector<vector<int>>>& s_info, vector<vector<int>>& segs, const string& filename, const INPUT_PROPERTY& input_property, INPUT_DATA& input_data, int debug){
    int cn_max = input_property.cn_max;

    // get CNs for all the bins or sites
    read_cn(s_info, filename, input_property.Ns, input_data.num_total_bins, cn_max, input_property.is_total, input_property.is_rcn, debug);

    // find location of segments (consecutive sites) first
    if(input_property.is_bin){
        // Read the input copy numbers while converting runs of variable bins into segments of constant cn values and group them by chromosome
        segs = get_invar_segs(s_info, input_property.Ns, input_data.num_total_bins, input_data.num_invar_bins, input_property.is_total, debug);
    }else{
        // Read the input copy numbers as they are and group them by chromosome
        segs = get_all_segs(s_info, input_property.Ns, input_data.num_total_bins, input_data.num_invar_bins, input_property.incl_all, input_property.is_total, debug);
    }

    input_data.seg_size = segs.size();

    int max_cn_state = cn_max;
    if(!input_property.is_total){   // haplotype-specific CN states in the rate matrix
        max_cn_state = (cn_max + 1) * (cn_max + 2) / 2 - 1;
    }

  return max_cn_state;
}



/************************* Functions called by tree building program *************************/
// only used in testing ML tree building
vector<vector<int>> read_data_var_regions(const string& filename, const INPUT_PROPERTY& input_property, INPUT_DATA& input_data, int debug){
  cout << "\nReading data and calculating CNA regions" << endl;

  vector<vector<vector<int>>> s_info;
  vector<vector<int>> segs;
  int max_cn_state = get_segs_cn(s_info, segs, filename, input_property, input_data, debug);

  vector<vector<int>> ret = group_segs(segs, s_info, input_property.Ns, max_cn_state, debug);

  return ret;
}



// // group observed copy numbers by chromosome for tree inference
// // return type: vobs[chr] = { {cn_sample1, cn_sample2, ...}, ...} }
// map<int, vector<vector<int>>> get_obs_vector_by_chr(map<int, vector<vector<int>>>& data, const int Ns){
//     map<int, vector<vector<int>>> vobs;
//     // Construct the CN matrix by chromosome
//     // Assume chromosomes in data are ordered numerically
//     // int total_chr = data.rbegin()->first;  // Some chromosomes got lost in the segment merging, so total_chr may not equal to data.size()
//     // int nchr = data.begin()->first;
//     // for(; nchr <= total_chr; nchr++){
//     for(auto dcn : data){
//         int nchr = dcn.first;
//         vector<vector<int>> obs_chr;
//         for(int nc = 0; nc < data[nchr].size(); ++nc){
//             vector<int> obs;
//             for(int i = 0; i < Ns; ++i){
//               obs.push_back(data[nchr][nc][i + 3]);
//             }
//             obs_chr.push_back(obs);
//         }
//         vobs[nchr] = obs_chr;
//     }
//     return vobs;
// }


/** 
 *  @brief Group observed copy numbers by chromosome for tree inference
 *  @param data Input data map where each key is a chromosome number and the value is a vector of vectors containing copy number information.
 *  @param Ns Number of samples.
 *  @tparam T Type of the copy number data (e.g., int, CN_CHANGE), i.e. vector<vector<CN_CHANGE>> or vector<vector<int>
 *  @return A map where each key is a chromosome number and the value is a vector of vectors containing copy number observations for each sample.
 *  return type: vobs[chr] = { {cn_sample1, cn_sample2, ...}, ...} } or vobs[chr] = { {cn_change_sample1, cn_change_sample2, ...}, ...} }
 */
void get_obs_vector_by_chr_state(const map<int, vector<vector<int>>>& data,  map<int, vector<vector<int>>>& vobs, int Ns){
    // Construct the CN matrix by chromosome
    // Assume chromosomes in data are ordered numerically
    // int total_chr = data.rbegin()->first;  // Some chromosomes got lost in the segment merging, so total_chr may not equal to data.size()
    // int nchr = data.begin()->first;
    // for(; nchr <= total_chr; nchr++){
    for(auto dcn : data){
        int nchr = dcn.first;
        vector<vector<int>> obs_chr;
        for(size_t nc = 0; nc < data.at(nchr).size(); ++nc){  // all sites on this chr
            vector<int> obs;
            // for copy number values
            for(int i = 0; i < Ns; ++i){
                obs.push_back(data.at(nchr).at(nc)[i + 3]);  // The first three elements in the vector are: chr, start, end; so copy numbers start from index 3
            }
            obs_chr.push_back(obs);
        }
        vobs[nchr] = obs_chr;
    }
}


/** 
 *  @brief Group observed copy number changes by chromosome for tree inference
 *  @param data Input data map where each key is a chromosome number and the value is a vector of vectors containing copy number information.
 *  @param Ns Number of samples.
 *  @tparam T Type of the copy number data (e.g., int, CN_CHANGE), i.e. vector<vector<CN_CHANGE>> or vector<vector<int>
 *  @return A map where each key is a chromosome number and the value is a vector of vectors containing copy number observations for each sample.
 *  return type: vobs[chr] = { {cn_sample1, cn_sample2, ...}, ...} } or vobs[chr] = { {cn_change_sample1, cn_change_sample2, ...}, ...} }
 */
void get_obs_vector_by_chr_change(const map<int, vector<vector<int>>>& data_change,  map<int, vector<vector<CN_CHANGE>>>& vobs_change, int Ns){
    // Construct the CN matrix by chromosome
    // Assume chromosomes in data are ordered numerically
    for(auto dcn : data_change){
        int nchr = dcn.first;
        vector<vector<CN_CHANGE>> obs_chr;
        for(size_t nc = 0; nc < data_change.at(nchr).size(); ++nc){  // all sites on this chr
            vector<CN_CHANGE> obs;
            // cout << "data[" << nchr << "][" << nc << "]: ";
            // print_vector<int>(data[nchr][nc]);            
            int j = 3;  // The first three elements in the vector are: chr, start, end; so copy numbers start from index 3
            for(int i = 0; i < Ns; i++){  // CN_CHANGE has 4 elements
                CN_CHANGE cn_change = {data_change.at(nchr).at(nc)[j], data_change.at(nchr).at(nc)[j + 1], data_change.at(nchr).at(nc)[j + 2], data_change.at(nchr).at(nc)[j + 3]};  // The first three elements in the vector are: chr, start, end; so copy numbers start from index 3
                j = j + 4; 
                obs.push_back(cn_change);
            }

            obs_chr.push_back(obs);
        }
        // will add new values when using []
        vobs_change[nchr] = obs_chr;
    }
}

/** 
 *  @brief Read data and group regions by chromosome, used in tree building for more comprehensive data processing  
 *  @param res Output parameter to store the grouped copy number information by chromosome.
 *  @param filename Input file name containing copy number data.
 *  @param input_property Structure containing properties of the input data.
 *  @param input_data Structure to store processed input data.
 *  @param seg_file Optional output file name to save segment information.
 *  @param debug Debug flag for verbose output.
 *  @return A map where each key is a chromosome number and the value is a vector of vectors containing copy number information for each segment on that chromosome.
*/
void read_data_var_regions_by_chr_state(map<int, vector<vector<int>>>& data, const string& filename, const INPUT_PROPERTY& input_property, INPUT_DATA& input_data, const string& seg_file, int debug){
    cout << "\nReading data and group regions by chromosome" << endl;
    vector<vector<vector<int>>> s_info; 
    vector<vector<int>> segs;

    int max_cn_state = get_segs_cn(s_info, segs, filename, input_property, input_data, debug);

    get_sample_ploidy(s_info, input_data.sample_avg_cn, input_data.sample_chr_cn, input_property.cn_max, input_property.is_total, debug);

    // get_sample_mcn(s_info, input_data.sample_max_cn, cn_max, input_property.is_total, debug);
    get_site_change(s_info, input_data.sample_avg_cn, input_data.sample_change_site, input_data.site_max_change, input_data.sample_max_cn, input_property.cn_max, input_property.is_total, debug);

    if(input_property.model == DECOMP){
        get_num_wgd(s_info, input_data.sample_avg_cn, input_data.sample_num_wgd, debug);
        get_chr_change(s_info, input_data.sample_avg_cn, input_data.sample_chr_cn, input_data.sample_change_chr, input_data.chr_max_change, input_property.cn_max, input_property.is_total, debug);
    }

    // combine segment locations and site-level CNs to get the final CN matrix
    group_segs_by_chr_state(segs, s_info, data, input_property.Ns, max_cn_state, seg_file, debug);
    if(debug){
        cout << "Reading input copy number state finished" << endl;
    }
}




/** 
 *  @brief Read data and group regions by chromosome, used in tree building for more comprehensive data processing  
 *  @param res Output parameter to store the grouped copy number change information by chromosome.
 *  @param filename Input file name containing copy number data.
 *  @param input_property Structure containing properties of the input data.
 *  @param input_data Structure to store processed input data, will be changed during processing.
 *  @param seg_file Optional output file name to save segment information.
 *  @param debug Debug flag for verbose output.
 *  @return A map where each key is a chromosome number and the value is a vector of vectors containing copy number information for each segment on that chromosome.
*/
void read_data_var_regions_by_chr_change(map<int, vector<vector<int>>>& data_change, const string& filename, const INPUT_PROPERTY& input_property, INPUT_DATA& input_data, const string& seg_file, int debug){
    cout << "\nReading data and group regions by chromosome" << endl;
    vector<vector<vector<int>>> s_info; 
    vector<vector<int>> segs;
    // convert copy numbers to copy number changes for decomposition model
    vector<vector<CN_CHANGE>> s_info_change;  

    int max_cn_state = get_segs_cn(s_info, segs, filename, input_property, input_data, debug);

    if(!input_property.is_rcn){
        get_sample_ploidy(s_info, input_data.sample_avg_cn, input_data.sample_chr_cn, input_property.cn_max, input_property.is_total, debug);

        // get_sample_mcn(s_info, input_data.sample_max_cn, cn_max, input_property.is_total, debug);
        get_site_change(s_info, input_data.sample_avg_cn, input_data.sample_change_site, input_data.site_max_change, input_data.sample_max_cn, input_property.cn_max, input_property.is_total, debug);

        get_num_wgd(s_info, input_data.sample_avg_cn, input_data.sample_num_wgd, debug);
        get_chr_change(s_info, input_data.sample_avg_cn, input_data.sample_chr_cn, input_data.sample_change_chr, input_data.chr_max_change, input_property.cn_max, input_property.is_total, debug);

        if(input_property.is_total){
            cout << "Converting total copy numbers to copy number changes for decomposition model" << endl;
            cn_to_decomposition(s_info, s_info_change, input_data, debug);
        }else{
            cout << "TODO: Converting haplotype-specific copy numbers to copy number changes for decomposition model" << endl;
        }
    }else{
        cout << "Converting relative copy numbers to copy number changes for decomposition model" << endl;
        rcn_to_decomposition(s_info, s_info_change, input_data.sample_change_site, input_data.site_max_change, input_data.sample_max_cn, debug);
    }
    
    group_segs_by_chr_change(segs, s_info_change, data_change, input_property.Ns, input_data.num_total_bins, seg_file, debug);
    
    cout << "=============Reading input copy number change finished==============\n" << endl;
}
