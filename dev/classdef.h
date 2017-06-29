#ifndef __CLASSDEF_H_
#define __CLASSDEF_H_

#include <stdlib.h>
#include <string>
#include <string.h>
#include <vector>
#include <map>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>


using namespace std;




class SNP {

 public:

  string id;
  int index;
  int dtss_bin; // distance to TSS, used by fastQTL


  // for grid estimating
  int grid_size;
  vector<double> log10_BF_vec;
  double *weightv;

  // possible imput
  double log10_BF;
  //
  double zval;
  //
  double beta;
  double se_beta;
  
  SNP(string snp_id, int snp_index){
    id = snp_id;
    index = snp_index;
  }




  // version with fixed BF input
  SNP(string snp_id, double snp_log10_BF, int snp_index){
    id = snp_id;
    log10_BF = snp_log10_BF;
    index = snp_index;
    grid_size = 1;
  }

  SNP(string snp_id, vector<double> snp_log10_BF_vec, int snp_index, double *wv){
    id = snp_id;
    log10_BF_vec  = snp_log10_BF_vec;
    index = snp_index;
    weightv = wv;
    grid_size = log10_BF_vec.size();
    if(grid_size==1){
      log10_BF = log10_BF_vec[0];
    }
  }
    

  double compute_log10_BF();

  void compute_log10_BF_vec_zval(vector<double> & grid);
  void compute_log10_BF_vec_bhat(vector<double> & grid);
    
 
};





class Locus {
  
 public:
  
  // possible gene information/annotation
  string id;
  
  vector<SNP> snpVec;

  vector<double> baseline_vec;
  
  vector<vector<double> > BF_wv_vec; 
  
  gsl_vector *prior_vec;
  gsl_vector *pip_vec;

  vector<double> wv_probv;
  
  double log10_lik; // log10 of marginal likelihood
  vector<double> log10_lik_wv; // for single_wv estimation
  vector<double> prior_pvec; // for single_wv estimation
  
  double fdr;
  
  
  Locus(string locus_id,  vector<SNP> & snpVec_){ id = locus_id;  snpVec = snpVec_; prior_vec = pip_vec =0;  };
  Locus(){};
  
  void get_BF_matrix();
  void EM_update();
  void compute_fdr();


};






class controller {

 public:
 
  controller(){
    p=kc=kd=dist_bin_level = 0;
    force_logistic = 0;
    dist_bin_size = -1;
    fastqtl_use_dtss = 0;
    dist_bin = 0;
    finish_em = 0;
    single_fuzzy_annot = 0;
    l1_lambda=l2_lambda=0;
    init_pi1 = 1e-3;
    print_avg = 0;
    print_lead_snp = 0;
    grid_size = 1;
    use_ash = 0;
    shrink_pi0 = 0.0;
  }

  // storage
  vector<Locus> locVec;
  
  
  int p; // number of loc-SNP pairs
  
  int kc; // number of continuous covariate
  int kd; // number of discrete covariate


  double dist_bin_size;


  map<string,int> loc_hash;
  map<string,int> snp_hash;
  vector<string> cvar_name_vec;
  vector<string> dvar_name_vec;



  
  gsl_vector_int *dist_bin;
  map<int, int> dtss_map;
  map<int, int> dtss_rmap;
  int dist_bin_level;
  double EM_thresh;
  
  gsl_vector_int *dlevel; // (kd+1) entry levels of each factor


  gsl_matrix *Xc;  // p x kc
  gsl_matrix_int *Xd; // p x kd
  
  
  gsl_vector *prior_vec;
  gsl_vector *pip_vec;
  gsl_vector *beta_vec;
  
  int ncoef;

  double final_log10_lik;
  

  double init_pi1;
  int nthread; 


  int finish_em;

  int fastqtl_use_dtss;
  int print_avg;
  int print_lead_snp;

  // for BSLMM grid
  double *grid_weight;
  int grid_size;

  vector<double> grid_vec;

  // use adaptive shrinkage for effect size estimation
  int use_ash;	

  double shrink_pi0;
  double pi0_ubd;


  void load_data(char *filename);    // load data with MatrixeQTL format -- default
  void load_data_zscore(char *filename); // load data with z-score/t-score
  void load_data_BF(char *filename); // load data with pre-computed log10 Bayes factors
  void load_data_fastqtl(char *filename); //load data with fastQTL format

  void load_data_BSLMM_BF(char *filename, int grid_size); // load data with BSLMM-precomputed BFs

  
  
  void load_map(char *gene_map_file, char *snp_map_file);
  void load_annotation(char *annot_file);
  int count_factor_level(int col);
  
  void simple_regression();
  void single_ct_regression();
  void single_probt_est();
  void single_probt_regression();
  int force_logistic;
  int single_fuzzy_annot;

  void update_grid_weight();
 
  
  double l1_lambda, l2_lambda; //shrinkage for enrich param est
  

  void init_params();
    
  void run_EM();

  
  void find_eGene(double thresh=0.05);
  void estimate();
  void dump_prior(char *path);
  void dump_pip(char *file);
  void dump_wv(char *file);
  
  vector<double> make_grid(double phi_min, double phi_max);

 private:
  double eval_likelihood(double x, int index);
  double fine_optimize_beta(int index, double est, double null_log10_lik, double &curr_log10_lik);
  



};




double log10_weighted_sum(vector<double> & val_vec, vector<double> & wts_vec);
double log10_weighted_sum(vector<double> &vec, double *wts);
double compute_log10_BF(double beta, double se_beta);
double compute_log10_BF(double z_score);
bool   rank_by_fdr (const Locus & loc1 , const Locus & loc2);
bool   sort_by_BF(const SNP &snp1, const SNP &snp2);  
int classify_dist_bin(int snp_pos, int tss, double bin_size = -1);
double map_bin_2_dist(int bin, double bin_size=-1);
#endif
