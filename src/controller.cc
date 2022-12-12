#include "classdef.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <omp.h>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include "logistic.h"
#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>

void controller::load_data(char *filename){



    ifstream dfile(filename, ios_base::in | ios_base::binary);
    boost::iostreams::filtering_istream in;

    in.push(boost::iostreams::gzip_decompressor());
    in.push(dfile);


    string line;
    istringstream ins;

    string curr_loc_id = "";

    vector<SNP> snp_vec;
    string loc_id;
    string snp_id;
    double beta;
    double se_beta;
    double t_val;
    int index_count = 0;
    int loc_count = 0;

    double median_zcount = 0;


    // for grid setting
    double min = 1; 
    double max = 0;; 
    while(getline(in,line)){

        ins.clear();
        ins.str(line);

        if(ins>>snp_id>>loc_id>>beta>>t_val){

            //printf("%s  %f\n", snp_id.c_str(), beta);
            if(curr_loc_id != loc_id){
                if(curr_loc_id != ""){

                    if(loc_hash.find(curr_loc_id)== loc_hash.end()){
                        Locus loc(curr_loc_id, snp_vec);
                        int pos = loc_hash.size();
                        loc_hash[curr_loc_id] = pos;
                        locVec.push_back(loc);	  
                        loc_count++;
                    }else{
                        // already defined, concate the existing snpVec
                        int pos = loc_hash[curr_loc_id];
                        for(int i=0;i<snp_vec.size();i++){
                            locVec[pos].snpVec.push_back(snp_vec[i]);
                        }	   
                    }

                    snp_vec.clear();

                }
                curr_loc_id = loc_id;

            }

            if(t_val == 0)
                t_val = 1e-8;

            se_beta = beta/t_val;
            SNP snp(snp_id, index_count);

            if(shrink_pi0>0 && pow(t_val,2)<= 0.4549){     
                median_zcount++;
            } 
            if(use_ash){
                if(se_beta<min)
                    min = se_beta;

                if(beta*beta - se_beta*se_beta > max){
                    max = pow(beta,2)-pow(se_beta,2);
                }

                snp.beta = beta;
                snp.se_beta = se_beta;
            }else{
                snp.log10_BF = compute_log10_BF(beta, se_beta);
                snp.grid_size = 1;
            }
            index_count++;
            snp_hash[snp_id] = 100;
            snp_vec.push_back(snp);

        }
    }

    dfile.close();

    if(loc_hash.find(curr_loc_id)== loc_hash.end()){
        Locus loc(curr_loc_id, snp_vec);
        int pos = loc_hash.size();
        loc_hash[curr_loc_id] = pos;
        locVec.push_back(loc);
        loc_count++;
    }else{
        // already defined, concate the existing snpVec                                               
        int pos = loc_hash[curr_loc_id];
        for(int i=0;i<snp_vec.size();i++){
            locVec[pos].snpVec.push_back(snp_vec[i]);
        }
    }


    if(use_ash){
        // compute BF vector
        double phi_min = min/10;
        double phi_max = 2*sqrt(max);
        if(phi_max<phi_min){
            phi_max = 8*phi_min;
        }

        grid_vec = make_grid(phi_min, phi_max);

        grid_size = grid_vec.size();

        grid_weight = new double[grid_size];
        for(int i=0;i<grid_size;i++){
            grid_weight[i] = 1.0/double(grid_size);
        }
        for(int i=0;i<locVec.size();i++){
            for(int j=0;j<locVec[i].snpVec.size();j++){
                locVec[i].snpVec[j].compute_log10_BF_vec_bhat(grid_vec);
                locVec[i].snpVec[j].weightv = grid_weight;
            }	
            locVec[i].get_BF_matrix(); 
        }	
    }

    p = index_count;

    if(shrink_pi0>0){
        if(set_pi0 > 0){
            pi0_ubd = set_pi0;
        }else{
            pi0_ubd = 2*median_zcount/p;
        }
    }	

    fprintf(stderr, "Read in %d loci, %d locus-SNP pairs ... \n",loc_count, p);

    prior_vec = gsl_vector_calloc(p);

}


void controller::load_data_fastqtl(char *filename){

    // fastQTL file contains dtss info, automatically parse and use dtss in analysis

    ifstream dfile(filename, ios_base::in | ios_base::binary);
    boost::iostreams::filtering_istream in;

    in.push(boost::iostreams::gzip_decompressor());
    in.push(dfile);


    string line;
    istringstream ins;

    string curr_loc_id = "";

    vector<SNP> snp_vec;
    string loc_id;
    string snp_id;

    double bhat;
    double sdbeta;
    double pval;
    double dtss;

    int index_count = 0;
    int loc_count = 0;

    map<int, int> bin_hash;

    double min = 1;
    double max = 0;

    double median_zcount = 0;


    while(getline(in,line)){

        ins.clear();
        ins.str(line);

        if(ins>> loc_id >> snp_id >> dtss >>pval >> bhat >> sdbeta){
            //printf("%s  %f\n", snp_id.c_str(), beta);                                                         
            if(curr_loc_id != loc_id){
                if(curr_loc_id != ""){

                    if(loc_hash.find(curr_loc_id)== loc_hash.end()){
                        Locus loc(curr_loc_id, snp_vec);
                        int pos = loc_hash.size();
                        loc_hash[curr_loc_id] = pos;
                        locVec.push_back(loc);
                        loc_count++;
                    }else{
                        // already defined, concate the existing snpVec                                              
                        int pos = loc_hash[curr_loc_id];
                        for(int i=0;i<snp_vec.size();i++){
                            locVec[pos].snpVec.push_back(snp_vec[i]);
                        }
                    }

                    snp_vec.clear();

                }
                curr_loc_id = loc_id;

            }

            //double log10_BF = compute_log10_BF(bhat, sdbeta);
            SNP snp(snp_id, index_count);
            if(shrink_pi0>0 && pow(bhat/sdbeta,2)<= 0.4549){
                median_zcount++;
            }
            if(use_ash){
                snp.beta = bhat;
                snp.se_beta = sdbeta;
                if(sdbeta<min){
                    min = sdbeta;
                }		
                if(bhat*bhat-sdbeta*sdbeta > max){
                    max = pow(bhat,2) - pow(sdbeta,2);
                }
            }else{
                snp.log10_BF = compute_log10_BF(bhat, sdbeta);
                snp.grid_size = 1;
            }	
            // handling dtss info
            if(fastqtl_use_dtss){
                int bin = classify_dist_bin(dtss, 0 , dist_bin_size);
                if(bin_hash.find(bin)==bin_hash.end()){
                    bin_hash[bin] = 0;
                }

                bin_hash[bin]++;
                snp.dtss_bin = bin;
            }

            index_count++;
            snp_hash[snp_id] = 100;
            snp_vec.push_back(snp);

        }
    }

    dfile.close();

    // record the last loc

    if(loc_hash.find(curr_loc_id)== loc_hash.end()){
        Locus loc(curr_loc_id, snp_vec);
        int pos = loc_hash.size();
        loc_hash[curr_loc_id] = pos;
        locVec.push_back(loc);
        loc_count++;
    }else{
        // already defined, concate the existing snpVec                                                       
        int pos = loc_hash[curr_loc_id];
        for(int i=0;i<snp_vec.size();i++){
            locVec[pos].snpVec.push_back(snp_vec[i]);
        }
    }

    if(use_ash){
        // compute BF vector
        double phi_min = min/10;
        double phi_max = 2*sqrt(max);
        if(phi_max<phi_min){
            phi_max = 8*phi_min;
        }	

        grid_vec = make_grid(phi_min, phi_max);
        grid_size = grid_vec.size();
        grid_weight = new double[grid_size];
        for(int i=0;i<grid_size;i++){
            grid_weight[i] = 1.0/double(grid_size);
        }

        for(int i=0;i<locVec.size();i++){
            for(int j=0;j<locVec[i].snpVec.size();j++){
                locVec[i].snpVec[j].compute_log10_BF_vec_bhat(grid_vec);
                locVec[i].snpVec[j].weightv = grid_weight;
            }		                    
            locVec[i].get_BF_matrix();
        }

    }

    p = index_count;

    if(shrink_pi0>0){
        if(set_pi0 > 0){
            pi0_ubd = set_pi0;
        }else{
            pi0_ubd = 2*median_zcount/p;
        }
    }	



    if(fastqtl_use_dtss){
        dtss_map[0] =0;
        dtss_rmap[0] = 0;
        int count = 1;

        dist_bin = gsl_vector_int_calloc(p);

        // re-map bin number to skip empty bins
        for (map<int,int>::iterator it=bin_hash.begin(); it!=bin_hash.end(); ++it){
            if(it->first==0)
                continue;
            dtss_map[it->first] = count;
            dtss_rmap[count] = it->first;
            count++;
        }
        dist_bin_level = count;

        for (int i=0;i<locVec.size();i++){

            for(int j=0;j<locVec[i].snpVec.size();j++){
                string snp_id = locVec[i].snpVec[j].id;
                gsl_vector_int_set(dist_bin,locVec[i].snpVec[j].index, dtss_map[locVec[i].snpVec[j].dtss_bin]);
            }
        }
    }


    fprintf(stderr, "Read in %d loci, %d locus-SNP pairs ... \n",loc_count, p);

    prior_vec = gsl_vector_calloc(p);



}








void controller::load_data_zscore(char *filename){



    ifstream dfile(filename, ios_base::in | ios_base::binary);
    boost::iostreams::filtering_istream in;

    in.push(boost::iostreams::gzip_decompressor());
    in.push(dfile);


    string line;
    istringstream ins;

    string curr_loc_id = "";

    vector<SNP> snp_vec;
    string loc_id;
    string snp_id;
    double z_val;
    int index_count = 0;
    int loc_count = 0;
    double max = 0;
    double median_zcount = 0;
    while(getline(in,line)){

        ins.clear();
        ins.str(line);

        if(ins>>snp_id>>loc_id>>z_val){

            //printf("%s  %f\n", snp_id.c_str(), z_val);

            if(curr_loc_id != loc_id){
                if(curr_loc_id != ""){

                    if(loc_hash.find(curr_loc_id)== loc_hash.end()){
                        Locus loc(curr_loc_id, snp_vec);
                        int pos = loc_hash.size();
                        loc_hash[curr_loc_id] = pos;
                        locVec.push_back(loc);	  
                        loc_count++;
                    }else{
                        // already defined, concate the existing snpVec
                        int pos = loc_hash[curr_loc_id];
                        for(int i=0;i<snp_vec.size();i++){
                            locVec[pos].snpVec.push_back(snp_vec[i]);
                        }	   
                    }

                    snp_vec.clear();

                }
                curr_loc_id = loc_id;

            }

            //double log10_BF = compute_log10_BF(z_val);
            //printf("%15s   %7.3f\n", snp_id.c_str(), log10_BF);
            SNP snp(snp_id, index_count);

            if(shrink_pi0>0 && pow(z_val,2)<= 0.4549){
                median_zcount++;
            }

            if(use_ash){	
                snp.zval = z_val;
                if(z_val*z_val - 1 > max){
                    max = z_val*z_val - 1;
                }	
            }else{
                snp.log10_BF = compute_log10_BF(z_val);
                snp.grid_size = 1;
            }
            index_count++;
            snp_hash[snp_id] = 100;
            snp_vec.push_back(snp);

        }
    }

    dfile.close();

    if(loc_hash.find(curr_loc_id)== loc_hash.end()){
        Locus loc(curr_loc_id, snp_vec);
        int pos = loc_hash.size();
        loc_hash[curr_loc_id] = pos;
        locVec.push_back(loc);
        loc_count++;
    }else{
        // already defined, concate the existing snpVec                                               
        int pos = loc_hash[curr_loc_id];
        for(int i=0;i<snp_vec.size();i++){
            locVec[pos].snpVec.push_back(snp_vec[i]);
        }
    }

    if(use_ash){ 
        // compute BF vector
        double phi_min = 1.0/10;
        double phi_max = 2*sqrt(max);
        if(phi_max<phi_min){
            phi_max = 8*phi_min;
        }

        grid_vec = make_grid(phi_min, phi_max);


        grid_size = grid_vec.size();
        grid_weight = new double[grid_size];
        for(int i=0;i<grid_size;i++){
            grid_weight[i] = 1.0/double(grid_size);
        }

        for(int i=0;i<locVec.size();i++){
            for(int j=0;j<locVec[i].snpVec.size();j++){
                locVec[i].snpVec[j].compute_log10_BF_vec_zval(grid_vec);
                locVec[i].snpVec[j].weightv = grid_weight;
            }
            locVec[i].get_BF_matrix();
        }	

    }

    p = index_count;

    if(shrink_pi0>0){
        if(set_pi0 > 0){
            pi0_ubd = set_pi0;
        }else{
            pi0_ubd = 2*median_zcount/p;
        }
    }	



    fprintf(stderr, "Read in %d loci, %d locus-SNP pairs ... \n", loc_count, p);

    prior_vec = gsl_vector_calloc(p);

}






// this function directly loads pre-computed log10 Bayes factors

void controller::load_data_BF(char *filename){



    ifstream dfile(filename, ios_base::in | ios_base::binary);
    boost::iostreams::filtering_istream in;

    in.push(boost::iostreams::gzip_decompressor());
    in.push(dfile);


    string line;
    istringstream ins;

    string curr_loc_id = "";

    vector<SNP> snp_vec;
    string loc_id;
    string snp_id;
    double log10_BF;
    int index_count = 0;
    int loc_count = 0;
    while(getline(in,line)){

        ins.clear();
        ins.str(line);

        if(ins>>snp_id>>loc_id>>log10_BF){

            //printf("%s  %f\n", snp_id.c_str(), log10_BF);

            if(curr_loc_id != loc_id){
                if(curr_loc_id != ""){

                    if(loc_hash.find(curr_loc_id)== loc_hash.end()){
                        Locus loc(curr_loc_id, snp_vec);
                        int pos = loc_hash.size();
                        loc_hash[curr_loc_id] = pos;
                        locVec.push_back(loc);	  
                        loc_count++;
                    }else{
                        // already defined, concate the existing snpVec
                        int pos = loc_hash[curr_loc_id];
                        for(int i=0;i<snp_vec.size();i++){
                            locVec[pos].snpVec.push_back(snp_vec[i]);
                        }	   
                    }

                    snp_vec.clear();

                }
                curr_loc_id = loc_id;

            }

            SNP snp(snp_id, log10_BF, index_count);
            index_count++;
            snp_hash[snp_id] = 100;
            snp_vec.push_back(snp);

        }
    }

    dfile.close();

    if(loc_hash.find(curr_loc_id)== loc_hash.end()){
        Locus loc(curr_loc_id, snp_vec);
        int pos = loc_hash.size();
        loc_hash[curr_loc_id] = pos;
        locVec.push_back(loc);
        loc_count++;
    }else{
        // already defined, concate the existing snpVec                                               
        int pos = loc_hash[curr_loc_id];
        for(int i=0;i<snp_vec.size();i++){
            locVec[pos].snpVec.push_back(snp_vec[i]);
        }
    }



    p = index_count;

    fprintf(stderr, "Read in %d loci, %d locus-SNP pairs ... \n",loc_count, p);

    prior_vec = gsl_vector_calloc(p);

}



// this function directly loads pre-computed log10 Bayes factors

void controller::load_data_BSLMM_BF(char *filename, int gsize){


    ifstream dfile(filename, ios_base::in | ios_base::binary);
    boost::iostreams::filtering_istream in;

    in.push(boost::iostreams::gzip_decompressor());
    in.push(dfile);


    string line;
    istringstream ins;

    string curr_loc_id = "";

    vector<SNP> snp_vec;
    vector<double> baseline_vec;
    string loc_id;
    string snp_id;
    double log10_BF;
    int index_count = 0;
    int loc_count = 0;


    // initialize grid_weight
    grid_size = gsize;
    grid_weight = new double[grid_size];
    for(int i=0;i<grid_size;i++){
        grid_weight[i] = 1.0/double(grid_size);
    }


    while(getline(in,line)){

        ins.clear();
        ins.str(line);

        if(ins>>snp_id>>loc_id){

            //printf("%s  %f\n", snp_id.c_str(), log10_BF);

            if(curr_loc_id != loc_id){
                if(curr_loc_id != ""){

                    if(loc_hash.find(curr_loc_id)== loc_hash.end()){
                        Locus loc(curr_loc_id, snp_vec);
                        if(baseline_vec.size()>0){
                            loc.baseline_vec = baseline_vec;
                        }
                        int pos = loc_hash.size();
                        loc_hash[curr_loc_id] = pos;
                        locVec.push_back(loc);
                        locVec[pos].get_BF_matrix();
                        loc_count++;
                    }else{
                        // already defined, concate the existing snpVec
                        int pos = loc_hash[curr_loc_id];
                        for(int i=0;i<snp_vec.size();i++){
                            locVec[pos].snpVec.push_back(snp_vec[i]);
                        }	   
                    }

                    snp_vec.clear();
                    baseline_vec.clear();
                }
                curr_loc_id = loc_id;

            }

            vector<double> bf_vec;
            double bf;
            for (int i=0;i<grid_size;i++){
                if(ins >> bf){
                    bf_vec.push_back(bf);
                }
            }
            if(snp_id == "BASELINE"){
                baseline_vec = bf_vec;
            }else{
                SNP snp(snp_id, bf_vec, index_count, grid_weight);
                index_count++;
                snp_hash[snp_id] = 100;
                snp_vec.push_back(snp);
            }
        }
    }

    dfile.close();

    if(loc_hash.find(curr_loc_id)== loc_hash.end()){
        Locus loc(curr_loc_id, snp_vec);
        if(baseline_vec.size()>0){
            loc.baseline_vec = baseline_vec;
        }
        int pos = loc_hash.size();
        loc_hash[curr_loc_id] = pos;
        locVec.push_back(loc);
        locVec[pos].get_BF_matrix();
        loc_count++;
    }else{
        // already defined, concate the existing snpVec                                               
        int pos = loc_hash[curr_loc_id];
        for(int i=0;i<snp_vec.size();i++){
            locVec[pos].snpVec.push_back(snp_vec[i]);
        }
    }



    p = index_count;

    fprintf(stderr, "Read in %d loci, %d locus-SNP pairs ... \n",loc_count, p);

    prior_vec = gsl_vector_calloc(p);

    /*
       for(int i=0;i<locVec.size();i++){
       for(int j=0;j<locVec[i].snpVec.size();j++){
       for(int k=0;k<grid_size;k++){
       printf("%7.2e  ", locVec[i].snpVec[j].log10_BF_vec[k]);
       printf("%7.2e    ", locVec[i].BF_wv_vec[k][j]);
       }
       }
       printf("\n");
       }
       */





}















void controller::load_map(char* gene_file, char *snp_file){

    if(fastqtl_use_dtss)
        return;

    if(strlen(gene_file)==0 || strlen(snp_file)==0){
        return;
    }



    map<string, int> gene_map;
    map<string, int> snp_map;

    ifstream gfile(gene_file, ios_base::in | ios_base::binary);
    boost::iostreams::filtering_istream in;

    in.push(boost::iostreams::gzip_decompressor());
    in.push(gfile);


    string line;
    istringstream ins;


    string gene_id;
    int chr;
    double tss1;
    double tss2;
    while(getline(in,line)){

        ins.clear();
        ins.str(line);

        if(ins>>gene_id>>chr>>tss1>>tss2){
            if(loc_hash.find(gene_id) != loc_hash.end()){
                gene_map[gene_id] = tss1;
            }
        }
    }
    gfile.close();


    ifstream sfile(snp_file, ios_base::in | ios_base::binary);
    boost::iostreams::filtering_istream snp_in;

    snp_in.push(boost::iostreams::gzip_decompressor());
    snp_in.push(sfile);

    istringstream snp_ins;


    string snp_id;
    double pos;

    while(getline(snp_in,line)){

        snp_ins.clear();
        snp_ins.str(line);

        if(snp_ins>>snp_id>>chr>>pos){
            if(snp_hash.find(snp_id)!=snp_hash.end()){
                snp_map[snp_id] = pos;
            }
        }
    }
    sfile.close();

    dist_bin = gsl_vector_int_calloc(p);
    map<int, int> bin_hash;

    for (int i=0;i<locVec.size();i++){
        string loc_id = locVec[i].id;
        for(int j=0;j<locVec[i].snpVec.size();j++){
            string snp_id = locVec[i].snpVec[j].id;
            int bin = classify_dist_bin(snp_map[snp_id], gene_map[loc_id],dist_bin_size);
            if(bin_hash.find(bin)==bin_hash.end()){
                bin_hash[bin] = 0;
            }
            bin_hash[bin]++;
            gsl_vector_int_set(dist_bin,locVec[i].snpVec[j].index, bin);
        }
    }


    dtss_map[0] =0;
    dtss_rmap[0] = 0;
    int count = 1;
    for (map<int,int>::iterator it=bin_hash.begin(); it!=bin_hash.end(); ++it){
        //printf("bin %d    %d \n", it->first, it->second);
        if(it->first==0)
            continue;
        dtss_map[it->first] = count;
        dtss_rmap[count] = it->first;
        //printf("%d   %f    %d\n",count, map_bin_2_dist(count,dist_bin_size),it->second);
        count++;
    }

    dist_bin_level = count;

    for(int i=0;i<p;i++){
        //printf("%d \n", gsl_vector_int_get(dist_bin,i));
        gsl_vector_int_set(dist_bin, i, dtss_map[gsl_vector_int_get(dist_bin,i)]);
    }
    //exit(1);

}





void controller::load_annotation(char* annot_file){


    map<string, vector<double> > annot_map;

    map<int, int> col2cat;
    map<int, int> col2cpos;
    map<int, int> col2dpos;
    int col_count = -1;


    if(strlen(annot_file)>0){




        ifstream afile(annot_file, ios_base::in | ios_base::binary);
        boost::iostreams::filtering_istream in;

        in.push(boost::iostreams::gzip_decompressor());
        in.push(afile);

        string line;
        istringstream ins;

        // parsing headers

        int count_d=0;
        int count_c=0;

        fprintf(stderr, "Loading annotations ... \n");


        while(getline(in,line)){
            ins.clear();
            ins.str(line);

            string token;

            while(ins>>token){

                if(col_count==-1){
                    if(token == "SNP" || token == "snp"){
                        col_count = 0;
                        continue;
                    }else{
                        break;
                    }
                }

                string cat = token.substr(token.size()-2, 2);
                // continuous
                if(cat == "_c" || cat =="_C"){
                    col2cat[col_count] = 1;
                    col2cpos[col_count] = kc;
                    string name = token.substr(0,token.size()-2);
                    cvar_name_vec.push_back(name);
                    kc++;
                } // discrete/categorical
                else{
                    col2cat[col_count] = 2;
                    col2dpos[col_count] = kd;
                    string name = token.substr(0,token.size()-2);
                    dvar_name_vec.push_back(name);
                    kd++;

                }

                col_count++;
            }


            if(col_count!=-1)
                break;

        }


        // read in data


        while(getline(in,line)){

            ins.clear();
            ins.str(line);
            string snp;


            ins>>snp;

            if(snp_hash[snp] != 100)
                continue;

            double val;
            vector<double> avec;
            while(ins>>val){
                avec.push_back(val);	
            }
            annot_map[snp] = avec;
        }
    }




    // memory allocation

    if(kc>0){
        Xc = gsl_matrix_calloc(p,kc);
    }


    int ncol = kd;
    if(dist_bin_level>0)
        ncol++;
    if(ncol>0){
        Xd = gsl_matrix_int_calloc(p,ncol);
        dlevel = gsl_vector_int_calloc(ncol);
    }


    if(kc+kd>0){
        for (int i=0;i<locVec.size();i++){
            for(int j=0;j<locVec[i].snpVec.size();j++){
                string snp_id = locVec[i].snpVec[j].id;
                int index = locVec[i].snpVec[j].index;

                vector<double> avec((kd+kc),0.0);
                if(annot_map.find(snp_id)!=annot_map.end())
                    avec = annot_map[snp_id];



                for(int k=0;k<avec.size();k++){

                    if(col2cat[k] == 1){
                        gsl_matrix_set(Xc,index,col2cpos[k],avec[k]);
                    }else{
                        gsl_matrix_int_set(Xd, index, col2dpos[k],int(avec[k]));
                    }

                }

            }
        }

        for(int i=0;i<kd;i++){
            int nl = count_factor_level(i);

            gsl_vector_int_set(dlevel, i,nl);
        }
    }


    if(dist_bin_level>0){
        gsl_matrix_int_set_col(Xd, ncol-1, dist_bin);
        gsl_vector_int_set(dlevel, ncol-1, dist_bin_level);
        dvar_name_vec.push_back(string("dtss"));
        kd++;
    }



    if(dist_bin != 0)
        gsl_vector_int_free(dist_bin);


}



int controller::count_factor_level(int col){

    map<int, int> rcd;
    for(int i=0;i<p;i++){
        int val = gsl_matrix_int_get(Xd,i,col);
        rcd[val] = 1;
    }

    return rcd.size();
}


vector<double> controller::make_grid(double phi_min, double phi_max){

    vector<double> gvec;
    gvec.push_back(phi_max); //phi_max
    double phi = phi_max;
    while(phi>phi_min){
        phi = phi/sqrt(2);
        gvec.push_back(phi); //phi
    }
    std::sort(gvec.begin(),gvec.end());

    return gvec;
}



void controller::init_params(){

    prior_vec = gsl_vector_calloc(p);
    pip_vec = gsl_vector_calloc(p);
    for(int i=0;i<p;i++){
        gsl_vector_set(prior_vec, i, init_pi1);
    }


    for(int i=0;i<locVec.size();i++){
        locVec[i].pip_vec = pip_vec;
        locVec[i].prior_vec = prior_vec;
    }




    ncoef = 0;
    for(int i=0;i<kd;i++){
        ncoef += gsl_vector_int_get(dlevel,i)-1;
    }

    ncoef += 1;
    ncoef += kc;



    beta_vec = gsl_vector_calloc(ncoef);

}


void controller::run_EM(){  

    // print header      	
    fprintf(stderr,"Starting EM ... \n");
    int count = 1;
    double last_log10_lik = -9999999;
    init_params();

    fprintf(stderr,"  Iter          loglik          Intercept    ");

    for(int i=0;i<dvar_name_vec.size();i++){
        string prefix = dvar_name_vec[i];
        int level = gsl_vector_int_get(dlevel,i);
        for(int j=1;j<level;j++){
            ostringstream stream;
            stream <<prefix<<"."<<j;
            string label = stream.str();
            fprintf(stderr, "%s\t", label.c_str());
        }

    }

    for(int i=0;i<cvar_name_vec.size();i++){
        fprintf(stderr, "%s\t", cvar_name_vec[i].c_str());
    }



    if(grid_size> 1 && !use_ash){
        fprintf(stderr,"  grid_weights (%d)\t",grid_size);
    }

    fprintf(stderr,"\n");

    int first_run = 1;

    // start iteration
    while(1){


        double curr_log10_lik = 0;

        //E-step
        for(int i=0;i<locVec.size();i++){
            locVec[i].EM_update();
            curr_log10_lik += locVec[i].log10_lik;
        }


        // M-step
        if(grid_size > 1){
            update_grid_weight();
        }

        if(first_run && single_fuzzy_annot == 1){
            //set shrinkage parameter
            single_probt_est();
            first_run = 0;
        }

        if(ncoef==1){
            simple_regression();
        }
        // only categrical annotations
        else if(kc==0 && kd!=0){
            if(kd == 1 && !force_logistic){
                single_ct_regression();
            }else{
                logistic_cat_fit(beta_vec, Xd, dlevel,pip_vec, l1_lambda,l2_lambda);
            }

            logistic_cat_pred(beta_vec, Xd, dlevel,prior_vec);
        }
        // only continuous annotations
        else if(kc!=0 && kd==0){

            logistic_cont_fit(beta_vec, Xc, pip_vec,l1_lambda, l2_lambda);      
            logistic_cont_pred(beta_vec,Xc, prior_vec);

        }
        // both continuous and categorical annotations
        else if(kc!=0 && kd !=0){
            logistic_mixed_fit(beta_vec, Xd, dlevel, Xc, pip_vec, l1_lambda,l2_lambda);
            logistic_mixed_pred(beta_vec,Xd, dlevel, Xc, prior_vec);
        }


        // output update to stderr

        fprintf(stderr,"%4d        %10.3f        ",count++,curr_log10_lik/log10(exp(1)));


        for(int i=0;i<ncoef;i++){
            fprintf(stderr, "%9.3f  ",gsl_vector_get(beta_vec,i));
        }


        if(grid_size>1 && !use_ash){
            for(int i=0;i<grid_size;i++){
                fprintf(stderr, "%9.3f  ",grid_weight[i]);
            }
        }

        fprintf(stderr,"\n");

        // convergence checking


        if(fabs(curr_log10_lik-last_log10_lik)<EM_thresh){
            final_log10_lik = curr_log10_lik;
            break;
        }

        last_log10_lik = curr_log10_lik;
    }


    double tp = 0;
    for(int i=0;i<p;i++){
        tp += gsl_vector_get(prior_vec,i);
    }
    if(print_avg)
        fprintf(stderr, "\nAvg. probability of association (%d unique variants): %9.3e\n", int(snp_hash.size()),tp/int(snp_hash.size()));


}




void controller::update_grid_weight(){

    if(grid_size<2)
        return;
    vector<double> pcount_vec(grid_size+1, 0);
    for(int i=0;i<locVec.size();i++){
        for(int j=0; j<grid_size+1;j++){
            // avoid numerical underflow

            pcount_vec[j] += locVec[i].wv_probv[j];
        }
    }



    double sum = 0;
    double w0 = pcount_vec[grid_size]/double(locVec.size());
    for(int i=0;i<grid_size;i++){
        grid_weight[i] = (pcount_vec[i]/double(locVec.size()))/(1-w0);
        if(grid_weight[i]<0)
            grid_weight[i] = 0;
        sum += grid_weight[i];
    }

    for(int i=0;i<grid_size;i++){
        grid_weight[i] = grid_weight[i]/sum;
    }


}



void controller::simple_regression(){

    double sum = 0;
    for(int i=0;i<p;i++){
        sum += gsl_vector_get(pip_vec,i);
    }
    double new_prior = sum/p;

    if(shrink_pi0 != 0){
        double pi0 = pi0_ubd*shrink_pi0 + (1-new_prior)*(1-shrink_pi0);
        new_prior = 1-pi0;
    }     

    for(int i=0;i<p;i++){
        gsl_vector_set(prior_vec,i,new_prior);
    }

    gsl_vector_set(beta_vec,0, log(new_prior/(1-new_prior)));
}



void controller::single_probt_regression(){
    /*
       for(int i=0;i<p;i++){

       double qv = gsl_matrix_get(Xc,i,0);
       double pv = gsl_vector_get(pip_vec,i);
       printf("%7.3e  %7.3e\n",pv, qv);
       }
       exit(0);
       */
    double b0 = gsl_vector_get(beta_vec,0);
    double b1 = gsl_vector_get(beta_vec,1);

    //gsl_vector_set(beta_vec,1, b1);
    double p1 = exp(b0+b1)/(1+exp(b0+b1));
    double p0 = exp(b0)/(1+exp(b0));


    for(int i=0;i<p;i++){
        double qv = gsl_matrix_get(Xc,i,0);
        //double new_prior = qv*p1 + (1-qv)*p0;
        double new_prior = exp(b0+b1*qv)/(1+exp(b0+b1*qv));
        gsl_vector_set(prior_vec,i,new_prior);   
    }

}



void controller::single_probt_est(){

    if(!(kc==0&&kd==1) && !(kd==0&&kc==1))
        return;

    double c00 = 0;
    double c11 = 0;
    double c10 = 0;
    double c01 = 0;

    for(int i=0;i<p;i++){
        double qv;
        if(kc==1)
            qv = gsl_matrix_get(Xc,i,0);
        if(kd==1)
            qv = double(gsl_matrix_int_get(Xd,i,0));

        double pv = gsl_vector_get(pip_vec,i);
        c11 += pv*qv;
        c10 += pv*(1-qv);
        c01 += (1-pv)*qv;
        c00 += (1-pv)*(1-qv);
    }

    double sd = sqrt( 1/c11 + 1/c10 + 1/c01 + 1/c00);
    if(l2_lambda==0){
        l2_lambda = sd;
        fprintf(stderr, "Applying adaptive L2 penalty = %.3f\n\n", l2_lambda);
        //fprintf(stderr, "Set L2 penalty = %.3f\n", sd);
        force_logistic = 1;
    }
    //fprintf (stderr, "CT sd = %.3f  %.1f \n",sd, c11);
}



void controller::single_ct_regression(){

    map<int,double> sum_pip;
    map<int,double> sum;

    int levels = gsl_vector_int_get(dlevel,0);

    for(int i=0;i<levels;i++){
        sum_pip[i] = sum[i] = 0;
    }

    for(int i=0;i<p;i++){
        int cat = gsl_matrix_int_get(Xd,i,0);
        sum_pip[cat] += gsl_vector_get(pip_vec,i);
        sum[cat] += 1;
    }


    for(int i=0;i<p;i++){
        int cat = gsl_matrix_int_get(Xd,i,0);
        gsl_vector_set(prior_vec,i,sum_pip[cat]/sum[cat]);
    }


    double baseline=0;
    for(int i=0;i<levels;i++){
        double new_prior = sum_pip[i]/sum[i];
        gsl_vector_set(beta_vec,i, log(new_prior/(1-new_prior))-baseline);
        if(i==0){
            baseline = log(new_prior/(1-new_prior));
        }

    }


}




// option 1: find egene


void controller::find_eGene(double fdr_thresh){

    if(!finish_em){
        run_EM();
        finish_em = 1;
    }

    for(int i=0;i<locVec.size();i++){
        locVec[i].compute_fdr();
    }

    std::sort(locVec.begin(), locVec.end(),rank_by_fdr);

    double rej = 1.0;
    double cpr = 0.0;
    int rej_decision = 1;
    int rej_count = 0;
    for(int i=0;i<locVec.size();i++){
        cpr += locVec[i].fdr;
        if(cpr/rej > fdr_thresh){
            rej_decision = 0;
        }
        printf("%5d  %20s    %9.3e    %d",int(rej), locVec[i].id.c_str(), locVec[i].fdr,rej_decision);
        rej++;
        if(rej_decision == 1){
            rej_count++;
        }

        if(print_lead_snp){
            std::sort(locVec[i].snpVec.begin(), locVec[i].snpVec.end(), sort_by_BF);
            printf("\t%20s", locVec[i].snpVec[0].id.c_str());
        }
        printf("\n");
    }

    fprintf(stderr,"\n\n    Total Loci: %d   Rejections: %d\n", int(locVec.size()), rej_count); 

}


// option 2: parameter estimation

void controller::estimate(){

    if(!finish_em){
        run_EM();
        finish_em = 1;
    }

    gsl_vector *saved_beta_vec = gsl_vector_calloc(ncoef);
    gsl_vector *saved_prior_vec = gsl_vector_calloc(p);
    gsl_vector_memcpy(saved_beta_vec, beta_vec);  
    gsl_vector_memcpy(saved_prior_vec, prior_vec);



    // CI for the intercept term
    double est = gsl_vector_get(beta_vec, 0);
    gsl_vector_set(beta_vec,0, 0.0);

    if(kc==0 && kd!=0){
        logistic_cat_pred(beta_vec, Xd, dlevel,prior_vec);
    }
    if(kc!=0 && kd==0){
        logistic_cont_pred(beta_vec,Xc, prior_vec);
    }
    if(kc!=0 && kd !=0){
        logistic_mixed_pred(beta_vec,Xd, dlevel, Xc, prior_vec);
    }

    if(kc==0 && kd ==0){

        for(int i=0;i<p;i++){
            gsl_vector_set(prior_vec,i,0.50);
        }
    }

    double null_log10_lik = 0;
    for(int k=0;k<locVec.size();k++){

        locVec[k].EM_update();
        null_log10_lik += locVec[k].log10_lik;

    }

    double diff = (final_log10_lik - null_log10_lik)/log10(exp(1));
    if(diff<1e-8){
        diff = 1e-8;
    }

    double sd = fabs(est)/sqrt(2*diff);
    printf("\n%25s  %9.3f     %9.3f  %9.3f\n", "Intercept", est, est-1.96*sd, est+1.96*sd);
    gsl_vector_memcpy(beta_vec, saved_beta_vec);
    gsl_vector_memcpy(prior_vec, saved_prior_vec);





    int index = 1; //start with 1, 0 always intercept



    for(int i=0;i<dvar_name_vec.size();i++){

        int level = gsl_vector_int_get(dlevel,i);
        string prefix = dvar_name_vec[i];
        for(int j=1;j<level;j++){
            ostringstream stream;
            if(prefix=="dtss"){
                int bin = dtss_rmap[j];
                stream<<"dtss."<<j<<"::["<<map_bin_2_dist(bin,dist_bin_size)<<"kb]";
            }else
                stream <<prefix<<"."<<j;
            string label = stream.str();
            double est = gsl_vector_get(beta_vec, index);
            gsl_vector_set(beta_vec,index, 0.0);

            if(kc==0 && kd!=0){
                logistic_cat_pred(beta_vec, Xd, dlevel,prior_vec);
            }
            if(kc!=0 && kd==0){
                logistic_cont_pred(beta_vec,Xc, prior_vec);
            }
            if(kc!=0 && kd !=0){
                logistic_mixed_pred(beta_vec,Xd, dlevel, Xc, prior_vec);
            }

            double null_log10_lik = 0;

            for(int k=0;k<locVec.size();k++){

                locVec[k].EM_update();
                null_log10_lik += locVec[k].log10_lik;

            }

            double diff = (final_log10_lik - null_log10_lik)/log10(exp(1));

            //printf("%f  %f  %f\n", null_log10_lik/log10(exp(1)), final_log10_lik/log10(exp(1)), diff);

            if(diff<0){
                double curr_log10_lik = final_log10_lik;
                est = fine_optimize_beta(index, est, null_log10_lik,  curr_log10_lik);
                diff = fabs(curr_log10_lik - null_log10_lik)/log10(exp(1));
            }

            if(diff<1e-8){
                diff = 1e-8;
            }
            double sd = fabs(est)/sqrt(2*diff);
            printf("%25s  %9.3f     %9.3f  %9.3f\n", label.c_str(), est, est-1.96*sd, est+1.96*sd);
            index++;
            // restore
            gsl_vector_memcpy(beta_vec, saved_beta_vec);  
            gsl_vector_memcpy(prior_vec, saved_prior_vec);
        }
    }


    for(int i=0;i<cvar_name_vec.size();i++){
        string label =  cvar_name_vec[i];
        double est = gsl_vector_get(beta_vec, index);
        gsl_vector_set(beta_vec,index, 0.0);


        if(kc==0 && kd!=0){
            logistic_cat_pred(beta_vec, Xd, dlevel,prior_vec);
        }
        if(kc!=0 && kd==0){
            logistic_cont_pred(beta_vec,Xc, prior_vec);
        }
        if(kc!=0 && kd !=0){
            logistic_mixed_pred(beta_vec,Xd, dlevel, Xc, prior_vec);
        }

        double null_log10_lik = 0;

        for(int k=0;k<locVec.size();k++){

            locVec[k].EM_update();
            null_log10_lik += locVec[k].log10_lik;

        }

        double diff = (final_log10_lik - null_log10_lik)/log10(exp(1));

        if(diff<0){
            double curr_log10_lik = final_log10_lik;
            est = fine_optimize_beta(index, est, null_log10_lik,  curr_log10_lik);
            diff = fabs(curr_log10_lik - null_log10_lik)/log10(exp(1));
        }

        if(diff<1e-8){
            diff = 1e-8;
        }
        double sd = fabs(est)/sqrt(2*diff);
        printf("%25s  %9.3f     %9.3f  %9.3f\n", label.c_str(), est, est-1.96*sd, est+1.96*sd);
        index++;
        // restore
        gsl_vector_memcpy(beta_vec, saved_beta_vec);
        gsl_vector_memcpy(prior_vec, saved_prior_vec);
    }


    // restore all
    for(int k=0;k<locVec.size();k++){
        locVec[k].EM_update();
    }

    gsl_vector_free(saved_beta_vec);
    gsl_vector_free(saved_prior_vec);

    if(grid_size>1){
        printf("\n\n");
        if(!use_ash){ 
            for(int i=0;i<grid_size;i++){
                printf("%25s %3d   %7.3e\n","grid_weight", (i+1), grid_weight[i]);
            }	
        }else{
            printf("%28s  %9s\n", "Effect size (phi)", "Weight");	    
            for(int i = 0;i<grid_size;i++){
                printf("%15s %9.3e       %9.3e\n", "",grid_vec[i], grid_weight[i]);
            }	
        }

        if(set_pi0 > 0){
            printf("\nfixed_pi0  %.3f  log10_likelihood  %9.3f\n", set_pi0, final_log10_lik);
        }


    }
}


void controller::dump_prior(char *prior_path){

    if(!finish_em){
        run_EM();
        finish_em = 1;
    }


    fprintf(stderr,"Output priors to directory %s...\n", prior_path);
    if(mkdir(prior_path, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH)!=0){
        fprintf(stderr, "Error: cannot create directory %s\n", prior_path);
        exit(1);
    }

    string dir = string(prior_path);
    for( int i=0; i<locVec.size();i++){
        string file_name = dir + "/" + locVec[i].id+".prior";
        FILE *fd = fopen(file_name.c_str(), "w");
        for(int j=0;j<locVec[i].snpVec.size();j++){
            int index = locVec[i].snpVec[j].index;
            string name = locVec[i].snpVec[j].id;
            fprintf(fd, "%s  %9.4e\n",name.c_str(), gsl_vector_get(prior_vec, index));
        }
        fclose(fd);
    }
}


void controller::dump_pip(char *file){
    if(!finish_em){
        run_EM();
        finish_em = 1;
    }
    fprintf(stderr,"Output PIPs to file %s: ...\n",file);
    FILE *fd = fopen(file, "w");
    for( int i=0; i<locVec.size();i++){
        string gname = locVec[i].id;
        for(int j=0;j<locVec[i].snpVec.size();j++){
            int index = locVec[i].snpVec[j].index;
            string name = locVec[i].snpVec[j].id;
            fprintf(fd, "%s\t%s\t%9.4e\t%9.4e\n",name.c_str(), gname.c_str(),gsl_vector_get(pip_vec, index) , gsl_vector_get(prior_vec, index));
        }
    }
    fclose(fd);
}



void controller::dump_wv(char *file){
    if(grid_size<2)
        return;

    if(!finish_em){
        run_EM();
        finish_em = 1;
    }

    fprintf(stderr,"Output grid weights to file %s: ...\n",file);
    FILE *fd = fopen(file, "w");
    for( int i=0; i<locVec.size();i++){
        locVec[i].compute_fdr();
        fprintf(fd, "%s\t%d\t\t%9.3e\t", locVec[i].id.c_str(), int(locVec[i].snpVec.size()), 1-locVec[i].fdr);
        for(int j=0;j<grid_size;j++){
            fprintf(fd, "\t%9.3e", locVec[i].wv_probv[j]);
        }
        fprintf(fd,"\n");
    }
}







double controller::fine_optimize_beta(int index, double est, double null_log10_lik, double &curr_log10_lik){

    double gr = (sqrt(5)-1)/2.0;

    /*
       double a = -fabs(est);
       double b = fabs(est);
       */

    double a = est;
    double b = 0;
    if(a > b){
        b = a;
        a = 0;
    }
    double c = b - gr*(b-a);
    double d = a + gr*(b-a);

    double thresh = 1e-3;
    if(fabs(c-d)<thresh){
        thresh = fabs(c-d)/2.0;
    }

    double fc;
    double fd;


    while((d-c)>thresh){

        fc = eval_likelihood(c, index);
        fd = eval_likelihood(d, index);


        //printf("%f %f %f %f   %f %f   %f\n",a, d
        if(fc > fd){
            /*
               if(fc > null_log10_lik){
               curr_log10_lik = fc;
               return c;
               }
               */
            b = d;
            d = c;
            c = b - gr*(b-a);
        }else{
            /*
               if(fd > null_log10_lik){
               curr_log10_lik = fd;
               return d;
               }
               */
            a = c;
            c = d;
            d = a + gr*(b-a);
        }
    }


    curr_log10_lik = fc;

    return (b+a)/2;



}



double controller::eval_likelihood(double x, int index){


    gsl_vector_set(beta_vec,index, x);

    if(kc==0 && kd!=0){
        logistic_cat_pred(beta_vec, Xd, dlevel,prior_vec);
    }
    if(kc!=0 && kd==0){
        logistic_cont_pred(beta_vec,Xc, prior_vec);
    }
    if(kc!=0 && kd !=0){
        logistic_mixed_pred(beta_vec,Xd, dlevel, Xc, prior_vec);
    }

    double log10_lik = 0;
    for(int k=0;k<locVec.size();k++){
        locVec[k].EM_update();
        log10_lik += locVec[k].log10_lik;
    }


    return log10_lik;
}





void Locus::get_BF_matrix(){
    if(snpVec[0].grid_size<2)
        return;

    for(int i = 0; i< snpVec[0].grid_size;i++){
        vector<double> BFv;
        for(int j=0;j<snpVec.size();j++){
            BFv.push_back(snpVec[j].log10_BF_vec[i]);
        }
        BF_wv_vec.push_back(BFv);
    }


}


void Locus::EM_update(){

    // compute log10_lik
    vector<double> BF_vec;  // snp-level avg. BF
    vector<double> p_vec;


    double locus_pi0 = 1;
    for(int i=0;i<snpVec.size();i++){
        double prior = gsl_vector_get(prior_vec, snpVec[i].index);
        snpVec[i].log10_BF = snpVec[i].compute_log10_BF();


        BF_vec.push_back(snpVec[i].log10_BF);
        p_vec.push_back(prior/(1-prior));
        locus_pi0 *= (1-prior);
    }

    if(locus_pi0 < 1e-100){
        locus_pi0 = 1e-100;
    }

    for(int i=0;i<snpVec.size();i++){
        p_vec[i] *= locus_pi0;
    }


    double baseline_val = 0;
    if(baseline_vec.size() > 0 ){
        baseline_val = log10_weighted_sum(baseline_vec, snpVec[0].weightv);
    }


    BF_vec.push_back(baseline_val);
    vector<double> snp_wv = p_vec; // wihtout null config prob
    p_vec.push_back(locus_pi0); // add null config prob

    log10_lik = log10_weighted_sum(BF_vec, p_vec);

    for(int i=0;i<snpVec.size();i++){
        double val = log10(p_vec[i]) + snpVec[i].log10_BF - log10_lik;
        double pip = pow(10, val);
        gsl_vector_set(pip_vec, snpVec[i].index, pip);
    }

    if(snpVec[0].grid_size>1){
        wv_probv.clear();
        double *weightv = snpVec[0].weightv;
        double sum = 0;
        for(int i=0;i<snpVec[0].grid_size;i++){
            vector<double> BFv = BF_wv_vec[i];
            double v = log10_weighted_sum(BFv, snp_wv);
            double w_val = pow(10, v-log10_lik)*weightv[i];
            if(w_val<0)
                w_val = 0;
            wv_probv.push_back(w_val);
            sum+= w_val;
        }
        //printf("%7.3f ", sum);
        double pip0 = pow(10, log10(locus_pi0) - log10_lik);
        //printf("%7.3f  %7.3f \n", pip0, pip0+sum);
        wv_probv.push_back(pip0);
    }

}







void Locus::compute_fdr(){

    double locus_pi0 = 1;

    // compute log10_lik
    vector<double> BF_vec;
    vector<double> p_vec;
    for(int i=0;i<snpVec.size();i++){
        double prior = gsl_vector_get(prior_vec, snpVec[i].index);
        BF_vec.push_back(snpVec[i].log10_BF);
        p_vec.push_back(prior/(1-prior));
        locus_pi0 *= (1-prior);
    }

    for(int i=0;i<snpVec.size();i++){
        p_vec[i] *= locus_pi0;
    }


    double baseline_val =0;
    if(baseline_vec.size() > 0 ){
        baseline_val = log10_weighted_sum(baseline_vec, snpVec[0].weightv);
    }

    BF_vec.push_back(baseline_val);
    p_vec.push_back(locus_pi0);

    log10_lik = log10_weighted_sum(BF_vec, p_vec);
    fdr =  pow(10,log10(locus_pi0)+baseline_val - log10_lik);

}


double SNP::compute_log10_BF(){

    if(grid_size==1){
        return log10_BF;
    }

    return log10_weighted_sum(log10_BF_vec, weightv);
}

void SNP::compute_log10_BF_vec_bhat(vector<double> &grid){

    if(se_beta == 0){
        log10_BF_vec = vector<double> (grid.size(),0);
        return;
    }   

    double z2 = pow((beta/se_beta), 2.0);
    double v2 = pow(se_beta, 2.0);
    vector<double> rstv;
    for(int i=0;i<grid.size();i++){
        double w2 = pow(grid[i],2.0);
        double val = 0.5*log(v2/(v2+w2)) + 0.5*z2*(w2/(v2+w2));
        log10_BF_vec.push_back(val/log(10));
    }

    grid_size = grid.size();

}



void SNP::compute_log10_BF_vec_zval(vector<double> &grid){

    double z2 = pow(zval, 2.0);
    vector<double> rstv;
    grid_size = grid.size();
    for(int i=0;i<grid_size;i++){
        double K = pow(grid[i],2);	  
        double val = -0.5*log(1+K) + 0.5*z2*K/(1+K);
        log10_BF_vec.push_back(val/log(10));
    }
}

double compute_log10_BF(double beta, double se_beta){

    if(se_beta == 0)
        return 0;
    //double phi[4] = {0.2,0.4,0.8, 1.6};
    double phi[4] = {1,4,16,25};
    int size = 4;
    double z2 = pow((beta/se_beta), 2.0);
    double v2 = pow(se_beta, 2.0);
    vector<double> rstv;
    vector<double> wtv(size,1.0/double(size));
    for(int i=0;i<size;i++){
        //double w2 = pow(se_beta*phi[i],2.0);
        double w2 = pow(phi[i],2.0);
        double val = 0.5*log(v2/(v2+w2)) + 0.5*z2*(w2/(v2+w2));
        rstv.push_back(val/log(10));
    }

    return log10_weighted_sum(rstv,wtv);
}


double compute_log10_BF(double z_score){

    double kv[4] = {1,4,16,25};
    int size = 4;
    double z2 = pow(z_score, 2.0);
    vector<double> rstv;
    vector<double> wtv(size,1.0/double(size));
    for(int i=0;i<size;i++){
        double val = 0.5*log(1/(1+kv[i])) + 0.5*z2*(kv[i]/(1+kv[i]));
        rstv.push_back(val/log(10));
    }

    return log10_weighted_sum(rstv,wtv);
}




double log10_weighted_sum(vector<double> &vec, vector<double> &wts){


    double max = vec[0];
    for(size_t i=0;i<vec.size();i++){
        if(vec[i]>max)
            max = vec[i];
    }
    double sum = 0;
    for(size_t i=0;i<vec.size();i++){
        sum += wts[i]*pow(10, (vec[i]-max));
    }

    return (max+log10(sum));
}


double log10_weighted_sum(vector<double> &vec, double *wts){


    double max = vec[0];
    for(size_t i=0;i<vec.size();i++){
        if(vec[i]>max)
            max = vec[i];
    }
    double sum = 0;
    for(size_t i=0;i<vec.size();i++){
        sum += wts[i]*pow(10, (vec[i]-max));
    }

    return (max+log10(sum));
}



bool rank_by_fdr (const Locus & loc1 , const Locus & loc2){
    return (loc1.fdr < loc2.fdr);
}


bool sort_by_BF(const SNP &snp1, const SNP &snp2){
    return snp1.log10_BF > snp2.log10_BF;
}





// default binning scheme
int classify_dist_bin(int pos, int tss, double bin_size){

    if(bin_size>0){
        return int((pos-tss)/bin_size);
    }
    // else default scheme
    double dist = (pos-tss)/1000.0;
    int sign = 1;
    if(dist < 0){
        sign = -1;
    }
    dist = fabs(dist);
    int bin;

    if(dist<0.5){
        bin = 0;
    }else if(dist<1){
        bin = 1;
    }else if(dist<2.5){
        bin = 2;
    }else if(dist<5){
        bin = 3;
    }else if(dist<10){
        bin = 4;
    }else if(dist <25){
        bin = 5;
    }else if(dist <50){
        bin = 6;
    }else if(dist <100){
        bin = 7;
    }else if(dist <250){
        bin =8;
    }else if(dist < 500){
        bin = 9;
    }else{
        bin = 10;
    }

    return bin*sign;

}



double map_bin_2_dist(int bin, double bin_size){

    // return the mid point of the bin
    if(bin_size>0){
        if(bin == 0){
            return 0;
        }
        if(bin>0){
            return 0.5*(bin+1)*bin_size/1000;
        }
        if(bin<0){
            return 0.5*(bin-1)*bin_size/1000;
        }
    }


    // else default scheme
    int sign = 1;
    if(bin < 0){
        sign = -1;
    }

    if(bin == 0){
        return 0;
    }
    if(abs(bin)==1){
        return sign*0.75;
    }
    if(abs(bin)==2){
        return sign*1.75;
    }
    if(abs(bin)==3){
        return sign*3.75;
    }
    if(abs(bin)==4){
        return sign*7.5;
    }
    if(abs(bin)==5){
        return sign*17.5;
    }
    if(abs(bin)==6){
        return sign*37.5;
    }
    if(abs(bin)==7){
        return sign*75;
    }
    if(abs(bin)==8){
        return sign*175;
    }
    if(abs(bin)==9){
        return sign*375;
    }

    return sign*750;


}
