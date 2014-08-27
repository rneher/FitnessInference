
/**
 * @file flusim.cpp
 * @brief simulate evolution of influenza with time resolved sampling and recorded genealogies
 * @author Richard Neher
 * @version 
 * @date 2013-12-22
 */

//g++ flusim.cpp -o flu -I/ebio/ag-neher/share/users/rneher/FFPopSim/ffpopsim/pkg/include -L/ebio/ag-neher/share/users/rneher/FFPopSim/ffpopsim/pkg/lib -lgsl -lgslcblas -lFFPopSim

//g++ flusim.cpp -o flu -I/home/richard/Projects/ffpopsim/pkg/include -L/home/richard/Projects/ffpopsim/pkg/lib -lgsl -lgslcblas -lFFPopSim

/* Include directives */
#include "ffpopsim_highd.h"
#include "multi_population.h"
#include <fstream>
#include <sstream>
#include <string>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>


// set up flu
int flu(int L, int N, int sample_size, int sampling_dt, double target_sigma, double mu, 
	double gamma_bene, double s_bene, double gamma, double s, double nflip, int flip_dt) {
    int gen_max = 20000, burn_in = 5000;
    gsl_rng* flurng = gsl_rng_alloc(RNG);

    // construct a directory to output the data
    stringstream fname;
    string today = "20140820";
    fname <<"../data_new/N_"<<N<<"_L_"<<L<<"_nflip_"<<nflip<<"_mu_"<<mu<<"_sdt_"<<sampling_dt<<"/"<<today;
    ofstream params((fname.str()+"_params.dat").c_str());
    ofstream stats((fname.str()+"_stats.dat").c_str());
    cout <<"opening file: "<<fname.str()+"_params.dat"<<endl;
    params <<"population size:\t"<<N<<"\nnumber of loci:\t"<<L
	   <<"\nmutation rate:\t"<<mu<<"\nbenefical gamma:\t"<<gamma_bene
	   <<"\ndeleterious gamma:\t"<<gamma<<"\neffect ratio(bene/dele):\t"<<s_bene/s
	   <<"\nbeneficial per gen:\t"<<nflip<<"\nfitness std:\t"<<target_sigma<<"\n";
    params.close();

    // set up the population and specify parameters
    haploid_highd pop(L);
    pop.set_mutation_rate(mu);
    // this is flu, hence no recombination
    pop.outcrossing_rate = 0;
    pop.crossover_rate = 0;
    pop.recombination_model = CROSSOVERS;
    vector <int> gen_loci;
    gen_loci.push_back(L/2);
    pop.track_locus_genealogy(gen_loci);

    vector <int> loci;
    int fitness_coeff_sign;
    for(int i=0; i< L; i++) {
	loci.assign(1, i);
	fitness_coeff_sign = -1;  // all initital fitness effects are negative, hence the population starts in an adapted state
	pop.add_fitness_coefficient(fitness_coeff_sign*gsl_ran_gamma(flurng, gamma, s), loci);		
	loci.clear();
    }

    pop.set_wildtype(N);		// start with a population of the right size in the adapted state

    stat_t fitstat;
    stat_t fitness_stat;
    double new_trait_weights[1];
    double pair_dist = 0;
    int actual_flips;
    // this is the initial sample to be used as outgroup in the analysis
    pop.tree_sample=1;
    for (int gi=0; gi<burn_in+gen_max; gi++){
	// flip loci, not every generation but every flipdt, this makes things 
	// faster as fitness doesn't need to be recomputed
	if (gi%flip_dt==0){
	    // number of flipped loci on average nflip*flip_dt, hence nflip per generation
	    actual_flips = gsl_ran_poisson(flurng, nflip*flip_dt);
	    for (int ni=0; ni<actual_flips; ni++){
		// flip in the first 500 position, mutations in the rest remain deleterious 
		loci.assign(1,gsl_rng_uniform_int(flurng, 500));
		fitness_coeff_sign = (pop.trait[0].get_additive_coefficient(loci[0])>0)?-1:1;
		pop.add_fitness_coefficient(fitness_coeff_sign*gsl_ran_gamma(flurng, gamma_bene, s_bene), loci);
		loci.clear();
	    }
	    // recalculate traits and fitness
	    if (actual_flips>0){
		pop.update_traits();
		pop.update_fitness();
	    }
	}
	pop.evolve(1);
	// sample only if post burn_in, outgroup has been written in the first iteration
	if (gi<burn_in)	pop.tree_sample=0;
	// sample at intervals sampling_dt 
	else if (gi%sampling_dt==0) pop.tree_sample=sample_size;
	else pop.tree_sample=0;

	// output of population statistics
	fitness_stat = pop.get_fitness_statistics();
	cout <<"gen: "<<gi<<" fitness std "<<sqrt(fitness_stat.variance)<<endl;
	if (gi%10 ==0){
	    pair_dist=0;
	    //calculate average pairwise distance between strains
	    for (int locus=0; locus<L; locus++){
		pair_dist+=pop.get_allele_frequency(locus)*(1-pop.get_allele_frequency(locus));
	    }
	    pair_dist*=2;
	    //write fitness std and pairwise distance to file
	    stats <<gi<<"\t"<<sqrt(fitness_stat.variance)<<"\t"<<pair_dist<<endl;
	}
	//rescale the trait weights to such that fitness_std deviation matches it
	new_trait_weights[0] = pop.get_trait_weight(0)*target_sigma/sqrt(fitness_stat.variance);
	// only do rescaling after half the burnin to allow diversity to build up before dividing by it
	if (gi>burn_in/2) pop.set_trait_weights(new_trait_weights);
    }

    // dump the tree of all sampled sequences to file 
    ofstream flu_tree((fname.str()+"_tree.nwk").c_str());
    rooted_tree sample_tree;
    pop.genealogy.trees[0].construct_subtree(pop.genealogy.trees[0].sampled_leafs, sample_tree);
    flu_tree <<sample_tree.print_newick()<<endl;
    flu_tree.close();

    // dump all sampled sequences to file
    ofstream flu_seqs((fname.str()+"_seqs.fa").c_str());
    flu_seqs <<pop.genealogy.trees[0].print_sequences()<<endl;
    flu_seqs.close();
    stats.close();
    return 0;
}


/* MAIN */
int main(int argc, char **argv){

    int status= 0;
    if (argc == 9) {
	int L = atoi(argv[1]);
	int N = atoi(argv[2]);
	int ssize = atoi(argv[3]);    // sample size
	int sdt = atoi(argv[4]);      // interval between samples
	double sigma = atof(argv[5]); // fitness standard deviation
	double nflip = atof(argv[6]); // number of effect reversals (deleterious mutation becoming beneficial) per generation
	double mu = atof(argv[7]);    //mutation rate
	double bene_dele_ratio = atof(argv[8]); // ratio of effects of beneficial and deleterious mutations
	status += flu(L,N,ssize,sdt,sigma, mu, 
		      2.0, 0.0025*bene_dele_ratio, 1.0, 0.0025, nflip, 1);
    }else{
	cout<<"Usage: "<<argv[0]<<endl;
	status = 1;
    } 
    return status;
}


