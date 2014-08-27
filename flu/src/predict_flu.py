#########################################################################################
#
# author: Richard Neher
# email: richard.neher@tuebingen.mpg.de
#
# Reference: Richard A. Neher, Colin A Russell, Boris I Shraiman. 
#            "Predicting evolution from the shape of genealogical trees"
#
#########################################################################################
#
# predict_flu.py
# This script provides influenza specific subclasses of the general fitness inference
# and prediction algorithm as well as a data set class of flu sequences. 
# 
# class flu_alignment inherits the alignement class from sequence_ranking.py
# and adds functionality to select different subsets of data based on the 
# dates and geographic regions
#
# class flu_ranking subclasses sequence_ranking and adds method to 
# attach an ad-hoc fitness shift to individual branches of the tree.
#


import sys
sys.path.append('../../prediction_src')
from sequence_ranking import alignment, sequence_ranking
from Bio import Phylo, AlignIO, Seq,SeqIO,SeqRecord
from Bio.Align import MultipleSeqAlignment, AlignInfo
from Bio.Alphabet import generic_dna, generic_protein
import pickle, gzip,time
import numpy as np
import tree_utils
from datetime import date
verbose = True

# antigenic sites in numbering starting at 0 (pos-1)
# positions as defined in Shih et al, PNAS 2007
HA1_antigenic_sites = [pos - 1  for pos in 
                       [122,124, 126, 131, 133, 135, 137,  142, 143, 144, 145, 146,     # epi A
                       155, 156, 157, 158, 159, 160, 163, 164, 186, 188,                # epi B
                        189, 190, 192, 193, 196, 197,                                   # epi Bcont
                       50, 53, 54, 275, 276, 278, 299, 307,                             # epi C
                       121, 172, 173, 174, 201, 207, 213, 217, 226, 227, 242, 244, 248]] #epi D

# cluster position in numbering starting at 0
cluster_positions = sorted([188,192,155,158,157, 154, 144])

def koel_mutation_branch_label(node):
    if len(node.get_terminals())>1:
        return branch_label(node, aa=True, display_positions = cluster_positions)
    else:
        return ''

def ofunc(fname,mode):
    '''
    unified function to open zipped and unzipped files
    parameters:
    fname   --  file name, if last two letters are gz, use gzip open
    mode    --  opening mode as in 'r', 'w'
    '''
    if fname[-2:]=='gz':
        return gzip.open(fname, mode)
    else:
        return open(fname, mode)


def coordinates_to_region(lat, lng):
    '''
    returns the region based on the geographic sectors defined by longitude and latitude
    argument:
    lat  -- latitude
    lng  -- longitude
    '''
    if lat>0:
        if lng<-20:
            return 'north america'
        elif lng>50:
            return 'asia'
        else:
            return 'europe'
    else:
        if lng<-20:
            return 'south america'
        else:
            return 'oceania'

##################################################
# ADDITIONAL RANKING FUNCTIONS#
# all additional ranking function have the mandatory arguments
# terminals and non_terminals
##################################################

def combined_ranking_internal(terminals, non_terminals):
    '''
    calculates a joint score of the mean_fitness and the estimated 
    expansion rate
    '''
    mean_fitness_std = np.std([node.mean_fitness for node in non_terminals])
    expansion_std = np.std([node.expansion_score for node in non_terminals])
    for node in non_terminals:
        node.expansion_fitness = node.mean_fitness/mean_fitness_std+\
                                 0.3*node.expansion_score/expansion_std        

def combined_ranking_external(terminals, non_terminals):
    '''
    calculates a joint score of the mean fitness and the time stamp of
    external nodes.
    '''
    mean_fitness_std = np.std([node.mean_fitness for node in terminals])
    time_std = np.std([node.date.toordinal() for node in terminals])
    for node in terminals:
        node.time_fitness = node.mean_fitness/mean_fitness_std+\
                                 node.date.toordinal()/time_std        


###########################################################
## flu_alignment class
##########################################################
class flu_alignment(alignment):
    '''
    influenza data set that selects sequences from an alignment based on an annotation file
    and selection criteria
    '''
    def __init__(self, aln_fname, outgroup, annotation, seq_names = None, criteria = None, 
                 cds=None, subsample_factor=1.0, collapse=False):
        '''
        parameters:
        aln_fname   --  filename with alignment
        outgroup    --  outgroup sequence  
        annotation  --  dictionary or panda DataFrame that holds spatio/temporal info
        cds         --  coding region
        '''
        self.aln_file_name = aln_fname
        self.annotation = annotation
        self.subsample_factor = subsample_factor
        if seq_names is not None:
            self.select_subset_from_names(seq_names)
        elif criteria is not None:
            self.select_subset(criteria)
        else:
            print "either names of criteria need to be given"

        # initialize the base class, allele frequencies and the tree will be calculated
        alignment.__init__(self,self.subset_aln, outgroup, cds, collapse=collapse)
        self.annotate_leafs()


    def select_subset_single(self, start_date, stop_date, regions=None, full_date=True):
        '''
        This function adds sequences to the target alignment based on a single date/region pair
        arguments:
        start_date  --  lower bound for date >=
        stop_date   --  upper bound for date <
        regions     --  list of regions as strings, optional
        full_date   --  require full date, default True
        '''
        seqs_to_keep = []
        # make all regions lower case to facilitate comparison
        if regions:
            tmp_regions = [region.lower() for region in regions]
        else:
            tmp_regions = []
        
        # loop over all terminal nodes and check dates and region
        for seqname, seq in self.sequence_lookup.iteritems():
            # check if sequences and annotation information is present
            # prune otherwise
            if seqname not in self.annotation:
                print 'no annotation', seqname
                continue
            if self.annotation[seqname]['date_info']=='no_date_info':
                print 'no date info', seqname
                continue

            # check date, prune otherwise
            if full_date:
                try:
                    date_ok = self.annotation[seqname]['date_info']=='full_date'
                except:
                    print 'no date info', seqname
                    date_ok = False
            else:
                date_ok=True

            if date_ok and start_date<=self.annotation[seqname]['date'] and \
                    stop_date>self.annotation[seqname]['date']:
                #passed date check, now check region
                in_region = True
                # check region, prune otherwise, if tmp_region list empty, go ahead
                if len(tmp_regions):
                    try:
                        region = self.annotation[seqname]['region'].lower()
                    except KeyError as e:
                        try:
                            region = coordinates_to_region(self.annotation[seqname]['lat'], 
                                                           self.annotation[seqname]['lng'])
                        except:
                            region = 'unknown'
                    if region not in tmp_regions:
                        in_region = False
                if in_region:
                    seqs_to_keep.append(self.sequence_lookup[seqname])
        return seqs_to_keep

    def select_subset_from_names(self, list_of_names):
        '''
        it creates an alignment with only the desired sequences
        parameters:
        list_of_names   --  list of names of all sequences to be kept
        '''
        if verbose:
            tmp_time = time.time()
            print "select_subset_and_subtree_from_names: Loading data..."

        with ofunc(self.aln_file_name, 'r') as infile:
            tmp_aln = AlignIO.read(infile, 'fasta')


        seqs_to_keep = set()
        self.subset_aln = MultipleSeqAlignment([])
        for seq in tmp_aln:
            if seq.name in list_of_names:
                self.subset_aln.append(seq)

        if verbose:
            print "done in",np.round(time.time()-tmp_time,2), 's'
            tmp_time = time.time()

    def select_subset(self, criteria):
        '''
        accepts a list of selection criteria and makes and alignment of all sequences
        that match at least one of the criteria
        parameters:
        criteria    --  list of criteria of form [start_date, stop_date, regions]
        '''
        if verbose:
            tmp_time = time.time()
            print "select_subset_and_subtree: Loading data..."

        with ofunc(self.aln_file_name, 'r') as infile:
            tmp_aln = AlignIO.read(infile, 'fasta')
        if verbose:
            print "loaded",len(tmp_aln),"sequences"
        self.sequence_lookup = {seq.name:seq for seq in tmp_aln}
        seqs_to_keep = set()

        # loop over different criteria 
        # produce a set of sequences to keep by union (present in at least one)
        for crit in criteria:
            start_date, stop_date, regions = crit[:3]
            tmp_seqs_to_keep = self.select_subset_single(start_date, stop_date, regions, full_date = True)
            print start_date,'to', stop_date,'in', regions,\
                ' found ', len(tmp_seqs_to_keep), 'leafs'
            if len(crit)>3 and crit[3]>0:
                self.subsample(tmp_seqs_to_keep, crit[3])
                print 'down-sampled to',len(tmp_seqs_to_keep)
            seqs_to_keep.update(tmp_seqs_to_keep)

        #print leafs_to_prune
        self.subset_aln = MultipleSeqAlignment([])
        for seq in seqs_to_keep:
            self.subset_aln.append(seq)
        if verbose:
            print "done in",np.round(time.time()-tmp_time,2), 's'
            tmp_time = time.time()


    def subsample(self, tmp_seqs_to_keep, nseqs):
        '''
        take lists, sample tmp_seqs down to size nseqs and add all
        not sampled nodes to the tmp_to_prune lists
        '''
        from random import sample as rd_sample
        tmp_nseqs = int(np.round(min(self.subsample_factor*len(tmp_seqs_to_keep), nseqs)))
        if len(tmp_seqs_to_keep)<=tmp_nseqs:
            return
        else:
            indices_to_pop = rd_sample(xrange(len(tmp_seqs_to_keep)), len(tmp_seqs_to_keep) - tmp_nseqs)
            indices_to_pop.sort(reverse=True)
            for ii in indices_to_pop:
                s_pop = tmp_seqs_to_keep.pop(ii)
   
    def annotate_leafs(self):
        '''
        annotate each leaf
        '''
        for node in self.T.get_terminals():
            if node.name in self.annotation:
                tree_utils.annotate_leaf(node, self.annotation[node.name])
            else:
                if verbose:
                    print "build_tree(): missing annotation for ",node.name

    def sampling_distribution(self, bins=20):
        all_dates = [node.date.toordinal() for node in self.T.get_terminals()]
        y,x = np.histogram(all_dates, bins=bins)
        return y,x,[date.fromordinal(int(b)) for b in x]


###########################################################################################
### flu prediction classs
##########################################################################################
        
class flu_ranking(sequence_ranking):
    def __init__(self, flu_data,
                 eps_branch_length=1e-5, boost = 0.0, time_bins=None, distance_scale = 3.0,
                 pseudo_count = 5, methods = ['mean_fitness', 'branch_length', 
                                              'depth', 'expansion_score'], D=0.5):
        sequence_ranking.__init__(self, flu_data, eps_branch_length = eps_branch_length, 
                                  time_bins=time_bins, pseudo_count = pseudo_count,
                                  methods= methods, D=D, distance_scale = distance_scale)
        self.display_positions =cluster_positions
        self.boost = boost
        if self.boost!=0.0:
            self.fitness_shift(pos = cluster_positions, dfit = self.boost)
    

    def fitness_shift(self, pos = [], dfit = 0.5):
        '''
        add mutations as labels to branches
        '''
        if self.data.protein:
            if len(pos)==0:
                pos = cluster_positions
            for node in self.terminals+self.non_terminals:
                node.fitness_shift = dfit*len([mut for mut in node.aa_mutations 
                                               if mut[0] in pos])
        else:
            print "No coding sequence provided --  can't calculate fitness shifts associated with amino acid substitutions"
            

######################
## quick test
######################
if __name__ == '__main__':
    flutype = 'H3N2'
    year = 2006
    prediction_regions = ['asia','north america']
    sample_size = 50
    aln_fname = '../data/'+flutype+'_HA1_all_years_filtered.fasta.gz'
    # load annotation
    with open('../data/'+flutype+'_annotations.pickle', 'r') as infile:
        annotation = pickle.load(infile)
    outgroup = SeqIO.read('../data/'+flutype+'_outgroup.fasta', 'fasta')

    # set up the filtering criteria and select sequences from the master alignment
    criteria=[(date(year-1, 5,1), date(year, 2,28), [region], sample_size)
            for region in prediction_regions]
    
    my_flu_alignment = flu_alignment(aln_fname, outgroup, annotation, criteria = criteria,
                                     cds = {'begin':0, 'end':987,'pad':0})

    # run the prediction
    prediction = flu_ranking(my_flu_alignment,boost = 0.5)
    top_seq = prediction.predict()
    print top_seq
    
    # plot the tree colored by the prediction
    tree_utils.plot_prediction_tree(prediction)

    # plot the distribution of sampling dates
    y,x,bin_names = my_flu_alignment.sampling_distribution(bins = 10)
    import matplotlib.pyplot as plt
    plt.figure()
    plt.plot(x[1:],y)
    plt.xticks(x[1:], map(str, [b for b in bin_names]), rotation = 30)
