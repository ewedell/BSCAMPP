import time, os, sys
import treeswift

from bscampp import get_logger, log_exception
from bscampp.configs import Configs
import bscampp.utils as utils

_LOG = get_logger(__name__)

'''
Function to read in the placement tree and alignment.
If query alignment is provided, will use the provided query instead of
the ones (potentially) included in the reference alignment
'''
def readData(workdir):
    t0 = time.perf_counter()
    _LOG.info('Reading in input data...')

    # (1) load reference tree
    tree = treeswift.read_tree_newick(Configs.tree_path)
    tree.resolve_polytomies()

    leaf_dict = tree.label_to_node(selection='leaves')
    # clean the leaf keys so that ' or " are not present
    ori_keys = list(leaf_dict.keys())
    for key in ori_keys:
        _node = leaf_dict[key]
        new_key = key.replace('\'', '')
        new_key = new_key.replace('\"', '')
        leaf_dict.pop(key)
        leaf_dict[new_key] = _node

    # (2) load reference alignment and query alignment (if provided) 
    if Configs.qaln_path is not None: 
        ref_dict = utils.read_data(Configs.aln_path)
        q_dict = utils.read_data(Configs.qaln_path)
        aln_path, qaln_path = Configs.aln_path, Configs.qaln_path
    else:
        aln_dict = utils.read_data(Configs.aln_path)
        ref_dict, q_dict = utils.seperate(aln_dict, leaf_dict)

        # after separating queries from the reference alignment, write
        # them to to TEMP/
        qaln_path = os.path.join(workdir, 'qaln.fa')
        write_fasta(temp_qaln_path, q_dict)
        
        aln_path = os.path.join(workdir, 'aln.fa')
        write_fasta(temp_aln_path, ref_dict)

    t1 = time.perf_counter()
    _LOG.info('Time to read in input data: {} seconds'.format(t1 - t0))
    return tree, aln_path, ref_dict, qaln_path, q_dict

'''
Function to get the closest leaf for each query sequence based on Hamming
distance
'''
def getClosestLeaves(aln_path, qaln_path, aln, qaln, workdir):
    t0 = time.perf_counter()
    _LOG.info('Computing closest leaves for query sequences...')
    query_votes_dict = dict()
    query_top_vote_dict = dict()
    
    tmp_output = os.path.join(workdir, 'closest.txt') 

    cmd = []
    if Configs.similarityflag:
        cmd.append(os.path.join(Configs.hamming_distance_dir, 'homology'))
    else:
        if fragment_flag == False:
            cmd.append(os.path.join(Configs.hamming_distance_dir, 'hamming'))
        else: 
            cmd.append(os.path.join(
                Configs.hamming_distance_dir, 'fragment_hamming'))
    cmd.extend([aln_path, str(len(aln)), qaln_path, str(len(qaln)),
        tmp_output, str(Configs.votes)])
    os.system(' '.join(cmd))

    # process closest leaves 
    unusable_queries = set()
    f = open(tmp_output)
    for line in f:
        line = line.strip()
        y = line.split(',')
        name = y.pop(0)
        for idx, taxon in enumerate(y):
            leaf, hamming = taxon.split(':')
            y[idx] = (leaf, int(hamming))

        y = sorted(y, key=lambda x: x[1])
        for idx, taxon in enumerate(y):
            y[idx] = taxon[0]

        if name.find(':') >= 0:
            name_list = name.split(":")
            name = name_list[0]
            ungapped_length = name_list[1]
            if y[0] == ungapped_length:
                _LOG.warning(f'Sequence {name}: no homologous sites found, '
                        'removed before placement.')
                unusable_queries.add(name)
        if name not in unusable_queries:
            query_votes_dict[name] = y
            query_top_vote_dict[name] = y[0]
    f.close()
    
    t1 = time.perf_counter()
    _LOG.info('Time to compute closest leaves: {} seconds'.format(t1 - t0)) 
    return query_votes_dict, query_top_vote_dict
