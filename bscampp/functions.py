import json, time, os, sys
import treeswift
from collections import defaultdict, Counter

from bscampp import get_logger, log_exception
from bscampp.configs import Configs
from bscampp.jobs import EPAngJob, TaxtasticJob, PplacerTaxtasticJob
from bscampp.utils import write_fasta
import bscampp.utils as utils

_LOG = get_logger(__name__)

'''
Function to read in the placement tree and alignment.
If query alignment is provided, will use the provided query instead of
the ones (potentially) included in the reference alignment
'''
def readData(workdir, dry_run=False):
    t0 = time.perf_counter()
    _LOG.info('Reading in input data...')

    if dry_run:
        return None, dict(), '', dict(), '', dict()

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
        write_fasta(qaln_path, q_dict)
        
        aln_path = os.path.join(workdir, 'aln.fa')
        write_fasta(aln_path, ref_dict)

    t1 = time.perf_counter()
    _LOG.info('Time to read in input data: {} seconds'.format(t1 - t0))
    return tree, leaf_dict, aln_path, ref_dict, qaln_path, q_dict

'''
Function to get the closest leaf for each query sequence based on Hamming
distance
'''
def getClosestLeaves(aln_path, qaln_path, aln, qaln, workdir, dry_run=False):
    t0 = time.perf_counter()
    _LOG.info('Computing closest leaves for query sequences...')

    if dry_run:
        return dict(), dict()

    query_votes_dict = dict()
    query_top_vote_dict = dict()
    tmp_output = os.path.join(workdir, 'closest.txt') 
    
    if Configs.subtreetype == "h":
        Configs.votes = Configs.subtreesize

    cmd = []
    if Configs.similarityflag:
        cmd.append(os.path.join(Configs.hamming_distance_dir, 'homology'))
    else:
        if Configs.fragmentflag == False:
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

'''
Function to assign queries to subtrees based on their votes
'''
def assignQueriesToSubtrees(query_votes_dict, query_top_vote_dict,
        tree, leaf_dict, dry_run=False):
    t0 = time.perf_counter()
    _LOG.info('Adding query votes to the placement tree...')

    if dry_run:
        return dict(), []

    # (1) go over the query votes and add them to corresponding leaves
    lf_votes = Counter()
    leaf_queries = dict()
    for name, y in query_votes_dict.items():
        lf_votes.update(y)
        for ind, leaf in enumerate(y):
            top_vote = False
            if ind == 0:
                top_vote = True
            if leaf not in leaf_queries:           
                leaf_queries[leaf] = {(name,top_vote)}
            else:
                leaf_queries[leaf].add((name,top_vote))

    subtree_dict = dict()
    subtree_leaf_label_dict = dict()
    most_common_index = 0
    
    # assign queries to subtrees, and remove them from the pool
    # repeat until all queries are assigned
    while len(query_votes_dict) > 0:
        _LOG.info("queries left to assign: {}".format(len(query_votes_dict)))
        (seed_label, node_votes) = lf_votes.most_common(
                most_common_index+1)[most_common_index]
        
        node_y = leaf_dict[seed_label]
        # extract [subtreesize] leaves
        labels = utils.subtree_nodes_with_edge_length(tree, node_y,
                Configs.subtreesize)
        subtree = tree.extract_tree_with(labels)
        label_set = set(labels)

        queries_by_subtree = set()
        subtree_query_set = set()

        # gather any other queries that can be used with this subtree
        for label in labels:
            leaf_queries_remove_set = set()
            if label in leaf_queries:
                    
                for leaf_query, top_vote in leaf_queries[label]:
                
                    if leaf_query not in query_votes_dict:
                        leaf_queries_remove_set.add((leaf_query, top_vote))
                        continue
                        
                    if top_vote:
                        subtree_query_set.add(leaf_query)
                        leaf_queries_remove_set.add((leaf_query, top_vote))
                    
                leaf_queries[label].difference_update(leaf_queries_remove_set)
        queries_by_subtree.update(subtree_query_set)

        if len(queries_by_subtree) > 0:
            subtree_dict[subtree] = (seed_label, queries_by_subtree)
            subtree_leaf_label_dict[subtree] = subtree.label_to_node(
                    selection='leaves')

        votes_b4 = len(list(lf_votes.elements()))
        for query in queries_by_subtree:
            if query in query_votes_dict:
                lf_votes.subtract(query_votes_dict[query])
                query_votes_dict.pop(query)

        if len(queries_by_subtree) == 0:
            # 10.27.2023 - Chengze Shen
            # >>> prevent going over the the total number of votes
            most_common_index += 1
        else:
            most_common_index = 0
            
    placed_query_list = []
    
    # reassign queries to the subtree minimizing total edge length 
    # from the query's top vote to the subtree's seedleaf
    new_subtree_dict = dict()
    for query, closest_label in query_top_vote_dict.items():
        best_subtree = None
        best_distance = 99999999999999999
        for subtree, value in subtree_dict.items():
            leaf_label_dict = subtree_leaf_label_dict[subtree]
            seed_label, _ = value
            if closest_label in leaf_label_dict:
                distance = subtree.distance_between(
                        leaf_label_dict[closest_label],
                        leaf_label_dict[seed_label])
                if distance < best_distance:
                   best_subtree = subtree
                   best_distance = distance
        if best_subtree in new_subtree_dict:
            new_subtree_dict[best_subtree].append(query)
        else:
            new_subtree_dict[best_subtree] = [query]

    t1 = time.perf_counter()
    _LOG.info('Time to assign queries to subtrees: {} seconds'.format(t1 - t0))
    return new_subtree_dict, placed_query_list


'''
Function to assign queries to subtrees as used in SCAMPP 
(subtrees are built using the nearest leaf as the seed sequence)
'''
def buildQuerySubtrees(query_votes_dict, query_top_vote_dict,
        tree, leaf_dict, dry_run=False):
    t0 = time.perf_counter()
    _LOG.info('(SCAMPP) Building query subtree for placement...')

    if dry_run:
        return dict(), []

    # (1) go over the query seed sequences to see if any queries use 
    # the same seed sequence (i.e. subtree)
    seed_queries = dict()
    for query, closest_leaf in query_top_vote_dict.items():
        if closest_leaf not in seed_queries:           
            seed_queries[closest_leaf] = [query]
        else:
            seed_queries[closest_leaf].append(query)

    new_subtree_dict = dict()
    # assign queries to subtrees, and remove them from the pool
    # repeat until all queries are assigned
    _total = 0
    for seed_label, queries in seed_queries.items():
        ####### additional logging for tracking progress 
        _total += 1
        if _total % 1000 == 0 or _total == len(seed_queries):
            _LOG.info(f"- Built {_total}/{len(seed_queries)} subtrees")

        node_y = leaf_dict[seed_label]
        # extract [subtreesize] leaves
        if Configs.subtreetype == "h":
            labels = query_votes_dict[queries[0]]
        elif Configs.subtreetype == "n":
            labels = utils.subtree_nodes(tree, node_y, Configs.subtreesize)
        else:
            labels = utils.subtree_nodes_with_edge_length(tree, node_y,
                Configs.subtreesize)
        subtree = tree.extract_tree_with(labels)
        new_subtree_dict[subtree] = queries

    placed_query_list = []

    t1 = time.perf_counter()
    _LOG.info('Time to assign queries to subtrees: {} seconds'.format(t1 - t0))
    return new_subtree_dict, placed_query_list

'''
Helper function to run a single placement task. Designed to use with
multiprocessing
'''
def placeOneSubtree():
    # TODO
    pass

'''
Function to perform placement of queries for each subtree
'''
def placeQueriesToSubtrees(tree, leaf_dict, new_subtree_dict, placed_query_list,
        aln, qaln, cmdline_args, workdir, pool, lock, dry_run=False):
    t0 = time.perf_counter()
    _LOG.info('Performing placement on each subtree...')

    if dry_run:
        return dict()

    # prepare to write an aggregated results to local
    jplace = dict()
    utils.add_edge_nbrs(tree)
    jplace["tree"] = utils.newick_edge_tokens(tree)
    placements = []

    # go over the dictionary of subtrees and their assigned queries
    # perform placement using either EPA-ng or pplacer
    final_subtree_count, total_subtrees_examined = 0, 0
    for subtree, query_list in new_subtree_dict.items():
        total_subtrees_examined += 1
        _LOG.info('- Subtree {}/{} with {} queries'.format(
            total_subtrees_examined, len(new_subtree_dict), len(query_list)))

        # empty subtree, continue
        if len(query_list) == 0:
            continue
        final_subtree_count += 1

        subtree_dir = os.path.join(workdir, f'subtree_{final_subtree_count}')
        if not os.path.isdir(subtree_dir):
            os.makedirs(subtree_dir)
        
        # name all temporary output files
        tmp_tree = os.path.join(subtree_dir, 'tree')
        tmp_aln = os.path.join(subtree_dir, f'subtree_{final_subtree_count}_aln.fa')
        tmp_qaln = os.path.join(subtree_dir, f'subtree_{final_subtree_count}_qaln.fa')
        tmp_output = os.path.join(subtree_dir,
                'subtree_{}_{}.jplace'.format(
                    final_subtree_count, Configs.placement_method))

        # extract corresponding ref sequences and queries
        tmp_leaf_dict = subtree.label_to_node(selection='leaves')
        if '' in tmp_leaf_dict:
            del tmp_leaf_dict['']
        tmp_ref_dict = {label : aln[label] for label in tmp_leaf_dict.keys()}
        tmp_q_dict = {name : qaln[name] for name in query_list}
        write_fasta(tmp_aln, tmp_ref_dict)
        write_fasta(tmp_qaln, tmp_q_dict)

        # process the subtree before placement
        subtree.resolve_polytomies()
        subtree.suppress_unifurcations()
        subtree.write_tree_newick(tmp_tree, hide_rooted_prefix=True)

        # 1.27.2025 - Chengze Shen
        # choose the placement method to run
        if Configs.placement_method == 'epa-ng':
            job = EPAngJob(path=Configs.epang_path,
                    info_path=Configs.info_path, tree_path=tmp_tree,
                    aln_path=tmp_aln, qaln_path=tmp_qaln,
                    outdir=subtree_dir, num_cpus=Configs.num_cpus)
            # for EPA-ng, ensure that outpath name is changed to the one we want
            _outpath = job.run(logging=f'subtree_{final_subtree_count}')
            os.system('mv {} {}'.format(_outpath, tmp_output))
        elif Configs.placement_method == 'pplacer':
            # build ref_pkg with info and tmp_tree and tmp_aln
            refpkg_dir = os.path.join(subtree_dir,
                    f'subtree_{final_subtree_count}.refpkg')
            taxit_job = TaxtasticJob(path=Configs.taxit_path,
                    outdir=refpkg_dir, name=f'subtree_{final_subtree_count}',
                    aln_path=tmp_aln, tree_path=tmp_tree,
                    info_path=Configs.info_path)
            _ = taxit_job.run()

            # run pplacer-taxtastic
            job = PplacerTaxtasticJob(path=Configs.pplacer_path,
                    refpkg_dir=refpkg_dir, model=Configs.model,
                    outpath=tmp_output, num_cpus=Configs.num_cpus,
                    qaln_path=tmp_qaln)
            tmp_output = job.run(logging=f'subtree_{final_subtree_count}')
        else:
            raise ValueError(
                    f"Placement method {Configs.placement_method} not recognized")

        # read in each placement result
        place_file = open(tmp_output, 'r')
        place_json = json.load(place_file)
        tgt = "n"
        if Configs.placement_method == 'pplacer':
            tgt = "nm"
        if len(place_json["placements"]) > 0:
            added_tree, edge_dict = utils.read_tree_newick_edge_tokens(
                    place_json["tree"])

            # obtain the fields for "p"
            fields = place_json["fields"]
            # set the fields in jplace accordingly
            if "fields" not in jplace:
                jplace["fields"] = fields
            field_to_idx = {field: i for i, field in enumerate(fields)}

            for tmp_place in place_json["placements"]:
                #print(tmp_place)
                placed_query_list.append(tmp_place[tgt][0])
                for i in range(len(tmp_place["p"])):
                    edge_num = tmp_place["p"][i][
                            field_to_idx['edge_num']]
                    edge_distal = tmp_place["p"][i][
                            field_to_idx['distal_length']]

                    right_n = edge_dict[str(edge_num)]
                    left_n = right_n.get_parent()

                    # left and right path_l and path_r are in added_tree
                    left, path_l = utils.find_closest(left_n, {left_n, right_n})
                    right, path_r = utils.find_closest(right_n, {left_n, right_n})

                    left = leaf_dict[left.get_label()]
                    right = leaf_dict[right.get_label()]
                    _, path = utils.find_closest(left, {left}, y=right)
                    # now left right and path are in tree

                    length = sum([x.get_edge_length() for x in path_l])+edge_distal
                    # left path length through subtree before placement node

                    target_edge = path[-1]

                    for j in range(len(path)):
                        length -= path[j].get_edge_length()
                        if length < 0:
                            target_edge = path[j]
                            break

                    #tmp_place["p"][i][field_to_idx['edge_num']] = 0

                    label = target_edge.get_label()

                    [taxon, target_edge_nbr] = label.split('%%',1)
                    tmp_place["p"][i][field_to_idx['distal_length']] = \
                            target_edge.get_edge_length()+length
                    tmp_place["p"][i][field_to_idx['edge_num']] = \
                            int(target_edge_nbr)

                placements.append(tmp_place.copy())
        place_file.close()
    _LOG.info(f'Final number of subtrees used: {final_subtree_count}')

    # prepare the output jplace to write
    jplace["placements"] = placements
    jplace["metadata"] = {"invocation": " ".join(cmdline_args)}
    jplace["version"] = 3
    #jplace["fields"] = ["distal_length", "edge_num", "like_weight_ratio", \
    #        "likelihood", "pendant_length"]

    t1 = time.perf_counter()
    _LOG.info('Time to place queries to subtrees: {} seconds'.format(t1 - t0))
    return jplace

'''
Function to write a given jplace object to local output
'''
def writeOutputJplace(output_jplace, dry_run=False):
    t0 = time.perf_counter()
    _LOG.info('Writing aggregated placements to local...')
    
    if dry_run:
        return

    outpath = os.path.join(Configs.outdir, Configs.outname)
    outf = open(outpath, 'w')
    json.dump(output_jplace, outf, sort_keys=True, indent=4)
    outf.close()

    t1 = time.perf_counter()
    _LOG.info('Time to build final jplace file: {} seconds'.format(t1 - t0))
