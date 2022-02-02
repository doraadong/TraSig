#!/usr/bin/env python

import os
import argparse
import pickle
import requests

import pandas as pd
import numpy as np
import h5py
import rpy2.robjects as robjects
from scipy import stats

"""
Optional parts skipped. See the tutorial for more details.

"""
if __name__ == '__main__':
    # parse command-line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', required=True, default='../input/',
                        help="string, folder to find inputs for trajectory inference")
    parser.add_argument('-o', '--output', required=True, default='../output/',
                        help="string, folder to save inputs for TraSig")
    parser.add_argument('-d', '--project', required=True, help="string, project name")
    parser.add_argument('-t', '--trajectoryFile', required=True,
                        default='../trajectory/output/output.h5', help="string, "
                                    "trajectory output file from dynverse, default ../trajectory/output/output.h5")
    parser.add_argument('-g', '--preprocess', required=True, help="string, preprocessing steps applied to the "
                                                                  "data / project, default 'None'", default="None")
    parser.add_argument('-b', '--modelName', required=True, help="string, name of the trajectory model")
    parser.add_argument('-e', '--otherIdentifier', required=False,
                        default="None", help="string, optional, other identifier for the output, default 'None'")
    parser.add_argument('-c', '--listType', required=False,
                        default='ligand_receptor', help="string, optional, interaction list type, "
                                                        "default ligand_receptor")
    parser.add_argument('-cp', '--pathLR', required=False,
                        default="../ligand_receptor_lists/ligand_receptor_FANTOM.pickle",
                        help="string, optional, path to the interaction list, "
                             "default ../ligand_receptor_lists/ligand_receptor_FANTOM.pickle")
    parser.add_argument('-y', '--genePairType', required=False,
                        default='interaction', help="string, optional, identifier for the type of genes to align, "
                                                    "e.g. interaction/cell_cycle, default interaction")
    parser.add_argument('-yp', '--pathAlign', required=False,
                        default="None", help="string, optional, path to the alignment genes list, "
                                             "set as 'None' if not doing alignment or using 'interaction' for alignment, "
                                             "default 'None'")

    args = parser.parse_args()
    print(args)

    # set parameters for data
    input_path = args.input
    output_path = args.output
    trajectory_filename = args.trajectoryFile
    project = args.project
    preprocess = args.preprocess
    model_name = args.modelName
    others = args.otherIdentifier
    list_type = args.listType
    lr_list_path = args.pathLR
    gene_pair_type = args.genePairType
    align_gene_path = args.pathAlign

    # specificy output names
    if preprocess != "None":
        _preprocess = f"_{preprocess}"
    else:
        _preprocess = ""

    if others == "None":
        condition = ""

    # naming output files
    suffix = f"{_preprocess}_{model_name}{condition}"

    ### Load expression and true labels
    filepath = f"{input_path}/{project}.rds"

    if os.path.exists(filepath):
        pass
    else:
        url = f"https://zenodo.org/record/1443566/files/real/silver/{project}.rds?download=1"
        r = requests.get(url)

        with open(filepath, 'wb') as f:
            f.write(r.content)

    ## Load expression
    filepath = f"{input_path}/{project}.rds"

    from rpy2.robjects import pandas2ri
    pandas2ri.activate()

    readRDS = robjects.r['readRDS']
    df = readRDS(filepath)
    # df = pandas2ri.rpy2py_dataframe(df)
    data_keys = list(df.names)

    cell_ids = df[data_keys.index('cell_ids')]
    expression = df[data_keys.index('expression')]
    genes = df[data_keys.index('feature_info')]['feature_id'].values

    N = len(cell_ids)  # number of cells
    G = len(genes)  # number of genes


    ##  Load true trajectory and labels

    # true trajectory
    milestones_true = df[data_keys.index('milestone_ids')]
    network_true = df[data_keys.index('milestone_network')]
    M_true = len(milestones_true)

    # add node index; node index consistent with index in 'milestone_ids'
    # will use node index to present node from now on
    network_true['idx_from'] = [list(milestones_true).index(i) for i in network_true['from']]
    network_true['idx_to'] = [list(milestones_true).index(i) for i in network_true['to']]

    membership_true = df[data_keys.index('milestone_percentages')]
    # assign cells to the most probable node
    assignment_true = membership_true[membership_true.groupby(['cell_id'])['percentage'].transform(max) == membership_true['percentage']]
    assignment_true.set_index('cell_id', inplace=True)
    assignment_true = assignment_true.reindex(cell_ids)
    clusters_true = [list(milestones_true).index(c) for c in assignment_true['milestone_id'].values]


    ### Load trajectory inference result
    f = h5py.File(trajectory_filename, 'r')

    key = 'data'

    # Get the HDF5 group
    group = f[key]

    _percentages = group['milestone_percentages']
    _network = group['milestone_network']
    _progressions = group['progressions']

    _cell_ids = list(_percentages['data']['cell_id'])
    _cell_ids = [i.decode('utf-8') for i in _cell_ids]
    estimated_percentages = pd.DataFrame(zip(_cell_ids, list(_percentages['data']['milestone_id']), list(_percentages['data']['percentage'])))
    estimated_percentages.columns = ['cell_id', 'milestone_id', 'percentage']

    _cell_ids = list(_progressions['data']['cell_id'])
    _cell_ids = [i.decode('utf-8') for i in _cell_ids]
    estimated_progressions = pd.DataFrame(zip(_cell_ids, list(_progressions['data']['from']), list(_progressions['data']['to']), list(_progressions['data']['percentage'])))
    estimated_progressions.columns = ['cell_id', 'from', 'to' ,'percentage']
    estimated_progressions = estimated_progressions.set_index("cell_id")
    estimated_progressions = estimated_progressions.reindex(assignment_true.index.values)  # assignment_true already reindexed by cell_ids

    estimated_network = pd.DataFrame(pd.DataFrame(zip(list(_network['data']['from']), list(_network['data']['to']), list(_network['data']['length']))))
    estimated_clusters = estimated_percentages.loc[estimated_percentages.groupby(["cell_id"])["percentage"].idxmax()].set_index('cell_id').reindex(cell_ids)
    estimated_clusters['milestone_id'] = [_c.decode("utf-8") for _c in estimated_clusters['milestone_id']]

    ### Prepare and save input for TraSig

    ## Save estimated cluster and progression time
    estimated_progressions['from'] = [i.decode('utf-8') for i in estimated_progressions['from']]
    estimated_progressions['to'] = [i.decode('utf-8') for i in estimated_progressions['to']]
    estimated_progressions['edge'] = estimated_progressions['from'] + '_' + estimated_progressions['to']

    # assign unique label (integer) to each edge
    edges = np.unique(estimated_progressions['edge'])

    edge2idx = {}
    for i, v in enumerate(edges):
        edge2idx[v] = i

    print(f"Edges and their new labels: {edge2idx}")

    estimated_progressions['idx_edge'] = estimated_progressions['edge'].replace(edge2idx)
    hid_var = {'cell_path': estimated_progressions['idx_edge'].values,
              'cell_time': estimated_progressions['percentage'].values,
              'cell_labels':assignment_true['milestone_id'].values}

    # save
    filename = f"{project}{_preprocess}_{model_name}_it2_hid_var.pickle"
    with open(os.path.join(output_path, filename), 'wb') as handle:
        pickle.dump(hid_var, handle, protocol=pickle.HIGHEST_PROTOCOL)

    ## Subsetting expression data
    ## filter expression data to keep only ligands and receptors (and genes used for alignment, if applicable)
    ## interaction list should be a pickle file wiht a python list of tuple: (ligand, receptor)
    ## alignment genes list should be a pickle file wiht a python list of string: gene
    # get interaction file
    with open(os.path.join(lr_list_path), 'rb') as handle:
        interaction_list = pickle.load(handle)

    ligands_receptors = np.unique([i[0] for i in interaction_list] + [i[1] for i in interaction_list])

    # get align genes
    if gene_pair_type != "None" and align_gene_path != "None":
        # get align genes (list of genes)
        with open(os.path.join(align_gene_path), 'rb') as handle:
            align_genes = pickle.load(handle)

        selected_genes = list(np.unique(list(ligands_receptors) + align_genes))
    else:
        align_genes = None
        selected_genes = ligands_receptors

    # get list of genes identified as ligand or receptor
    genes_upper = [g.upper() for g in genes]
    kepted_genes = list(set(genes_upper).intersection(set(selected_genes)))

    df = pd.DataFrame(expression)
    df.columns = genes_upper
    df.index = cell_ids

    df_sub = df[kepted_genes]

    # save filtered expression
    if align_genes == None:
        filename = f"{project}{_preprocess}_{list_type}.txt"
    else:
        filename = f"{project}{_preprocess}_{list_type}_{gene_pair_type}.txt"

    data_file = os.path.join(output_path, filename)
    df_sub.to_csv(data_file)

    # save filtered interactions (list of (ligand, receptor) that are expressed)
    filtered_interactions = []
    for i, j in interaction_list:
        if i in kepted_genes and j in kepted_genes:
            filtered_interactions.append((i, j))

    filename = f"{list_type}_{project}{_preprocess}.pickle"
    with open(os.path.join(output_path, filename), 'wb') as handle:
        pickle.dump(filtered_interactions, handle, protocol=pickle.HIGHEST_PROTOCOL)

    # save filtered align genes
    if align_genes != None:
        align_genes_filtered = list(set(kepted_genes).intersection(set(align_genes)))
        if len(align_genes_filtered) == 0:
            raise ValueError( f"No align genes expressed in the dataset! "
                              f"Consider using another set of genes for alignment!")
        else:
            filename = f"align_{gene_pair_type}_filtered_{project}{_preprocess}.pickle"
            with open(os.path.join(output_path, filename), 'wb') as handle:
                pickle.dump(align_genes_filtered, handle, protocol=pickle.HIGHEST_PROTOCOL)


    ## Save correspondence from sampling time to paths
    cell_ori_time = np.repeat(0, N)  # put all cells at time 0 if sampling time unknow

    unique_days = np.unique(cell_ori_time)
    sorted_days = list(np.sort(unique_days))
    cell_paths = np.unique(hid_var["cell_path"])

    sampleT2path = dict.fromkeys(range(len(sorted_days)))  # use index of sorted sampling time as key
    for k, v in sampleT2path.items():
        sampleT2path[k] = []

    for i, cur_path in enumerate(cell_paths):
        # print("current path (edge)", cur_path)

        # get data corresponding to a path
        condition = hid_var["cell_path"] == cur_path
        cur_days = np.array(cell_ori_time)[condition]

        # get the sampling time for the majority cells
        mode, count = stats.mode(cur_days)
        # print(
        #     f"Sampling time for the majority of cells: {mode[0]}, making {round(float(count[0]) / len(cur_days), 2)}% percent")
        cur_sampleT = mode[0]

        # will use index instead of input time
        sampleT2path[sorted_days.index(cur_sampleT)].append(cur_path)

    # save the dictionary
    filename = 'sampling_time_per_path_' + project + suffix + '.pickle'

    with open(os.path.join(output_path, filename), 'wb') as handle:
        pickle.dump(sampleT2path, handle, protocol=pickle.HIGHEST_PROTOCOL)