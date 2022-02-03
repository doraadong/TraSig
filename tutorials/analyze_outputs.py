#!/usr/bin/env python

import pickle
import os
import argparse

import numpy as np
import pandas as pd

# load packages required for analysis
import statsmodels.api as sm
import statsmodels as sm
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

from trasig.utils import str2bool

if __name__ == '__main__':
    # parse command-line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', required=True, default='../input/',
                        help="string, folder to find TraSig's inputs")
    parser.add_argument('-o', '--output', required=True, default='../output/',
                        help="string, folder to find TraSig's outputs")
    parser.add_argument('-d', '--project', required=True, help="string, project name")
    parser.add_argument('-g', '--preprocess', required=True, help="string, preprocessing steps applied to the "
                                                                  "data / project, default None", default="None")
    parser.add_argument('-b', '--modelName', required=True, help="string, name of the trajectory model")
    parser.add_argument('-t', '--listType', required=False,
                        default='ligand_receptor', help="string, optional, "
                                                        "interaction list type, default ligand_receptor")
    parser.add_argument('-p', '--otherIdentifier', required=False,
                        default="None", help="string, optional, other identifier for the output, default None")
    parser.add_argument('-l', '--nLap', required=False, default=20, help="integer, optional, "
                                                                         "sliding window size, default 20")
    parser.add_argument('-m', '--metric', required=False, default='dot', help="string, optional, "
                                                                              "scoring metric, default dot")
    parser.add_argument('-z', '--nan2zero', required=False, type=str2bool,
                        default=True, help="boolean, optional, if treat nan as zero, default True")
    parser.add_argument('-n', '--numPerms', required=False,
                        default=10000, help="integer, optional, number of permutations, default 10000")
    parser.add_argument('-s', '--startingTreatment', required=False,
                        default="smallerWindow", help="string, optional, way to treat values at the beginning of an "
                                                      "edge with sliding window size smaller than nLap, "
                                                      "None/parent/discard/smallerWindow, default smallerWindow, "
                                                      "need to provide an extra input 'path_info.pickle' "
                                                      "for 'parent' option")
    # the following specifiy alignment related arguements
    parser.add_argument('-a', '--alignType', required=False,
                        default='unaligned', help="string, optional, how to align edges, "
                                                  "options: unaligned/aligned-fixed/aligned-specific, "
                                                  "default unaligned")
    parser.add_argument('-y', '--genePairType', required=False,
                        default='interaction', help="string, optional, identifier for the type of genes to align, "
                                                    "e.g. interaction/cell_cycle, default interaction")
    parser.add_argument('-f', '--smooth', required=False,
                        default=1, help="float, optional, smoothing parameter for splines, default 1")
    parser.add_argument('-v', '--overlap', required=False,
                        default=0.5, help="float, optional, overlap threshold for alignment, default 0.5")
    parser.add_argument('-r', '--rate', required=False,
                        default=1, help="integer, optional, sampling rate for aligned time points, default 1")
    parser.add_argument('-e', '--errorType', required=False, default="cosine",
                        help="string, optional, type of distance metric for alignment (MSE, cosine or corr), "
                             "default cosine")
    parser.add_argument('-k', '--aRate', required=False,
                        default=0.05, help="float, optional, rate to sample parameter a for alignment, default 0.05")
    parser.add_argument('-j', '--bRate', required=False,
                        default=2.5, help="float, optional, rate to sample parameter b for alignment, default 2.5")

    args = parser.parse_args()
    print(args)

    # set parameters for data
    input_path = args.input
    output_path = args.output
    project = args.project
    preprocess = args.preprocess
    model_name = args.modelName
    list_type = args.listType
    others = args.otherIdentifier

    if preprocess != "None":
        _preprocess = f"_{preprocess}"
    else:
        _preprocess = ""

    if others == "None":
        others = ""

    # set parameters for calculating metrics
    n_lap = int(args.nLap)
    metrics = [args.metric]
    nan2zero = args.nan2zero
    num_perms = int(args.numPerms)
    startingTreatment = args.startingTreatment

    if startingTreatment != "None":
        _startingTreatment = f"_{startingTreatment}"
    else:
        _startingTreatment = ""

    # set parameters for alignment
    align_type = args.alignType

    if align_type != "unaligned":
        SMOOTH_PARAMETER = float(args.smooth)
        OVERLAP_THRESHOLD = float(args.overlap)
        SAMPLING_RATE = int(args.rate)
        gene_pair_type = args.genePairType
        error_type = args.errorType
        a_rate = float(args.aRate)
        b_rate = float(args.bRate)

    else:
        pass

    ### load inputs
    suffix = f"{project}_{list_type}{_preprocess}_{model_name}"
    suffix = f"{suffix}{_startingTreatment}_nlap_{n_lap}{others}"

    # set parameters for alignment
    if align_type != "unaligned":
        suffix = f"{suffix}_{gene_pair_type}_{SMOOTH_PARAMETER}_{error_type}_{a_rate}_{b_rate}"
    else:
        pass

    # get interaction file (list of (ligand, receptor/target))
    filename = f"{list_type}_{project}{_preprocess}.pickle"
    with open(os.path.join(input_path, filename), 'rb') as handle:
        interaction_list = pickle.load(handle)

    # load expression data
    filename = f"{project}{_preprocess}_{list_type}.txt"
    if align_type != "unaligned" and gene_pair_type != "interaction":
        filename = f"{project}{_preprocess}_{list_type}_{gene_pair_type}.txt"
    print("Load: ", filename)

    data_file = os.path.join(input_path, filename)
    df = pd.read_csv(data_file, index_col=0)
    cell_exps = df.values
    gene_names = list(df.columns.values)  # assume unique

    # (optional) load corresponding between sampling time and path
    filename = f"sampling_time_per_path_{project}{_preprocess}_{model_name}.pickle"
    with open(os.path.join(input_path, filename), 'rb') as handle:
        time2path = pickle.load(handle)

    path2time = dict()
    for k, ps in time2path.items():
        for p in ps:
            path2time[p] = k

    # load path & time assignment
    # original assignment
    hid_var_file = f"{project}{_preprocess}_{model_name}_it2_hid_var.pickle"
    with open(os.path.join(input_path, hid_var_file), 'rb') as handle:
        hid_var = pickle.load(handle, encoding="latin1")

    unique_paths = np.unique(hid_var["cell_path"])
    all_times = [round(i, 2) for i in np.arange(0, 1.01, 0.01)]  # all possible labels for cell time
    cell_paths_o = hid_var["cell_path"]
    cell_times_o = hid_var["cell_time"]

    # put time (from 0 to 1) into 101 bins
    cell_times_o = np.round(cell_times_o, 2)

    ### load outputs
    # load the scores on the original data
    _n = 0

    _columns = dict.fromkeys(metrics)
    for m in metrics:
        _columns[m] = []

    _columns.update({'pair': [], 'gene_pair_id': []})

    # load results
    filename = f"{suffix}_metrics_{_n}.pickle"
    data_file = os.path.join(output_path, filename)

    with open(data_file, 'rb') as handle:
        results = pickle.load(handle)

    for pair, mets in results.items():
        for m in metrics:
            _columns[m] += list(mets[m])

        _columns['pair'] += list(np.repeat(pair, len(mets[m])))
        _columns['gene_pair_id'] += list(range(len(mets[m])))

    df = pd.DataFrame(_columns)
    num_pairs = len(results[pair][m])

    # load permutation results
    suffix = f"{suffix}_{int(np.log10(num_perms))}"
    if align_type != "unaligned":
        suffix = f"{suffix}_{align_type}"
    else:
        pass

    child_suffix = f"{suffix}_{metrics[0]}"

    filename = f"{suffix}_permutation_results.pickle"
    data_file = os.path.join(output_path, filename)

    with open(data_file, 'rb') as handle:
        pair2counts = pickle.load(handle)

    # turn to p-values
    for pair, _ in pair2counts.items():
        for m in metrics:
            pair2counts[pair][m] = (pair2counts[pair][m] + 1) / (num_perms + 1)

    # add to the dataframe
    _columns = dict.fromkeys(metrics)
    for m in metrics:
        _columns[m] = []

    for pair, counts in pair2counts.items():
        for m in metrics:
            _columns[m] += list(counts[m])

    for m in metrics:
        df[f"{m}_p"] = _columns[m]

    # add ligand target info
    df['ligand'] = [interaction_list[int(i)][0] for i in df['gene_pair_id']]
    df['target'] = [interaction_list[int(i)][1] for i in df['gene_pair_id']]
    ligand_list = np.unique(df['ligand'])

    # add more info about cell clusters
    df['sender'] = [i.split('_')[0] for i in df['pair']]
    df['receiver'] = [i.split('_')[1] for i in df['pair']]
    df['sender'] = df['sender'].astype('int')
    df['receiver'] = df['receiver'].astype('int')
    df['time-sender'] = [path2time[i] for i in df['sender']]
    df['time-receiver'] = [path2time[i] for i in df['receiver']]

    ## label clusters using the true labels of the majority of cells (for plotting)
    # build path2label
    unique_days = np.unique(hid_var['cell_labels'])
    cell_paths = np.unique(hid_var["cell_path"])

    _dict = dict.fromkeys(range(len(cell_paths)))

    for i, cur_path in enumerate(cell_paths):
        # print("------current path", cur_path)

        # get data corresponding to a path
        condition = hid_var["cell_path"] == cur_path
        cur_labels = hid_var['cell_labels'][condition]
        try:
            cur_labels = [i.decode('UTF-8') for i in cur_labels]
        except AttributeError:
            pass

        # get the sampling time for the majority cells
        mode, count = stats.mode(cur_labels)

        major_percent = round(float(count[0]) / len(cur_labels), 2)
        # print(mode[0], major_percent)

        cur_label = mode[0]

        # add more labels if cells of the major cell type make less than 90% of the whole population
        if major_percent < 0.9:
            cur_label += '(' + str(major_percent) + ')'

            labels, counts = np.unique(cur_labels, return_counts=True)
            sorted_counts, idxs = np.unique(counts, return_index=True)
            #  print(zip(sorted_counts, labels[idxs]))

            count = 0
            while major_percent < 0.9:
                # add more labels until major_percent >= 0.9
                add_counts = sorted_counts[::-1][1 + count]
                _add_percent = round(add_counts / len(cur_labels), 2)

                major_percent += _add_percent
                # print(major_percent)

                cur_label += '\n '
                cur_label += labels[idxs][::-1][1 + count]
                cur_label += '(' + str(round(_add_percent, 2)) + ')'

                count += 1

        _dict[cur_path] = cur_label

    path2label = _dict


    ## Adjust p-values for multiple comparisons
    _p = df['dot_p'].values.copy()

    for pair in results.keys():
        condition = np.where(df['pair'] == pair)[0]
        adjusted = sm.stats.multitest.fdrcorrection(df['dot_p'].values[condition])
        _p[condition] = adjusted[1]

    df['dot_p_adjusted'] = _p

    ### Infer interactions among cell clusters (edges)
    # reset output path to save analysis results
    output_path = f"{output_path}/analysis"

    if not os.path.exists(output_path):
        os.makedirs(output_path)

    print(f"Analysis outputs to be saved at {output_path}")

    df_pool = pd.DataFrame(list(set(df['pair'])))
    df_pool.columns = ['pair']

    df_pool['sender'] = [i.split('_')[0] for i in df_pool['pair']]
    df_pool['receiver'] = [i.split('_')[1] for i in df_pool['pair']]
    df_pool['sender'] = df_pool['sender'].astype('int')
    df_pool['receiver'] = df_pool['receiver'].astype('int')
    df_pool['time-sender'] = [path2time[i] for i in df_pool['sender']]
    df_pool['time-receiver'] = [path2time[i] for i in df_pool['receiver']]
    # df_pool = df_pool[df_pool['time-sender'] == df_pool['time-receiver']]  # if only keep pairs sampled at the same time

    ## Calculate summary score over all ligand-receptor pairs
    cutoff = 0.05
    name_p = 'dot_p_adjusted'

    _counts = []
    for p in df_pool['pair']:
        condition = df['pair'] == p
        _counts.append((df[condition][name_p] < cutoff).sum())

    df_pool['counts'] = _counts

    # subset only contains significant pairs
    condition = df[name_p] < cutoff
    df_sig = df[condition].copy()
    df_sig.reset_index(inplace=True)

    # order clusters (edges / paths) by sampling time
    path_order_time = []
    for k, v in time2path.items():
        path_order_time = path_order_time + v

    df_pool['sender'] = pd.Categorical(df_pool['sender'], path_order_time)
    df_pool['receiver'] = pd.Categorical(df_pool['receiver'], path_order_time)
    df_pool.sort_values(['sender', 'receiver'], inplace=True)

    _vmin = min(df_pool['counts'])
    _vmax = max(df_pool['counts'])

    method = "TraSig"

    metric = 'counts'
    _center_value = _vmin
    # plot only pairs sampled at the same time
    df_plot = df_pool[df_pool['time-sender'] == df_pool['time-receiver']].pivot(index='sender', columns='receiver', values=metric)

    # sort by column names
    df_plot = df_plot.sort_index(axis=1)

    # sns.set_style("white")
    plt.figure(figsize=(5, 5))
    sns.set_context("paper", font_scale=2)

    ax = sns.heatmap(df_plot.values, xticklabels=True, yticklabels=True,
                     vmin=_vmin, vmax=_vmax, center=_center_value, cmap="RdBu_r")
    plt.xticks(rotation=90)
    plt.ylabel("Sender")
    plt.xlabel('Receiver')
    ax.set_xticklabels(df_plot.index.values)
    ax.set_yticklabels(df_plot.index.values)

    if 'tf' in model_name:
        _traj = 'CSHMM'
        plt.title(f"{method} using \n output from {_traj}")
    else:
        _traj = model_name.split('_')[1].capitalize()
        plt.title(f"{method} using \n output from {_traj}")

    filename = f"{child_suffix}_summary_scores.png"
    plt.savefig(os.path.join(output_path, filename), bbox_inches="tight", dpi=300, format="png")
    filename = f"{child_suffix}_summary_scores.eps"
    plt.savefig(os.path.join(output_path, filename), bbox_inches="tight", dpi=300, format="eps")

    # save summary score
    df_pool['sender-label'] = df_pool['sender'].replace(path2label)
    df_pool['receiver-label'] = df_pool['receiver'].replace(path2label)

    cols_order = ['sender', 'sender-label', 'receiver', 'receiver-label', 'counts']
    df_out = df_pool[cols_order].copy()

    filename = f"{child_suffix}_summary_score.csv"
    df_out.to_csv(os.path.join(output_path, filename), index=False)

    ## Save significant ligand-receptor pairs
    df_sig['sender-label'] = df_sig['sender'].replace(path2label)
    df_sig['receiver-label'] = df_sig['receiver'].replace(path2label)

    # sort ligand-receptors in each cluster pair by their scores
    _dfs = []
    pairs_ts = np.unique(df_sig['pair'])
    for pair in pairs_ts:
        condition = df_sig['pair'] == pair
        _dfs.append(df_sig[condition].sort_values('dot', ascending=False))

    df_sorted = pd.concat(_dfs)
    cols_order = ['pair', 'time-sender', 'sender', 'sender-label', 'time-receiver', 'receiver', 'receiver-label',
                  'ligand', 'target', 'dot', 'dot_p', 'dot_p_adjusted']
    df_sorted = df_sorted[cols_order]

    df_sorted.columns = ['interaction pair', 'sender sampling time', 'sender', 'sender-label', 'receiver sampling time',
                         'receiver', 'receiver-label', 'ligand', 'target', 'score', 'score p-value',
                         'score p-value adjusted']

    filename = f"{child_suffix}_significant_pairs.csv"
    df_sorted.to_csv(os.path.join(output_path, filename), index=False)