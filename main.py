#!/usr/bin/env python

"""

Updaets log:
10-21-21: for the expression (filtered by ligand-receptor pairs), instead of having "_lr.txt"; use
                 "_{list_type}.txt"
9-17-21: add alignment using Dynamic Time Warping (DTW) option
4-19-21: add startingTreatment option

"""

import argparse
import pickle
import os
import gc
import copy

import numpy as np
import bottleneck as bn
import pandas as pd

from trasig.trasig import trasig
from trasig.alignment import DTW
from trasig.utils import MPR, str2bool, getPairedPaths

# multi-threading function
def process_per_permutation(ns, save=False):
    """

    Given expressions and a particular cluster assignment, calculate the metrics.

    Global: results, cell_paths_o, cell_times_o, cell_exps, unique_paths, time2path,
    path2time, all_times, metrics, n_lap, nan2zero, gene_names

    Parameters
    ----------
    ns: list or generator
        permutation indexes
    save: boolean
        if save scoring result for the current permutation

    Returns
    -------
    _counts: dictionary
        pair of paths / edges (string) as keys, counts as values; count representing number of permutations
        as least as extreme as the observed value

    """
    saveCounts = True

    try:
        _counts = copy.deepcopy(pair2counts)
    except NameError:  # when the true results not exist
        saveCounts = False
        _counts = None

    for _n in ns:

        if _n != 0:  # run for the permutation
            cell_paths, cell_times = trasig.generate_permutation(_n)
        else:  # run for the original data
            cell_paths, cell_times = cell_paths_o, cell_times_o

        # path to idx
        path2idx = dict.fromkeys(unique_paths)

        # get sum expression & number of cells per bin (time)
        exp = []
        count = []

        for i, cur_path in enumerate(unique_paths):
            # get data corresponding to a path
            condition = cell_paths == cur_path
            cur_exp = cell_exps[condition]
            cur_time = cell_times[condition]

            # inspect how many cells assigned to each location(time on path)
            df_exp = pd.DataFrame(cur_exp[:, 0])
            df_exp['time'] = cur_time  # all values has 2 digits in the dataframe
            count_bin = df_exp.groupby('time').count()
            count_bin.reset_index(inplace=True)

            # get sum expression per gene per bin/time
            df_exp = pd.DataFrame(cur_exp)
            df_exp['time'] = cur_time  # all values has 2 digits in the dataframe
            sum_bin = df_exp.groupby('time').sum()
            sum_bin.reset_index(inplace=True)

            # get all available bins
            df_at = pd.DataFrame(all_times)  # all values has 2 digits in the dataframe
            df_at.columns = ['time']

            # get full set of expression for all bins; assign np.nan if no cells
            sum_full = sum_bin.merge(df_at, on='time', how='right')
            sum_full.sort_values('time', inplace=True)
            sum_full.reset_index(inplace=True, drop=True)
            sum_full.drop(['time'], axis=1, inplace=True)

            # get full set of counts per bins; assign np.nan if no cells
            count_full = count_bin.merge(df_at, on='time', how='right')
            count_full.sort_values('time', inplace=True)
            count_full.reset_index(inplace=True, drop=True)
            count_full.drop(['time'], axis=1, inplace=True)

            exp.append(sum_full.values[:, None])
            count.append(count_full.iloc[:, 0].values[:, None])
            path2idx[cur_path] = i

        exp = np.concatenate(exp, axis=1)  # bin size x number of paths x genes
        count = np.concatenate(count, axis=1)  # bin size x number of paths

        # deal with initial several values (with window size smaller than n_lap)
        if startingTreatment == "parent":

            n_bins = exp.shape[0]
            n_genes = exp.shape[2]

            filename = f"path_info_{project}{_preprocess}_{model_name}.pickle"
            with open(os.path.join(input_path, filename), 'rb') as handle:
                path_info = pickle.load(handle, encoding="latin1")

            exps_parent = dict.fromkeys(unique_paths)
            counts_parent = dict.fromkeys(unique_paths)

            # for each path, find its parent's expression
            for i, path in enumerate(path_info):

                try:
                    exp_parent = exp[:, path2idx[path['parent_path']], :]  # use idx instead; not path"ID"
                    count_parent = count[:, path2idx[path['parent_path']]]
                except KeyError:
                    exp_parent = np.zeros((n_bins, n_genes))
                    count_parent = np.zeros((n_bins))

                exps_parent[path['ID']] = exp_parent
                counts_parent[path['ID']] = count_parent

            # append parent's expression to the path
            exps_with_parent = []
            counts_with_parent = []

            for i, cur_path in enumerate(unique_paths):
                cur_exp = exp[:, i, :]
                exp_parent = exps_parent[cur_path]

                exp_with_parent = np.concatenate([exp_parent, cur_exp], axis=0)

                cur_count = count[:, i]
                count_parent = counts_parent[cur_path]

                count_with_parent = np.concatenate([count_parent, cur_count], axis=0)

                exps_with_parent.append(exp_with_parent[:, None])
                counts_with_parent.append(count_with_parent[:, None])

            exp = np.concatenate(exps_with_parent, axis=1)  # bin size x number of paths x genes
            count = np.concatenate(counts_with_parent, axis=1)  # bin size x number of paths x genes

        elif startingTreatment == "smallerWindow":
            if n_lap > 1:
                smaller_n_lap = int(n_lap / 2)
            else:
                smaller_n_lap = 1

            # calculate move mean for a smaller window size
            _exp_sum = bn.move_sum(exp, window=smaller_n_lap, min_count=1, axis=0)
            _count_sum = bn.move_sum(count, window=smaller_n_lap, min_count=1, axis=0)
            exp_smaller = _exp_sum / _count_sum[:, :, None]

            if nan2zero:
                exp_smaller = np.nan_to_num(exp_smaller)

        else:
            pass

        del cell_paths
        del cell_times
        del cur_exp
        del cur_time
        del df_exp
        del count_bin
        del sum_bin
        del df_at
        del sum_full
        del count_full
        # collect unused variables
        gc.collect()

        # get expressions over a window
        # still hold the space for the first few elements where window size is not large enough
        # min_count=1: still calculate for the first few elements even though the length is small than n_lap
        # also, np.nan is ignored as long as 1 element is not nan
        _exp_sum = bn.move_sum(exp, window=n_lap, min_count=1, axis=0)
        _count_sum = bn.move_sum(count, window=n_lap, min_count=1, axis=0)
        exp = _exp_sum / _count_sum[:, :, None]

        if nan2zero:
            exp = np.nan_to_num(exp)

        if startingTreatment == "parent":
            exp = exp[n_bins:, :, :]  # remove parenting values

            del exps_with_parent
            del counts_with_parent
            del exp_with_parent
            del count_with_parent
            del exps_parent
            del counts_parent
            gc.collect()

        elif startingTreatment == "discard":
            exp = exp[n_lap:, :, :]

        elif startingTreatment == "smallerWindow":
            exp = exp[n_lap:, :, :]
            exp = np.concatenate([exp_smaller[smaller_n_lap:n_lap, :, :], exp], axis=0)

            del exp_smaller
            gc.collect()
        else:
            pass

        del _exp_sum
        del _count_sum
        # collect unused variables
        gc.collect()

        if align_type == "unaligned":

            # build pair-up sender and receiver expressions
            _sender = []
            _receiver = []

            for i in interaction_list:
                ligand = i[0]
                receptor = i[1]

                _sender.append(exp[:, :, gene_names.index(ligand)][:, :, None])
                _receiver.append(exp[:, :, gene_names.index(receptor)][:, :, None])

            _sender = np.concatenate(_sender, axis=2)
            _receiver = np.concatenate(_receiver, axis=2)

            del exp
            # collect unused variables
            gc.collect()

            # split into paths (TBD: faster if not split; check enisum)
            sender_exp = {}
            receiver_exp = {}

            # for each path, prepare a sender and a receiver expressions
            for i, cur_path in enumerate(unique_paths):
                sender_exp[cur_path] = _sender[:, i, :]
                receiver_exp[cur_path] = _receiver[:, i, :]

            # changed to process all pairs of paths?
            all_path_results = {}
            for cur_path in unique_paths:
                i = trasig.process_per_path(cur_path, sender_exp, receiver_exp)
                all_path_results.update(i)

            del sender_exp
            del receiver_exp
            del _sender
            del _receiver

            # collect unused variables
            gc.collect()

        else:
            observed_exist = False
            errors_path = os.path.join(output_path, f"{suffix}_0_errors.pickle")
            if os.path.isfile(errors_path):
               observed_exist = True

            if observed_exist and (_n == 0 or align_type == "aligned-fixed"):

                # build splines
                T = exp.shape[0]

                splines = dtw.get_splines(exp, T, gene_names, unique_paths)

                del exp
                # collect unused variables
                gc.collect()

                # read in the best alignment for the observed sample
                errors_path = os.path.join(output_path, f"{suffix}_0_errors.pickle")
                params_path = os.path.join(output_path, f"{suffix}_0_params.pickle")

                with open(errors_path, 'rb') as handle:
                    errors = pickle.load(handle)

                with open(params_path, 'rb') as handle:
                    align_params = pickle.load(handle)

            else:
                # conduct alignment
                if _n == 0:
                    save_errors = True
                    save_parameters = True

                    # save
                    errors_path = os.path.join(output_path, f"{suffix}_{_n}_errors.pickle")
                    params_path = os.path.join(output_path, f"{suffix}_{_n}_params.pickle")
                else:
                    save_errors = False
                    save_parameters = False

                    errors_path = None
                    params_path = None


                T = exp.shape[0]
                G = len(gene_names)

                splines = dtw.get_splines(exp, T, gene_names, unique_paths)

                del exp
                # collect unused variables
                gc.collect()

                T_REF_MIN = 0
                T_REF_MAX = T

                timepoints_aligned_ref, timepoints_aligned_sample, mask, align_params = \
                    dtw.get_params(G, T_REF_MIN, T_REF_MAX, save_parameters=save_parameters, params_path=params_path)

                errors = dtw.get_alignment_error(unique_paths, time2path, path2time, gene_pairs,
                                                 timepoints_aligned_ref, timepoints_aligned_sample,
                                                 mask, splines, save_errors=save_errors, errors_path=errors_path)

                del timepoints_aligned_ref
                del timepoints_aligned_sample
                del mask

                # collect unused variables
                gc.collect()


            # calculate score
            S, R = dtw.evaluate_splines_paths(time2path, path2time, errors, align_params, unique_paths,
                                              interaction_list, splines)

            # cap at 0 (remove negative values' impact on scores)
            S[S < 0] = 0
            R[R < 0] = 0

            all_path_results = trasig.process_all_paths(S, R)

            del splines
            del errors
            del align_params
            del S
            del R

            # collect unused variables
            gc.collect()

        # save results
        if save:
            filename = f"{suffix}_metrics_{_n}.pickle"
            data_file = os.path.join(output_path, filename)
            with open(data_file, 'wb') as handle:
                pickle.dump(all_path_results, handle, protocol=pickle.HIGHEST_PROTOCOL)

        # update counts
        if saveCounts:
            for pair, mets in all_path_results.items():
                for m in metrics:
                    _add = np.array((mets[m] >= results[pair][m]) * 1)
                    _counts[pair][m] += _add

        del all_path_results

        # collect unused variables
        gc.collect()

    return _counts


if __name__ == '__main__':
    # parse command-line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', required=True, default='../input/', help="string, folder to find inputs")
    parser.add_argument('-o', '--output', required=True, default='../output/', help="string, folder to put outputs")
    parser.add_argument('-d', '--project', required=True, help="string, project name")
    parser.add_argument('-g', '--preprocess', required=True, help="string, preprocessing steps applied to the "
                                                                  "data / project, default None", default="None")
    parser.add_argument('-b', '--modelName', required=True, help="string, name of the trajectory model")
    parser.add_argument('-t', '--listType', required=False,
                        default='ligand_receptor', help="string, optional, "
                                                        "interaction list type, default ligand_receptor")
    parser.add_argument('-l', '--nLap', required=False, default=20, help="integer, optional, "
                                                                         "sliding window size, default 20")
    parser.add_argument('-m', '--metric', required=False, default='dot', help="string, optional, "
                                                                              "scoring metric, default dot")
    parser.add_argument('-z', '--nan2zero', required=False, type=str2bool,
                        default=True, help="boolean, optional, if treat nan as zero, default True")
    parser.add_argument('-n', '--numPerms', required=False,
                        default=10000, help="integer, optional, number of permutations, default 10000")
    parser.add_argument('-p', '--multiProcess', required=False, type=str2bool,
                        default=True, help="boolean, optional, if use multi-processing, default True")
    parser.add_argument('-c', '--ncores', required=False, default=4, help="integer, optional, number of cores to use "
                                                                          "for multi-processing, default 4")
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
    
    if preprocess != "None":
        _preprocess = f"_{preprocess}"
    else:
        _preprocess = ""

    suffix = f"{project}_{list_type}{_preprocess}_{model_name}"

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

    suffix = f"{suffix}{_startingTreatment}_nlap_{n_lap}"

    # set parameters for multi-processing
    multi_processing = args.multiProcess
    ncores = int(args.ncores)

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

    # initialize trasig
    trasig = trasig(metrics, time2path, path2time, cell_times_o, cell_paths_o, unique_paths)

    # create output folder if not exist yet
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    # initialize alignment (if using alignment)
    if align_type != "unaligned":
        dtw = DTW(SMOOTH_PARAMETER, OVERLAP_THRESHOLD, SAMPLING_RATE, gene_pair_type, error_type, a_rate, b_rate)

        # get gene pairs for alignment
        if gene_pair_type == "interaction":
            gene_pairs = interaction_list
        else:
            # get alignment genes list
            filename = f"align_{gene_pair_type}_filtered_{project}{_preprocess}.pickle"
            with open(os.path.join(input_path, filename), 'rb') as handle:
                align_genes = pickle.load(handle)

            if len(align_genes) == 0:
                raise ValueError(f"No align genes available! Consider using 'interaction' or another set of "
                                 f"genes for alignment!")

            gene_pairs = list(zip(align_genes, align_genes))

    else:
        pass

    # load the result from the original data
    _n = 0

    filename = f"{suffix}_metrics_{_n}.pickle"
    data_file = os.path.join(output_path, filename)
    # check if file exists if not, run original
    if not os.path.isfile(data_file):
        _ = process_per_permutation([_n], save=True)

    # load results
    with open(data_file, 'rb') as handle:
        results = pickle.load(handle)

    # initialize to calculate p-value
    pair2counts = dict.fromkeys(results.keys())
    for pair, counts in pair2counts.items():
        pair2counts[pair] = dict.fromkeys(metrics)
        for m in metrics:
            pair2counts[pair][m] = np.repeat(0, len(results[pair][m]))

    # if just parallel among different paths, then slower than not parallel at all
    if multi_processing:
        # sequential due to too large memory consumption given the return variables
        _quotient = int(num_perms // (ncores - 1))
        _remainder = int(num_perms % (ncores - 1))

        # the remainder part normally contains very few jobs compared to the others
        # actually only use ncores - 1
        _groups = []

        idx_start = 1
        for i in range(ncores - 1):
            _groups.append(range(idx_start, _quotient + idx_start))
            idx_start += _quotient

        _groups.append(range(idx_start, _remainder + idx_start))

        print(f"Permutations grouped to different cores: {_groups}")

        MPRWork = MPR(process_per_permutation, _groups, ncores)
        Res = MPRWork.poolwork()
        del MPRWork  # delete the multi-threading worker to avoid memory leak

        # update the global variable
        for i in range(len(Res)):
            _results = Res[i]
            if len(_groups[i]) > 0:  # skips groups range(1,1)
                for pair, counts in pair2counts.items():
                    for m in metrics:
                        pair2counts[pair][m] += _results[pair][m]
        del _results
        del Res
    else:

        for _n in range(num_perms):
            print(_n)
            process_per_permutation(_n)

    # save
    suffix = f"{suffix}_{int(np.log10(num_perms))}"

    if align_type != "unaligned":
        suffix = f"{suffix}_{align_type}"
    else:
        pass

    filename = f"{suffix}_permutation_results.pickle"
    data_file = os.path.join(output_path, filename)

    with open(data_file, 'wb') as handle:
        pickle.dump(pair2counts, handle, protocol=pickle.HIGHEST_PROTOCOL)