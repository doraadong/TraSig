#!/usr/bin/env python

from collections import defaultdict

import numpy as np
from scipy import interpolate
import pickle

from trasig.utils import getPairedPaths
"""

Modified based on Jun Ding's alignment code for gene expression profiles and Jose Lugo-Martinez's alignment code. 

Ref: Ruiz-Perez, D., Lugo-Martinez, J., Bourguignon, N., Mathee, K., Lerner, B., Bar-Joseph, Z., 
& Narasimhan, G. (2021). Dynamic Bayesian Networks for Integrating Multi-omics Time Series Microbiome Data. 
Msystems, 6(2), e01105-20.

"""


def warpFunction(a, b, s, warpType='linear'):
    if warpType == 'exponential':
        return np.exp((s - b) / a)
    else:
        return (s - b) / a


def warpFunctionInverse(a, b, t, warpType='linear'):
    if warpType == 'exponential':
        return a * np.log(t) + b
    else:
        return (a * t) + b

class DTW:

    def __init__(self, SMOOTH_PARAMETER, OVERLAP_THRESHOLD, SAMPLING_RATE, gene_pair_type, error_type, a_rate, b_rate):
        self.SMOOTH_PARAMETER = SMOOTH_PARAMETER
        self.OVERLAP_THRESHOLD = OVERLAP_THRESHOLD
        self.SAMPLING_RATE = SAMPLING_RATE
        self.gene_pair_type = gene_pair_type
        self.error_type = error_type
        self.a_rate = a_rate
        self.b_rate = b_rate

        self.A_RANGE = [0.5, 2.0]
        self.B_RANGE = [-50, 50]

    def get_splines(self, exp, T, gene_names, unique_paths, degree=3, save_splines=False, splines_path=None):
        """

        Get splines per edge and per gene.

        Parameters
        ----------
        exp: numpy array,
            number of time points x number of edges x number of genes, expression values
        T: integer,
            number of time points
        gene_names: list of strings,
            gene symbols
        unique_paths: numpy array,
            number of edges x 1
        degree: integer
            degree of splines, default 3
        save_splines: boolean
            if save the fitted spline parameters, default False
        splines_path: string
            path to save the spline parameters, default None

        Returns
        -------
        splines: dictionary with edge as the 1st level key and gene as the 2nd level,
            fitted splines parameters (tuple)

        """
        timepoints = np.arange(T)

        splines = defaultdict(dict)
        for i, cur_path in enumerate(unique_paths):
            for j, cur_gene in enumerate(gene_names):
                _exp = exp[:, i, j]
                tck = interpolate.splrep(timepoints, _exp, k=degree, s=self.SMOOTH_PARAMETER)
                splines[cur_path][cur_gene] = tck

        if save_splines:
            # save fitted splines
            with open(splines_path, 'wb') as handle:
                pickle.dump(splines, handle, protocol=pickle.HIGHEST_PROTOCOL)

        return splines


    def get_params(self, G, T_REF_MIN, T_REF_MAX, save_parameters=False, params_path=None):

        """

        Get all possible combinations of parameter values for the alignment.

        Parameters
        ----------
        G: integer
            number of genes
        T_REF_MIN: integer
            minimum time in the reference
        T_REF_MAX: integer
            maximum time in the reference
        save_parameters: boolean
            if save the candidate values for parameters, default False
        params_path: string
            path to save the alignment parameters, default None

        Returns
        -------
        timepoints_aligned_ref: numpy array, number of aligned time points x number of combinations of parameters
            aligned time points in reference
        timepoints_aligned_sample: numpy array, number of aligned time points x number of combinations of parameters
            aligned time points warped to sample
        mask: numpy array, number of aligned time points x number of combinations of parameters
            binary matrix with 0 for padded zeros, 1 otherwise
        align_params: list
            each element is a list [a, b, alpha, beta, overlap]

        """

        ## get alignment parameters

        # the minimum time in the sample edge same as the reference edge (also for the max)
        T_SAMPLE_MIN = T_REF_MIN
        T_SAMPLE_MAX = T_REF_MAX

        max_length = 0  # may be larger than the original time points given different sampling rate
        align_params = []
        for a in np.arange(self.A_RANGE[0], self.A_RANGE[1] + self.a_rate,  self.a_rate):
            for b in np.arange(self.B_RANGE[0], self.B_RANGE[1] + self.b_rate, self.b_rate):
                alpha = max(T_REF_MIN, warpFunctionInverse(a, b, T_SAMPLE_MIN))
                beta = min(T_REF_MAX, warpFunctionInverse(a, b, T_SAMPLE_MAX))
                overlap = (beta - alpha) / (T_REF_MAX - T_REF_MIN)
                if overlap > self.OVERLAP_THRESHOLD and alpha < beta:

                    _aligned = np.arange(alpha, (beta + self.SAMPLING_RATE), self.SAMPLING_RATE)

                    if len(_aligned) > max_length:
                        max_length = len(_aligned)

                    align_params.append([a, b, alpha, beta, overlap])

        print(f"Max length of the aligned time points is {max_length}")

        # save alignment parameters
        if save_parameters:
            with open(params_path, 'wb') as handle:
                pickle.dump(align_params, handle, protocol=pickle.HIGHEST_PROTOCOL)

        ## get all possible time points
        timepoints_aligned_ref = []
        timepoints_aligned_sample = []

        mask = []

        for a in np.arange(self.A_RANGE[0], self.A_RANGE[1] + self.a_rate,  self.a_rate):
            for b in np.arange(self.B_RANGE[0], self.B_RANGE[1] + self.b_rate, self.b_rate):
                alpha = max(T_REF_MIN, warpFunctionInverse(a, b, T_SAMPLE_MIN))
                beta = min(T_REF_MAX, warpFunctionInverse(a, b, T_SAMPLE_MAX))
                overlap = (beta - alpha) / (T_REF_MAX - T_REF_MIN)
                if overlap > self.OVERLAP_THRESHOLD and alpha < beta:
                    _ref = np.arange(alpha, (beta + self.SAMPLING_RATE), self.SAMPLING_RATE)
                    _sample = warpFunction(a, b, _ref)

                    pad_length = max_length - len(_ref)
                    _ref_p = np.pad(_ref, (0, pad_length), 'constant')
                    _sample_p = np.pad(_sample, (0, pad_length), 'constant')

                    _mask = np.array(list(np.ones(len(_ref))) + list(np.zeros(pad_length)))

                    timepoints_aligned_ref.append(_ref_p)
                    timepoints_aligned_sample.append(_sample_p)
                    mask.append(_mask)

        timepoints_aligned_ref = np.array(timepoints_aligned_ref).T
        timepoints_aligned_sample = np.array(timepoints_aligned_sample).T
        mask = np.array(mask).T

        return timepoints_aligned_ref, timepoints_aligned_sample, mask, align_params

    def get_alignment_error(self, unique_paths, time2path, path2time, gene_pairs, timepoints_aligned_ref,
                            timepoints_aligned_sample, mask, splines, save_errors=False, errors_path=None):
        """

        Find alignment errors for all pairs of edges.

        Parameters
        ----------
        unique_paths: numpy array,
            number of edges x 1
        time2path: dictionary
            time as key and list of edges / paths as value
        path2time: dictionary
            edge / path as key and time as value
        gene_pairs: list of tuples
            pairs of genes used for alignment
        timepoints_aligned_ref: numpy array, number of aligned time points x number of combinations of parameters
            aligned time points in reference
        timepoints_aligned_sample: numpy array, number of aligned time points x number of combinations of parameters
            aligned time points warped to sample
        mask: numpy array, number of aligned time points x number of combinations of parameters
            binary matrix with 0 for padded zeros, 1 otherwise
        splines: dictionary with edge as the 1st level key and gene as the 2nd level,
            fitted splines parameters (tuple)
        save_errors: boolean
            if save the errors (dictionary, each value as number of combinations of parameters x 1) for all pairs of
            edges, default False
        errors_path: string
            path to save the errors, default None


        Returns
        -------
        errors: dictionary
            f"{cur_path}_{other_path}" as key and a numpy array (umber of combinations of parameters x 1) as value,
            alignment errors for all pairs of edges

        """
        errors = defaultdict(dict)

        for cur_path in unique_paths:
            paired_paths = getPairedPaths(cur_path, time2path, path2time)
            for other_path in paired_paths:
                errors[f"{cur_path}_{other_path}"] = self.get_alignment_error_pairwise(cur_path, other_path, gene_pairs,
                                                                                       timepoints_aligned_ref,
                                                                                       timepoints_aligned_sample, mask,
                                                                                       splines)
        if save_errors:
            # save alignment errors
            with open(errors_path, 'wb') as handle:
                pickle.dump(errors, handle, protocol=pickle.HIGHEST_PROTOCOL)

        return errors

    def get_alignment_error_pairwise(self, ref_path, sample_path, gene_pairs, timepoints_aligned_ref,
                                     timepoints_aligned_sample, mask, splines, save_pairwise_errors=False,
                                     pairwise_errors_path=None):

        """

        Get alignment errors for a pair of edges.


        Parameters
        ----------
        ref_path: string
            name of the reference edge
        sample_path: string
            name of the sample edge
        gene_pairs: list of tuples
            pairs of genes used for alignment
        timepoints_aligned_ref: numpy array, number of aligned time points x number of combinations of parameters
            aligned time points in reference
        timepoints_aligned_sample: numpy array, number of aligned time points x number of combinations of parameters
            aligned time points warped to sample
        mask: numpy array, number of aligned time points x number of combinations of parameters
            binary matrix with 0 for padded zeros, 1 otherwise
        splines: dictionary with edge as the 1st level key and gene as the 2nd level,
            fitted splines parameters (tuple)
        save_pairwise_errors: boolean
            if save the errors (number of gene pairs x number of combinations of parameters)
            between the current edge pairs, default False
        pairwise_errors_path: string
            path to save the pairwise errors, default None


        Returns
        -------
        errors: numpy array, number of combinations of aligned parameters x 1
            overall error for each parameter combination


        Notes
        -----
        Currently not consider weights for different gene pairs.

        """

        num_loops = len(gene_pairs)
        error_per_genepair = defaultdict(dict)

        # TBD: parellel?
        for i in range(num_loops):
            try:
                cur_pair = gene_pairs[i]
            except IndexError:
                print(i)

            gene_a, gene_b = cur_pair

            _ref = interpolate.splev(timepoints_aligned_ref, splines[ref_path][gene_a])
            _sample = interpolate.splev(timepoints_aligned_sample, splines[sample_path][gene_b])

            _ref *= mask
            _sample *= mask

            if self.error_type == "MSE":
                _error = ((_ref - _sample) ** 2).sum(axis=0)
            elif self.error_type == "cosine":
                _nominator = np.einsum('ij,ij->j', _ref, _sample)
                _denominator = np.linalg.norm(_ref, axis=0) * np.linalg.norm(_sample, axis=0)
                _error = 1 - _nominator / _denominator
            elif self.error_type == "corr":
                _ref = _ref - np.mean(_ref, axis=0)
                _sample = _sample - np.mean(_sample, axis=0)

                _nominator = np.einsum('ij,ij->j', _ref, _sample)
                _denominator = np.linalg.norm(_ref, axis=0) * np.linalg.norm(_sample, axis=0)
                _error = 1 - _nominator / _denominator
            else:
                raise NotImplementedError(f"Now only support MSE, corr or cosine.")

            error_per_genepair[i] = _error

        if save_pairwise_errors:
            # save alignment errors
            with open(pairwise_errors_path, 'wb') as handle:
                pickle.dump(error_per_genepair, handle, protocol=pickle.HIGHEST_PROTOCOL)

        errors_total = np.nan_to_num(np.array(list(error_per_genepair.values())))
        errors = errors_total.sum(axis=0)

        return errors


    def evaluate_splines_paths(self, time2path, path2time, errors, align_params, unique_paths, interaction_list,
                               splines, NUM_SCORING = 100):

        """

        Estimate expression values using fitted splines.

        Parameters
        ----------
        time2path: dictionary
            time as key and list of edges / paths as value
        path2time: dictionary
            edge / path as key and time as value
        errors: dictionary
            f"{cur_path}_{other_path}" as key and a numpy array (umber of combinations of parameters x 1) as value,
            alignment errors for all pairs of edges
        align_params: list
            all combinations of possible parameter values, each element is a list [a, b, alpha, beta, overlap]
        unique_paths: numpy array,
            number of edges x 1
        interaction_list: list of tuples
            pairs of genes used for assessing interactions
        splines: dictionary with edge as the 1st level key and gene as the 2nd level,
            fitted splines parameters (tuple)
        NUM_SCORING: integer
            number of samples to take from the aligned time line

        Returns
        -------
        SS: numpy array, number of time points x number of pair of edges x number of interaction pairs
            values in the sender
        RR: numpy array, number of time points x number of pair of edges x number of interaction pairs
            values in the receiver

        """

        _RR = []
        _SS = []

        for cur_path in unique_paths:
            paired_paths = getPairedPaths(cur_path, time2path, path2time)

            for other_path in paired_paths:

                align_param = align_params[np.argmin(errors[f"{cur_path}_{other_path}"])]

                a = align_param[0]
                b = align_param[1]

                start_ref = align_param[2]
                end_ref = align_param[3]

                time_ref = np.linspace(start=start_ref, stop=end_ref, num=NUM_SCORING)
                time_sample = warpFunction(a, b, time_ref)

                s = []
                r = []

                # for each ligand - receptor pair
                for i in interaction_list:
                    ligand = i[0]
                    receptor = i[1]

                    spl = splines[cur_path][ligand]
                    _s = interpolate.splev(time_ref, spl)

                    spl = splines[other_path][receptor]
                    _r = interpolate.splev(time_sample, spl)

                    s.append(_s[:, None])
                    r.append(_r[:, None])

                S = np.concatenate(s, axis=1)
                R = np.concatenate(r, axis=1)

                _RR.append(R[:, None])
                _SS.append(S[:, None])

        RR = np.concatenate(_RR, axis=1)
        SS = np.concatenate(_SS, axis=1)

        return SS, RR
