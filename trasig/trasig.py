
import numpy as np

from trasig.utils import getPairedPaths

class trasig:

    """

    Updates log:
    9-17-21: Add functions calculating scores for all path pairs altogehter (instead of process each pair at a time);
            these functions related to calculating the scores using aligned profiles

    """
    def __init__(self, metrics,  time2path, path2time, cell_times_o, cell_paths_o, unique_paths):
        self.metrics = metrics
        self.time2path = time2path
        self.path2time = path2time
        self.cell_times_o = cell_times_o
        self.cell_paths_o = cell_paths_o
        self.unique_paths = unique_paths

    def paired_dot(self, S, R, normalized=True):
        """
        Calculate dot product for each column in S and its paired column
        (same index) in R.

        Parameters
        ----------
            S: numpy array, number of time points x number of interaction pairs
                values in the sender
            R: numpy array, number of time points x number of interaction pairs
                values in the receiver
            normalized: boolean
                normalize by number of time points

        """
        if normalized:
            return np.einsum('ij,ij->j', S, R) / S.shape[0]
        else:
            return np.einsum('ij,ij->j', S, R)

    def paired_cov(self, S, R, normalized=True):
        """
        Calculate covariance estimate for each column in S and its paired column
        (same index) in R.

        Parameters:
        -----------
            normalized: boolean, normalize by number of rows (samples); biased
                        estimate for population covariance

        """
        _S = S - np.mean(S, axis=0)
        _R = R - np.mean(R, axis=0)
        return self.paired_dot(_S, _R, normalized)

    def paired_dot_paths(self, S, R, normalized=True):
        """
        Calculate dot product for each element defined by both the 2nd and 3rd axes in S and its paired element
        (same index) in R.

        Parameters
        ----------
            S: numpy array, number of time points x number of pair of edges x number of interaction pairs
                values in the sender
            R: numpy array, number of time points x number of pair of edges x number of interaction pairs
                values in the receiver
            normalized: boolean
                normalize by number of time points

        """
        if normalized:
            return np.einsum('ijk,ijk->jk', S, R) / S.shape[0]
        else:
            return np.einsum('ijk,ijk->jk', S, R)


    def paired_cov_paths(self, S, R, normalized=True):
        """
        Calculate covariance estimate for each element defined by both the 2nd and 3rd axes in S and its paired element
        (same index) in R.

        Parameters:
        -----------
            normalized: boolean, normalize by number of rows (samples); biased
                        estimate for population covariance

        """
        _S = S - np.mean(S, axis=0)
        _R = R - np.mean(R, axis=0)
        return self.paired_dot_paths(_S, _R, normalized)


    def generate_permutation(self, i, permute_time=True):
        """
        Generate ith permutation by permutating the cluster labels and the cell time assignment while
        keeping the original distribution of time conditional on label and the marginal distribution
        of labels (i.e. the fraction of each time in a label and the fraction of a label) unchanged.

        """
        np.random.seed(i)

        # generate a permutated cell assignment (path & time)
        permutated_paths = np.random.permutation(self.cell_paths_o)
        permutated_times = self.cell_times_o  # initialized as the original

        if permute_time:
            for i, cur_path in enumerate(self.unique_paths):
                # get original cell time assignment for the path
                cur_times = self.cell_times_o[self.cell_paths_o == cur_path]
                permutated_times[permutated_paths == cur_path] = np.random.permutation(cur_times)

        return permutated_paths, permutated_times

    def process_all_paths(self, S, R, round_values=True, round_decimals=10):

        """

        Parameters
        ----------
            S: numpy array, number of time points x number of pair of edges x number of interaction pairs
                values in the sender
            R: numpy array, number of time points x number of pair of edges x number of interaction pairs
                values in the receiver
            round_values: boolean
                if round values to some decimals (default 10); this is recommended given that output from bottleneck
                move_sum sometimes gives negative values very closed to 0 (e.g. -4e-17); these values will make the dot
                product also negative and thus make value >= 0 (if real dot product being 0) as false

        Returns
        -------
            all_path_results: dictionary of 2 layers
                1st layer indexed by cur_path_other_path, 2nd layer indexed
                by metric.

        """
        results = {}
        for metric in self.metrics:
            if metric == 'dot':
                result = self.paired_dot_paths(S, R)
            elif metric == 'cov':
                result = self.paired_cov_paths(S, R)
            else:
                raise NotImplementedError(f"Now only support dot or cov.")

            if round_values:
                results[metric] = np.round(result, round_decimals)
            else:
                results[metric] = result

        all_path_results = {}
        i = 0
        for cur_path in self.unique_paths:
            paired_paths = getPairedPaths(cur_path, self.time2path, self.path2time)
            for other_path in paired_paths:
                results_per_pair = dict.fromkeys(self.metrics)
                for metric in self.metrics:
                    results_per_pair[metric] = results[metric][i]
                all_path_results[f"{cur_path}_{other_path}"] = results_per_pair
                i += 1

        return all_path_results

    def process_per_path(self, cur_path, sender_exp, receiver_exp, round_values=True, round_decimals=10):
        """
        Parameters
        ----------
            round_values: boolean
                if round values to some decimals (default 10); this is recommended given that output from bottleneck
                move_sum sometimes gives negative values very closed to 0 (e.g. -4e-17); these values will make the dot
                product also negative and thus make value >= 0 (if real dot product being 0) as false

        Returns
        -------
            all_results: dictionary of 2 layers
                1st layer indexed by cur_path_other_path, 2nd layer indexed
                by metric.


        """
        S = sender_exp[cur_path]
        paired_paths = getPairedPaths(cur_path, self.time2path, self.path2time)

        all_results = {}
        for other_path in paired_paths:
            R = receiver_exp[other_path]

            results = {}
            for metric in self.metrics:
                if metric == 'dot':
                    result = self.paired_dot(S, R)
                elif metric == 'cov':
                    result = self.paired_cov(S, R)
                else:
                    raise NotImplementedError(f"Now only support dot or cov.")

                if round_values:
                    results[metric] = np.round(result, round_decimals)
                else:
                    results[metric] = result

            all_results[f"{cur_path}_{other_path}"] = results

        return all_results




