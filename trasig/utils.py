"""

Helper classes / functions.

9-17-21: Move getPairedPaths to utils.

"""


def getPairedPaths(cur_path, time2path, path2time, same_time=True):
    """
    Get edges from the same time or (closed time if same_time = False).

    Parameters
    ----------
    cur_path: integer
        name of the current edge / path
    time2path: dictionary
        time as key and list of edges / paths as value
    path2time: dictionary
        edge / path as key and time as value
    same_time: boolean
        if only select edge / path belonging to the same time


    Returns
    -------
    paired_paths: list
        edge / paths closed in time to the current edge / path

    """

    time = path2time[cur_path]

    paired_paths = []
    for k, v in time2path.items():
        if same_time:
            if k == time:
                paired_paths += time2path[k]
        else:
            if k <= time + 1 and k >= time - 1:  # include 1 sampling time before and after
                paired_paths += time2path[k]

    return paired_paths


from multiprocessing import Pool
import argparse

class MPR:
    """
    Adapted from:

    Author: Jun Ding
    Project: SCDIFF2
    Ref: Ding, J., Aronow, B. J., Kaminski, N., Kitzmiller, J., Whitsett, J. A., & Bar-Joseph, Z.
    (2018). Reconstructing differentiation networks and their regulation from time series
    single-cell expression data. Genome research, 28(3), 383-395.

    """

    def __init__(self, Func, DataList, ncores=None):
        self.Func = Func
        self.DataList = DataList
        self.ncores = ncores

    def poolwork(self):
        pool = Pool(processes=self.ncores, maxtasksperchild=1)
        Res = pool.map_async(self.Func, self.DataList)
        pool.close()
        pool.join()
        Res = Res.get()
        del pool
        return Res



def str2bool(v):
    """
    Helper to pass boolean arguements.

    Extracted from: https://stackoverflow.com/questions/15008758/parsing-boolean-values-with-argparse
    Author: @Maxim

    """

    if isinstance(v, bool):
       return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')