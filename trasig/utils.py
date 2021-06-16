# ------ helper classes / funcitons --------
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
    Helper to pass arguements.

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