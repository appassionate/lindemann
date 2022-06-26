import os
import MDAnalysis as mda
from MDAnalysis.analysis.base import (AnalysisBase,
                                      AnalysisFromFunction,
                                      analysis_class)

from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis

from MDAnalysis.lib.distances import distance_array

import numpy as np
# import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import ticker

from typing import Type
import numbers
from MDAnalysis.lib.log import ProgressBar


import logging
logger = logging.getLogger(__name__)

# import need correct

from MDAnalysis.lib.util import (fixedwidth_bins, iterable, asiterable,
                                 deprecate,)
import warnings
from MDAnalysis.analysis.density import Density



class ParallelAnalysisBase(AnalysisBase):

    def __init__(self, atomgroup, verbose=True):
        trajectory = atomgroup.universe.trajectory

        #np.seterr(divide='ignore', invalid='ignore')
        super(ParallelAnalysisBase, self).__init__(trajectory, verbose=verbose)

        # trajectory value initial
        self.ag = atomgroup
        self.ag_natoms = len(self.ag)


        #parallel value initial
        self.para = None
        self._para_region = None
        #MORE NEED TO CORRECT


    def _parallel_init(self, *args, **kwargs):

        start = self._para_region.start
        stop = self._para_region.stop
        step = self._para_region.step
        self._setup_frames(self._trajectory, start, stop, step)
        self._prepare()

    def run(self, start=None, stop=None, step=None, verbose=None):

        #self._trajectory._reopen()
        if verbose == True:
            print(" ", end='')
        super().run(start, stop, step, verbose)### will be problem for conclude operation
        
        if self.para:
            block_result = self._para_block_result()
            if block_result == None:
                raise ValueError("in parallel, block result has not been defined or no data output!")
            #logger.info("block_anal finished.")
            return block_result

    def _para_block_result(self):
        # data need to be transformed, which is usually relative to values in prepare() method.

        return None


    def _parallel_conclude(self, rawdata):

        raise NotImplementedError("ParallelAnalysisBase is father class")


 

class LindemannAna(ParallelAnalysisBase):

    def __init__(self, atomgroup, verbose=True, cell=None, last_only=False):
        
        np.seterr(divide='ignore', invalid='ignore',) #might affect global variables

        self.cell = cell
        self.last_only = last_only
        
        super(LindemannAna, self).__init__(atomgroup, verbose=verbose)
        
    def _prepare(self):

        natoms = self.ag_natoms
        
        self.array_mean = np.zeros((natoms, natoms))
        self.array_var = np.zeros((natoms, natoms))
        
        self.series_lindex = np.empty((0, natoms, natoms), dtype=float)
        self.series_mean = np.empty((0, natoms, natoms), dtype=float)
        self.series_var = np.empty((0, natoms, natoms), dtype=float)
        
    def _single_frame(self):
        
        tmp = self.ag.positions.copy()
        array_distance = distance_array(tmp, tmp, box=self.cell)

        #array_distance = distance_array(self.ag.positions, self.ag.positions, box=self.cell)
        array_mean = self.array_mean
        array_var = self.array_var

        for i in range(self.ag_natoms):
            for j in range(i + 1, self.ag_natoms):
                xn = array_distance[i, j]
                mean = array_mean[i, j]
                var = array_var[i, j]
                delta = xn - mean

                # update mean
                self.array_mean[i, j] = mean + delta / (self._frame_index + 1)
                # update variance
                array_var[i, j] = var + delta * (xn - array_mean[i, j])

        for i in range(self.ag_natoms):
            for j in range(i + 1, self.ag_natoms):
                self.array_mean[j, i] = self.array_mean[i, j]
                self.array_var[j, i] = self.array_var[i, j]
        
        if not self.last_only:
            
            self.series_mean = np.r_[self.series_mean, [array_mean]]
            self.series_var = np.r_[self.series_var, [array_var]]
        
        self.array_mean = array_mean
        self.array_var = array_var
                
    def _conclude(self):
        
        #lindemann_indices = self.series_lindex
        nframes = len(self.series_var)
        
        #for last_only
        if self.last_only:
            self.last_atom_lindex = np.nanmean(np.sqrt(self.array_var/self.n_frames)/self.array_mean, axis=1)
            self.last_frame_lindex = np.mean(self.last_atom_lindex)
            self.results = (self.last_frame_lindex, self.last_atom_lindex)


        ## process var and mean to "raw" lindex
        self.series_lindex = np.divide(np.sqrt(np.divide(self.series_var, np.arange(1,nframes+1).reshape(nframes,1,1))), self.series_mean)
        
        #process "raw" series lindex
        self.per_atom_lindex = np.array([np.nanmean(i, axis=1) for i in self.series_lindex])
        self.per_frame_lindex = np.array([np.mean(np.nanmean(i, axis=1)) for i in self.series_lindex])
        
        self.results = (self.per_frame_lindex, self.per_atom_lindex)
        

    def _para_block_result(self):
        #which is correspond to _parallel_conclude rawdata
        if self.last_only:
            return self.array_mean, self.array_var, (self.start, self.stop, self.step)
        else:
            return self.series_mean, self.series_var
    
    
    def _last_only_parallel_conclude(self, rawdata, *args, **kwargs):
        
        mean_results = ([ block_result[0] for block_result in rawdata])
        var_results = ([ block_result[1] for block_result in rawdata])
        range_mark_results = ([ block_result[2] for block_result in rawdata])
        
        n_atoms = len(mean_results[0])
        n_block = len(range_mark_results)

        sample_lens = []
        for range_mark in range_mark_results: 
            block = range(range_mark[0], range_mark[1],range_mark[2])
            sample_lens.append(len(block))
        n_samples = sum(sample_lens)

        mean_results = ([ block_result[0] for block_result in rawdata])
        var_results = ([ block_result[1] for block_result in rawdata])

        M = np.zeros([n_atoms,n_atoms])
        V = np.zeros([n_atoms,n_atoms])

        for j in range(n_atoms):
            for k in range(j+1, n_atoms):
                for i in range(n_block):
                    M[j,k] = M[j,k] + mean_results[i][j,k]*sample_lens[i]/n_samples
                for i in range(n_block):
                    V[j,k] = V[j,k] + var_results[i][j,k] +sample_lens[i]*(M[j,k] - mean_results[i][j,k])**2

        for i in range(n_atoms):
            for j in range(i+1, n_atoms):
                M[j, i] = M[i, j]
                V[j, i]  = V[i, j]
        
        self.array_mean = M
        self.array_var = V
        

        #conclude
        self._conclude()
    
    
    #@jit(nopython=True)
    def _parallel_conclude(self, rawdata, *args, **kwargs):
        # future for map results process,might be static, without it, the parallel_exec will be interrupt.
        # conclude section: align the value of mean and var to serial execution
        
        ###new update###
        import time
        time_1 = time.time()
        
        self._parallel_init()

        if self.last_only:
            self._last_only_parallel_conclude(rawdata)
            return True
        
        
        mean_raw_data = [i[0] for i in rawdata]
        var_raw_data = [i[1] for i in rawdata]
        
        whole_true_var = [] #using numpy latter
        whole_true_mean = []

        for num, (slice_var, slice_mean) in enumerate(zip(var_raw_data, mean_raw_data)):

            latter_tail_mean = [ _slice_mean[-1] for _slice_mean in mean_raw_data[:num]]
            latter_tail_var = [ _slice_var[-1] for _slice_var in var_raw_data[:num]]
            latter_slice_frames = [ len(_slice_mean) for _slice_mean in mean_raw_data[:num]]
            latter_frames_sum = int(np.sum(latter_slice_frames))

            for k in range(len(slice_mean)):

                cur_frame_nums = latter_frames_sum + (k+1)
                sample_len = np.array(latter_slice_frames + [k+1])
                ratio_mean = np.array(latter_slice_frames + [k+1])/cur_frame_nums
                #calc mean
                tail_mean = np.array(latter_tail_mean + [slice_mean[k]])
                true_mean = [tail_mean[idx]*ratio_mean[idx] for idx in range(len(tail_mean))]
                true_mean = np.sum(true_mean, axis=0)
                
                #calc var
                tail_var = np.array(latter_tail_var + [slice_var[k]])   
                tail_unknown_factors = np.r_[np.square(true_mean - tail_mean)[:num], [np.square(true_mean - slice_mean[k])]]
                #val_1 =  np.array([tail_var[idx]*1 for idx in range(len(tail_var))])
                val_1 = tail_var
                val_1 = np.sum(val_1, axis=0)
                val_2 =  np.array([tail_unknown_factors[idx]*sample_len[idx] for idx in range(len(tail_unknown_factors))])
                val_2 = sum(val_2, )
                true_var = val_1 + val_2 #sum to var
                
                #add to list
                whole_true_mean.append(true_mean)
                whole_true_var.append(true_var)
        

        whole_true_var = np.array(whole_true_var)  #or true mean .r_?
        whole_true_mean = np.array(whole_true_mean)
        
        self.series_mean = whole_true_mean
        self.series_var = whole_true_var

        self._conclude()

        time_2 = time.time()
        delta_time = time_2 - time_1
        
        return "FINISH PARA CONCLUDE, USING {:.4f} s".format(delta_time)
