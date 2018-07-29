# -*- coding: utf-8 -*-
"""
Execute CEA calculation
"""

import numpy as np
from scipy import interpolate
import pandas as pd
import matplotlib.pyplot as plt
import os, sys, glob, shutil
import re, copy
import tqdm
from subprocess import*
import warnings
from cea_exe import CEA_execute


class Read_datset:
    """
    Read_datset(fld_path, fexten="csv")
    
    Class to read and interpolate datasets calculated parameter with respect to every O/F and Pc
    
    Parameter
    ---------
    self.fld_path: string
        folder path containing dataset files

    self.fexten: string
        Defalut-value = "csv".
        Extension of dataset file. Default is CSV file.

    Class variable
    --------
    self.of: ndarray
        oxidizer to fuel ratio
    
    self.Pc: ndarray
        chamber pressure [MPa]

    """
    
    def __init__(self, fld_path, fexten="csv"):
        self.fld_path = fld_path
        self.fexten = fexten
        if os.path.exists(self.fld_path):
            flist = self.get_flist()
#            print(flist)
            init_fpath = os.path.join(self.fld_path, flist[0] +"."+self.fexten)
            dataframe = pd.read_csv(init_fpath, header=0, index_col=0, comment="#")
            self.of = np.asarray([float(i) for i in dataframe.index])
            self.Pc = np.asarray([float(i) for i in dataframe.columns])
        else:
            print("There is no such a dataset file/n{}".format(self.fld_path))
    
    def _read_csv_(self, param_name):
        """
        Read a csv-type-dataset files
        
        Parameter
        ---------
        param_name: string
            Parameter name which is a dataset file name \n
            e.g. "CSTAR", "GAMMAs", "T_c", "Cp_c"
        
        Return
        ------
        array: 2-ndarray, array.shape -> (of.size(), Pc.size())
            values of a calculated parameter with respect to every O/F and Pc
        """
        fpath = os.path.join(self.fld_path, param_name+"."+self.fexten)
        
        if os.path.exists(fpath):
            dataframe = pd.read_csv(fpath, header=0, index_col=0, comment="#")
            array = np.asarray(dataframe)
        else:
            print("There is no such a parameter \"{}\"\n".format(param_name))
            print("Please select a parameter from below list\n")
            print(self.get_flist())
            sys.exit(1)
        return(array)

    def get_flist(self):
        """
        Get dataset files list
        
        Return
        -------
        file_list: list
            list of csv files
        """
        split =  lambda r: os.path.splitext(r)[0] # get file name without extention
        file_list = [os.path.basename(split(r))  for r in glob.glob(self.fld_path + "/*.{}".format(self.fexten))]
        return(file_list)

    def gen_func(self, param_name, extraporate="linear"):
        """
        Generate function of calculated parameter with respect to O/F and Pc.
        Return a value after reading csv file and interpolate the data.
        
        Parameter
        ---------
        param_name: string
            Parameter name which is a dataset file name \n
            e.g. "CSTAR", "GAMMAs", "T_c", "Cp_c"
        
        extraporate: string; optional
            "exp": using f(x)=theta*numpy.exp(x-p)+q to extraporate the region out of database
            "exp2"
            "ln"
            "inverse"
            "power"
            "linear" Default; using linear function to extraporate

        Return
        ------
        func: function(of, Pc)
            A function which return a interpolated value (array-like)
        """
        array = self._read_csv_(param_name)
        func_interp = interpolate.interp2d(self.of, self.Pc, array.T, kind="cubic", bounds_error=False)
        def func(of, Pc):
            """Function to do interpolation and linear-extrapolation about the assigned database.
            Extrapolation is avalable only when assigned O/F is out of data-base range
            
            Parameter
            -----------
            of: float,
                O/F
            Pc: float,
                Chamber Pressure [Pa]
            
            Return
            ----------
            val: float
                interpolated or extrapolated value
            """
            Pc = Pc*1.0e-6
            cstr_array = func_interp(self.of, Pc)
            def extrapfunc_exp(of, a, b,diff, ddiff):
                theta = ddiff/diff
                p = np.log(diff/theta)/theta + a
                q = b - np.exp(theta*(a-p))
                val = np.exp(theta*(a-p))+q
                return(val)
            def extrapfunc_exp2(of, a,b, diff, ddiff):
                p = a - ddiff/diff
                theta = np.power(diff,2)/(ddiff*np.exp(ddiff/diff))
                q = b - theta*np.exp(a-p)
                val = theta*np.exp(a-p) + q
                return(val)
            def extrapfunc_ln(of, a,b, diff, ddiff):
                p = diff/ddiff + a
                theta = -np.power(diff, 2.0)/ddiff
                q = b - theta*np.log(a-p)
                val = theta*np.log(of - p) + q
                return(val)
            def extrapfunc_inverse(of, a,b, diff, ddiff):
                p = 2*diff/ddiff + a
                theta = -diff*np.power(a-p, 2.0)
                q = b - theta/(a-p)
                val = theta/(of-p) + q
                return(val)
            def extrapfunc_power(of, a, b, diff, ddiff, dddiff):
                phi = (dddiff*diff-2*np.power(ddiff,2))/(dddiff*diff-np.power(ddiff,2))
                p = a - (phi-1)*diff/ddiff
                theta = diff/(phi*np.power(a-p,phi-1))
                q = b - theta*np.power(a-p,phi)
                val = theta*np.power(of-p,phi)+q
                return(val)
            def extrapfunc_linear(of, a, b, diff):
                val= diff*(of-a) + b
                return(val)
            
            if of<self.of.min(): #when assigned O/F is smaller than minimum O/F of database
                diff_begin = (-3*cstr_array[0] +4*cstr_array[1] -cstr_array[2])/(2*(self.of[1]-self.of[0]))
                ddiff_begin = (2*cstr_array[0] -5*cstr_array[1] + 4*cstr_array[2] -cstr_array[3])/np.power((self.of[1]-self.of[0]),2.0)
#                dddiff_begin = 
                a = self.of.min()
                b = cstr_array[0]
                if extraporate=="exp":
                    val = extrapfunc_exp(of, a, b, diff_begin, ddiff_begin)
                elif extraporate=="exp2":
                    val = extrapfunc_exp2(of, a, b, diff_begin, ddiff_begin)
                elif extraporate=="ln":
                    val = extrapfunc_ln(of, a, b, diff_begin, ddiff_begin)
                elif extraporate=="inverse":
                    val = extrapfunc_inverse(of, a, b, diff_begin, ddiff_begin)
                elif extrapfunc_power=="power":
                    val = extrapfunc_power(of, a, b, diff_begin, ddiff_begin, dddiff_begin)
                elif extraporate=="linear":
                    val = extrapfunc_linear(of, a, b, diff_begin)
            elif self.of.max()<of: #when assigned O/F is larger than maximum O/F of database
                diff_end = (cstr_array[len(cstr_array)-3] -4* cstr_array[len(cstr_array)-2] +3*cstr_array[len(cstr_array)-1])/(2*(self.of[len(self.of)-1]-self.of[len(self.of)-2]))
                ddiff_end = (-2*cstr_array[len(cstr_array)-4] +4*cstr_array[len(cstr_array)-3] -5*cstr_array[len(cstr_array)-2] +2*cstr_array[len(cstr_array)-1])/np.power(self.of[len(self.of)-1]-self.of[len(self.of)-2], 2.0)
#                dddiff_end = 
                a = self.of.max()
                b = cstr_array[len(cstr_array)-1]
                if extraporate=="exp":
                    val = extrapfunc_exp(of, a, b, diff_end, ddiff_end)
                elif extraporate=="exp2":
                    val = extrapfunc_exp2(of, a, b, diff_end, ddiff_end)
                elif extraporate=="ln":
                    val = extrapfunc_ln(of, a, b, diff_end, ddiff_end)
                elif extraporate=="inverse":
                    val = extrapfunc_inverse(of, a, b, diff_end, ddiff_end)
                elif extrapfunc_power=="power":
                    val = extrapfunc_power(of, a, b, diff_end, ddiff_end, dddiff_end)
                elif extraporate=="linear":
                    val = extrapfunc_linear(of, a, b, diff_end)
            else: #when assigned O/F is with in the range of O/F
                    val = func_interp(of, Pc)[0]
            return(val)
        
        return(func)
        
    def plot(self, param_name, pickup_num):    
        """
        Plot graph about relationship of param to of and Pc
        
        Parameters
        ----------
        param_name: string
            Parameter name which is a dataset file name \n
            e.g. "CSTAR", "GAMMAs", "T_c", "Cp_c"

        pickup_num: int
            The number of "Pc" picked up to draw a graph 
        """
        if pickup_num > len(self.Pc):
            pass
        else:
            Pc_nlm = (self.Pc[-1]-self.Pc[0])/(pickup_num-1)
            pick_idx = lambda i: np.abs(np.asarray(self.Pc)-(Pc_nlm*i + self.Pc[0])).argmin()
            idx = [pick_idx(i) for i in range(pickup_num)]

        array = self._read_csv_(param_name)        
        plt.rcParams["font.family"] = "Times New Roman"
        plt.rcParams["font.size"] = 17
        fig = plt.figure(figsize=(8,6))
        ax = fig.add_subplot(111)
        for i in idx:
            ax.plot(self.of, array[:,i], label=r"$P_c$ = {} MPa".format(self.Pc[i]))
        ax.legend(loc="best", fontsize=16)
        ax.set_xlabel(r"$O/F$")
        ax.set_ylabel("${}$".format(param_name))


if __name__ == "__main__":

    dbfld_path = os.path.join("cea_db", "sample", "csv_database")
    inst2 = Read_datset(dbfld_path)
    func = inst2.gen_func("CSTAR", extraporate="linear")
    func(100, 1.0e+6)
    of_range = np.arange(0.01, 100, 0.1)
    plt.plot(of_range, np.array([func(of, 1.0e+6) for of in of_range]))
