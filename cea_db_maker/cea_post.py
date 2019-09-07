# -*- coding: utf-8 -*-
"""
Read and utilize the result of CEA calculation
"""

import numpy as np
from scipy import interpolate
import pandas as pd
import matplotlib.pyplot as plt
import os, sys, glob


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
            def cea_exe(of, Pc):
                """ Using single execute of CEA instead of not using extraporation 
                
                Parameter
                --------
                of: float
                    O/F
                Pc: float
                    chamber pressure [MPa]
                """
                pass
            
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
                if extraporate == False:
                    val = cea_exe(of, Pc)
                elif extraporate=="exp":
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
                if extraporate == False:
                    val = cea_exe(of, Pc)
                elif extraporate=="exp":
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
        
    def plot(self, param_name, of_range, Pc_plot):    
        """
        Plot graph about relationship of param to of and Pc
        
        Parameters
        ----------
        param_name: string
            Parameter name which is a dataset file name \n
            e.g. "CSTAR", "GAMMAs", "T_c", "Cp_c"
        
        of_range: list
            [minimum of, maximum of]
        
        Pc_plot: 1-ndarray
            This array contain the chamber pressure [Pa] which you want to plot a graph
        """
        func = self.gen_func(param_name, extraporate="linear")
        plt.rcParams["font.family"] = "Times New Roman"
        plt.rcParams["font.size"] = 17
        fig = plt.figure(figsize=(8,6))
        ax = fig.add_subplot(111)
        for Pc in Pc_plot:
            ax.plot(self.of, np.array([func(of,Pc) for of in self.of]), label=r"$P_c$ = {} MPa".format(round(Pc*1.0e-6, 2)))
        ax.legend(loc="best", fontsize=16)
        ax.set_xlabel(r"$O/F$")
        ax.set_ylabel("${}$".format(param_name))
        ax.set_xlim(of_range[0], of_range[1])
        imgf_name = os.path.join(self.fld_path, param_name + ".jpg")
        plt.savefig(imgf_name)


if __name__ == "__main__":

    while(True):
        fld_name = input("Input Folder Name (e.g. \"O2+PMMA\")\n>> ")
        dbfld_path = os.path.join("cea_db", fld_name, "csv_database")
        if os.path.exists(dbfld_path):
            inst = Read_datset(dbfld_path)
            param_name = input("\nInput parameter name (That is same as csv file name.)\n>> ")
            print("\n\nPlease input the range of O/F where you want to plot.\ne.g. If the range is 0.5 to 5.0 \n0.5 5.0")
            of_range = list(map(lambda x: float(x) ,input(">> ").split()))
            print("\n\nPlease input the range of Pc where you want to plot.\nRange: 0.2 ~ 99.99, Minimum interval: 00.1 MPa\ne.g. If the range is 0.5 to 5.1 MPa and the interval is 0.1 MPa\n0.5 5.0 0.1")
            Pc_range = list(map(lambda x: float(x)*1.0e+6 ,input(">> ").split()))
            Pc_plot = np.arange(Pc_range[0], Pc_range[1], Pc_range[2])
            inst.plot(param_name, of_range, Pc_plot)
            break
        else:
            print("There is no such a dataset file/n{}".format(dbfld_path))

#func_cgr = inst.gen_func("MoleFraction/C(gr)")
#func_c2h4 = inst.gen_func("MoleFraction/C2H4")
#func_c2h6 = inst.gen_func("MoleFraction/C2H6")
#func_ch4 = inst.gen_func("MoleFraction/CH4")
#func_co = inst.gen_func("MoleFraction/CO")
#func_co2 = inst.gen_func("MoleFraction/CO2")
#func_cooh = inst.gen_func("MoleFraction/COOH")
#func_h = inst.gen_func("MoleFraction/H")
#func_h2 = inst.gen_func("MoleFraction/H2")
#func_h2o = inst.gen_func("MoleFraction/H2O")
#func_hco = inst.gen_func("MoleFraction/HCO")
#func_o = inst.gen_func("MoleFraction/O")
#func_o2 = inst.gen_func("MoleFraction/O2")
#func_oh = inst.gen_func("MoleFraction/OH")
#of_range = np.arange(0.1, 7.0, 0.1)
#Pc = 1.0e+6
#plt.xlabel("O/F")
#plt.ylabel("Mole fraction")
#plt.plot(of_range, np.array([func_cgr(of, Pc) for of in of_range]), label="C(gr)")
#plt.plot(of_range, np.array([func_c2h4(of, Pc) for of in of_range]), label="C2H4")
#plt.plot(of_range, np.array([func_c2h6(of, Pc) for of in of_range]), label="C2H6")
#plt.plot(of_range, np.array([func_ch4(of, Pc) for of in of_range]), label="CH4")
#plt.plot(of_range, np.array([func_co(of, Pc) for of in of_range]), label="CO")
#plt.plot(of_range, np.array([func_co2(of, Pc) for of in of_range]), label="CO2")
#plt.plot(of_range, np.array([func_cooh(of, Pc) for of in of_range]), label="COOH")
#plt.plot(of_range, np.array([func_h(of, Pc) for of in of_range]), label="H")
#plt.plot(of_range, np.array([func_h2(of, Pc) for of in of_range]), label="H2")
#plt.plot(of_range, np.array([func_h2o(of, Pc) for of in of_range]), label="H2O")
#plt.plot(of_range, np.array([func_hco(of, Pc) for of in of_range]), label="HCO")
#plt.plot(of_range, np.array([func_o(of, Pc) for of in of_range]), label="O")
#plt.plot(of_range, np.array([func_o2(of, Pc) for of in of_range]), label="O2")
#plt.plot(of_range, np.array([func_oh(of, Pc) for of in of_range]), label="OH")
#plt.legend(fontsize=10, loc="upper right")