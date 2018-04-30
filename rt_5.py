# -*- coding: utf-8 -*-
"""
Created on Sat Apr 21 18:03:51 2018

@author: T.J.-LAB-PC
"""

import numpy as np
import pandas as pd
from scipy import integrate
from scipy import optimize
import rt_1
from tqdm import tqdm

class Main:
    def __init__(self, ex_df, func_cstr, func_gamma, input_param):
        self.ex_df = ex_df
        self.func_cstr = func_cstr
        self.func_gamma = func_gamma
        self.input_param = input_param
        self.of_init = rt_1.main(self.ex_df, self.input_param).of
        self.anl_df = pd.DataFrame([], index=self.ex_df.index)
        self.counter_eta_iterat = 0
        
    def do_iterat(self, maxiter=10):
#        self.iterat_newton_eta(maxiter, eta_init=0.7)
        self.iterat_brentq_eta(maxiter)

    def iterat_brentq_eta(self, maxiter, eta_min=0.01, eta_max=1.0):
        """
        Function to calculate O/F using Newton-Raphson method; in this case, using Secant method
        
        Parameter
        ---------
        eta: float
            c* efficiency
        
        Return
        ---------
        df: pandas.DataFrame(index, columns)
            index: float
                time
            columns: string
                "of": O/F
                "mp": fuel mass flow rate, [kg/s]
        """
        self.counter_eta_iterat = 1
        self.eta = optimize.brentq(self.func_Mf, eta_min, eta_max, maxiter=maxiter, xtol=1.0e-3, full_output=False)
        return(self.eta)

    def iterat_newton_eta(self, maxiter, eta_init=1.0):
        """
        Function to calculate O/F using Newton-Raphson method; in this case, using Secant method
        
        Parameter
        ---------
        eta: float
            c* efficiency
        
        Return
        ---------
        df: pandas.DataFrame(index, columns)
            index: float
                time
            columns: string
                "of": O/F
                "mp": fuel mass flow rate, [kg/s]
        """
        self.counter_eta_iterat = 1
        self.eta = optimize.newton(self.func_Mf, eta_init, maxiter=maxiter, tol=1.0e-2)
        return(self.eta)


    def iterat_newton_of(self, eta):
        """
        Function to calculate O/F using Newton-Raphson method; in this case, using Secant method
        
        Parameter
        ---------
        eta: float
            c* efficiency
        
        ex_df: pandas.DataFrame
            which must include the parameters: "Pc",chamber pressure; "mox", oxidizer mass flow rate.
            And also its index must indicate time [s]
            
        func_cstr: function(of, Pc)
            A function which return a interpolated value (array-like)    
    
        Return
        ---------
        df: pandas.DataFrame(index, columns)
            index: float
                time
            columns: string
                "of": O/F
                "mp": fuel mass flow rate, [kg/s]
        """
        of = np.array([])
        j=0
        if self.counter_eta_iterat == 1:
            string = "st"
        elif self.counter_eta_iterat == 2:
            string = "nd"
        elif self.counter_eta_iterat == 3:
            string = "rd"
        else:
            string = "th"
        print("{}".format(self.counter_eta_iterat)+ string + " iteration")
        print("eta = {}".format(eta))
        for i in tqdm(ex_df.index):
            func = lambda of: of - eta*self.ex_df.mox[i]*func_cstr(of, self.ex_df.Pc[i])[0]
            tmp = optimize.newton(func, self.of_init[i], maxiter=10, tol=1.0e-5)
            of = np.append(of, tmp)
            j +=1
        self.anl_df["of"] = of
        plt.plot(self.anl_df.of, label="{} interation".format(self.counter_eta_iterat))
        plt.xlabel("Time [s]")
        plt.ylabel("O/F [-]")
        plt.legend()
        plt.show()
        return(of)
    
    
    def func_Mf(self, eta):
        """
        Function to calculate O/F using Newton-Raphson method; in this case, using Secant method
        
        Parameter
        ---------
        eta: float
            c* efficiency

        Return
        ---------
        diff: float
            difference between fuel consumption obtained by experiment and iterative calculation
        """
        of = self.iterat_newton_of(eta)
        mf = np.array(self.ex_df.mox)/of
        self.anl_df["mf"] = mf
        self.Mf = integrate.simps(mf, np.array(self.ex_df.index))
        self.diff_Mf = self.Mf - self.input_param["Mf"]
        print("Difference of Mf = {} [kg]\n\n".format(self.diff_Mf))
        self.counter_eta_iterat += 1
        return(self.diff_Mf)


if __name__ == "__main__":
    import RockCombstAnly
    import matplotlib.pyplot as plt
    inst = RockCombstAnly.Cui_input()
    db_of = RockCombstAnly.RT(inst).of
    db_Pc = RockCombstAnly.RT(inst).Pc
    ex_df = RockCombstAnly.RT(inst).ex_df
    func_cstr = RockCombstAnly.RT(inst).cstr
    func_gamma = RockCombstAnly.RT(inst).gamma
    input_param = RockCombstAnly.RT(inst).input_param
    
    result = Main(ex_df, func_cstr, func_gamma, input_param)
#    plt.plot(result.index, result.of)