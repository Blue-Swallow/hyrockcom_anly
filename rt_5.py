# -*- coding: utf-8 -*-
"""
Created on Sat Apr 21 18:03:51 2018

@author: T.J.-LAB-PC
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import integrate
from scipy import optimize
import rt_1
from tqdm import tqdm
import os

class Main:
    def __init__(self, ex_df, func_cstr, func_gamma, input_param, plot=False):
        self.ex_df = ex_df
        self.func_cstr = func_cstr
        self.func_gamma = func_gamma
        self.input_param = input_param
        self.plot = plot
        self.of_init = rt_1.Main(self.ex_df, self.input_param).of
        self.anl_df = pd.DataFrame([], index=self.ex_df.index)
        self.counter_eta_iterat = 0
        
    def execute_RT(self, maxiter=30):
#        self.iterat_brentq_eta(maxiter, eta_min=0.5, eta_max=2.0)
        self.iterat_newton_eta(maxiter, eta_init=0.7)
        return(pd.concat([self.ex_df, self.anl_df], axis=1))
        
    def iterat_brentq_eta(self, maxiter, eta_min=0.5, eta_max=2.0):
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
        self.eta = optimize.brentq(self.func_error_Mf, eta_min, eta_max, maxiter=maxiter, xtol=1.0e-3, full_output=False)
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
        self.eta = optimize.newton(self.func_error_Mf, eta_init, maxiter=maxiter, tol=1.0e-2)
        return(self.eta)

    def func_error_Mf(self, eta):
        """
        Function to calculate the error of fuel mass consumpiton [%]
        
        Parameter
        ---------
        eta: float
            c* efficiency

        Return
        ---------
        error: float
            error of fuel mass consumption [%]
        """
        Mf = self.func_Mf(eta)
        diff_Mf = self.input_param["Mf"] - Mf
        error = diff_Mf/Mf
        print("Error of Mf = {} [%]\n\n".format(error*1.0e+2))
        self.anl_df["eta"] = np.array([eta for i in self.anl_df.index])
        return(error)

    def func_Mf(self, eta):
        """
        Function to calculate fuel mass consumption
        
        Parameter
        ---------
        eta: float
            c* efficiency

        Return
        ---------
        Mf: float
            Fuel mass consumption [kg]
        """
        of = self.iterat_newton_of(eta)
        mf = np.array(self.ex_df.mox)/of
        self.anl_df["mf"] = mf
        self.Mf = integrate.simps(mf, np.array(self.ex_df.index))

        self.counter_eta_iterat += 1
        return(self.Mf)
    

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
        for i in tqdm(self.ex_df.index):
            of_init = self.of_init.where(self.of_init>0, other=1.0e-2)
            try:
                tmp = optimize.newton(self.func_error_eq14, of_init[i], maxiter=100, tol=1.0e-5, args=(i, eta))
            except:
                tmp = optimize.brentq(self.func_error_eq14, 1.0e-3, self.of_init.max(), maxiter=100, xtol=1.0e-5, args=(i, eta))                
            of = np.append(of, tmp)
            j +=1
        self.anl_df["of"] = of
        if self.plot:
            plt.plot(self.anl_df.of, label="{} iteration".format(self.counter_eta_iterat))
            plt.xlabel("Time [s]")
            plt.ylabel("O/F [-]")
            plt.legend()
            plt.show()
        return(of)

    
    def func_error_of(self, of, t, eta):
        """ Return the error between assigned O/F and calculated O/F
        """
        of_cal = self.func_of(of, t, eta)
        diff = of_cal - of
        error = diff/of_cal
        return(error)

    def func_of(self, of, t, eta):
        """ Return: Calculated O/F; float
        """
        Pc = self.ex_df.Pc[t]
        mox = self.ex_df.mox[t]
        At = np.pi*np.power(self.input_param["Dt"], 2)/4
        cstr = self.func_cstr(of, Pc)
        of = eta*mox*cstr/(Pc*At - mox*eta*cstr)
        return(of)
    
    def func_error_eq14(self, of, t, eta):
        """ Return the error of Eq.14
        """
        left_eq14 = self.func_left_eq14(of,t,eta)
        right_eq14 = self.func_right_eq14(t)
        diff = left_eq14 - right_eq14
        error = diff/right_eq14
        return(error)

    def func_left_eq14(self, of, t, eta):
        """
        Left-hand side value of Eq.(14), Nagata et.al., "Evaluations of Data Reduction Methods for Hybrid Rockets", 65th IAC, 2014.
        
        Return
        -----------
        left_val: float
            eta *cstr_th *(1+1/(O/F)) 
        """
        Pc = self.ex_df.Pc[t]
        cstr = self.func_cstr(of, Pc)
        left_val = eta*cstr*(1 + 1/of)
        return(left_val)

    def func_right_eq14(self, t):
        """
        Right-hand side value of Eq.(14), Nagata et.al., "Evaluations of Data Reduction Methods for Hybrid Rockets", 65th IAC, 2014.

        Return
        -----------
        right_val: float
            Pc*At/mox 
        """
        Pc = self.ex_df.Pc[t]
        mox = self.ex_df.mox[t]
        At = np.pi*np.power(self.input_param["Dt"], 2)/4
        right_val = Pc*At/mox
        return(right_val)



if __name__ == "__main__":
    import RockCombstAnly
    inst = RockCombstAnly.Cui_input()
    db_of = RockCombstAnly.RT(inst).of
    db_Pc = RockCombstAnly.RT(inst).Pc
    ex_df = RockCombstAnly.RT(inst).ex_df
    func_cstr = RockCombstAnly.RT(inst).cstr
    func_gamma = RockCombstAnly.RT(inst).gamma
    input_param = RockCombstAnly.RT(inst).input_param
    
    result = Main(ex_df, func_cstr, func_gamma, input_param)
    df = result.execute_RT(maxiter=20)
    df.to_csv(os.path.join(inst.ex_path, "result.csv"))
    

