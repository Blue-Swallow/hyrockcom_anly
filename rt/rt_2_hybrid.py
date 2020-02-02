# -*- coding: utf-8 -*-
"""
Created on Sat Apr 21 18:03:51 2018

@author: T.J.-LAB-PC
"""

import numpy as np
import pandas as pd
import warnings
import matplotlib.pyplot as plt
from scipy import integrate
from scipy import optimize
from . import rt_1
from tqdm import tqdm
import os

class Main:
    def __init__(self, ex_df, func_cstr, func_gamma, input_param, plot=False):
        self.ex_df = ex_df
        self.func_cstr = func_cstr
        self.func_gamma = func_gamma
        self.input_param = input_param
        self.plot = plot
        self.of_rt1 = rt_1.Main(self.ex_df, self.input_param).of
        self.anl_df = pd.DataFrame([], index=self.ex_df.index)
        self.counter_eta_iterat = 0
        
    def execute_RT(self, maxiter=30, eta_init=0.8, filter_level=10.0):
        dx = 1.0e-3 # dO/F
        x_init = 1.0 # initial x when newton iteration is conducted.
        x_min = 1.0e-2
        x_max = 1.0e+2
        dydx = lambda x: (self.func_left_eq14_ms(x+dx) - self.func_left_eq14_ms(x-dx))/(2*dx) # delivative of the left hand side of Eq.14
        """ plot part to confirm the shape of dydx """
        # x_array = np.arange(0.1, 3.0, 1e-2)
        # plt.plot(x_array, np.array([dydx(x) for x in x_array]))
        # plt.xlim(0, 3.0)
        # plt.ylim(-20, 20)
        # plt.grid()
        # plt.xlabel(r"$O/F$ [-]")
        # plt.ylabel(r"$ \frac{d}{dO/F} c^*(1+ \frac{1}{O/F})$")
        # plt.savefig("dydx.png", dpi=400)
        xx_0 = optimize.newton(dydx, x_init) # seek the solution of dydx; namely the peak of concave and convex of LHS of Eq.14
        try:
            xx_1 = optimize.brentq(dydx, x_min, xx_0 -dx)
        except ValueError:
            xx_1 = xx_0
        if round(xx_0, 5) == round(xx_1, 5):
            xx_1 = optimize.brentq(dydx, xx_0 +dx, x_max)
        if xx_0 > xx_1:
            xx_tmp = xx_0
            xx_0 = xx_1
            xx_1 = xx_tmp
        self.filter_level = filter_level
        self.ms_upper_end = self.func_left_eq14_ms(xx_1)
        self.ms_bottom_end = self.func_left_eq14_ms(xx_0)
        self.ms_left_end = self.iterat_newton_of_ms(self.ms_upper_end, x_min, xx_0)
        self.ms_right_end = self.iterat_newton_of_ms(self.ms_bottom_end, xx_1, x_max)

        try:
            self.iterat_newton_eta(maxiter, eta_init=eta_init)
        except:
            self.iterat_brentq_eta(maxiter, eta_min=0.1, eta_max=2.0)
        self.anl_df["gamma"] = np.array([self.func_gamma(self.anl_df.of[t], self.ex_df.Pc[t]) for t in self.anl_df.index])
        self.anl_df["cstr_th"] = np.array([self.func_cstr(self.anl_df.of[t], self.ex_df.Pc[t]) for t in self.anl_df.index])
        self.anl_df["cstr_ex"] = self.ex_df.Pc*np.pi/4*np.power(self.input_param["Dt"], 2)/(self.ex_df.mox + self.anl_df.mf)
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
        self.eta = optimize.newton(self.func_error_Mf, eta_init, maxiter=maxiter, tol=1.0e-3)
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
        print("Error of Mf = {} [%]\n".format((diff_Mf/self.input_param["Mf"])*1.0e+2))
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
        switch = np.array([])
        warnings.filterwarnings("error")
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
        filter_center = eta*(self.ms_upper_end + self.ms_bottom_end)/2
        filter_width = (self.ms_upper_end - self.ms_bottom_end)*self.filter_level
        filter_upper_end = filter_center + filter_width/2
        filter_bottom_end = filter_center - filter_width/2
        for i in tqdm(self.ex_df.index):
            of_bound_max = 1.0e+3 # maximum value of O/F boundary at O/F iteration 
            of_init = self.of_rt1.where(self.of_rt1>0, other=of_bound_max)
            self.of_init = of_init
            rhs = self.ex_df.Pc[i]*(np.pi*np.power(self.input_param["Dt"], 2)/4)/self.ex_df.mox[i]
            try:
                if (filter_upper_end > rhs) and (rhs > filter_bottom_end): # using RT-5
                    tmp = optimize.newton(self.func_error_eq14_rt5, of_init[i], maxiter=100, tol=1.0e-5, args=(i, eta))
                    switch_tmp = "RT-5"
                else: # using RT-2
                    tmp = optimize.newton(self.func_error_eq14_rt2, of_init[i], maxiter=100, tol=1.0e-5, args=(i, eta))
                    switch_tmp = "RT-2"
            except:
                try:
                    if (filter_upper_end > rhs) and (rhs > filter_bottom_end): # using RT-5         
                        tmp = optimize.brentq(self.func_error_eq14_rt5, 1.0e-3, of_bound_max, maxiter=100, xtol=1.0e-5, args=(i, eta))
                        switch_tmp = "RT-5"
                    else: # using RT-2
                        tmp = optimize.brentq(self.func_error_eq14_rt2, 1.0e-3, of_bound_max, maxiter=100, xtol=1.0e-5, args=(i, eta))
                        switch_tmp = "RT-2"
                except ValueError:
                    tmp = of_bound_max
            of = np.append(of, tmp)
            switch = np.append(switch, switch_tmp)
            j +=1
        self.anl_df["of"] = of
        self.anl_df["switch"] = switch
        if self.plot:
            plt.plot(self.anl_df.of, label="{} iteration".format(self.counter_eta_iterat))
            plt.xlabel("Time [s]")
            plt.ylabel("O/F [-]")
            plt.legend()
            plt.show()
        warnings.resetwarnings()
        return(of)

    
    def iterat_newton_of_ms(self, rhs, x_left, x_right):
        """
        Function to calculate O/F. This function is specified to be applid to multi-solution patch.
        
        Parameter
        ---------
        rhs: float
            R.H.S. of Eq.14
        x_left: float
            left end of the seek region of O/F
        x_right: float
            right end of the seek region of O/F
           
        Return
        ---------
        of: float
            O/F, solution of c*(1+1/O/F) - R.H.S. = 0
        """
        warnings.filterwarnings("error")
        of = optimize.brentq(self.func_error_eq14_ms, x_left, x_right, maxiter=100, xtol=1.0e-5, args=(rhs,))
        warnings.resetwarnings()
        return(of)
    
    def func_error_eq14_rt2(self, of, t, eta):
        """ Return the error of Eq.14
        """
        left_eq14 = self.func_left_eq14_rt2(of,t,eta)
        right_eq14 = self.func_right_eq14(t)
        diff = left_eq14 - right_eq14
        error = diff/right_eq14
        return(error)

    def func_error_eq14_rt5(self, of, t, eta):
        """ Return the error of Eq.14
        """
        left_eq14 = self.func_left_eq14_rt5(of,t,eta)
        right_eq14 = self.func_right_eq14(t)
        diff = left_eq14 - right_eq14
        error = diff/right_eq14
        return(error)

    def func_error_eq14_ms(self, of, rhs):
        """ Return the error of Eq.14
        """
        left_eq14 = self.func_left_eq14_ms(of)
        right_eq14 = rhs
        diff = left_eq14 - right_eq14
        error = diff/right_eq14
        return(error)

    def func_left_eq14_rt2(self, of, t, eta):
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

    def func_left_eq14_rt5(self, of, t, eta):
        """
        Left-hand side value of Eq.(14), Nagata et.al., "Evaluations of Data Reduction Methods for Hybrid Rockets", 65th IAC, 2014.
        
        Return
        -----------
        left_val: float
            eta *cstr_th *(1+1/(O/F)) 
        """
        Pc = self.ex_df.Pc[t]
        cstr = self.func_cstr(self.of_init[t], Pc)
        left_val = eta*cstr*(1 + 1/of)
        return(left_val)

    def func_left_eq14_ms(self, of):
        """
        Left-hand side value of Eq.(14), Nagata et.al., "Evaluations of Data Reduction Methods for Hybrid Rockets", 65th IAC, 2014.
        This function is specified to be applid to multi-solution patch.

        Return
        -----------
        left_val: float
            cstr_th *(1+1/(O/F)) 
        """
        Pc = 1.0e+6 #[Pa] # temporary Pc which is used to calculate delivative of c*(1+1/O/F). Pc does not affect the shape of cstr
        cstr = self.func_cstr(of, Pc)
        left_val = cstr*(1 + 1/of)
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
    import HyRockCom_Anly_cui_tmp
    inst = HyRockCom_Anly_cui_tmp.Cui_input()
    db_of = HyRockCom_Anly_cui_tmp.RT(inst).of
    db_Pc = HyRockCom_Anly_cui_tmp.RT(inst).Pc
    ex_df = HyRockCom_Anly_cui_tmp.RT(inst).ex_df
    func_cstr = HyRockCom_Anly_cui_tmp.RT(inst).cstr
    func_gamma = HyRockCom_Anly_cui_tmp.RT(inst).gamma
    input_param = HyRockCom_Anly_cui_tmp.RT(inst).input_param
    
    result = Main(ex_df, func_cstr, func_gamma, input_param)
    df = result.execute_RT(maxiter=20, filter_level=10.0)
    df.to_csv(os.path.join(inst.ex_path, "result.csv"))
    

