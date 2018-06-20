# -*- coding: utf-8 -*-
"""
Created on Sat May 12 21:42:28 2018

@author: TJ
"""
import numpy as np
import pandas as pd
import warnings
import matplotlib.pyplot as plt
from scipy import integrate
from scipy import optimize
from tqdm import tqdm
import os
import rt_1
import rt_5

class Main:
    def __init__(self, ex_df, func_cstr, func_gamma, input_param, plot=False):
        self.ex_df = ex_df
        self.func_cstr = func_cstr
        self.func_gamma = func_gamma
        self.input_param = input_param
        self.plot = plot
        self.func_Pe = class_Pe(self.func_gamma, self.input_param["eps"]).gen_func()
#        self.of_init = rt_1.Main(self.ex_df, self.input_param).of
        self.of_init = rt_5.Main(self.ex_df, self.func_cstr, self.func_gamma, self.input_param, ).execute_RT().of
        self.anl_df = pd.DataFrame([], index=self.ex_df.index)
        self.counter_lmbd_iterat = 0
        self.Pa = 0.1013*1.0e+6 #Atmospheric pressure [Pa]


    def execute_RT(self, maxiter=30, lmbd_init=0.9):
        lmbd = self.iterat_lmbd(maxiter=maxiter, lmbd_init=lmbd_init, lmbd_min=0.1, lmbd_max=2.0)
        self.anl_df["lambda"] = np.array([lmbd for i in self.anl_df.index])
        self.anl_df["Pe"] = np.array([self.func_Pe(self.anl_df.of[t], self.ex_df.Pc[t]) for t in self.anl_df.index])
        self.anl_df["Ve"] = np.array([func_Ve(self.anl_df.of[t], self.ex_df.Pc[t], self.input_param["eps"], self.func_cstr, self.func_gamma) for t in self.anl_df.index])
        self.cal_eta()
        return(pd.concat([self.ex_df, self.anl_df], axis=1))

    def iterat_lmbd(self, maxiter=50, lmbd_init=0.75, lmbd_min=0.5, lmbd_max=2.0):
        """
        Function to calculate nozzle coefficient
        
        Parameter
        ---------
        maxiter: int
            the number of maimum iteration
        
        Return
        ---------
        self.lmbd: float
            nozzle discharge coefficient
        """
        try:
            self.counter_lmbd_iterat = 1
            self.lmbd = optimize.newton(self.func_error_Mf, lmbd_init, maxiter=maxiter, tol=1.0e-3)
        except:
            self.counter_lmbd_iterat = 1
            self.lmbd = optimize.brentq(self.func_error_Mf, lmbd_min, lmbd_max, maxiter=maxiter, xtol=1.0e-3, full_output=False)
        return(self.lmbd)

    def func_error_Mf(self, lmbd):
        Mf_cal = self.func_Mf_cal(lmbd)
        Mf_ex = self.func_Mf_ex()
        diff = Mf_cal - Mf_ex
        error = diff/Mf_cal
        print("Error of Mf = {} [%]\n".format(diff/Mf_ex*1.0e+2))
        self.error_Mf = error
        self.counter_lmbd_iterat += 1
        return(error)

    def func_Mf_cal(self, lmbd):
        """
        Function to calculate difference of fuel consumption
        
        Parameter
        ---------
        lmbd: float
            nozzle discharge coefficient

        Return
        ---------
        diff: float
            difference between fuel consumption obtained by experiment and iterative calculation
        """
        of = self.iterat_of(lmbd)
        mf = np.array(self.ex_df.mox)/of
        for i in np.where(mf<0): #when mf<0, zero-value is inserted to mf
            mf[i] = 0.0
        self.anl_df["mf"] = mf
        Mf_cal = integrate.simps(mf, np.array(self.ex_df.index))
        return(Mf_cal)

    def func_Mf_ex(self):
        val = self.input_param["Mf"]
        return(val)

    def iterat_of(self, lmbd):
        """
        Function to calculate O/F using Newton-Raphson method; in this case, using Secant method
        
        Parameter
        ---------
        lmbd: float
            nozzle discharge coefficient 
    
        Return
        ---------
        df: pandas.DataFrame(index, columns)
            index: float
                time
            columns: string
                "of": O/F
                "mp": fuel mass flow rate, [kg/s]
        """
        self.of = np.array([])
        warnings.filterwarnings("error")        
        j=0
        if self.counter_lmbd_iterat == 1:
            string = "st"
        elif self.counter_lmbd_iterat == 2:
            string = "nd"
        elif self.counter_lmbd_iterat == 3:
            string = "rd"
        else:
            string = "th"
        print("{}".format(self.counter_lmbd_iterat)+ string + " iteration")
        print("lambda = {}".format(lmbd))
        eps = self.input_param["eps"]
        Ae = np.power(self.input_param["Dt"], 2)*np.pi/4 * eps
        self.of_init = self.of_init.where(self.of_init>0, other=1.0e-2)
        for i in tqdm(self.ex_df.index):
            try:
                tmp = optimize.newton(self.func_error_eq9,  self.of_init[i], maxiter=100, tol=1.0e-3, args=(i, lmbd, eps, Ae))
            except:
                tmp = optimize.brentq(self.func_error_eq9, 1.0e-2, 300, maxiter=100, xtol=1.0e-3, args=(i, lmbd, eps, Ae))
            self.of = np.append(self.of, tmp)
            j +=1
        self.anl_df["of"] = self.of
        if self.plot:
            plt.plot(self.anl_df.of, label="{} iteration".format(self.counter_lmbd_iterat))
            plt.xlabel("Time [s]")
            plt.ylabel("O/F [-]")
            plt.legend()
            plt.show()
        warnings.resetwarnings()
        return(self.of)

    def func_error_of(self, of, time, lmbd, eps, Ae):
        Ve = func_Ve(of, self.ex_df.Pc[time], eps, self.func_cstr, self.func_gamma)
        Pe = self.func_Pe(of, self.ex_df.Pc[time])
        mox = self.ex_df.mox[time]
        F = self.ex_df.F[time]
        of_cal = lmbd*mox*Ve / (F - (Pe-self.Pa)*Ae - lmbd*Ve*mox)
        error_of = (of - of_cal)/of_cal
        self.error_of = error_of
        return(self.error_of)
            
    def func_error_eq9(self, of, t, lmbd,eps,Ae):
        """ Return the error of Eq.14
        """
        left_eq9 = self.func_left_eq9(t)
        right_eq9 = self.func_right_eq9(of,t,lmbd,eps,Ae)
        diff = right_eq9 - left_eq9
        error = diff/left_eq9
        return(error)

    def func_left_eq9(self, time):
        """
        Left-hand side value of Eq.(9), Nagata et.al., "Evaluations of Data Reduction Methods for Hybrid Rockets", 65th IAC, 2014.
        
        Return
        -----------
        left_val: float
            F : thrust ;[N] 
        """
        F = self.ex_df.F[time]
        left_val = F
        return(left_val)

    def func_right_eq9(self, of, time, lmbd, eps, Ae):
        """
        Right-hand side value of Eq.(14), Nagata et.al., "Evaluations of Data Reduction Methods for Hybrid Rockets", 65th IAC, 2014.

        Return
        -----------
        right_val: float
            lambda*Ve*mox*(1+1/of) + (Pe-Pa)*Ae
        """
        Ve = func_Ve(of, self.ex_df.Pc[time], eps, self.func_cstr, self.func_gamma)
        mox = self.ex_df.mox[time]
        left_term = lmbd*Ve*mox*(1+1/of)
        Pe = self.func_Pe(of, self.ex_df.Pc[time])
        right_term = (Pe-self.Pa)*Ae
        right_val = left_term + right_term
        return(right_val)
    
    def cal_eta(self):
        eta = np.array([])
        cstr_ex = np.array([])
        for time in self.anl_df.index:
            Pc = self.ex_df.Pc[time]
            cstr_th = self.func_cstr(self.anl_df.of[time], Pc)
            At = np.pi*np.power(self.input_param["Dt"], 2.0)/4
            tmp_cstr_ex = Pc*At/(self.ex_df.mox[time]+self.anl_df.mf[time])
            tmp_eta = tmp_cstr_ex/cstr_th
            cstr_ex = np.append(cstr_ex, tmp_cstr_ex)
            eta = np.append(eta, tmp_eta)
        self.anl_df["cstr_ex"]=cstr_ex
        self.anl_df["eta"]=eta


class class_Pe():
    """Class to calculate Nozzle exit pressure, Pe [Pa].
    Parameter
    ----------
    func_gamma: function
        function generated by the gen_func() method in execute.Read_datset module
    """
    def __init__(self, func_gamma, eps):
        self.func_gamma = func_gamma
        self.eps = eps
        
    def gen_func(self):
        func = self.iterat_Pe
        return(func)

    def iterat_Pe(self, of, Pc):
        """
        of: float
            O/F
        Pc: float
            Chamber pressure [Pa]
        eps: float
            Nozzle expansion ratio
        """
        if of<=0:
            of = 1.0e-2
        try:
            Pe = optimize.newton(self.func_error_eps, Pc/2, maxiter=100, tol=1.0e+3, args=(of, Pc))
        except:
            Pe = optimize.brentq(self.func_error_eps, 1, Pc/2, maxiter=100, xtol=1.0e+3, full_output=False, args=(of, Pc))
        self.Pe=Pe
        return(self.Pe)
    
    def func_error_eps(self, Pe, of , Pc):
        eps_cal = self.func_eps_cal(Pe, of, Pc)
        diff = eps_cal - self.eps
        error = diff/eps_cal
        return(error)
    
    def func_eps_cal(self, Pe, of, Pc):
        gam = self.func_gamma(of, Pc)
        eps_cal = np.power((2/(gam+1)), 1/(gam-1)) * np.power(Pc/Pe, 1/gam) / np.sqrt((gam+1)/(gam-1)*(1-np.power(Pe/Pc, (gam-1)/gam)))
        return(eps_cal)
    


def func_Ve(of, Pc, eps, func_cstr, func_gamma):
    if of<=0:
        of = 1.0e-3
    func_Pe = class_Pe(func_gamma, eps).gen_func()
    gam = func_gamma(of, Pc)
    SON_c = func_cstr(of, Pc)*gam*np.sqrt(np.power(2/(gam+1), (gam+1)/(gam-1)))
    Pe = func_Pe(of, Pc)
    Ve = np.sqrt(2/(gam-1)*(1-np.power(Pe/Pc, (gam-1)/gam))) * SON_c
    return(Ve)

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
    df = result.execute_RT()
    df.to_csv(os.path.join(inst.ex_path, "result.csv"))
