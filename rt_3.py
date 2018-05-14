# -*- coding: utf-8 -*-
"""
Created on Sat May 12 21:42:28 2018

@author: TJ
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import integrate
from scipy import optimize
from tqdm import tqdm
import rt_1

class Main:
    def __init__(self, ex_df, func_cstr, func_gamma, input_param):
        self.ex_df = ex_df
        self.func_cstr = func_cstr
        self.func_gamma = func_gamma
        self.input_param = input_param
        self.of_init = rt_1.main(self.ex_df, self.input_param).of
        self.anl_df = pd.DataFrame([], index=self.ex_df.index)
        self.counter_lmbd_iterat = 0
        self.Pa = 0.1013*1.0e+6 #Atmospheric pressure [Pa]

    def do_iterat(self, maxiter=10):
        self.iterat_brentq_lmbd(maxiter)

    def iterat_brentq_lmbd(self, maxiter, lmbd_min=0.01, lmbd_max=1.0):
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
        self.counter_lmbd_iterat = 1
        self.lmbd = optimize.brentq(self.func_Mf, lmbd_min, lmbd_max, maxiter=maxiter, xtol=1.0e-3, full_output=False)
        return(self.lmbd)


    def func_Mf(self, lmbd):
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
        of = self.iterat_newton_of(lmbd)
        mf = np.array(self.ex_df.mox)/of
        self.anl_df["mf"] = mf
        self.Mf = integrate.simps(mf, np.array(self.ex_df.index))
        self.diff_Mf = self.Mf - self.input_param["Mf"]
        print("Difference of Mf = {} [kg]\n\n".format(self.diff_Mf))
        self.counter_lmbd_iterat += 1
        return(self.diff_Mf)


    def iterat_newton_of(self, lmbd):
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
        for i in tqdm(self.ex_df.index):
#==============================================================================
#             func = lambda of: lmbd*self.ex_df.mox[i]*func_Ve(of, self.ex_df.Pc[i], eps, self.func_cstr, self.func_gamma)\
#             / (self.ex_df.F[i] - (func_Pe(of, self.ex_df.Pc[i], eps, self.func_gamma)-self.Pa)*Ae\
#                - lmbd*func_Ve(of, self.ex_df.Pc[i], eps, self.func_cstr, self.func_gamma)*self.ex_df.mox[i])
#==============================================================================
            tmp = optimize.newton(self.func_of, self.of_init[i], maxiter=10, tol=1.0e-5, args=(i, lmbd, eps, Ae))
            self.of = np.append(self.of, tmp)
            j +=1
        self.anl_df["of"] = self.of
        plt.plot(self.anl_df.of, label="{} iteration".format(self.counter_lmbd_iterat))
        plt.xlabel("Time [s]")
        plt.ylabel("O/F [-]")
        plt.legend()
        plt.show()
        return(self.of)

    def func_of(self, of, time, lmbd, eps, Ae):
        Ve = func_Ve(of, self.ex_df.Pc[time], eps, self.func_cstr, self.func_gamma)
        Pe = func_Pe(of, self.ex_df.Pc[time], eps, self.func_gamma)
        mox = self.ex_df.mox[time]
        F = self.ex_df.F[time]
        of_cal = lmbd*mox*Ve / (F - (Pe-self.Pa)*Ae - lmbd*Ve*mox)
        self.diff_of = of - of_cal
        return(self.diff_of)
    
    
def func_Pe(of, Pc, eps, func_gamma):
    gam = func_gamma(of, Pc*1.0e-6)[0]
    func = lambda Pe: eps - np.power((2/(gam+1)), 1/(gam-1)) * np.power(Pc/Pe, 1/gam) / np.sqrt((gam+1)/(gam-1)*(1-np.power(Pe/Pc, (gam-1)/gam)))   
    Pe = optimize.brentq(func, 1, Pc/2, maxiter=10, xtol=1.0e+3, full_output=False)        
    return(Pe)

def func_Ve(of, Pc, eps, func_cstr, func_gamma):
    gam = func_gamma(of, Pc*1.0e-6)[0]
    SON_c = func_cstr(of, Pc*1.0e-6)[0]*gam*np.sqrt(np.power(2/(gam+1), (gam+1)/(gam-1)))
    Pe = func_Pe(of, Pc, eps, func_gamma)
    Ve = np.sqrt(2/(gam-1)*(1-np.power(Pe/Pc, (gam-1)/gam))) * SON_c
    return(Ve)
    
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
#    result.do_iterat()
    result.iterat_newton_of(0.9)
