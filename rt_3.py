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
#        tmp_list = [integrate.simps(self.ex_df.mox, self.ex_df.index)/self.input_param["Mf"] for i in ex_df.index]
#        self.of_init = pd.Series(tmp_list, index=ex_df.index)
        self.anl_df = pd.DataFrame([], index=self.ex_df.index)
        self.counter_lmbd_iterat = 0
        self.Pa = 0.1013*1.0e+6 #Atmospheric pressure [Pa]
#        self.iterat_lmbd(maxiter=20)
#        self.cal_eta()

    def execute_RT(self):
        lmbd = self.iterat_lmbd(maxiter=10, lmbd_init=0.75, lmbd_min=0.5, lmbd_max=2.0)
        self.anl_df["lambda"] = np.array([lmbd for i in self.anl_df.index])
        self.anl_df["eta"] = self.cal_eta()
        self.anl_df["Pe"] = np.array([func_Pe(self.anl_df.of[t], self.ex_df.Pc[t], self.input_param["eps"], self.func_gamma).output() for t in self.anl_df.index])
        self.anl_df["Ve"] = np.array([func_Ve(self.anl_df.of[t], self.ex_df.Pc[t], self.input_param["eps"], self.func_cstr, self.func_gamma) for t in self.anl_df.index])
        return(self.anl_df)

    def iterat_lmbd(self, maxiter=10, lmbd_init=0.75, lmbd_min=0.5, lmbd_max=2.0):
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
        try:
            self.lmbd = optimize.newton(self.func_error_Mf, lmbd_init, maxiter=maxiter, tol=1.0e-3)
        except:
            self.lmbd = optimize.brentq(self.func_error_Mf, lmbd_min, lmbd_max, maxiter=maxiter, xtol=1.0e-3, full_output=False)
        return(self.lmbd)

    def func_error_Mf(self, lmbd):
        Mf_cal = self.func_Mf_cal(lmbd)
        Mf_ex = self.func_Mf_ex()
        diff = Mf_cal - Mf_ex
        error = diff/Mf_cal
        print("Difference of Mf = {} [kg]\n\n".format(diff))
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
        of = self.iterat_newton_of(lmbd)
        mf = np.array(self.ex_df.mox)/of
        for i in np.where(mf<0): #when mf<0, zero-value is inserted to mf
            mf[i] = 0.0
        self.anl_df["mf"] = mf
        Mf_cal = integrate.simps(mf, np.array(self.ex_df.index))
        return(Mf_cal)

    def func_Mf_ex(self):
        val = self.input_param["Mf"]
        return(val)

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
        of_init = self.of_init.where(self.of_init<=0, 0.1e-3)
        for i in tqdm(self.ex_df.index):
            try:
#                tmp = optimize.newton(self.func_error_of, of_init[i], maxiter=100, tol=1.0e-3, args=(i, lmbd, eps, Ae))
                tmp = optimize.newton(self.func_error_eq9, of_init[i], maxiter=100, tol=1.0e-3, args=(i, lmbd, eps, Ae))
            except:
#                tmp = optimize.brentq(self.func_error_of, -1.0e-3, 1.0e+3, maxiter=100, xtol=1.0e-3, args=(i, lmbd, eps, Ae))
                tmp = optimize.brentq(self.func_error_eq9, 1.0e-3, 1.0e+3, maxiter=100, xtol=1.0e-3, args=(i, lmbd, eps, Ae))
            self.of = np.append(self.of, tmp)
            j +=1
        self.anl_df["of"] = self.of
        plt.plot(self.anl_df.of, label="{} iteration".format(self.counter_lmbd_iterat))
        plt.xlabel("Time [s]")
        plt.ylabel("O/F [-]")
        plt.legend()
        plt.show()
        return(self.of)

    def func_error_of(self, of, time, lmbd, eps, Ae):
        Ve = func_Ve(of, self.ex_df.Pc[time], eps, self.func_cstr, self.func_gamma)
        Pe = func_Pe(of, self.ex_df.Pc[time], eps, self.func_gamma).output()
        mox = self.ex_df.mox[time]
        F = self.ex_df.F[time]
        of_cal = lmbd*mox*Ve / (F - (Pe-self.Pa)*Ae - lmbd*Ve*mox)
        error_of = (of - of_cal)/of_cal
#        error_of = of - of_cal
        self.error_of = error_of
        return(self.error_of)
            
    def func_error_eq9(self, of, t, lmbd,eps,Ae):
        """ Return the error of Eq.14
        """
        left_eq9 = self.func_left_eq9(t)
        right_eq9 = self.func_right_eq9(of,t,lmbd,eps,Ae)
        diff = right_eq9 - left_eq9
        error = diff/left_eq9
#        error = diff/right_eq9
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
        Pe = func_Pe(of, self.ex_df.Pc[time], eps, self.func_gamma).output()
        right_term = (Pe-self.Pa)*Ae
        right_val = left_term + right_term
        return(right_val)
    
    def cal_eta(self):
        eta = np.array([])
        for time in self.anl_df.index:
            Pc = self.ex_df.Pc[time]
            cstr_th = self.func_cstr(self.anl_df.of[time], Pc)
            At = np.pi*np.power(self.input_param["Dt"], 2.0)/4
            cstr_ex = Pc*At/(self.ex_df.mox[time]+self.anl_df.mf[time])
            tmp = cstr_ex/cstr_th
            eta = np.append(eta, tmp)
        self.anl_df["eta"]=eta

class func_Pe():
    def __init__(self, of, Pc, eps, func_gamma):
        self.of = of
        self.Pc = Pc
        self.eps = eps
        self.func_gamma = func_gamma
        self.Pe = self.iterat_Pe()
    
    def output(self):
        return(self.Pe)
    
    def iterat_Pe(self):
        try:
            Pe = optimize.newton(self.func_error_eps, self.Pc/2, maxiter=100, tol=1.0e+3)
        except:
            Pe = optimize.brentq(self.func_error_eps, 1, self.Pc/2, maxiter=100, xtol=1.0e+3, full_output=False)
        return(Pe)
    
    def func_error_eps(self, Pe):
        eps_cal = self.func_eps_cal(Pe)
        eps_real = self.func_eps_real()
        diff = eps_cal - eps_real
        error = diff/eps_cal
        return(error)
    
    def func_eps_cal(self, Pe):
        gam = self.func_gamma(self.of, self.Pc)
        eps_cal = np.power((2/(gam+1)), 1/(gam-1)) * np.power(self.Pc/Pe, 1/gam) / np.sqrt((gam+1)/(gam-1)*(1-np.power(Pe/self.Pc, (gam-1)/gam)))
        return(eps_cal)
    
    def func_eps_real(self):
        return(self.eps)


def func_Ve(of, Pc, eps, func_cstr, func_gamma):
    gam = func_gamma(of, Pc)
    SON_c = func_cstr(of, Pc)*gam*np.sqrt(np.power(2/(gam+1), (gam+1)/(gam-1)))
    Pe = func_Pe(of, Pc, eps, func_gamma).output()
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
    df = result.execute_RT()
#    result.iterat_lmbd(maxiter=20)
#    result.cal_eta()
#    result.iterat_newton_of(0.9)
    
#    time = 0.0
#    lmbd = 0.75
#    eps = 2.0
#    Ae = inst.input_param["eps"]*np.power(inst.input_param["Dt"],2.0)*np.pi/4
#    of_range = np.arange(-5, 10, 0.1)
#    plt.plot(of_range, [result.func_error_of(x, time, lmbd, eps, Ae) for x in of_range])
#    plt.plot(of_range, [func_Pe(x, result.ex_df.Pc[time], 2.0, func_gamma) for x in of_range])
#    optimize.newton(result.func_error_of, result.of_init[time], maxiter=100, tol=1.0e-3, args=(time, lmbd, eps, Ae))
#    optimize.newton(result.func_error_of, 35, maxiter=100, tol=1.0e-3, args=(time, lmbd, eps, Ae))
#    optimize.brentq(result.func_error_of, -1.0e+3, 1.0e+3, maxiter=50, xtol=1.0e-3, args=(time, lmbd, eps, Ae))
