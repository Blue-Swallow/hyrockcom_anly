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
# from . import rt_5
from . import rt_1


class Main:
    def __init__(self, ex_df, func_cstr, func_gamma, input_param, plot=False):
        self.ex_df = ex_df
        self.func_cstr = func_cstr
        self.func_gamma = func_gamma
        self.input_param = input_param
        self.plot = plot
        self.func_Pe = class_Pe(self.func_gamma, self.input_param["eps"]).gen_func()
        self.of_init = rt_1.Main(self.ex_df, self.input_param).of
        # self.of_init = rt_5.Main(self.ex_df, self.func_cstr, self.func_gamma, self.input_param, ).execute_RT().of
        self.anl_df = pd.DataFrame([], index=self.ex_df.index)
        self.counter_lmbd_iterat = 0
        self.Pa = 0.1013*1.0e+6 #Atmospheric pressure [Pa]


    def execute_RT(self, maxiter=30, lmbd_init=0.9, filter_level=10.0):
        self.filter_level = filter_level
        self.filter_level = 1.0
        lmbd = self.iterat_lmbd(maxiter=maxiter, lmbd_init=lmbd_init, lmbd_min=0.1, lmbd_max=2.0)
        print("\nNow calculating the ramained parameters and output files. Please wait.")
        self.anl_df["lambda"] = np.array([lmbd for i in self.anl_df.index])
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            self.anl_df["Pe"] = np.array([self.func_Pe(self.anl_df.of[t], self.ex_df.Pc[t]) for t in self.anl_df.index])
            self.anl_df["Ve"] = np.array([func_Ve(self.anl_df.of[t], self.ex_df.Pc[t], self.input_param["eps"], self.func_cstr, self.func_gamma) for t in self.anl_df.index])
        self.anl_df["CF"] = self.anl_df.Ve*self.ex_df.mox*(1+1/self.anl_df.of)/(self.ex_df.Pc*np.power(self.input_param["Dt"], 2)*np.pi/4)
        self.anl_df["gamma"] = np.array([self.func_gamma(self.anl_df.of[t], self.ex_df.Pc[t]) for t in self.anl_df.index])
        self.anl_df["cstr_th"] = np.array([self.func_cstr(self.anl_df.of[t], self.ex_df.Pc[t]) for t in self.anl_df.index])
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
        warnings.filterwarnings("ignore")    
        self.of = np.array([])
        self.switch = np.array([])   
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
        dx = 1.0e-3 # dO/F
        x_init = 1.0 # initial x when newton iteration is conducted.
        x_min = 1.0e-2
        x_max = 1.0e+2
        eps = self.input_param["eps"]
        Ae = np.power(self.input_param["Dt"], 2)*np.pi/4 * eps
        self.of_init = self.of_init.where(self.of_init>0, other=1.0e-2)
        for i in tqdm(self.ex_df.index):
            dydx = lambda x: (self.func_right_eq9(x+dx, i, lmbd, eps, Ae) - self.func_right_eq9(x-dx, i, lmbd, eps, Ae))/(2*dx) # delivative of the left hand side of Eq.14
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
            ms_upper_end = self.func_right_eq9(xx_1, i, lmbd, eps, Ae) # multisolution region upper end
            if ms_upper_end <0: # if upper end of multi solution region less than zero, skip following caluculation.
                ms_bottom_end = self.func_right_eq9(xx_0, i, lmbd, eps, Ae) # multisolution region bottom end
                ms_left_end = self.iterat_newton_of_ms(ms_upper_end, x_min, xx_0, i, lmbd, eps, Ae)
                ms_right_end = self.iterat_newton_of_ms(ms_bottom_end, xx_1, x_max, i, lmbd, eps, Ae)
                filter_center = (ms_upper_end + ms_bottom_end)/2
                filter_width = (ms_upper_end - ms_bottom_end)*self.filter_level
                self.filter_upper_end = filter_center + filter_width/2
                self.filter_bottom_end = filter_center - filter_width/2
                self.filter_left_end = self.iterat_newton_of_ms(self.filter_upper_end, x_min, xx_0, i, lmbd, eps, Ae)
                try:
                    self.filter_right_end = self.iterat_newton_of_ms(self.filter_bottom_end, xx_1, x_max, i, lmbd, eps, Ae)
                except:
                    # ############ [ TEMP PLOT ] #################
                    print("\nWarning Warning\n")
                    tmp_x_array = np.arange(x_min, 35+dx/2, 0.1)
                    plt.plot(tmp_x_array, np.array([self.func_right_eq9(x, i, lmbd, eps, Ae) for x in tmp_x_array]), color="r")
                    plt.plot(tmp_x_array, np.array([self.filter_bottom_end for x in tmp_x_array]), color="b")
                    plt.ylim(self.filter_bottom_end-filter_width, self.filter_upper_end)
                    plt.xlabel("O/F")
                    plt.ylabel("Thrust $F$ [N]")
                    plt.savefig("tmp.png", dpi=400)
                    import sys
                    sys.exit()
                    
                    # ############ [ TEMP PLOT ] #################
            else:
                self.upper_end = 0
                self.bottom_end = 0
            of_bound_max = 1.0e+3 # maximum value of O/F boundary at O/F iteration 
            of_init = self.of_init.where(self.of_init>0, other=of_bound_max)
            self.of_init = of_init
            lhs = self.ex_df.F[i]
            try:
                if (self.filter_upper_end > lhs) and (lhs > self.filter_bottom_end): # using Patch
                    tmp = optimize.newton(self.func_error_eq9_patch,  self.of_init[i], maxiter=100, tol=1.0e-3, args=(i, lmbd, eps, Ae))
                    switch_tmp = "ON"
                else: # using RT-3
                    tmp = optimize.newton(self.func_error_eq9,  self.of_init[i], maxiter=100, tol=1.0e-3, args=(i, lmbd, eps, Ae))
                    switch_tmp = "OFF"
            except:
                try:
                    if (self.filter_upper_end > lhs) and (lhs > self.filter_bottom_end): # using Patch
                        tmp = optimize.brentq(self.func_error_eq9_patch, 1.0e-2, of_bound_max, maxiter=100, xtol=1.0e-3, args=(i, lmbd, eps, Ae))
                        switch_tmp = "ON"
                    else: # using RT-3
                        tmp = optimize.brentq(self.func_error_eq9, 1.0e-2, of_bound_max, maxiter=100, xtol=1.0e-3, args=(i, lmbd, eps, Ae))
                        switch_tmp = "OFF"
                except ValueError:
                    # if O/F exceeds of_bound_max and stop at optimization, maximum value of O/F boundary is assined as O/F
                    tmp = of_bound_max
            self.of = np.append(self.of, tmp)
            self.switch = np.append(self.switch, switch_tmp)
            j +=1
        self.anl_df["of"] = self.of
        self.anl_df["patch"] = self.switch
        if self.plot:
            plt.plot(self.anl_df.of, label="{} iteration".format(self.counter_lmbd_iterat))
            plt.xlabel("Time [s]")
            plt.ylabel("O/F [-]")
            plt.legend()
            plt.show()
        warnings.resetwarnings()
        return(self.of)

    def iterat_newton_of_ms(self, lhs, x_left, x_right, t, lmbd,eps,Ae):    
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
        of = optimize.brentq(self.func_error_eq9_ms, x_left, x_right, maxiter=100, xtol=1.0e-5, args=(lhs,t,lmbd,eps,Ae))
        return(of)

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
        """ Return the error of Eq.9
        """
        left_eq9 = self.func_left_eq9(t)
        right_eq9 = self.func_right_eq9(of,t,lmbd,eps,Ae)
        diff = right_eq9 - left_eq9
        error = diff/left_eq9
        return(error)

    def func_error_eq9_patch(self, of, t, lmbd, eps, Ae):
        """ Return the error of Eq.14
        """
        right_eq9 = self.func_right_eq9_patch(of)
        left_eq9 = self.func_left_eq9(t)
        diff = right_eq9 - left_eq9
        error = diff/left_eq9
        return(error)

    def func_error_eq9_ms(self, of, lhs, t, lmbd, eps, Ae):
        """ Return the error of Eq.14
        """
        right_eq9 = self.func_right_eq9(of, t, lmbd, eps, Ae)
        left_eq9 = lhs
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
        Right-hand side value of Eq.(9), Nagata et.al., "Evaluations of Data Reduction Methods for Hybrid Rockets", 65th IAC, 2014.

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

    def func_right_eq14_patch(self, of):
        """
        Left-hand side value of Eq.(14), Nagata et.al., "Evaluations of Data Reduction Methods for Hybrid Rockets", 65th IAC, 2014.
        
        Return
        -----------
        left_val: float
            eta *cstr_th *(1+1/(O/F)) 
        """
        inclination = (self.filter_bottom_end - self.filter_upper_end)/(self.filter_right_end - self.filter_left_end)
        section = self.filter_upper_end - inclination
        right_val = inclination *of + section
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
        warnings.filterwarnings("ignore")
        if of<=0:
            of = 1.0e-2
        try:
            Pe = optimize.newton(self.func_error_eps, Pc/2, maxiter=100, tol=1.0e+3, args=(of, Pc))
        except:
            Pe = optimize.brentq(self.func_error_eps, 1, Pc/2, maxiter=100, xtol=1.0e+3, full_output=False, args=(of, Pc))
        self.Pe=Pe
        warnings.resetwarnings()
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
#    import sys
#    sys.path.append(os.getcwd())
    import HyRockCom_Anly_cui
    inst = HyRockCom_Anly_cui.Cui_input()
    db_of = HyRockCom_Anly_cui.RT(inst).of
    db_Pc = HyRockCom_Anly_cui.RT(inst).Pc
    ex_df = HyRockCom_Anly_cui.RT(inst).ex_df
    func_cstr = HyRockCom_Anly_cui.RT(inst).cstr
    func_gamma = HyRockCom_Anly_cui.RT(inst).gamma
    input_param = HyRockCom_Anly_cui.RT(inst).input_param
    
    result = Main(ex_df, func_cstr, func_gamma, input_param)
    df = result.execute_RT()
    df.to_csv(os.path.join(inst.ex_path, "result.csv"))
