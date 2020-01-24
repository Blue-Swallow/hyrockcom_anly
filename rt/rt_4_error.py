# -*- coding: utf-8 -*-
"""
Error Analysis for RT-4

@author: T.J.
"""

import copy
import pandas as pd
import numpy as np
from . import rt_4

def main(anl_df, ex_df, input_param, func_cstr, func_gamma):
    """Functioni of error analisys for RT3
    
    Parameters
    ----------
    anl_df : pandas.DataFrame
        Containing Rt-3 calculation result
    
    ex_df: pandas.DataFrame
        Containing experimental data
    
    input_param: dict
        Containing experimental condition and error of conditioni

    func_cstr: function
        Function of specific exhaust velociy generated by the gen_func() method in cea_post.Read_datset module

    func_gamma: function
        Function of specific heat generated by the gen_func() method in cea_post.Read_datset module

    Returns
    -------
    anl_df : pandas.DataFrame
        Containg error analisys result as well as Rt-3 result
    """
    print("\nNow analyzing chamber pressure error...")
    delta_Pc = ex_df.Pc*1.1
    ex_df_Pc = ex_df.copy()
    ex_df_Pc.Pc = delta_Pc
    anl_df_dPc = rt_4.Main(ex_df_Pc, func_cstr, func_gamma, input_param).execute_RT()
    dof_dPc = (anl_df_dPc.of - anl_df.of)/(delta_Pc-ex_df.Pc)
    deta_dPc = (anl_df_dPc.eta - anl_df.eta)/(delta_Pc-ex_df.Pc)
    dlmbd_dPc = (anl_df_dPc["lambda"] -anl_df["lambda"])/(delta_Pc-ex_df.Pc)
    print("\nNow analyzing oxidizer mass flow rate error...")
    delta_mox = ex_df.mox*1.1
    ex_df_mox = ex_df.copy()
    ex_df_mox.mox = delta_mox
    anl_df_dmox = rt_4.Main(ex_df_mox, func_cstr, func_gamma, input_param).execute_RT()
    dof_dmox = (anl_df_dmox.of - anl_df.of)/(delta_mox - ex_df.mox)
    deta_dmox = (anl_df_dmox.eta -anl_df.eta)/(delta_mox - ex_df.mox)
    dlmbd_dmox = (anl_df_dmox["lambda"] - anl_df["lambda"])/(delta_mox - ex_df.mox)
    print("\nNow analyzing thrust error...")
    delta_F = ex_df.F*1.1
    ex_df_F = ex_df.copy()
    ex_df_F.F = delta_F
    anl_df_dF = rt_4.Main(ex_df_F, func_cstr, func_gamma, input_param).execute_RT()
    dof_dF = (anl_df_dF.of - anl_df.of)/(delta_F - ex_df.F)
    deta_dF = (anl_df_dF.eta -anl_df.F)/(delta_F - ex_df.F)
    dlmbd_dF = (anl_df_dF["lambda"] - anl_df["lambda"])/(delta_F - ex_df.F)
    print("\nNow analyzing fuel mass consumption error...")
    delta_Mf = input_param["Mf"]*1.1
    input_param_Mf = copy.deepcopy(input_param)
    input_param_Mf["Mf"] = delta_Mf
    anl_df_Mf = rt_4.Main(ex_df, func_cstr, func_gamma, input_param_Mf).execute_RT()
    dof_dMf = (anl_df_Mf.of - anl_df.of)/(delta_Mf - input_param["Mf"])
    deta_dMf = (anl_df_Mf.eta - anl_df.eta)/(delta_Mf - input_param["Mf"])
    dlmbd_dMf = (anl_df_Mf["lambda"] - anl_df["lambda"])/(delta_Mf - input_param["Mf"])
    print("\nNow analyzing nozzle throat diameter error...")
    delta_Dt = input_param["Dt"]*1.1
    input_param_Dt = copy.deepcopy(input_param)
    input_param_Dt["Dt"] = delta_Dt
    anl_df_Dt = rt_4.Main(ex_df, func_cstr, func_gamma, input_param_Dt).execute_RT()
    dof_dDt = np.pi*(anl_df_Dt.of - anl_df.of)/(delta_Dt - input_param["Dt"])
    deta_dDt = (anl_df_Dt.eta - anl_df.eta)/(delta_Dt - input_param["Dt"])
    dlmbd_dDt = (anl_df_Dt["lambda"] - anl_df["lambda"])/(delta_Dt - input_param["Dt"])
    
    dPc = input_param["dPc"]
    dmox = input_param["dmox"]
    dF = input_param["dF"]
    dMf = input_param["dMf"]
    dDt = input_param["dDt"]
    dof = np.sqrt(np.power(dof_dPc,2)*np.power(dPc,2) + np.power(dof_dMf,2)*np.power(dMf,2) + np.power(dof_dmox,2)*np.power(dmox,2) + np.power(dof_dF,2)*np.power(dF,2) + np.power(dof_dDt,2)*np.power(dDt,2))
    anl_df["dof"] = dof                
    anl_df["dmf"] = np.sqrt(np.power(1/anl_df.of*dmox, 2) + np.power(-anl_df.mox/np.power(anl_df.of,2)*dof, 2))
    anl_df["deta"] = np.sqrt(np.power(deta_dPc*dPc, 2) + np.power(deta_dmox*dmox, 2) + np.power(deta_dMf*dMf, 2) + np.power(deta_dDt*dDt, 2) + np.power(deta_dF*dF, 2))
    anl_df["dlambda"] = np.sqrt(np.power(dlmbd_dPc*dPc, 2) + np.power(dlmbd_dmox*dmox, 2) + np.power(dlmbd_dMf*dMf, 2) + np.power(dlmbd_dDt*dDt, 2) + np.power(dlmbd_dF*dF, 2))
    return anl_df
