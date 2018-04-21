# -*- coding: utf-8 -*-
"""
Created on Sat Apr 21 18:03:51 2018

@author: T.J.-LAB-PC
"""

import numpy as np
import pandas as pd
from scipy import integrate
from rt_1 import func_cstr_ave


def main(ex_df, db_of, db_Pc, func_cstr, func_gamma, input_param):
    """
    Function to calculate O/F and time-averaged characteristic exhaust velocity
    
    Parameter
    ----------
    ex_df: pandas.DataFrame
        which must include the parameters: "Pc",chamber pressure; "mox", oxidizer mass flow rate.
        And also its index must indicate time [s]
    input_param: dict
        which must include the parameters: "Mf", fuel consumption; "Dt", nozzle throat diameter
    
    Return
    -------
    anl_df: pandas.DataFrame
        which include the calcualted results of O/F and averaged c*
    """
    cstr_ave = func_cstr_ave(ex_df, input_param)
    At = np.power(input_param["Dt"], 2.0)*np.pi/4
    Pc = np.array(ex_df.Pc)
    mox = np.array(ex_df.mox)
    of = (cstr_ave*mox)/(Pc*At-cstr_ave*mox)
    dic = {"of": of,
           "cstr_ave": [cstr_ave for i in ex_df.index]}
    anl_df = pd.DataFrame(dic, index=ex_df.index)
    return(anl_df)
    

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
    
    result = main(ex_df, db_of, db_Pc, func_cstr, func_gamma, input_param)
    plt.plot(result.index, result.of)