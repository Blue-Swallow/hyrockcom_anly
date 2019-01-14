# -*- coding: utf-8 -*-
"""
Created on Sat Apr 14 13:40:16 2018

@author: T.J.-LAB-PC
"""

import os, sys
import copy
import json
import numpy as np
import pandas as pd
import cea_exe
import cea_pre
import cea_post
import rt_1
import rt_3
import rt_5


class RT():
    """
    Class executing re-construction technique and acquiering calculated results.
    
    Parameters
    ---------
    inst: class
        instance generated by "Cui_input()" class
        
    Class variable
    --------------
    self.of: 1-d ndarray
        list of oxidizer to fuel ratio calculated by CEA
    
    self.Pc: 1-d ndarray
        list of chamber pressure by CEA
        
    self.ex_df: pandas.DataFrame
        data frame of experimet data

    self.cstr: function
        function of characteristic exhaust velocity: self.cstr(of, Pc)
        
    self.gamma: function
        function of specific heat ratio at the end of chamber: self.gamma(of, Pc)
    """
    def __init__(self, inst):
        self.inst = inst
        self.of = inst.cea_db.of #range of O/F in the database calculated by CEA
        self.Pc = inst.cea_db.Pc #range of Chamber Pressure in the database calculated by CEA
        self.ex_df = inst.ex_df #Data Frame of experiment data
        self.input_param = inst.input_param
        self.cstr = inst.cea_db.gen_func("CSTAR") #data-base of characteristics exhaust velocity
        self.gamma = inst.cea_db.gen_func("GAMMAs_c") #data-base of specific heat ratio at chamber
        
    def call_rt(self):
        """ Execute RT and return the analized data as data frame
        
        Return
        -------
        anl_df: DataFrame
            analized data with using Re-construction technique
        """
        print("\nNow executing RT calculation. Please wait.")
        if self.input_param["mode"] == 1:
            anl_df = rt_1.Main(self.ex_df, self.input_param)
            if self.input_param["mode_error"] == "y":
                print("\nNow analyzing chamber pressure error...")
                delta_Pc = self.ex_df.Pc*1.1
                ex_df_Pc = self.ex_df.copy()
                ex_df_Pc.Pc = delta_Pc
                anl_df_dPc = rt_1.Main(ex_df_Pc, self.input_param)
                dof_dPc = (anl_df_dPc.of - anl_df.of)/(delta_Pc-self.ex_df.Pc)
                print("\nNow analyzing oxidizer mass flow rate error...")
                delta_mox = self.ex_df.mox*1.1
                ex_df_mox = self.ex_df.copy()
                ex_df_mox.mox = delta_mox
                anl_df_dmox = rt_1.Main(ex_df_mox, self.input_param)
                dof_dmox = (anl_df_dmox.of - anl_df.of)/(delta_mox - self.ex_df.mox)
                print("\nNow analyzing fuel mass consumption error...")
                delta_Mf = self.input_param["Mf"]*1.1
                input_param_Mf = copy.deepcopy(self.input_param)
                input_param_Mf["Mf"] = delta_Mf
                anl_df_Mf = rt_1.Main(self.ex_df, input_param_Mf)
                dof_dMf = (anl_df_Mf.of - anl_df.of)/(delta_Mf - self.input_param["Mf"])
                print("\nNow analyzing nozzle throat diameter error...")
                delta_Dt = self.input_param["Dt"]*1.1
                input_param_Dt = copy.deepcopy(self.input_param)
                input_param_Dt["Dt"] = delta_Dt
                anl_df_Dt = rt_3.Main(self.ex_df, input_param_Dt)
                dof_dDt = np.pi*(anl_df_Dt.of - anl_df.of)/(delta_Dt - self.input_param["Dt"])
                
                dPc = self.input_param["dPc"]
                dmox = self.input_param["dmox"]
                dMf = self.input_param["dMf"]
                dDt = self.input_param["dDt"]
                dof = np.sqrt(np.power(dof_dPc,2)*np.power(dPc,2) + np.power(dof_dMf,2)*np.power(dMf,2) + np.power(dof_dmox,2)*np.power(dmox,2) + np.power(dof_dDt,2)*np.power(dDt,2))
                anl_df["dof"] = dof
                of_tmp = anl_df.of + dof
                mf_tmp = anl_df.mox/of_tmp
                anl_df["dmf"] = np.abs(mf_tmp - anl_df["mf"])

#        if self.input_param["mode"] == 2:
#            anl_df = rt_2.main(self.ex_df, self.of, self.Pc, self.cstr, self.gamma, self.input_param)
        if self.input_param["mode"] == 3:
            anl_df = rt_3.Main(self.ex_df, self.cstr, self.gamma, self.input_param).execute_RT()
            if self.input_param["mode_error"] == "y":
                print("\nNow analyzing chamber pressure error...")
                delta_Pc = self.ex_df.Pc*1.1
                ex_df_Pc = self.ex_df.copy()
                ex_df_Pc.Pc = delta_Pc
                anl_df_dPc = rt_3.Main(ex_df_Pc, self.cstr, self.gamma, self.input_param).execute_RT()
                dof_dPc = (anl_df_dPc.of - anl_df.of)/(delta_Pc-self.ex_df.Pc)
                print("\nNow analyzing oxidizer mass flow rate error...")
                delta_mox = self.ex_df.mox*1.1
                ex_df_mox = self.ex_df.copy()
                ex_df_mox.mox = delta_mox
                anl_df_dmox = rt_3.Main(ex_df_mox, self.cstr, self.gamma, self.input_param).execute_RT()
                dof_dmox = (anl_df_dmox.of - anl_df.of)/(delta_mox - self.ex_df.mox)
                print("\nNow analyzing thrust error...")
                delta_F = self.ex_df.F*1.1
                ex_df_F = self.ex_df.copy()
                ex_df_F.F = delta_F
                anl_df_dF = rt_3.Main(ex_df_F, self.cstr, self.gamma, self.input_param).execute_RT()
                dof_dF = (anl_df_dF.of - anl_df.of)/(delta_F - self.ex_df.F)
                print("\nNow analyzing fuel mass consumption error...")
                delta_Mf = self.input_param["Mf"]*1.1
                input_param_Mf = copy.deepcopy(self.input_param)
                input_param_Mf["Mf"] = delta_Mf
                anl_df_Mf = rt_3.Main(self.ex_df, self.cstr, self.gamma, input_param_Mf).execute_RT()
                dof_dMf = (anl_df_Mf.of - anl_df.of)/(delta_Mf - self.input_param["Mf"])
                print("\nNow analyzing nozzle throat diameter error...")
                delta_Dt = self.input_param["Dt"]*1.1
                input_param_Dt = copy.deepcopy(self.input_param)
                input_param_Dt["Dt"] = delta_Dt
                anl_df_Dt = rt_3.Main(self.ex_df, self.cstr, self.gamma, input_param_Dt).execute_RT()
                dof_dDt = np.pi*(anl_df_Dt.of - anl_df.of)/(delta_Dt - self.input_param["Dt"])
                
                dPc = self.input_param["dPc"]
                dmox = self.input_param["dmox"]
                dF = self.input_param["dF"]
                dMf = self.input_param["dMf"]
                dDt = self.input_param["dDt"]
                dof = np.sqrt(np.power(dof_dPc,2)*np.power(dPc,2) + np.power(dof_dMf,2)*np.power(dMf,2) + np.power(dof_dmox,2)*np.power(dmox,2) + np.power(dof_dF,2)*np.power(dF,2) + np.power(dof_dDt,2)*np.power(dDt,2))
                anl_df["dof"] = dof                
                anl_df["dof"] = dof
                of_tmp = anl_df.of + dof
                mf_tmp = anl_df.mox/of_tmp
                anl_df["dmf"] = np.abs(mf_tmp - anl_df["mf"])
#        if self.input_param["mode"] == 4:
#            anl_df = rt_4.main(self.ex_df, self.cstr, self.gamma, self.input_param)
                
        if self.input_param["mode"] == 5:
            anl_df = rt_5.Main(self.ex_df, self.cstr, self.gamma, self.input_param).execute_RT()
            if self.input_param["mode_error"] == "y":
                print("\nNow analyzing chamber pressure error...")
                delta_Pc = self.ex_df.Pc*1.1
                ex_df_Pc = self.ex_df.copy()
                ex_df_Pc.Pc = delta_Pc
                anl_df_dPc = rt_5.Main(ex_df_Pc, self.cstr, self.gamma, self.input_param).execute_RT()
                dof_dPc = (anl_df_dPc.of - anl_df.of)/(delta_Pc-self.ex_df.Pc)
                print("\nNow analyzing oxidizer mass flow rate error...")
                delta_mox = self.ex_df.mox*1.1
                ex_df_mox = self.ex_df.copy()
                ex_df_mox.mox = delta_mox
                anl_df_dmox = rt_5.Main(ex_df_mox, self.cstr, self.gamma, self.input_param).execute_RT()
                dof_dmox = (anl_df_dmox.of - anl_df.of)/(delta_mox - self.ex_df.mox)
                print("\nNow analyzing fuel mass consumption error...")
                delta_Mf = self.input_param["Mf"]*1.1
                input_param_Mf = copy.deepcopy(self.input_param)
                input_param_Mf["Mf"] = delta_Mf
                anl_df_Mf = rt_5.Main(self.ex_df, self.cstr, self.gamma, input_param_Mf).execute_RT()
                dof_dMf = (anl_df_Mf.of - anl_df.of)/(delta_Mf - self.input_param["Mf"])
                print("\nNow analyzing nozzle throat diameter error...")
                delta_Dt = self.input_param["Dt"]*1.1
                input_param_Dt = copy.deepcopy(self.input_param)
                input_param_Dt["Dt"] = delta_Dt
                anl_df_Dt = rt_5.Main(self.ex_df, self.cstr, self.gamma, input_param_Dt).execute_RT()
                dof_dDt = np.pi*(anl_df_Dt.of - anl_df.of)/(delta_Dt - self.input_param["Dt"])
                
                dPc = self.input_param["dPc"]
                dmox = self.input_param["dmox"]
                dMf = self.input_param["dMf"]
                dDt = self.input_param["dDt"]
                dof = np.sqrt(np.power(dof_dPc,2)*np.power(dPc,2) + np.power(dof_dMf,2)*np.power(dMf,2) + np.power(dof_dmox,2)*np.power(dmox,2) + np.power(dof_dDt,2)*np.power(dDt,2))
                anl_df["dof"] = dof                
                anl_df["dof"] = dof
                of_tmp = anl_df.of + dof
                mf_tmp = anl_df.mox/of_tmp
                anl_df["dmf"] = np.abs(mf_tmp - anl_df["mf"])
                
        self.anl_df = anl_df
        return(self.anl_df)
        

class Cui_input():
    """
    Class to attract information through CUI to generate .inp file
    
    Class variable
    ----------
    self._langlist_: list ["jp","en",...]
        Contain initial two characters of language name.
    
    self._sntns_df_ : dict {outer_key: innerdict, ...}
        outer_key: string\n
            it is representative of the questionaire type. \n
        innerdict: string\n
            it is a questionaier sentense.

    self._lang_: string
        Selected language from the "_langlist_"
        
    self.ex_path: string
        Path containing experiment data
    
    self.ex_file: string
        The name of experimant data file (.csv)
        
    self.ex_param: dict {param_name: symbol}
        Parameter dictionary
        param_name: string; e.g., time [s], Oxidizer mass flow rate [kg/s], ...
        symbol: string; e.g., t, mox, ... 
        
    self.ex_df: DataFrame
        Data frame of measured parameter in an experiment
    
    self.input_param: dict{key, value}
        key: str, "mode"
        value: int, 1, 3, 4, 5
        the value is reperesentative number,
        1: RT-1, Re-construction technique 1; using average c*
        3: RT-3, Re-construction technique 3; using constant nozzle discharge coefficient lambda1 
        4: RT-4, Re-construction technique 4; using constant nozzle discharge coefficient lambda2 
        5: RT-5, Re-construction technique 5; using constant c* efficiency
        
        key: str, "Dt"
        value: float, nozzle throat diameter [m]
        
        key: str, "eps"
        value: float, nozzle expansion ratio [-]
        
        key: str, "Mf"
        value: float, fuel consumption [kg]
    
    self.cea_path: string
        Path containing the results of cea calculation
   
        
    """
    _langlist_ = ["jp","en"]
    _tmp_ = dict()
    _tmp_["oxid"] = {"jp": "\n\n計算オプション(0~2)を選択してください．\n例: 0: 全域平衡計算\n    1: 燃焼器内のみ平衡計算\n    2: スロートまで平衡計算",
                     "en": "\n\nPlease select option (0-2) of calculation.\ne.g.: 0: equilibrium during expansion\n      1: frozen after the end of chamber\n      2: frozen after nozzle throat"}

    def __init__(self):
        self._sntns_df_ = pd.DataFrame([], columns=self._langlist_)
        for i in self._tmp_:
            self._sntns_df_ = self._sntns_df_.append(pd.DataFrame(self._tmp_[i], index=[i]))
#        self._inp_lang_()
        self._get_expath_()
        self.input_param = dict()
        flag, cond_dict = self._select_cond_()
        if flag:
            self.input_param = cond_dict
        else:
            self._select_mode_()
            self._input_nozzle_()
            self._input_eps_()
            self._input_consump_()
            flag_error = self._select_error_analysis_()
            if flag_error:
                self._input_error_Pc_()
                self._input_error_mox_()
                self._input_error_Mf_()
                self._input_error_Dt_()
                if self.input_param["mode"] == 3:
                    self._input_error_F_()
            else:
                pass
        self._get_ceapath_()
        self.cea_db = cea_post.Read_datset(self.cea_path)
        self._cond_out_()

    def _inp_lang_(self):
        """
        Select user language
        """
        print("Please select language.\n{}".format(self._langlist_))
        lang = input()
        if lang in self._langlist_:
            self._lang_ = lang
        else:
            print("There is no such language set!")

    def _get_expath_(self):
        """
        Get the experiment folder path, file name and data and, create template file to contain an experimental data
            ex_path: string, folder path 
            ex_file: string, file name
            ex_df: data frame of experintal data
            ex_param: dictionary of experimental data
        """
        cadir = os.path.dirname(os.path.abspath(__file__))
        while(True):
            foldername = input("\nInput a Experiment Name (The Name of Folder Containing Experiment Data in \"exp_dat\" folder) \n>>")
            self.ex_path = os.path.join(cadir, "exp_dat" ,foldername)
            if os.path.exists(self.ex_path):
                self.ex_file = "ex_dat.csv"
                file_path = os.path.join(self.ex_path, self.ex_file)
                p_name = ("time [s]", "Oxidizer mass flow rate [g/s]", "Thrust [N]", "Chamber pressure [MPaG]")
                symbol = ("t", "mox", "F", "Pc")
                self.ex_param = dict(zip(p_name, symbol))
                if os.path.exists(file_path):
                    self.ex_df = pd.read_csv(file_path,header=1, index_col=0)
                    self.ex_df.mox = self.ex_df.mox * 1.0e-3 #convert [g/s] to [kg/s]
                    self.ex_df.Pc = self.ex_df.Pc * 1.0e+6 + 0.1013e+6 #convert [MPaG] to [Pa]
                    break
                else: # create template file
                    print("\nThere is no such a experiment data\n{}".format(file_path))
                    flag = input("\nDo you want to make a template file ?\n  y/n ?\n>>")
                    if flag == "y":
                        df = pd.DataFrame(self.ex_param, index=[0], columns=p_name)
                        df.to_csv(file_path, index= False)
                        print("\nComplete to generate a template file. \nPlease input the experimental data.")
                        print("Please aboid using any negative values and data that has vivration with time")
                    elif flag == "n":
                        sys.exit()
            else:
                print("There is no such a Folder\n{}".format(self.ex_path))           
    
    def _select_cond_(self):
        """
        Select whether using values contained in "cond.json" file
        
        Return
        -------
        flag: bool
            True: using the values
            False: don't use the values
        cond_dict: dict
            dictionary wihch contains the condition values
        """
        cond_dict = dict([])
        cond_path = os.path.join(self.ex_path, "cond.json")
        if os.path.exists(cond_path):
            cond_json = open(cond_path, "r")
            cond_dict = json.load(cond_json)
            while(True):
                print("\nDo you want to use following values contained in \"cond.json\"?")
                print(cond_dict)
                char = input("  y/n ?\n>>")
                if char == "y":
                    flag = True
                    break
                elif char == "n":
                    flag = False
                    break
                else:
                    pass
        else:
            flag = False
        return(flag, cond_dict)
        
    def _select_error_analysis_(self):
        """
        Select whether execute error analysys or not
        
        Return
        ------
        flag: bool
            True: execute error analisys
            False: don't execute error analisys
        """
        while(True):
            print("\nDo you want to execute error analisys?")
            char = input("  y/n ?\n>>")
            if char == "y":
                flag = True
                self.input_param["mode_error"] = "y"
                break
            elif char == "n":
                flag = False
                self.input_param["mode_error"] = "n"
                break
            else:
                pass
        return(flag)
        
           
    def _select_mode_(self):
        """
        Select a calculation mode; RT-1,2,3,4,5,...
        """
        while(True):
            print("\nSelect calculation mode.")
            mode = {1: "RT-1",
#                    2: "RT-2",
                    3: "RT-3",
#                    4: "RT-4",
                    5: "RT-5"}
            inp = int(input(" 1: RT-1; assuming c* is constant\n"+\
#                  "2: RT-2; assuming c* efficiency is constant\n"+\
                  " 3: RT-3; assuming nozzle discharge coefficient is constant; lambda1\n"+\
#                  "4: RT-4; assuming nozzle discharge coefficient is constant; lambda2\n"+\
                  " 5: RT-5; assuming c* efficiency is constant. RT-2 improved with initial O/F\n>>"))
            if inp in mode.keys():
                self.input_param["mode"] = inp
                break
            else:
                print("There is no such a mode \"{}\"\n".format(inp))
                
    def _input_nozzle_(self):
        """
        Input the nozzle diameter [mm]
        """
#        print("Please input nozzle throat diameter [mm]\n>>")
        self.input_param["Dt"] = float(input("\nInput nozzle throat diameter [mm]\n>>"))*1.0e-3

    def _input_eps_(self):
        """
        Input the nozzle expansion ratio [-]
        """
#        print("Please input nozzle expansion ratio")
        self.input_param["eps"] = float(input("\nInput nozzle expansion ratio\n>>"))
        
    def _input_consump_(self):
        """
        Input the fuel consumption [g]
        """
#        print("Please input fuel consumption [g]")
        self.input_param["Mf"] = float(input("\nInput fuel consumption [g]\n>>"))*1.0e-3
    
    def _input_error_Pc_(self):
        """
        Input the error of pressure sensor which measures chamber pressure [Pa]
        """
        self.input_param["dPc"] = float(input("\nInput error of pressure sensor wihch measures chamber pressue [MPa]\n>>"))*1.0e+6

    def _input_error_mox_(self):
        """
        Input the error of oxidizer pass flow rate [kg/s]
        """
        self.input_param["dmox"] = float(input("\nInput error of oxidizer mass flow rate [g/s]\n>>"))*1.0e-3

    def _input_error_F_(self):
        """
        Input the error of load cell which measures thrust [N]
        """
        self.input_param["dF"] = float(input("\nInput error of load cell which measures thrust [N]\n>>"))

    def _input_error_Mf_(self):
        """
        Input the error of total fuel mass consumption [kg]
        """
        self.input_param["dMf"] = float(input("\nInput error of fuel mass consumption [g]\n>>"))*1.0e-3

    def _input_error_Dt_(self):
        """
        Input the error of nozzle throat diamter [m]
        """
        self.input_param["dDt"] = float(input("\nInput error of nozzle throat diameter [mm]\n>>"))*1.0e-3

    
    def _get_ceapath_(self):
        """
        Return the folder path cantaining the results of cea calculation.
            cea.path: string, folder path containing "out" folder
        """
        cadir = os.path.dirname(os.path.abspath(__file__))
        while(True):
            foldername = input("\nInput the Folder Name Containing Results of CEA in \"cea_db\" folder \n>>")
            self.cea_path = os.path.join(cadir, "cea_db", foldername, "csv_database")
            if os.path.exists(self.cea_path):
                break
            else:
                print("\nThere is no such a dataset folder/n{}".format(self.cea_path))
                flag = input("\nDo you want to make a dataset of CEA result ?\n  y/n \n>>")
                if flag == "y":
                    print("\nInput some information to generate \".inp\" files to execute CEA.")
                    generate_class = cea_pre.Cui_input()
                    fld_path = generate_class.gen_all()
                    print("\nNow doing CEA calculation and output csv type data-base files. Please wait.")
                    execute_class = cea_exe.CEA_execute(fld_path=fld_path)
                    execute_class.all_exe()
                    break
                elif flag == "n":
                    pass
                
    def _cond_out_(self):
        """
        Output the input value as condition file: "cond.json"
        """
        output_fname = open(os.path.join(self.ex_path, "cond.json"), "w")
        json.dump(self.input_param, output_fname)


if __name__ == "__main__":
    inst = Cui_input()
    df = RT(inst).call_rt()
    df.to_csv(os.path.join(inst.ex_path, "result.csv"))
    

    
