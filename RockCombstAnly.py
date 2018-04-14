# -*- coding: utf-8 -*-
"""
Created on Sat Apr 14 13:40:16 2018

@author: T.J.-LAB-PC
"""

import os, sys
import numpy as np
import pandas as pd
import execute as ex


class Cui_input():
    """
    Class to attract information through CUI to generate .inp file
    
    Parameters
    ----------
    langlist: list ["jp","en",...]
        Contain initial two characters of language name.
    
    sntns_df : dict {outer_key: innerdict, ...}
        outer_key: string\n
            it is representative of the questionaire type. \n
        innerdict: string\n
            it is a questionaier sentense.

    lang: string
        Selected language from the "langlist"
        
    ex_path: string
        Path containing experiment data
    
    ex_file: string
        experimant data file (.csv)
        
    mode: int
        Reperesentative number,
        1: RT-1, Re-construction technique 1 
        2: RT-2, Re-construction technique 2
        3: RT-3, Re-construction technique 3 
        4: RT-4, Re-construction technique 4
        5: RT-5, Re-construction technique 5
    
    cea_path: string
        Path containing the results of cea calculation
        
    
        
    """
    langlist = ["jp","en"]
    _tmp_ = dict()
    _tmp_["oxid"] = {"jp": "\n\n計算オプション(0~2)を選択してください．\n例: 0: 全域平衡計算\n    1: 燃焼器内のみ平衡計算\n    2: スロートまで平衡計算",
                     "en": "\n\nPlease select option (0-2) of calculation.\ne.g.: 0: equilibrium during expansion\n      1: frozen after the end of chamber\n      2: frozen after nozzle throat"}

    def __init__(self):
        self.sntns_df = pd.DataFrame([], columns=self.langlist)
        for i in self._tmp_:
            self.sntns_df = self.sntns_df.append(pd.DataFrame(self._tmp_[i], index=[i]))
#        self._inp_lang_()
        self._get_expath_()
        self._get_ceapath_()

    def _inp_lang_(self):
        """
        Select user language
        """
        print("Please select language.\n{}".format(self.langlist))
        lang = input()
        if lang in self.langlist:
            self.lang = lang
        else:
            print("There is no such language set!")

    def _get_expath_(self):
        """
        Return the folder path cantaining experiment data.
        """
        cadir = os.path.dirname(os.path.abspath(__file__))
        while(True):
            print("\nInput a Experiment Name (The Name of Folder Containing Experiment Data) >>")
            foldername = input()
            self.ex_path = os.path.join(cadir, foldername)
            if os.path.exists(self.ex_path):
                filename = "ex_dat.csv"
                self.ex_file = os.path.join(self.ex_path, filename)
                if os.path.exists(self.ex_file):
                    break
                else:
                    print("\nThere is no such a experiment data/n{}".format(self.ex_path))
                    print("\nDo you want to make a template file ?\ny/n ?")
                    flag = input()
                    if flag == "y":
                        dic = {"t":"[s]",
                               "mox":"[kg/s]",
                               "Pc":"[MPaG]",
                               "F":"[N]"}
                        df = pd.DataFrame(dic, index=[0])
                        df.to_csv(self.ex_file, index= False)
                    elif flag == "n":
                        sys.exit()
            else:
                print("\nThere is no such a Folder\n{}".format(self.ex_path))           


    def _get_ceapath_(self):
        """
        Return the folder path cantaining the results of cea calculation.
        """
        cadir = os.path.dirname(os.path.abspath(__file__))
        while(True):
            print("\nInput the Folder Name Containing Results of CEA >>")
            foldername = input()
            self.cea_path = os.path.join(cadir, foldername)
            self.cea_path = os.path.join(self.cea_path, "out")
            if os.path.exists(self.cea_path):
                break
            else:
                print("There is no such a dataset folder/n{}".format(self.cea_path))

# =============================================================================
#     def _inp_oxid_(self):
#         """
#         Input oxidizer species
#         """
#         print(self.sntns_df[self.lang]["oxid"])
#         oxid = input()
#         self.oxid = oxid
# =============================================================================

    

    
if __name__ == "__main__":
    test = Cui_input()

    
