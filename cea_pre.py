# -*- coding: utf-8 -*-
"""
Generate NASA-CEA input file "*.inp"
"""

import os
import numpy as np
import pandas as pd
from tqdm import tqdm

cond_name = "cond.txt"




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
    
    option: int
        represented number
        0: equilibrium during expansion
        1: frozen after the end of chamber
        2: frozen after nozzle throat
    
    oxid: string
        must use a symbol written in NASA RP-1311-P2 app.B
        
    fuel: string
        fuel name. It isn't necessary to select from NASA RP-1311-P2 app.B
    
    o_itemp: float
        Oxidizer initial temperature [K]

    f_itemp: float
        Fuel initial temperature [K]

    f_enthlpy: float
        Fuel satndard enthalpy of formation [kJ/mol]
    
    f_elem: dict, {"A": num, ...}
        A: string, symbol of contained element in fuel
        num: int, the number of element contained in 1 mol fuel [mol]
        
    eps: float
        Nozzle area ratio, Ae/At.
    """
    
    langlist = ["jp","en"]
    _tmp_ = dict()
#    _tmp_["oxid"] = {"jp": "酸化剤の種類と質量分率(%)を入力してください.\n*記号はNASA RP-1311-P2 app.B に準拠\n\n例)\nO2(L) 80, N20 20\n",
#                     "en": "Please input oxidizer name and its mass fraction (%).\n*Concerning spiecies name, please refer \"NASA RP-1311-P2 app.B.\"\n\ne.g.\nO2(L) 80, N20 20"}
    _tmp_["option"] = {"jp": "\n\n計算オプション(0~2)を選択してください．\n例: 0: 全域平衡計算\n    1: 燃焼器内のみ平衡計算\n    2: スロートまで平衡計算",
                     "en": "\n\nPlease select option (0-2) of calculation.\ne.g.: 0: equilibrium during expansion\n      1: frozen after the end of chamber\n      2: frozen after nozzle throat"}

    _tmp_["oxid"] = {"jp": "\n\n酸化剤の種類を入力してください.\n*記号はNASA RP-1311-P2 app.B に準拠\n例: O2(L)",
                     "en": "\n\nPlease input oxidizer name.\n*Concerning spiecies name, please refer \"NASA RP-1311-P2 app.B.\"\ne.g.: O2(L)"}

    _tmp_["fuel"] = {"jp": "\n\n燃料の名前を入力してください\n例: PMMA",
                     "en": "\n\nPlease input fuel name.\"\ne.g.: PMMA"}

    _tmp_["o_wt"] = {"jp": "\n\n酸化剤の質量分率[%]を入力してください",
                     "en": "\n\nPlease input oxidizer mass fraction [%]"}
    
    _tmp_["f_wt"] = {"jp": "\n\n燃料の質量分率[%]を入力してください",
                     "en": "\n\nPlease input fuel mass fraction [%]"}
    
    _tmp_["o_itemp"] = {"jp": "\n\n酸化剤の初期温度[K]を入力してください",
                     "en": "\n\nPlease input initial oxidizer temperature [K]"}
    
    _tmp_["f_itemp"] = {"jp": "\n\n燃料の初期温度[K]を入力してください",
                     "en": "\n\nPlease input initial fuel temperature [K]"}

    _tmp_["o_enthlpy"] = {"jp": "\n\n酸化剤の標準生成エンタルピ[kJ/mol]を入力してください",
                     "en": "\n\nPlease input standard enthalpy of formation [kJ/mol] respect to oxidizer"}
    
    _tmp_["f_enthlpy"] = {"jp": "\n\n燃料の標準生成エンタルピ[kJ/mol]を入力してください",
                     "en": "\n\nPlease input standard enthalpy of formation [kJ/mol] respect to fuel"}

    _tmp_["o_elem"] = {"jp": "\n\n1molの酸化剤に含まれる元素とそのmol数を入力してください.\n例: N 2 O 1",
                     "en": "\n\nPlease input the element symbols and its mol number contained per 1 mol of oxidizer\ne.g.: N 2 O 1"}
    
    _tmp_["f_elem"] = {"jp": "\n\n1molの燃料に含まれる元素とそのmol数を入力してください.\n例: C 5 H 2 O 6",
                     "en": "\n\nPlease input the element symbols and its mol number contained per 1 mol of fuel\ne.g.: C 5 H 2 O 6"}

    _tmp_["eps"] = {"jp": "\n\n開口比Ae/Atを入力してください.",
                    "en": "\n\nPlease input the area ratio, Ae/At."}
    
    _tmp_["of"] = {"jp": "\n\n計算するO/Fの範囲を入力してください.\n許容範囲：0.01~99.99 , 最小刻み幅：0.01\n例) 0.5~10 を 0.1毎に計算する場合.\n0.5　10.1　0.1",
                "en": "\n\nPlease input the range of O/F where you want to calculate.\nRange: 0.01 ~ 99.99, Minimum interval: 00.1\ne.g. If the range is 0.5 to 10 and the interval is 0.1\n0.5 10.1 0.1"}

    _tmp_["Pc"] = {"jp": "\n\n計算する燃焼室圧力[MPa]の範囲を入力してください.\n許容範囲：0.2~100 MPa, 最小刻み幅：0.01　MPa\n例) 0.5~5.0 MPa を 0.1 MPa毎に計算する場合.\n0.5　5.1　0.1",
                "en": "\n\nPlease input the range of Chamber pressure [MPa] where you want to calculate.\nRange: 0.2~100 MPa, Minimum interval: 0.01 MPa\ne.g. If the range is 0.5 to 5.1 MPa and the interval is 0.1 MPa\n0.5 5.0 0.1"}

    _tmp_["cont_react"] = {"jp": "\n化学種の入力を続けますか? \"y/n\"",
                "en": "\nDo you want to continue inputting reactant information? \"y/n\""}

    _tmp_["cont_enthlpy"] = {"jp": "\n標準生成エンタルピと構成元素を手動入力しますか? \"y/n\"",
                "en": "\nDo you want to manually input standard enthalpy and constitution of element? \"y/n\""}

    _tmp_["confirm"] = {"jp": "\n入力した内容は正確ですか? \"y/n\"",
                "en": "\nIs the inputted data correct? \"y/n\""}


    def __init__(self):
        self.sntns_df = pd.DataFrame([], columns=self.langlist)
        for i in self._tmp_:
            self.sntns_df = self.sntns_df.append(pd.DataFrame(self._tmp_[i], index=[i]))
        self._inp_lang_()
        self._inp_option_()
        self.list_oxid = self._inp_react_("oxid")
        self.list_fuel = self._inp_react_("fuel")
#        self._inp_oxid_()
#        self._inp_fuel_()
#        self._inp_o_itemp_()
#        self._inp_f_itemp_()
#        self._inp_f_enthlpy_()
#        self._inp_f_elem_()
        self._inp_eps_()
        self._inp_of_()
        self._inp_Pc_()

    def _inp_lang_(self):
        """
        Select user language
        """
        print("Please select language.\n{}".format(self.langlist))
        lang = input(">> ")
        if lang in self.langlist:
            self.lang = lang
        else:
            print("There is no such language set!")

    def _inp_option_(self):
        """
        Select a option of calculation: whether equilibrium or frozen composition.
        """
        print(self.sntns_df[self.lang]["option"])
        option = int(input(">> "))
        if option == 0:
            option = "equilibrium"
        elif option == 1:
            option = "frozen nfz=1" # frozen composition after the end of chamber
        elif option == 2:
            option = "frozen nfz=2" # frozen composition after nozzle throat
        else:
            print("Please re-confirm input integer!")
        self.option = option

    def _inp_react_(self, ident):
        while(True):
            list_react = []
            while(True):
                name = self._inp_name_(ident)
                wt = self._inp_wt_(ident)
                temp = self._inp_temp_(ident)
                print(self.sntns_df[self.lang]["cont_enthlpy"])
                flag = input(">> ")
                if flag == "y":
                    enthalpy = self._inp_enthlpy_(ident)
                    elem = self._inp_elem_(ident)
                else:
                    enthalpy = ""
                    elem = ""
                list_react.append({"name": name,
                                  "wt": wt,
                                  "temp": temp,
                                  "h": enthalpy,
                                  "elem": elem})
                print(self.sntns_df[self.lang]["cont_react"])
                flag = input(">> ")
                if flag == "n":
                    break
            print(self.sntns_df[self.lang]["confirm"])
            print(list_react)
            flag = input(">> ")
            if flag == "y":
                break
        return(list_react)

# =============================================================================
#     def _inp_oxid_(self):
#         """
#         Input oxidizer species
#         """
#         print(self.sntns_df[self.lang]["oxid"])
#         oxid = input()
# #        if oxid in self.langlist:
# #            return(lang)
# #        else:
# #            print("There is no such species in themo.inp!")
#         self.oxid = oxid
# =============================================================================

    def _inp_name_(self, ident):
        """
        Input fuel species
        """
        if ident == "fuel":
            print(self.sntns_df[self.lang]["fuel"])
        elif ident == "oxid":
            print(self.sntns_df[self.lang]["oxid"])
#        fuel = input()
        name = input(">> ")
#        if fuel in self.langlist:
#            return(lang)
#        else:
#            print("There is no such species in themo.inp!")
#        self.fuel = fuel
        return(name)

# =============================================================================
#     def _inp_o_itemp_(self):
#         """
#         Input oxidizer initial temperature
#         """
#         print(self.sntns_df[self.lang]["o_itemp"])
#         o_itemp = float(input())
# #        if fuel in self.langlist:
# #            return(lang)
# #        else:
# #            print("There is no such species in themo.inp!")
#         self.o_itemp = o_itemp
# =============================================================================

    def _inp_wt_(self, ident):
        """
        Input weight fraction
        """
#        print(self.sntns_df[self.lang]["wt"])
        if ident == "fuel":
            print(self.sntns_df[self.lang]["f_wt"])
        elif ident == "oxid":
            print(self.sntns_df[self.lang]["o_wt"])
#        print("Please input weight fractioin")
#        f_itemp = float(input())
        wt = float(input(">> "))
#        if fuel in self.langlist:
#            return(lang)
#        else:
#            print("There is no such species in themo.inp!")
#        self.f_itemp = f_itemp
        return(wt)

    def _inp_temp_(self, ident):
        """
        Input fuel initial temperature
        """
        if ident == "fuel":
            print(self.sntns_df[self.lang]["f_itemp"])
        elif ident == "oxid":
            print(self.sntns_df[self.lang]["o_itemp"])
#        f_itemp = float(input())
        temp = float(input(">> "))
#        if fuel in self.langlist:
#            return(lang)
#        else:
#            print("There is no such species in themo.inp!")
#        self.f_itemp = f_itemp
        return(temp)
    
    def _inp_enthlpy_(self, ident):
        """
        Input fuel standard enthalpy of formation
        """
        if ident == "fuel":
            print(self.sntns_df[self.lang]["f_enthlpy"])
        elif ident == "oxid":
            print(self.sntns_df[self.lang]["o_enthlpy"])
#        f_enthlpy = float(input())
        enthalpy = float(input(">> "))
#        if fuel in self.langlist:
#            return(lang)
#        else:
#            print("There is no such species in themo.inp!")
#        self.f_enthlpy = f_enthlpy
        return(enthalpy)
        
    def _inp_elem_(self, ident):
        """
        Input element contained in fuel and its mol number contained in per 1 mol
        """
        if ident == "fuel":
            print(self.sntns_df[self.lang]["f_elem"])
        elif ident == "oxid":
            print(self.sntns_df[self.lang]["o_elem"])
#        self.f_elem = {}
        elem = input(">> ")
#        elem = input().split()
#        for i in range(len(tmp)):
#            if i%2==1:
#                self.f_elem[tmp[i-1]] = float(tmp[i])
        return(elem)
                
    def _inp_eps_(self):
        """
        Input nozzle-area ratio
        """      
        print(self.sntns_df[self.lang]["eps"])
        self.eps = float(input(">> "))
#        if fuel in self.langlist:
#            return(lang)
#        else:
#            print("There is no such species in themo.inp!")
        
    def _inp_of_(self):
        """
        Input calculation range of O/F
        """      
        print(self.sntns_df[self.lang]["of"])
        self.of = list(map(lambda x: float(x) ,input(">> ").split()))
        

    def _inp_Pc_(self):
        """
        Input calculation range of chamber pressure, Pc.
        """      
        print(self.sntns_df[self.lang]["Pc"])
        self.Pc = list(map(lambda x: float(x) ,input(">> ").split()))

    def _getpath_(self):
        """
        Return the folder path which will cantain cea files: .inp, .out and csv cea-database
        """
        cadir = os.path.dirname(os.path.abspath(__file__))
        foldername = input("Input a Case Name (Folder Name)\n>>")
        path = os.path.join(cadir, "cea_db", foldername)
        return(path)


    def gen_all(self):
        """
        Generate input file with respect to every condition
        
        Parameters
        ----------
        path: string
            Folder path where this function will make "inp" floder storing ".inp"
    
        func: function    
        -----------------
        (oxid, fuel, dh, Pc, of, n, elem) = func(path)
         
        Parameters
            path: string\n
        Returns
            oxid: string, e.g."oxid=O2(L) wt=100 t,k=90.15"  \n
            fuel: string, e.g."fuel=PMMA wt=100 t,k=298.15"  \n
            dh: float  \n
            Pc: list, [start, end, interval], each element type is float  \n
            of: list, [start, end, interval], each element type is float  \n
            
        """
        path = os.path.join(self._getpath_(), "inp")
#        oxid,fuel,dh,Pc,of,n,elem  = read_cond(path)
        of = np.arange(self.of[0], self.of[1], self.of[2])
        Pc = np.arange(self.Pc[0], self.Pc[1], self.Pc[2])
#        oxid = "oxid={} wt={} t,k={}".format(self.oxid, 100, self.o_itemp)
#        fuel = "fuel={} wt={} t,k={}".format(self.fuel, 100, self.f_itemp)
    
#        elem_input = ""
#        for i in self.f_elem:
#            elem_input = elem_input+"{} {} ".format(i, self.f_elem[i])
        for i in tqdm(range(np.size(Pc))):
            for j in range(np.size(of)):
#                make_inp(path, self.option, of[j], Pc[i], oxid, fuel, self.f_enthlpy, elem_input, self.eps)
                make_inp(path, self.option, of[j], Pc[i], self.list_oxid, self.list_fuel, self.eps)
        return(path)



def make_inp(path, option, of, Pc, list_oxid, list_fuel, eps, n=""):
    """
    Write information in input file, "*.inp".
    
    Parameter
    ---------
    path : string
        Path at which "*.inp" file is saved
    option: string
        Calculation option, wtheher using equilibrium composition or frozen composition.
    of: float,
        O/F
    Pc: float,
        Camberpressure, [MPa]
    list_oxid: list,
        The list has an information about oxidizer as dict type; dict{"name":name, "wt":weight fraction, "temp":initial temperature, "h":enthalpy, "elem"element}
    list_fuel: list
        The list has an information about fuel as dict type; dict{"name":name, "wt":weight fraction, "temp":initial temperature, "h":enthalpy, "elem"element}
    eps: float,
        Area ratio of nozzle throat & exit, Ae/At
    n: int
        polimerization number
    """

#    if len(n) == 0:
#        fld_name = "inp"
#    else:
#        fld_name = os.path.join("inp", "n={}".fromat(n))
        
#    if os.path.exists(os.path.join(path, fld_name)):
#        pass
#    else:
#        print("{},  {}".format(path, fld_name))
#        os.makedirs(os.path.join(path, fld_name))
    if os.path.exists(path):
        pass
    else:
#        print("{}".format(path))
        os.makedirs(path)
    num_round = int(2) #the number of decimal places in "Pc" & "of"
#    inp_fname = "Pc_{:0^5}__of_{:0^5}.inp".format(round(Pc,num_round), round(of,num_round)) #.inp file name, e.g. "Pc=1.00_of=6.00.inp"
    inp_fname = "Pc_{:0>5.2f}__of_{:0>5.2f}.inp".format(round(Pc,num_round), round(of,num_round)) #.inp file name, e.g. "Pc=1.00_of=6.00.inp"
    file = open(os.path.join(path,inp_fname), "w")

    Pc = Pc * 10    #Pc:Chamber pressure [bar]
    prob = "case={} o/f={} rocket {} tcest,k=3800 p,bar={} sup,ae/at={}".format(inp_fname, round(of,4), option, round(Pc,4), round(eps,4))
    oxid = ""
    for i in range(len(list_oxid)):
        if len(str(list_oxid[i]["h"]))==0:
            oxid = oxid + "\toxid={} wt={} t,k={} \n".format(list_oxid[i]["name"], list_oxid[i]["wt"], list_oxid[i]["temp"])
        else:
            oxid = oxid + "\toxid={} wt={} t,k={} h,kj/mol={} {} \n".format(list_oxid[i]["name"], list_oxid[i]["wt"], list_oxid[i]["temp"], list_oxid[i]["h"], list_oxid[i]["elem"])
    fuel = ""
    for i in range(len(list_fuel)):
        if len(str(list_fuel[i]["h"]))==0:
            fuel = fuel + "\tfuel={} wt={} t,k={} \n".format(list_fuel[i]["name"], list_fuel[i]["wt"], list_fuel[i]["temp"])
        else:
            fuel = fuel + "\tfuel={} wt={} t,k={} h,kj/mol={} {} \n".format(list_fuel[i]["name"], list_fuel[i]["wt"], list_fuel[i]["temp"], list_fuel[i]["h"], list_fuel[i]["elem"])
#    outp = "siunits short"
    outp = "transport"
#    file.write("prob\n\t{}\nreact\n\t{}\n\t{}\noutput\t{}\nend\n".format(prob,oxid,fuel,outp))
    file.write("prob\n\t{}\nreact\n{}{}output\t{}\nend\n".format(prob,oxid,fuel,outp))
    file.close()



if __name__ == "__main__":
#    path = _getpath_()
#    test = gen_all(path)

    myclass = Cui_input()
#    path = "D:\T.J\Github\HybridRocketCombustionSim\Develop\RockCombustSim"
#    option = "frozen nfz=2"
#    of = 1.0
#    Pc = 1.0
#    list_oxid = myclass.list_oxid
#    list_fuel = myclass.list_fuel
#    eps = 1.0
#    make_inp(path, option, of, Pc, list_oxid, list_fuel, eps)
    myclass.gen_all()
