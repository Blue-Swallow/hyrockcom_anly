# -*- coding: utf-8 -*-
"""
Execute CEA calculation
"""

import numpy as np
from scipy import interpolate
import pandas as pd
import matplotlib.pyplot as plt
import os, sys, glob, shutil, platform
import re, copy
import tqdm
from subprocess import*
import warnings
import time
from cea_pre import make_inp_name, make_inp


class CEA_execute:
    """
    Class to excecute CEA calculation
    """

    def __init__(self, fld_path=None):
        """"        
        Parameters
        ----------
        fld_path : string, optional
            Folder's path of data-base, by default None
        """
        self.fld_path = fld_path
        self.platform = platform.system()
        pass
    
    def _getpath_(self):
        """
        Return the folders path
        
        Return
        ------
        cadir: string
            Folder path containing this code file
            Correspond with the path containing "FACE2.exe"

        inpfld_path: string,
            Folder's path containing input files, "*.inp" 

        outfld_path: string,
            Folder's path to contain output files, "*.out"
            
        dbfld_path: string
            Folder's path to contain csv database file, "*.csv"
        """
        
        cadir = os.path.dirname(os.path.abspath(__file__))
        if self.fld_path is None:
            input_path = input("Input Folder Name (e.g. \"O2+PMMA\")\n>>")
            self.fld_path = os.path.join(cadir, "cea_db", input_path)
#        print("Input Polymerization Number. If you did't assign it, please input \"n\" \n(Directory structure is like \"O2+PMMA/n=100\")")
#        num = input()
        num = "n"
        if num=="n":
#            global outfld_path
            inpfld_path = os.path.join(self.fld_path, "inp")
            outfld_path = os.path.join(self.fld_path, "out")
            dbfld_path = os.path.join(self.fld_path, "csv_database")
        else:
            inpfld_path = os.path.join(self.fld_path, "inp_n={}".format(num))
            outfld_path = os.path.join(self.fld_path, "out_n={}".format(num))
            dbfld_path = os.path.join(self.fld_path + "csv_database_n={}".format(num))
        if os.path.exists(inpfld_path):
            if os.path.exists(outfld_path):
                pass
            else:
                os.mkdir(outfld_path) #make output folder
            if os.path.exists(dbfld_path):
                pass
            else:
                os.mkdir(dbfld_path) #make output folder
            # if os.path.exists(os.path.join(dbfld_path, "MoleFraction")):
            #     pass
            # else:
            #     os.mkdir(os.path.join(dbfld_path, "MoleFraction")) #make output folder of mole_fraction
        else:
            sys.exit("There is no such a directory, \n\"{}\"".format(self.inpfld_path))
        return(cadir, inpfld_path, outfld_path, dbfld_path)


    def _csv_out_(self, dbfld_path, of, Pc, val_dict, point):
        """
        Write out calculattion results in csv-files
        
        Parameters
        ----------
        dbfld_path: string
            folder path where csv-files is outputted 
        of: list
            contains O/F values
        Pc: list,
            contains Pc values
        val_dict: dict, {key: value}
            key: string, parameter name \n
            value: 2-ndarray, parameter of row is of, parameter of column is Pc
        point: string
            identifer of file name : e.g. if point="c" -> file name = "xxx_c"
        """
        if len(point)==0:
            for i in val_dict:
                dir = os.path.dirname(os.path.join(dbfld_path, i + ".csv"))
                if not os.path.exists(dir):
                    os.mkdir(dir)
                df = pd.DataFrame(val_dict[i], index=of, columns=Pc).sort_index().sort_index(axis=1)
                df.to_csv(os.path.join(dbfld_path, i) + ".csv")
        else:
            for i in val_dict:
                df = pd.DataFrame(val_dict[i], index=of, columns=Pc).sort_index().sort_index(axis=1)
                df.to_csv(os.path.join(dbfld_path, i)+"_"+point+".csv")

    def single_exe(self, cea_dirpath, inp_fname):
        """
        One-time CEA execution
        
        Parameter
        --------
        cea_dirpath: string
            CEA directory path
        inp_fname : string
            Input file name with out ".inp" extension \n
            It is required to put ".inp" file in the same directory with "FCEA2.exe"
        """
        #cea_fname : Name and case of CEA input-file & output-file
        os.chdir(cea_dirpath)
        if self.platform == "Windows":
            cea_path = os.path.join(cea_dirpath, "FCEA2.exe")
        elif self.platform == "Linux":
            cea_path = os.path.join(cea_dirpath, "FCEA2_linux")
        elif self.platform == "Darwin":
            cea_path = os.path.join(cea_dirpath, "FCEA2_mac")
        else:
            print("Sorry, this program does not support \"{}\" system.".format(self.platform))
            sys.exit()
        command = os.path.join(cea_dirpath,inp_fname) + "\n"
        p = Popen(cea_path, stdin=PIPE, stdout=PIPE)
        p.communicate(input=bytes(command,"utf-8"))
#        time.sleep(0.1)
        p.wait()
        os.chdir("..")
        return

      
    def all_exe(self):
        """
        Execute CEA calculation with respect to all conditon: Pc & o/f
        
        Return
        ------
        cea_out: tuple, (therm_param, rocket_param)
            therm_param: dict, Thermodynamics parameters  \n
            rocket_param: dict, Rocket parameters
        """
        cadir, inpfld_path, outfld_path, dbfld_path = self._getpath_()
        cea_dirpath = os.path.join(cadir, "cea")
        split =  lambda r: os.path.splitext(r)[0] # get file name without extention
        inp_list = [os.path.basename(split(r))  for r in glob.glob(inpfld_path + "/*.inp")]        

        num_point = 0               
        for i, fname in enumerate(tqdm.tqdm(inp_list)):
            shutil.copy(os.path.join(inpfld_path,fname+".inp"), os.path.join(cadir,"cea","tmp.inp"))
            self.single_exe(cea_dirpath, "tmp")
            shutil.copy(os.path.join(cea_dirpath, "tmp.out"), os.path.join(outfld_path, fname+".out"))
            cond, therm, trans, rock, mole = Read_output("cea").read_out("tmp")

            if cond["O/F"] == 0.0 or type(cond["O/F"]) is not float:
                cond["O/F"] = float(re.sub("Pc_.....__of_" ,"", fname))
            if cond["Pc"] == 0.0 or type(cond["Pc"]) is not float:
                cond["Pc"] = float(re.sub("Pc_", "", re.sub("__of_....." ,"", fname)))
            if len(mole) != 0:
                if len(mole[list(mole.keys())[0]]) > num_point:
                    num_point =  len(mole[list(mole.keys())[0]])
            therm.update(trans) #combine dict "therm" and dict "trans"

            #initialize each container array
            if i ==0: 
                of = []
                Pc = [] 
                value_c = copy.deepcopy(therm)
                value_t = copy.deepcopy(therm)
                value_e = copy.deepcopy(therm)
                for j in therm:
                    value_c[j] = np.zeros((0,0), float)
                    value_t[j] = np.zeros((0,0), float)
                    value_e[j] = np.zeros((0,0), float)                    
                value_rock = copy.deepcopy(rock)
                for j in rock:
                    value_rock[j] = np.zeros((0,0), float)
                value_mole = copy.deepcopy(mole)
                keys_mole = list(mole.keys())
                for j in mole:
                    for k in range(mole[j].__len__()):
                        value_mole[j][k] = np.zeros((0,0), float)

            list_only = [x for x in keys_mole if x not in list(mole.keys())]

            #extend row of array when o/f is renewed
            if cond["O/F"] not in of and type(cond["O/F"]) is float:
                of.append(cond["O/F"])
                for j in therm:
                    value_c[j] = np.append(value_c[j], np.zeros((1,value_c[j].shape[1]), float), axis=0)
                    value_t[j] = np.append(value_t[j], np.zeros((1,value_t[j].shape[1]), float), axis=0)
                    value_e[j] = np.append(value_e[j], np.zeros((1,value_e[j].shape[1]), float), axis=0)
                for j in rock:
                    value_rock[j] = np.append(value_rock[j], np.zeros((1,value_rock[j].shape[1]), float), axis=0)
                for j in mole:
                    if j not in keys_mole:
                        value_mole[j] = [np.array([]) for k in range(mole[j].__len__())]
                        for k in range(mole[j].__len__()):
                            value_mole[j][k] = np.zeros((len(of), len(Pc)), float)
                    else:
                        for k in range(mole[j].__len__()):
                            value_mole[j][k] = np.append(value_mole[j][k], np.zeros((1,value_mole[j][k].shape[1]), float), axis=0)
                for i in list_only:
                    for k in range(mole[j].__len__()):
                        value_mole[i][k] = np.append(value_mole[i][k], np.zeros((1,value_mole[i][k].shape[1]), float), axis=0)
            #extend column of array when Pc is renewed
            if cond["Pc"] not in Pc and type(cond["Pc"]) is float:
                Pc.append(cond["Pc"])
                for j in therm:
                    value_c[j] = np.append(value_c[j], np.zeros((value_c[j].shape[0],1), float), axis=1)
                    value_t[j] = np.append(value_t[j], np.zeros((value_t[j].shape[0],1), float), axis=1)
                    value_e[j] = np.append(value_e[j], np.zeros((value_e[j].shape[0],1), float), axis=1)
                for j in rock:
                    value_rock[j] = np.append(value_rock[j], np.zeros((value_rock[j].shape[0],1), float), axis=1)
                for j in mole:
                    if j not in keys_mole:
                        value_mole[j] = [np.array([]) for k in range(mole[j].__len__())]
                        for k in range(mole[j].__len__()):
                            value_mole[j][k] = np.zeros((len(of), len(Pc)), float)
                    else:
                        for k in range(mole[j].__len__()):
                            value_mole[j][k] = np.append(value_mole[j][k], np.zeros((value_mole[j][k].shape[0],1), float), axis=1)
                for i in list_only:
                    # for k in range(mole[j].__len__()):
                    for k in range(num_point):
                        value_mole[i][k] = np.append(value_mole[i][k], np.zeros((value_mole[i][k].shape[0],1), float), axis=1)

            p = of.index(cond["O/F"])
            q = Pc.index(cond["Pc"])
            for j in therm:
                #Substitute each themodynamic value
                if len(therm[j]) == 0:
                    value_c[j][p,q] = np.nan
                    value_t[j][p,q] = np.nan
                    value_e[j][p,q] = np.nan
                else:
                    value_c[j][p,q] = therm[j][0]
                    value_t[j][p,q] = therm[j][1]
                    value_e[j][p,q] = therm[j][2]
            for j in rock:
                #Substitute each rocket-parameter value
                if len(rock[j]) == 0:
                    value_rock[j][p,q] = np.nan
                else:    
                    value_rock[j][p,q] = rock[j][1]
            for j in mole:
                #Substitute each mole fraction
                if j not in keys_mole:
                    value_mole[j] = [np.array([]) for k in range(mole[j].__len__())]
                    for k in range(mole[j].__len__()):
                        value_mole[j][k] = np.zeros((len(of), len(Pc)), float)
                        value_mole[j][k][p,q] = mole[j][k]
                    keys_mole.append(j)
                else:
                    for k in range(mole[j].__len__()):
                        value_mole[j][k][p,q] = mole[j][k]
                        # if there is not molecular in the dict of mole, following operation input 0.0 in to database
                        for i in list_only:
                            value_mole[i][k][p,q] = 0.0
        
        # exchange the content of "value_mole"
        tmp_mole = [np.nan for i in range(value_mole[j].__len__())]
        for i in range(value_mole[j].__len__()):
            tmp_dic = {}
            for j in value_mole:
                tmp_dic[j] = value_mole[j][i]
            tmp_mole[i] = tmp_dic
        value_mole = tmp_mole

        self._csv_out_(dbfld_path, of, Pc, value_c, point="c") #write out in csv-file
        self._csv_out_(dbfld_path, of, Pc, value_t, point="t") #write out in csv-file
        self._csv_out_(dbfld_path, of, Pc, value_e, point="e") #write out in csv-file
        self._csv_out_(dbfld_path, of, Pc, value_rock, point="") #write out in csv-file
        point_list_mole = ["Chamber", "Throat", "Exit"]
        for i in range(value_mole.__len__()):
            self._csv_out_(os.path.join(dbfld_path, "MoleFraction@" + point_list_mole[i]), of, Pc, value_mole[i], point="") #write mole-fraction out in csv-file
            
        return(of, Pc, value_c, value_t, value_e, value_rock, value_mole)



class CEA_onetime_execute(CEA_execute):
    """
    Class for executing CEA in onetime
    base: CEA_execute
    """
    def __init__(self, fld_path=None):
        super().__init__(fld_path=fld_path)
        cadir = os.path.dirname(os.path.abspath(__file__))
        self.cea_dirpath = os.path.join(cadir, "cea")

    def onetime_exe_name(self, option, list_species, Pc, eps):
        """
        Execute CEA in onetime; preparing ".inp", execute CEA and read ".out" file,
        when you assign several species and its fraction of amount.

        Parameter
        ---------
        option: string
            Calculation option, wtheher using equilibrium composition or frozen composition.
            "equilibrium", "frozen nfz=1" or "frozen nfz=2"
        list_species: list of dictionary,
            The list has an information about chemical species as dict type; dict{"name":name, "wt":weight fraction, "temp":initial temperature, "h":enthalpy, "elem"element}
        Pc: float,
            Camberpressure, [Pa]
        eps: float,
            Area ratio of nozzle throat & exit, Ae/At
        
        Return
        ----------
        output: dictionary of dictionary
            {"cond": dictionary of calculating condition,
            "therm": dictionary of thermal property,
            "trans": dictionary of transport property,
            "rock": dictionary of rocket property,
            "mole": dictionary of mole fraction}
        """
        Pc = Pc*1e-6    # convert unit Pa to MPa
        make_inp_name(self.cea_dirpath, option, list_species, Pc, eps, fname="tmp")
        output = self._exe_post_process()
        return output

    def onetime_exe_of(self, option, of, Pc, list_oxid, list_fuel, eps):
        """
        Execute CEA in onetime; preparing ".inp", execute CEA and read ".out" file,
        when you assign only oxid and fuel species, and O/F.

        Parameter
        ---------
        option: string
            Calculation option, wtheher using equilibrium composition or frozen composition.
            "equilibrium", "frozen nfz=1" or "frozen nfz=2"
        of: float,
            O/F
        Pc: float,
            Camberpressure, [Pa]
        list_oxid: list,
            The list has an information about oxidizer as dict type; dict{"name":name, "wt":weight fraction, "temp":initial temperature, "h":enthalpy, "elem"element}
        list_fuel: list
            The list has an information about fuel as dict type; dict{"name":name, "wt":weight fraction, "temp":initial temperature, "h":enthalpy, "elem"element}
        eps: float,
            Area ratio of nozzle throat & exit, Ae/At
        
        Return
        ----------
        output: dictionary of dictionary
            {"cond": dictionary of calculating condition,
            "therm": dictionary of thermal property,
            "trans": dictionary of transport property,
            "rock": dictionary of rocket property,
            "mole": dictionary of mole fraction}
        """
        Pc = Pc*1e-6    # convert unit Pa to MPa
        make_inp(self.cea_dirpath, option, of, Pc, list_oxid, list_fuel, eps, fname="tmp")
        output = self._exe_post_process()
        return output

    def _exe_post_process(self):
        self.single_exe(self.cea_dirpath, "tmp")
        cond, therm, trans, rock, mole = Read_output("cea").read_out("tmp")
        output = {"cond": cond,
                  "therm": therm,
                  "trans": trans,
                  "rock": rock,
                  "mole": mole
                  }
        return output


class Read_output:
    """ Read contens in one ".out" file
    
    Parameter
    ---------
    fld_path: string
        folder path which contains .out file
    """
    
    def __init__(self, fld_path):
        self.fld_path = fld_path
    """
    Class to read ".out" file
    """
    cond_param = ["O/F", "Pc", "PHI"]
    therm_param = ["P", "T", "RHO", "H", "U", "G", "S", "M", "Cp", "GAMMAs", "SON", "MACH"]
    rock_param  = ["CSTAR", "CF", "Ivac", "Isp"]
    trans_param = ["VISC", "CONDUCTIVITY", "PRANDTL"]
    
    def _vextract_(self, str_list):
        """
        Extract calculated value from splitted data-list containing raw-data string
        """
        # extract = [str_list[i] for i in range(len(str_list)) if (str_list[i].replace(".","")).replace("-","").isdigit()]
        extract = []
        for i in range(len(str_list)):     # extract only values and eliminate chracters.
            tmp = (str_list[i].replace(".", "")).replace("-", "")
            if tmp.isdigit(): # choose only numerical value
                extract.append(str_list[i])
            elif re.search("\*\*\*", tmp):    # change "*******" to np.nan
                extract.append(np.nan)
            else:
                pass
        val_list = []
        for i in extract:
            if i is np.nan:
                #contain Nan value to val_list
                val_list.append(i)
            elif re.search("-.$", i):
                #calculate exponents, which are sometimes included in "RHO"-row include 
                flag = re.search("-.$", i)
                exp = flag.group()
                base = i.replace(exp, "")
                val_list.append(float(base)*np.power(10, float(exp)))
            elif (float(i)==0.0 and len(i)==1):
                #eliminate eponent "0", which are sometimes included in "RHO"-row include 
                pass
            else:
                val_list.append(float(i))
        return(val_list)
                
    
    def read_out(self, cea_fname):
        """
        Read a ".out" file
        
        Parameters
        ----------
        cea_fname: sting
            file name of ".out" file with out ".out" extension.
        
        Return
        ------
        cond_param: dict {key: elem}
            key: string, a name of condition parameter: e.g. "O/F", "Pc", "PHI"
            elem: list, values of condition parameter
        
        therm_param: dict {key: elem}
            key: string, a name of thermodynamic parameter: e.g. "GAMMAs", "T", "RHO", etc...
            elem: list [c, t, e],
                c: float, a value in chamber
                t: float, a value at throat
                e: float, a value at the end of nozzle

        trans_param: dict {key: elem}
            key: string, a name of transport parameter: e.g. "VISC", "CONDUCTIVITY", "PLANDTL", etc...
            elem: list [t, e]  *in this case "t" and "e" values are equivalent
                t: float, a value at the throat
                e: float, a value at the end of nozzle        

        rock_param: dict {key: elem}
            key: string, a name of rocket parameter: e.g. "CSTR", "Isp", "CF", etc...
            elem: list [t, e]  *in this case "t" and "e" values are equivalent
                t: float, a value at the throat
                e: float, a value at the end of nozzle
        mole_fraction: dict {key: elem}
            key: string, a name of chemical species
            elem: list of float, [c, t, e] *in this case "c" is at chamber, "t" is at throat and "e" is at nozzle exit            
        """
        out_fpath = os.path.join(self.fld_path, cea_fname + ".out")
        file = open(out_fpath,"r")
        
#        cond_param = ["O/F", "Pc", "PHI"]
        cond_param = self.cond_param
        emp_list = ["" for i in range(len(cond_param))]
        cond_param = dict(zip(cond_param, emp_list))
        
#        therm_param = ["P", "T", "RHO", "H", "U", "G", "S", "M", "Cp", "GAMMAs", "SON", "MACH"]
        therm_param = self.therm_param
        emp_list = ["" for i in range(len(therm_param))]
        therm_param = dict(zip(therm_param, emp_list))
        
#        rock_param  = ["Ae/At", "CSTAR", "CF", "Ivac", "Isp"]
        rock_param = self.rock_param
        emp_list = ["" for i in range(len(rock_param))]
        rock_param = dict(zip(rock_param, emp_list))
        
#        rock_param  = ["Ae/At", "CSTAR", "CF", "Ivac", "Isp"]
        trans_param = self.trans_param
        emp_list = ["" for i in range(len(trans_param))]
        trans_param = dict(zip(trans_param, emp_list))
        
        mole_fraction = {}
    
        line = str("null")
        flag_cp = False
        count_trans = 0
        flag_mole = False
        while line:
            line = file.readline()
            warnings.filterwarnings("ignore") # ignore Further Warnings about "empty-string"
            dat = re.split("[\s=]+",line)
            del(dat[0])
            if (len(dat) >= 1):
                del(dat[-1])
            if(len(dat)==0): # empty line
                pass
            else: # not-empty line
                dat_head = dat[0].split(",")[0]
                if(dat_head in cond_param and len(dat)>3):
                    tmp = self._vextract_(dat)
                    cond_param["O/F"] = tmp[0]
                    cond_param["PHI"] = tmp[3]
                elif(dat_head == "Pc"):
                    cond_param["Pc"] = round(float(dat[1])*1.0e-1, 4)
                elif(dat_head in therm_param): #line containing therm_param
                    if (dat_head!="Cp"):
                        therm_param[dat_head] = self._vextract_(dat)
                    elif (dat_head=="Cp" and flag_cp==False):
                        therm_param[dat_head] = self._vextract_(dat)
                        flag_cp = True
                    if (dat_head == "P"):
                        # cond_param["Pc"] = round(therm_param[dat_head][0] *1.0e-1, 4)
                        therm_param[dat_head] = [round(i*1.0e-1, 4) for i in therm_param[dat_head]]
                elif(dat_head in rock_param): #line containing rock_param
                    rock_param[dat_head] = self._vextract_(dat)
                elif(dat_head in trans_param):
                    if (dat_head=="VISC"):
                        trans_param[dat_head] = self._vextract_(dat)
                    elif(count_trans < 3):
                        trans_param[dat_head] = self._vextract_(dat)
                        count_trans += 1
                elif(dat_head == "MOLE"):
                    flag_mole = True
                    continue
                elif(dat_head == "*"):
                    flag_mole = False
                elif(flag_mole):
                    for i in range(len(dat)):
                        tmp_dat = dat[i].replace(".", "")
                        if tmp_dat.isdecimal():
                            tmp_fraction.append(float(dat[i]))
                        else:
                            tmp_fraction = []
                            key = dat[i].strip("*")
                        mole_fraction[key] = tmp_fraction
        file.close()  
        return(cond_param, therm_param, trans_param, rock_param, mole_fraction)


if __name__ == "__main__":
# Following Part is Normal Code for Using This Program
    inst = CEA_execute()
    of, Pc, value_c, value_t, value_e, value_rock, value_mole = inst.all_exe()

# Following Part is for debugging the method of single CEA execute
    # inst2 = CEA_onetime_execute()
    # OPTIOIN = "frozen nfz=2"
    # LIST_SPECIES = [
    #                     {"name": "O2",
    #                     "wt": 50,
    #                     "temp": 290,
    #                     "h": "",
    #                     "elem": ""
    #                     },
    #                     {"name": "N2",
    #                     "wt": 25,
    #                     "temp": 290,
    #                     "h": "",
    #                     "elem": ""
    #                     },
    #                     {"name": "PE",
    #                     "wt": 25,
    #                     "temp": 290,
    #                     "h": -54.2,
    #                     "elem": "C 2 H 4"
    #                     }
    #                 ]
    # PC = 1.0e+6 # [Pa]
    # EPS = 1.0
    # output = inst2.onetime_exe_name(OPTIOIN, LIST_SPECIES, PC, EPS)
    # print(output)

# Following Part is for debugging the method of reading .out file.
#    fld_path = "cea"
#    fname = 'debug'
#    Read = Read_output(fld_path)
#    result = Read.read_out(fname)
#    cond, therm, trans, rock, mole = result
#    print(cond, therm, trans, rock, mole)


