# -*- coding: utf-8 -*-
"""
Execute CEA calculation
"""

import numpy as np
from scipy import interpolate
import pandas as pd
import matplotlib.pyplot as plt
import os, sys, glob, shutil
import re, copy
import tqdm
from subprocess import*
import warnings
import time


class CEA_execute:
    """
    Class to excecute CEA calculation
    """

    def __init__(self, fld_path=None):
        self.fld_path = fld_path
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
        else:
            sys.exit("There is no such a directory, \n\"{}\"".format(self.fld_path))
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
                df = pd.DataFrame(val_dict[i], index=of, columns=Pc)
                df.to_csv(os.path.join(dbfld_path, i) + ".csv")
        else:
            for i in val_dict:
                df = pd.DataFrame(val_dict[i], index=of, columns=Pc)
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
        os.chdir("cea")
        cea_path = os.path.join(cea_dirpath, "FCEA2.exe")
        command = os.path.join(cea_dirpath,inp_fname) + "\n"
        p = Popen(cea_path, stdin=PIPE, stdout=PIPE)
        p.communicate(input=bytes(command,"utf-8"))
#        time.sleep(0.1)
        p.wait()
        os.chdir("..")
        return

    def onetime_exe(self):
        """
        Execute CEA in onetime; preparing ".inp", execute CEA and read ".out" file
        """
        pass
       
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
          
               
        for i, fname in enumerate(tqdm.tqdm(inp_list)):
            shutil.copy(os.path.join(inpfld_path,fname+".inp"), os.path.join(cadir,"cea","tmp.inp"))
            self.single_exe(cea_dirpath, "tmp")
            shutil.copy(os.path.join(cea_dirpath, "tmp.out"), os.path.join(outfld_path, fname+".out"))
            cond, therm, trans, rock = Read_output("cea").read_out("tmp")
            
            therm.update(trans) #combine dict "therm" and dict "trans"

            if i ==0: 
                #initialize
                of = []
                Pc = [] 
                value_c = copy.deepcopy(therm)
                value_t = copy.deepcopy(therm)
                value_e = copy.deepcopy(therm)
                for j in therm:
                    value_c[j] = np.empty((0,0), float)
                    value_t[j] = np.empty((0,0), float)
                    value_e[j] = np.empty((0,0), float)                    
                value_rock = copy.deepcopy(rock)
                for j in rock:
                    value_rock[j] = np.empty((0,0), float)
                
            if cond["O/F"] not in of:
                #extend row of array when o/f is renewed
                of.append(cond["O/F"])
                for j in therm:
                    value_c[j] = np.append(value_c[j], np.empty((1,value_c[j].shape[1]), float), axis=0)
                    value_t[j] = np.append(value_t[j], np.empty((1,value_t[j].shape[1]), float), axis=0)
                    value_e[j] = np.append(value_e[j], np.empty((1,value_e[j].shape[1]), float), axis=0)
                for j in rock:
                    value_rock[j] = np.append(value_rock[j], np.empty((1,value_rock[j].shape[1]), float), axis=0)

            if cond["Pc"] not in Pc:
                #extend column of array when Pc is renewed
                Pc.append(cond["Pc"])
                for j in therm:
                    value_c[j] = np.append(value_c[j], np.empty((value_c[j].shape[0],1), float), axis=1)
                    value_t[j] = np.append(value_t[j], np.empty((value_t[j].shape[0],1), float), axis=1)
                    value_e[j] = np.append(value_e[j], np.empty((value_e[j].shape[0],1), float), axis=1)
                for j in rock:
                    value_rock[j] = np.append(value_rock[j], np.empty((value_rock[j].shape[0],1), float), axis=1)

            p = of.index(cond["O/F"])
            q = Pc.index(cond["Pc"])
            for j in therm:
                #Substitute each themodynamic value
                value_c[j][p,q] = therm[j][0]
                value_t[j][p,q] = therm[j][1]
                value_e[j][p,q] = therm[j][2]
            for j in rock:
                #Substitute each rocket-parameter value
                value_rock[j][p,q] = rock[j][1]
                
        self._csv_out_(dbfld_path, of, Pc, value_c, point="c") #write out in csv-file
        self._csv_out_(dbfld_path, of, Pc, value_t, point="t") #write out in csv-file
        self._csv_out_(dbfld_path, of, Pc, value_e, point="e") #write out in csv-file
        self._csv_out_(dbfld_path, of, Pc, value_rock, point="") #write out in csv-file
            
        return(of, Pc, value_c, value_t, value_e, value_rock)


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
        extract = [str_list[i] for i in range(len(str_list)) if (str_list[i].replace(".","")).replace("-","").isdigit()]
        val_list = []
        for i in extract:
            if re.search("-.$", i):
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
    
        line = str("null")
        flag_cp = False
        count_trans = 0
        flag_mole = False
        while line:
            line = file.readline()
            warnings.filterwarnings("ignore") # ignore Further Warnings about "empty-string"
            dat = re.split("[\s=]*",line)
            del(dat[0])
            if(len(dat)==0): # empty line
                pass
            else: # not-empty line
                dat_head = dat[0].split(",")[0]
                if(dat_head in cond_param and len(dat)>3):
                    tmp = self._vextract_(dat)
                    cond_param["O/F"] = tmp[0]
                    cond_param["PHI"] = tmp[3]
                elif(dat_head in therm_param): #line containing therm_param
                    if (dat_head!="Cp"):
                        therm_param[dat_head] = self._vextract_(dat)
                    elif (dat_head=="Cp" and flag_cp==False):
                        therm_param[dat_head] = self._vextract_(dat)
                        flag_cp = True
                    if (dat_head == "P"):
                        cond_param["Pc"] = round(therm_param[dat_head][0] *1.0e-1, 4)
                        therm_param[dat_head] = [round(i*1.0e-1, 4) for i in therm_param[dat_head]]
                elif(dat_head in rock_param): #line containing rock_param
                    rock_param[dat_head] = self._vextract_(dat)
                elif(dat_head in trans_param):
                    if (dat_head=="VISC"):
                        trans_param[dat_head] = self._vextract_(dat)
                    elif(count_trans < 3):
                        trans_param[dat_head] = self._vextract_(dat)
                        count_trans += 1
#                elif(dat_head == "MOLE"):
#                    flag = True
#                elif
                    
        file.close()

    #    therm_ntpl = collections.namedtuple("thermval",["c","t","e"])    
    #    rock_ntpl = collections.namedtuple("rockval",["t","e"])
    #    P = therm_ntpl(c=therm_param["P"][0], t=therm_param["P"][1], e=therm_param["P"][2])
    #    T = therm_ntpl(c=therm_param["T"][0], t=therm_param["T"][1], e=therm_param["T"][2])
    #    rho = therm_ntpl(c=therm_param["RHO"][0], t=therm_param["RHO"][1], e=therm_param["RHO"][2])
    #    H = therm_ntpl(c=therm_param["H"][0], t=therm_param["H"][1], e=therm_param["H"][2])
    #    U = therm_ntpl(c=therm_param["U"][0], t=therm_param["U"][1], e=therm_param["U"][2])
    #    G = therm_ntpl(c=therm_param["G"][0], t=therm_param["G"][1], e=therm_param["G"][2])
    #    S = therm_ntpl(c=therm_param["S"][0], t=therm_param["S"][1], e=therm_param["S"][2])
    #    M = therm_ntpl(c=therm_param["M"][0], t=therm_param["M"][1], e=therm_param["M"][2])
    #    Cp = therm_ntpl(c=therm_param["Cp"][0], t=therm_param["Cp"][1], e=therm_param["Cp"][2])
    #    gamma = therm_ntpl(c=therm_param["GAMMAs"][0], t=therm_param["GAMMAs"][1], e=therm_param["GAMMAs"][2])
    #    SON = therm_ntpl(c=therm_param["SON"][0], t=therm_param["SON"][1], e=therm_param["SON"][2])
    #    MACH = therm_ntpl(c=therm_param["MACH"][0], t=therm_param["MACH"][1], e=therm_param["MACH"][2])
    #    Eps = rock_ntpl(t=rock_param["Ae/At"][0], e=rock_param["Ae/At"][1])
    #    cstr = rock_ntpl(t=rock_param["CSTAR"][0], e=rock_param["CSTAR"][1])
    #    CF = rock_ntpl(t=rock_param["CF"][0], e=rock_param["CF"][1])
    #    Ispvac = rock_ntpl(t=rock_param["Ivac"][0], e=rock_param["Ivac"][1])
    #    Isp = rock_ntpl(t=rock_param["Isp"][0], e=rock_param["Isp"][1])    
    
        return(cond_param, therm_param, trans_param, rock_param)


if __name__ == "__main__":
    inst = CEA_execute()
    of, Pc, value_c, value_t, value_e, value_rock = inst.all_exe()

