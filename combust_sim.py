# -*- coding: utf-8 -*-
"""
Created on Sat Aug  4 07:22:22 2018

@author: TJ
"""

class Gen_excond_table:
    def __init__(self, d, N, Df, eta, rho_f, Rm, Tox, C1, C2, cea_path):
        self.a = 1 - N*np.power(d, 2)/np.power(Df, 2)
        self.Df = Df*1.0e-3
        self.eta = eta
        self.rho_f = rho_f
        self.Rm = Rm
        self.Tox = Tox
        self.C1 = C1
        self.C2 = C2
        self.cea_path = cea_path
        self.mox_range, self.Dt_range = self._input_range_() #generate the calculation range of mox[g/s] and Dt [mm]
#        self.table = self.gen_table()

    def _input_range_(self):
        print("\n\nInput the range of mox [g/s], oxidizer mass flow rate, where you want to generate the table." )
        print("\ne.g. If the range is 10.0 to 20.0 g/s and the interval is 1.0 g/s\n10.0 20.0 1.0")
        tmp = list(map(lambda x: float(x) ,input().split()))
        mox_range = np.arange(tmp[0], tmp[1], tmp[2])*1.0e-3 # generate range and convert [g/s] to [kg/s]

        print("\n\nInput the range of Dt [mm], nozzle throat diameter, where you want to generate the table." )
        print("\ne.g. If the range is 5.0 to 10.0 mm and the interval is 1.0 mm\n5.0 10.0 1.0")
        tmp = list(map(lambda x: float(x) ,input().split()))
        Dt_range = np.arange(tmp[0], tmp[1], tmp[2])*1.0e-3 # generate range and convert [mm] to [m]
        return(mox_range, Dt_range)

    def gen_table(self):
        df = pd.DataFrame({}, index=self.mox_range)
        for Dt in tqdm(self.Dt_range, desc="Dt loop", leave=True):
            Pc = np.array([])
            for mox in tqdm(self.mox_range, desc="mox loop", leave=False):
                tmp = self.converge_Pc(mox, Dt, Pc_init=1.0e+6)
                Pc = np.append(Pc, tmp)
            df[Dt] = Pc # insert calculated Pc to data frame, df.
        return(df)
        

    def converge_Pc(self, mox, Dt, Pc_init=1.0e+6):

        try:
            Pc = optimize.newton(self._iterat_func_Pc_, Pc_init, maxiter=10, tol=1.0e-3, args=(mox, Dt))
        except:
            Pc = optimize.brentq(self._iterat_func_Pc_, 1.0e+4, 10e+6, maxiter=50, xtol=1.0e-3, args=(mox, Dt), full_output=False)
        return(Pc)
    
    def _iterat_func_Pc_(self, Pc, mox, Dt):
        Vox = func_Vox(Pc, mox, self.Rm, self.Tox, self.Df, self.a)
        Vf = func_Vf(Vox, Pc, self.C1, self.C2, n=1.0)
        of = func_of(Vox, Vf, self.rho_f, self.Df, self.a)
        func_cstr = gen_func_cstr(self.cea_path)
        Pc_cal = func_Pc_cal(of, Pc, mox, func_cstr, Dt, self.eta)
        diff = Pc_cal - Pc
        error = diff/Pc_cal
        return(error)

def func_Vox(Pc, mox, Rm, Tox, Df, a):
    Af = np.pi*np.power(Df, 2)/4
    Vox = mox*Rm*Tox/(Pc*Af*(1-a)) #[m/s]
    return(Vox)

def func_Vf(Vox, Pc, C1, C2, n=1.0):
    Vf = (C1/Vox + C2)*np.power(Pc, n) #[m/s]
    return(Vf)
    
def func_of(Vox, Vf, rho_f, Df, a):
    Af = np.pi*np.power(Df, 2)/4
    of = Vox/(rho_f*Af*a*Vf)
    return(of)

def gen_func_cstr(cea_path):
    func = Read_datset(cea_path).gen_func("CSTAR")
    return(func)

def func_Pc_cal(of, Pc, mox, func_cstr, Dt, eta):
    cstr = func_cstr(of,Pc*1.0e-6)[0]
    At = np.pi*np.power(Dt, 2)/4
    Pc_cal = eta*cstr*mox*(1 + 1/of)/At
    return(Pc_cal)
   

if __name__ == "__main__":
#    cea_path = Cui_input().cea_path
    d = 0.3 #[mm]
    N = 433 #[個]
    Df = 38 #[mm]
    eta = 0.9 #[-]
    rho_f = 1191 #[kg/m^3]
    Rm = 259.8 #[J/kg/K]
    Tox = 280 #[K]
    C1 = 9.34e-8 # SI-unit
    C2 = 2.46e-9 # SI-unit  
    instance = Gen_excond_table(d, N, Df, eta, rho_f, Rm, Tox, C1, C2, cea_path)
    
    table = instance.gen_table()
    table.index = np.round(table.index*1.0e+3, 3) # convert the unit of mox to [g/s]
    table = table*1.0e-6 #convert the unit of Pc to [MPa]
    table.columns = np.round(table.columns*1.0e+3, 3) #convert the unit of Dt to [mm]
    
# =============================================================================
#     import matplotlib.pylab as plt
#     mox = 10.0e-3 #[kg/s]
#     Dt = 7.0e-3 #[mm]
#     Pc_range = np.arange(0.2, 1.0, 0.1)*1.0e+6
#     error = [instance._iterat_func_Pc_(x, mox, Dt) for x in Pc_range]
#     plt.plot(Pc_range*1.0e-6, error)
# =============================================================================
