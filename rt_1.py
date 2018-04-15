# -*- coding: utf-8 -*-
"""
Created on Sat Apr 14 13:21:33 2018

@author: T.J.-LAB-PC
"""

def main(ex_df, db_of, db_Pc, func_cstr, func_gamma):
    pass

if __name__ == "__main__":
    import RockCombstAnly
    inst = RockCombstAnly.Cui_input()
    db_of = RockCombstAnly.RT(inst).of
    db_Pc = RockCombstAnly.RT(inst).Pc
    ex_df = RockCombstAnly.RT(inst).ex_df
    func_cstr = RockCombstAnly.RT(inst).cstr
    func_gamma = RockCombstAnly.RT(inst).gamma
    
    main(ex_df, db_of, db_Pc, func_cstr, func_gamma)