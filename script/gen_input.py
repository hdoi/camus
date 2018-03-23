#!/usr/bin/env python
""" import my lib """
import sys
import lib
import copy
import filler_gen
sys.path.append('.')


def convert_chi_to_a(chi):
    """ convert chi parameter to aij parameter """
    for i in range(len(chi)):
        chi[i][2] = 25.0 + chi[i][2] / 0.306
    return chi

if __name__ == "__main__":
    argvs = sys.argv
    argc = len(argvs)
    if (argc < 2):
        print "./gen_input.py monomer_model/calc_setting"
        quit()
    exec(open(argvs[1], 'r').read())
    # read aij
    exec(open(aij_file, 'r').read())
    # read monomer_list.dat
    exec(open(monomer_file, 'r').read())
    try:
        for ii in range(len(filler_monomer)):
            monomers.append(filler_monomer[ii])
    except NameError:
        pass

    ####### setting #######
    for ii in range(len(ratio_list)):
        if len(total_num_list) == 1:
            total_num = total_num_list[0]
        else:
            total_num = total_num_list[ii]
        if len(step_list) == 1:
            step = step_list[0]
        else:
            step = step_list[ii]

        # gen file name
        name_tail = name_tail_list[ii]
        file_list = copy.deepcopy(default_file)
        for ij in enumerate(default_file):
            file_list[ij[1]] = file_list[ij[1]] + name_tail

        [total_num, monomer_ratio, num_of_mol] = lib.calc_part_num(
            total_num, ratio_list[ii], monomers)
        print '[total_num , monomer_ratio, num_mol] = ', total_num, monomer_ratio, num_of_mol

        try:
            lipid_condition
        except NameError:
            lipid_condition = []
        total_num = lib.aij_part_bond_file_gen(aij, monomers, num_of_mol, total_num, density, file_list, phys_param[
                                               'use_random_initial_coord'], lipid_condition)

        [phys_param['x'], phys_param['y'], phys_param['z']
         ] = lib.calc_box_length(total_num, density, [1, 1, 1])
        input_param = [total_num, step, dump_ene, dump, file_list, phys_param]
        lib.gen_input_file(input_param)
