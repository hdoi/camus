#!/usr/bin/env python
import copy
import random
import math


def calc_part_num(total_num, monomer_ratio, monomers):
    total_ratio = 0
    for ii in monomer_ratio.items():
        if ii[1] != 'solvent':
            total_ratio += ii[1]
        if ii[1] == 'solvent':
            solvent_name = ii[0]
    solvent_ratio = 100 - total_ratio
    monomer_ratio[solvent_name] = solvent_ratio
    print 'monomer ratio => ',  monomer_ratio

    monomer_dict = {}
    for ii in range(len(monomers)):
        monomer_dict[monomers[ii]['name']] = monomers[ii]

    # calc real particle index
    nmol_list = {}
    nmol_counter = []
    for ratio in monomer_ratio.items():
        target_particle_num = int(total_num * ratio[1] / 100.0)
        monomer_part_num = len(monomer_dict[ratio[0]]['particle'])
        nmol = target_particle_num / monomer_part_num
        nmol_list[ratio[0]] = nmol
        nmol_counter.append(nmol * monomer_part_num)

    sum_nmol_counter = sum(nmol_counter)
    monomer_ratio_real = []
    for ii in range(len(monomer_ratio)):
        monomer_ratio_real.append(
            100 * nmol_counter[ii] / float(sum_nmol_counter))

    return [total_num, monomer_ratio, nmol_list]


def calc_box_length(part_num, density, aspect_ratio):
    """ part_num, density, aspect_ratio = [x, y, z] -> x:y:z """
    box_length = (float(part_num) / float(density) /
                  (aspect_ratio[0] * aspect_ratio[1] * aspect_ratio[2]))**(1.0 / 3.0)
    return [box_length * aspect_ratio[0], box_length * aspect_ratio[1], box_length * aspect_ratio[2]]


def gen_input_file(input_param):
    if input_param[2] == 0:
        input_param[2] = input_param[1] / 1000
    if input_param[2] == 0:
        input_param[2] = 1
    if input_param[3] == 0:
        input_param[3] = input_param[1] / 100
    if input_param[3] == 0:
        input_param[3] = 1

    file_list = input_param[4]
    phys_param = input_param[5]
    f = open(file_list['filename'], 'w')
    print >> f, "&input1"
    print >> f, "  natom                    = ", input_param[0]
    print >> f, "  nstep_total              = ", input_param[1]
    print >> f, "  dump_ene_step            = ", input_param[2]
    print >> f, "  dump_step                = ", input_param[3]
    print >> f, "&end"
    print >> f, "&data_file"
    print >> f, "  init_pos_file            = ", "'" + \
        file_list['init_pos_file'] + "'"
    print >> f, "  init_vel_file            = ", "'" + \
        file_list['init_vel_file'] + "'"
    print >> f, "  aij_file                 = ", "'" + \
        file_list['aij_file'] + "'"
    print >> f, "  particle_file            = ", "'" + \
        file_list['particle_file'] + "'"
    print >> f, "  bond_file                = ", "'" + \
        file_list['bond_file'] + "'"
    print >> f, "  xyz_file                 = ", "'" + \
        file_list['xyz_file'] + "'"
    print >> f, "  ene_file                 = ", "'" + \
        file_list['ene_file'] + "'"
    print >> f, "  vel_file                 = ", "'" + \
        file_list['vel_file'] + "'"
    print >> f, "  pos_file                 = ", "'" + \
        file_list['pos_file'] + "'"
    print >> f, "  input_dump_file          = ", "'" + \
        file_list['input_dump_file'] + "'"
    print >> f, "&end"
    print >> f, "&phys_param"
    print >> f, "  box_length(1)            = ", phys_param['x']
    print >> f, "  box_length(2)            = ", phys_param['y']
    print >> f, "  box_length(3)            = ", phys_param['z']
    print >> f, "  dt                       = ", phys_param['dt']
    print >> f, "  kb                       = ", phys_param['kb']
    print >> f, "  temperature              = ", phys_param['temperature']
    print >> f, "  lambda                   = ", phys_param['lambda']
    print >> f, "  rand_seed                = ", phys_param['rand_seed']
    print >> f, "&end"
    print >> f, "&dpd_parameter"
    print >> f, "  gamma                    = ", phys_param['gamma']
    print >> f, "  cutoff                   = ", phys_param['cutoff']
    print >> f, "&end"
    print >> f, "&calc_flags"
    print >> f, "  use_random_initial_coord = ", phys_param[
        'use_random_initial_coord']
    print >> f, "&end"
    return


def gen_aij_file(aij, filename):
    f = open(filename, 'w')
    aij_name = []
    for ii in range(len(aij)):
        for ij in range(2):
            if (aij[ii][ij] in aij_name) == False:
                aij_name.append(aij[ii][ij])

    # init
    aij_data = []
    for ii in range(len(aij_name)):
        temp = []
        for ij in range(len(aij_name)):
            temp.append(25.0)
        aij_data.append(temp)

    #
    for ii in range(len(aij)):
        part1 = aij_name.index(aij[ii][0])
        part2 = aij_name.index(aij[ii][1])
        aij_data[part1][part2] = aij[ii][2]
        aij_data[part2][part1] = aij[ii][2]

    print >>f, """# 1: comment, 2: "particle type1" "particle type2" aij"""
    # output aij_file
    for ii in range(len(aij_data)):
        for ij in range(ii, len(aij_data[0])):
            print >>f, ii + 1, ij + 1, aij_data[ii][ij]
    return aij_name


def gen_part_file(part_file, part_list):
    # output part file
    f_part = open(part_file, 'w')
    print >> f_part, "# particle index, particle type, fix or move (0-> fix, 1->move)"
    for ii in range(len(part_list)):
        print >> f_part, ii + 1,  part_list[ii][0] + 1, part_list[ii][1]
    return


def lipid_membrane(pos, part_list, sim_box, aij_name, lipid_part, mol_list):
    print lipid_part
    half_z = 0.0
    for ij in range(len(lipid_part)):
        hydrophilic = lipid_part[ij]['hydrophil']
        hydrophobic = lipid_part[ij]['hydrophobic']
        lipid_name = lipid_part[ij]['name']
        direction = 1
        for ii in range(len(mol_list)):
            if mol_list[ii]['name'] == lipid_name:
                now_mol = mol_list[ii]
                for ij in range(now_mol['start'], now_mol['end'] + 1):
                    part_type = part_list[ij][0]
                    part_name = aij_name[part_type]
                    if hydrophilic.count(part_name) != 0:
                        pos[ij][2] = half_z + (2.0) * direction
                    elif hydrophobic.count(part_name) != 0:
                        pos[ij][2] = half_z + (1.0) * direction
                direction *= -1

    return pos

def calc_center(coord):
    center = [0.0, 0.0, 0.0]
    for i in range(len(coord)):
        center += coord[i]
    center = center / float(len(coord))
    return center

def trans_move(orig_coord, mv_vec):
    coord = copy.deepcopy(orig_coord)
    for i in range(len(coord)):
        coord[i] += mv_vec
    return coord

def move_to_origin(coord):
    center = calc_center(coord)
    coord = trans_move(coord, -center)
    return coord

def rand_rotate(coord):
    [a, b, c] = random_rad()
    coord = rotate_xyz(coord, a, b, c)
    return coord

def random_rad():
    return [2.0 * math.pi * random.random(),  2.0 * math.pi * random.random(), 2.0 * math.pi * random.random()]

def dot3(a,b):
    c = [ 
            [ sum( [a[0][i]*b[i][j] for i in range(len(b))] ) for j in range(len(b)) ], 
            [ sum( [a[1][i]*b[i][j] for i in range(len(b))] ) for j in range(len(b)) ], 
            [ sum( [a[2][i]*b[i][j] for i in range(len(b))] ) for j in range(len(b)) ] ]
    return c

def dot2(a,b):
    c = [ sum( [a[i]*b[i][j] for i in range(len(b))] ) for j in range(len(b)) ]
    return c

def rotate_xyz(carte, px, py, pz):
    Rx =    [[1, 0, 0],
             [0, math.cos(px), math.sin(px)],
             [0, -math.sin(px), math.cos(px)]]
    Ry =    [[math.cos(py), 0, -math.sin(py)],
             [0, 1, 0],
             [math.sin(py), 0, math.cos(py)]]
    Rz =    [[math.cos(pz), math.sin(pz), 0],
             [-math.sin(pz), math.cos(pz), 0],
             [0, 0, 1]]
    
    te = dot3( dot3(Rz, Ry),  Rx)
    te = list(map( list, zip(*te)))
    return dot2( carte, te)

def calc_col_radii(filler_pos):
    collision_radii = 0.0
    for ii in range(len(filler_pos)):
        dist =  filler_pos[ii][0]*filler_pos[ii][0] + filler_pos[ii][1]*filler_pos[ii][1] + filler_pos[ii][2]*filler_pos[ii][2]
        dist = dist**(1/2.0)
        if collision_radii < dist :
            collision_radii = dist
    return collision_radii

def calc_dist(coord):
    dist = 0.0
    for ii in range(3):
        dist += coord[ii]*coord[ii]
    dist = dist**(1/2.0)
    return dist

def calc_not_collision_position(filler_center_radii, sim_box, collision_radii):
    collision_flag = True
    dpos = [0,0,0]
    now_count = 0
    while collision_flag == True:
        collision_flag = False
        for ii in range(len(filler_center_radii)):
            dpos = [0,0,0]
            for ij in range(3):
                drand = sim_box[ij] * (random.random() - 0.5) - filler_center_radii[ii][0][ij]
                dpos[ij] = drand - round(drand/sim_box[ij])*sim_box[ij]
            dist = calc_dist(dpos)
            if dist < filler_center_radii[ii][1]+collision_radii:
                collision_flag = True
        now_count += 1
        if now_count == 1000:
            break
    rand_vec = dpos
    return rand_vec

def gen_pos(part_list, sim_box, bond_list, aij_name):
    # gen pos
    pos = []
    filler_center_radii = []
    for ii in range(len(part_list)):
        if ii % 1000 == 0:
            print 'now mol =>', ii, '/', len(part_list)
        if len(part_list[ii]) >= 3: # filler
            if ii == part_list[ii][3]:
                filler_pos = []
                for ij in range(part_list[ii][3], part_list[ii][4]+1):
                    filler_pos.append(part_list[ij][2])
                if len(filler_center_radii) != 0:
                    filler_pos = rand_rotate(move_to_origin(filler_pos))
                else:
                    filler_pos = move_to_origin(filler_pos)
                col_radii = calc_col_radii(filler_pos)
                [x,y,z] = calc_not_collision_position(filler_center_radii, sim_box, col_radii)
                filler_center_radii.append([[x,y,z],col_radii])
                for ij in range(len(filler_pos)):
                    pos.append([filler_pos[ij][0]+x,filler_pos[ij][1]+y,filler_pos[ij][2]+z])
        else:
            [x,y,z] = calc_not_collision_position(filler_center_radii, sim_box, 0.0)
            pos.append([
                sim_box[0] * (random.random() - 0.5),
                sim_box[1] * (random.random() - 0.5),
                sim_box[2] * (random.random() - 0.5)])
            
    
    # same molecule grouping
    for ii in range(len(bond_list)):
        if len(bond_list[ii]) == 2 or len(bond_list[ii]) == 4:  # [atom ID1, atom ID2]
            atom1 = bond_list[ii][0]
            atom2 = bond_list[ii][1]
        elif len(bond_list[ii]) == 6:  # [harmonic or morse, atom ID1, atom ID2, length, spring const, spring const2]
            atom1 = bond_list[ii][1]
            atom2 = bond_list[ii][2]
        if len(part_list[atom1]) == 2: # filler
          for ij in range(3):
              pos[atom2][ij] = pos[atom1][ij]

    return pos


def gen_init_pos_file(init_pos_file, pos):
    # output part file
    f_part = open(init_pos_file, 'w')
    for ii in range(len(pos)):
        print >> f_part, ii + 1, pos[ii][0], pos[ii][1], pos[ii][2]
    return


def gen_init_vel_file(init_vel_file, part_list):
    # output part file
    f_part = open(init_vel_file, 'w')
    for ii in range(len(part_list)):
        print >> f_part, ii + 1, random.random() - 0.5, random.random() - 0.5, random.random() - 0.5
    return


def gen_bond_file(bond_file, bond_list):
    # output bond file
    f_bond = open(bond_file, 'w')
    print >> f_bond, "# bond type (1-> harmonic, 2-> morse), particle index A, particle index B, length, spring const, spring const2 ( morse bond only ), cutoff"
    for ii in range(len(bond_list)):
        if len(bond_list[ii]) == 2:  # harmonic bond [atom ID1, atom ID2]
            print >> f_bond, 1, bond_list[ii][0] + \
                1, bond_list[ii][1] + 1, 0.86, 4.0, 0.0, 10000.0
        elif len(bond_list[ii]) == 4:  # harmonic bond [atom ID1, atom ID2, length, spring const]
            print >> f_bond, 1, bond_list[ii][0] + \
                1, bond_list[ii][1] + 1, bond_list[ii][2], bond_list[ii][3], 0.0, 10000.0
        elif len(bond_list[ii]) == 6:  # [harmonic or morse, atom ID1, atom ID2, length, spring const, spring const2]
            print >> f_bond, \
                bond_list[ii][0], bond_list[ii][1] + 1, bond_list[ii][2] + 1, bond_list[ii][3], bond_list[ii][4], bond_list[ii][5], 10000.0
        elif len(bond_list[ii]) == 7:  # [harmonic or morse, atom ID1, atom ID2, length, spring const, spring const2]
            print >> f_bond, \
                bond_list[ii][0], bond_list[ii][1] + 1, bond_list[ii][2] + 1, bond_list[ii][3], bond_list[ii][4], bond_list[ii][5], bond_list[ii][6]
    return


def gen_particle(total_num, monomers, num_of_mol, aij_name):
    monomer_dict = {}
    for ii in range(len(monomers)):
        monomer_dict[monomers[ii]['name']] = monomers[ii]

    # nmol_list = [] # [10, 10]
    # nmol_name = [] # ['molA', 'molB']

    part_list = []
    bond_list = []
    mol_list = []
    part_counter = 0
    for ii in enumerate(num_of_mol):
        mol_name = ii[1]
        nmol = num_of_mol[mol_name]
        for imol in range(nmol):
            if part_counter % 100000 == 0:
                print 'now particle => ', part_counter
            # particle section
            mol_particle = monomer_dict[mol_name]['particle']
            natom = len(mol_particle)
            for ij in range(natom):
                iatom = aij_name.index(mol_particle[ij])
                ifix = monomer_dict[mol_name]['move'][ij]
                try:
                    coord = monomer_dict[mol_name]['coord'][ij]
                    part_list.append([iatom, ifix, coord, part_counter, part_counter+natom-1])
                except:
                    part_list.append([iatom, ifix])

            # molecule list section
            mol_list.append({'name': mol_name, 'start': part_counter, 'end': part_counter + natom - 1})

            # bond section
            bond = monomer_dict[mol_name]['bond']
            for ij in range(len(bond)):
                if len(bond[ij]) == 2:  # [atom ID1, atom ID2]
                    temp_bond = [bond[ij][0] + part_counter,
                                 bond[ij][1] + part_counter]
                elif len(bond[ij]) == 4:  # [atom ID1, atom ID2, length, spring const]
                    temp_bond = copy.deepcopy(bond[ij])
                    temp_bond[0] += part_counter
                    temp_bond[1] += part_counter
                elif len(bond[ij]) == 6:  # [harmonic or morse, atom ID1, atom ID2, length, spring const, spring const2]
                    temp_bond = copy.deepcopy(bond[ij])
                    temp_bond[1] += part_counter
                    temp_bond[2] += part_counter
                elif len(bond[ij]) == 7:  # [harmonic or morse, atom ID1, atom ID2, length, spring const, spring const2, cutoff]
                    temp_bond = copy.deepcopy(bond[ij])
                    temp_bond[1] += part_counter
                    temp_bond[2] += part_counter

                bond_list.append(temp_bond)
            part_counter += len(monomer_dict[mol_name]['particle'])
    return [part_counter, part_list, bond_list, mol_list]


def aij_part_bond_file_gen(aij, monomers, num_of_mol, total_num, density, file_list, use_random_init_coord, lipid_condition):

    # aij file
    aij.sort()
    aij_name = gen_aij_file(aij, file_list['aij_file'])

    # particle file
    [total_num, part_list, bond_list, mol_list] = gen_particle(
        total_num, monomers, num_of_mol, aij_name)

    # part file
    gen_part_file(file_list['particle_file'], part_list)
    # bond file
    gen_bond_file(file_list['bond_file'], bond_list)

    # use_random_init_coord : '.true.'-> init pos and vel with camus
    # use_random_init_coord : '.false.' -> init pos and vel with this script
    if use_random_init_coord == '.false.':
        # gen init coord
        sim_box = calc_box_length(total_num, density, [1, 1, 1])
        pos = gen_pos(part_list, sim_box, bond_list, aij_name)
        if len(lipid_condition) != 0:
            pos = lipid_membrane(pos, part_list, sim_box, aij_name,
                                 lipid_condition, mol_list)
        gen_init_pos_file(file_list['init_pos_file'], pos)

        # gen init vel
        gen_init_vel_file(file_list['init_vel_file'], part_list)
    return total_num


def convert_3let1let_1let3let(name):
    amino_acid_data = {
        "Ala": "A",  "A": "Ala",
        "Arg": "R",  "R": "Arg",
        "Asn": "N",  "N": "Asn",
        "Asp": "D",  "D": "Asp",
        "Cys": "C",  "C": "Cys",
        "Glu": "E",  "E": "Glu",
        "Gln": "Q",  "Q": "Gln",
        "Gly": "G",  "G": "Gly",
        "His": "H",  "H": "His",
        "Ile": "I",  "I": "Ile",
        "Leu": "L",  "L": "Leu",
        "Lys": "K",  "K": "Lys",
        "Met": "M",  "M": "Met",
        "Phe": "F",  "F": "Phe",
        "Pro": "P",  "P": "Pro",
        "Ser": "S",  "S": "Ser",
        "Thr": "T",  "T": "Thr",
        "Trp": "W",  "W": "Trp",
        "Tyr": "Y",  "Y": "Tyr",
        "Val": "V",  "V": "Val"}
    return amino_acid_data[name]

