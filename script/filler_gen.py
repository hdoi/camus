#!/usr/bin/env python
import copy

def scale_for_density(filler, density):
    # calc density
    factor = (4.0/density)**(1.0/3.0)
    for ix in range(len(filler)):
        for iy in range(len(filler[0])):
            for iz in range(len(filler[0][0])):
                for ii in range(len(filler[0][0][0]['coord'])):
                    for ij in range(3):
                        filler[ix][iy][iz]['coord'][ii][ij] = factor*filler[ix][iy][iz]['coord'][ii][ij]
    return filler

def delete_particle(delete_part_list, filler, bond_list):
    # particle
    for ix in range(len(filler)):
        for iy in range(len(filler[0])):
            for iz in range(len(filler[0][0])):
                for ipart in range(len(filler[ix][iy][iz]['num'])):
                    for idel in range(len(delete_part_list)):
                        if filler[ix][iy][iz]['num'][ipart] == delete_part_list[idel]:
                            filler[ix][iy][iz]['num'][ipart] = "del"
                            filler[ix][iy][iz]['coord'][ipart][0] = "del"
                            filler[ix][iy][iz]['coord'][ipart][1] = "del"
                            filler[ix][iy][iz]['coord'][ipart][2] = "del"
                for ibond in range(len(filler[ix][iy][iz]['bond'])):
                    for idel in range(len(delete_part_list)):
                        if filler[ix][iy][iz]['bond'][ibond][0] == delete_part_list[idel]:
                            filler[ix][iy][iz]['bond'][ibond] = ['del','del']
                        if filler[ix][iy][iz]['bond'][ibond][1] == delete_part_list[idel]:
                            filler[ix][iy][iz]['bond'][ibond] = ['del','del']
    # bond
    for ibond in range(len(bond_list)):
        for idel in range(len(delete_part_list)):
            if bond_list[ibond][0] == delete_part_list[idel] :
                bond_list[ibond] = ['del', 'del']
            if bond_list[ibond][1] == delete_part_list[idel] :
                bond_list[ibond] = ['del', 'del']

    # delete in filler
    for ix in range(len(filler)):
        for iy in range(len(filler[0])):
            for iz in range(len(filler[0][0])):
                while 'del' in filler[ix][iy][iz]['num']: filler[ix][iy][iz]['num'].remove('del')
                
                for ipart in range(len(filler[ix][iy][iz]['coord'])):
                    while 'del' in filler[ix][iy][iz]['coord'][ipart]: filler[ix][iy][iz]['coord'][ipart].remove('del')
                while [] in filler[ix][iy][iz]['coord']: filler[ix][iy][iz]['coord'].remove([])
                
                for ibond in range(len(filler[ix][iy][iz]['bond'])):
                    while 'del' in filler[ix][iy][iz]['bond'][ibond]: filler[ix][iy][iz]['bond'][ibond].remove('del')
                while [] in filler[ix][iy][iz]['bond']: filler[ix][iy][iz]['bond'].remove([])
    
    # delete in bond_list
    for ibond in range(len(bond_list)):
        while 'del' in bond_list[ibond]: bond_list[ibond].remove('del')
    while [] in bond_list: bond_list.remove([])
    
    return [filler, bond_list]

def cutting_up_for_cylinder_tube(filler, bond_list, inner_radii, outer_radii):
    # calc center
    center = [0,0,0]
    counter = 0
    for ix in range(len(filler)):
        for iy in range(len(filler[0])):
            for iz in range(len(filler[0][0])):
                for ii in range(len(filler[0][0][0]['coord'])):
                    counter += 1
                    center[0] += filler[ix][iy][iz]['coord'][ii][0]
                    center[1] += filler[ix][iy][iz]['coord'][ii][1]
                    center[2] += filler[ix][iy][iz]['coord'][ii][2]
    center[0] = center[0] / float(counter)
    center[1] = center[1] / float(counter)
    center[2] = center[2] / float(counter)

    delete_part_list = []
    for ix in range(len(filler)):
        for iy in range(len(filler[0])):
            for iz in range(len(filler[0][0])):
                for ii in range(len(filler[0][0][0]['coord'])):
                    dist = 0.0
                    for ij in range(2):
                        dr = filler[ix][iy][iz]['coord'][ii][ij] - center[ij]
                        dist += dr*dr
                    dist = dist**(1/2.0)
                    if dist < inner_radii or dist > outer_radii :
                        delete_part_list.append(filler[ix][iy][iz]['num'][ii])
    return delete_part_list

def cutting_up_for_cylinder(filler, bond_list, radius):
    # calc center
    center = [0,0,0]
    counter = 0
    for ix in range(len(filler)):
        for iy in range(len(filler[0])):
            for iz in range(len(filler[0][0])):
                for ii in range(len(filler[0][0][0]['coord'])):
                    counter += 1
                    center[0] += filler[ix][iy][iz]['coord'][ii][0]
                    center[1] += filler[ix][iy][iz]['coord'][ii][1]
                    center[2] += filler[ix][iy][iz]['coord'][ii][2]
    center[0] = center[0] / float(counter)
    center[1] = center[1] / float(counter)
    center[2] = center[2] / float(counter)

    delete_part_list = []
    for ix in range(len(filler)):
        for iy in range(len(filler[0])):
            for iz in range(len(filler[0][0])):
                for ii in range(len(filler[0][0][0]['coord'])):
                    dist = 0.0
                    for ij in range(2):
                        dr = filler[ix][iy][iz]['coord'][ii][ij] - center[ij]
                        dist += dr*dr
                    dist = dist**(1/2.0)
                    if dist > radius :
                        delete_part_list.append(filler[ix][iy][iz]['num'][ii])
    return delete_part_list

def cutting_up_for_sphere(filler, bond_list, radius):
    # calc center
    center = [0,0,0]
    counter = 0
    for ix in range(len(filler)):
        for iy in range(len(filler[0])):
            for iz in range(len(filler[0][0])):
                for ii in range(len(filler[0][0][0]['coord'])):
                    counter += 1
                    center[0] += filler[ix][iy][iz]['coord'][ii][0]
                    center[1] += filler[ix][iy][iz]['coord'][ii][1]
                    center[2] += filler[ix][iy][iz]['coord'][ii][2]
    center[0] = center[0] / float(counter)
    center[1] = center[1] / float(counter)
    center[2] = center[2] / float(counter)

    delete_part_list = []
    for ix in range(len(filler)):
        for iy in range(len(filler[0])):
            for iz in range(len(filler[0][0])):
                for ii in range(len(filler[0][0][0]['coord'])):
                    dist = 0.0
                    for ij in range(3):
                        dr = filler[ix][iy][iz]['coord'][ii][ij] - center[ij]
                        dist += dr*dr
                    dist = dist**(1/2.0)
                    if dist > radius :
                        delete_part_list.append(filler[ix][iy][iz]['num'][ii])
    return delete_part_list

def gen_block(cell):
    x = cell[0]
    y = cell[1]
    z = cell[2]

    unit_cell_pos = []
    unit_cell_pos.append([0.0, 0.0, 0.0])
    unit_cell_pos.append([0.5, 0.5, 0.0])
    unit_cell_pos.append([0.5, 0.0, 0.5])
    unit_cell_pos.append([0.0, 0.5, 0.5])
    unit_cell_bond = [[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]]
    unit_cell_num = [0, 1, 2, 3]
    unit_cell = {'coord': unit_cell_pos,
                 'bond': unit_cell_bond, 'num': unit_cell_num}
    
    filler = []
    counter = 0
    for ix in range(x):
        tempy = []
        for iy in range(y):
            tempz = []
            for iz in range(z):
                temp = copy.deepcopy(unit_cell)
                for ii in range(len(unit_cell_pos)):
                    temp['coord'][ii][0] = unit_cell_pos[ii][0] + ix
                    temp['coord'][ii][1] = unit_cell_pos[ii][1] + iy
                    temp['coord'][ii][2] = unit_cell_pos[ii][2] + iz
                    temp['num'][ii] += counter
                for ii in range(len(unit_cell['bond'])):
                    temp['bond'][ii][0] += counter
                    temp['bond'][ii][1] += counter
                counter += len(unit_cell_pos)
                tempz.append(temp)
            tempy.append(tempz)
        filler.append(tempy)
    
    # bond generate
    bond_list = []
    for ix in range(x):
        for iy in range(y):
            for iz in range(z):
                for ix2 in range(ix, ix + 2):
                    if ix2 == x:
                        continue
                    for iy2 in range(iy, iy + 2):
                        if iy2 == y:
                            continue
                        for iz2 in range(iz, iz + 2):
                            if iz2 == z:
                                continue
                            if [ix, iy, iz] == [ix2, iy2, iz2]:
                                continue
                            bond_list.extend(make_bond(filler[ix][iy][iz], filler[ix2][iy2][iz2]))
    return [filler, bond_list]

def make_bond(filler1, filler2):
    # create bond
    # bond length 0.7071
    new_bond = []
    ref_dist = (0.5**2 + 0.5**2)**(0.5)
    for ii in range(len(filler1['coord'])):
        for ij in range(len(filler2['coord'])):
            dist = 0.0
            for ik in range(3):
                dist += (filler1['coord'][ii][ik] - filler2['coord'][ij][ik])**2
            dist = dist**0.5
            if dist == ref_dist:
                new_bond.append([filler1['num'][ii], filler2['num'][ij]])
    # filler = {'num': [4, 5, 6, 7], 'coord': [[0.0, 0.0, 1.0], [0.5, 0.5, 1.0], [0.5, 0.0, 1.5], [0.0, 0.5, 1.5]], 'bond': [[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]]}
    return new_bond

def calc_new_num(cell, filler, bond_list, id_list):
    x = cell[0]
    y = cell[1]
    z = cell[2]
    for ix in range(x):
        for iy in range(y):
            for iz in range(z):
                if filler[ix][iy][iz]['num'] != []:
                    for ii in range(len(filler[ix][iy][iz]['num'])):
                        filler[ix][iy][iz]['num'][ii] =  id_list.index(filler[ix][iy][iz]['num'][ii])
                    for ii in range(len(filler[ix][iy][iz]['bond'])):
                        if len(filler[ix][iy][iz]['bond'][ii]) != 0:
                            filler[ix][iy][iz]['bond'][ii][0] = id_list.index(filler[ix][iy][iz]['bond'][ii][0])
                            filler[ix][iy][iz]['bond'][ii][1] = id_list.index(filler[ix][iy][iz]['bond'][ii][1])
    for ii in range(len(bond_list)):
        bond_list[ii][0] = id_list.index( bond_list[ii][0] )
        bond_list[ii][1] = id_list.index( bond_list[ii][1] )
    return [filler, bond_list]

def make_id_list(cell, filler):
    x = cell[0]
    y = cell[1]
    z = cell[2]
    id_list = []
    for ix in range(x):
        for iy in range(y):
            for iz in range(z):
                if filler[ix][iy][iz]['num'] != []:
                    for ii in range(len(filler[ix][iy][iz]['num'])):
                        id_list.append(filler[ix][iy][iz]['num'][ii])
    return id_list

def convert_to_monomers(name, cell, filler, bond_list, movable, spring, particle_type):
    filler_monomer = {'name':name,'particle':[],'move':[],'bond':[],'coord':[]}
    id_list = make_id_list(cell, filler)
    # calc new atom id
    [filler, bond_list] = calc_new_num(cell, filler,bond_list, id_list)
    
    x = cell[0]
    y = cell[1]
    z = cell[2]
    for ix in range(x):
        for iy in range(y):
            for iz in range(z):
                # coord copy
                if filler[ix][iy][iz]['num'] != []:
                    for ik in range(len(filler[ix][iy][iz]['num'])):
                        filler_monomer['coord'].append( filler[ix][iy][iz]['coord'][ik] )
                        filler_monomer['move'].append(movable)
                        filler_monomer['particle'].append(particle_type)
                # bond copy
                if len(filler[ix][iy][iz]['bond']) != []:
                    for ik in range(len(filler[ix][iy][iz]['bond'])):
                        filler_monomer['bond'].append(
                        [filler[ix][iy][iz]['bond'][ik][0],filler[ix][iy][iz]['bond'][ik][1],
                        spring[0],spring[1]])
    # extra bond copy
    for ii in range(len(bond_list)):
        filler_monomer['bond'].append([bond_list[ii][0],bond_list[ii][1],spring[0],spring[1]])
    return filler_monomer

def output_monomer_xyz(filler_monomer):
    print len(filler_monomer['coord'])
    print filler_monomer['name']
    for ii in range(len(filler_monomer['coord'])):
        print "He", filler_monomer['coord'][ii][0], filler_monomer['coord'][ii][1], filler_monomer['coord'][ii][2]
    return

def gen_filler_monomer(cell, density, radius, shape, name, movable, spring, particle_type):
    if shape == 'sphere':
        radii = radius[0]
    if shape == 'cylinder':
        radii = radius[0]
    if shape == 'cylinder_tube':
        inner_radii = radius[0]
        outer_radii = radius[1]

    [filler, bond_list] = gen_block(cell)
    filler = scale_for_density(filler, density)

    if shape == 'sphere':
        delete_part_list = cutting_up_for_sphere(filler, bond_list, radii)
    if shape == 'cylinder':
        delete_part_list = cutting_up_for_cylinder(filler, bond_list, radii)
    if shape == 'cylinder_tube':
        delete_part_list = cutting_up_for_cylinder_tube(filler, bond_list, inner_radii, outer_radii)
    if shape != 'block':
        [filler, bond_list] = delete_particle(delete_part_list, filler, bond_list)
    filler_monomer = convert_to_monomers(name, cell, filler, bond_list, movable, spring, particle_type)
    return filler_monomer

if __name__ == "__main__":
    density = 3.0
    cell = [2,2,2]
    radius = [1,2]
    
    for shape in ['block', 'sphere', 'cylinder', 'cylinder_tube']:
    #for shape in ['block']:
        name = shape
        movable = 0 # 0 is fixed. 1 is movable.
        spring = [0.84, 4] # spring length, spring constant
        particle_type = 'X'
        filler_monomer = gen_filler_monomer(cell, density, radius, shape, name, movable, spring, particle_type)
        output_monomer_xyz(filler_monomer)
