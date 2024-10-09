from particle import Particle

import settings
settings.init()  

from ast import literal_eval
import numpy as np
import random
from random import gauss
import os
import math
import sys

from espressomd import checkpointing
# import espressomd.virtual_sites

class Polymer(Particle):

    def __init__(self):
        print("A new Polymer has been created!")

    def load(self):
        pass

    def bake(self, sb, length, N, periodicity):
        """
        [creates N polymers of length (with bonds)]

        """

        fene = sb.param._max_bond_ext

        self.ids = []
        self.positions = []
        self.bonds = []

        def make_rand_unit_vector(dims):
            vec = [random.uniform(-1, 1) for i in range(dims)]
            mag = sum(x**2 for x in vec) ** .5
            return [x/mag for x in vec]

        def is_distance_at_least_one(points, point):
            def euclidean_distance(p1, p2):
                return math.sqrt((p1[0] - p2[0])**2 + (p1[1] - p2[1])**2 + (p1[2] - p2[2])**2)
            
            for p in points:
                if euclidean_distance(p, point) < 1:
                    return False
            return True

        if periodicity == 3:

            #random
            for i in range(N):
                #get a random number in range of box_len
                x = random.uniform(0, sb.box_len[0])
                y = random.uniform(0, sb.box_len[1])
                z = random.uniform(20, sb.box_len[2]-20)

                self.positions.append([x, y, z])
                for j in range(length-1):
                    
                    #generate random vector and multiply it by fene
                    done = False
                    while not done:
                        vec = make_rand_unit_vector(3)
                        # vec = [x*fene*0.5 for x in vec]

                        x += vec[0]
                        y += vec[1]
                        z += vec[2]

                        if z > 10 and z < sb.box_len[2]-10 and is_distance_at_least_one(self.positions, [x, y, z]):
                            done = True

                    self.positions.append([x, y, z])



            # #aligned in x direction
            # xyzs = []
            # for y_i in range(int(sb.box_len[1]/(fene*2))):
            #     for z_i in range(int((sb.box_len[2])/(fene*1.5))):
            #         for x_i in range(int(sb.box_len[0]/((length*1.5)*fene))):
            #             xyzs.append([x_i, y_i, z_i])

            # print(len(xyzs))
                    
            # while len(self.positions) < (N-1)*length:

            #     #pickup randmly and uniformly from yzs and remove it from yzs
            #     xyz = random.choice(xyzs)
            #     xyzs.remove(xyz)
            #     x_i = xyz[0]
            #     y_i = xyz[1]
            #     z_i = xyz[2]

            #     for i in range(length):
            #         self.positions.append([x_i*((length*1.5)*fene)+fene*i*0.7, y_i*2*fene+fene, z_i*1.5*fene+fene])



        elif periodicity == 2:
            posible_starts = []
            for y_i in range(int(sb.box_len[1]/(fene*2.9))):
                for z_i in range(int((sb.box_len[2]-4*sb.param._gap)/(fene*2))):
                    #generate random number between 0.8 and 1.2
                    rnd_x = random.uniform(0.8, 1.2)
                    for x_i in range(int(sb.box_len[0]/((length*rnd_x)*fene))):
                        posible_starts.append((x_i, y_i,z_i, rnd_x))

            #take random element from possible_starts
            for i in range(N):
                x_i, y_i, z_i, rnd_x = random.choice(posible_starts)
                posible_starts.remove((x_i, y_i,z_i, rnd_x))

                for i in range(length):
                    self.positions.append([x_i*((length*rnd_x)*fene)+fene*i*0.6 + (length*rnd_x)*fene*0.5, 
                                           y_i*2.9*fene+fene, 
                                           z_i*2*fene+fene + sb.param._gap])

            # for x_i in range(int(sb.box_len[0]/((length*2)*fene))):
            
            #     if len(self.positions) < N*length:
            #         for i in range(length):
            #             self.positions.append([x_i*((length*2)*fene)+fene*i*0.7, y_i*2.9*fene+fene, z_i*2*fene+fene + sb.param._gap])
            #     else:
            #         break
        
        for p_id in range(N):
            for bead in range(length):
                self.ids.append(p_id*length+bead)
                if bead != 0:
                    self.bonds.append([(0, p_id*length+bead-1)])
                else:
                    self.bonds.append([])


    # def add_calc_cm(self, sb, new_continue_sim):

    #     system = sb.s 
    #     p = sb.param

    #     tot_mass = p._mass_m + p._mass_gel

    #     def calc_cm(self):
    #         cm_mnp = system.analysis.center_of_mass(p._mnp_type)
    #         cm_bead = system.analysis.center_of_mass(p._p_type)

    #         cm = [0]*3

    #         #TODO: calculate CM more adequate and do it for many particles

    #         for idx, (mnp, b) in enumerate(zip(cm_mnp, cm_bead)):
    #             cm[idx] = (mnp*p._mass_m + b*p._mass_gel)/tot_mass

    #         return cm

    #     sb.calc_cm = calc_cm.__get__(sb)


    def initialize(self, sb, new_continue_sim):

        def make_rand_vector(dims):
            vec = [gauss(0, 1) for i in range(dims)]
            mag = sum(x**2 for x in vec) ** .5
            return [x/mag for x in vec]

        system = sb.s 
        p = sb.param

        self.new_ids = []

        if new_continue_sim == "new":

            #creates gel (with bonds) in system

            max_id = system.part.highest_particle_id

            if max_id == -1:

                system.non_bonded_inter[p._p_type, p._p_type].lennard_jones.set_params(
                    epsilon=p._lj_eps, sigma=p._lj_sigma, cutoff=p._lj_cut, shift=p._lj_shift)

                system.non_bonded_inter[p._p_type, p._wall_type].wca.set_params(
                    epsilon=p._lj_eps*2, sigma=2*p._lj_sigma
                )

                #how to calculate CM of this particle

                tot_mass = p._mass_gel

                # def calc_cm(self):
                #     cm_bead = system.analysis.center_of_mass(p._p_type)

                #     cm = [0]*3

                #     #TODO: calculate CM more adequate and do it for many particles

                #     for idx, (mnp, b) in enumerate(zip(cm_mnp, cm_bead)):
                #         cm[idx] = (mnp*p._mass_m + b*p._mass_gel)/tot_mass

                #     return cm

                # sb.calc_cm = calc_cm.__get__(sb)


            for id, pos in enumerate(self.positions):
                sys_id = id + max_id + 1
                system.part.by_id(sys_id).pos = (pos[0], pos[1], pos[2])
                self.new_ids.append(sys_id)
            
                system.part.by_id(sys_id).type = 0
                system.part.by_id(sys_id).rinertia = np.ones(3) * p._momI_gel
                system.part.by_id(sys_id).mass = p._mass_gel
                system.part.by_id(sys_id).dip = (0.0, 0.0, 0.0)
                system.part.by_id(sys_id).rotation = (1, 1, 1)
                # system.part.by_id(sys_id).quat = [ 0.5, -0.5,  0.5, -0.5]
                
                for b in self.bonds[id]:
                    if len(b) == 2:
                        system.part.by_id(sys_id).add_bond((b[0], b[1]))

            #add LJ(WCA)


