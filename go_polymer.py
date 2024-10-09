#Playground for naughty kids
from simulation_box import SimulationBox
from particle.polymer import Polymer
import espressomd
import numpy as np
import functools
import sys 
import io
import os
import time


# For utf-8!
sys.stdout = io.TextIOWrapper(sys.stdout.detach(), encoding="utf-8")
sys.stderr = io.TextIOWrapper(sys.stderr.detach(), encoding="utf-8")

# Immediate output!
print = functools.partial(print, flush=True)

sb = SimulationBox()
polymer = Polymer()

max_ts = 0.00000001

if sb.new_continue_sim == "new":

    if sb.param._N_part != 0:
        sb.add_polymers(polymer, periodicity = 2)

    if sb.thermostat == 'lb':
        sb.init_thermostat(zero_kbT_B = False, zero_kbT_lb = False, act_on_virtual=False)
        sb.add_walls_and_shear_flow(shear_flow=True)
    elif sb.thermostat == 'lang':
        sb.init_thermostat_lang()
    elif sb.thermostat == 'brownian':
        sb.init_thermostat_brownian()

    print("Adding particles..")

    
    sb.print_parameters()

    if sb.thermostat == 'lb' and sb.param._N_part != 0:
        sb.init_part_vel_form_interp_lb()

    sb.s.time = 0
    sb.init_vtf()
    sb.write_parameters_to_json()

elif sb.new_continue_sim == "continue":
    sb.additional_afterload_init(polymer)



sb.init_data_containers_and_clock()


print("Running..")
while(sb.num_of_steps - sb.pos_dip_id > 0):
    print(sb.pos_dip_id, sb.s.time)

    sb.param._timesteps = 15#500

    if sb.hmf_freq != 0:
        sb.int_rot_field()
    else:
        sb.run(sb.param._timesteps)
        sb.pos_dip_id += 1

        try:
            sb.write_sim_data(sb.pos_dip_id, checkpoint_step=10)
        except:
            print(str(sys.exc_info()))
            pid = int(str(sys.exc_info()[1]).split(" ")[-1].split(",")[0])
            pid1 = int(str(sys.exc_info()[1]).split(" ")[-2].split(",")[0])

            p = sb.s.part.by_id(pid1)
            p.delete_bond(p.bonds[0])

            with open(sb.project_path + "/brocken_bonds", "a") as file:
                file.write(str(sb.s.time)+ " " + str(pid) + " " + str(pid1)+ "\n")

        sb.be_faster(max_ts = max_ts)


print("Done 👍")
