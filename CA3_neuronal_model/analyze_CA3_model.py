import numpy as np
from neuron import h
import neuron
from typing import Tuple
def record_all_init()-> Tuple[h.Vector, h.Vector,h.Vector]:
     vs=h.Vector()
     vs.record(h.apical[14](0.5)._ref_v)
     Ca_rec= h.Vector()
     Ca_rec.record(h.apical[14](0.5)._ref_cai)
     t_rec= h.Vector()
     t_rec.record(h._ref_t)
     return vs, Ca_rec, t_rec
# Load the simulation file
h.load_file("run_model_CA3.hoc")
# Set up the simulation parameters
h.tstop = 15000  # Set simulation time (in ms)
h.dt = 0.025  # Set time step (in ms)
BK_con_vals=np.array([0,0.0005,0.001,0.005,0.01,0.015,0.02,0.025,0.03]) # mh/cm^2
# iterate n (length of BK_con_vals) times
# have a matrix v_gatherer to store v_rec_temp
v_gatherer= []
t_gatherer= []
Ca_gatherer=[]
for bk_con in BK_con_vals:
    print(f"Running simulation for BK_con = {bk_con} mh/cm^2")
    v_rec_temp,Ca_rec_temp,t_temp= record_all_init()
    # Update the BK_con value in the hoc file
    neuron.hoc.execute(f"bk_con={bk_con}")
    neuron.hoc.execute("update_bk_con()")
    # Run the simulation
    h.run()
    # Store the voltage and time data in the gatherer matrices
    v_gatherer.append(np.array(v_rec_temp))
    t_gatherer.append(np.array(t_temp))
    Ca_gatherer.append(np.array(Ca_rec_temp))
# Plotting voltage traces for different BK_con values
import matplotlib.pyplot as plt
plt.figure(figsize=(10,6))
for i, bk_con in enumerate(BK_con_vals):
    plt.plot(t_gatherer[i], v_gatherer[i],
             label=f'BK_con={bk_con:.4f} mho/cm²')

plt.title('Voltage Traces for Different BK Conductance Values')
plt.xlabel('Time (ms)')
plt.ylabel('Membrane Voltage (mV)')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()
# Plotting Calcium traces
plt.figure(figsize=(10,6))
for i, bk_con in enumerate(BK_con_vals):
    plt.plot(t_gatherer[i],Ca_gatherer[i],
             label=f'BK_con={bk_con:.4f} mho/cm²')
plt.title('Intracellular calcium traces for Different BK Conductance Values')
plt.xlabel('Time (ms)')
plt.ylabel('Membrane Voltage (mV)')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()
##########################################
from scipy.io import savemat
savemat("CA3_BK_results.mat", {
    "BK_con_vals": BK_con_vals,
   "v_gatherer": np.array(v_gatherer, dtype=object),
    "t_gatherer": np.array(t_gatherer, dtype=object),
   "Ca_gatherer": np.array(Ca_gatherer, dtype=object)
})
print("Saved BK_con_vals, v_gatherer, t_gatherer, and Ca_gatherer to CA3_BK_results.mat")