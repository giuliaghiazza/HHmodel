import matplotlib.pyplot as plt
import numpy as np
import math
    
#Parameters for potassium and sodium
def alpha_n(v):
    v = v*1000 #in Volt
    a = 0.01 * (-v-55)/( math.exp((-v-55)/10.0) -1) * 1000 #in per S
    return a
    
def beta_n(v):
    v = v*1000
    b = 0.125*math.exp((-v-65)/80.0) * 1000
    return b
    
def alpha_m(v):
    v = v*1000
    am = 0.1 * (-v-40)/( math.exp((-v-40)/10.0) -1) * 1000
    return am
    
def beta_m(v):
    v = v*1000
    bm = 4*math.exp((-v-65)/18.0) * 1000
    return bm
    
def alpha_h(v):
    v = v*1000
    ah = 0.07*math.exp((-v-65)/20.0) * 1000
    return ah
    
def beta_h(v):
    v = v*1000
    bh = 1/( math.exp((-v-35)/10.0) +1) * 1000
    return bh
    
    
# Values
dt     = 10E-6; 
Cm     = 100E-12;   
v_init = -80E-3;   
    
# HH 
# Initial conditions
n = alpha_n(v_init)/(alpha_n(v_init)+beta_n(v_init))          
m = alpha_m(v_init)/(alpha_m(v_init)+beta_m(v_init))
h = alpha_h(v_init)/(alpha_h(v_init)+beta_h(v_init))
    
Gkmax = 1E-6      # Maximum potassium conductance
Gnamax= 7E-6      # Maximum sodium conductance
Gleak  = 5E-9      # 5 nS conductance
    
Ek     = -80E-3    # K equilibrium potential 
Ena    = 40E-3     # Na equilibrium potential 
Eleak  = -70E-3    # Reversal potential of -60 mV
    
# Injected Current step
current_magnitude = 200E-12; 

i_inj = np.concatenate( (np.zeros([round(0.2/dt),1]),
                            current_magnitude*np.ones([round(0.3/dt), 1]),
                            np.zeros([round(0.5/dt), 1])) )
v_out = np.zeros(np.size(i_inj))
    
for t in range(np.size(v_out)):
    if t == 1:
        v_out[t] = v_init; #Initial condition
    else:
        #Calculate how the particles change
        dn = (alpha_n(v_out[t-1]) * (1 - n) - beta_n(v_out[t-1]) * n) * dt
        n = n + dn
    
        dm = (alpha_m(v_out[t-1]) * (1 - m) - beta_m(v_out[t-1]) * m) * dt
        m = m + dm
    
        dh = (alpha_h(v_out[t-1]) * (1 - h) - beta_h(v_out[t-1]) * h) * dt
        h = h + dh
        
        # Potassium channel conductance
        Gk = Gkmax*n*n*n*n 
        # Potassium current 
        i_k = Gk * (v_out[t-1] - Ek)
    
        # Sodium channel conductance      
        Gna = Gnamax*m*m*m*h
        # Sodium current
        i_na = Gna * (v_out[t-1] - Ena)
        # Ion channels
        i_leak = Gleak * (v_out[t-1] - Eleak)  
    
        i_total = i_inj[t] - i_leak - i_k - i_na
    
        # Capacitor equation
        dv = (i_total/Cm) * dt
        # Add dv to last known voltage
        v_out[t] = v_out[t-1] + dv      
    
#Make the graph
t_vec = np.linspace(0, dt*np.size(v_out), np.size(v_out))
plt.plot(t_vec, v_out)
plt.xlabel('Time (s)')
plt.ylabel('Membrane Potential (V)')
plt.show()