# %%
# Importing the data
from pyapep import simsep
import matplotlib.pyplot as plt
import numpy as np

# %%
# H2, N2, CH4, CO, CO2
# K1 (mol/kg )
K1 = np.array([8.519, 6.774, 13.69, 9.614, 8.984])

# K2(mol/kg/K *1E-3)
K2 = np.array([-13.57, -11.57, -31.67, -19.27, -9.867])*1E-3

# K3(1/kPa 1E-6)
K3 = np.array([5.651, 3.778, 38.36, 13.95, 2.266])*1E-6

# K4 (K)
K4 = np.array([640, 1650, 1083, 1496, 3130])

# K5 (dim.less *1E-2)
K5 = np.array([95.91, 76.74, 66.51, -14.61, 36.22])*1E-2

# K6 (K)
K6 = [9.776, 78.83, 97.8, 407, 454.4]

def isomix_cand1(P,T):
    qms_list = []
    Bs_list = []
    n_list = []
    # qms = K1 + K2*T (mol/kg)
    # Bs = K3*exp(K4/T) (kPa ^(-1) )
    # n =K5 + K6/T  (dimensionless)
    T_av = np.mean(T)
    for ii in range(len(P)):
        qms_tmp = K1[ii] + K2[ii]*T_av
        Bs_tmp = K3[ii]*np.exp(K4[ii]/T_av)
        n_tmp = K5[ii]+K6[ii]/T_av
        qms_list.append(qms_tmp)
        Bs_list.append(Bs_tmp)
        n_list.append(n_tmp)
    #print('qms_list', qms_list)
    #print('Bs_list', Bs_list)
    #print('n_list', n_list)
    q_return = []
    for ii, pp in enumerate(P):
        P_kPa = pp*100
        P_kPa[P_kPa < 0] = 0
        try:
            n_inv = 1/n_list[ii]
            numer = qms_list[ii]*(Bs_list[ii]*P_kPa)**n_inv
            denom = 1+ (Bs_list[ii]*P_kPa)**n_inv
        except:
            print(n_list[ii])
        q_tmp = numer/denom
        q_return.append(q_tmp)
    return q_return

P_list = []
P_list.append(1*np.ones([100,]))
P_list.append(1*np.ones([100,]))
P_list.append(0.5*np.ones([100,]))
P_list.append(0.5*np.ones([100,]))
P_list.append(0.5*np.ones([100,]))
T_test = 293

isomix_cand1(P_list, T_test)

# %%
# Test the isotherm model
P_ran = np.linspace(0,10, 51)
P_ran_input = [P_ran, ]*5
T_ran = 308*np.ones([51,])
qH2, qN2, qCH4, qCO, qCO2 = isomix_cand1(P_ran_input, T_ran)
    
# %%
# Test result H2
plt.figure()
plt.plot(P_ran, qH2, 'ko')
plt.xlabel('pressure (bar)')
plt.ylabel('uptake (mol/kg)')
plt.title('H2')

# %%
# Test result N2
plt.figure()
plt.plot(P_ran, qN2, 'ko' )
plt.xlabel('pressure (bar)')
plt.ylabel('uptake (mol/kg)')
plt.title('N2')

# %%
# Test result CH4
plt.figure()
plt.plot(P_ran, qCH4, 'ko')
plt.xlabel('pressure (bar)')
plt.ylabel('uptake (mol/kg)')
plt.title('CH4')

# %%
# Test result CO
plt.figure()
plt.plot(P_ran, qCO, 'ko')
plt.xlabel('pressure (bar)')
plt.ylabel('uptake (mol/kg)')
plt.title('CO')

# %%
# Test result CO2
plt.figure()
plt.plot(P_ran, qCO2, 'ko')
plt.xlabel('pressure (bar)')
plt.ylabel('uptake (mol/kg)')
plt.title('CO2')

# %%
L = 1
D = 0.2
N = 26
A_cros = D**2/4*np.pi
epsi = 0.4      # m^3/m^3
rho_s = 1000    # kg/m^3
dp = 0.02       # m (particle diameter)

n_comp = 5

# %%
# Column Geometry

A_cros = D**2/4*np.pi
epsi = 0.4      # m^3/m^3
# %%
# Define a column    
c1 = simsep.column(L, A_cros,n_comp, N, E_balance=False)

# %%
# Adsorbent info
c1.adsorbent_info(isomix_cand1, epsi, dp, rho_s,)

# %%
# Gas properties
# H2, N2, CH4, CO, CO2
Mw = [2, 28, 14, 28, 44]
mu = [1.81E-5,]*5 # Pa sec (visocisty of gas)

c1.gas_prop_info(Mw,mu)
# %%
##### PRESSURIZATION #####
# %%
# Mass transfer
k_MTC = [2E-3, ]*4 + [2E-3]
#k_MTC = [5E-5, ]*4 + [2E-3]
D_disp = [1E-7, ] *5
a_surf = 1
c1.mass_trans_info(k_MTC, a_surf, D_disp)

# %%
# Boundary conditions
P_out = 1     # Ignore This Value
P_in = 8      # Important (during pressurization)
T_in = 300
y_in = [0.05, 0.55, 0.05, 0.05, 0.30] # H2, N2, CH4, CO, CO2
Cv_in = 1E-4        # m^3/sec/bar
Cv_out = 0          # m^3/sec/bar
u_feed = 0.01            # m/s
Q_in = u_feed*A_cros*epsi  # volumetric flowrate
c1.boundaryC_info(P_out, P_in, T_in, y_in, Cv_in,Cv_out,Q_in, 
                    assigned_v_option = True, 
                    foward_flow_direction=True)

# %%
# Initial conditions
P_init = 1.5*np.ones([N,])

y_init = [1.00*np.ones([N,]),
          0.00*np.ones([N,]),
          0.00*np.ones([N,]),
          0.00*np.ones([N,]),
          0.00*np.ones([N,]),]
P_part = [1.00*np.ones([N,])*P_init,
          0.00*np.ones([N,])*P_init,
          0.00*np.ones([N,])*P_init,
          0.00*np.ones([N,])*P_init,
          0.00*np.ones([N,])*P_init]
Tg_init = 300*np.ones([N,])
Ts_init = 300*np.ones([N,])
q1,q2,q3,q4,q5 = isomix_cand1(P_part, Tg_init)
q_init = [q1,q2,q3,q4,q5]

c1.initialC_info(P_init, Tg_init, Ts_init, y_init, q_init )
print(c1)

# %%
# Run
y_pr, z_pr, t_pr = c1.run_ma_linvel(100,50)
# %%
# Graph of 1st component
c1.Graph(5, 4, 
         yaxis_label= 'Gas Conc. CO2 (mol/m$^3$)')
# %%
c1.Graph(5, 9, 
         yaxis_label= 'Solid uptake CO2 (mol/kg)')
# %%
c1.Graph(5,2,
         yaxis_label= 'Gas Conc. Comp. 3 (mol/m$^3$)')
# %%
c1.Graph_P(5)

# %%
# Change the initial values
c1.next_init()
# %%
##### Adsorption Step #####
P_out = 7.9     
P_in = 8      # Ignore this
T_in = 300
y_in = [0.05, 0.55, 0.05, 0.05, 0.30] # H2, N2, CH4, CO, CO2
Cv_in = 5E-4        # m^3/sec/bar
Cv_out = 5E-1          # m^3/sec/bar
u_feed = 0.015            # m/s
Q_in = u_feed*A_cros*epsi  # volumetric flowrate
c1.boundaryC_info(P_out, P_in, T_in, y_in, Cv_in,Cv_out,Q_in, 
                    assigned_v_option = True, 
                    foward_flow_direction=True)
# Mass transfer
#k_MTC = [5E-5, 5E-5, 5E-5, 5E-5, 2E-3]
#c1.mass_trans_info(k_MTC, a_surf, D_disp)
#y_ad, z_ad, t_ad = c1.run_ma_linvel(30,20)
y_ad, z_ad, t_ad = c1.run_ma_linvel(300,25)

# %%
# Graph of 1st component
c1.Graph(30, 4, 
         yaxis_label= 'Gas Conc. Comp. 1 (mol/m$^3$)')
# %%
c1.Graph(30, 9, 
         yaxis_label= 'Gas Conc. Comp. 2 (mol/m$^3$)')
# %%
c1.Graph(20,2,
         yaxis_label= 'Gas Conc. Comp. 3 (mol/m$^3$)')
# %%
c1.Graph_P(50)

# %%
# Changing initial values
c1.next_init()
# %% 
########### Blowdown ##########
P_out = 0.5     
P_in = 7.9      
T_in = 300
y_in = [0.05, 0.55, 0.05, 0.05, 0.30] # H2, N2, CH4, CO, CO2
Cv_in = 0        # m^3/sec/bar
Cv_out = 1E-4          # m^3/sec/bar
u_feed = 0.015            # m/s
Q_in = u_feed*A_cros*epsi  # volumetric flowrate
c1.boundaryC_info(P_out, P_in, T_in, y_in, Cv_in,Cv_out,Q_in, 
                    assigned_v_option = True, 
                    foward_flow_direction=False)
# Mass transfer
k_MTC = [5E-5, 5E-5, 5E-5, 5E-5, 5E-3]
#D_disp = [1E-7, ]*5
#a_surf = 1
c1.mass_trans_info(k_MTC, a_surf, D_disp)
y_bl, z_bl, t_bl = c1.run_ma_linvel(100,50)
# %%
# Graph of 1st component
c1.Graph(10, 4, 
         yaxis_label= 'Gas Conc. CO2 (mol/m$^3$)')
# %%
c1.Graph(10, 9, 
         yaxis_label= 'Solid uptake CO2 (mol/kg)')
# %%
c1.Graph_P(10)
# %%
c1.next_init()
# %%
# %% 
########### Purge step ###########
P_out = 1.4     
P_in = 1.6      
T_in = 300
y_in = [1.00, 0.0, 0.0, 0.0, 0.0] # H2, N2, CH4, CO, CO2
Cv_in = 5E-1        # m^3/sec/bar
Cv_out = 5E-1          # m^3/sec/bar
u_feed = 0.002            # m/s
Q_in = u_feed*A_cros*epsi  # volumetric flowrate
c1.boundaryC_info(P_out, P_in, T_in, y_in, Cv_in,Cv_out,Q_in, 
                    assigned_v_option = True, 
                    foward_flow_direction=False)
k_MTC = [5E-5, 5E-5, 5E-5, 5E-5, 2E-3]
c1.mass_trans_info(k_MTC, a_surf, D_disp)
y_pu, z_pu, t_pu = c1.run_ma_linvel(300,50)
# %%
c1.Graph(30, 4, 
         yaxis_label= 'Gas Conc. Comp. 1 (mol/m$^3$)')
# %%
c1.Graph(30, 9, 
         yaxis_label= 'Gas Conc. Comp. 2 (mol/m$^3$)')
# %%
