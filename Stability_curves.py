import matplotlib.pyplot as plt

############################################################################ SEMF CONSTANTS ###########################################################################

a_v= 15.75                                        # Volume coefficient
a_s= 17.8                                         # Surface coefficient
a_c= 0.711                                        # Coulomb coefficient
a_a= 23.7                                         # Asymmetry coefficient
a_p= 34.0                                         # Pairing coefficient

############################################################################# DEFINING SEMF ###########################################################################

def binding_energy(Z,N):

    A= Z+N                                        # Calculates mass number

    if Z % 2 == 0 and N % 2 == 0:                 # Even condition

        pairing= a_p / A**(3 / 4)

    elif Z % 2 == 1 and N % 2 == 1:               # Odd condition
        pairing= -a_p / A**(3 / 4)

    else:
        pairing= 0.0

    return (a_v * A - a_s * A**(2 / 3) - a_c * Z * (Z - 1) / A**(1 / 3) - a_a * (A - 2 * Z)**2 / A + pairing)           # SEMF equation


#################################################################### BETA STABILITY PREDICTION ########################################################################

min_A= 16                                        # Sets where SEMF comparison begins
max_A= 300                                       # Sets largest mass number to calculate

beta_N= []                                       # Stores neutron numbers
beta_Z= []                                       # Stores proton numbers

for A in range(min_A, max_A + 1):
      
    best_binding= None
    best_Z= None

    for Z in range(1, A):
        N= A - Z
        BE= binding_energy(Z, N)

        if best_binding is None or BE > best_binding:
            best_binding= BE
            best_Z= Z

    beta_N.append(A - best_Z)
    beta_Z.append(best_Z)


###################################################################### DRIPLINE PREDICTION #############################################################################

max_Z= 118
neutron_drip_N= []
neutron_drip_Z= []
proton_drip_N= []
proton_drip_Z= []

for Z in range(1, max_Z + 1):
    was_bound = False
    for N in range(max(1, min_A - Z), max_A - Z):
        added_neutron_energy = binding_energy(Z, N + 1) - binding_energy(Z, N)                      # Neutron drip line
        if added_neutron_energy > 0:
            was_bound = True
        elif was_bound:
            neutron_drip_N.append(N)
            neutron_drip_Z.append(Z)
            break


for N in range(1, max_A):
    was_bound = False
    for Z in range(max(1, min_A - N), min(max_Z, max_A - N)):
        added_proton_energy = binding_energy(Z + 1, N) - binding_energy(Z, N)                       # Proton drip line
        if added_proton_energy > 0:
            was_bound = True
        elif was_bound:
            proton_drip_Z.append(Z)
            proton_drip_N.append(N)
            break

plot_min_A = 20
neutron_N_plot = []
neutron_Z_plot = []

for N, Z in zip(neutron_drip_N, neutron_drip_Z):
    if N + Z >= plot_min_A:
        neutron_N_plot.append(N)
        neutron_Z_plot.append(Z)        


######################################################################### PLOTTING #####################################################################################

plt.figure(figsize=(18, 10))
plt.plot(beta_N, beta_Z, label="SEMF beta-stability prediction")
plt.plot(neutron_N_plot, neutron_Z_plot, label="SEMF neutron drip line")
plt.plot(proton_drip_N, proton_drip_Z, label="SEMF proton drip line")
plt.plot(range(1, max_Z + 1), range(1, max_Z + 1), "--", label="N = Z")

plt.xlabel("Neutron number, N")
plt.ylabel("Proton number, Z")
plt.title("SEMF-predicted beta stability and drip lines")
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()

plt.savefig("SEMF_stability_and_driplines.png", dpi=200)
plt.show()

