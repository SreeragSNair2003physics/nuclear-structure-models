import numpy as np
import matplotlib.pyplot as plt

############################################################################### CONSTANTS ############################################################################

hbar_omega= 1.0                    # Oscillator energy unit
max_shell= 6                       # Energy shell
states=[]

########################################################################## DEFINING OSCILLATOR ########################################################################

delta_values = np.linspace(-0.2, 0.2, 121)                                                             # Deformation values

for N in range(max_shell + 1):
    for n_z in range(N + 1):
        n_perp = N - n_z
        states.append((N, n_z, n_perp)) 


for N, n_z, n_perp in states:
    energy = hbar_omega * (N + 1.5 - delta_values * (2*n_z - n_perp) / 3)                              # First order anisotropic oscillator energy equation
    plt.plot(delta_values, energy / hbar_omega, color=f"C{N}", label=f"N={N}" if n_z == 0 else None)


######################################################################## PLOTTING #####################################################################################

plt.xlabel("Deformation, δ")
plt.ylabel(r"Energy, $E/(\hbar\omega_0)$")
plt.title("First-order axial deformation splitting of 3D HO levels")
plt.grid(True)
plt.legend(title="Spherical shell")
plt.figtext(0.5, 0.01, "Small-deformation approximation; spin-orbit and Nilsson l² terms not included.", ha="center", fontsize=9)
plt.tight_layout(rect=(0, 0.05, 0.82, 1))
plt.savefig("anisotropic_ho_first_order.png", dpi=300, bbox_inches="tight")
plt.show()