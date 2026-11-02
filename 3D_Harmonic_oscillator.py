import matplotlib.pyplot as plt

############################################################################### CONSTANTS ############################################################################

hbar_omega= 1.0
max_shell= 6

########################################################################## DEFINING OSCILLATOR ########################################################################

shell_numbers= []
shell_energies= []
shell_degeneracies= []

for N in range(max_shell + 1):

    energy= ( N + 1.5) * hbar_omega
    shell_numbers.append(N)
    shell_energies.append(energy)

    degeneracy= (N + 1) * (N + 2)
    shell_degeneracies.append(degeneracy)
    print(f"\nShell N={N}: energy = {energy:.2f}, states = {degeneracy}")

particles_filled = 0
print("\nSpherical HO shell closures:")

for N, degeneracy in zip(shell_numbers, shell_degeneracies):
    particles_filled += degeneracy
    print(f"\nAfter shell N={N}: {particles_filled}")


orbital_letters = "spdfghi"

for N in shell_numbers:
    print(f"\nOrbitals in shell N={N}:")
    for l in range(N, -1, -2):
        n_r = (N - l) // 2
        orbital = f"{n_r}{orbital_letters[l]}"
        print(f"{orbital}: l={l}")    

######################################################################## PLOTTING #####################################################################################

plt.figure(figsize=(14, 7))

for N, energy, degeneracy in zip(shell_numbers, shell_energies, shell_degeneracies):
    plt.hlines(energy, N - 0.3, N + 0.3, linewidth=3)
    plt.text(N, energy + 0.12, f"g={degeneracy}", ha="center")

plt.xlabel("Oscillator shell number, N")
plt.ylabel(r"Energy, $E/(\hbar\omega)$")
plt.title("3D Harmonic Oscillator Energy Levels")
plt.grid(True)
plt.savefig("3D_HO_energy_levels.png", dpi=300)
plt.show()    