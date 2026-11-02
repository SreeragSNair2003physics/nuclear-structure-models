import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import linear_sum_assignment

############################################################################## CONSTANTS ################################################################################

hbaromegazero= 1.0
N_max= 7
Kappa= [0.12, 0.12, 0.105, 0.09, 0.09, 0.07, 0.065, 0.06]
mu= [0.0, 0.0, 0.0, 0.24, 0.18, 0.30, 0.30, 0.30]
delta_vals = np.linspace(-0.3, 0.3, 121)

########################################################################## BASIS GENERATION  ###########################################################################

def basisgenerator(Nmax):

    basis= []
    
    for N in range(0, Nmax+1):

        L_min= 0 if N % 2 == 0 else 1

        for L in range(N,L_min-1,-2):

            for Lambda in range(-L, L+1):

                J= L + 0.5 

                for Omega in np.arange(-J, J+1):

                    Sigma= Omega - Lambda

                    if abs(abs(Sigma) - 0.5)<= 1e-8: 

                        basis.append((N, L, Lambda, Sigma)) 

    return basis 




######################################################################### DEFINING HAMILTONIAN #########################################################################

def Hamiltonian(basis, kappa, delta) :
   
    basis_size= len(basis)
    H=np.zeros([basis_size, basis_size])
    for i, state_i in enumerate(basis):
        for j, state_j in enumerate(basis):
            N_i,L_i,Lambda_i,Sigma_i = state_i 
            N_j,L_j,Lambda_j,Sigma_j = state_j 



            H_ij= 0.0

            if (N_i==N_j) and (L_i==L_j) and (Lambda_i==Lambda_j) and abs(Sigma_i-Sigma_j)<1e-8:

                H_ij += N_i + 1.5
                H_ij += - Kappa[N_i]* mu[N_i]*(L_i*(L_i+1) - N_i*(N_i+3)/2)


            if (N_i==N_j) and (L_i==L_j):                                                                     

                if (Lambda_j== Lambda_i+1) and abs(Sigma_j - (Sigma_i-1)) < 1e-8:
                    ldots= (1/2) * np.sqrt((L_i - Lambda_i) * (L_i + Lambda_i + 1))
                    
                
                elif(Lambda_j== Lambda_i-1) and abs(Sigma_j - (Sigma_i+1)) < 1e-8:
                    ldots= (1/2) * np.sqrt((L_i + Lambda_i) * (L_i - Lambda_i + 1))
                
                elif(Lambda_j== Lambda_i) and abs(Sigma_j - (Sigma_i)) < 1e-8:
                    ldots= Lambda_i*Sigma_i
                
                else: 
                    ldots=0.0

                H_ij += -2*Kappa[N_i] * ldots




            r2=0.0

            if (N_i==N_j) and (Lambda_i==Lambda_j) and abs(Sigma_i-Sigma_j)<1e-8:                             
                if (L_j ==L_i-2):
                    r2= np.sqrt((N_i - L_i + 2) * (N_i + L_i + 1))
                
                elif(L_j==L_i):
                    r2= (N_i + (3/2))
                
                elif(L_j == L_i+2):
                    r2= np.sqrt((N_i - L_i) * (N_i + L_i + 3))



            Y20 = 0.0
            if (N_i==N_j) and abs(Sigma_i-Sigma_j)<1e-8 and Lambda_i==Lambda_j:
                l, lp, m = L_i, L_j, Lambda_i
                if lp == l:
                    Y20 = np.sqrt(5/(4*np.pi))*(l*(l+1)-3*m**2)/((2*l-1)*(2*l+3))
                elif lp == l+2:
                    Y20 = np.sqrt(5/(16*np.pi))*3*np.sqrt(((l+1)**2-m**2)*((l+2)**2-m**2)) \
                        /((2*l+3)*np.sqrt((2*l+1)*(2*l+5)))
                elif lp == l-2:
                    l = lp; 
                    Y20 = np.sqrt(5/(16*np.pi))*3*np.sqrt(((l+1)**2-m**2)*((l+2)**2-m**2)) \
                        /((2*l+3)*np.sqrt((2*l+1)*(2*l+5)))



            H_delta= -1 * delta * hbaromegazero * (4/3) * np.sqrt((np.pi)/5) * r2 * Y20   
            
            H_ij+= H_delta  


            H[i,j]= H_ij  

    return H
            


def evalue_calculation(kappa, delta, basis):
    H = Hamiltonian(basis, kappa, delta)
    keys = [(N, int(round(2*(Lam+Sig)))) for (N, L, Lam, Sig) in basis]
    out = {}
    for key in sorted(set(keys)):
        if key[1] <= 0:          
            continue
        idx = [i for i, k in enumerate(keys) if k == key]
        out[key] = np.linalg.eigvalsh(H[np.ix_(idx, idx)])
    return out


def block_eigensystems(kappa, delta, basis):
    H = Hamiltonian(basis, kappa, delta)
    keys = [(N, int(round(2 * (Lam + Sig))) )
            for N, L, Lam, Sig in basis]
    systems = {}

    for key in sorted(set(keys)):
        if key[1] <= 0:
            continue
        idx = [i for i, k in enumerate(keys) if k == key]
        block = H[np.ix_(idx, idx)]
        energies, vectors = np.linalg.eigh(block)
        systems[key] = (idx, energies, vectors)

    return systems

############################################################################### PLOTTING  ##############################################################################

basis = basisgenerator(N_max)

curves = {}
for delta in delta_vals:
    for key, ev in evalue_calculation(Kappa, delta, basis).items():
        curves.setdefault(key, []).append(ev)

particle_numbers = [8, 20, 28, 50, 82, 126, 184]
gap_curves = {p: [] for p in particle_numbers}

for d in range(len(delta_vals)):
    levels = []
    for energy_list in curves.values():
        levels.extend(np.repeat(energy_list[d], 2))
    levels = np.sort(levels)

    for p in particle_numbers:
        gap_curves[p].append(levels[p] - levels[p - 1])        


spherical_eigenvalues = []
for values in evalue_calculation(Kappa, 0.0, basis).values():
    spherical_eigenvalues.extend(values)

spherical_eigenvalues.sort()
spherical_levels = []
for energy in spherical_eigenvalues:
    if not spherical_levels or abs(energy - spherical_levels[-1][0]) > 1e-7:
        spherical_levels.append([energy, 2])
    else:
        spherical_levels[-1][1] += 2



magic_targets = {8, 20, 28, 50}
magic_positions = {}
filled = 0

for i, (energy, degeneracy) in enumerate(spherical_levels[:-1]):
    filled += degeneracy
    if filled in magic_targets:
        next_energy = spherical_levels[i + 1][0]
        magic_positions[filled] = (energy + next_energy) / 2

print("Magic-number circle positions:", magic_positions)



orbital_letters = "spdfghij"
spherical_labels = []

for N in range(N_max + 1):
    for L in range(N, -1, -2):
        nr = (N - L) // 2 + 1
        j_values = [2 * L + 1] if L == 0 else [2 * L + 1, 2 * L - 1]

        for j2 in j_values:
            ls_value = L / 2 if j2 == 2 * L + 1 else -(L + 1) / 2
            l2_value = L * (L + 1) - N * (N + 3) / 2
            energy = N + 1.5 - 2 * Kappa[N] * ls_value - Kappa[N] * mu[N] * l2_value
            label = f"{nr}{orbital_letters[L]}$_{{{j2}/2}}$"
            spherical_labels.append((energy, label))



H_spherical = Hamiltonian(basis, Kappa, 0.0)
H_at_unit_deformation = Hamiltonian(basis, Kappa, 1.0)
deformation_operator = H_at_unit_deformation - H_spherical

basis_keys = [(N, int(round(2 * (Lam + Sig))))
              for N, L, Lam, Sig in basis]
asymptotic_labels = {}
asymptotic_vectors = {}

for key in sorted(set(basis_keys)):
    if key[1] <= 0:
        continue
    N, om2 = key
    idx = [i for i, k in enumerate(basis_keys) if k == key]

    block_states = [basis[i] for i in idx]
    deformation_block = deformation_operator[np.ix_(idx, idx)]
    labels_for_block = []
    vectors_for_block = []
    sectors = sorted(set((state[2], state[3]) for state in block_states))

    for Lambda, Sigma in sectors:
        sector_idx = [p for p, state in enumerate(block_states)
                      if state[2] == Lambda and state[3] == Sigma]
        sector_matrix = deformation_block[np.ix_(sector_idx, sector_idx)]
        slopes, sector_vectors = np.linalg.eigh(sector_matrix)
        for q, slope in enumerate(slopes):
            nz = int(round(N / 3 - slope))
            label = f"{om2}/2[{N}{nz}{Lambda}]"
            ref_vector = np.zeros(len(block_states))
            ref_vector[sector_idx] = sector_vectors[:, q]
            labels_for_block.append((slope, label))
            vectors_for_block.append(ref_vector)

    asymptotic_labels[key] = labels_for_block
    asymptotic_vectors[key] = vectors_for_block


# Labels are read off at a deformation far outside the plot (delta = 3), where every level is a
# pure |N nz Lambda Sigma> state, then attached to the curve with the same index in its block.
H_asym = Hamiltonian(basis, Kappa, 3.0)
curve_labels = {}
for key in asymptotic_labels:
    idx = [i for i, k in enumerate(basis_keys) if k == key]
    _, vectors_asym = np.linalg.eigh(H_asym[np.ix_(idx, idx)])
    reference_vectors = np.column_stack(asymptotic_vectors[key])
    overlaps = np.abs(reference_vectors.T @ vectors_asym) ** 2
    label_idx, level_idx = linear_sum_assignment(-overlaps)
    names = [lab for _, lab in asymptotic_labels[key]]
    curve_labels[key] = [None] * len(idx)
    for a, b in zip(label_idx, level_idx):
        curve_labels[key][b] = names[a]

fig, ax = plt.subplots(figsize=(8.5, 10.5), dpi=100)

ax.axvline(0.0, color="k", lw=0.8, zorder=1)

for (N, om2), energy_list in curves.items():
    energies = np.array(energy_list)
    style = "-" if N % 2 == 0 else (0, (2.5, 1.2))

    for level in range(energies.shape[1]):
        ax.plot(
            delta_vals,
            energies[:, level],
            linestyle=style,
            color="k",
            lw=0.9,
            zorder=2
        )


for number, energy in magic_positions.items():
    ax.text(
        0.0, energy, str(number),
        ha="center", va="center", fontsize=8,
        bbox=dict(boxstyle="circle,pad=0.25",
                  facecolor="white", edgecolor="black")
    )




visible = [(e, l) for e, l in sorted(spherical_labels) if 2.0 <= e <= 5.4]
for i, (energy, label) in enumerate(visible):
    side = -1 if i % 2 == 0 else 1            # alternate left / right of the axis
    ax.text(side * 0.01, energy + (0.08 if side < 0 else -0.07), label,
            ha="right" if side < 0 else "left", va="center",
            fontsize=10, fontweight="bold", zorder=5,
            bbox=dict(boxstyle="square,pad=0.08", fc="white", ec="none"))

ax.set_xlim(-0.3, 0.3)
ax.set_ylim(2.0, 5.4)
ax.set_xticks(np.arange(-0.3, 0.31, 0.1))
ax.set_yticks(np.arange(2.0, 5.01, 0.5))

ax.set_xlabel(r"$\varepsilon_2$", fontsize=13)
ax.set_ylabel(r"$E_{\rm s.p.}\ (\hbar\omega_0)$", fontsize=13)
plt.tight_layout()


fig.canvas.draw()
renderer = fig.canvas.get_renderer()
ax_box = ax.get_window_extent(renderer)
E_all = {key: np.array(v) for key, v in curves.items()}
curve_list = [(key, k) for key, E in E_all.items() for k in range(E.shape[1])]
pts, ids = [], []
for cid, (key, k) in enumerate(curve_list):
    c = E_all[key][:, k]
    vis = (c > 2.0) & (c < 5.4)
    pts.append(ax.transData.transform(np.column_stack([delta_vals[vis], c[vis]])))
    ids.append(np.full(vis.sum(), cid))
pts, ids = np.vstack(pts), np.concatenate(ids)
placed = []

def try_place(x, curve, cid, text, allow_line_overlap):
    y = float(np.interp(x, delta_vals, curve))
    p0 = ax.transData.transform((x - 0.005, np.interp(x - 0.005, delta_vals, curve)))
    p1 = ax.transData.transform((x + 0.005, np.interp(x + 0.005, delta_vals, curve)))
    ang = np.degrees(np.arctan2(p1[1] - p0[1], p1[0] - p0[0]))
    t = ax.text(x, y, text, fontsize=7, rotation=ang, rotation_mode="anchor",
                ha="center", va="bottom", zorder=3)
    bb = t.get_window_extent(renderer).expanded(1.02, 1.10)
    ok = (bb.x0 > ax_box.x0 and bb.x1 < ax_box.x1 and bb.y0 > ax_box.y0 and bb.y1 < ax_box.y1
          and not any(bb.overlaps(o) for o in placed))
    if ok and not allow_line_overlap:
        hit = ((pts[:, 0] > bb.x0) & (pts[:, 0] < bb.x1) & (pts[:, 1] > bb.y0)
               & (pts[:, 1] < bb.y1) & (ids != cid))
        ok = not hit.any()
    if ok:
        placed.append(bb)
        return True
    t.remove()
    return False

candidates = sorted(np.arange(0.05, 0.295, 0.005), key=lambda v: abs(v - 0.17))
for cid, (key, k) in enumerate(curve_list):
    curve, text = E_all[key][:, k], curve_labels[key][k]
    for side in (-1, +1):
        for allow in (False, True):           # first a clean spot, then accept line overlap
            done = False
            for ax_x in candidates:
                x = side * ax_x
                y = float(np.interp(x, delta_vals, curve))
                if 2.04 < y < 5.36 and try_place(x, curve, cid, text, allow):
                    done = True
                    break
            if done:
                break

plt.savefig("nilsson_diagram.png", dpi=150)
plt.show()

fig, ax = plt.subplots(figsize=(9, 6))

for p, gaps in gap_curves.items():
    ax.plot(delta_vals, gaps, label=f"p = {p}")

ax.set_xlabel(r"Deformation, $\varepsilon_2$")
ax.set_ylabel(r"Shell gap, $\Delta E/(\hbar\omega_0)$")
ax.set_title("Shell-gap evolution with deformation")
ax.grid(True, alpha=0.3)
ax.legend(title="Particles filled", ncol=2)
plt.tight_layout()
plt.savefig("Deformation Dependence of Neutron Shell Gaps", dpi=150)
plt.show()