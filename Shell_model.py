import matplotlib.pyplot as plt

def orbital_label(l):
    ls = ['s', 'p', 'd', 'f', 'g', 'h']
    return ls[l] if l < len(ls) else f"l={l}"

def plot_shell_model(C=-0.1, D= 0.0225, N_max=5):
    fig, ax = plt.subplots(figsize=(12, 18))
    ax.set_title("Shell Model Energy Levels", fontsize=20)
    ax.set_ylabel("Energy", fontsize=16)
    ax.set_xticks([])
    ax.grid(True, linestyle='--', alpha=0.3)

    for N in range(N_max + 1):
        base_E = N + 1.5
        ax.plot([0, 0.75], [base_E, base_E], color='black')
        ax.text(0.2, base_E + 0.1, f"N={N}", fontsize=12, weight='bold')

        for l in range(N % 2, N + 1, 2):
            n = (N + 2 - l) // 2
            E = base_E + D * l * (l + 1)
            E_up = E + C * l / 2
            E_low = E - C * (l + 1) / 2
            label = orbital_label(l)

            # Main orbital
            ax.plot([1, 1.75], [E, E], color='green')
            ax.text(1.3, E + 0.05, f"{label} (l={l})", fontsize=12)

            # Spin-orbit split
            ax.plot([2.25, 3], [E_up, E_up], color='blue')
            ax.text(2.3, E_up - 0.07, f"{n}{label}{int(2*l+1)}/2", fontsize=10)

            if l > 0:
                ax.plot([2.25, 3], [E_low, E_low], color='red')
                ax.text(3.05, E_low - 0.1, f"{n}{label}{int(2*l-1)}/2", fontsize=10)

            # Connectors
            ax.plot([0.75, 1], [base_E, E], '--', color='gray')
            ax.plot([1.75, 2.25], [E, E_up], '--', color='gray')
            if l > 0:
                ax.plot([1.75, 2.25], [E, E_low], '--', color='gray')

    plt.savefig("Shell model of nucleus.png", dpi=150)
    plt.show()

plot_shell_model()
plt.xticks(fontsize=12)  # X-tick font size
plt.yticks(fontsize=12)  # Y-tick font size

