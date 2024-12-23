# Chargement des fichiers
fichiers = [
    '1er_sans_aimant.csv', '2eme_sans_aimant.csv',
    '1er_avec_aimant.csv', '2eme_avec_aimant.csv',
    '3eme_avec_aimant.csv', '4eme_avec_aimant.csv'
]
donnees = charger_donnees(fichiers)

# Renommage des colonnes et mise à l'échelle
renommer_colonnes(donnees)
temps_image(donnees)

# Constante L
L = 17 * 28.246

# Calcul des résultats
resultats = {}
for i, df in enumerate(donnees):
    b, H = calculer_angles(df, L)
    phi_eq, theta_eq = calculer_barycentre(df)
    x_eq, y_eq = calculer_coordonnees_equilibre(b, L, H, theta_eq, phi_eq)
    b, H = nouvelles_coord(df, L, x_eq, y_eq)
    phi_eq, theta_eq = calculer_barycentre(df)
    x_eq, y_eq = calculer_coordonnees_equilibre(b, L, H, theta_eq, phi_eq)
    x_eq, y_eq = x_eq / (28.246 * 100), y_eq / (28.246 * 100)
    resultats[f"df_{i + 1}"] = {
        'b': b, 'H': H, 'phi_eq': phi_eq,
        'theta_eq': theta_eq, 'x_eq': x_eq, 'y_eq': y_eq
    }
    passage_m(df)
    calculer_vitesses(df)

# Affichage des résultats
for cle, res in resultats.items():
    print(f"{cle}: b = {res['b'] / (28.246 * 100):.3f}, H = {res['H'] / (28.246 * 100):.3f}")
    print(f"{cle}: x_eq = {res['x_eq']:.3f}, y_eq = {res['y_eq']:.3f}, phi_eq = {res['phi_eq']:.3f}, theta_eq = {res['theta_eq']:.3f}")

# Calcul des pulsations
resultats_pulsations = {}
for i, df in enumerate(donnees[:2]):
    try:
        T0, w0 = calculer_pulsation(df, colonne_angle='phi', temps='T')
        resultats_pulsations[f"df_{i + 1}"] = {'T0': T0, 'w0': w0}
    except ValueError as e:
        print(f"Erreur pour df_{i + 1}: {e}")

for cle, res in resultats_pulsations.items():
    print(f"{cle}: T0 = {res['T0']:.3f} s, w0 = {res['w0']:.3f} rad/s")

# Calcul de l'énergie cinétique
for df in donnees:
    calculer_energie_cinetique(df)

# Sous-graphes pour les oscillations et l'énergie cinétique
def tracer_graphes_separes(donnees, resultats):
    # Tracer des oscillations (Phi et Theta)
    fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(12, 12))
    for i, df in enumerate(donnees[:2]):
        axes[i, 0].plot(df['T'], df['phi'], label=f"Oscillations Phi - df_{i + 1}")
        axes[i, 0].set_title(f'Oscillations de Phi - df_{i + 1}')
        axes[i, 0].set_xlabel("Temps (s)")
        axes[i, 0].set_ylabel("Phi")
        axes[i, 0].grid()

        axes[i, 1].plot(df['T'], df['theta'], label=f"Oscillations Theta - df_{i + 1}")
        axes[i, 1].set_title(f'Oscillations de Theta - df_{i + 1}')
        axes[i, 1].set_xlabel("Temps (s)")
        axes[i, 1].set_ylabel("Theta")
        axes[i, 1].grid()

    plt.tight_layout()
    plt.show()

    fig, axes = plt.subplots(nrows=4, ncols=2, figsize=(12, 24))
    for i, df in enumerate(donnees[2:]):
        idx = i + 2
        axes[i, 0].plot(df['T'], df['phi'], label=f"Oscillations Phi - df_{idx + 1}")
        axes[i, 0].set_title(f'Oscillations de Phi - df_{idx + 1}')
        axes[i, 0].set_xlabel("Temps (s)")
        axes[i, 0].set_ylabel("Phi")
        axes[i, 0].grid()

        axes[i, 1].plot(df['T'], df['theta'], label=f"Oscillations Theta - df_{idx + 1}")
        axes[i, 1].set_title(f'Oscillations de Theta - df_{idx + 1}')
        axes[i, 1].set_xlabel("Temps (s)")
        axes[i, 1].set_ylabel("Theta")
        axes[i, 1].grid()

    plt.tight_layout()
    plt.show()

    # Portraits de phase
    fig, axes = plt.subplots(2, 1, figsize=(10, 12))
    for i, df in enumerate(donnees[:2]):
        phi = df['phi']
        theta = df['theta']
        axes[i].scatter(phi, theta, label=f"Portrait de phase - df_{i + 1}", s=5)
        pt_eq = [resultats[f"df_{i + 1}"]['phi_eq'], resultats[f"df_{i + 1}"]['theta_eq']]
        axes[i].scatter(pt_eq[0], pt_eq[1], label=f"Point d'équilibre - df_{i + 1}", color='red', s=5)

        axes[i].set_xlabel("Phi")
        axes[i].set_ylabel("Theta")
        axes[i].set_title(f"Portrait de phase - df_{i + 1}")
        axes[i].legend()
        axes[i].grid()

    plt.tight_layout()
    plt.show()

    fig, axes = plt.subplots(4, 1, figsize=(10, 24))
    for i, df in enumerate(donnees[2:]):
        idx = i + 2
        phi = df['phi']
        theta = df['theta']
        axes[i].scatter(phi, theta, label=f"Portrait de phase - df_{idx + 1}", s=5)

        axes[i].set_xlabel("Phi")
        axes[i].set_ylabel("Theta")
        axes[i].set_title(f"Portrait de phase - df_{idx + 1}")
        axes[i].legend()
        axes[i].grid()

    plt.tight_layout()
    plt.show()

    # Ajustement exponentiel pour la dissipation de l'énergie
    fig, axes = plt.subplots(2, 1, figsize=(10, 12))
    for i, df in enumerate(donnees[:2]):
        temps = df['T']
        energie_cinetique = df['E_c']
        axes[i].plot(temps, energie_cinetique, label=f"Énergie cinétique - df_{i + 1}")

        parametres_initiaux = [energie_cinetique.iloc[0], 0.1]
        parametres_opt, _ = curve_fit(modele_exponentiel, temps, energie_cinetique, p0=parametres_initiaux)
        gamma = parametres_opt[1]
        b = 2 * gamma
        print(f"df_{i + 1}: gamma = {gamma:.4f}, coefficient de frottement b = {b:.4f} N·s/m")
        axes[i].plot(temps, modele_exponentiel(temps, *parametres_opt), label=f"Ajustement exponentiel - df_{i + 1}", linestyle="--")

        axes[i].set_xlabel("Temps (s)")
        axes[i].set_ylabel("Énergie cinétique")
        axes[i].set_title(f"Variation de l'énergie cinétique - df_{i + 1}")
        axes[i].legend()
        axes[i].grid()

    plt.tight_layout()
    plt.show()

    fig, axes = plt.subplots(4, 1, figsize=(10, 24))
    for i, df in enumerate(donnees[2:]):
        idx = i + 2
        temps = df['T']
        energie_cinetique = df['E_c']
        axes[i].plot(temps, energie_cinetique, label=f"Énergie cinétique - df_{idx + 1}")

        axes[i].set_xlabel("Temps (s)")
        axes[i].set_ylabel("Énergie cinétique")
        axes[i].set_title(f"Variation de l'énergie cinétique - df_{idx + 1}")
        axes[i].legend()
        axes[i].grid()

    plt.tight_layout()
    plt.show()

tracer_graphes_separes(donnees, resultats)


# ---------------------------- Calcul de l'énergie potentielle ----------------------------

# 12. Calcul de l'énergie potentielle pour chaque DataFrame
for df in donnees[:2]:
    df = calculer_energie_potentielle(df,L)
    
# Calcul de l'énergie totale pour les deux premiers DataFrames
for i, df in enumerate(donnees[:2]):
    df['E_totale'] = df['E_c'] + df['E_p']  # Somme de l'énergie cinétique et potentielle
    
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

for i, ax in enumerate(axes):
    df = donnees[i]
    ax.plot(df['T'], df['E_totale'], label=f"Énergie totale - df_{i+1}", color='red')
    ax.set_xlabel("Temps (T)")
    ax.set_ylabel("Énergie totale (E_totale)")
    ax.set_title(f"Énergie totale vs Temps - df_{i+1}")
    ax.legend()
    ax.grid()
    if i < 2:
        parametres_initiaux = [df['E_totale'].iloc[0], 0.1]
        parametres_opt, _ = curve_fit(modele_exponentiel, df['T'], df['E_totale'], p0=parametres_initiaux)
        gamma = parametres_opt[1]
        b=2*gamma
        print(f"df_{i+1}: gamma = {gamma:.4f}, coefficient de frottement b = {b:.4f} N·s/m")
        axes[i].plot(df['T'], modele_exponentiel(df['T'], *parametres_opt), label=f"Ajustement exponentiel - df_{i+1}", linestyle="--")

plt.tight_layout()
plt.show()
