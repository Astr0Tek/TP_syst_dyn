import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import find_peaks
from scipy.optimize import curve_fit

# ---------------------------- Fonctions de chargement et traitement des données ----------------------------

# Fonction pour charger les données
def charger_donnees(fichiers):
    return [pd.read_csv(f, skiprows=2) for f in fichiers]

# Fonction pour renommer les colonnes
def renommer_colonnes(donnees):
    for df in donnees:
        df.columns = ['Etiquette', 'T', 'X', 'Y']

#Fonction pour mettre le temps en s
def temps_image(donnees):
    for df in donnees:
        df['T']=df['T']/100


# ---------------------------- Fonctions de calcul dynamique ----------------------------

# Fonction pour calculer les angles et les paramètres dynamiques
def calculer_angles(df, L):
    df['r'] = np.sqrt(df['X']**2 + df['Y']**2)
    masque_petits_angles = np.abs(df['Y']) / L < 0.27
    indices_valides = np.where(masque_petits_angles)[0]

    if len(indices_valides) < 2:
        raise ValueError("Pas assez de points respectant l'approximation des petits angles.")

    idx1, idx2 = indices_valides[0], indices_valides[-1]
    r1, r2 = df['r'].iloc[idx1], df['r'].iloc[idx2]
    y1, y2 = df['Y'].iloc[idx1], df['Y'].iloc[idx2]

    A1 = (1 - np.sqrt(1 - y1**2 / L**2)) * L**2
    A2 = (1 - np.sqrt(1 - y2**2 / L**2)) * L**2
    B = (r1 * y2) / (r2 * y1)
    H = abs((A1 * B - A2) / (1 - B))
    b = (r2 * (H + A1)) / y2

    df['A'] = df['r'] / b
    facteur_B = H / L + 1
    df['theta'] = np.arcsin(df['A'] * (facteur_B + np.sqrt(2 * facteur_B**2 - df['A']**2 * facteur_B**2 - df['A']**2 - 1)) / (df['A']**2 + 1))
    df['phi'] = np.arctan2(df['Y'], df['X'])

    return b, H

# Fonction pour calculer les vitesses
def calculer_vitesses(df):
    df['v_theta'] = df['r'] * np.gradient(df['theta'], df['T'])
    df['v_phi'] = df['r'] * np.sin(df['theta']) * np.gradient(df['phi'], df['T'])

# Fonction pour calculer l'énergie cinétique
def calculer_energie_cinetique(df, masse=1):
    df['E_c'] = 0.5 * masse * (df['v_theta']**2 + df['v_phi']**2)

# Fonction pour calculer l'énergie potentielle gravitationnelle
def calculer_energie_potentielle(df, masse=1, g=9.81):
    df['E_p'] = masse * g * df['Y']
    return df

# ---------------------------- Fonctions d'analyse des oscillations ----------------------------

# Fonction pour calculer la période T0 et la pulsation w0
def calculer_pulsation(df, colonne_angle='phi', temps='T'):
    angles = df[colonne_angle]
    temps = df[temps]
    
    indices_pics, _ = find_peaks(angles)
    
    if len(indices_pics) < 2:
        raise ValueError("Pas assez de pics pour déterminer la période.")
    
    temps_pics = temps.iloc[indices_pics].values
    differences_temps = np.diff(temps_pics)
    T0 = np.mean(differences_temps)
    w0 = 2 * np.pi / T0
    
    return T0, w0

# ---------------------------- Fonctions de calcul du barycentre et des coordonnées d'équilibre ----------------------------

# Fonction pour calculer le barycentre
def calculer_barycentre(df):
    poids = np.abs(df['v_phi'])
    barycentre_phi = np.average(df['phi'], weights=poids)
    barycentre_theta = np.average(df['theta'], weights=poids)
    return barycentre_phi, barycentre_theta

# Fonction pour calculer les coordonnées d'équilibre
def calculer_coordonnees_equilibre(b, L, H, theta_eq, phi_eq):
    A = b * (L * np.sin(theta_eq)) / (H + (1 - np.cos(theta_eq)) * L)
    x_eq = A / np.sqrt(1 + np.tan(phi_eq)**2)
    y_eq = x_eq * np.tan(phi_eq)
    return x_eq, y_eq

# Fonction pour calculer la distance entre les points d'équilibre
def calculer_distance(resultats):
    coordonnees = [(res['x_eq'], res['y_eq']) for res in resultats.values()]
    return np.sqrt((coordonnees[1][0] - coordonnees[0][0])**2 + (coordonnees[1][1] - coordonnees[0][1])**2)

# ---------------------------- Modèle exponentiel et ajustement ----------------------------

# Fonction exponentielle pour ajuster l'énergie au temps
def modele_exponentiel(t, E0, gamma):
    return E0 * np.exp(-gamma * t)

# ---------------------------- Chargement des données et calculs ----------------------------

# 1. Chargement des données
fichiers = ['1er_sans_aimant.csv', '2eme_sans_aimant.csv', 
            '1er_avec_aimant.csv', '2eme_avec_aimant.csv', 
            '3eme_avec_aimant.csv', '4eme_avec_aimant.csv']
donnees = charger_donnees(fichiers)

# 2. Renommage des colonnes + mise à la bonne échelle de temps
renommer_colonnes(donnees)
temps_image(donnees)

# 3. Constante L
L = 17*28.246

# 4. Calcul des angles et des vitesses
for df in donnees:
    try:
        b, H = calculer_angles(df, L)
        calculer_vitesses(df)
    except ValueError as e:
        print(f"Erreur pour df_{donnees.index(df)+1}: {e}")

# ---------------------------- Analyse du barycentre et des coordonnées d'équilibre ----------------------------

# 5. Calcul du barycentre et des coordonnées d'équilibre pour les deux premiers DataFrames
resultats = {}
for i, df in enumerate(donnees[:2]):
    try:
        phi_eq, theta_eq = calculer_barycentre(df)
        x_eq, y_eq = calculer_coordonnees_equilibre(b, L, H, theta_eq, phi_eq)
        resultats[f"df_{i+1}"] = {'b': b, 'H': H, 'phi_eq': phi_eq, 'theta_eq': theta_eq, 'x_eq': x_eq, 'y_eq': y_eq}
    except ValueError as e:
        print(f"Erreur pour df_{i+1}: {e}")

# 6. Calcul de la distance entre les points d'équilibre
distance = calculer_distance(resultats)

# Affichage des résultats des coordonnées d'équilibre
for cle, res in resultats.items():
    print(f"{cle}: x_eq = {res['x_eq']:.3f}, y_eq = {res['y_eq']:.3f}, phi_eq = {res['phi_eq']:.3f}, theta_eq = {res['theta_eq']:.3f}")
print(f"Distance entre les points d'équilibre : {distance:.3f}")

# ---------------------------- Analyse de la pulsation ----------------------------

# 7. Calcul de la période T0 et de la pulsation w0 pour les deux premiers DataFrames
resultats_pulsations = {}
for i, df in enumerate(donnees[:2]):
    try:
        T0, w0 = calculer_pulsation(df, colonne_angle='phi', temps='T')
        resultats_pulsations[f"df_{i+1}"] = {'T0': T0, 'w0': w0}
    except ValueError as e:
        print(f"Erreur pour df_{i+1}: {e}")

# Affichage des résultats de la pulsation
for cle, res in resultats_pulsations.items():
    print(f"{cle}: T0 = {res['T0']:.3f} s, w0 = {res['w0']:.3f} rad/s")

# ---------------------------- Calcul de l'énergie et affichage des résultats ----------------------------

# 8. Calcul de l'énergie cinétique pour tous les DataFrames
for df in donnees:
    calculer_energie_cinetique(df)

# 9. Création de sous-graphes pour les oscillations et l'énergie cinétique
fig, axes = plt.subplots(nrows=len(donnees), ncols=2, figsize=(12, 6 * len(donnees)))
for i, df in enumerate(donnees):
    axes[i, 0].plot(df['T'], df['phi'], label=f"Oscillations Phi - df_{i+1}")
    axes[i, 0].set_title(f'Oscillations de Phi - df_{i+1}')
    axes[i, 0].set_xlabel("Temps (s)")
    axes[i, 0].set_ylabel("Phi")
    axes[i, 0].grid()
    
    axes[i, 1].plot(df['T'], df['theta'], label=f"Oscillations Theta - df_{i+1}")
    axes[i, 1].set_title(f'Oscillations de Theta - df_{i+1}')
    axes[i, 1].set_xlabel("Temps (s)")
    axes[i, 1].set_ylabel("Theta")
    axes[i, 1].grid()

plt.tight_layout()
plt.show()

# 10. Portrait de phase
fig, axes = plt.subplots(len(donnees), 1, figsize=(10, 6 * len(donnees)))
for i, df in enumerate(donnees):
    phi = df['phi']
    theta = df['theta']
    axes[i].scatter(phi, theta, label=f"Portrait de phase - df_{i+1}")
    
    if i < 2:
        pt_eq = [resultats[f"df_{i+1}"]['phi_eq'],resultats[f"df_{i+1}"]['theta_eq']]
        axes[i].scatter(pt_eq[0], pt_eq[1], label=f"Point d'équilibre - df_{i+1}", color='red')
    
    axes[i].set_xlabel("Phi")
    axes[i].set_ylabel("Theta")
    axes[i].set_title(f"Portrait de phase - df_{i+1}")
    axes[i].legend()
    axes[i].grid()

plt.tight_layout()
plt.show()

# ---------------------------- Ajustement exponentiel pour la dissipation de l'énergie ----------------------------

# 11. Ajustement exponentiel pour la dissipation de l'énergie
fig, axes = plt.subplots(len(donnees), 1, figsize=(10, 6 * len(donnees)))
for i, df in enumerate(donnees):
    temps = df['T']
    energie_cinetique = df['E_c']
    axes[i].plot(temps, energie_cinetique, label=f"Énergie cinétique - df_{i+1}")
    
    if i < 2:
        parametres_initiaux = [energie_cinetique.iloc[0], 0.1]
        parametres_opt, _ = curve_fit(modele_exponentiel, temps, energie_cinetique, p0=parametres_initiaux)
        gamma = parametres_opt[1]
        print(f"Taux de dissipation de l'énergie (gamma) pour df_{i+1}: {gamma:.4f} s^-1")
        axes[i].plot(temps, modele_exponentiel(temps, *parametres_opt), label=f"Ajustement exponentiel - df_{i+1}", linestyle="--")
    
    axes[i].set_xlabel("Temps (s)")
    axes[i].set_ylabel("Énergie cinétique")
    axes[i].set_title(f"Variation de l'énergie cinétique - df_{i+1}")
    axes[i].legend()
    axes[i].grid()

plt.tight_layout()
plt.show()

# ---------------------------- Calcul de l'énergie potentielle ----------------------------

# 12. Calcul de l'énergie potentielle pour chaque DataFrame
for df in donnees[:2]:
    df = calculer_energie_potentielle(df)
    
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

plt.tight_layout()
plt.show()
