import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import find_peaks
from scipy.optimize import curve_fit

# Fonction pour charger les données
def charger_donnees(fichiers):
    return [pd.read_csv(f, skiprows=2) for f in fichiers]

# Fonction pour renommer les colonnes
def renommer_colonnes(donnees):
    for df in donnees:
        df.columns = ['Etiquette', 'T', 'X', 'Y']

# Fonction pour mettre le temps en s
def temps_image(donnees):
    for df in donnees:
        df['T'] = df['T']
        
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
L = 0.17

x_eq,y_eq=np.average(donnees[1]['X']),np.average(donnees[1]['Y'])

for df in donnees:
    df['x']=(df['X']-x_eq)
    df['y']=(df['Y']-y_eq)
    df['r']=np.sqrt(df['x']**2 + df['y']**2)
    
#x_eq,y_eq=x_eq/(100*28.246),y_eq/(100*28.246)
    
# Fonction pour calculer les vitesses
def vitesses(df):
    df['v_x'] = np.gradient(df['x'], df['T']) *100
    df['v_y'] = np.gradient(df['y'], df['T'])  *100
    df['v_theta'] = df['r'] * np.gradient(df['theta'], df['T'])*100
    df['v_phi'] = df['r'] * np.sin(df['theta']) * np.gradient(df['phi'], df['T'])*100

def e_c(df, masse=1):
    df['E_c'] = 0.5 * masse * (df['v_x']**2 +df['v_y']**2 +((np.gradient(df['z'],df['T'])) *100)**2)

def modele_exponentiel(t, E0, gamma):
    return E0 * np.exp(-gamma * t)

def e_t(df,masse=1):
    df['E_t']= df['E_c'] +masse*9.81*df['z']

b=220.23
    
H=0.083975
    
plt.show()

def angles(df, L,b,H):
    df['A'] = df['r']/b
    B = H / L + 1
    df['phi'] = np.arctan(df['y']/(df['x']))
    df['theta'] =np.sign(df['x'])* np.arcsin((df['A']*(B-np.sqrt(df['A']**2 +1 -df['A']**2 *B**2)))/(df['A']**2 +1))

    
for df in donnees:
    angles(df,L,b,H)


for df in donnees:
    df['x']=L*np.sin(df['theta'])*np.cos(df['phi'])
    df['y']=L*np.sin(df['theta'])*np.sin(df['phi'])
    df['z']=L*(1-np.cos(df['theta']))
    
for df in donnees:
    vitesses(df)
    e_c(df)
    e_t(df)

def tracer_graphes_separes(donnees):
    
    # Tracer oscillations de Phi
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    for i, df in enumerate(donnees[:2]):
        axes[i].plot(df['T'] / 100, df['phi'], label=f"Oscillations Phi - df_{i + 1}")
        axes[i].set_title(f'Oscillations de Phi - df_{i + 1}')
        axes[i].set_xlabel("Temps (s)")
        axes[i].set_ylabel("Phi (rad)")
        axes[i].grid()
        axes[i].legend()
    plt.tight_layout()
    plt.savefig(f'Oscillations de Phi - df_{i + 1}.jpeg')
    plt.show()

    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    axes = axes.flatten()
    for i, df in enumerate(donnees[2:]):
        idx = i + 2
        axes[i].plot(df['T'] / 100, df['phi'], label=f"Oscillations Phi - df_{idx + 1}")
        axes[i].set_title(f'Oscillations de Phi - df_{idx + 1}')
        axes[i].set_xlabel("Temps (s)")
        axes[i].set_ylabel("Phi (rad)")
        axes[i].grid()
        axes[i].legend()
    plt.tight_layout()
    plt.savefig(f'Oscillations de Phi - df_{idx + 1}.jpeg')
    plt.show()

    # Tracer oscillations de Theta
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    for i, df in enumerate(donnees[:2]):
        axes[i].plot(df['T'] / 100, df['theta'], label=f"Oscillations Theta - df_{i + 1}")
        axes[i].set_title(f'Oscillations de Theta - df_{i + 1}')
        axes[i].set_xlabel("Temps (s)")
        axes[i].set_ylabel("Theta (rad)")
        axes[i].grid()
        axes[i].legend()
    plt.tight_layout()
    plt.savefig(f'Oscillations de Theta - df_{i + 1}.jpeg')
    plt.show()

    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    axes = axes.flatten()
    for i, df in enumerate(donnees[2:]):
        idx = i + 2
        axes[i].plot(df['T'] / 100, df['theta'], label=f"Oscillations Theta - df_{idx + 1}")
        axes[i].set_title(f'Oscillations de Theta - df_{idx + 1}')
        axes[i].set_xlabel("Temps (s)")
        axes[i].set_ylabel("Theta (rad)")
        axes[i].grid()
        axes[i].legend()
    plt.tight_layout()
    plt.savefig(f'Oscillations de Theta - df_{idx + 1}.jpeg')
    plt.show()

    # Portraits de phase
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    for i, df in enumerate(donnees[:2]):
        phi = df['phi']
        theta = df['theta']
        axes[i].scatter(phi, theta, label=f"Portrait de phase - df_{i + 1}", s=5)
        axes[i].set_title(f"Portrait de phase - df_{i + 1}")
        axes[i].set_xlabel("Phi (rad)")
        axes[i].set_ylabel("Theta (rad)")
        axes[i].legend()
        axes[i].grid()
    plt.tight_layout()
    plt.savefig(f"Portrait de phase - df_{i + 1}.jpeg")
    plt.show()

    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    axes = axes.flatten()
    for i, df in enumerate(donnees[2:]):
        idx = i + 2
        phi = df['phi']
        theta = df['theta']
        axes[i].scatter(phi, theta, label=f"Portrait de phase - df_{idx + 1}", s=5)
        axes[i].set_title(f"Portrait de phase - df_{idx + 1}")
        axes[i].set_xlabel("Phi (rad)")
        axes[i].set_ylabel("Theta (rad)")
        axes[i].legend()
        axes[i].grid()
    plt.tight_layout()
    plt.savefig(f"Portrait de phase - df_{idx + 1}.jpeg")
    plt.show()
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    for i, df in enumerate(donnees[:2]):
        x = df['x']
        y = df['y']
        axes[i].scatter(x, y, label=f"Portrait de phase - df_{i + 1}", s=5)
        axes[i].set_title(f"Portrait de phase (x,y) - df_{i + 1}")
        axes[i].set_xlabel("x (m)")
        axes[i].set_ylabel("y (m)")
        axes[i].legend()
        axes[i].grid()
    plt.tight_layout()
    plt.savefig(f"Portrait de phase (x,y) - df_{i + 1}.jpeg")
    plt.show()

    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    axes = axes.flatten()
    for i, df in enumerate(donnees[2:]):
        idx = i + 2
        x = df['x']/(100*28.246)
        y = df['y']/(100*28.246)
        axes[i].scatter(x, y, label=f"Portrait de phase - df_{idx + 1}", s=5)
        axes[i].set_title(f"Portrait de phase (x,y) - df_{idx + 1}")
        axes[i].set_xlabel("x (m)")
        axes[i].set_ylabel("y (m)")
        axes[i].legend()
        axes[i].grid()
    plt.tight_layout()
    plt.savefig(f"Portrait de phase (x,y) - df_{idx + 1}.jpeg")
    plt.show()


    # Variation de l'énergie cinétique
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    for i, df in enumerate(donnees[:2]):
        temps = df['T'] / 100
        energie_cinetique = df['E_c']
        axes[i].plot(temps, energie_cinetique, label=f"Énergie cinétique - df_{i + 1}")
        axes[i].set_title(f"Variation de l'énergie cinétique - df_{i + 1}")
        axes[i].set_xlabel("Temps (s)")
        axes[i].set_ylabel("Énergie cinétique (J/kg)")
        axes[i].legend()
        axes[i].grid()
    plt.tight_layout()
    plt.savefig(f"Variation de l'énergie cinétique - df_{i + 1}.jpeg")
    plt.show()

    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    axes = axes.flatten()
    for i, df in enumerate(donnees[2:]):
        idx = i + 2
        temps = df['T'] / 100
        energie_cinetique = df['E_c']
        axes[i].plot(temps, energie_cinetique, label=f"Énergie cinétique - df_{idx + 1}")
        axes[i].set_title(f"Variation de l'énergie cinétique - df_{idx + 1}")
        axes[i].set_xlabel("Temps (s)")
        axes[i].set_ylabel("Énergie cinétique (J/kg)")
        axes[i].legend()
        axes[i].grid()
    plt.tight_layout()
    plt.savefig(f"Variation de l'énergie cinétique - df_{idx + 1}.jpeg")
    plt.show()
    
    # Energie mécanique
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    for i, df in enumerate(donnees[:2]):
        temps = df['T'] / 100
        energie_cinetique = df['E_t']
        axes[i].plot(temps, energie_cinetique, label=f"Énergie mécanique - df_{i + 1}")
        parametres_initiaux = [df['E_t'].iloc[0], 0.1]
        parametres_opt, _ = curve_fit(modele_exponentiel, df['T']/100, df['E_t'], p0=parametres_initiaux)
        gamma = parametres_opt[1]
        b=2*gamma
        print(f"df_{i+1}: gamma = {gamma:.4f}, coefficient de frottement b = {b:.4f} N·s/m")
        axes[i].plot(df['T']/100, modele_exponentiel(df['T']/100, *parametres_opt), label=f"Ajustement exponentiel - df_{i+1}", linestyle="--")
        axes[i].set_title(f"Energie mécanique - df_{i + 1}")
        axes[i].set_xlabel("Temps (s)")
        axes[i].set_ylabel("Énergie mécanique (J/kg)")
        axes[i].legend()
        axes[i].grid()
    plt.tight_layout()
    plt.savefig(f"Energie mécanique - df_{i + 1}.jpeg")
    plt.show()

    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    axes = axes.flatten()
    for i, df in enumerate(donnees[2:]):
        idx = i + 2
        temps = df['T'] / 100
        energie_cinetique = df['E_t']
        axes[i].plot(temps, energie_cinetique, label=f"Énergie mécanique - df_{idx + 1}")
        axes[i].set_title(f"Energie mécanique - df_{idx + 1}")
        axes[i].set_xlabel("Temps (s)")
        axes[i].set_ylabel("Énergie mécanique (J/kg)")
        axes[i].legend()
        axes[i].grid()
    plt.tight_layout()
    plt.savefig(f"Energie mécanique - df_{idx + 1}.jpeg")
    plt.show()


tracer_graphes_separes(donnees)

def calculer_pulsation(df, colonne_angle='theta', temps='T'):
    angles = df[colonne_angle]
    temps = df[temps]/100

    indices_pics, _ = find_peaks(angles)

    if len(indices_pics) < 2:
        raise ValueError("Pas assez de pics pour déterminer la période.")

    temps_pics = temps.iloc[indices_pics].values
    differences_temps = np.diff(temps_pics)
    T0 = np.mean(differences_temps)
    w0 = 2 * np.pi / T0

    return T0, w0

resultats_pulsations = {}
for i, df in enumerate(donnees[:2]):
    try:
        T0, w0 = calculer_pulsation(df, colonne_angle='theta', temps='T')
        resultats_pulsations[f"df_{i + 1}"] = {'T0': T0, 'w0': w0}
    except ValueError as e:
        print(f"Erreur pour df_{i + 1}: {e}")

for cle, res in resultats_pulsations.items():
    print(f"{cle}: T0 = {res['T0']:.3f} s, w0 = {res['w0']:.3f} rad/s")
    
from mpl_toolkits.mplot3d import Axes3D

def espace_phases_3d_groupes(donnees, cols=['x', 'y', 'z']):
    """
    Trace les espaces des phases 3D dans deux fenêtres :
    - Une fenêtre avec une ligne et deux colonnes pour les 2 premiers DataFrames.
    - Une fenêtre avec deux lignes et deux colonnes pour les 4 derniers DataFrames.
    
    Args:
        donnees (list): Liste de DataFrames contenant les données de phase.
        cols (list): Liste des colonnes à utiliser pour l'espace des phases, ici ['x', 'y', 'z'].
    """
    # Première fenêtre (1 ligne, 2 colonnes pour les 2 premiers DataFrames)
    fig1, axes1 = plt.subplots(1, 2, figsize=(16, 8), subplot_kw={'projection': '3d'})
    
    # Tracer pour les 2 premiers DataFrames
    for i, df in enumerate(donnees[:2]):  # Les 2 premiers DataFrames
        x = df[cols[0]]
        y = df[cols[1]]
        z = df[cols[2]]
        axes1[i].scatter(x, y, z, label=f"Trajectoire - df_{i+1}", lw=0.8, s=5)
        axes1[i].set_title(f"Espace des phases 3D - df_{i+1}")
        axes1[i].set_xlabel(f"{cols[0]} (m)")
        axes1[i].set_ylabel(f"{cols[1]} (m)")
        axes1[i].set_zlabel(f"{cols[2]} (m)")
        axes1[i].legend()
    fig1.savefig(f"Espace des phases 3D - df_{i+1}.jpeg")
    # Deuxième fenêtre (2 lignes, 2 colonnes pour les 4 derniers DataFrames)
    fig2, axes2 = plt.subplots(2, 2, figsize=(16, 12), subplot_kw={'projection': '3d'})
    axes2 = axes2.flatten()  # Aplatir pour itérer facilement sur 4 axes
    
    # Tracer pour les 4 derniers DataFrames
    for i, df in enumerate(donnees[2:]):  # Les 4 derniers DataFrames
        x = df[cols[0]]
        y = df[cols[1]]
        z = df[cols[2]]
        axes2[i].scatter(x, y, z, label=f"Trajectoire - df_{i+3}", lw=0.8, s=5)
        axes2[i].set_title(f"Espace des phases 3D - df_{i+3}")
        axes2[i].set_xlabel(f"{cols[0]} (m)")
        axes2[i].set_ylabel(f"{cols[1]} (m)")
        axes2[i].set_zlabel(f"{cols[2]} (m)")
        axes2[i].legend()
    fig2.savefig(f"Espace des phases 3D - df_{i+3}.jpeg")
    # Ajuster les espacements pour une meilleure présentation
    plt.tight_layout()
    plt.show()

# Tracer les espaces des phases pour les 6 DataFrames dans la disposition demandée
espace_phases_3d_groupes(donnees)

