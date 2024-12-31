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
        b=gamma/2
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


def calculer_pulsation_et_amortissement(df, colonne_angle='theta', temps='T'):
    # Extraction des angles et du temps
    angles = df[colonne_angle]
    temps = df[temps] / 100  # Conversion en secondes

    # Détection des pics
    indices_pics, _ = find_peaks(angles)

    # Vérification du nombre de pics
    if len(indices_pics) < 2:
        raise ValueError("Pas assez de pics pour déterminer la période et l'amortissement.")

    # Extraction des temps des pics
    temps_pics = temps.iloc[indices_pics].values

    # Calcul des différences de temps entre les pics
    differences_temps = np.diff(temps_pics)

    # Calcul de la période moyenne
    T0 = np.mean(differences_temps)  # Pas de correction pour T0/2 ici
    w0 = 2 * np.pi / T0

    # Calcul de l'amortissement
    amplitudes = np.abs(angles.iloc[indices_pics].values)  # Prendre la valeur absolue des amplitudes
    if len(amplitudes) < 2:
        raise ValueError("Pas assez d'amplitudes pour calculer l'amortissement.")

    ln_ratios = np.log(amplitudes[:-1] / amplitudes[1:])
    b = np.abs(np.mean(ln_ratios) / T0)

    return T0, w0, b

# Utilisation
resultats = {}

for i, df in enumerate(donnees[:2]):  # Limité aux 2 premiers DataFrames
    try:
        # Calcul T0, w0 et b en une seule fonction
        T0, w0, b = calculer_pulsation_et_amortissement(df, colonne_angle='theta', temps='T')
        
        # Stockage des résultats
        resultats[f"df_{i + 1}"] = {'T0': T0, 'w0': w0, 'b': b}
    except ValueError as e:
        print(f"Erreur pour df_{i + 1}: {e}")

# Affichage des résultats
for cle, res in resultats.items():
    print(f"{cle}: T0 = {res['T0']:.3f} s, w0 = {res['w0']:.3f} rad/s, b = {res['b']:.3f} N·s/m")



    
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

def distance_aux_attracteurs(x, y, attracteurs):
    """
    Calcule la distance minimale d'un point (x, y) à l'un des attracteurs.
    """
    distances = [np.sqrt((x - ax)**2 + (y - ay)**2) for ax, ay in attracteurs]
    return min(distances)

def visualiser_bassins_attraction_avec_couleur(donnees, attracteurs):
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    
    for i, df in enumerate(donnees[2:]):  # Trajectoires des 4 derniers DataFrames (chaotiques)
        x = df['x']
        y = df['y']
        
        # Calcul des distances à l'attracteur le plus proche pour chaque point
        distances = [distance_aux_attracteurs(xi, yi, attracteurs) for xi, yi in zip(x, y)]
        
        # Normalisation des distances pour l'échelle de couleur
        norm = plt.Normalize(min(distances), max(distances))
        cmap = plt.get_cmap("coolwarm")
        
        # Trace les trajectoires en fonction des distances avec une échelle de couleur
        ax = axes[i // 2, i % 2]  # Positionner les graphiques sur une grille 2x2
        scatter = ax.scatter(x, y, c=distances, cmap=cmap, norm=norm, s=10)
        ax.set_title(f"Portrait de phase (x,y) - df_{i+3}")
        ax.set_xlabel("x (m)")
        ax.set_ylabel("y (m)")
        ax.grid()
        
        # Ajouter la barre de couleur
        fig.colorbar(scatter, ax=ax, label="Distance aux attracteurs")
    plt.savefig(f"Portrait de phase (x,y) - df_{i+3}.jpeg")
    plt.tight_layout()
    plt.show()

# Liste des attracteurs (en coordonnées x, y)
attracteurs = [(-0.03, 0), (0.015, -0.026), (0.015, 0.026)]

# Visualiser les bassins d'attraction pour les 4 derniers DataFrames avec couleurs en fonction de la distance
visualiser_bassins_attraction_avec_couleur(donnees, attracteurs)


def visualiser_bassins_attraction_combinés(donnees, attracteurs, colormap="viridis", gamma=0.5):
    """
    Visualise les bassins d'attraction combinés pour les 4 derniers DataFrames
    avec des couleurs nuancées en fonction des distances aux attracteurs.
    
    Args:
        donnees (list): Liste de DataFrames contenant les colonnes x, y.
        attracteurs (list): Liste de tuples (x, y) représentant les attracteurs.
        colormap (str): Nom du colormap à utiliser.
        gamma (float): Facteur d'ajustement non linéaire pour accentuer les nuances.
    """
    # Créer une figure pour tout afficher dans un seul graphique
    plt.figure(figsize=(10, 8))
    
    # Initialiser les listes pour stocker toutes les positions (x, y) et distances
    x_combined = []
    y_combined = []
    distances_combined = []
    
    for df in donnees[2:]:  # Trajectoires des 4 derniers DataFrames (chaotiques)
        x_combined.extend(df['x'])
        y_combined.extend(df['y'])
        
        # Calcul des distances à l'attracteur le plus proche pour chaque point
        for xi, yi in zip(df['x'], df['y']):
            distances_combined.append(distance_aux_attracteurs2(xi, yi, attracteurs))
    
    # Appliquer une transformation gamma pour mieux nuancer les couleurs
    distances_combined = np.power(distances_combined, gamma)
    
    # Normalisation des distances pour l'échelle de couleur
    norm = plt.Normalize(vmin=min(distances_combined), vmax=max(distances_combined))
    
    # Trace les trajectoires en fonction des distances avec une échelle de couleur
    scatter = plt.scatter(x_combined, y_combined, c=distances_combined, cmap=colormap, norm=norm, s=10)
    plt.title("Bassins d'attraction combinés")
    plt.xlabel("x (m)")
    plt.ylabel("y (m)")
    plt.grid()
    
    # Ajouter la barre de couleur
    plt.colorbar(scatter, label="Distance (ajustée) aux attracteurs")
    plt.savefig("Bassins d'attraction combinés.jpeg")
    plt.tight_layout()
    plt.show()

def distance_aux_attracteurs2(x, y, attracteurs):
    """
    Calcule la distance minimale entre un point (x, y) et une liste d'attracteurs.
    
    Args:
        x (float): Coordonnée x du point.
        y (float): Coordonnée y du point.
        attracteurs (list): Liste de tuples (x, y) représentant les attracteurs.
    
    Returns:
        float: Distance minimale au plus proche attracteur.
    """
    return min(np.sqrt((ax - x)**2 + (ay - y)**2) for ax, ay in attracteurs)

# Exemple d'appel avec une nuance accentuée et un colormap différent
visualiser_bassins_attraction_combinés(donnees, attracteurs, colormap="plasma", gamma=0.8)

