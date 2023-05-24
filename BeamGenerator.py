# Generates ions based on a specified beam matrix
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
sns.set()

# Set random seed to ensure repeatability
rng = np.random.default_rng(1)

# Generate phase pts
def getPhasePts(mat, N, plotPhase = False, dist='Gaussian', name=None, axis_label=None, color='skyblue'):
    if dist == 'Gaussian':
        A = rng.multivariate_normal([0, 0], 1/6 * mat, N)
        
        pos = A[:, 0]
        mom = A[:, 1]
                
    if plotPhase:
        plt.scatter(pos[:] * 1e3, mom[:] * 1e3, label=name, color = color)

    return pos, mom       # Generated phase points

""" Elinas input"""
r11 = 6.25000000e-06
r12 = 0
r22 = 4.00e-06
r33 = 6.25000000e-06
r22 = 4.00e-06
r33 = 6.25e-06
r34 = 0
r44 = 4.00e-06


"""
# Cosy Result for beam 1:
r11 = 0.4915759459186747e-4
r12 = -0.8691129663384992e-4
r22 = 0.1663742549074718e-3
r33 = 0.1668938718754474e-4
r34 = -0.4689750558581079e-4
r44 = 0.1692309423420552e-3

# Cosy Result for beam 2:
r11 = 0.2360717709791814e-4
r12 = -0.4179106819195326e-4
r22 = 0.1004558089040632e-3
r33 = 0.1850465441691150e-4
r34 = -0.3635787204502700e-4
r44 = 0.1052102625894840e-3


# Cosy Result for beam 3:
r11 = 0.1807039611590347e-4
r12 = -0.1908721197332630e-4
r22 = 0.2154472200927816e-4
r33 = 0.1807039611590347e-4
r34 = -0.1908721197332630e-4
r44 = 0.2154472200927816e-4
"""

# Longitudinal beam matrix
r55 = 0.021025          # Bunch length (seconds (?))
r56 = 0.00723832        # Coupling [l - E]
r66 = 0.00249575        # Energy spread (percent)

meanE = 52000   # Specify mean energy (eV)
mass = 39 # a.m.u.

r = np.zeros((6, 6))

# X
r[0, 0] = r11
r[0, 1] = r12
r[1, 1] = r22
r[1, 0] = r12

# Y
r[2, 2] = r33
r[2, 3] = r34
r[3, 3] = r44
r[3, 2] = r34

# l
r[4, 4] = r55
r[4, 5] = r56
r[5, 5] = r66
r[5, 4] = r56

# Number of ions to generate:
N = 350

plt.figure(figsize=(28/2.54, 14/2.54))
plt.suptitle("Input beam for PDT simulations")
plt.subplot(121)

# Generate and plot x and y phase portraits
x, px = getPhasePts(np.array(r[:2, :2]), N, plotEig=False, plotPhase=True, name='X', color='skyblue')
y, py = getPhasePts(np.array(r[2:4, 2:4]), N, plotEig=False, plotPhase=True, name='Y', color='indianred')
plt.ylabel('$p_x$ [mrad]')
plt.xlabel('$x$ [mm]')

plt.legend()
plt.subplot(122)

# Generate and plot longitudinal phase portraits
l, dE = getPhasePts(np.array(r[4:, 4:]), N, plotEig=False, plotPhase=True, color='skyblue')

plt.ylabel('$E/\Delta E$')
plt.xlabel('$\ell$ [mm]')

plt.savefig('BeamGeneration/GeneratedBeam.png')
plt.savefig('BeamGeneration/GeneratedBeam.eps')

# Translate result into readable input for simion
E = meanE * (1 + dE)
v0 = np.sqrt(2 * E / mass * 9.64853322e7)
vx = px * v0
vy = py * v0
vz = np.sqrt(v0**2 - vx**2 - vy**2)

result = np.array([x, y, l, vx, vy, vz, E]).T
np.savetxt('BeamGeneration/GeneratedParticles.csv', result, delimiter=';', header="x (m); y (m); l (m); vx (m/s); vy (m/s); vz (m/s); E")

