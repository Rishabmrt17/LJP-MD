#!/usr/bin/env python
# coding: utf-8

# In[5902]:


#import modules
import numpy as np
import matplotlib.pyplot as plt


# In[5903]:


def run_md(N_atoms, r, forces, T, m, size, step, delta_step, dt):
    """
    MD program using velocity verlet algorithm
    N_atoms = number_of_atoms
    r = distance_between_particles
    forces = forces_on_the_particles
    T = Boltzmann_temperature
    m = mass_of_particle
    size = size_of_the_box
    step = number_of_steps
    delta_step = frequency_in_steps
    dt = time_step
    """
        # Main MD loop
          # open trajectory file
    file = open("traj_MD.xyz", 'w')
        # initialize positions
    positions = initialize_positions(N_atoms, size)
    
        # initialize velocities
    velocities = initialize_velocities(m, kB, T)
    
    for steps in range(step):
    
        # Propagate Positions
        update_positions(positions, velocities, forces, dt, size)
        
        # Propagate Velocities
        update_velocities(velocities, forces, dt)
        
        if (step%delta_step==0):
            write_trajectory_frame(positions, file, step)
    
    # close trajectory file
    file.close()
    return(file)

# In[5904]:


#################Sub-Routines#################

# initialize positions
def initialize_positions(N_atoms, size):
    """Initialize positions"""
    positions = np.random.rand(N_atoms, 3) * size
    
    return positions


# In[5905]:


N_atoms = 2
size = 10
positions = np.random.rand(N_atoms, 3) * size
print(positions)


# In[5906]:


def lj_potential(r, epsilon, sigma):
    """Lennard-Jones potential function"""
    potential = 4 * epsilon * ((sigma / r)**12 - (sigma / r)**6)
    return potential


# In[5907]:


def compute_forces(positions, f, size):
    """
    Compute forces
    epsilon = welldepth of the particles
    sigma = distance at which potential is zero
    r = distance between two particles
    """
    r6 = (sigma / r) ** 6
    r12 = r6 ** 2
    forces = -24 * epsilon * (2 * r12 - r6) / r
    return forces


# In[5908]:


sigma = 1
epsilon = 1
r = 1
forces = -24 * epsilon * (2 * r12 - r6) / r
print(forces)


# In[5909]:


# initialize velocities
def initialize_velocities(m, kB, T):
    """
    kB = boltzmann_temp
    """
    velocities = np.random.normal(loc=0, scale=1, size=(N_atoms, 3)) * np.sqrt(kB*T/m)
    
    return velocities


# In[5910]:


N_atoms = 2
kB = 0.08314
T = 300
m = 12
velocities = np.random.normal(loc=0, scale=np.sqrt(T), size=(N_atoms,3)) * np.sqrt(kB*T/m)
print(velocities)


# In[5911]:


# Propagate Positions
def update_positions(positions, velocities, dt, forces, size):
    """
    positions = particle_positions
    velocities = particle_velocities
    dt = time_step
    forces = force_on_particles
    size = size_of_box
    """
    
    positions += velocities * dt + 0.5 * forces * dt**2
    
     # wrap into central box (box is from 0 to size in each dimension)
    for i in range(N):
        for j in range(N):
            if positions[i,j] < 0:
                positions[i,j] += size
            elif positions[i,j] > size:
                positions[i,j] -= size
    return positions


# In[5912]:


dt = 0.01
size = 10
positions += velocities * dt + 0.5 * forces * dt**2


# In[5913]:


print(positions)


# In[5914]:


# Propagate Velocities
def update_velocities(velocities, forces, dt):
    """
    velocities = particle_velocities
    forces = forces
    dt = time_step
    
    """
    
    velocities += forces * dt
    return velocities


# In[5915]:


dt = 0.01
velocities += forces * dt
print(velocities)


# In[5916]:


# Compute Energy
def kinetic_energy(N_atoms, kB, T):
    """kB - Boltzmann constant"""
    
    Energy = 1.5*N_atoms*kB*T
    return Energy


# In[5917]:


N_atoms = 2
kB = 0.08314
T = 300
Energy = 1.5*N_atoms*kB*T
print(Energy)


# In[5918]:


# Trajectory frame        
def write_trajectory_frame(positions, file_pointer, step):
    """
    positions = particle_positions
    file_pointer = trajectory_file_pointer
    step = step_number
    """
    for i in range(N_atoms):
        file_pointer.write("%10.5f %10.5f %10.5f\n" % ( positions[i,0],  positions[i,1], positions[i,2]))


# In[5919]:


#test sim
N_atoms = 2
T = 300
m = 12
size = 10
step = 100
delta_step = 10
dt = 0.01
sim = run_md(N_atoms, r, forces, T, m, size, step, delta_step, dt)


# In[5920]:


print(sim)


# In[5921]:


print(file)


# In[5922]:


def compute_radial_distribution(positions, size, N_atoms, num_bins):
    """Compute the radial distribution function (g(r)) for a given set of coordinates."""
    rdf = np.zeros(num_bins)
    bin_width = size / (2 * num_bins)
    positions = np.random.rand(N_atoms,3)
    for i in range(N_atoms - 1):
        for j in range(i + 1, N_atoms):
            dx = positions[i, 0] - positions[j, 0]
            dy = positions[i, 1] - positions[j, 1]

            # Apply periodic boundary conditions
            dx -= size * np.round(dx / size)
            dy -= size * np.round(dy / size)

            r = np.sqrt(dx ** 2 + dy ** 2)
            
            # Assign particles to bins and increment the respective bin count
            bin_index = int(r / bin_width)
            if bin_index < num_bins:
                rdf[bin_index] += 2  # Increase count for pair (i, j)
    
    # Normalize the radial distribution function
    density = N_atoms / (size ** 3)
    shell_volume = (4 / 3) * np.pi * (bin_width ** 3) * np.arange(1, num_bins + 1) ** 3
    normalization = density * shell_volume
    rdf /= normalization
    
    return rdf


# In[5927]:


# Simulation parameters
N_atoms = 2
num_bins = 100
size = 10

rdf = compute_radial_distribution(positions, size, N_atoms, num_bins)
    
r = np.linspace(0, size / 2, num_bins)


# In[5928]:


print(rdf)


# In[5929]:


print(r)


# In[5930]:


plt.plot(r,rdf)





