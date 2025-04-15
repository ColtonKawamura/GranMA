# Theoretical Background

## Fundamentals of Granular Media

Granular materials are collections of discrete solid particles that interact through contact forces. Unlike conventional solids, liquids, or gases, granular materials exhibit complex behaviors that emerge from the collective interactions of many particles. Key characteristics include:

- **Discrete Nature**: Individual particles with finite size
- **Energy Dissipation**: Through friction and inelastic collisions
- **Force Chains**: Heterogeneous stress distribution through particle networks
- **History-Dependent Behavior**: Response depends on preparation and loading history

## Acoustic Wave Propagation in Granular Media

### Wave Types

In granular materials, we observe several types of waves:

1. **Compressional (P) Waves**: Particle motion parallel to the direction of wave propagation
2. **Shear (S) Waves**: Particle motion perpendicular to the direction of wave propagation
3. **Rotational Waves**: Associated with particle rotations
4. **Coupled Modes**: Complex combinations of the above

### Effective Medium Theory

For long-wavelength acoustic waves (where wavelength >> particle size), the granular medium can be approximated as a continuous effective medium with:

- Effective elastic moduli dependent on:
  - Coordination number (average number of contacts per particle)
  - Contact stiffness (derived from Hertzian or other contact models)
  - Confining pressure
  - Particle arrangement (fabric)

### Non-Linear Effects

Granular materials exhibit pronounced non-linear acoustic behaviors:

- **Pressure-Dependent Wave Speed**: $v \propto P^{1/6}$ for Hertzian contacts
- **Amplitude-Dependent Attenuation**: Higher amplitudes lead to increased attenuation
- **Frequency-Dependent Dispersion**: Wave speed varies with frequency
- **Scattering**: Multiple scattering from grain-scale heterogeneities

## Contact Mechanics

### Hertzian Contact Model

For spherical particles, the Hertzian contact model provides:

- Normal force: $F_n = \frac{4}{3}E^*\sqrt{R^*}\delta^{3/2}$
- Where:
  - $E^*$ is the effective elastic modulus: $\frac{1}{E^*} = \frac{1-\nu_1^2}{E_1} + \frac{1-\nu_2^2}{E_2}$
  - $R^*$ is the effective radius: $\frac{1}{R^*} = \frac{1}{R_1} + \frac{1}{R_2}$
  - $\delta$ is the overlap between particles

### Tangential Contact and Friction

The tangential force at contacts involves:

- Tangential stiffness: $k_t \propto k_n$
- Coulomb friction criterion: $|F_t| \leq \mu |F_n|$
- History-dependent tangential displacement

## Numerical Implementation

### Discrete Element Method (DEM)

Our simulations use the Discrete Element Method which:

1. Tracks positions, velocities, and orientations of all particles
2. Computes contact forces based on particle overlaps
3. Integrates Newton's equations of motion to update particle states
4. Handles boundary conditions and external forcing

### Molecular Dynamics Approach

Time integration is performed using velocity Verlet or similar algorithms with:

- Adaptive timestep selection for numerical stability
- Damping to account for energy dissipation
- Specialized routines for efficient contact detection

## References

1. WIP!