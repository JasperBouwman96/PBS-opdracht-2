#ifndef DYNAMICS_H_
#define DYNAMICS_H_

/* Update velocities for one time step using velocities at half the time step*/
void update_positions(struct Parameters*, struct Nbrlist*, struct Vectors*);

/* Update velocities for half a time step using forces. This function returns the kinetic energy. */
double update_velocities_half_dt(struct Parameters*, struct Nbrlist*, struct Vectors*);

/* Apply boundary conditions. In this of periodic BCs case particles are put back in the box */
void boundary_conditions(struct Parameters*, struct Vectors*);

/* Change velocities according by Berendsen thermostatting */
void berendsen_thermostat(struct Parameters*, struct Vectors*, double);

#endif /* DYNAMICS_H_ */
