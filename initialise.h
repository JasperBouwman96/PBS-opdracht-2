#ifndef INITIALISE_H_
#define INITIALISE_H_

/* Initialize all the variables of an MD simulation. */
void initialise_variables(struct Parameters*, struct Vectors*, size_t*, double*);

/* Initialize positions of all particles */
void initialise_positions(struct Parameters*, struct Vectors*);

/* Initialize the velocities of all particles */
void initialise_velocities(struct Parameters*, struct Vectors*);

#endif /* INITIALISE_H_ */
