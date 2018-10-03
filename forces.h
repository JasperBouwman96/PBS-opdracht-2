#ifndef FORCES_H_
#define FORCES_H_

/* Compute forces on particles using the pairs in a neighbor list.
This function returns the total potential energy of the system. */
double calculate_forces(struct Parameters*, struct Nbrlist*, struct Vectors*);

#endif /* FORCES_H_ */
