#ifndef MEMORY_H_
#define MEMORY_H_

/* Allocate the arrays in 'vectors' needed to store information of all particles */
void alloc_vectors(struct Vectors *vectors, size_t N);

/* Free the arrays in 'vectors' */
void free_vectors(struct Vectors *vectors);

/* Allocate all variables needed in the MD simulation */
void alloc_memory(struct Parameters*, struct Vectors*, struct Nbrlist*);

/* Free the memory allocated by alloc_memory */
void free_memory(struct Vectors *, struct Nbrlist*);

#endif /* MEMORY_H_ */
