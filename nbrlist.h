#ifndef NBRLIST_H_
#define NBRLIST_H_

#include <stddef.h>

/* Allocate arrays needed to store the cell-linked-list data*/
void alloc_celllist(struct Parameters*, struct Celllist*);

/* Free arrays used for the cell-linked-list */
void free_celllist(struct Celllist*);

/* Build the cell-linked-list */
void build_celllist(struct Parameters*, struct Vectors*, struct Celllist*);

/* Allocate arrays needed to store the neighbor list */
void alloc_nbrlist(struct Parameters*, struct Nbrlist*);

/* Free arrays use to store the neighbor list */
void free_nbrlist(struct Nbrlist*);

/* Build the neighbor list */
void build_nbrlist(struct Parameters*, struct Vectors*, struct Nbrlist*);

/* Update the connecting vectors of the neighbor list and if needed rebuild it.*/
int update_nbrlist(struct Parameters*, struct Vectors*, struct Nbrlist*);

#endif /* NBRLIST */
