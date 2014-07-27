/*
 *   Binary heap code from
 *     https://gist.github.com/martinkunev/1365481
 *   modified for use by Wray Buntine
 * 
 *   Heap stored in data[] vector with root at data[0]
 *   The heap_remove() function requires keeping a parallel
 *   inverted index to locate the element to remove.
 *
 *   WARNING: Requires C99 compatible compiler
 */

#ifndef __HEAP_H
#define __HEAP_H

#include <stdint.h>

struct heap_s
{
  uint32_t size;             // size of allocated array
  uint32_t count;            // number of elements
  uint32_t *data;            // Array with the elements
  uint32_t *lookup;          // where to find entry
  /*
   *   comparison function; return true if k1 is
   *   better or as good as k2
   */
  int (*cmp)(uint32_t k1, uint32_t k2, void *par);  
  void *par;                 //  auxiliary data for cmp()
};

/*
 *   all four return non-zero on error,
 *   usually out of memory
 */
int heap_zero(struct heap_s *h,
              int (*cmp)(uint32_t k1, uint32_t k2, void *par),
              void *par);
int heap_init(struct heap_s *h,
              uint32_t *data,  uint32_t count,
              int (*cmp)(uint32_t k1, uint32_t k2, void *par),
              void *par);
int heap_push(struct heap_s *h, uint32_t value);
int heap_pop(struct heap_s *h);
int heap_remove(struct heap_s *h, uint32_t value);

void heap_free(struct heap_s *h);
void heap_print(struct heap_s *h);

// Returns the biggest element in the heap
#define heap_front(h) (*(h)->data)

#endif
