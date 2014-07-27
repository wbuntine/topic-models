/*
 *   Binary heap code from
 *     https://gist.github.com/martinkunev/1365481
 *   modified for use by Wray Buntine
 * 
 *   Heap stored in data[] vector with root at data[0]
 */

#include <unistd.h>
#include <stdlib.h>
#include <assert.h>

#include "heap.h"

static int base_size=4;

// Prepares the heap for use
int heap_zero(struct heap_s *h,
               int (*cmp)(uint32_t k1, uint32_t k2, void *par),
               void *par)
{
  *h = (struct heap_s){
    .size = base_size,
    .count = 0,
    .data = malloc(sizeof(h->data[0]) * base_size),
    .lookup = malloc(sizeof(h->data[0]) * base_size),
    .cmp = cmp,
    .par = par
  };
  if (!h->data || !h->lookup) return 1;
  return 0;
}

int heap_init(struct heap_s *h,
              uint32_t *data,  uint32_t count,
              int (*cmp)(uint32_t k1, uint32_t k2, void *par),
              void *par)
{
  uint32_t item, index, swap, other;
  uint32_t temp;
  
  *h = (struct heap_s){
    .size = count,
    .count = count,
    .data = data,
    .lookup = malloc(sizeof(h->data[0]) * count),
    .cmp = cmp,
    .par = par
  };
  if (!h->lookup) return 1;
  /*
   *  Heapifies a non-empty array
   */
  // Move every non-leaf element to the right position in its subtree
  item = (count >> 1) - 1;
  while (1) {
    // Find the position of the current element in its subtree
    temp = data[item];
    for (index = item; 1; index = swap) {
      // Find the child to swap with
      swap = (index << 1) + 1;
      // If there are no children, the current element is positioned
      if (swap >= count) 
        break; 
      other = swap + 1;
      if ((other < count) && 
          cmp(data[other], data[swap], par)) swap = other;
      // If bigger child is smaller than or equal to parent, is reordered
      if ( cmp(temp, data[swap], par) ) 
        break;    
      data[index] = data[swap];
    }
    if (index != item) data[index] = temp;
    if (!item) break;
    --item;
  }
  /*
   *  now build lookup
   */
  for (item=0; item<count; item++) {
    assert(data[item]<count);
    h->lookup[data[item]] = item;
  }
  return 0;
}

// Inserts element to the heap
int heap_push(struct heap_s *h, uint32_t value)
{
  uint32_t index, parent;
  assert(value<h->count);
  
  // Resize the heap if it is too small to hold all the data
  if (h->count == h->size)
    {
      h->size <<= 1;
      h->data = realloc(h->data, sizeof(h->data[0]) * h->size);
      if (!h->data) return 1; // Exit if the memory allocation fails
    }
  
  // Find out where to put the element and put it
  for(index = h->count++; index; index = parent) {
    parent = (index - 1) >> 1;
    if ( h->cmp(h->data[parent], value, h->par) )
      break;
    h->data[index] = h->data[parent];
    h->lookup[h->data[index]] = index;
  }
  h->data[index] = value;
  h->lookup[value] = index;
  return 0;
}

void heap_free(struct heap_s *h) {
  free(h->lookup);
  free(h->data);
}

/*
 *   Removes a particular element from the heap, 
 *      - find using .lookup
 *      - ripple up to root (pretending is max)
 *      - then do a heap_pop()
 */
int heap_remove(struct heap_s *h, uint32_t value) 
{
  uint32_t index, parent;
  uint32_t temp;
  
  if ( h->lookup[value]>=h->size ) 
    return 1;
  
  // Find out where to put the element and put it
  for (index = h->lookup[value]; index; index = parent) {
    parent = (index - 1) >> 1;
    h->data[index] = h->data[parent];
    h->lookup[h->data[index]] = index;
  }
  h->lookup[value] = UINT32_MAX;
  heap_pop(h);
  return 0;
}

// Removes the root element from the heap
int heap_pop(struct heap_s *h)
{
  uint32_t index, swap, other;
  // Remove the biggest element
  uint32_t temp = h->data[--h->count];
    
  h->lookup[h->data[0]] = UINT32_MAX;

  // Resize the heap if it's consuming too much memory
  if ((h->count <= (h->size >> 2)) && (h->size > base_size)) {
    h->size >>= 1;
    h->data = realloc(h->data, sizeof(h->data[0]) * h->size);
    if (!h->data) return 1; // Exit if the memory allocation fails
  }
  
  // Reorder the elements
  for(index = 0; 1; index = swap) {
    // Find the child to swap with
    swap = (index << 1) + 1;
    // If there are no children, the heap is reordered
    if (swap >= h->count) 
      break; 
    other = swap + 1;
    if ((other < h->count) 
        && h->cmp(h->data[other], h->data[swap], h->par)) 
      swap = other;
    // If the bigger child is less than or equal to its parent, 
    // the heap is reordered
    if ( h->cmp(temp, h->data[swap], h->par)  )
      break;                
    h->data[index] = h->data[swap];
    h->lookup[h->data[index]] = index;
  }
  h->data[index] = temp;
  h->lookup[temp] = index;
  return 0;
}

#ifdef HEAP_MAIN
#include <stdio.h>

void heap_print(struct heap_s *h) {
  int i;
  int p2=2;
  printf("Heap count=%u, size=%u\n1: ", h->count, h->size);
  for (i=0; i<h->count; i++) {
    printf(" %u", h->data[i]);
    if ( i+1==p2-1 && i<h->count-1 ) {
      p2 <<= 1;
      printf("\n%d: ",p2/2);
    }
  }
  printf("\n");
  printf("Lookup: ");
  for (i=0; i<h->size; i++) {
    printf(" %u", h->lookup[i]);
  }
  printf("\n");
}

#define DIM 20

int icmp(uint32_t k1, uint32_t k2, void *par) {
  if ( k2<=k1 )
    return 1;
  return 0;
}

int main(int argc, char* argv[]) {
  struct heap_s heap;
  uint32_t data[DIM];
  int i,j;

  for (i=0; i<DIM; i++)
    data[i] = i;
  for (i=0; i<DIM; i++) {
    uint32_t swap;
    j = rand()%DIM;
    if ( i!=j ) {
      swap = data[i];
      data[i] = data[j];
      data[j] = swap;
    }
  }

  printf("data:");
  for (i=0; i<DIM; i++)
    printf(" %u", data[i]);
  printf("\n");

  heap_init(&heap, &data[0], DIM, icmp, NULL);
  heap_print(&heap);

  printf("top: %u\n", heap_front(&heap));
  heap_pop(&heap);
  printf("top: %u\n", heap_front(&heap));
  heap_pop(&heap);
  printf("top: %u\n", heap_front(&heap));
  heap_print(&heap);
  printf("removing %u\n",data[5]);
  heap_remove(&heap,data[5]);
  heap_print(&heap);
  printf("removing %u\n",data[6]);
  heap_remove(&heap,data[6]);
  heap_print(&heap);
  printf("removing %u\n",data[12]);
  heap_remove(&heap,data[12]);
  heap_print(&heap);

  return 0;
}

#endif
