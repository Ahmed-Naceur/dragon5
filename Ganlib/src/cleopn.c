
/*****************************************/
/*             CLE-2000 API              */
/*    LIFO stack utility for CLE-2000    */
/*     AUTHOR: A. Hebert ; 27/07/10      */
/*****************************************/

#include <stdlib.h>
#include <string.h>
#include "cle2000.h"

void cleopn(lifo **my_lifo)
{
   (*my_lifo) = (lifo *) malloc(sizeof(lifo));
   (*my_lifo)->nup = 0;
   (*my_lifo)->root = NULL;
   (*my_lifo)->node = NULL;
}
lifo_node * clepop(lifo **my_lifo)
{
   lifo_node *my_node;
   if ((*my_lifo)->nup == 0) return NULL;
   my_node = (*my_lifo)->node;
   (*my_lifo)->node = my_node->daughter;
   my_node->daughter = NULL;
   (*my_lifo)->nup--;
   return my_node;
}
void clepush(lifo **my_lifo, lifo_node *my_node)
{
   lifo_node *daughter_node;
   if ((*my_lifo)->nup == 0) {
      (*my_lifo)->root = my_node;
      daughter_node = NULL;
   } else {
      daughter_node = (*my_lifo)->node;
   }
   (*my_lifo)->node = my_node;
   (*my_lifo)->node->daughter = daughter_node;
   (*my_lifo)->nup++;
}
int_32 clecls(lifo **my_lifo)
{
   if ((*my_lifo)->nup != 0) return -1;
   free(*my_lifo);
   (*my_lifo) = NULL;
   return 0;
}
lifo_node * clenode(lifo **my_lifo, const char *name)
{
   lifo_node *my_node;
   my_node = (*my_lifo)->node;
   if (my_node == NULL) return NULL;
   while (my_node->daughter != NULL) {
      if (strcmp(my_node->name, name) == 0) return my_node;
      my_node = my_node->daughter;
   }
   if (strcmp(my_node->name, name) == 0) return my_node;
   return NULL;
}
lifo_node * clepos(lifo **my_lifo, int_32 ipos)
{
   lifo_node *my_node;
   int_32 iloop;
   if ((ipos > (*my_lifo)->nup - 1) || (ipos < 0)) return NULL;
   my_node = (*my_lifo)->node;
   for (iloop = 0; iloop < (*my_lifo)->nup - ipos - 1; ++iloop) {
      my_node = my_node->daughter;
   }
   return my_node;
}
void clelib(lifo **my_lifo)
{
   lifo_node *my_node;
   int_32 iloop;
   printf("\n lifo content:\n node     type   name........\n");
   for (iloop = 0; iloop < (*my_lifo)->nup; ++iloop) {
      my_node = clepos(my_lifo, iloop);
      printf(" %4d    (%4d)  %s\n", (int)iloop, (int)my_node->type, my_node->name);
   }
}
