/* File: stack.h
* Desc: interface of a stack ADT using an opaque type.
* Cristina Cifuentes
* 13 Aug 1997
*/

#ifndef _STACK_H_
#define _STACK_H_

typedef struct _Stack *Stack;

Stack	stackNew(void);
int	stackEmpty(Stack);
void	stackPush(Stack, void *);
void   *stackPop(Stack);
void    stackFree(Stack *);

#endif
