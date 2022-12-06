/* File: stack.c
* Desc: implementation of the Stack interface using opaque types.
* Cristina Cifuentes
* 13 Aug 1997
*/

#define _CRT_SECURE_NO_DEPRECATE
#define _CRT_NONSTDC_NO_DEPRECATE
#include <assert.h>
#include <stdlib.h>
#include "stack.h"

struct _Stack {
	int count;			/* # elements */
	struct elem {
		void        *x;		/* actual element */
		struct elem *link;	/* next element of the stack */
	} *head;
};


Stack stackNew(void)
{
	Stack stk;

	stk = (Stack)malloc(sizeof(struct _Stack));
	stk->count = 0;
	stk->head = NULL;
	return stk;
}


int stackEmpty(Stack stk)
{
	assert(stk);
	return (stk->count == 0);
}


void stackPush(Stack stk, void *x)
{
	struct elem *t;

	assert(stk);
	t = (struct elem *)malloc(sizeof(struct elem));
	t->x = x;
	t->link = stk->head;
	stk->head = t;
	stk->count++;
}


void *stackPop(Stack stk)
{
	void *x;
	struct elem *t;

	assert(stk);
	assert(stk->count > 0);
	t = stk->head;
	stk->head = t->link;
	stk->count--;
	x = t->x;
	free(t);
	t = NULL;
	return x;
}


void stackFree(Stack *stk)
{
	struct elem *t, *u;

	assert(stk && *stk);
	for (t = (*stk)->head; t; t = u)
	{
		u = t->link;
		free(t);
		t = NULL;
	}
	free(*stk);
	*stk = NULL;
}


