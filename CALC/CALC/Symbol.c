/* symbol.c */

#define _CRT_SECURE_NO_DEPRECATE
#define _CRT_NONSTDC_NO_DEPRECATE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "integer.h"
#include "calc.h"
#include "fun.h"
#include "stack.h"
#ifdef _WIN32
#include "ytab.h"
#else
#include "y.tab.h"
#endif

static Symbol *symlist = NULL; /* symbol table */


Argument createArg(void *arg, int type, int defined)
{
	Symbol *sym;
	Argument Arg;
	if (!(Arg = malloc(sizeof(struct _Argument))))
	{
		printf("Not enough memory to allocate argument\n");
		exit(1);
	}
	sym = (Symbol *)arg;
	Arg->defined = defined;
	switch (type) {
	case NUM:
		Arg->type = NUM;
		Arg->u.num = arg;

		break;
	case VARADR:
		Arg->type = VARADR;
		Arg->u.varAdr = &(sym->u.symval);

		break;
	case ARR:
		Arg->type = ARR;
		/* in this case a symbol has been passed to createArg. it's ok */
		Arg->u.array = sym->u.symarr;
		break;
	case ARRADR:
		Arg->type = ARRADR;
		Arg->u.arrayAdr = &(sym->u.symarr);
		break;
	case POLY:
		Arg->type = POLY;
		Arg->u.poly = arg;
		break;
	}
	return Arg;
}

void freeArg(Argument Arg)
{

	switch ((Arg)->type) {
	case NUM:
		FREEMPI((Arg)->u.num);
		break;
	case POLY:
		DELETEPI((Arg)->u.poly);
		break;
	}
	free(Arg);
}


Symbol *lookup(char *s, int typ)  /* find s of type  typ  in symbol table */
{
	Symbol *sp;

	for (sp = symlist; sp != NULL; sp = sp->next)
		if (strcmp(sp->name, s) == 0 && (sp->type == typ))
			return sp;
	return NULL;  /* NULL ==> not found */
}

/* add function to symlist. s = name, t=type */
Symbol *installFunc(const char *s, int t, const int *argTypes)
{
	Symbol *sp;

	sp = (Symbol *)mmalloc(sizeof(Symbol));
	sp->name = (char *)mmalloc(1 + strlen(s));
	strcpy(sp->name, s);
	sp->type = t;
	sp->argTypes = (int *)argTypes;
	sp->next = symlist;  /* put at front of list */
	symlist = sp;
	return sp;
}


Symbol *install(char *s, int t)  /* install s in symbol table */
{
	Symbol *sp;

	sp = (Symbol *)mmalloc(sizeof(Symbol));
	sp->name = (char *)mmalloc(1 + strlen(s));
	strcpy(sp->name, s);
	sp->type = t;
	sp->argTypes = NULL;
	sp->next = symlist;  /* put at front of list */
	symlist = sp;
	return sp;
}

void clean_symtab()
/* deallocates all memory allocated during the CALC session. */
{
	Symbol *tmp;
	int typ;


	while (symlist)
	{
		tmp = symlist->next;
		typ = symlist->type;
		if (typ == VAR && symlist->u.symval != NULL)
			FREEMPI(symlist->u.symval);
		if (typ == ARRAY)
		{
			FREEMPIA(symlist->u.symarr);
		}
		if (typ == POLYVAR && symlist->u.sympval != NULL) {
			DELETEPI(symlist->u.sympval);
		}
		ffree((char*)(symlist->name), 1 + strlen(symlist->name));
		ffree((char*)symlist, sizeof(Symbol));
		symlist = tmp;
	}
}
