#ifndef YYERRCODE
#define YYERRCODE 256
#endif

#define NUMBER 257
#define POL 258
#define ARGSTACK 259
#define ANARRAY 260
#define ARRAY 261
#define POLYVAR 262
#define VAR 263
#define POLYTERM 264
#define BLTIN 265
#define BLTINV 266
#define BLTINP 267
#define UNDEF 268
#define UNARYMINUS 269
typedef union {
	MPI *val;  /* actual value */
	POLYI pol;
	Stack argStack;
	MPIA arr;
	Symbol *sym; /* symbol table pointer */
} YYSTYPE;
extern YYSTYPE yylval;

