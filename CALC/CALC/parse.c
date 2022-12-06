#ifndef lint
static char const
yyrcsid[] = "$FreeBSD: src/usr.bin/yacc/skeleton.c,v 1.28 2000/01/17 02:04:06 bde Exp $";
#endif
#include <stdlib.h>
#define YYBYACC 1
#define YYMAJOR 1
#define YYMINOR 9
#define YYLEX yylex()
#define YYEMPTY -1
#define yyclearin (yychar=(YYEMPTY))
#define yyerrok (yyerrflag=0)
#define YYRECOVERING() (yyerrflag!=0)
static int yygrowstack();
#define YYPREFIX "yy"
// #line 11 "parse.y"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "integer.h"
#include "fun.h"
#include "calc.h"
#include "stack.h"

#define POLYID 'X'
union {
	MPI *mpi;
	POLYI polyi;
} retval;



int  rettype; /* 0 for nothing, 1 for MPI, 2 for POLYI, 3 for failure of
			  * function */

int yyparse(void);
int yylex(void);


// #line 38 "parse.y"
typedef union {
	MPI *val;  /* actual value */
	POLYI pol;
	Stack argStack;
	MPIA arr;
	Symbol *sym; /* symbol table pointer */
} YYSTYPE;
// #line 51 "y.tab.c"
#define YYERRCODE 256
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
const short yylhs[] = { -1,
0,    0,    0,    0,    0,    0,    0,    0,    0,    4,
4,    5,    5,    5,    5,    5,    5,    5,    5,    5,
5,    5,    5,    5,    5,    5,    5,    5,    5,    5,
5,    5,    6,    6,    7,    7,    2,    1,    1,    1,
1,    1,    1,    1,    1,    1,    1,    1,    1,    1,
1,    1,    3,    3,    9,    9,    8,    8,    8,    8,
8,    8,    8,    8,    8,    8,
};
const short yylen[] = { 2,
0,    2,    3,    3,    3,    3,    3,    3,    3,    3,
3,    1,    1,    2,    4,    3,    1,    3,    3,    3,
2,    3,    3,    3,    3,    3,    3,    3,    3,    3,
3,    2,    6,    6,    2,    3,    3,    1,    1,    4,
4,    4,    1,    3,    3,    3,    3,    3,    3,    2,
3,    2,    2,    3,    2,    2,    2,    2,    4,    5,
3,    4,    3,    3,    5,    6,
};
const short yydefred[] = { 1,
0,    0,   38,    0,    0,    0,    0,    0,    0,    0,
0,    2,    0,    0,    0,    0,    0,    0,    0,    7,
0,    0,    0,    0,    0,    0,   52,   53,   32,    0,
0,   43,   17,    0,    0,    0,    0,    0,    0,    0,
0,    0,    0,    5,    3,    6,    9,    0,    0,    0,
0,    0,    0,    8,    0,    4,    0,    0,    0,    0,
0,    0,    0,    0,    0,    0,   55,    0,    0,    0,
56,    0,   51,   31,    0,    0,    0,    0,    0,    0,
0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
0,    0,    0,    0,    0,    0,   41,    0,    0,    0,
0,   57,    0,   58,    0,    0,   40,    0,    0,    0,
0,    0,   61,   63,   64,   42,    0,   34,    0,    0,
59,    0,   62,    0,   35,   65,    0,   60,   36,   66,
};
const short yydgoto[] = { 1,
69,   32,   16,   33,   59,   19,  118,   71,   27,
};
const short yysindex[] = { 0,
-10,   -1,    0,  -75,  -29,  -43,  -60,    6,    6,    6,
34,    0,   34,    5,   41,   50,   78,  469,   83,    0,
16,   34,   34,   34,   34,  -37,    0,    0,    0,    4,
-166,    0,    0,   64,   21,  368,   19,   34,   34,   34,
34,   34,   34,    0,    0,    0,    0,   34,   34,   34,
34,   34,   34,    0,   34,    0,   39,   49,  491,  102,
491,   60,  102,  491,  -92,   23,    0, -256,  -16,  478,
0,   34,    0,    0,   34,  -23,  228,  -23,  228,  -92,
-27,  -92,  -92,  -92,  -23,  228,  -23,  228,  -92,  -27,
-27,  -27,  -92,   82,   -7,   59,    0,   25,   31,  -21,
-28,    0,  -28,    0,   91,  -92,    0,   34,   34,   28,
33,  -28,    0,    0,    0,    0,   38,    0,  102,  -28,
0,   32,    0,   34,    0,    0,  -28,    0,    0,    0,
};
const short yyrindex[] = { 0,
0,    0,    0,    0,  445,   -4,  157,    0,    0,    0,
0,    0,    0,    0,  113,    0,  500,    0,    0,    0,
0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
169,    0,    0,  454,    0,    0,  534,    0,    0,    0,
0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
0,    0,    0,    0,    0,    0,  120,    0,    0,    0,
43,    0,  313,  107,  546,    0,    0,    0,    0,    0,
0,    0,    0,    0,    0,  125,   67,  147,  130,  338,
606,  392,  407,  433,  294,  319,  644,  525,  555,  615,
627,  636,  570,    0,    0,   12,    0,    0,    0,    0,
0,    0,    0,    0,    0,  591,    0,    0,    0,    0,
0,    0,    0,    0,    0,    0,    0,    0,  121,    0,
0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
};
const short yygindex[] = { 0,
348,  136,    0,  140,  292,    0,   22,  589,   35,
};
#define YYTABLESIZE 716
const short yytable[] = { 12,
68,   43,   13,   67,   99,   39,  100,   11,   20,   68,
23,   13,   55,   42,   44,   21,   11,   24,   40,  113,
42,   42,  112,   41,  102,   40,   38,  101,   39,   13,
41,   22,   39,   25,   11,   39,   39,   39,   39,   39,
39,   42,   39,   28,   29,   26,   40,   38,   42,   39,
45,   41,   10,   42,   42,   13,   42,   42,   42,   46,
11,   73,   40,   38,   13,   39,   53,   41,  121,   11,
43,  120,  128,   13,   42,  127,   24,   43,   11,   40,
38,  124,   39,   10,   41,   42,   10,   47,   39,   39,
40,   38,   56,   39,   72,   41,   42,   37,   43,   95,
97,   40,   38,   55,   39,   42,   41,   24,   57,   24,
24,   24,   75,   98,   43,  108,   11,  110,   42,  109,
39,  111,  107,   40,   38,  122,   39,   42,   41,   54,
33,   43,   40,   38,   44,   39,   15,   41,   42,   26,
17,   96,   43,   40,   38,  129,   39,   11,   41,   43,
11,    0,    0,   43,   43,   43,   45,   43,    0,   43,
0,    0,  125,    0,   44,   44,   13,   44,   44,   44,
26,   37,   26,   26,   26,   43,    0,    0,   50,    0,
0,    0,    0,  116,   43,    0,   45,   45,    0,   45,
45,   45,    0,   13,    0,   43,   13,   13,   13,   13,
13,   13,    0,   13,    0,   50,   43,    0,   50,   50,
50,   50,   50,   50,    0,   50,    0,   44,    0,    3,
0,    0,    0,   66,    5,    6,    7,    8,    3,   10,
0,    0,   66,    5,    6,    7,    8,    0,   10,   45,
37,    0,    0,    0,    0,    2,    3,   37,    0,   44,
4,    5,    6,    7,    8,    9,   10,    0,    0,   39,
0,   50,   50,    0,   52,    0,    0,   55,   37,   50,
0,   45,    3,    0,   51,   42,   30,    5,    6,    7,
8,    3,   10,    0,   37,   30,    5,    6,    7,    8,
3,   10,   18,   50,   30,    5,    6,    7,    8,    0,
10,   37,   34,   23,   36,    0,    0,    0,    0,    0,
0,    0,   37,   61,    0,   64,    0,   70,    0,    0,
0,   53,   37,   37,    0,    0,    0,    0,   18,   77,
79,   81,    0,   23,   23,    0,   23,   23,   23,   86,
88,   90,   91,   92,    0,   37,    0,   46,   14,    0,
0,    0,   37,   37,   37,    0,   37,    0,   31,   18,
35,   18,   18,   18,    0,   37,    0,    0,   58,   60,
62,   63,   65,    0,   46,    0,   43,   46,   46,   46,
46,   46,   46,    0,   46,   76,   78,   80,   82,   83,
84,    0,   70,    0,   70,   85,   87,   89,   60,   60,
93,   47,   94,   70,   52,   37,    0,   55,   74,   50,
48,   70,   49,    0,   51,    0,   48,    0,   70,  105,
0,    0,  106,    0,    0,    0,    0,    0,   47,    0,
46,   47,   47,   47,   47,   47,   47,   37,   47,    0,
0,    0,   49,   48,    0,  105,   48,   48,   48,   48,
48,   48,    0,   48,   12,  117,  119,    0,    0,    0,
0,   53,   46,   21,    0,    0,    0,    0,    0,   49,
0,  117,   49,   49,   49,   49,   49,   49,   54,   49,
0,   12,    0,    0,   47,   12,   12,   12,   12,   12,
21,   12,    0,    0,   21,   21,   21,   21,   21,   48,
21,    0,    0,    0,    0,   52,    0,    0,   55,    0,
50,   48,    0,   49,   52,   51,   47,   55,  104,   50,
48,  103,   49,    0,   51,   49,    0,   52,    0,    0,
55,   48,   50,   48,   22,   49,   17,   51,   12,   17,
0,   17,   17,   14,   17,    0,   17,   21,    0,    0,
0,    0,    0,    0,    0,   16,    0,   49,    0,    0,
0,    0,   53,    0,   28,   22,    0,   22,   22,   22,
14,   53,    0,   14,   14,   14,   14,   14,   14,   30,
14,    0,   16,    0,   53,   16,   16,   16,   16,   16,
16,   28,   16,   17,   28,   28,   28,   28,   28,   28,
15,   28,    0,    0,    0,    0,   30,    0,    0,   30,
30,   30,   30,   30,   30,   27,   30,    0,    0,    0,
0,    0,    0,    0,   29,    0,    0,   15,    0,    0,
15,   15,   15,   15,   15,   15,   20,   15,    0,    0,
0,    0,   27,    0,    0,   19,   27,   27,   27,   27,
27,   29,   27,   25,    0,   29,   29,   29,   29,   29,
0,   29,    0,   20,    0,    0,    0,   20,   20,   20,
20,   20,   19,   20,    0,    0,   19,   19,   19,   19,
19,    0,   19,   25,   25,    0,   25,   25,   25,  114,
0,  115,    0,    0,    0,    0,    0,    0,    0,    0,
123,    0,    0,    0,    0,    0,    0,    0,  126,    0,
0,    0,    0,    0,    0,  130,
};
const short yycheck[] = { 10,
38,   94,   40,   41,  261,   10,  263,   45,   10,   38,
40,   40,   40,   37,   10,   91,   45,   61,   42,   41,
37,   10,   44,   47,   41,   42,   43,   44,   45,   40,
47,   61,   37,   94,   45,   40,   41,   42,   43,   44,
45,   37,   47,    9,   10,   40,   42,   43,   37,   45,
10,   47,   10,   42,   43,   40,   45,   37,   47,   10,
45,   41,   42,   43,   40,   45,   94,   47,   41,   45,
94,   44,   41,   40,   37,   44,   10,   94,   45,   42,
43,   44,   45,   41,   47,   37,   44,   10,   93,   94,
42,   43,   10,   45,   91,   47,   37,  264,   94,   61,
41,   42,   43,   40,   45,   94,   47,   41,   93,   43,
44,   45,   94,   91,   94,  123,   10,   93,   37,   61,
125,   91,   41,   42,   43,   93,   45,   37,   47,   10,
10,   94,   42,   43,   10,   45,    1,   47,   37,   10,
1,   93,   94,   42,   43,  124,   45,   41,   47,   37,
44,   -1,   -1,   94,   42,   43,   10,   45,   -1,   47,
-1,   -1,  125,   -1,   40,   41,   10,   43,   44,   45,
41,  264,   43,   44,   45,   94,   -1,   -1,   10,   -1,
-1,   -1,   -1,   93,   94,   -1,   40,   41,   -1,   43,
44,   45,   -1,   37,   -1,   94,   40,   41,   42,   43,
44,   45,   -1,   47,   -1,   37,   94,   -1,   40,   41,
42,   43,   44,   45,   -1,   47,   -1,   93,   -1,  257,
-1,   -1,   -1,  261,  262,  263,  264,  265,  257,  267,
-1,   -1,  261,  262,  263,  264,  265,   -1,  267,   93,
264,   -1,   -1,   -1,   -1,  256,  257,  264,   -1,  125,
261,  262,  263,  264,  265,  266,  267,   -1,   -1,  264,
-1,   93,   94,   -1,   37,   -1,   -1,   40,  264,   42,
-1,  125,  257,   -1,   47,  264,  261,  262,  263,  264,
265,  257,  267,   -1,  264,  261,  262,  263,  264,  265,
257,  267,    1,  125,  261,  262,  263,  264,  265,   -1,
267,  264,   11,   10,   13,   -1,   -1,   -1,   -1,   -1,
-1,   -1,  264,   22,   -1,   24,   -1,   26,   -1,   -1,
-1,   94,   10,  264,   -1,   -1,   -1,   -1,   10,   38,
39,   40,   -1,   40,   41,   -1,   43,   44,   45,   48,
49,   50,   51,   52,   -1,  264,   -1,   10,    1,   -1,
-1,   -1,   40,   41,  264,   -1,   44,   -1,   11,   41,
13,   43,   44,   45,   -1,  264,   -1,   -1,   21,   22,
23,   24,   25,   -1,   37,   -1,  264,   40,   41,   42,
43,   44,   45,   -1,   47,   38,   39,   40,   41,   42,
43,   -1,  101,   -1,  103,   48,   49,   50,   51,   52,
53,   10,   55,  112,   37,   93,   -1,   40,   41,   42,
43,  120,   45,   -1,   47,   -1,   10,   -1,  127,   72,
-1,   -1,   75,   -1,   -1,   -1,   -1,   -1,   37,   -1,
93,   40,   41,   42,   43,   44,   45,  125,   47,   -1,
-1,   -1,   10,   37,   -1,   98,   40,   41,   42,   43,
44,   45,   -1,   47,   10,  108,  109,   -1,   -1,   -1,
-1,   94,  125,   10,   -1,   -1,   -1,   -1,   -1,   37,
-1,  124,   40,   41,   42,   43,   44,   45,   10,   47,
-1,   37,   -1,   -1,   93,   41,   42,   43,   44,   45,
37,   47,   -1,   -1,   41,   42,   43,   44,   45,   93,
47,   -1,   -1,   -1,   -1,   37,   -1,   -1,   40,   -1,
42,   43,   -1,   45,   37,   47,  125,   40,   41,   42,
43,   44,   45,   -1,   47,   93,   -1,   37,   -1,   -1,
40,  125,   42,   43,   10,   45,   37,   47,   94,   40,
-1,   42,   43,   10,   45,   -1,   47,   94,   -1,   -1,
-1,   -1,   -1,   -1,   -1,   10,   -1,  125,   -1,   -1,
-1,   -1,   94,   -1,   10,   41,   -1,   43,   44,   45,
37,   94,   -1,   40,   41,   42,   43,   44,   45,   10,
47,   -1,   37,   -1,   94,   40,   41,   42,   43,   44,
45,   37,   47,   94,   40,   41,   42,   43,   44,   45,
10,   47,   -1,   -1,   -1,   -1,   37,   -1,   -1,   40,
41,   42,   43,   44,   45,   10,   47,   -1,   -1,   -1,
-1,   -1,   -1,   -1,   10,   -1,   -1,   37,   -1,   -1,
40,   41,   42,   43,   44,   45,   10,   47,   -1,   -1,
-1,   -1,   37,   -1,   -1,   10,   41,   42,   43,   44,
45,   37,   47,   10,   -1,   41,   42,   43,   44,   45,
-1,   47,   -1,   37,   -1,   -1,   -1,   41,   42,   43,
44,   45,   37,   47,   -1,   -1,   41,   42,   43,   44,
45,   -1,   47,   40,   41,   -1,   43,   44,   45,  101,
-1,  103,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
112,   -1,   -1,   -1,   -1,   -1,   -1,   -1,  120,   -1,
-1,   -1,   -1,   -1,   -1,  127,
};
#define YYFINAL 1
#ifndef YYDEBUG
#define YYDEBUG 0
#endif
#define YYMAXTOKEN 269
#if YYDEBUG
const char * const yyname[] = {
	"end-of-file",0,0,0,0,0,0,0,0,0,"'\\n'",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,"'%'","'&'",0,"'('","')'","'*'","'+'","','","'-'",0,"'/'",0,0,0,0,0,
	0,0,0,0,0,0,0,0,"'='",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	"'['",0,"']'","'^'",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	"'{'",0,"'}'",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,"NUMBER","POL","ARGSTACK","ANARRAY","ARRAY",
	"POLYVAR","VAR","POLYTERM","BLTIN","BLTINV","BLTINP","UNDEF","UNARYMINUS",
};
const char * const yyrule[] = {
	"$accept : list",
	"list :",
	"list : list '\\n'",
	"list : list varassign '\\n'",
	"list : list arrayassign '\\n'",
	"list : list expr '\\n'",
	"list : list VOID '\\n'",
	"list : list error '\\n'",
	"list : list polyexpr '\\n'",
	"list : list polyassign '\\n'",
	"polyassign : POLYVAR '=' polyexpr",
	"polyassign : VAR '=' polyexpr",
	"polyexpr : POLYVAR",
	"polyexpr : POLYTERM",
	"polyexpr : expr POLYTERM",
	"polyexpr : expr POLYTERM '^' expr",
	"polyexpr : POLYTERM '^' expr",
	"polyexpr : polyassign",
	"polyexpr : polyexpr '+' polyexpr",
	"polyexpr : polyexpr '%' polyexpr",
	"polyexpr : polyexpr '/' polyexpr",
	"polyexpr : '-' polyexpr",
	"polyexpr : polyexpr '-' polyexpr",
	"polyexpr : polyexpr '+' expr",
	"polyexpr : expr '+' polyexpr",
	"polyexpr : polyexpr '-' expr",
	"polyexpr : expr '-' polyexpr",
	"polyexpr : expr '*' polyexpr",
	"polyexpr : polyexpr '*' expr",
	"polyexpr : polyexpr '*' polyexpr",
	"polyexpr : polyexpr '^' expr",
	"polyexpr : '(' polyexpr ')'",
	"polyexpr : BLTINP arglist",
	"arrayassign : ARRAY '[' expr ']' '=' expr",
	"arrayassign : ARRAY '[' ']' '=' '{' arraylist",
	"arraylist : expr '}'",
	"arraylist : expr ',' arraylist",
	"varassign : VAR '=' expr",
	"expr : NUMBER",
	"expr : VAR",
	"expr : polyexpr '(' expr ')'",
	"expr : POLYVAR '(' expr ')'",
	"expr : ARRAY '[' expr ']'",
	"expr : varassign",
	"expr : expr '+' expr",
	"expr : expr '-' expr",
	"expr : expr '*' expr",
	"expr : expr '/' expr",
	"expr : expr '%' expr",
	"expr : expr '^' expr",
	"expr : '-' expr",
	"expr : '(' expr ')'",
	"expr : BLTIN arglist",
	"VOID : BLTINV arglist",
	"VOID : ARRAY '[' ']'",
	"arglist : '(' ')'",
	"arglist : '(' arguments",
	"arguments : expr ')'",
	"arguments : polyexpr ')'",
	"arguments : ARRAY '[' ']' ')'",
	"arguments : '&' ARRAY '[' ']' ')'",
	"arguments : '&' VAR ')'",
	"arguments : '&' VAR ',' arguments",
	"arguments : expr ',' arguments",
	"arguments : polyexpr ',' arguments",
	"arguments : ARRAY '[' ']' ',' arguments",
	"arguments : '&' ARRAY '[' ']' ',' arguments",
};
#endif
#if YYDEBUG
#include <stdio.h>
#endif
#ifdef YYSTACKSIZE
#undef YYMAXDEPTH
#define YYMAXDEPTH YYSTACKSIZE
#else
#ifdef YYMAXDEPTH
#define YYSTACKSIZE YYMAXDEPTH
#else
#define YYSTACKSIZE 10000
#define YYMAXDEPTH 10000
#endif
#endif
#define YYINITSTACKSIZE 200
int yydebug;
int yynerrs;
int yyerrflag;
int yychar;
short *yyssp;
YYSTYPE *yyvsp;
YYSTYPE yyval;
YYSTYPE yylval;
short *yyss;
short *yysslim;
YYSTYPE *yyvs;
int yystacksize;
// #line 567 "parse.y"


static char *p;

static Stack argStack = NULL;


void Parse(char *s)

{
	p = s;
	argStack = stackNew();
	yyparse();
	switch (rettype) {
	case RET_MPI:
		FREEMPI(retval.mpi);
		break;
	case RET_POLY:
		DELETEPI(retval.polyi);
		break;

	}
	rettype = RET_NULL;
	stackFree(&argStack);


}

/* Checks if the arguments contained in the source stack match the arguments
* of the function whose information is stored in the supplied symbol.
* If the arguments match then the function returns a variable stack that
* can be sent straight to a wrapper function (see note at top of file)
* If the arguments do not match then an error message is printed, and the
* function returns NULL.
* Regardless of whether the match is successful or not, all the arguments on
* the globar variable argStack are popped, leaving an empty stack.
*/

Stack checkArgs(Symbol *s)
{
	Stack varStack = stackNew();
	Stack tmpStack = stackNew();
	int *argTypes = s->argTypes;
	int nArgs = argTypes[0];
	int quitFlag = 0;
	int i;
	Argument Arg = NULL;
	for (i = 1; ; i++) {
		if (i > nArgs && stackEmpty(argStack))
			break;
		if (!stackEmpty(argStack))
			Arg = stackPop(argStack);
		else
			quitFlag = 1; /* This means there are too few arguments sent to this
						  * function */


						  /* oh no! we quit if this is true */
		if (quitFlag == 1 || Arg->type != argTypes[i]) { /* short circuit eval! */
			int j;
			/* pop all entries off argument stack */
			if (!quitFlag) { /* only free these if the reason we are exiting is
							 * that there are too few Arguments */
				freeArg(Arg);
				while (!stackEmpty(argStack)) {
					Arg = stackPop(argStack);
					freeArg(Arg);
				}
			}
			/* pop all entries off stack to be returned
			* It is helpful that i contains a value one greater than the number
			* of elements on the stack to be returned. */
			while (!stackEmpty(tmpStack)) {
				Arg = stackPop(tmpStack);
				switch (Arg->type) {
				case NUM:
					FREEMPI(Arg->u.num);
					break;
				case POLY:
					DELETEPI(Arg->u.poly);
					break;
					/* we do not free arrays, array address or variable addresses */
				}
			}

			printf("This function has the format: %s(", s->name);
			if (nArgs > 0) {
				for (j = 1; j <= nArgs; j++) {
					switch (argTypes[j]) {
					case NUM:
						printf("number");
						break;
					case VARADR:
						printf("&var");
						break;
					case ARR:
						printf("array[]");
						break;
					case ARRADR:
						printf("&array[]");
						break;
					case POLY:
						printf("poly");
						break;
					}
					if (j != nArgs)
						printf(", ");
					else
						printf(")\n");
				}
			}
			else
				printf(")\n");

			stackFree(&tmpStack);
			return NULL;
		}
		stackPush(tmpStack, Arg); /* put it back in the stack. in reverse order */

	}

	while (!stackEmpty(tmpStack)) {
		Arg = stackPop(tmpStack);
		switch (Arg->type) {
		case NUM:
			stackPush(varStack, Arg->u.num);
			break;
		case VARADR:
			if (Arg->defined == 1) /* if the variable has previously been assigned.
								   */
				FREEMPI(*(Arg->u.varAdr));
			stackPush(varStack, Arg->u.varAdr);
			break;
		case ARR:
			stackPush(varStack, Arg->u.array);
			break;
		case ARRADR:
			if (Arg->defined == 1) /* if the variable has previously been assigned.
								   */
				FREEMPIA(*(Arg->u.arrayAdr));

			stackPush(varStack, Arg->u.arrayAdr);
			break;
		case POLY:
			stackPush(varStack, Arg->u.poly);
			break;
		}
	}
	stackFree(&tmpStack);
	return varStack;
}

int yylex()
{
	MPI *Temp;
	char c;
	int typ;

	while ((c = *p++) == ' ' || c == '\t')
		;
	if (c == '\0')
		return 0;
	if (isdigit((int)c))
	{
		yylval.val = CHANGE((USL)(c - '0'));
		while (isdigit((int)*p))
		{
			Temp = yylval.val;
			yylval.val = MULT_I(yylval.val, 10L);
			FREEMPI(Temp);
			Temp = yylval.val;
			yylval.val = ADD0_I(yylval.val, (USL)(*p++ - '0'));
			FREEMPI(Temp);
		}
		return NUMBER;
	}
	p--;
	if (isalpha((int)c))
	{
		Symbol *s;
		char sbuf[100], *tmp = sbuf;

		do {
			*tmp++ = *p;
			p++;
			c = (*p);
		} while (c && isalnum((int)c));
		*tmp = '\0';
		if (c == '[')
			typ = ARRAY;
		else if ((c == '^' && *(tmp - 1) == POLYID && (tmp - 1) == sbuf) ||
			(*(p - 1) == POLYID && strlen(sbuf) == 1)) {
			return POLYTERM; /* don't want to install this */
		}
		else
			typ = VAR;
		s = lookup(sbuf, typ);
		if (s == NULL)
			s = lookup(sbuf, BLTIN);
		if (s == NULL)
			s = lookup(sbuf, BLTINV);
		if (s == NULL)
			s = lookup(sbuf, BLTINP);
		if (s == NULL)
			s = lookup(sbuf, POLYVAR);
		if (s == NULL)
			s = install(sbuf, UNDEF);
		yylval.sym = s;
		if (s->type == UNDEF)
			return(typ);
		else
			return (s->type);

	}
	p++;
	return (int)c; /* returns +, -, ^, *, / etc */
}

void yyerror(s)
char *s;
{
	warning(s, "");
}
// #line 619 "y.tab.c"
/* allocate initial stack or double stack size, up to YYMAXDEPTH */
static int yygrowstack()
{
	int newsize;
	ptrdiff_t i;
	short *newss;
	YYSTYPE *newvs;

	if ((newsize = yystacksize) == 0)
		newsize = YYINITSTACKSIZE;
	else if (newsize >= YYMAXDEPTH)
		return -1;
	else if ((newsize *= 2) > YYMAXDEPTH)
		newsize = YYMAXDEPTH;
	i = yyssp - yyss;
	newss = yyss ? (short *)realloc(yyss, newsize * sizeof *newss) :
		(short *)malloc(newsize * sizeof *newss);
	if (newss == NULL)
		return -1;
	yyss = newss;
	yyssp = newss + i;
	newvs = yyvs ? (YYSTYPE *)realloc(yyvs, newsize * sizeof *newvs) :
		(YYSTYPE *)malloc(newsize * sizeof *newvs);
	if (newvs == NULL)
		return -1;
	yyvs = newvs;
	yyvsp = newvs + i;
	yystacksize = newsize;
	yysslim = yyss + newsize - 1;
	return 0;
}

#define YYABORT goto yyabort
#define YYREJECT goto yyabort
#define YYACCEPT goto yyaccept
#define YYERROR goto yyerrlab

#ifndef YYPARSE_PARAM
#if defined(__cplusplus) || __STDC__
#define YYPARSE_PARAM_ARG void
#define YYPARSE_PARAM_DECL
#else	/* ! ANSI-C/C++ */
#define YYPARSE_PARAM_ARG
#define YYPARSE_PARAM_DECL
#endif	/* ANSI-C/C++ */
#else	/* YYPARSE_PARAM */
#ifndef YYPARSE_PARAM_TYPE
#define YYPARSE_PARAM_TYPE void *
#endif
#if defined(__cplusplus) || __STDC__
#define YYPARSE_PARAM_ARG YYPARSE_PARAM_TYPE YYPARSE_PARAM
#define YYPARSE_PARAM_DECL
#else	/* ! ANSI-C/C++ */
#define YYPARSE_PARAM_ARG YYPARSE_PARAM
#define YYPARSE_PARAM_DECL YYPARSE_PARAM_TYPE YYPARSE_PARAM;
#endif	/* ANSI-C/C++ */
#endif	/* ! YYPARSE_PARAM */

int
yyparse(YYPARSE_PARAM_ARG)
YYPARSE_PARAM_DECL
{
	register int yym, yyn, yystate;
#if YYDEBUG
	register const char *yys;

	if ((yys = getenv("YYDEBUG")))
	{
		yyn = *yys;
		if (yyn >= '0' && yyn <= '9')
			yydebug = yyn - '0';
	}
#endif

	yynerrs = 0;
	yyerrflag = 0;
	yychar = (-1);

	if (yyss == NULL && yygrowstack()) goto yyoverflow;
	yyssp = yyss;
	yyvsp = yyvs;
	*yyssp = yystate = 0;

yyloop:
	if ((yyn = yydefred[yystate])) goto yyreduce;
	if (yychar < 0)
	{
		if ((yychar = yylex()) < 0) yychar = 0;
#if YYDEBUG
		if (yydebug)
		{
			yys = 0;
			if (yychar <= YYMAXTOKEN) yys = yyname[yychar];
			if (!yys) yys = "illegal-symbol";
			printf("%sdebug: state %d, reading %d (%s)\n",
				YYPREFIX, yystate, yychar, yys);
		}
#endif
	}
	if ((yyn = yysindex[yystate]) && (yyn += yychar) >= 0 &&
		yyn <= YYTABLESIZE && yycheck[yyn] == yychar)
	{
#if YYDEBUG
		if (yydebug)
			printf("%sdebug: state %d, shifting to state %d\n",
				YYPREFIX, yystate, yytable[yyn]);
#endif
		if (yyssp >= yysslim && yygrowstack())
		{
			goto yyoverflow;
		}
		*++yyssp = yystate = yytable[yyn];
		*++yyvsp = yylval;
		yychar = (-1);
		if (yyerrflag > 0)  --yyerrflag;
		goto yyloop;
	}
	if ((yyn = yyrindex[yystate]) && (yyn += yychar) >= 0 &&
		yyn <= YYTABLESIZE && yycheck[yyn] == yychar)
	{
		yyn = yytable[yyn];
		goto yyreduce;
	}
	if (yyerrflag) goto yyinrecovery;
#if defined(lint) || defined(__GNUC__)
	goto yynewerror;
#endif
yynewerror:
	yyerror("syntax error");
#if defined(lint) || defined(__GNUC__)
	goto yyerrlab;
#endif
yyerrlab:
	++yynerrs;
yyinrecovery:
	if (yyerrflag < 3)
	{
		yyerrflag = 3;
		for (;;)
		{
			if ((yyn = yysindex[*yyssp]) && (yyn += YYERRCODE) >= 0 &&
				yyn <= YYTABLESIZE && yycheck[yyn] == YYERRCODE)
			{
#if YYDEBUG
				if (yydebug)
					printf("%sdebug: state %d, error recovery shifting\
 to state %d\n", YYPREFIX, *yyssp, yytable[yyn]);
#endif
				if (yyssp >= yysslim && yygrowstack())
				{
					goto yyoverflow;
				}
				*++yyssp = yystate = yytable[yyn];
				*++yyvsp = yylval;
				goto yyloop;
			}
			else
			{
#if YYDEBUG
				if (yydebug)
					printf("%sdebug: error recovery discarding state %d\n",
						YYPREFIX, *yyssp);
#endif
				if (yyssp <= yyss) goto yyabort;
				--yyssp;
				--yyvsp;
			}
		}
	}
	else
	{
		if (yychar == 0) goto yyabort;
#if YYDEBUG
		if (yydebug)
		{
			yys = 0;
			if (yychar <= YYMAXTOKEN) yys = yyname[yychar];
			if (!yys) yys = "illegal-symbol";
			printf("%sdebug: state %d, error recovery discards token %d (%s)\n",
				YYPREFIX, yystate, yychar, yys);
		}
#endif
		yychar = (-1);
		goto yyloop;
	}
yyreduce:
#if YYDEBUG
	if (yydebug)
		printf("%sdebug: state %d, reducing by rule %d (%s)\n",
			YYPREFIX, yystate, yyn, yyrule[yyn]);
#endif
	yym = yylen[yyn];
	yyval = yyvsp[1 - yym];
	switch (yyn)
	{
	case 3:
// #line 64 "parse.y"
	{
		retval.mpi = yyvsp[-1].val;
		rettype = RET_MPI;
	}
	break;
	case 4:
// #line 69 "parse.y"
	{
		rettype = RET_MPIA;
	}
	break;
	case 5:
// #line 73 "parse.y"
	{

		retval.mpi = yyvsp[-1].val;
		if (!(rettype == FUNC_FAIL)) {
			PRINTI(yyvsp[-1].val);
			printf("\n");
		}
		rettype = RET_MPI;

	}
	break;
	case 6:
// #line 84 "parse.y"
	{
		rettype = 0;
	}
	break;
	case 7:
// #line 88 "parse.y"
	{
		yyerrok;
	}
	break;
	case 8:
// #line 92 "parse.y"
	{
		retval.polyi = yyvsp[-1].pol;
		if (!(rettype == FUNC_FAIL)) {
			PRINTPI(yyvsp[-1].pol);
			printf("\n");
		}
		rettype = RET_POLY;
	}
	break;
	case 9:
// #line 101 "parse.y"
	{
		retval.polyi = yyvsp[-1].pol;
		rettype = RET_POLY;
	}
	break;
	case 10:
// #line 108 "parse.y"
	{
		DELETEPI(yyvsp[-2].sym->u.sympval);
		yyvsp[-2].sym->u.sympval = yyvsp[0].pol;
		yyval.pol = COPYPI(yyvsp[0].pol);
	}
	break;
	case 11:
// #line 114 "parse.y"
	{
		if (yyvsp[-2].sym->type == VAR)
			FREEMPI(yyvsp[-2].sym->u.symval);
		yyvsp[-2].sym->u.sympval = yyvsp[0].pol;
		yyvsp[-2].sym->type = POLYVAR;
		yyval.pol = COPYPI(yyvsp[0].pol);
	}
	break;
	case 12:
// #line 123 "parse.y"
	{
		if (yyvsp[0].sym->type == UNDEF)
			execerror("is an undefined variable", yyvsp[0].sym->name);
		else if (yyvsp[0].sym->u.sympval != (POLYI)NULL)
			yyval.pol = COPYPI(yyvsp[0].sym->u.sympval);
		else
			yyval.pol = (POLYI)NULL;
	}
	break;
	case 13:
// #line 132 "parse.y"
	{
		MPI *ONE;
		ONE = ONEI();
		yyval.pol = NULL;
		PINSERTPI(1, ONE, &yyval.pol, 1);
		FREEMPI(ONE);
	}
	break;
	case 14:
// #line 140 "parse.y"
	{
		yyval.pol = NULL;
		PINSERTPI(1, yyvsp[-1].val, &yyval.pol, 1);
		FREEMPI(yyvsp[-1].val);
	}
	break;
	case 15:
// #line 146 "parse.y"
	{
		yyval.pol = NULL;
		PINSERTPI(CONVERTI(yyvsp[0].val), yyvsp[-3].val, &yyval.pol, 1);
		FREEMPI(yyvsp[-3].val);
		FREEMPI(yyvsp[0].val);

	}
	break;
	case 16:
// #line 154 "parse.y"
	{
		MPI *O;
		O = ONEI();
		yyval.pol = NULL;
		PINSERTPI(CONVERTI(yyvsp[0].val), O, &yyval.pol, 1);
		FREEMPI(O);
		FREEMPI(yyvsp[0].val);
	}
	break;
	case 17:
// #line 163 "parse.y"
	{
	}
	break;
	case 18:
// #line 166 "parse.y"
	{
		yyval.pol = ADDPI(yyvsp[-2].pol, yyvsp[0].pol);
		DELETEPI(yyvsp[-2].pol);
		DELETEPI(yyvsp[0].pol);
	}
	break;
	case 19:
// #line 173 "parse.y"
	{
		yyval.pol = MODPI(yyvsp[-2].pol, yyvsp[0].pol);
		DELETEPI(yyvsp[-2].pol);
		DELETEPI(yyvsp[0].pol);

	}
	break;
	case 20:
// #line 180 "parse.y"
	{
		yyval.pol = DIVPI(yyvsp[-2].pol, yyvsp[0].pol);
		DELETEPI(yyvsp[-2].pol);
		DELETEPI(yyvsp[0].pol);
	}
	break;
	case 21:
// #line 186 "parse.y"
	{
		POLYI Z;
		Z = ZEROPI(); /* returns zero polynomial */
		yyval.pol = SUBPI(Z, yyvsp[0].pol);
		DELETEPI(Z);
		DELETEPI(yyvsp[0].pol);
	}
	break;
	case 22:
// #line 194 "parse.y"
	{
		yyval.pol = SUBPI(yyvsp[-2].pol, yyvsp[0].pol);
		DELETEPI(yyvsp[-2].pol);
		DELETEPI(yyvsp[0].pol);
	}
	break;
	case 23:
// #line 200 "parse.y"
	{
		POLYI O, P;
		O = ONEPI();
		P = SCALARPI(yyvsp[0].val, O);
		yyval.pol = ADDPI(yyvsp[-2].pol, P);
		DELETEPI(O);
		DELETEPI(P);
		DELETEPI(yyvsp[-2].pol);
		FREEMPI(yyvsp[0].val);
	}
	break;
	case 24:
// #line 211 "parse.y"
	{
		POLYI O, P;
		O = ONEPI();
		P = SCALARPI(yyvsp[-2].val, O);
		yyval.pol = ADDPI(yyvsp[0].pol, P);
		DELETEPI(O);
		DELETEPI(P);
		DELETEPI(yyvsp[0].pol);
		FREEMPI(yyvsp[-2].val);
	}
	break;
	case 25:
// #line 223 "parse.y"
	{
		POLYI O, P;
		O = ONEPI();
		P = SCALARPI(yyvsp[0].val, O);
		yyval.pol = SUBPI(yyvsp[-2].pol, P);
		DELETEPI(O);
		DELETEPI(P);
		DELETEPI(yyvsp[-2].pol);
		FREEMPI(yyvsp[0].val);
	}
	break;
	case 26:
// #line 234 "parse.y"
	{
		POLYI O, P;
		O = ONEPI();
		P = SCALARPI(yyvsp[-2].val, O);
		yyval.pol = SUBPI(P, yyvsp[0].pol);
		DELETEPI(O);
		DELETEPI(P);
		DELETEPI(yyvsp[0].pol);
		FREEMPI(yyvsp[-2].val);
	}
	break;
	case 27:
// #line 245 "parse.y"
	{
		yyval.pol = SCALARPI(yyvsp[-2].val, yyvsp[0].pol);
		FREEMPI(yyvsp[-2].val);
		DELETEPI(yyvsp[0].pol);
	}
	break;
	case 28:
// #line 251 "parse.y"
	{
		yyval.pol = SCALARPI(yyvsp[0].val, yyvsp[-2].pol);
		FREEMPI(yyvsp[0].val);
		DELETEPI(yyvsp[-2].pol);
	}
	break;
	case 29:
// #line 257 "parse.y"
	{
		yyval.pol = MULTPI(yyvsp[-2].pol, yyvsp[0].pol);
		DELETEPI(yyvsp[-2].pol);
		DELETEPI(yyvsp[0].pol);
	}
	break;
	case 30:
// #line 263 "parse.y"
	{
		yyval.pol = POWERPI(yyvsp[-2].pol, CONVERTI(yyvsp[0].val));
		DELETEPI(yyvsp[-2].pol);
		FREEMPI(yyvsp[0].val);
	}
	break;
	case 31:
// #line 269 "parse.y"
	{
		yyval.pol = yyvsp[-1].pol;
	}
	break;
	case 32:
// #line 273 "parse.y"
	{
		Stack varStack;
		if ((varStack = checkArgs(yyvsp[-1].sym))) {
			yyval.pol = (*(yyvsp[-1].sym->u.ptrp))(varStack);
			stackFree(&varStack);
		}
		else {
			rettype = FUNC_FAIL; /* return failure of function */
			yyval.pol = NULL;
		}
	}
	break;
	case 33:
// #line 290 "parse.y"
	{
		unsigned n = CONVERTI(yyvsp[-3].val);
		if (yyvsp[-5].sym->type == UNDEF) {
			yyvsp[-5].sym->u.symarr = BUILDMPIA();
			yyvsp[-5].sym->type = ARRAY;
		}
		ADD_TO_MPIA(yyvsp[-5].sym->u.symarr, yyvsp[0].val, n);
		FREEMPI(yyvsp[-3].val);
		FREEMPI(yyvsp[0].val);

	}
	break;
	case 34:
// #line 302 "parse.y"
	{
		if (yyvsp[-5].sym->type != UNDEF)
			FREEMPIA(yyvsp[-5].sym->u.symarr);
		yyvsp[-5].sym->u.symarr = yyvsp[0].arr;
		yyvsp[-5].sym->type = ARRAY;

	}
	break;
	case 35:
// #line 310 "parse.y"
	{
		yyval.arr = BUILDMPIA();
		ADD_TO_MPIA(yyval.arr, yyvsp[-1].val, 0);
		FREEMPI(yyvsp[-1].val);
	}
	break;
	case 36:
// #line 316 "parse.y"
	{
		MPIA_INSERT(yyvsp[0].arr, yyvsp[-2].val, 0);
		FREEMPI(yyvsp[-2].val);
		yyval.arr = yyvsp[0].arr;
	}
	break;
	case 37:
// #line 323 "parse.y"
	{
		if (yyvsp[-2].sym->type != UNDEF && yyvsp[-2].sym->u.symval != NULL)
			FREEMPI(yyvsp[-2].sym->u.symval);
		yyvsp[-2].sym->u.symval = yyvsp[0].val;
		if (yyvsp[-2].sym->u.symval != (MPI *)NULL)
			yyval.val = COPYI(yyvsp[-2].sym->u.symval);
		else
			yyval.val = (MPI *)NULL;
		yyvsp[-2].sym->type = VAR;
	}
	break;
	case 38:
// #line 334 "parse.y"
	{
		yyval.val = yyvsp[0].val;
	}
	break;
	case 39:
// #line 338 "parse.y"
	{
		if (yyvsp[0].sym->type == UNDEF)
			execerror("is an undefined variable", yyvsp[0].sym->name);
		else
			yyval.val = COPYI(yyvsp[0].sym->u.symval);
	}
	break;
	case 40:
// #line 350 "parse.y"
	{
		yyval.val = VALPI(yyvsp[-3].pol, yyvsp[-1].val);
		FREEMPI(yyvsp[-1].val);
		DELETEPI(yyvsp[-3].pol);
	}
	break;
	case 41:
// #line 356 "parse.y"
	{
		yyval.val = VALPI(yyvsp[-3].sym->u.sympval, yyvsp[-1].val);
		FREEMPI(yyvsp[-1].val);
	}
	break;
	case 42:
// #line 361 "parse.y"
	{
		unsigned long ind;
		if (yyvsp[-3].sym->type == UNDEF)
		{
			FREEMPI(yyvsp[-1].val);
			execerror("[] is an undefined array", yyvsp[-3].sym->name);
		}
		ind = CONVERTI(yyvsp[-1].val);
		if (ind >= yyvsp[-3].sym->u.symarr->size)
		{
			FREEMPI(yyvsp[-1].val);
			execerror("array is too small", yyvsp[-3].sym->name);
		}
		yyval.val = COPYI(yyvsp[-3].sym->u.symarr->A[ind]);
		FREEMPI(yyvsp[-1].val);
	}
	break;
	case 43:
// #line 378 "parse.y"
	{
	}
	break;
	case 44:
// #line 381 "parse.y"
	{
		yyval.val = ADDI(yyvsp[-2].val, yyvsp[0].val);
		FREEMPI(yyvsp[-2].val);
		FREEMPI(yyvsp[0].val);
	}
	break;
	case 45:
// #line 387 "parse.y"
	{
		yyval.val = SUBI(yyvsp[-2].val, yyvsp[0].val);
		FREEMPI(yyvsp[-2].val);
		FREEMPI(yyvsp[0].val);
	}
	break;
	case 46:
// #line 393 "parse.y"
	{
		yyval.val = MULTI(yyvsp[-2].val, yyvsp[0].val);
		FREEMPI(yyvsp[-2].val);
		FREEMPI(yyvsp[0].val);
	}
	break;
	case 47:
// #line 399 "parse.y"
	{
		if ((yyvsp[0].val)->S <= 0)
		{
			FREEMPI(yyvsp[-2].val);
			FREEMPI(yyvsp[0].val);
			execerror(" divisor <= 0", "");
		}
		yyval.val = INTI(yyvsp[-2].val, yyvsp[0].val);
		FREEMPI(yyvsp[-2].val);
		FREEMPI(yyvsp[0].val);
	}
	break;
	case 48:
// #line 411 "parse.y"
	{
		if ((yyvsp[0].val)->S <= 0)
		{
			FREEMPI(yyvsp[-2].val);
			FREEMPI(yyvsp[0].val);
			execerror(" divisor <= 0", "");
		}
		yyval.val = MOD(yyvsp[-2].val, yyvsp[0].val);
		FREEMPI(yyvsp[-2].val);
		FREEMPI(yyvsp[0].val);
	}
	break;
	case 49:
// #line 423 "parse.y"
	{
		if ((yyvsp[0].val)->S < 0)
		{
			FREEMPI(yyvsp[-2].val);
			FREEMPI(yyvsp[0].val);
			execerror("negative exponent", "");
		}
		if ((yyvsp[0].val)->D > 0)
		{
			FREEMPI(yyvsp[-2].val);
			FREEMPI(yyvsp[0].val);
			execerror("exponent >= R0", "");
		}
		yyval.val = POWERI(yyvsp[-2].val, (unsigned int)(CONVERTI(yyvsp[0].val)));
		FREEMPI(yyvsp[-2].val);
		FREEMPI(yyvsp[0].val);
	}
	break;
	case 50:
// #line 441 "parse.y"
	{
		yyval.val = MINUSI(yyvsp[0].val);
		FREEMPI(yyvsp[0].val);
	}
	break;
	case 51:
// #line 446 "parse.y"
	{
		yyval.val = yyvsp[-1].val;
	}
	break;
	case 52:
// #line 450 "parse.y"
	{
		Stack varStack;
		if ((varStack = checkArgs(yyvsp[-1].sym))) {
			yyval.val = (*(yyvsp[-1].sym->u.ptr))(varStack);
			stackFree(&varStack);
		}
		else {
			yyval.val = NULL;
			rettype = FUNC_FAIL;
		}
	}
	break;
	case 53:
// #line 463 "parse.y"
	{
		Stack varStack;
		if ((varStack = checkArgs(yyvsp[-1].sym))) {
			(*(yyvsp[-1].sym->u.ptrv))(varStack);
			stackFree(&varStack);
		}
		else {
			rettype = FUNC_FAIL;
		}
	}
	break;
	case 54:
// #line 473 "parse.y"
	{
		if (yyvsp[-2].sym->type != UNDEF) {
			PRINTIA(yyvsp[-2].sym->u.symarr);
		}
		else {
			execerror("undefined array", yyvsp[-2].sym->name);
		}
	}
	break;
	case 55:
// #line 482 "parse.y"
	{
	}
	break;
	case 56:
// #line 485 "parse.y"
	{
	}
	break;
	case 57:
// #line 490 "parse.y"
	{
		stackPush(argStack, createArg(yyvsp[-1].val, NUM, 0));
	}
	break;
	case 58:
// #line 494 "parse.y"
	{
		stackPush(argStack, createArg(yyvsp[-1].pol, POLY, 0));
	}
	break;
	case 59:
// #line 498 "parse.y"
	{
		if (yyvsp[-3].sym->type != UNDEF) {
			stackPush(argStack, createArg(yyvsp[-3].sym, ARR, 0));
		}
		else {
			execerror("[] is an undefined array", yyvsp[-3].sym->name);
			rettype = FUNC_FAIL;
		}
	}
	break;
	case 60:
// #line 507 "parse.y"
	{

		/* S.Seefried.  Normally one would free an old array if it
		* existed.  I have deferred this until I am sure all arguments
		* on the command line are valid for the function in question.
		* See checkArgs.  It is here that I free. */
		if (yyvsp[-3].sym->type != UNDEF)
			stackPush(argStack, createArg(yyvsp[-3].sym, ARRADR, 1));
		else
			stackPush(argStack, createArg(yyvsp[-3].sym, ARRADR, 0));
		yyvsp[-3].sym->type = ARRAY;
	}
	break;
	case 61:
// #line 520 "parse.y"
	{

		if (yyvsp[-1].sym->type != UNDEF)
			stackPush(argStack, createArg(yyvsp[-1].sym, VARADR, 1));
		else
			stackPush(argStack, createArg(yyvsp[-1].sym, VARADR, 0));
		yyvsp[-1].sym->type = VAR;
	}
	break;
	case 62:
// #line 530 "parse.y"
	{
		if (yyvsp[-2].sym->type != UNDEF)
			stackPush(argStack, createArg(yyvsp[-2].sym, VARADR, 1));
		else
			stackPush(argStack, createArg(yyvsp[-2].sym, VARADR, 0));
		yyvsp[-2].sym->type = VAR;
	}
	break;
	case 63:
// #line 538 "parse.y"
	{
		stackPush(argStack, createArg(yyvsp[-2].val, NUM, 0));
	}
	break;
	case 64:
// #line 542 "parse.y"
	{
		stackPush(argStack, createArg(yyvsp[-2].pol, POLY, 0));
	}
	break;
	case 65:
// #line 546 "parse.y"
	{
		if (yyvsp[-4].sym->type != UNDEF) {
			stackPush(argStack, createArg(yyvsp[-4].sym, ARR, 0));
		}
		else {
			rettype = FUNC_FAIL;
		}
	}
	break;
	case 66:
// #line 554 "parse.y"
	{
		/* S.Seefried.  Normally one would free an old array if it
		* existed.  I have deferred this until I am sure all arguments
		* on the command line are valid for the function in question.
		* See checkArgs.  It is here that I free. */

		if (yyvsp[-4].sym->type != UNDEF)
			stackPush(argStack, createArg(yyvsp[-4].sym, ARRADR, 1));
		else
			stackPush(argStack, createArg(yyvsp[-4].sym, ARRADR, 0));
		yyvsp[-4].sym->type = ARRAY;
	}
	break;
// #line 1422 "y.tab.c"
	}
	yyssp -= yym;
	yystate = *yyssp;
	yyvsp -= yym;
	yym = yylhs[yyn];
	if (yystate == 0 && yym == 0)
	{
#if YYDEBUG
		if (yydebug)
			printf("%sdebug: after reduction, shifting from state 0 to\
 state %d\n", YYPREFIX, YYFINAL);
#endif
		yystate = YYFINAL;
		*++yyssp = YYFINAL;
		*++yyvsp = yyval;
		if (yychar < 0)
		{
			if ((yychar = yylex()) < 0) yychar = 0;
#if YYDEBUG
			if (yydebug)
			{
				yys = 0;
				if (yychar <= YYMAXTOKEN) yys = yyname[yychar];
				if (!yys) yys = "illegal-symbol";
				printf("%sdebug: state %d, reading %d (%s)\n",
					YYPREFIX, YYFINAL, yychar, yys);
			}
#endif
		}
		if (yychar == 0) goto yyaccept;
		goto yyloop;
	}
	if ((yyn = yygindex[yym]) && (yyn += yystate) >= 0 &&
		yyn <= YYTABLESIZE && yycheck[yyn] == yystate)
		yystate = yytable[yyn];
	else
		yystate = yydgoto[yym];
#if YYDEBUG
	if (yydebug)
		printf("%sdebug: after reduction, shifting from state %d \
to state %d\n", YYPREFIX, *yyssp, yystate);
#endif
	if (yyssp >= yysslim && yygrowstack())
	{
		goto yyoverflow;
	}
	*++yyssp = yystate;
	*++yyvsp = yyval;
	goto yyloop;
yyoverflow:
	yyerror("yacc stack overflow");
yyabort:
	return (1);
yyaccept:
	return (0);
}
