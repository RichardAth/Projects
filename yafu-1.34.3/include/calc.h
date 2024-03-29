#pragma once
/*----------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Ben Buhrow. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	

Some parts of the code (and also this header), included in this 
distribution have been reused from other sources. In particular I 
have benefitted greatly from the work of Jason Papadopoulos's msieve @ 
www.boo.net/~jasonp, Scott Contini's mpqs implementation, and Tom St. 
Denis Tom's Fast Math library.  Many thanks to their kind donation of 
code to the public domain.
       				   --bbuhrow@gmail.com 11/24/09
----------------------------------------------------------------------*/

#include "yafu.h"
#include "arith.h"
#include "factor.h"

//symbols in calc
#define EOE 1
#define IMM 2
#define NUM 3
#define OP 4
#define RP 5
#define LP 6
#define CH 7
#define AMBIG 8
#define COMMA 9
#define SPACE 10

//operator associativity
#define RIGHT 1
#define LEFT 0

#define NUM_FUNC 67

//arbitrary precision calculator
//void testcalc(void);
//void handle_singleop(char *arg1, int op);
//int single_op(char s);
//int dual_op(char s);
int preprocess(str_t *str);
//int get_el_type(char s);
//int processIMM(int opcode, str_t *str);
int calc(str_t *str, fact_obj_t *fobj);
//int isInt(char s);

//static int getAssoc(const char *s);
//int processOP(char *s, str_t *n1, str_t *n2);

int getFunc(const char *s, int *nargs);
//int getArgs(str_t *args, str_t *in, int num);
int new_uvar(const char *name, mpz_t data);
int set_uvar(const char *name, const mpz_t data);
int get_uvar(const char *name, mpz_t data);
static void free_uvars();
int invalid_dest(const char *dest);
//int invalid_num(char *num);
//int calc2(str_t *in);
static char** tokenize(const char *in, int *token_types, int *num_tokens);
//static int get_el_type2(char s);
static int is_new_token(int el_type, int el_type2);
void calc_finalize();
int calc_init();

//user variables
uvars_t uvars;
