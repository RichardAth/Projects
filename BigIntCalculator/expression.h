#pragma once

#define COPYRIGHT_SPANISH "Hecho por Darío Alpern. Actualizado el 18 de febrero de 2018."
#define COPYRIGHT_ENGLISH "Written by Dario Alpern. Last updated on 18 February 2018."

#ifdef __EMSCRIPTEN__
int stamp(void);
#endif
void databack(char *data);

extern char inputString[1000000];
extern char output[300000];
extern BigInteger valueX;
extern int counterC;
extern int expressionNbr;


#define FACTORIZATION_FUNCTIONS 1
