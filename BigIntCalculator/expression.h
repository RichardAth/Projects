#pragma once

#define COPYRIGHT_SPANISH "Hecho por Darío Alpern. Actualizado el 17 de setiembre de 2022."
#define COPYRIGHT_ENGLISH "Written by Dario Alpern. Last updated on 17 September 2022."

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
