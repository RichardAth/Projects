#pragma once
extern int lang;

void fcubesText(char *input, int groupLen);
void fsquaresText(char *input, int groupLen);
void contfracText(char *input, int groupLen);
void dilogText(char *baseText, char *powerText, char *modText, int groupLen);
void gaussianText(char *valueText, int doFactorization);
void ecmFrontText(char *tofactorText, int doFactorization, char *knownFactors, char *output);
//void polyFactText(char *modText, char *polyText, int groupLen);
extern void quadmodText(char *quadrText, char *linearText, char *constText,
	char *modText, int groupLength);