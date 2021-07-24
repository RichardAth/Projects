/* calculators/test.c
@alpertron alpertron Fixed calculation when adding new relation on SIQS algorithm
242 lines (240 sloc) 6.52 KB */
#define _CRT_SECURE_NO_DEPRECATE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <locale.h>
#include <time.h>
#include <ctype.h>
#include "bignbr.h"
#include "highlevel.h"
#include "factor.h"
#include "batch.h"
#include "expression.h"


int factors[5000];
//struct sFactors astFactors[1000];
clock_t beginTime;
extern int number[MAX_LEN];
extern int nbrLimbs;
extern int lang, groupLen;
extern limb TestNbr[MAX_LEN];



static limb Product[32];
char input[10000];
extern char tofactorDec[30000];

static BigInteger dividend, divisor, quotient;
static char helpmessage[] = 
"0 \texit\n"
"1 \tfsquaresText   (user supplied input)\n"
"2 \tfcubesText     (arg1, 6)\n"
"3 \tmultiply & squareRoot\n"
"4 \tBigIntDivide\n"
"5 \tFcubesText   (series of values)\n"
"6 \tcontfracText(arg1 arg2 arg3)\n"
"7 \tfsquaresText (arg1, 6)\n"
"8 \tModInvBigNbr\n"
"10\tEmcFrontText (user supplied input)\n"
"11	\tDilogText\n"
"12 \tgaussianText\n"
"13 \temcFrontText (predefined values)\n"
"14 \tDec2Bin & Bin2Dec\n"
"15 \tdoWork\n"
"16 \tquadmodText\n"
"17 \tChange options\n";

static char ayuda[] =
"0 \tsalida\n"
"1 \tfsquaresText   (entrada proporcionada por el usuario)\n"
"2 \tfcubesText     (arg1, 6)\n"
"3 \tmultiply & squareRoot\n"
"4 \tBigIntDivide\n"
"5 \tFcubesText   (serie de valores)\n"
"6 \tcontfracText(arg1 arg2 arg3)\n"
"7 \tfsquaresText (arg1, 6)\n"
"8 \tModInvBigNbr\n"
"10\tEmcFrontText (entrada proporcionada por el usuario)\n"
"11	\tDilogText\n"
"12 \tgaussianText\n"
"13 \temcFrontText (valores predefinidos)\n"
"14 \tDec2Bin & Bin2Dec\n"
"15 \tdoWork\n"
"16 \tquadmodText\n"
"17 \tCambiar opciones\n";

extern void doWork(void);

#ifdef _DEBUG
	bool debug = true;
#else
	bool debug = false;
#endif

int main(int argc, char *argv[])
{
	limb Factor1[] = { 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x00, 0x00, 0x00 };
	limb Factor2[] = { 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x00, 0x00, 0x00 };
	limb Factor3[] = { 29504, 29490, 19798, 633, 181, 0, 0, 0, 0, 0 };
	limb Factor4[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
	limb Factor5[] = { 32767, 32767, 32767, 32767, 0, 0, 0, 0 };
	limb Factor6[] = { 32767, 32767, 32767, 32767, 0, 0, 0, 0 };
	char expr[] = "123456789012345";
	int len, i;
	char ibuff[80] = { 0 };
	char *ptrInput;
	int debug_code=1;
	int performFactorization = true;

	setlocale(LC_ALL, "en-GB");	
	char banner[] = "Compiled on "  __DATE__ " at " __TIME__ "\n";
	printf("%s", banner);
	while (debug_code != 0) {
		if (lang)
			printf ("ingrese el número de prueba o AYUDA\n");
		else
			printf (" enter test number or HELP\n");
		gets_s(ibuff, sizeof(ibuff));
		if (strcmp(ibuff, "HELP") == 0 || strcmp(ibuff, "help") == 0
			|| strcmp(ibuff, "AYUDA") == 0 || strcmp(ibuff, "ayuda") == 0) {
			if (lang)
				printf(ayuda);
			else
				printf(helpmessage);
			printf("enter test number\n");
			gets_s(ibuff, sizeof(ibuff));
		}
		while (!isdigit(ibuff[0])) {
			printf("input must be a number!\n");
			gets_s(ibuff, sizeof(ibuff));
		}
		debug_code = atoi(ibuff);  // returns 0 if input not valid
		clock_t beginTime = clock();


		/*	1	fsquaresText   
			2	fcubesText     (arg1, 6)
			3	multiply & squareRoot
			4	BigIntDivide
			5	FcubesText   (series of values)
			6	contfracText  (arg1 arg2 arg3)
			7	fsquaresText (arg1, 6)
			8	ModInvBigNbr
			9	polyFactText  // not working
			10  Factorisation
			11	DilogText
			12 gaussianText
			13 emcFrontText  predefined values
			14 Dec2Bin & Bin2Dec
			15 doWork
			16 quadmodText
		*/
		switch (debug_code) {
		case 0: break;
		case 1: {		//#if DEBUG_CODE == 1
			//  fsquaresText("n(10^32)", 6);
			fsquaresText(NULL, 6);
			printf("%s\n", output);
			break;
		}
		case 2: {		//#elif DEBUG_CODE == 2
			fcubesText(argv[1], 6);
			printf("%s\n", output);
			break;
		}
		case 3: {		//#elif DEBUG_CODE == 3
			for (i = 0; i < 20; i++)
			{
				Factor1[i].x = 7 - i;
				Factor2[i].x = 1 + i;
				Product[i].x = 21 + i;
			}
			multiply(Factor5, Factor6, Product, 4, NULL);
			for (i = 0; i < 32; i++)
			{
				printf("%04X ", Product[i].x);
			}
			printf("\n");
			Factor1[0].x = 0;
			Factor1[1].x = 0;
			Factor1[2].x = 0;
			Factor1[3].x = 0;
			Factor1[4].x = 728;
			Factor1[5].x = 32767;
			Factor1[6].x = 32767;
			squareRoot(&Factor1[2], Factor2, 4, &len);
			//squareRoot(Factor1, Factor2, 7, &len);
			break;
		}
		case 4: {		//#elif DEBUG_CODE == 4
			dividend.limbs[0].x = 0;
			dividend.limbs[1].x = 0;
			dividend.limbs[2].x = 0;
			dividend.limbs[3].x = 0x20;
			dividend.nbrLimbs = 4;
			dividend.sign = SIGN_POSITIVE;
			divisor.limbs[0].x = 0;
			divisor.limbs[1].x = 0x05;
			divisor.nbrLimbs = 2;
			dividend.sign = SIGN_NEGATIVE;
			BigIntDivide(&dividend, &divisor, &quotient);
			break;
		}
		case 5: {		//#elif DEBUG_CODE == 5
			strcpy(expr, "123456789012345");
			for (i = sizeof(expr) - 1; i >= 0; i -= 3)
			{
				expr[i] = 0;  // reduce length of string
				fcubesText(expr, 6);
				printf("%s\n\n", output);
			}
			break;
		}
		case 6: {		//#elif DEBUG_CODE == 6
			if (argc != 4)
			{
				printf("num delta den\n");
				break;
			}
			ptrInput = input;
			strcpy(ptrInput, argv[1]);
			ptrInput += strlen(ptrInput) + 1;
			strcpy(ptrInput, argv[2]);
			ptrInput += strlen(ptrInput) + 1;
			strcpy(ptrInput, argv[3]);
			ptrInput += strlen(ptrInput) + 1;
			contfracText(input, 20000);
			printf("%s\n", output);
			break;
		}
		case 7: {		//#elif DEBUG_CODE == 7
			if (argc != 2)
			{
				printf("num\n");
				return 0;
			}
			ptrInput = input;
			strcpy(ptrInput, argv[1]);
			fsquaresText(input, 6);
			printf("%s\n", output);
			break;
		}
		case 8: {		//#elif DEBUG_CODE == 8
			;
			enum eExprErr rc;
			int NumberLength;
			BigInteger num, mod, inv;

			if (argc != 3)
			{
				printf("num mod\n");
				break;
			}
			rc = ComputeExpression(argv[1], 1, &num);  // num = arg1
			if (rc != EXPR_OK)
			{
				textError(output, rc);
			}
			else
			{
				rc = ComputeExpression(argv[2], 1, &mod);  // mod = arg2
				if (rc != EXPR_OK)
				{
					textError(output, rc);
				}
				else
				{
					NumberLength = mod.nbrLimbs;
					if (num.nbrLimbs < mod.nbrLimbs)
					{
						memset(&num.limbs[num.nbrLimbs], 0, (mod.nbrLimbs - num.nbrLimbs) * sizeof(limb));
					}
					memcpy(TestNbr, mod.limbs, NumberLength * sizeof(limb));
					TestNbr[NumberLength].x = 0;
					GetMontgomeryParms(NumberLength);
					ModInvBigNbr(num.limbs, inv.limbs, mod.limbs, NumberLength);
					Bin2Dec(inv.limbs, output, NumberLength, 200);  // convert number to text
				}
			}
			printf("%s", output);
			break;
		}
				//case 9: {			//#elif DEBUG_CODE == 9
				//	if (argc != 3)
				//	{
				//		printf("modulus polynomial\n");
				//		return 0;
				//	}
				//	polyFactText(argv[1], argv[2], 6);
				//	printf("%s\n", output);
				//	break;
				//}
		case 10: {
			output[0] = '\0';
			valuesProcessed = 0;
			ecmFrontText(NULL, performFactorization, NULL, output);
			printf("%s\n", output);
			break;
		}
		case 11: {			//#elif DEBUG_CODE == 11
			if (argc != 4)
			{
				printf("base power modulus\n");
				break;
			}
			dilogText(argv[1], argv[2], argv[3], 6);
			printf("%s\n", output);
			break;
		}
		case 12: {			//#elif DEBUG_CODE == 12
			if (argc != 3)
			{
				printf("value factorize\n");
				break;
			}
			gaussianText(argv[1], argv[2][0]);
			printf("%s\n", output);
			break;
		}
		case 13: {			//#elif DEBUG_CODE == 13
			lang = 0;      // set Language (0=English or 1=Spanish)
			hexadecimal = 0;
			if (argc == 3)
			{
				ecmFrontText(argv[1], 1, argv[2], output);
				printf("%s\n", output);
			}
			else if (argc == 2 || !debug)
			{
#ifdef _DEBUG
				char *ptrKnownFactors = strchr(argv[1], '=');
#else
				char *ptrKnownFactors = NULL;
#endif
				char text[100];
#if 1
				strcpy(text, "10**45+572");
				//    strcpy(text, "x=10**45+572;x=x+1;c<1000;x");
				ecmFrontText(text, 1, ptrKnownFactors, output);
				printf("%s\n", output);

				output[0] = '\0';
				ecmFrontText(NULL, 0, NULL, output);
				printf("%s\n", output);

				output[0] = '\0';
				valuesProcessed = 0;
				strcpy(text, "10**45+573");
				//    strcpy(text, "x=1;x=x+1;x<1001;c");
				ecmFrontText(text, 1, NULL, output);
				printf("%s\n", output);

				output[0] = '\0';
				ecmFrontText(NULL, 0, NULL, output);
				printf("%s\n", output);
				break;
#endif
				if (ptrKnownFactors)
				{                          // There is equal sign.
					*ptrKnownFactors = 0;    // Replace equal sign by string terminator.
					ptrKnownFactors++;
				}
				sprintf(text, "%s\n", argv[1]);
				ecmFrontText(text, 0, ptrKnownFactors, output);
				printf("%s\n", output);
			}
			else
			{
				printf("value [known factors]\n");
			}
			break;
		}
		case 14: {		//#elif DEBUG_CODE == 14

			limb internalNotation[100];
			static int bitGroups;
			char textInput[500];
			char textOutput[500];
			textInput[0] = '1';
			memset(&textInput[1], '0', 150);
			textInput[151] = 0;
			Dec2Bin(textInput, internalNotation, 151, &bitGroups);
			Bin2Dec(internalNotation, textOutput, 17, 6);
			textOutput[200] = 0;
			printf("result is: %s\n", textOutput);
			break;
		}
		case 15: {	
			/*note: inputstring consists of a series of null-terminated strings*/
			memcpy(inputString, "6,-2,,1234,00102^1042+1"
				"\0"
				"2^1042+1=5^1(0)*16673^1(0)*627186185377^1(16673)*131294792925870751515684960383613518415615538737991528767912593379854404518341858118366491474959205710499826133822402120149306175263402700301^1(16673)*6864797660130609714981900799081393217269435300143305409394463459185543183397652346775704046543201000705776033378429553397612687501667381169885775070966579201^1(2)"
				"\0\0"
				"222"
				"\0", 
				393 + 1);
			doWork();
			printf("%s\n", output);
			break;
		}
		case 16: {		// #elif DEBUG_CODE == 16
			quadmodText("1", "0", "-41", "5^10", 6);
			//quadmodText("7", "3", "5", "77", 6);
			//quadmodText("8", "3", "7", "16", 6);
			//quadmodText("1", "1", "-42", "10000", 6);
			//quadmodText("1", "0", "-41", "625", 6);
			//quadmodText("1", "1", "0", "112856782", 6);
			//  quadmodText("1", "1", "0", "56428391", 6);
			//  quadmodText("1", "1", "0", "2", 6);
			printf("%s\n", output);
			break;
		}
		case 17: {
			printf("Enter option: options are: \n S (Spanish)  \t-S (English)\n"
				"P (Pretty Print \t-P (normal print)\n"
				"O (Only evaluate i.e. compute value of expr but no factorisation etc\n"
				"-O (evaluate expr, then do factorisation etc\n");
			int ic = getchar();
			switch (toupper(ic)) {
				case 'S':
					lang = 1;		// change language to Spanish
					break;
				case 'P':
					prettyprint = true;
					break;
				case 'O':
					performFactorization = false;
					break;
				case '-': {
					int ic = getchar();
					switch (toupper(ic)) {
					case 'S':
						lang = 0;		// change language to English
						break;
					case 'P':
						prettyprint = false;
						break;
					case 'O':
						performFactorization = true;
						break;
					default:
						printf("invalid option\n");
						break;
					}
					break;
				}
				default:
					printf("invalid option\n");
					break;
			}
			while (ic != '\n')
				ic = getchar();
			break;
		}

		default:
			fprintf(stderr, "** invalid test number selected \n");
		}
	}
	return EXIT_SUCCESS;
}
