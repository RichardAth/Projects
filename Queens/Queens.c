// This program is intended to solve the Eight Queens Puzzle.
// The program starts by placing a queen in the 1st row, then placing a queen
// somewhere in the next row that is not blocked by any of the queens already
// in place. Then moving to the next row, and so on until the last row is reached.
// If no queen can be placed in a given row the program backtracks to the 
// previous row and moves the queen on that row to the right. If there are no more 
// squares available on a row to move the queen to, the queen is removed and the  
// program moves back up to the previous row. 

// the convention used is that row 0, column 0 is the top left corner of the board.

// the board size is defined by B_SIZE. To change the board size, alter this macro
// and recompile the program.

#include <stdio.h>
#include <stdbool.h>
#include <minmax.h>
#include <string.h>
#include <stdlib.h>
#include <windows.h>

#define B_SIZE 8		// N.B. If B_SIZE is > 13 there are more than 10,000 solutions!

int SolCnt = 0;
HANDLE stanout;

// This function checks whether there is already a queen in row 'row'.
// If so it returns the value 'false'. If the row is available it returns 'true'
_inline bool RowIsFree(const bool board[][B_SIZE], const int row) {

	for (int i = 0; i < B_SIZE; i++) {
		if (board[row][i]) {
#ifdef debug
			if (SolCnt > 0) printf("RowIsFree false row=%d, i=%d,\n",
				row, i);
#endif
			return false;
		}
	}
	return true;  // if we drop through to here the row is free
}

// This function checks whether there is already a queen in column 'col'.
// If so it returns the value 'false'. If the column is available it returns 'true'
_inline bool ColIsFree(const bool board[][B_SIZE], const int col) {

	for (int i = 0; i < B_SIZE; i++) {
		if (board[i][col]) {
#ifdef debug
			if (SolCnt > 0) printf("ColIsFree false col=%d, i=%d,\n",
				col, i);
#endif
			return false;
		}
	}
	return true; // if we drop through to here the column is free
}

// This function checks whether there is already a queen in the diagonal top left to
// bottom right which passes throug the square specified by row and col.
// If so it returns the value 'false'. If the diagonal is available it returns 'true'
bool LDiagIsFree(const bool board[][B_SIZE], const int row, const int col) {
	int x, sRow, sCol;

	x = min(row, col);
	sRow = row - x;   // either sRow, sCol, or both are zero
	sCol = col - x;

	while ((sRow < B_SIZE) && (sCol < B_SIZE)) {
		if (board[sRow] [sCol]) {
#ifdef debug
			if (SolCnt > 0) printf("LdiagIsFree false row=%d, col=%d, sRow=%d, sCol=%d\n",
				row, col, sRow, sCol);
#endif
			return false;
		}
		sRow++;
		sCol++;
	}
	return true;   // if we drop through to here the diagonal is free
}

// This function checks whether there is already a queen in the diagonal top right to
// bottom left which passes throug the square specified by row and col.
// If so it returns the value 'false'. If the diagonal is available it returns 'true'
bool RDiagIsFree(const bool board[][B_SIZE], const int row, const int col) {
	int x, sRow, xCol;

	xCol = (B_SIZE - 1) - col;
	x = min(row, xCol);
	sRow = row - x;   // either sRow, xCol, or both are zero
	xCol = xCol - x;

	while ((sRow < B_SIZE) && (xCol < B_SIZE)) {
		if (board[sRow] [(B_SIZE - 1) - xCol]) {
#ifdef debug
			if (SolCnt > 0) printf("RdiagIsFree false row=%d, col=%d, sRow=%d, sCol=%d\n",
				row, col, sRow, (B_SIZE - 1) - xCol);
#endif
			return false;
		}
		sRow++;
		xCol++;
	}
	return true;   // if we drop through to here the diagonal is free
}

// this function checks whether the square specified by row and col is free
// by checking whether there is already a queen in that row, column, or either
// diagonal. If there is it returns the value 'false'.
// This function also checks whether the square is blocked in SymChk. If it is then 
// the solution would be a rotation or reflection of an earlier solution, so in
// this case also, the function returns 'false', otherwise it returns 'true'.
bool SpaceIsFree(const bool board[][B_SIZE], const int row, const int col, 
	const bool SymChk[][B_SIZE]) {

#ifdef debug
	if (SolCnt > 0) printf("SpaceIsFree row=%d, col=%d\n", row, col);
#endif
	if (ColIsFree(board, col) &&
//		RowIsFree(board, row) &&   // don't need to check the row because the way
		                           // that queens are placed onto the board ensures 
		                           // that there is never more than 1 queen in one row 
		LDiagIsFree(board, row, col) &&
		RDiagIsFree(board, row, col) 
		&& (!SymChk[row][col])
        ) {
#ifdef debug
		if (SolCnt > 0) printf("SpaceIsFree row=%d, col=%d, line=%d, true\n", row, col, __LINE__);
#endif
		return true;
	}
	else {
#ifdef debug
		if (SolCnt > 0) printf("SpaceIsFree row=%d, col=%d, false\n", row, col);
#endif
		return false;
	}
}

// print the 'board'. A 'Q' represents a queen.  
void printB (const bool board[][B_SIZE]) {
	BOOL r;

	// print chessboard pattern
	for (int i = 0; i < B_SIZE; i++) {
		for (int j = 0; j < B_SIZE; j++) {
			if (!((i + j) % 2)) {
				r = SetConsoleTextAttribute(stanout,
					BACKGROUND_BLUE | BACKGROUND_GREEN | BACKGROUND_RED
					| BACKGROUND_INTENSITY);   // black text on white background.
			}
			else {
				r = SetConsoleTextAttribute(stanout,
					FOREGROUND_BLUE | FOREGROUND_GREEN | FOREGROUND_RED |
					FOREGROUND_INTENSITY);	   // white text on black background.
			}
			if (!r) {
				fprintf(stderr, "SetConsoleTextAttribute failed %d at line %d\n", GetLastError(), __LINE__);
				return;
			}
			if (board[i][j]) printf("Q ");
			else printf("  ");
		}
		r = SetConsoleTextAttribute(stanout,
			BACKGROUND_BLUE | BACKGROUND_GREEN | BACKGROUND_RED
			| BACKGROUND_INTENSITY);   // black text on white background.
		if (!r) {
			fprintf(stderr, "SetConsoleTextAttribute failed %d at line %d\n", GetLastError(), __LINE__);
			return;
		}
		printf("\n");
	}

	r = SetConsoleTextAttribute(stanout,
		BACKGROUND_BLUE | BACKGROUND_GREEN | BACKGROUND_RED
		| BACKGROUND_INTENSITY);   // black text on white background.
	if (!r) {
		fprintf(stderr, "SetConsoleTextAttribute failed %d at line %d\n", GetLastError(), __LINE__);
	}

	/*for (int r = 0; r < B_SIZE; r++) {
		for (int c = 0; c < B_SIZE; c++) {
			if (board[r][c]) printf("Q ");
			else printf(". ");
		}
		printf("\n");
	}*/
	return;
}

// print the 'board'. The board is stored in compact form.
void printC(const unsigned char board[]) {
	char pbuff[B_SIZE * 2 + 1];
	printf("\n");
	for (int i = 0; i < B_SIZE; i++) {

		// set up pbuff as blanks
		for (int j = 0; j < B_SIZE; j++) {
			pbuff[j * 2] = '.';
			pbuff[j * 2 + 1] = ' ';
		}

		pbuff[board[i] * 2] = 'Q';  // mark the position of the queen
		pbuff[B_SIZE * 2] = '\0';	// ensure string is null-terminated
		printf(pbuff); printf("\n");
	}
	for (int i = 0; i < B_SIZE; i++) {
		printf("%d,", board[i]);	// print the board array.
	}
	printf("\n");
}

// Check whether or not solution soln is in the list SolLst
_inline bool checklist(const unsigned char SolList[][B_SIZE], 
	const unsigned char soln[], const int SolCnt) {
	for (int i = 0; i < SolCnt; i++) {
		if (memcmp(SolList[i], soln, B_SIZE) == 0) {
			printf("Duplicate solution matches sol. %d ", i + 1);
			return true;
		}
	}
	return false;
}

// reflect the solution along the diagonal from top left to bottom right.
_inline void reflect(const unsigned char charS1[], unsigned char chars2[]) {
	for (int i = 0; i < B_SIZE; i++) {
		chars2[charS1[i]] = i;
	}
}

// rotate the solution 90 degrees clockwise.
_inline void rotate (const unsigned char charS1[], unsigned char chars2[]) {
	for (int i = 0; i < B_SIZE; i++) {
		chars2[charS1[i]] = B_SIZE - 1 -i;
	}
}

// Ckeck whether the solution on the board is a rotation or reflection of a 
// solution found previously. If so return 'true' otherwise return 'false'
// Only one reflection needs to be checked, on a diagonal from top left to
// bottom right. Other reflections have already been eliminated by another method.
bool DupChk(const bool board[][B_SIZE]) {

// for reasons of efficiency the solution on the board is converted
// to a compact representation. An 8X8 array is reduced to an 8 byte
// list. A board size > 256X256 could not be handled because the maximum
// column number that can be stored in 1 byte is 255.

#define MaxNoOfSols 10000	// N.B. If B_SIZE is > 13 there are more than 10,000 solutions!
	static unsigned char SolList[MaxNoOfSols][B_SIZE];  // stores all solutions found so far
														// in compact form. 
	static int SolCnt = 0;				// number of solutions in SolList
	unsigned char CurSol[B_SIZE], CurSol90[B_SIZE], CurSol180[B_SIZE], 
		CurSol270[B_SIZE], CurSolR[B_SIZE];

	// convert current board layout to compact form, store it in CurSol. Assume  
	// there is only 1 queen per row (which has to be true for a valid solution)
	for (int i = 0; i < B_SIZE; i++) 
		for (int j = 0; j < B_SIZE; j++ ) {
			if (board[i][j]) CurSol[i] = j;
	}

	rotate(CurSol, CurSol90);	// rotate board 90 degrees clockwise
	if (checklist(SolList, CurSol90, SolCnt)) {
		printf("rotated 90 degrees\n");
		return true;  // return 'true' if match found
	}

	rotate(CurSol90, CurSol180);	// rotate board 90 degrees clockwise
	if (checklist(SolList, CurSol180, SolCnt)) {
		printf("rotated 180 degrees\n");
		return true;  // return 'true' if match found
	}
	rotate(CurSol180, CurSol270);	// rotate board 90 degrees clockwise
	if (checklist(SolList, CurSol270, SolCnt)) {
		printf("rotated 270 degrees\n");
		return true;  // return 'true' if match found
	}

	reflect(CurSol, CurSolR);		// reflect board around diagonal.
	if (checklist(SolList, CurSolR, SolCnt)) {
		printf("diagonal reflection\n");
		return true;  // return 'true' if match found
	}
	if (SolCnt >= MaxNoOfSols) {
		printf("** too many solutions to process!! Abort at line %d\n", __LINE__);
		abort();   // SolList is full
	}

	// the current solution does not match any previous solution, so add it
	// to the list of solutions.
	memcpy(SolList[SolCnt], CurSol, B_SIZE);
	SolCnt++;
	return false;		// current solution does NOT match any previous one.
}

// set certain entries to 'true' to represent squares where a queen must not
// be placed in order to avoid solutions which are reflections or rotations
// of a previous solution. This suppresses most reflections and rotations but
// not all:

// If there is a queen in the top LH corner a duplicate solution reflected along 
// the top left to bottom R corner will not be suppressed.
// If there is not a queen in the top LH corner, reflections and rotations of 
// the reflections will be suppressed, but some rotations of some solutions will
// not be suppressed.
// For an 8x8 board there are 92 solutions, but only 12 when rotations and 
// reflections are suppressed. This check reduces the 92 to 18, i.e. it allows 6
// duplicate solutions.
void setSym(const int x, bool SymChk[][B_SIZE]) {
	static int xsave = -1;

	if (xsave == x) return;
	xsave = x;
	SymChk[0]             [x]              = true; // top row
	SymChk[0]             [B_SIZE - 1 - x] = true;
	SymChk[B_SIZE - 1]    [x]              = true; // bottom row
	SymChk[B_SIZE - 1]    [B_SIZE - 1 - x] = true;
	SymChk[x]             [0]              = true; // leftmost column
	SymChk[B_SIZE - 1 - x][0]              = true;
	SymChk[x]             [B_SIZE - 1]     = true; // rightmost column
	SymChk[B_SIZE - 1 - x][B_SIZE - 1]     = true;
#ifdef debug
	printf("setSym line %d, x=%d\n", __LINE__, x);
	for (int r = 0; r < B_SIZE; r++) {
		for (int c = 0; c < B_SIZE; c++) {
			if (SymChk[r][c]) printf("Q ");
			else printf(". ");
		}
		printf("\n");

	}
	printf("press any key\n");
	getc(stdin);
#endif 
	return;
}

/* see https://docs.microsoft.com/en-us/windows/console/setcurrentconsolefontex */
void setFontSize(int FontSize) {
	CONSOLE_FONT_INFOEX info = { 0 };
	info.cbSize = sizeof(info);
	info.dwFontSize.Y = FontSize; // leave X as zero
	info.FontWeight = FW_NORMAL;
	wcscpy(info.FaceName, L"Lucida Console");
	SetCurrentConsoleFontEx(GetStdHandle(STD_OUTPUT_HANDLE), FALSE, &info);
}

// function to clear the screen and reset cursor to top.
// use     hConsole = GetStdHandle (STD_OUTPUT_HANDLE) ;
// for the stdout console window
int clearScreen(HANDLE hConsole)
{
	COORD coordScreen = { 0, 0 };    // home for the cursor 
	DWORD cCharsWritten;
	DWORD dwConSize;
	CONSOLE_SCREEN_BUFFER_INFO csbi;

	// Get the current text attribute & screen size.

	if (!GetConsoleScreenBufferInfo(hConsole, &csbi))
	{
		fprintf(stderr, "** GetConsoleScreenBufferInfo failed with %d!\n", GetLastError());
		Beep(750, 1000);
		return EXIT_FAILURE;
	}

	// Get the number of character cells in the current buffer. 

	dwConSize = csbi.dwSize.X * csbi.dwSize.Y;

	// Fill the entire screen with blanks.

	if (!FillConsoleOutputCharacter(hConsole,       // Handle to console screen buffer 
		(TCHAR) ' ',     // Character to write to the buffer
		dwConSize,       // Number of cells to write 
		coordScreen,     // Coordinates of first cell 
		&cCharsWritten))// Receive number of characters written
	{
		fprintf(stderr, "** FillConsoleOutputCharacter failed! with %d\n", GetLastError());
		Beep(750, 1000);
		return  EXIT_FAILURE;
	}

	// Set the buffer's attributes accordingly.

	if (!FillConsoleOutputAttribute(hConsole,         // Handle to console screen buffer 
		csbi.wAttributes, // Character attributes to use
		dwConSize,        // Number of cells to set attribute 
		coordScreen,      // Coordinates of first cell 
		&cCharsWritten)) // Receive number of characters written
	{
		fprintf(stderr, "**  FillConsoleOutputAttribute failed with %d!\n", GetLastError());
		Beep(750, 1000);
		return  EXIT_FAILURE;
	}

	// Put the cursor at its home coordinates.

	SetConsoleCursorPosition(hConsole, coordScreen);
	return EXIT_SUCCESS;
}

int main(int argc, char ** argv)
{
	bool board [B_SIZE][B_SIZE];   // represents the chessboard. An array element is 'true'
	                               // if there is a queen in the corresponding position on
	                               // the board.
	bool SymChk [B_SIZE][B_SIZE];  // used to avoid solutions which are rotations or reflections
	                               // of a solution found earlier
	int bt[B_SIZE] = { 0 };   // this array is used when backtracking.
	int row = 0, col = 0;
	BOOL r;
 	system("PAUSE");

	setFontSize(16); /* set font and font size */

	// clear the board.
	for (int r = 0; r < B_SIZE; r++) {
		for (int c = 0; c < B_SIZE; c++) {
			board  [r][c] = false;
			SymChk [r][c] = false;
		}
	}

	stanout = GetStdHandle(STD_OUTPUT_HANDLE);
	if (stanout == INVALID_HANDLE_VALUE)
	{
		fprintf(stderr, "GetStdHandle failed with %d at line %d\n", GetLastError(), __LINE__);
		return EXIT_FAILURE;
	}

	r = SetConsoleTextAttribute(stanout,
		BACKGROUND_BLUE | BACKGROUND_GREEN | BACKGROUND_RED
		| BACKGROUND_INTENSITY);   // black text on white background.

	if (!r) {
		fprintf(stderr, "SetConsoleTextAttribute failed %d at line %d\n", GetLastError(), __LINE__);
		return;
	}


	clearScreen(stanout);

#ifdef debug
	printf("line %d\n", __LINE__);
#endif
	for (row = 0; row < B_SIZE; row++) {
		for (col = bt[row]; col < B_SIZE; col++) {
			if (SpaceIsFree(board, row, col, SymChk)) {
				board[row][col] = true;
				bt [row] = col+1;   // save position for backtracking.
#ifdef debug
				printf("line %d, row=%d, col=%d, bt=", __LINE__, row, col);
				for (int i = 0; i < B_SIZE; i++) printf(" %d, ", bt[i]) ;
				printf("\n");
#endif
				break;
			}
		}
#ifdef debug
		printf("line %d, row=%d, col=%d\n", __LINE__, row, col);
#endif

		if (col >= B_SIZE) {   // if failed to place a queen in this row
#ifdef debug
			printf("** line %d fail row %d\n", __LINE__, row);
#endif
			if (row == 0) break;   // row doesn't go -ve. 'break' jumps out 
								   // of for loop; exits from 'main'

#ifdef debug
			if (SolCnt > 0) {
				printB(board);
				printf("** fail line %d, row %d, bt=", __LINE__, row);
				for (int i = 0; i < B_SIZE; i++) printf(" %d, ", bt[i]);
				printf("  press any key\n");
				cc = getc(stdin);
			}
#endif
			board[row - 1][bt[row - 1] - 1] = false;  // remove queen in row above
			bt[row] = 0;  // start again from 1st col in this row
			row -= 2;     // adjust the for-loop counter to move back 1 row 
						  // instead of forward 1 row 
			              // (remember that for loop increments by 1)

		}
		if (row == (B_SIZE - 1)) { // is this the last row?

			// once we have a queen in the last row we have found a solution
			// so let's print it
			if (!DupChk(board)) {
				SolCnt++;
				printf("Solution %d\n", SolCnt);
				printB(board);
			}

			// next, we backtrack i.e remove the queens from the 2 bottom rows
			// and continue looking for more solutions.
			board[row][col] = false;
			bt[row] = 0;
			board[row - 1][bt[row - 1] - 1] = false;
			row -= 2;
#ifdef debug
			printf("line %d, row=%d, col=%d  ", __LINE__, row, col);
			printf("press any key to continue\n");
			cc = getc(stdin);
#endif
		}
		if (bt[0] >1) setSym(bt[0]-2, SymChk);  // set indicators to suppress
		                                        // reflections and rotations.
	}
	system("PAUSE");
	return 0;
}