/*
This file is part of Alpertron Calculators.
Copyright 2015 Dario Alejandro Alpern
Alpertron Calculators is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
Alpertron Calculators is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with Alpertron Calculators.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "pch.h"

#pragma warning(disable : 4996)

/* get clock time in 1/10th of a second */
double tenths(void) {
	clock_t now = std::clock();
	return (double)now / (CLOCKS_PER_SEC / 10); ;
}

/* Convert seconds to days, hours, minutes and seconds.
& move *pptrText past added text (normally 14 chars)
output format is "nd nnh nnm nns"*/
void GetDHMS(char **pptrText, int seconds) {
	char *ptrText = *pptrText;
	
	auto len = sprintf(ptrText, "%dd %2dh %2dm %2ds \n", seconds / 86400, 
		(seconds / 3600) % 24, (seconds / 60) % 60, seconds % 60);
	*pptrText += len-1;    // pointer points at null terminator
}

/* Convert tenths to days, hours, minutes and seconds.tenths.
& move *pptrText past added text (normally 16 chars)
output format is "nd nnh nnm nn.ns"*/
void GetDHMSt(char **pptrText, int tenths) {
	char *ptrText = *pptrText;
	int seconds = tenths / 10;

	auto len = sprintf(ptrText, "%dd %2dh %2dm %2d.%ds \n", seconds / 86400,
		(seconds / 3600) % 24, (seconds / 60) % 60, seconds % 60,
		tenths%10);
	*pptrText += len - 1;    // pointer points at null terminator
}
