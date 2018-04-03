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

#include <string>

/* equivalent to:
_ui64toa((unsigned long long)nbr, *ptrOutput, 10));
*ptroutput += (strlen(*ptrOutput)); */
void int2dec(char **pOutput, long long nbr)
{
	char *ptrOutput = *pOutput;
	int significantZero = 0;
	unsigned long long div = 1000000000000000000;
	unsigned long long value = (unsigned long long)nbr;
	while (div > 0)
	{
		long long digit;

		digit = value / div;
		if (digit > 0 || significantZero != 0)
		{
			significantZero = 1;
			*ptrOutput++ = (char)(digit + (int)'0');
		}
		value %= div;
		div /= 10;
	}
	if (significantZero == 0)
	{
		*ptrOutput++ = '0';
	}
	*ptrOutput = '\0';   // null terminate O/P
	*pOutput = ptrOutput;
}
