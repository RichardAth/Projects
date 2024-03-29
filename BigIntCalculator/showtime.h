#pragma once
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
extern double originalTenthSecond;
extern int oldTimeElapsed;
double tenths(void);

void GetDHMS(char **pptrText, int seconds);
void GetDHMSt(char **pptrText, int tenths);

const char * myTime(void);    /* get hh:mm:ss */
const char* myTimePW(void);   /* get hh:mm:ss.msec */
const char* myTimeP();        /* get hh:mm:ss.msec */
const char* myDateTime();     /* Www Mmm dd hh:mm:ss yyyy */
