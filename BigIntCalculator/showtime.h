#pragma once

extern double originalTenthSecond;
extern int oldTimeElapsed;

double tenths(void);

void GetDHMS(char **pptrText, int seconds);
void GetDHMSt(char **pptrText, int tenths);

//void showElapsedTime(char **pptrOutput);
