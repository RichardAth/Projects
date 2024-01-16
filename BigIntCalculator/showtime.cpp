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
#include <chrono>
#include <sstream>
#include <iomanip>

#pragma warning(disable : 4996)

/* get clock time in 1/10th of a second */
double tenths(void) {
	clock_t now = std::clock();
	return (double)now / (CLOCKS_PER_SEC / 10); ;
}

/* Convert interval in seconds to days, hours, minutes and seconds.
& move *pptrText past added text (normally 14 chars)
output format is "nd nnh nnm nns"*/
void GetDHMS(char **pptrText, int seconds) {
	char *ptrText = *pptrText;
	
	auto len = std::sprintf(ptrText, "%dd %2dh %2dm %2ds \n", seconds / 86400, 
		(seconds / 3600) % 24, (seconds / 60) % 60, seconds % 60);
	*pptrText += len-1;    // pointer points at null terminator
}

/* Convert interval in tenths of a second to days, hours, minutes and 
seconds.tenths, & move *pptrText past added text (normally 16 chars)
output format is "nd nnh nnm nn.ns"*/
void GetDHMSt(char **pptrText, int tenths) {
	char *ptrText = *pptrText;
	int seconds = tenths / 10;

	auto len = std::sprintf(ptrText, "%dd %2dh %2dm %2d.%ds \n", seconds / 86400,
		(seconds / 3600) % 24, (seconds / 60) % 60, seconds % 60,
		tenths%10);
	*pptrText += len - 1;    // pointer points at null terminator
}

/* get current time of day in format hh:mm:ss */
const char* myTime(void) {
	static char timestamp[10];   // time in format hh:mm:ss
	struct tm newtime;

	const time_t current = std::time(NULL);  // time as seconds elapsed since midnight, January 1, 1970
	localtime_s(&newtime, &current);    // convert time to tm structure
	/* convert time to hh:mm:ss */
	std::strftime(timestamp, sizeof(timestamp), "%H:%M:%S", &newtime);
	return timestamp;
}

/* get current time of day in format hh:mm:ss.msec (windows-only version)*/
const char* myTimePW(void) {
	static char timestamp[15];   // time in format hh:mm:ss.msec
	FILETIME SystemTimeAsFileTime;
	SYSTEMTIME   sysTime;

	GetSystemTimePreciseAsFileTime(&SystemTimeAsFileTime);

	if (FileTimeToSystemTime(&SystemTimeAsFileTime, &sysTime)) {
		sprintf_s(timestamp, sizeof(timestamp), "%02d:%02d:%02d.%03d",
			sysTime.wHour, sysTime.wMinute, sysTime.wSecond, sysTime.wMilliseconds);
		return timestamp;
	}
	else {
		ErrorDisp("__FUNCTION__");
		return nullptr;
	}
}

/* this version does not rely on any windows-specific functions */
/* see https://stackoverflow.com/questions/16077299/how-to-print-current-time-with-milliseconds-using-c-c11/66291527#66291527
*/
const char* myTimeP()
{
	using namespace std::chrono;
	using clock = system_clock;

	static char timestamp[15];   // time in format hh:mm:ss.msec
	const auto current_time_point{ clock::now() };  /* get current time */

	/* convert to time_t */
	const time_t current_time{ clock::to_time_t(current_time_point) };

	/* convert to tm struct */
	const tm current_localtime{ *std::localtime(&current_time) };

	/* convert to time interval since start of epoch (1 Jan 1970) */
	const auto current_time_since_epoch{ current_time_point.time_since_epoch() };
	/* get number of milliseconds */
	const auto current_milliseconds{ duration_cast<milliseconds> (current_time_since_epoch).count() % 1000 };

	size_t offset = std::strftime(timestamp, sizeof(timestamp), "%T", 
		&current_localtime);  /* get HH:MM:SS into timestamp */
	/* append milliseconds */
	sprintf_s(timestamp + offset, sizeof(timestamp) - offset, ".%03lld", current_milliseconds);
	return timestamp;
}
