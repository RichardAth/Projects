#include "pch.h"
#include <wtypes.h>
#include <winver.h>
/*
*** Don't forget! *****
add to Project Properties->Linker->Input->Additional Dependencies -> Add Version.lib
*/

#include <windef.h>
#include <winbase.h>
#include <shlwapi.h>

/* returns version info from .exe file, also date & time file last modified */
void VersionInfo(const LPCSTR path, int ver[4], std::string &modified) {
	DWORD handle;
	std::string FileInfo;
	int len;

	VS_FIXEDFILEINFO *fixedinfo;  /* pointer to root block info */

//typedef struct tagVS_FIXEDFILEINFO
//{
//  DWORD   dwSignature;            /* e.g. 0xfeef04bd */
//	DWORD   dwStrucVersion;         /* e.g. 0x00000042 = "0.42" */
//	DWORD   dwFileVersionMS;        /* e.g. 0x00030075 = "3.75" */
//	DWORD   dwFileVersionLS;        /* e.g. 0x00000031 = "0.31" */
//	DWORD   dwProductVersionMS;     /* e.g. 0x00030010 = "3.10" */
//	DWORD   dwProductVersionLS;     /* e.g. 0x00000031 = "0.31" */
//	DWORD   dwFileFlagsMask;        /* = 0x3F for version "0.42" */
//	DWORD   dwFileFlags;            /* e.g. VFF_DEBUG | VFF_PRERELEASE */
//	DWORD   dwFileOS;               /* e.g. VOS_DOS_WINDOWS16 */
//	DWORD   dwFileType;             /* e.g. VFT_DRIVER */
//	DWORD   dwFileSubtype;          /* e.g. VFT2_DRV_KEYBOARD */
//	DWORD   dwFileDateMS;           /* e.g. 0 */
//	DWORD   dwFileDateLS;           /* e.g. 0 */
//} VS_FIXEDFILEINFO;

	auto BufLen = GetFileVersionInfoSizeA(path, &handle); /* get size needed for data */
	if (BufLen == 0) {
		ErrorDisp(__FUNCTION__);
		return;
	}

	FileInfo.resize(BufLen + 1);  // allocate storage 
	/* get version info into FileInfo */
	bool rv = GetFileVersionInfoA(path, handle, BufLen, (LPVOID)FileInfo.data());
	if (!rv) {
		std::cerr << "can't get version info \n";
		return;
	}

	/* read root block */
	rv = VerQueryValueA(FileInfo.data(), "\\", (LPVOID *)&fixedinfo, (PUINT) &len);
	if (!rv) {
		std::cerr << "can't get version info root block \n";
		return;
	}

	/* copy 4 version digits */
	ver[0] = fixedinfo->dwFileVersionMS >> 16;
	ver[1] = fixedinfo->dwFileVersionMS & 0xffff;
	ver[2] = fixedinfo->dwFileVersionLS >> 16;
	ver[3] = fixedinfo->dwFileVersionLS & 0xffff;

	/*  copy file date */
	struct _stat64 fileStat;
	struct tm ftimetm;
	
	int err = _stat64(path, &fileStat);
	if (err == 0) {
		modified.resize(31);
		auto ftime = fileStat.st_mtime;  // time last modified in time_t format
		localtime_s(&ftimetm, &ftime);   // convert to tm format
		/* convert to dd/mm/yyyy  at hh:mm:ss */
		std::strftime(&modified[0], 30, "%a %d/%m/%Y at %T", &ftimetm);
	}
	else
		std::cout << path << " not found \n";
}

#define PACKVERSION(major,minor) MAKELONG(minor,major)
typedef HRESULT(CALLBACK* DLLGETVERSIONPROC2)(DLLVERSIONINFO2*);

/* get the version number of the dll specified in DllName */
static bool GetVersion(LPCTSTR lpszDllName, DWORD * Major, DWORD * Minor, 
	DWORD * Build, DWORD * platform) {
	HINSTANCE hinstDll;
	DWORD dwVersion = 0;
	DLLVERSIONINFO2 dvi;
	HRESULT hr;
	DLLGETVERSIONPROC2 pDllGetVersion = nullptr;

	// For security purposes, LoadLibrary should be provided with a fully qualified 
	// path to the DLL. The lpszDllName variable should be tested to ensure that it 
	// is a fully qualified path before it is used. 
	hinstDll = LoadLibrary(lpszDllName);

	if (hinstDll) {
		pDllGetVersion = (DLLGETVERSIONPROC2)GetProcAddress(hinstDll, "DllGetVersion");

		// Because some DLLs might not implement this function, you must test for 
		// it explicitly. Depending on the particular DLL, the lack of a DllGetVersion 
		// function can be a useful indicator of the version. 

		if (pDllGetVersion) {
			ZeroMemory(&dvi, sizeof(dvi));
			dvi.info1.cbSize = sizeof(dvi);
			hr = (*pDllGetVersion)(&dvi);
			if (SUCCEEDED(hr)) 	{
				dwVersion = PACKVERSION(dvi.info1.dwMajorVersion, dvi.info1.dwMinorVersion);
				*Major = dvi.info1.dwMajorVersion;
				*Minor = dvi.info1.dwMinorVersion;
				*Build = dvi.info1.dwBuildNumber;
				*platform = dvi.info1.dwPlatformID;
			}
		}
		FreeLibrary(hinstDll);
	}
	return (pDllGetVersion != nullptr);
}


/* get version of ComCtl32.dll */
/* specifying just the file name, not the full path, ensures that the version 
from the manifest is used. */
const LPCTSTR lpszDllName = L"ComCtl32.dll";
// const LPCTSTR lpszDllName = L"C:\\Windows\\System32\\ComCtl32.dll";
// const LPCTSTR lpszDllName = L"C:/Windows/WinSxS/x86_microsoft.windows.common-controls_6595b64144ccf1df_6.0.19041.1110_none_a8625c1886757984/comctl32.dll";
/*                           "C:\Windows\WinSxS\x86_microsoft.windows.common-controls_6595b64144ccf1df_6.0.19041.1110_none_a8625c1886757984\comctl32.dll" */
DWORD getComCtlVer(void) {
	DWORD Major, Minor, build, platform;
	DWORD version;
	if (GetVersion(lpszDllName, &Major, &Minor, &build, &platform)) {
		if (verbose >= 1)
			printf_s("ComCtl32.dll version = %u.%u.%u \n", Major, Minor, build);

		version = PACKVERSION(Major, Minor);
		return version;
	}
	else return 0;
}