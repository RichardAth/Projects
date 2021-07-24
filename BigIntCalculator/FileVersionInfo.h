/********************************************************************
*
* Copyright (C) 1999-2000 Sven Wiegand
* Copyright (C) 2000-2001 ToolsCenter
*
* This file is free software; you can redistribute it and/or
* modify, but leave the headers intact and do not remove any
* copyrights from the source.
*
* If you have further questions, suggestions or bug fixes, visit
* our homepage
*
*    http://www.ToolsCenter.org
*
********************************************************************/



#include <winver.h>
#include <atltime.h>
#include "pch.h"

typedef unsigned short      WORD;
typedef unsigned long       DWORD;

class AFX_EXT_CLASS CFileVersionInfo
{
	// construction/destruction
public:
	CFileVersionInfo();
	virtual ~CFileVersionInfo();

	// operations
public:
	bool Create(HMODULE hModule = NULL);
	bool Create(LPCTSTR lpszFileName);

	// attribute operations
public:
	WORD GetFileVersion(int nIndex) const;
	WORD GetProductVersion(int nIndex) const;
	DWORD GetFileFlagsMask() const;
	DWORD GetFileFlags() const;
	DWORD GetFileOs() const;
	DWORD GetFileType() const;
	DWORD GetFileSubtype() const;
	CTime GetFileDate() const;

	CString GetCompanyName() const;
	CString GetFileDescription() const;
	CString GetFileVersion() const;
	CString GetInternalName() const;
	CString GetLegalCopyright() const;
	CString GetOriginalFileName() const;
	CString GetProductName() const;
	CString GetProductVersion() const;
	CString GetComments() const;
	CString GetLegalTrademarks() const;
	CString GetPrivateBuild() const;
	CString GetSpecialBuild() const;

	// implementation helpers
protected:
	virtual void Reset();
	bool GetTranslationId(LPVOID lpData, UINT unBlockSize, WORD wLangId, DWORD &dwId, bool bPrimaryEnough = FALSE);

	// attributes
private:
	VS_FIXEDFILEINFO m_FileInfo;

	CString m_strCompanyName;
	CString m_strFileDescription;
	CString m_strFileVersion;
	CString m_strInternalName;
	CString m_strLegalCopyright;
	CString m_strOriginalFileName;
	CString m_strProductName;
	CString m_strProductVersion;
	CString m_strComments;
	CString m_strLegalTrademarks;
	CString m_strPrivateBuild;
	CString m_strSpecialBuild;
};

