﻿#include "pch.h"
#include <fcntl.h>
#include <io.h>
#include "resource.h"

struct HelpDialogResp{
    int radiobutton;
    int Ok_Cancel;
};

extern HANDLE hConsole;   /* STD_OUTPUT_HANDLE */

static HelpDialogResp hResp;  /* user choices in the help dialog box are saved here */

/* process messages from dialog box for Help topics
returns TRUE(1) or FALSE(0) */
static INT_PTR helpDialogAct(HWND DiBoxHandle,
    UINT message,
    WPARAM additionalInfo1,
    LPARAM additionalInfo2) {

    int wpHi = HIWORD(additionalInfo1);     /* 0 or control notification code */
    int wpLo = LOWORD(additionalInfo1);     /* menu or control identifier */

    int button = 0;  /* set according to radio button selected */
 
    switch (message) {

    case WM_DESTROY:        /* 0x02 */
    case  WM_MOVE:          /* 0x03 */
    case WM_ACTIVATE:       /* 0x06 */
    case WM_SETFOCUS:       /* 0x07 */
    case WM_KILLFOCUS:      /* 0x08 */
    case WM_PAINT:          /* 0x0f */
    case WM_ERASEBKGND:     /* 0x14 */
    case WM_SHOWWINDOW:     /* 0x18 */
    case WM_ACTIVATEAPP:    /* ox1c */
    case WM_SETCURSOR:      /* 0x20 */
    case WM_MOUSEACTIVATE:  /* 0x21 */
    case WM_GETMINMAXINFO:  /* 0x24 */
    case WM_SETFONT:        /* 0x30 */
    case WM_WINDOWPOSCHANGING:  /* 0x46 */
    case WM_WINDOWPOSCHANGED:  /* 0x47 */
    case WM_NOTIFY:            /* 0x4e */
    case WM_HELP:              /* ox53 */
    case WM_NOTIFYFORMAT:      /* 0x55 */
    case WM_GETICON:           /* 0x7f */
    case WM_NCDESTROY:         /* 0x82 */
    case WM_NCHITTEST:     /* 0x84 */
    case WM_NCPAINT:       /* 0x85 */
    case WM_NCACTIVATE:    /* 0x86 */
    case 0x90:             /* can't find any documentation for this */
    case WM_NCMOUSEMOVE:   /* 0xa0 */
    case WM_NCLBUTTONDOWN: /* 0xa1 */
    case 0xae:             /* can't find any documentation for this */
        return FALSE;

    case WM_INITDIALOG:    /* 0x110 */ {
        /* select 1st radio button */
        BOOL rv = CheckRadioButton(DiBoxHandle, general, pari, general);
        if (rv == FALSE) {
            ErrorDisp(__FUNCTION__);
        }
        /* return true so system sets focus to 1st control */
        return TRUE;
    }

    case WM_COMMAND: /* 0x111 control selected by user */
        switch (wpLo) {  /* switch according to control selected */
        case general:
        case function:
        case expression:
        case Other:
        case YAFU:
        case Test:
        case Msieve:
        case Background:
        case Loop:
        case QMES:
        case pari:
            button = wpLo - general;  /* general -> 0, function -> 1 etc. */
            hResp.radiobutton = button;
            //if (verbose >0)
            //    std::cout << "button = " << buttText[button] << '\n';
            return FALSE;

        case IDOK:       /* OK         */
        case IDCANCEL:   /* cancel */
            hResp.Ok_Cancel = wpLo;
            EndDialog(DiBoxHandle, wpLo);  /* close dialog box */
            return TRUE;

        default:  /* unknown control*/
            std::cout << "HelpDialog WM_COMMAND  wpHi = " << wpHi << " wpLo = " << wpLo
                << " info2 = " << additionalInfo2 << '\n';
        }
        return FALSE;

    case WM_SYSCOMMAND:             /* 0x112 */
    case WM_CHANGEUISTATE:          /* 0x127 */
    case WM_UPDATEUISTATE:          /* 0x128 */
    case WM_QUERYUISTATE:           /* 0x129 */
    case WM_CTLCOLORBTN:            /* 0x135 */
    case WM_CTLCOLORDLG:            /* 0x136 */
    case WM_CTLCOLORSTATIC:         /* 0x138 */
    case WM_MOUSEFIRST:             /* 0x200 */
    case WM_LBUTTONDOWN:            /* 0x201 */
    case WM_LBUTTONUP:              /* 0x202 */
    case WM_CAPTURECHANGED:         /* 0x215 */
    case WM_MOVING:                 /* 0x216 */
    case WM_DEVICECHANGE:           /* 0X219 */
    case WM_ENTERSIZEMOVE:          /* 0x231*/
    case WM_EXITSIZEMOVE:           /* 0x232 */
    case WM_IME_SETCONTEXT:         /* 0x281 */
    case WM_IME_NOTIFY:             /* 0x282 */
    case WM_NCMOUSELEAVE:           /* 0x2a2 */
    case WM_PRINTCLIENT:            /* 0x318 */
    case WM_DWMNCRENDERINGCHANGED:  /* 0x31f */
    case WM_USER:                   /* 0x400 */
        return FALSE;

    default:  /* unexpected message type */
        printf_s("HelpDialog msg = %x Info1 = %lld info2 =%lld \n", message, additionalInfo1,
            additionalInfo2);
        return FALSE;
    }

    return FALSE;
}

/* display a dialog box so that the user can select a help topic */
static long long helpdiag(void) {
    hResp.radiobutton = 0;  /* set default value */
    auto rv = DialogBoxParamW(GetModuleHandle(nullptr), MAKEINTRESOURCE(help_dialog),
        handConsole, helpDialogAct, (LPARAM)99L);
    if (rv != IDOK && rv != IDCANCEL) {
        std::cout << "rv = " << rv << '\n';
        ErrorDisp(__FUNCTION__);
        return rv;
    }

    return rv;
}


/*search the docfile for the required help topic. The topic can be specified in
the command parameter. If this is not present the topic can be selected via a dialog box;
the required entry is specified by hResp.radiobutton.
Just search the file for the heading, and print everything until the next
heading is found.
If the help file is not found it is possible to change the path to access the file.
This change is retained in the BigIntCalculator.ini file.
The help file is not limited to ASCII characters. Most UTF-8 characters can be printed. */
void helpfunc(const std::vector<std::string>& command)
{
    FILE* doc;
    std::string nstr;
    wchar_t wstr[1025];
    bool printtopic = false;
    bool line1 = true;
    std::string expr = " ";
    char* newpathC;
    std::string newpath;
    const unsigned char BOM[] = { 0xEF, 0xBB, 0xBF };   /* byte order marker for UTF-8 */
    bool UTF8 = false;          /* set true if UTF8 BOM found */
    std::string helptopic;
    int lineCount = 0;

    /* list of topics that can be selected (usually via the help dialog box). 
       If new topics are addded to the help file the dialog box and associated 
       functions should also be updated. */
    const char buttText[][20] = { "HELP", "FUNCTION", "EXPRESSION", "OTHER",
            "YAFU", "TEST", "MSIEVE", "BACKGROUND", "LOOP", "QMES",
            "PARI", "AYUDA" };

    if (command.size() > 1) {
        helptopic = command[1];
    }
    else
        helptopic.clear();

    if (lang == 0) {     /* if English select one from multiple topics  */
        if (helptopic.empty()) {
            helpdiag();  /* result saved in global hResp */
            if (verbose > 0 && hResp.Ok_Cancel == IDOK) {
                std::cout << "button = " << buttText[hResp.radiobutton] << '\n';
            }
            if (hResp.Ok_Cancel == IDCANCEL) {
                std::cout << "help cancelled \n";
                return;
            }
            helptopic = buttText[hResp.radiobutton];
        }
    }
    else   /* español. Spanish language help is not divided into topics. */
        if (helptopic.empty())
            helptopic = buttText[11];  /* AYUDA */
retry:
    //open the doc file and search for a matching topic
    errno_t ecode = fopen_s(&doc, helpFilePath.data(), "r");
    if (ecode != 0) {
        /* failed to open the help file*/
        std::perror("fopen error: help file not found");

        while (std::toupper(expr[0]) != 'Y') {
            std::cout << "Do you want to search for the help file? (Y/N) \n";
            std::getline(std::cin, expr);
            if (std::toupper(expr[0]) == 'N')
                return;
        }
        newpathC = getFileName("Text\0*.TXT\0\0", handConsole);
        if (newpathC == NULL) {
            std::cout << "command cancelled \n";
            return;
        }
        else {
            helpFilePath = newpathC; /* copy new path for doc file to permanent storage */
            writeIni();       /* update the .ini file*/
            goto retry;       /* let's try to open the file again */
        }
    }

    /* doc file has been opened successfully */
    /* change stdout mode to support unicode (need to use wprintf, not printf)
    this allows multi-byte charactes to be printed */
    std::fflush(stdout);
    int oldmode = _setmode(_fileno(stdout), _O_U8TEXT);
    if (verbose > 0)
        wprintf_s(L"searching for help on '%S'\n", helptopic.c_str());

    /* exit this loop when reached EOF or the next topic after the one required
    is reached */
    while (!std::feof(doc)) {
        nstr.resize(1024);
        char* rv = std::fgets(nstr.data(), (int)nstr.size(), doc);   //read a line
        if (rv == NULL)
            if (std::feof(doc)) {
                break;
            }
            else {
                std::perror("fgets error reading Help file");
                break;
            }

        if (line1) {
            line1 = false;     /* only do this on 1st line */
            int d = (unsigned char)nstr[0] - BOM[0];    /* BOM has to be declared as unsigned char */
            int d1 = (unsigned char)nstr[1] - BOM[1];
            int d2 = (unsigned char)nstr[2] - BOM[2];
            if (d == 0 && d1 == 0 && d2 == 0) {
                nstr.erase(0, 3);  /* remove BOM */
                UTF8 = true;
            }
        }

        //is this a header?
        size_t lastchar = strlen(nstr.c_str())-1;  /* index of \n char*/
        if ((nstr.front() == '[') && (nstr[lastchar-1] == ']')) {
            if (_stricmp(helptopic.c_str(), "LIST")== 0) {
                wprintf_s(L"%S", nstr.c_str()); /* print header*/
            }
            else {
                if (printtopic)
                    break;  /* we have reached the start of the next topic, so exit
                               Only print 1 topic per help command */

                               //does it match our topic?
                nstr.resize(lastchar-1);  /* remove ']' and '\n'  */
               
                if (_stricmp(helptopic.c_str(), nstr.c_str()+ 1) == 0)
                    /* we get a match if the topic between [ and ] = helptopic */
                    printtopic = true;   /* we have found the required topic*/
                lineCount = 0;    /* reset line count */
            }
        }
        else {  /* not a header line */
            if (printtopic) {
                /* print only if within the required topic */
                if (lineCount > 27) {
                    std::wcout << L"** More (y/n) ? ";  /* pause every 27 lines of output */
                    std::getline(std::cin, expr);
                    expr[0] = std::toupper(expr[0]);
                    if (expr[0] != 'Y')
                        break;
                    else
                        lineCount = 0;
                }
                /* convert ascii or utf-8 to wide chars, print line */
                auto rv = MultiByteToWideChar(CP_UTF8, 0, nstr.c_str(), 
                    (int)nstr.size(), wstr, sizeof(wstr)/2);
                assert(rv > 0);
                wprintf_s(L"%s", wstr);
                lineCount++;
            }
        }
    }

    if (std::feof(doc)) {
        if (printtopic || (_stricmp(helptopic.c_str(), "LIST") == 0))
            std::wprintf(L"\n");   /* contrary to the POSIX standard, the last line of the file
                            may not end with newline */
        else
            if (lang)
                wprintf_s(L"Ayuda para %S no encontrado \n", helptopic.c_str());
            else
                wprintf_s(L"Help for %S not found \n", helptopic.c_str());
    }
    /* change stdout back to normal */
    std::fflush(stdout);
    _setmode(_fileno(stdout), oldmode);
    std::fclose(doc);    /* close doc file */
    return;
}