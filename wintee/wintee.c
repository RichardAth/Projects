/* ***** BEGIN LICENSE BLOCK *****
 * Version: MPL 1.1
 *
 * The contents of this file are subject to the Mozilla Public License Version
 * 1.1 (the "License"); you may not use this file except in compliance with
 * the License. You may obtain a copy of the License at
 * http://www.mozilla.org/MPL/
 *
 * Software distributed under the License is distributed on an "AS IS" basis,
 * WITHOUT WARRANTY OF ANY KIND, either express or implied. See the License
 * for the specific language governing rights and limitations under the
 * License.
 *
 * The Original Code is Point2 Inc. code.
 *
 * The Initial Developer of the Original Code is Ryan Buhl
 * Portions created by the Initial Developer are Copyright (C) 2008
 * the Initial Developer. All Rights Reserved.
 *
 * ***** END LICENSE BLOCK ***** */

 /*
	 WinTee.c
	 Original code written by Ryan Buhl (rsbuhl@gmail.com)

	 This application approximates the tee command in *nix enviroments which
	 allow the user to pipe standard input to standard output and one or more
	 files.

	 The text of the usage display was borrowed from the usage text displayed by
	 the GNU coreutils command tee from *nix.  No GNU code was used within this
	 application.  Any similarities between this and other tee type commands is
	 coincidence and strictly a case of "great minds think alike".
 */

//#define _CRT_SECURE_NO_WARNINGS
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <io.h>
#include <fcntl.h>
#include <sys\stat.h>
#include "queue.h"
#include <dos.h>
#include <string.h>
#include <windows.h>

#define MAX_BUF_SIZE 1024

#undef FALSE                   /* remove definitions buried somewhere in windows.h */
#undef TRUE
#define FALSE 0
#define TRUE !FALSE

#define VERSION "1.0.1"
#define AUTHORS "Ryan Buhl"
#define AUTHOR_CONTACT "rbuhl@point2.com"

int open_outfile(int argc, char*argv[]);
int get_opt(int argc, char *argv[]);

int output_to_stdout = FALSE;
int append_files = FALSE;
int ignore_interrupts = FALSE;
int show_help = FALSE;

/* Need list to keep track of the file handles opened for output so we can
	close them when we're done */

TAILQ_HEAD(tailhead, entry) head;
struct tailhead *headp;  /* Tail Queue head */
struct entry             /* Tail Queue entry node */
{
	int outfile;
	TAILQ_ENTRY(entry) entries; /* Tail Queue */
};

/*
	Function to display usage information to the user so they know how this
		application works.
*/
void display_help(void) {
	printf("Usage: wtee [OPTION]... [FILE]...\n");
	printf("Copy standard input to each FILE, and also to standard output\n\n");
	printf("-a \t\t append to the given FILEs, do not overwrite\n");
	printf("-i \t\t ignore interrupt signals\n");
	printf("-?, --help \t display this help and exit\n");
	printf("--version \t output version information and exit\n\n");
	printf("If a FILE is -, copy again to standard output\n\n");

	printf("Report bugs to %s\n", AUTHOR_CONTACT);
}

/**
	Function to display version and bug contact information
*/
void display_version(void) {
	printf("wtee %s\n", VERSION);
	printf("There is NO WARRANTY, to the extent permitted by law.\n");
	printf("\nWritten by: %s\n", AUTHORS);
}

/**
	Function to parse the command line and search for expected options prefixed
		with a '-' character.

	@param argc Count of paramenters contained within the argv parameter.
	@param argv Pointer to character arrays containing the command line.

	@returns 0 if the command line contains only valid options; non-zero if an invalid option is found.

	@note This function is *not* the equivelant of the *nix function getopt.
	@note Even though this function returns non-zero to indicate that an
		  invalid command argument was given, a non-zero return value does not
		  nescesarily indicate that the function failed.
*/
int get_opt(int argc, char *argv[]) {
	char **base_arg;
	int i, result = 0, processed_parameter = FALSE;
	char *p;

	if ((NULL == argv) || (NULL == argv[0]))
		return 0;

	base_arg = argv;

	for (i = 0; i < argc; i++)
	{
		p = base_arg[i];  // copy pointer to parameter
		switch (*p)      // switch on 1st character of parameter
		{
		default:
			break;      // parameter does not begin with - Ignore it for now.
		case '-':
			++p;
			while (*p) 	{
				switch (*p) {
				default:
					result = 1;
					break;
				case 'a':
					append_files = TRUE;
					processed_parameter = TRUE;
					break;
				//case 'i':
				//	disable(); /* disable interrupts */
				//	processed_parameter = TRUE;
				//	break;
				case '?':
					display_help();
					exit(0);
					break;
				case '-':
					p++;
					if (!_stricmp(p, "version"))
					{
						display_version();
						exit(0);    // note: in this case the function does not return; the program stops
					}
					else if (!_stricmp(p, "help"))
					{
						display_help();
						exit(0);
					}
					else
						result = 1; // -- is not followed by a valid parameter
				}
				p++;

			}
			continue;
		}
	}

	return result;
}

/**
	Function to open and add streams to the stream list.

	@param argc Number of arguments contained within argv.
	@param argv Pointer to list of stream names to open.

	@returns number of output files opened.
*/
int open_outfile(int argc, char *argv[]) {
	char **base_v = argv, *p;
	int fHandle;            // file handle
	errno_t errnum;
	struct entry *EntryP;
	int list_count = 0;

	argc--; base_v++;

	while (argc--) {
		p = *base_v;   // set p to point to a command parameter
		if (*p == '-') {
			if (1 != strlen(p))
				base_v++;  // ignore parameter of form "-....." for now
			else
				output_to_stdout = TRUE;
			continue;
		}

		if (!append_files)
			//fHandle = _open(*base_v++, O_CREAT | O_WRONLY, S_IREAD | S_IWRITE);
			errno = _sopen_s(&fHandle, *base_v++, 
				O_CREAT | O_WRONLY,   // oflag:  create file if it doesn't already exist, write only
				_SH_DENYNO,           // shflag: permits read and write access
				S_IREAD | S_IWRITE);  // pmode:  reading & writing permitted
		else
			//fHandle = _open(*base_v++, O_CREAT | O_APPEND | O_WRONLY, S_IREAD | S_IWRITE);
			errnum = _sopen_s(&fHandle, *base_v++, 
				O_CREAT | O_APPEND | O_WRONLY, 
				_SH_DENYNO, 
				S_IREAD | S_IWRITE );

		if (-1 < fHandle) {
			/* add to the list */
			if (NULL != (EntryP = malloc(sizeof(struct entry)))) {
				/* have a list node so update and move on */
				EntryP->outfile = fHandle;
				TAILQ_INSERT_TAIL(&head, EntryP, entries);
				list_count++;
			}
			else {
				/* no memory, so close the file */
				_close(fHandle);
			}
		}
		else {
			perror("file open failure");
		}
	}

	return list_count;
}

/**
	Function that writes a given number of bytes from a buffer to a given stream.

	@param file Stream to write the buffer to.
	@param buffer Pointer to the buffer to write to the stream.
	@param bytes_to_write Number of bytes to write from the buffer.

	@returns 0 if there were no problems writing to the stream; Non-zero if errors occured.
*/
int write_to_file(int file, char *buffer, int bytes_to_write) {
	int bytes_written, total_written = 0;
	char *p = buffer;

	while (0 < (bytes_written = _write(file, p, bytes_to_write))) {
		/* short write */
		bytes_to_write -= bytes_written;
		p += bytes_written;
		total_written += bytes_written;
	}

	if (bytes_written == -1)
		return -1;

	return bytes_to_write - total_written;
}

int main(int argc, char *argv[]) {
	char buf[MAX_BUF_SIZE + 1] = { 0 };
	int bytes_read;
	int STDIN, STDOUT;
	//int redirect;
	struct entry *np;

	/* increase process priority */
	SetPriorityClass(GetCurrentProcess(), ABOVE_NORMAL_PRIORITY_CLASS);

	get_opt(argc, argv);

	/* Init the head */
	TAILQ_INIT(&head);

	STDIN = _fileno(stdin);
	STDOUT = _fileno(stdout);

	open_outfile(argc, argv);
	/* while input, output */
	while (0 < (bytes_read = _read(STDIN, buf, MAX_BUF_SIZE)))
	{
		for (np = head.tqh_first; np != NULL; np = np->entries.tqe_next) {
			int rv = write_to_file(np->outfile, buf, bytes_read);
			if (rv != 0) {
				fprintf(stderr, "*** wintee error writing to file \n");
				abort();
			}
		}

		/* now write to the display */
		write_to_file(STDOUT, buf, bytes_read);

		/* a filename of '-' indicates that the output should go to stdout again */
		if (output_to_stdout)
			write_to_file(STDOUT, buf, bytes_read);
	}

	/* Clean up the outfile queue */
	while (head.tqh_first != NULL) {
		_close(head.tqh_first->outfile);
		TAILQ_REMOVE(&head, head.tqh_first, entries);
	}

	return 0;
}