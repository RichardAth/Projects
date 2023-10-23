/* vim: set et ts=4
 *
 * Copyright (C) 2015-2021 the json-parser authors  All rights reserved.
 * https://github.com/json-parser/json-parser
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *   notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *   notice, this list of conditions and the following disclaimer in the
 *   documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR AND CONTRIBUTORS ``AS IS'' AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHOR OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
 * OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE.
 */

#ifdef _MSC_VER
#ifndef _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
#endif
#endif



#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <mpir.h>
#include "boost/multiprecision/gmp.hpp" 
#include <vector>
typedef boost::multiprecision::mpz_int Znum;

#define ZT(a) a.backend().data()  /* access mpz_t within a Znum (Boost mpz_int)*/

#include "json.h"

 /*
  * Test for json.c
  *
  * Compile (static linking) with
  *         gcc -o test_json -I.. test_json.c ../json.c -lm
  *
  * Compile (dynamic linking) with
  *         gcc -o test_json -I.. test_json.c -lm -ljsonparser
  *
  * USAGE: ./test_json <json_file>
  */

/* forward declarations */
static void process_value(json_value* value, int depth);
static void process_value_s(json_value* value, int depth, const char* name,
    const int index);

/* indent according to depth */
static void print_depth_shift(int depth)
{
    int j;
    for (j = 0; j < depth; j++) {
        printf("  ");
    }
}

Znum ToBeFactored;
std::vector <Znum> factors;


/* check whether the product of all the factors is equal to the number to be
factored */
bool sanityCheck(const Znum &ToBeFactored, const std::vector<Znum> &factors) {
    Znum residue = ToBeFactored;
    Znum remainder;

    for (auto f : factors) {
        /* divide by each of the factors in turn, checking the remainder each time */
        mpz_fdiv_qr(ZT(residue), ZT(remainder), ZT(residue), ZT(f));
        assert(remainder == 0);
    }
    if (residue != 1) {
        gmp_printf("factors OK but residue after removing factors  = %Zd \n",
            ZT(residue));
        return false;
    }
    else {
        printf("all factors found \n");
        return true;
    }
}

/* save specified values. At the moment, we are only looking for the number 
to be factored and the factors found. These are decimal mumbers saved as text 
strings */
static void save_value(json_value* value, const int index, int depth) {
    json_type t = value->type;
    char* sp;
    Znum numValue;
    int rv;

    assert(t == json_string); /* at the moment, only strings can be processed. */
    sp = value->u.string.ptr;
    rv = mpz_set_str(ZT(numValue), sp, 10);
    if (rv != 0) {
        fprintf(stderr, "invalid value: not a decimal number: %s \n", sp);
        return;
    }
    switch (index) {
    case 0:
        gmp_printf("number to be factored = %Zd \n", ZT(numValue));
        ToBeFactored = numValue;
        break;
    case 1:
        gmp_printf("factor = %Zd \n", ZT(numValue));
        factors.push_back(numValue);
        break;

    default:
        abort();   /* WTF? */
    }
}

static void process_object_s(json_value* value, int depth, const char *name,
    const int index)
{
    int length, x;
    if (value == NULL) {
        return;
    }
    length = value->u.object.length;
    for (x = 0; x < length; x++) {
        if (strcmp(value->u.object.values[x].name, name) != 0)
            continue;
        //print_depth_shift(depth);
        //printf("object[%d].name = %s\n", x, value->u.object.values[x].name);
        process_value_s(value->u.object.values[x].value, depth + 1, name, index);
    }
}

static void process_object(json_value* value, int depth)
{
    int length, x;
    if (value == NULL) {
        return;
    }
    length = value->u.object.length;
    for (x = 0; x < length; x++) {
        print_depth_shift(depth);
        printf("object[%d].name = %s\n", x, value->u.object.values[x].name);
        process_value(value->u.object.values[x].value, depth + 1);
    }
}

static void process_array(json_value* value, int depth)
{
    int length, x;
    if (value == NULL) {
        return;
    }
    length = value->u.array.length;
    printf("array\n");
    for (x = 0; x < length; x++) {
        process_value(value->u.array.values[x], depth);
    }
}

static void process_array_s(json_value* value, int depth, const char *name,
    const int index)
{
    int length, x;
    if (value == NULL) {
        return;
    }
    length = value->u.array.length;
    //printf("array\n");
    for (x = 0; x < length; x++) {
        process_value_s(value->u.array.values[x], depth, name, index);
    }
}

static void process_value_s(json_value* value, int depth, const char* name,
    const int index) {
    if (value == NULL) {
        return;
    }
    switch (value->type) {
    case json_none:
    case json_null:
        break;
    case json_integer:
    case json_double:
    case json_string:
    case json_boolean:
        save_value(value, index, depth);
        break;

    case json_object:
        process_object_s(value, depth + 1, name, index);
        break;
    case json_array:
        process_array_s(value, depth + 1, name, index);
        break;

    default:
        printf("json: invalid type\n");
        break;
    }
}

static void process_value(json_value* value, int depth)
{
    if (value == NULL) {
        return;
    }
    if (value->type != json_object) {
        print_depth_shift(depth);
    }
    switch (value->type) {
    case json_none:
        printf("none\n");
        break;
    case json_null:
        printf("null\n");
        break;
    case json_object:
        process_object(value, depth + 1);
        break;
    case json_array:
        process_array(value, depth + 1);
        break;
    case json_integer:
        printf("int: %10ld\n", (long)value->u.integer);
        break;
    case json_double:
        printf("double: %f\n", value->u.dbl);
        break;
    case json_string:
        printf("string: %s\n", value->u.string.ptr);
        break;
    case json_boolean:
        printf("bool: %d\n", value->u.boolean);
        break;
    }
}

/* returns -1 for read error, +1 for parsing error, otherwise
treats each line as a separate json record; parses it and prints 
selected fields of last record in a decoded form. */
int process_file_s(FILE* fp, int* counter,
    const char* name_list[], const int name_list_size) {
    char* string_p = NULL;
    char buffer[4096] = "" ;
    json_char* json = NULL;
    json_value* value = NULL;

    *counter = 0;  /* counts number of records in file */
    for (;;) {
        string_p = fgets(buffer, sizeof(buffer), fp); /* get 1 record */
        if (string_p != NULL) {
            //printf("%s \n", buffer);
            json = (json_char*)buffer;
            json_value_free(value);  /* avoid memory leakage */
            value = json_parse(json, strlen(buffer));
            if (value == NULL) {
                fprintf(stderr, "Unable to parse data\n");
                return(1);
            }
             (*counter)++;  /* count number of records */
        }
        else {   /* error or end-of-file */
            if (feof(fp)) {
#ifdef _DEBUG
                process_value(value, 0);
#endif
                /* process last record read */
                factors.clear();
                for (int index = 0; index < name_list_size; index++)
                    process_value_s(value, 0, name_list[index], index);
                json_value_free(value);  /* avoid memory leakage */
                if (!factors.empty())
                    sanityCheck(ToBeFactored, factors);
                break;    /* normal exit */
            }
            else {
                perror("fgets error");
                return(-1);
            }
        }
    }
    return(0);  /* normal exit */
}

/* returns -1 for read error, +1 for parsing error, otherwise
treats each line as a separate json record; prints it, parses it
and prints it in a decoded form. */
int process_file(FILE *fp, int *counter) {
    char* string_p;
    char buffer[4096];
    json_char* json;
    json_value* value;

    *counter = 0;

    for (;;) {
        string_p = fgets(buffer, sizeof(buffer), fp);
        if (string_p != NULL) {
            //printf("%s \n", buffer);
            json = (json_char*)buffer;
            value = json_parse(json, sizeof(buffer));

            if (value == NULL) {
                fprintf(stderr, "Unable to parse data\n");
                return(1);
            }

            process_value(value, 0);
            json_value_free(value);  /* avoid memory leakage */
            (*counter)++;  /* count number of records */
        }
        else {   /* error or end-of-file */
            if (feof(fp))
                break;    /* normal exit */
            else {
                perror("fgets error");
                return(-1);
            }
        }
    }

    return(0);


    //file_contents = (char*)malloc(filestatus.st_size);
    //if (file_contents == NULL) {
    //    fprintf(stderr, "Memory error: unable to allocate %d bytes\n", file_size);
    //    return 1;
    //}
    //f_rt = fread(file_contents, 1, file_size, fp);
    ///* f_rt = number of elements read */
    //if (f_rt < file_size) {
    //    fprintf(stderr, "Unable to read content of %s\n", filename);
    //    fprintf(stderr, "Expected %d bytes, got %lld bytes \n", file_size, f_rt);
    //    /*fclose(fp);
    //    free(file_contents);
    //    return 1;*/
    //}
    //fclose(fp);
    //file_contents[f_rt] = '\0';  /* null terminate string*/
    //printf("%s\n", file_contents);

    //printf("--------------------------------\n\n");

    //json = (json_char*)file_contents;

    //value = json_parse(json, file_size);

    //if (value == NULL) {
    //    fprintf(stderr, "Unable to parse data\n");
    //    free(file_contents);
    //    exit(1);
    //}

    //process_value(value, 0);

    //json_value_free(value);
    //free(file_contents);
    //return 0;
}

int main(int argc, char** argv) {
    char* filename;
    FILE* fp;
    struct stat filestatus;
    int file_size;
    int counter = 0;
    int rv;
    const char *name_list[] = { "input-decimal", "factors-prime" };

    if (argc != 2) {
        fprintf(stderr, "%s <file_json>\n", argv[0]);
        return 1;
    }
    filename = argv[1];

    if (stat(filename, &filestatus) != 0) {
        fprintf(stderr, "File %s not found\n", filename);
        return 1;
    }
    file_size = filestatus.st_size;
        fp = fopen(filename, "rt");
    if (fp == NULL) {
        perror(NULL);
        fprintf(stderr, "Unable to open %s\n", filename);
        fclose(fp);
        return 1;
    }
    rv = process_file_s(fp, &counter, name_list, 2);
    printf("file contains %d records, size = %.1f Kb\n", counter,
        (double)file_size/1024.0) ;
    fclose(fp);
    return rv;
}
