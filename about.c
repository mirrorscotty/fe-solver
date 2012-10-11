/**
 * @file about.c
 * This file just handles stuff with the about dialog for the program.
 */


#include <string.h>
#include "about.h"

#define LINEBRK strcat(result, "<br/>")

/**
 * @brief Make the text to in the about dialog.
 * @return A character pointer to the text.
 */
char *genAboutString()
{
    char *progName = "Finite Element Solver";
    char *version = "0.06";
    char *sourceURL = "https://github.com/mirrorscotty/fe-solver/";

    char *copyrightYear = "2012";
    char *author = "Alex Griessman";
    char *email = "agriessm@purdue.edu";

    char *license =
            "Redistribution and use in source and binary forms, with or without "
            "modification, are permitted provided that the following conditions are met: "
            "<ul>"
            "<li>Redistributions of source code must retain the above copyright "
                 "notice, this list of conditions and the following disclaimer.</li>"
            "<li>Redistributions in binary form must reproduce the above copyright "
                  "notice, this list of conditions and the following disclaimer in the "
                  "documentation and/or other materials provided with the distribution.</li>"
            "<li>Neither the name of Purdue University nor the "
                  "names of its contributors may be used to endorse or promote products "
                  "derived from this software without specific prior written permission.</li>"
            "</ul>"

            "THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS \"AS IS\" AND "
            "ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED "
            "WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE "
            "DISCLAIMED. IN NO EVENT SHALL Alex Griessman BE LIABLE FOR ANY "
            "DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES "
            "(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; "
            "LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND "
            "ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT "
            "(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS "
            "SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.";


    /* Code to compile the above strings into something useful! */
    char* result;
    int length = strlen(progName)
            +strlen(version)
            +strlen(sourceURL)
            +strlen(copyrightYear)
            +strlen(author)
            +strlen(email)
            +strlen(lisence)
            +400; /* Extra wiggle room */

    result = (char*) calloc(length, sizeof(char));
    strcpy(result, "<b>");
    strcat(result, progName);
    strcat(result, " Version ");
    strcat(result, version);
    strcat(result, "</b> ");
    LINEBRK;
    strcat(result, "The source code for this program may be downloaded from <a href=\"");
    strcat(result, sourceURL);
    strcat(result, "\">");
    strcat(result, sourceURL);
    strcat(result, "</a>");
    LINEBRK;
    LINEBRK;
    strcat(result, "Copyright (c) ");
    strcat(result, copyrightYear);
    strcat(result, ", ");
    strcat(result, author);
    strcat(result, " (<a href=\"");
    strcat(result, email);
    strcat(result, "\">");
    strcat(result, email);
    strcat(result, "</a>)");
    LINEBRK;
    strcat(result, "All rights reserved.");
    LINEBRK;
    strcat(result, "<hr/>");
    LINEBRK;
    strcat(result, license);

    return result;
}
