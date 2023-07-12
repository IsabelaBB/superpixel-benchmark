#include "Utils.h"

//=============================================================================
// Void
//=============================================================================
void printError(const char* function_name, const char* message, ...)
{
    va_list args;
    char full_msg[4096];

    va_start(args, message);
    vsprintf(full_msg, message, args);
    va_end(args);

    fprintf(stderr, "\nError in %s:\n%s!\n", function_name, full_msg);
    fflush(stdout);
    exit(-1);
}

void printWarning(const char *function_name, const char *message, ...)
{
    va_list args;
    char full_msg[4096];

    va_start(args, message);
    vsprintf(full_msg, message, args);
    va_end(args);

    fprintf(stdout, "\nWarning in %s:\n%s!\n", function_name, full_msg);
}

/*
Find a string key in arguments and storage in argument variable
Input: argv and argc received in main, stringKey to match, argument variable to store match value
Output: argument variable with the match value or '-'
*/
char *parseArgs(char *argv[], int argc, const char *stringKey){
    
    for(int i=1; i < argc; i++){
        // if both strings are identical
        if(strcmp(argv[i], stringKey) == 0){
            return argv[i+1];
        }
    }
    return "-";
}