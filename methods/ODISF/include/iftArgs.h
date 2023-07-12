/*****************************************************************************\
* iftArgs.h
*
* AUTHOR  : Felipe Belem
* DATE    : 2021-02-24
* LICENSE : MIT License
* EMAIL   : felipe.belem@ic.unicamp.br
\*****************************************************************************/
#ifndef IFT_ARGS_H
#define IFT_ARGS_H

#ifdef __cplusplus
extern "C" {
#endif

#include "ift.h"

/*****************************************************************************\
*
*                                   STRUCTS
*
\*****************************************************************************/
/*
  Stores the references of the arguments given in the command line interface,
  simplifying the process of parsing its content for most practical standalone
  applications.
*/
typedef struct ift_args iftArgs;

/*****************************************************************************\
*
*                               PUBLIC FUNCTIONS
*
\*****************************************************************************/
//===========================================================================//
// CONSTRUCTOR & DESTRUCTOR
//===========================================================================//
/*
  Creates a new iftArgs object in which references to the array of arguments
  given in parameter. The user must also provide the exact number of arguments
  within the given array.
*/
iftArgs *iftCreateArgs
(const int argc, const char **argv);

/*
  Deallocates the object memory and sets it to NULL.
*/
void iftDestroyArgs
(iftArgs **args);

//===========================================================================//
// GETTER
//===========================================================================//
/*
  Gets the original string value associated to the argument token provided in 
  parameter. The token MUST NOT contain the default prefix.
*/
const char *iftGetArg
(const iftArgs *args, const char *token);

//===========================================================================//
// VERIFIERS
//===========================================================================//
/*  
  Verifies if the argument token exists in the container. The token MUST NOT
  contain the default prefix.
*/
bool iftExistArg
(const iftArgs *args, const char *token);

/*
  Verifies if the existing argument token has a value associated to it. The 
  token MUST NOT contain the default prefix.
*/
bool iftHasArgVal
(const iftArgs *args, const char *token);

#ifdef __cplusplus
}
#endif

#endif