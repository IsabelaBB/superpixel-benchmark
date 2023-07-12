/*****************************************************************************\
* iftArgs.c
*
* AUTHOR  : Felipe Belem
* DATE    : 2021-02-24
* LICENSE : MIT License
* EMAIL   : felipe.belem@ic.unicamp.br
\*****************************************************************************/
#include "iftArgs.h"

/*****************************************************************************\
*
*                                 DEFINITIONS
*
\*****************************************************************************/
#define _IFT_ARGS_PREFIX "--" // Default prefix for any token

/*****************************************************************************\
*
*                                   STRUCTS
*
\*****************************************************************************/
struct ift_args
{
  int argc; // Number of arguments given by the user
  const char **REF_ARGV; // Reference to the array of argument values
};


/*****************************************************************************\
*
*                              PRIVATE FUNCTIONS
*
\*****************************************************************************/
//===========================================================================//
// GENERAL & AUXILIARY
//===========================================================================//
/*
  Gets the argument index within the container, if it exists. If not, a 
  negative number (-1) indicates the inexistence of such token in the array.
*/
int _iftGetArgIdx
(const iftArgs *args, const char *token)
{
  #ifdef IFT_DEBUG //--------------------------------------------------------//
  assert(args != NULL);
  assert(token != NULL);
  #endif //------------------------------------------------------------------//
  bool found;
  int i, idx;
  char *full_token;

  full_token = iftConcatStrings(2, _IFT_ARGS_PREFIX, token);

  i = 0; found = false;
  while(i < args->argc && found == false)
  {
    found = iftCompareStrings(args->REF_ARGV[i], full_token);
    
    if(found != true) ++i;
  }
  free(full_token);

  if(i == args->argc) idx = -1;
  else idx = i;

  return idx;
}

/*****************************************************************************\
*
*                               PUBLIC FUNCTIONS
*
\*****************************************************************************/
//===========================================================================//
// CONSTRUCTOR & DESTRUCTOR
//===========================================================================//
iftArgs *iftCreateArgs
(const int argc, const char **argv)
{
  #ifdef IFT_DEBUG //--------------------------------------------------------//
  assert(argc >= 0);
  #endif //------------------------------------------------------------------//
  iftArgs *args;

  args = malloc(sizeof(iftArgs));
  assert(args != NULL);

  args->argc = argc;
  args->REF_ARGV = argv;

  return args;
}

void iftDestroyArgs
(iftArgs **args)
{
  #ifdef IFT_DEBUG //--------------------------------------------------------//
  assert(args != NULL && *args != NULL);
  #endif //------------------------------------------------------------------//
  free(*args);

  *args = NULL;
}

//===========================================================================//
// GETTER
//===========================================================================//
inline const char *iftGetArg
(const iftArgs *args, const char *token)
{
  #ifdef IFT_DEBUG //--------------------------------------------------------//
  assert(args != NULL);
  assert(token != NULL);
  assert(iftExistArg(args, token));
  assert(iftHasArgVal(args, token));
  #endif //------------------------------------------------------------------//

  return args->REF_ARGV[_iftGetArgIdx(args, token) + 1];
}


//===========================================================================//
// VERIFIERS
//===========================================================================//
inline bool iftExistArg
(const iftArgs *args, const char *token)
{
  #ifdef IFT_DEBUG //--------------------------------------------------------//
  assert(args != NULL);
  assert(token != NULL);
  #endif //------------------------------------------------------------------//
  
  return _iftGetArgIdx(args, token) >= 0;
}

bool iftHasArgVal
(const iftArgs *args, const char *token)
{
  #ifdef IFT_DEBUG //--------------------------------------------------------//
  assert(args != NULL);
  assert(token != NULL);
  assert(iftExistArg(args, token) == true);
  #endif //------------------------------------------------------------------//
  const char *VAL;
  bool has_val;
  int idx;

  idx = _iftGetArgIdx(args, token);
  VAL = args->REF_ARGV[idx + 1]; // The subsequent must be the value
  
  // If it starts with the prefix, it is a token, not a value
  if(idx == args->argc - 1 || iftStartsWith(VAL, _IFT_ARGS_PREFIX) == true) 
  {
    has_val = false;
  }
  else has_val = true;

  return has_val;
}