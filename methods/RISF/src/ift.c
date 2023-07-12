#include "ift.h"

// ---------- iftGQueue.c start 

iftGQueue *iftCreateGQueue(int nbuckets, int nelems, int *value)
{
    iftGQueue *Q=NULL;
    
    Q = (iftGQueue *) iftAlloc(1, sizeof(iftGQueue));
    
    if (Q != NULL)
    {
        Q->C.first = (int *)iftAlloc((nbuckets+1), sizeof(int));
        Q->C.last  = (int *)iftAlloc((nbuckets+1), sizeof(int));
        Q->C.nbuckets = nbuckets;
        if ( (Q->C.first != NULL) && (Q->C.last != NULL) )
        {
            Q->L.elem = (iftGQNode *)iftAlloc(nelems, sizeof(iftGQNode));
            Q->L.nelems = nelems;
            Q->L.value   = value;
            if (Q->L.elem != NULL)
            {
                iftResetGQueue(Q);
            }
            else
                iftError(MSG_MEMORY_ALLOC_ERROR, "iftCreateGQueue");
        }
        else
            iftError(MSG_MEMORY_ALLOC_ERROR, "iftCreateGQueue");
    }
    else
        iftError(MSG_MEMORY_ALLOC_ERROR, "iftCreateGQueue");
    
    /* default */
    
    iftSetTieBreak(Q,FIFOBREAK);
    iftSetRemovalPolicy(Q,MINVALUE);
    
    return(Q);
}

void iftDestroyGQueue(iftGQueue **Q)
{
    iftGQueue *aux=*Q;
    
    if (aux != NULL)
    {
        if (aux->C.first != NULL) iftFree(aux->C.first);
        if (aux->C.last  != NULL) iftFree(aux->C.last);
        if (aux->L.elem  != NULL) iftFree(aux->L.elem);
        iftFree(aux);
        *Q = NULL;
    }
}

int iftEmptyGQueue(iftGQueue *Q)
{
    int last,current;
    
    if (Q->C.removal_policy == MINVALUE)
        current=iftSafeMod(Q->C.minvalue,Q->C.nbuckets);
    else
        current=Q->C.nbuckets - 1 - (iftSafeMod(Q->C.maxvalue,Q->C.nbuckets));
    
    if (Q->C.first[current] != IFT_NIL)
        return 0;
    
    last = current;
    
    current = iftSafeMod(current + 1, Q->C.nbuckets);
    
    while ((Q->C.first[current] == IFT_NIL) && (current != last))
    {
        current = iftSafeMod(current + 1, Q->C.nbuckets);
    }
    
    if (Q->C.first[current] == IFT_NIL)
    {
        if (Q->C.first[Q->C.nbuckets] == IFT_NIL)
        {
            //Changed by Falcao and Nikolas
            // iftResetGQueue(Q);
            return(1);
        }
    }
    
    return (0);
}

void iftInsertGQueue(iftGQueue **Q, int elem)
{
    int bucket,minvalue=(*Q)->C.minvalue,maxvalue=(*Q)->C.maxvalue;
    
    if (((*Q)->L.value[elem] == IFT_INFINITY_INT) || ((*Q)->L.value[elem] == IFT_INFINITY_INT_NEG))
        bucket=(*Q)->C.nbuckets;
    else
    {
        if ((*Q)->L.value[elem] < minvalue)
            minvalue = (*Q)->L.value[elem];
        if ((*Q)->L.value[elem] > maxvalue)
            maxvalue = (*Q)->L.value[elem];
        if ((maxvalue-minvalue) > ((*Q)->C.nbuckets-1))
        {
            (*Q) = iftGrowGQueue(Q,2*(maxvalue-minvalue)+1);
            fprintf(stdout,"Warning: Doubling queue size\n");
        }
        if ((*Q)->C.removal_policy==MINVALUE)
        {
            bucket=iftSafeMod((*Q)->L.value[elem],(*Q)->C.nbuckets);
        }
        else
        {
            bucket=(*Q)->C.nbuckets-1-(iftSafeMod((*Q)->L.value[elem],(*Q)->C.nbuckets));
        }
        (*Q)->C.minvalue = minvalue;
        (*Q)->C.maxvalue = maxvalue;
    }
    if ((*Q)->C.first[bucket] == IFT_NIL)
    {
        (*Q)->C.first[bucket]   = elem;
        (*Q)->L.elem[elem].prev = IFT_NIL;
    }
    else
    {
        (*Q)->L.elem[(*Q)->C.last[bucket]].next = elem;
        (*Q)->L.elem[elem].prev = (*Q)->C.last[bucket];
    }
    
    (*Q)->C.last[bucket]     = elem;
    (*Q)->L.elem[elem].next  = IFT_NIL;
    (*Q)->L.elem[elem].color = IFT_GRAY;
}

int iftRemoveGQueue(iftGQueue *Q)
{
    int elem= IFT_NIL, next, prev;
    int last, current;
    
    if (Q->C.removal_policy==MINVALUE)
        current=iftSafeMod(Q->C.minvalue,Q->C.nbuckets);
    else
        current=Q->C.nbuckets-1-iftSafeMod(Q->C.maxvalue,Q->C.nbuckets);
    
    /** moves to next element **/
    
    if (Q->C.first[current] == IFT_NIL)
    {
        last = current;
        
        current = iftSafeMod(current + 1, Q->C.nbuckets);
        
        while ((Q->C.first[current] == IFT_NIL) && (current != last))
        {
            current = iftSafeMod(current + 1, Q->C.nbuckets);
        }
        
        if (Q->C.first[current] != IFT_NIL)
        {
            if (Q->C.removal_policy==MINVALUE)
                Q->C.minvalue = Q->L.value[Q->C.first[current]];
            else
                Q->C.maxvalue = Q->L.value[Q->C.first[current]];
        }
        else
        {
            if (Q->C.first[Q->C.nbuckets] != IFT_NIL)
            {
                current = Q->C.nbuckets;
                if (Q->C.removal_policy==MINVALUE)
                    Q->C.minvalue = Q->L.value[Q->C.first[current]];
                else
                    Q->C.maxvalue = Q->L.value[Q->C.first[current]];
            }
            else
            {
                iftError("iftGQueue is empty\n", "iftRemoveGQueue");
            }
        }
    }
    
    if (Q->C.tiebreak == LIFOBREAK)
    {
        elem = Q->C.last[current];
        prev = Q->L.elem[elem].prev;
        if (prev == IFT_NIL)           /* there was a single element in the list */
        {
            Q->C.last[current] = Q->C.first[current]  = IFT_NIL;
        }
        else
        {
            Q->C.last[current]   = prev;
            Q->L.elem[prev].next = IFT_NIL;
        }
    }
    else   /* Assume FIFO policy for breaking ties */
    {
        elem = Q->C.first[current];
        next = Q->L.elem[elem].next;
        if (next == IFT_NIL)           /* there was a single element in the list */
        {
            Q->C.first[current] = Q->C.last[current]  = IFT_NIL;
        }
        else
        {
            Q->C.first[current] = next;
            Q->L.elem[next].prev = IFT_NIL;
        }
    }
    
    Q->L.elem[elem].color = IFT_BLACK;
    
    return elem;
}

void iftRemoveGQueueElem(iftGQueue *Q, int elem)
{
    int prev,next,bucket;
    
    if ((Q->L.value[elem] == IFT_INFINITY_INT) || (Q->L.value[elem] == IFT_INFINITY_INT_NEG))
        bucket = Q->C.nbuckets;
    else
    {
        if (Q->C.removal_policy == MINVALUE)
            bucket = iftSafeMod(Q->L.value[elem],Q->C.nbuckets);
        else
            bucket = Q->C.nbuckets-1-iftSafeMod(Q->L.value[elem],Q->C.nbuckets);
    }
    
    prev = Q->L.elem[elem].prev;
    next = Q->L.elem[elem].next;
    
    /* if elem is the first element */
    if (Q->C.first[bucket] == elem)
    {
        Q->C.first[bucket] = next;
        if (next == IFT_NIL) /* elem is also the last one */
            Q->C.last[bucket] = IFT_NIL;
        else
            Q->L.elem[next].prev = IFT_NIL;
    }
    else    /* elem is in the middle or it is the last */
    {
        Q->L.elem[prev].next = next;
        if (next == IFT_NIL) /* if it is the last */
            Q->C.last[bucket] = prev;
        else
            Q->L.elem[next].prev = prev;
    }
    
    Q->L.elem[elem].color = IFT_WHITE;
    
}

void iftResetGQueue(iftGQueue *Q)
{
    int i;
    
    Q->C.minvalue = IFT_INFINITY_INT;
    Q->C.maxvalue = IFT_INFINITY_INT_NEG;
    /* No need for that, since the programmer might have changed them  */
    //    iftSetTieBreak(Q,FIFOBREAK);
    //    iftSetRemovalPolicy(Q,MINVALUE);
    for (i=0; i < Q->C.nbuckets+1; i++)
        Q->C.first[i]=Q->C.last[i]= IFT_NIL;
    
    for (i=0; i < Q->L.nelems; i++)
    {
        Q->L.elem[i].next =  Q->L.elem[i].prev = IFT_NIL;
        Q->L.elem[i].color = IFT_WHITE;
    }
    
}

iftGQueue *iftGrowGQueue(iftGQueue **Q, int nbuckets)
{
    iftGQueue *Q1=iftCreateGQueue(nbuckets,(*Q)->L.nelems,(*Q)->L.value);
    int i,bucket;
    
    Q1->C.minvalue  = (*Q)->C.minvalue;
    Q1->C.maxvalue  = (*Q)->C.maxvalue;
    Q1->C.tiebreak = (*Q)->C.tiebreak;
    Q1->C.removal_policy = (*Q)->C.removal_policy;
    for (i=0; i<(*Q)->C.nbuckets; i++)
        if ((*Q)->C.first[i] != IFT_NIL)
        {
            bucket = iftSafeMod((*Q)->L.value[(*Q)->C.first[i]],Q1->C.nbuckets);
            Q1->C.first[bucket] = (*Q)->C.first[i];
            Q1->C.last[bucket]  = (*Q)->C.last[i];
        }
    if ((*Q)->C.first[(*Q)->C.nbuckets] != IFT_NIL)
    {
        bucket = Q1->C.nbuckets;
        Q1->C.first[bucket] = (*Q)->C.first[(*Q)->C.nbuckets];
        Q1->C.last[bucket]  = (*Q)->C.last[(*Q)->C.nbuckets];
    }
    
    for (i=0; i < (*Q)->L.nelems; i++)
        Q1->L.elem[i]  = (*Q)->L.elem[i];
    
    iftDestroyGQueue(Q);
    return(Q1);
}

// ---------- iftGQueue.c end 
// ---------- iftBasicDataTypes.c start 

iftGVal iftInitBoolGVal(bool val) 
{
    iftGVal gv;
    gv.type = IFT_BOOL_TYPE;
    gv.bool_val = val;
    return gv;
}

iftGVal iftInitCharGVal(char val) 
{
    iftGVal gv;
    gv.type = IFT_CHAR_TYPE;
    gv.char_val = val;
    return gv;
}

iftGVal iftInitUCharGVal(uchar val) 
{
    iftGVal gv;
    gv.type = IFT_UCHAR_TYPE;
    gv.uchar_val = val;
    return gv;
}

iftGVal iftInitStrGVal(const char *val) 
{
    char *str = (char *) calloc(strlen(val)+1, sizeof(char));
    strcpy(str, val);
    iftGVal gv;
    gv.type = IFT_STR_TYPE;
    gv.str_val = str;
    return gv;
}

iftGVal iftInitLongGVal(long val) 
{
    iftGVal gv;
    gv.type = IFT_LONG_TYPE;
    gv.long_val = val;
    return gv;
}

iftGVal iftInitULongGVal(ulong val) 
{
    iftGVal gv;
    gv.type = IFT_ULONG_TYPE;
    gv.ulong_val = val;
    return gv;
}

iftGVal iftInitDblGVal(double val) 
{
    iftGVal gv;
    gv.type = IFT_DBL_TYPE;
    gv.dbl_val = val;
    return gv;
}

iftGVal iftInitIntArrayGVal(iftIntArray *array) 
{
    iftGVal gv;
    gv.type = IFT_INT_ARRAY_TYPE;
    gv.int_array_val = array;
    return gv;
}

iftGVal iftInitDblArrayGVal(iftDblArray *array) 
{
    iftGVal gv;
    gv.type = IFT_DBL_ARRAY_TYPE;
    gv.dbl_array_val = array;
    return gv;
}

iftGVal iftInitStrArrayGVal(iftStrArray *array) 
{
    iftGVal gv;
    gv.type = IFT_STR_ARRAY_TYPE;
    gv.str_array_val = array;
    return gv;
}

iftGVal iftInitIntMatrixGVal(iftIntMatrix *mat) 
{
    iftGVal gv;
    gv.type = IFT_INT_MATRIX_TYPE;
    gv.int_matrix_val = mat;
    return gv;
}

iftGVal iftInitDblMatrixGVal(iftMatrix *mat) 
{
    iftGVal gv;
    gv.type = IFT_DBL_MATRIX_TYPE;
    gv.dbl_matrix_val = mat;
    return gv;
}

iftGVal iftInitStrMatrixGVal(iftStrMatrix *mat) 
{
    iftGVal gv;
    gv.type = IFT_STR_MATRIX_TYPE;
    gv.str_matrix_val = mat;
    return gv;
}

iftGVal iftInitDictGVal(iftDict *dict) 
{
    iftGVal gv;
    gv.type = IFT_DICT_TYPE;
    gv.dict_val = dict;
    return gv;
}

iftGVal iftInitPtrGVal(void *val) 
{
    iftGVal gv;
    gv.type = IFT_PTR_TYPE;
    gv.ptr_val = val;
    return gv;
}

void iftFreeGVal(iftGVal gv) 
{
    if (gv.type == IFT_STR_TYPE)
        free(gv.str_val);
    else if (gv.type == IFT_INT_ARRAY_TYPE)
        iftDestroyIntArray(&gv.int_array_val);
    else if (gv.type == IFT_DBL_ARRAY_TYPE)
        iftDestroyDblArray(&gv.dbl_array_val);
    else if (gv.type == IFT_STR_ARRAY_TYPE)
        iftDestroyStrArray(&gv.str_array_val);
    else if (gv.type == IFT_INT_MATRIX_TYPE)
        iftDestroyIntMatrix(&gv.int_matrix_val);
    else if (gv.type == IFT_DBL_MATRIX_TYPE)
        iftDestroyMatrix(&gv.dbl_matrix_val);
    else if (gv.type == IFT_STR_MATRIX_TYPE)
        iftDestroyStrMatrix(&gv.str_matrix_val);

    gv.ptr_val = NULL;
}

const char *iftCDataTypeToString(iftCDataType datatype) 
{
    switch(datatype) {
        case IFT_UNTYPED:
            return "untyped";
        case IFT_BOOL_TYPE:
            return "boolean";
        case IFT_CHAR_TYPE:
            return "char";
        case IFT_UCHAR_TYPE:
            return "unsigned char";
        case IFT_STR_TYPE:
            return "string (char*)";
        case IFT_INT_TYPE:
            return "int";
        case IFT_UINT_TYPE:
            return "unsigned int";
        case IFT_LONG_TYPE:
            return "long";
        case IFT_ULONG_TYPE:
            return "unsigned long";
        case IFT_FLT_TYPE:
            return "float";
        case IFT_DBL_TYPE:
            return "double";
        case IFT_INT_ARRAY_TYPE:
            return "iftIntArray";
        case IFT_DBL_ARRAY_TYPE:
            return "iftDblArray";
        case IFT_STR_ARRAY_TYPE:
            return "iftStrArray";
        case IFT_INT_MATRIX_TYPE:
            return "iftIntMatrix";
        case IFT_DBL_MATRIX_TYPE:
            return "iftMatrix";
        case IFT_STR_MATRIX_TYPE:
            return "iftStrMatrix";
        case IFT_DICT_TYPE:
            return "iftDict*";
        case IFT_PTR_TYPE:
            return "void*";
        default:
            return "unknown type"; // just to avoid compilation warnings
    }
}

char *iftGValToString(iftGVal gval) 
{
    char *str = (char *) calloc(IFT_STR_DEFAULT_SIZE, sizeof(char));

    switch(gval.type) {
        case IFT_BOOL_TYPE:
            sprintf(str, "%s", gval.bool_val ? "true" : "false");
            break;
        case IFT_CHAR_TYPE:
            sprintf(str, "\'%c\'", gval.char_val);
            break;
        case IFT_UCHAR_TYPE:
            sprintf(str, "\'%c\'", gval.char_val);
            break;
        case IFT_STR_TYPE:
            sprintf(str, "\"%s\"", gval.str_val);
            break;
        case IFT_LONG_TYPE:
            sprintf(str, "%ld", gval.long_val);
            break;
        case IFT_ULONG_TYPE:
            sprintf(str, "%lu", gval.ulong_val);
            break;
        case IFT_DBL_TYPE:
            sprintf(str, "%lf", gval.dbl_val);
            break;
        case IFT_PTR_TYPE:
            if (gval.ptr_val == 0x0) {
                sprintf(str, "null");
            }
            else sprintf(str, "%p", gval.ptr_val);
            break;
        default:
            strcpy(str, "");
            break;
    }

    return str;
}

bool iftCompareGVal(iftGVal val1, iftGVal val2) 
{
    if (val1.type == val2.type) {
        switch(val1.type) {
            case IFT_BOOL_TYPE:
                return (val1.bool_val == val2.bool_val);
            case IFT_CHAR_TYPE:
                return (val1.char_val == val2.char_val);
            case IFT_UCHAR_TYPE:
                return (val1.uchar_val == val2.uchar_val);
            case IFT_STR_TYPE:
                return (strcmp(val1.str_val, val2.str_val) == 0);
            case IFT_LONG_TYPE:
                return (val1.long_val == val2.long_val);
            case IFT_ULONG_TYPE:
                return (val1.ulong_val == val2.ulong_val);
            case IFT_DBL_TYPE:
                return (val1.dbl_val == val2.dbl_val);
            case IFT_PTR_TYPE:
                return (val1.ptr_val == val2.ptr_val);
            default:
                return false; // IFT_UNTYPED
        }
    }
    else return false;
}

iftGVal iftCopyGVal(iftGVal gval) 
{
    iftGVal out;
    out.type = gval.type;

    switch (gval.type) {
        case IFT_BOOL_TYPE:
            out = iftInitBoolGVal(gval.bool_val);
            break;
        case IFT_CHAR_TYPE:
            out = iftInitCharGVal(gval.char_val);
            break;
        case IFT_UCHAR_TYPE:
            out = iftInitUCharGVal(gval.uchar_val);
            break;
        case IFT_STR_TYPE:
            out = iftInitStrGVal(gval.str_val);
            break;
        case IFT_INT_TYPE:
        case IFT_LONG_TYPE:
            out = iftInitLongGVal(gval.long_val);
            break;
        case IFT_UINT_TYPE:
        case IFT_ULONG_TYPE:
            out = iftInitULongGVal(gval.ulong_val);
            break;
        case IFT_FLT_TYPE:
        case IFT_DBL_TYPE:
            out = iftInitDblGVal(gval.dbl_val);
            break;
        case IFT_INT_ARRAY_TYPE:
            out.int_array_val = iftCreateIntArray(gval.int_array_val->n);
            iftCopyIntArray(out.int_array_val->val, gval.int_array_val->val, gval.int_array_val->n);
            break;
        case IFT_DBL_ARRAY_TYPE:
            out.dbl_array_val = iftCopyDblArray(gval.dbl_array_val->val, gval.dbl_array_val->n);
            break;
        case IFT_STR_ARRAY_TYPE:
            out.str_array_val = iftCopyStrArray(gval.str_array_val->val, gval.str_array_val->n);
            break;
        case IFT_INT_MATRIX_TYPE:
            out.int_matrix_val = iftCopyIntMatrix(gval.int_matrix_val->val, gval.int_matrix_val->nrows, gval.int_matrix_val->ncols);
            break;
        case IFT_DBL_MATRIX_TYPE:
            out.dbl_matrix_val = iftCopyMatrix(gval.dbl_matrix_val);
            break;
        case IFT_STR_MATRIX_TYPE:
            out.str_matrix_val = iftCopyStrMatrix(gval.str_matrix_val->val, gval.str_matrix_val->nrows, gval.str_matrix_val->ncols);
            break;
        case IFT_DICT_TYPE:
            out.dict_val = iftCopyDict(gval.dict_val);
            break;
        case IFT_UNTYPED:
        case IFT_PTR_TYPE:
        default:
            out.ptr_val = gval.ptr_val;
            break;
    }

    return out;
}

bool iftGetBoolVal(iftGVal gval) 
{
    return gval.bool_val;
}

char iftGetCharVal(iftGVal gval) 
{
    return gval.char_val;
}

uchar iftGetUCharVal(iftGVal gval) 
{
    return gval.uchar_val;
}

char *iftGetStrVal(iftGVal gval) 
{
    return gval.str_val;
}

const char *iftGetConstStrVal(iftGVal gval) 
{
    return gval.str_val;
}

long iftGetLongVal(iftGVal gval) 
{
    return gval.long_val;
}

ulong iftGetULongVal(iftGVal gval) 
{
    return gval.ulong_val;
}

double iftGetDblVal(iftGVal gval) 
{
    return gval.dbl_val;
}

iftIntArray *iftGetIntArrayVal(iftGVal gval) 
{
    return gval.int_array_val;
}

iftDblArray *iftGetDblArrayVal(iftGVal gval) 
{
    return gval.dbl_array_val;
}

iftStrArray *iftGetStrArrayVal(iftGVal gval) 
{
    return gval.str_array_val;
}

iftIntMatrix *iftGetIntMatrixVal(iftGVal gval)
{
    return gval.int_matrix_val;
}

iftStrMatrix *iftGetStrMatrixVal(iftGVal gval) 
{
    return gval.str_matrix_val;
}

iftDict *iftGetDictVal(iftGVal gval) 
{
    return gval.dict_val;
}

void *iftGetPtrVal(iftGVal gval) 
{
    return gval.ptr_val;
}

void iftCopyVoxel(iftVoxel *src, iftVoxel *dst) 
{
    (*dst).x = (*src).x;
    (*dst).y = (*src).y;
    (*dst).z = (*src).z;
}

// ---------- iftBasicDataTypes.c end
// ---------- iftIntArray.c start

iftIntArray *iftCreateIntArray(long n) 
{
    iftIntArray *iarr = (iftIntArray*) iftAlloc(1, sizeof(iftIntArray));
    
    iarr->n = n;
    iarr->val = iftAllocIntArray(n);
    
    return iarr;
}

void iftDestroyIntArray(iftIntArray **iarr) 
{
    if (iarr != NULL && *iarr != NULL) {
        iftIntArray *iarr_aux = *iarr;
        
        if (iarr_aux->val != NULL)
            iftFree(iarr_aux->val);
        iftFree(iarr_aux);
        *iarr = NULL;
    }
}

void iftShuffleIntArray(int* array, int n) 
{
    // Start from the last element and swap one by one. We don't
    // need to run for the first element that's why i > 0
    int j;
    for (int i = n-1; i > 0; i--)
    {
        // Pick a random index from 0 to i
        j = rand() % (i+1);
        // Swap arr[i] with the element at random index
        iftSwap(array[i],array[j]);
    }
}

// ---------- iftIntArray.c end
// ---------- iftFloatArray.c start

iftFloatArray *iftCreateFloatArray(long n) 
{
    iftFloatArray *darr = (iftFloatArray *) iftAlloc(1, sizeof(iftFloatArray));
    
    darr->n   = n;
    darr->val = iftAllocFloatArray(n);
    
    return darr;
}

void iftDestroyFloatArray(iftFloatArray **darr) 
{
    if (darr != NULL && *darr != NULL) {
        iftFloatArray *darr_aux = *darr;
        
        if (darr_aux->val != NULL)
            iftFree(darr_aux->val);
        iftFree(darr_aux);
        *darr = NULL;
    }
}

// ---------- iftFloatArray.c end 
// ---------- iftDblArray.c start 

iftDblArray *iftCreateDblArray(long n) 
{
    iftDblArray *darr = (iftDblArray *) iftAlloc(1, sizeof(iftDblArray));
    
    darr->n   = n;
    darr->val = iftAllocDoubleArray(n);
    
    return darr;
}

iftDblArray *iftCopyDblArray(const double* array, long n) 
{
    iftDblArray* out = iftCreateDblArray(n);
    
    memcpy(out->val, array, n*sizeof(double));
    
    return out;
}

void iftDestroyDblArray(iftDblArray **darr) 
{
    
    if(darr == NULL){
        return;
    }
    if(*darr == NULL){
        return;
    }
    
    
    if (darr != NULL && *darr != NULL) {
        iftDblArray *darr_aux = *darr;
        
        if (darr_aux->val != NULL)
            iftFree(darr_aux->val);
        iftFree(darr_aux);
        *darr = NULL;
        
    }
}

// ---------- iftDblArray.c end 
// ---------- iftStrArray.c start 

iftStrArray *iftCreateStrArray(long n) 
{
    iftStrArray *sarr = (iftStrArray *) iftAlloc(1, sizeof(iftStrArray));
    
    sarr->n   = n;
    sarr->val = (char**) iftAlloc(n, sizeof(char*));
    for (int i = 0; i < n; i++) {
        sarr->val[i] = iftAllocCharArray(2048);
    }
    
    return sarr;
}

iftStrArray *iftCopyStrArray(char **src_arr, long n) 
{
    iftStrArray *copy = iftCreateStrArray(n);
    
    for (int i = 0; i < copy->n; i++) {
        iftFree(copy->val[i]);
        copy->val[i] = iftCopyString(src_arr[i]);
    }
    
    return copy;
}

void iftDestroyStrArray(iftStrArray **sarr) 
{
    if (sarr != NULL && *sarr != NULL) {
        iftStrArray *sarr_aux = *sarr;
        
        if (sarr_aux->val != NULL) {
            for (int i = 0; i < sarr_aux->n; i++) {
                iftFree(sarr_aux->val[i]);
            }
            iftFree(sarr_aux->val);
        }
        iftFree(sarr_aux);
        *sarr = NULL;
    }
}

// ---------- iftStrArray.c end 
// ---------- iftIntMatrix.c start

iftIntMatrix *iftCreateIntMatrix(int ncols, int nrows) 
{
    iftIntMatrix *iM = (iftIntMatrix *) iftAlloc(1, sizeof(iftIntMatrix));
    
    iM->ncols = ncols;
    iM->nrows = nrows;
    
    iM->tbrow = iftAllocIntArray(nrows);
    for (int r = 0; r < nrows; r++) {
        iM->tbrow[r] = r * ncols;
    }
    
    iM->n = ncols * nrows;
    iM->val = iftAllocIntArray(iM->n);
    
    return iM;
}

iftIntMatrix *iftCopyIntMatrix(int *val, int nrows, int ncols) 
{
    iftIntMatrix *mat = iftCreateIntMatrix(ncols, nrows);
    
    for (int i = 0; i < mat->n; i++)
        mat->val[i] = val[i];
    
    return mat;
}

void iftDestroyIntMatrix(iftIntMatrix **iM) 
{
    iftIntMatrix *aux = *iM;
    
    if (aux != NULL) {
        if (aux->val != NULL)
            iftFree(aux->val);
        iftFree(aux->tbrow);
        iftFree(aux);
        *iM = NULL;
    }
}

// ---------- iftIntMatrix.c end 
// ---------- iftCommon.c start 

int iftFArgmax(const float *x, int n)
{
    int i, best_i = IFT_NIL;
    float x_max = IFT_INFINITY_FLT_NEG;

    for(i = 0; i < n; i++)
    {
        if(x[i] > x_max)
        {
            best_i = i;
            x_max = x[i];
        }
    }

    return best_i;
}

float iftFeatDistance(float *A, float *B, int n)
{
    float dist=0.0;
    int    i;
    for (i=0; i < n; i++)
        dist += (A[i]-B[i])*(A[i]-B[i]);

    return(sqrt(dist));
}

int iftRandomInteger (int low, int high)
{
    int k;
    double d;

    d = (double) rand () / ((double) RAND_MAX + 0.5);
    k =  iftMin((int)(d * (high - low + 1.0) ) + low, high);

    return k;
}

void iftSkipComments(FILE *fp)
{
    //skip for comments
    while (fgetc(fp) == '#') {
        while (fgetc(fp) != '\n');
    }
    fseek(fp,-1,SEEK_CUR);

}

double iftLog(double val, double base) 
{
    return (log(val) / log(base));
}

long iftNormalizationValue(long maxval) 
{
    long norm_val = 1;
    
    if (maxval < 0)
        iftError("Input value %ld < 0", "iftNormalizationValue", maxval);
    else if (maxval <= 1)
        norm_val = 1;
    else if (maxval <= 255)
        norm_val = 255;
    else if (maxval <= 4095)
        norm_val = 4095;
    else if (maxval <= 65535)
        norm_val = 65535;
    else if (maxval <= 4294967295)
        norm_val = 4294967295;
    else iftError("Invalid maxval number %ld with number of bits > 32. It only supports values within [0, 2ˆn_bits -1], " \
                  "where n_bits in {1, 8, 12, 16, 32}", "iftNormalizationValue", maxval);
    
    return norm_val;
}

int iftAlmostZero(double x)
{
    return (x >= -IFT_EPSILON) && (x <= IFT_EPSILON);
}

void iftRandomSeed(unsigned int seed)
{

    srand(seed);

}

inline int iftSafeMod(int a, int n)
{
    int r = a % n;

    return (r >= 0) ? r : n+r;
}

bool iftIsPrime(long n) 
{
    if (n < 0)
        iftError("Number is Negative (not Natural): %ld... Try Natural Numbers", "iftIsPrime", n);

    if (n <= 1)
        return false;
    else
    if ((n == 2) || (n == 3))
        return true;
    else
    if (((n % 2) == 0) || ((n % 3) == 0))
        return false;

    long sqrt_of_n = (long) sqrt(n); // floor of sqrt(n)

    for (long d = 5; d <= sqrt_of_n; d = d + 6)
        if (((n % d) == 0) || ((n % (d+2)) == 0))
            return false;

    return true;
}

void iftNormalizeFeatures(float *feats, int nelems)
{
    int i;
    float maximum = IFT_INFINITY_FLT_NEG;
    for (i = 0; i < nelems; ++i)
    {
        maximum = iftMax(maximum, feats[i]);
    }

    for (i = 0; i < nelems; ++i)
    {
        feats[i]/=maximum;
    }
}

timer *iftTic()
{
    timer *tic=NULL;
    tic = (timer *)iftAlloc(1, sizeof(timer));
    gettimeofday(tic,NULL);
    return(tic);
}

timer *iftToc()
{
    timer *toc=NULL;
    toc = (timer *)iftAlloc(1, sizeof(timer));
    gettimeofday(toc,NULL);
    return(toc);
}

float iftCompTime(timer *tic, timer *toc)
{
    float t=0.0;
    if ((tic!=NULL)&&(toc!=NULL)){
        t = (toc->tv_sec-tic->tv_sec)*1000.0 +
            (toc->tv_usec-tic->tv_usec)*0.001;
        iftFree(tic);iftFree(toc);
    }
    return(t);
}

// ---------- iftCommon.c end
// ---------- iftAdjacency.c start 

iftAdjRel  *iftCreateAdjRel(int n)
{
  iftAdjRel *A=(iftAdjRel *)iftAlloc(1,sizeof(iftAdjRel ));

  A->dx = (int *)iftAllocIntArray(n);
  A->dy = (int *)iftAllocIntArray(n);
  A->dz = (int *)iftAllocIntArray(n);
  A->n  = n;

  return(A);
}

void     iftDestroyAdjRel(iftAdjRel **A)
{
  iftAdjRel *aux = *A;

  if (aux != NULL){
    if (aux->dx != NULL) iftFree(aux->dx);
    if (aux->dy != NULL) iftFree(aux->dy);
    if (aux->dz != NULL) iftFree(aux->dz);
    iftFree(aux);
    *A = NULL;
  }
}

iftAdjRel *iftSpheric(float r) 
{
    int r0 = (int) r;
    float r2 = (int) (r*r + 0.5);

    int n = 0;
    for (int dz = -r0; dz <= r0; dz++)
        for (int dy = -r0; dy <= r0; dy++)
            for (int dx =- r0; dx <= r0; dx++)
                if ( ((dx*dx) + (dy*dy) + (dz*dz)) <= r2)
                    n++;


    iftAdjRel *A = iftCreateAdjRel(n);

    int i = 0;
    int i0 = 0;
    for (int dz = -r0; dz <= r0; dz++)
        for (int dy = -r0; dy <= r0; dy++)
            for (int dx = -r0; dx <= r0; dx++)
                if ( ((dx*dx) + (dy*dy) + (dz*dz)) <= r2) {
                    A->dx[i] = dx;
                    A->dy[i] = dy;
                    A->dz[i] = dz;
                
                if ((dx == 0) && (dy == 0) && (dz == 0))
                    i0 = i;
                i++;
            }

    // shift to right and place central voxel at first
    for (int i = i0; i > 0; i--) {
        int dx = A->dx[i];
        int dy = A->dy[i];
        int dz = A->dz[i];
        A->dx[i] = A->dx[i-1];
        A->dy[i] = A->dy[i-1];
        A->dz[i] = A->dz[i-1];
        A->dx[i-1] = dx;
        A->dy[i-1] = dy;
        A->dz[i-1] = dz;
    }


    // sort by radius, so the 6 closest neighbors will come first
    float *dr = iftAllocFloatArray(A->n);
    for (int i = 0; i < A->n; i++)
        dr[i] = A->dx[i]*A->dx[i] + A->dy[i]*A->dy[i] + A->dz[i]*A->dz[i];

    iftIntArray *idxs = iftIntRange(0, A->n-1, 1);
    iftFQuickSort(dr, idxs->val, 0, A->n-1, IFT_INCREASING);
    iftAdjRel *Asort = iftCreateAdjRel(A->n);

    for (int i = 0; i < A->n; i++) {
        int idx = idxs->val[i];
        Asort->dx[i] = A->dx[idx];
        Asort->dy[i] = A->dy[idx];
        Asort->dz[i] = A->dz[idx];
    }

    iftFree(dr);
    iftDestroyIntArray(&idxs);
    iftDestroyAdjRel(&A);

    return Asort;
}

iftAdjRel *iftCircular(float r)
{
    int r0 = (int) r;
    float r2 = (int) (r*r + 0.5);
    
    int n = 0;
    for (int dy = -r0; dy <= r0; dy++)
        for (int dx = -r0; dx <= r0; dx++)
            if (((dx*dx) + (dy*dy)) <= r2)
                n++;

    iftAdjRel *A = iftCreateAdjRel(n);
    int i = 0;
    int i0 = 0;
    for (int dy = -r0; dy <= r0; dy++)
        for (int dx = -r0; dx <= r0; dx++)
            if (((dx*dx) + (dy*dy)) <= r2) {
                A->dx[i] = dx;
                A->dy[i] = dy;
                A->dz[i] = 0;

                if ((dx==0) && (dy==0))
                    i0 = i;
                i++;
            }

    // shift to right and place central pixel at first
    for (int i = i0; i > 0; i--) {
        int dx     = A->dx[i];
        int dy     = A->dy[i];
        A->dx[i]   = A->dx[i-1];
        A->dy[i]   = A->dy[i-1];
        A->dx[i-1] = dx;
        A->dy[i-1] = dy;
    }


    // sort by radius, so the 4 closest neighbors will come first
    float *dr = iftAllocFloatArray(A->n);
    for (int i = 0; i < A->n; i++)
        dr[i] = A->dx[i]*A->dx[i] + A->dy[i]*A->dy[i];


    iftIntArray *idxs = iftIntRange(0, A->n-1, 1);
    iftFQuickSort(dr, idxs->val, 0, A->n-1, IFT_INCREASING);
    iftAdjRel *Asort = iftCreateAdjRel(A->n);

    for (int i = 0; i < A->n; i++) {
        int idx = idxs->val[i];
        Asort->dx[i] = A->dx[idx];
        Asort->dy[i] = A->dy[idx];
        Asort->dz[i] = A->dz[idx];
    }

    iftFree(dr);
    iftDestroyIntArray(&idxs);
    iftDestroyAdjRel(&A);

    return Asort;
}

iftAdjRel *iftCopyAdjacency(const iftAdjRel *A) 
{
  iftAdjRel *B = iftCreateAdjRel(A->n);
  int i;
  for (i=0; i < A->n; i++) {
    B->dx[i] = A->dx[i];
    B->dy[i] = A->dy[i];
    B->dz[i] = A->dz[i];
  }

  return(B);
}

inline iftVoxel iftGetAdjacentVoxel(const iftAdjRel *A, iftVoxel u, int adj)
{
  iftVoxel v;

  v.x = u.x + A->dx[adj];
  v.y = u.y + A->dy[adj];
  v.z = u.z + A->dz[adj];

  return(v);
}

void iftMaxAdjShifts(const iftAdjRel *A, int *dx, int *dy, int *dz)
{
  int i, d[3];

  *dx = *dy = *dz = 0.0;

  for (i=0; i < A->n; i++) {
    d[0] = abs(A->dx[i]);
    d[1] = abs(A->dy[i]);
    d[2] = abs(A->dz[i]);
    if (*dx < d[0]) *dx = d[0];
    if (*dy < d[1]) *dy = d[1];
    if (*dz < d[2]) *dz = d[2];
  }

}

iftAdjRel *iftAdjacencyBoundaries(const iftAdjRel *A, const iftAdjRel *B_in) 
{
    bool is_3D_adj = false;
    for (int i = 0; i < A->n; i++) {
        if (A->dz[i] != 0) {
            is_3D_adj = true;
            break;
        }
    }

    iftAdjRel *B = NULL;
    if (B_in)
        B = iftCopyAdjacency(B_in);
    else
        B = (is_3D_adj) ? iftSpheric(1.0) : iftCircular(1.0);

    // get the maximum distance between the center of A and the other displacement vectors
    float max_dist = sqrtf(A->dx[A->n-1]*A->dx[A->n-1] + A->dy[A->n-1]*A->dy[A->n-1] + A->dz[A->n-1]*A->dz[A->n-1]);

    int xsize, ysize, zsize;
    xsize = ysize = zsize = (3 * ceil(max_dist)); // find out a given size to create an image that fits A
    if (!is_3D_adj)
        zsize = 1;

    // create "an image" and centralize A on it
    int ***M = iftAlloc(zsize, sizeof(int**));
    for (int z = 0; z < zsize; z++) {
        M[z] = iftAlloc(ysize, sizeof(int*));
        
        for (int y = 0; y < ysize; y++)
            M[z][y] = iftAlloc(xsize, sizeof(int));
    }

    iftVoxel center = {xsize/2, ysize/2, zsize/2};

    // label the adjacent on the "image"
    for (int i = 0; i < A->n; i++) {
        iftVoxel u;
        u.x = center.x + A->dx[i];
        u.y = center.y + A->dy[i];
        u.z = center.z + A->dz[i];

        M[u.z][u.y][u.x] = 1;
    }
    

    iftList *boundaries = iftCreateList();

    for (int i = 0; i < A->n; i++) {
        iftVoxel u;
        u.x = center.x + A->dx[i];
        u.y = center.y + A->dy[i];
        u.z = center.z + A->dz[i];

        for (int j = 1; j < B->n; j++) {
            // adjacent voxel of v
            iftVoxel v;
            v.x = u.x + B->dx[j];
            v.y = u.y + B->dy[j];
            v.z = u.z + B->dz[j];

            // if v is out of domain or it has label != of u
            if ((v.x < 0 || v.x >= xsize) || (v.y < 0 || v.y >= ysize) || (v.z < 0 || v.z >= zsize) ||
                (M[u.z][u.y][u.x] != M[v.z][v.y][v.x])) {
                iftInsertListIntoTail(boundaries, i);
                break;
            }
        }
    }

    iftAdjRel *Abound = iftCreateAdjRel(boundaries->n);

    int i = 0;
    while (!iftIsEmptyList(boundaries)) {
        int j = iftRemoveListTail(boundaries);
        Abound->dx[i] = A->dx[j];
        Abound->dy[i] = A->dy[j];
        Abound->dz[i] = A->dz[j];
        i++;
    }
    iftDestroyList(&boundaries);


    for (int z = 0; z < zsize; z++) {
        for (int y = 0; y < ysize; y++)
            iftFree(M[z][y]);
        iftFree(M[z]);
    }        
    iftFree(M);

    if (B_in)
        iftDestroyAdjRel(&B);

    return Abound;
}

// ---------- iftAdjacency.c end
// ---------- iftCSV.c start 

bool _iftHasCSVHeader(const char *csv_pathname, char separator) 
{
    bool has_header = false;
    
    FILE *fp = fopen(csv_pathname, "rb");
    if (fp == NULL)
        iftError(MSG_FILE_OPEN_ERROR, "_iftHasCSVHeader", csv_pathname);
    
    // reads the first line for checking if it is a CSV header
    char *line = iftGetLine(fp);
    if (line != NULL) {
        bool is_first_line_ok = true;
        
        // if all columns of first line have letters in the beginning of the string, the row can be the header
        char strSeparator[2] = {separator, '\0'};
        iftSList *SL = iftSplitString(line, strSeparator);
        while (!iftIsSListEmpty(SL)) {
            char *column = iftRemoveSListHead(SL);
            
            if (!iftRegexMatch(column, "^[a-zA-Z]+.*$", separator)) {
                is_first_line_ok = false;
                iftFree(column);
                break;
            }
            iftFree(column);
        }
        iftDestroySList(&SL);
        
        if (is_first_line_ok) {
            iftFree(line);
            
            line = iftGetLine(fp);
            if (line != NULL) {
                iftSList *SL = iftSplitString(line, strSeparator);
                iftFree(line);
                
                while (SL->n != 0) {
                    char *column = iftRemoveSListHead(SL);
                    
                    // if at least one column of the second row is a number (integer or real)
                    // the first row is a header
                    if (iftRegexMatch(column, "^[0-9]+(.[0-9]+)?$", separator)) {
                        iftFree(column);
                        has_header = true;
                        break;
                    }
                    iftFree(column);
                }
                iftDestroySList(&SL);
            }
        }
    }
    fclose(fp);
    
    return has_header;
}

void _iftCountNumOfRowsAndColsFromCSVFile(const char *csv_pathname, long *nrows, long *ncols, char separator) 
{
    char strSeparator[2] = {separator, '\0'};
    
    FILE *fp = fopen(csv_pathname, "rb");
    if (fp == NULL)
        iftError(MSG_FILE_OPEN_ERROR, "_iftCountNumOfRowsAndColsFromCSVFile", csv_pathname);
    
    *nrows = 0;
    *ncols = 0;
    
    // reads the first line from the file to get the number of cols
    iftSList *SL = NULL;
    char *line = iftGetLine(fp);
    
    // gets the number of columns from the first line, because such number must be the same for
    // the entire csv
    if (line != NULL) {
        SL = iftSplitString(line, strSeparator);
        (*nrows)++;
        *ncols = SL->n;
        iftDestroySList(&SL);
    }
    
    iftFree(line);
    line = iftGetLine(fp);
    
    // gets each line of the file
    while (line != NULL) {
        SL = iftSplitString(line, strSeparator);
        
        if (*ncols != SL->n)
            iftError("Number of Columns is different in the lines: %d - %d",
                     "_iftCountNumOfRowsAndColsFromCSVFile", *ncols, SL->n);
        
        iftDestroySList(&SL);
        (*nrows)++;
        
        iftFree(line);
        line = iftGetLine(fp);
    }
    
    fclose(fp);
}

iftCSV *_iftCreateCSVWithoutStringAllocation(long nrows, long ncols) 
{
    iftCSV *csv = (iftCSV*) iftAlloc(1, sizeof(iftCSV));
    
    csv->nrows = nrows;
    csv->ncols = ncols;
    
    // allocates the CSV string matrix
    csv->data = (char***) iftAlloc(nrows, sizeof(char**));
    for (long i = 0; i < nrows; i++) {
        csv->data[i] = (char**) iftAlloc(ncols, sizeof(char*));
    }
    
    return csv;
}

iftCSV *iftCreateCSV(long nrows, long ncols) 
{
    iftCSV *csv = (iftCSV*) iftAlloc(1, sizeof(iftCSV));
    csv->header = NULL;
    
    csv->nrows = nrows;
    csv->ncols = ncols;
    
    // allocates the CSV string matrix
    csv->data = (char***) iftAlloc(nrows, sizeof(char**));
    for (long i = 0; i < nrows; i++) {
        csv->data[i] = (char**) iftAlloc(ncols, sizeof(char*));
        
        for (long j = 0; j < ncols; j++) {
            csv->data[i][j] = iftAllocCharArray(IFT_STR_DEFAULT_SIZE);
        }
    }
    
    return csv;
}

void iftDestroyCSV(iftCSV **csv) 
{
    iftCSV *csv_aux = *csv;
    
    if (csv_aux != NULL) {
        if (csv_aux->data != NULL) {
            // deallocates the CSV string matrix
            for (long i = 0; i < csv_aux->nrows; i++) {
                if (csv_aux->data[i] != NULL)
                    for (long j = 0; j < csv_aux->ncols; j++) {
                        iftFree(csv_aux->data[i][j]);
                    }
                iftFree(csv_aux->data[i]);
            }
        }
        iftFree(csv_aux->data);
        
        if (csv_aux->header != NULL) {
            for (int c = 0; c < csv_aux->ncols; c++)
                iftFree(csv_aux->header[c]);
            iftFree(csv_aux->header);
        }
        iftFree(csv_aux);
        *csv = NULL;
    }
}

// ---------- iftCSV.c end
// ---------- iftColor.c start 

iftColorTable *iftCreateColorTable(int n_colors) 
{
  
  if (n_colors <= 0)
      iftError("Invalid num of colors: %d <= 0", "iftCreateColorTable", n_colors);
  
    iftColorTable *ctb = (iftColorTable*) iftAlloc(1, sizeof(iftColorTable));
    ctb->ncolors       = n_colors;
    ctb->color         = (iftColor*) iftAlloc(n_colors, sizeof(iftColor));
    if (ctb->color == NULL)
        iftError("Cannot allocate Color Table", "iftCreateColorTable");

    int S = 255;

    //iftRandomSeed(time(NULL));
    iftRandomSeed(7);

    for (int c = 0, h = iftRandomInteger(0,359); c < n_colors; c++, h = iftRandomInteger(0,359)) {
      ctb->color[c].val[0] = h;
      ctb->color[c].val[1] = S;
      ctb->color[c].val[2] = 255;
      if (c%10 == 0) {
    	S = S - 10;
    	if (S < 20)
    	  S = 255;
      }
      ctb->color[c] = iftRGBtoYCbCr(iftHSVtoRGB(ctb->color[c], 255), 255);
    }

    return ctb;
}

void iftDestroyColorTable(iftColorTable **ctb)
{
    iftColorTable *aux=*ctb;

    if (aux != NULL) {
        iftFree(aux->color);
        iftFree(aux);
        *ctb = NULL;
    }

}

iftColor iftRGBtoYCbCrBT2020(iftColor cin, const int rgbBitDepth, const int yCbCrBitDepth)
{
    int minLum, minChr, quantLum, quantChr;
    iftColor cout;

    switch (yCbCrBitDepth) {
        case 8:
            minLum = 16; // 16 * 2^(bitDepth-8)
            minChr = 128; // 128 * 2^(bitDepth-8)
            quantLum = 219.0; // 219 * 2^(bitDepth-8)
            quantChr = 224.0;  // 224 * 2^(bitDepth-8)
            break;
        case 10:
            minLum = 64; // 16 * 2^(bitDepth-8)
            minChr = 512; // 128 * 2^(bitDepth-8)
            quantLum = 876; // 219 * 2^(bitDepth-8)
            quantChr = 896;  // 224 * 2^(bitDepth-8)
            break;
        case 12:
            minLum = 256; // 16 * 2^(bitDepth-8)
            minChr = 2048; // 128 * 2^(bitDepth-8)
            quantLum = 3504; // 219 * 2^(bitDepth-8)
            quantChr = 3584;  // 224 * 2^(bitDepth-8)
            break;
        case 16:
            minLum = 4096; // 16 * 2^(bitDepth-8)
            minChr = 32768; // 128 * 2^(bitDepth-8)
            quantLum = 56064.0; // 219 * 2^(bitDepth-8)
            quantChr = 57344.0;  // 224 * 2^(bitDepth-8)
            break;
        default:
            iftError("Bit depth not specified in BT.2020", "iftRGBtoYCbCrBT2020");
            cout.val[0] = cout.val[1] = cout.val[2] = 0;
            return cout;
    }

    double maxRgbValue = (double) ((1 << rgbBitDepth) - 1);
    double r = cin.val[0] / maxRgbValue;
    double g = cin.val[1] / maxRgbValue;
    double b = cin.val[2] / maxRgbValue;

    double y = 0.2627 * r + 0.6780 * g + 0.0593 * b;
    double cb = (b - y) / 1.8814;
    double cr = (r - y) / 1.4746;

    // clip luminance to [0..1] and chrominance to [-0.5..0.5]
    if (y < 0.0) y = 0.0;
    else if (y > 1.0) y = 1.0;
    if (cb < -0.5) cb = -0.5;
    else if (cb > 0.5) cb = 0.5;
    if (cr < -0.5) cr = -0.5;
    else if (cr > 0.5) cr = 0.5;

    // perform quantization
    cout.val[0] = (int) (y * quantLum) + minLum;
    cout.val[1] = (int) (cb * quantChr) + minChr;
    cout.val[2] = (int) (cr * quantChr) + minChr;

    return cout;
}

iftColor iftRGBtoYCbCr(iftColor cin, int normalization_value)
{
    iftColor cout;
    float a = (16.0/255.0)*(float)normalization_value;
    float b = (128.0/255.0)*(float)normalization_value;

    cout.val[0]=(int)(0.256789062*(float)cin.val[0]+
                      0.504128906*(float)cin.val[1]+
                      0.09790625*(float)cin.val[2]+a);
    cout.val[1]=(int)(-0.148222656*(float)cin.val[0]+
                      -0.290992187*(float)cin.val[1]+
                      0.439214844*(float)cin.val[2]+b);
    cout.val[2]=(int)(0.439214844*(float)cin.val[0]+
                      -0.367789063*(float)cin.val[1]+
                      -0.071425781*(float)cin.val[2]+b);

    for(int i=0; i < 3; i++) {
        if (cout.val[i] < 0) cout.val[i] = 0;
        if (cout.val[i] > normalization_value) cout.val[i] = normalization_value;
    }

    return(cout);
}

iftColor iftYCbCrtoRGB(iftColor cin, int normalization_value)
{
    iftColor cout;
    float a = (16.0/255.0)*(float)normalization_value;
    float b = (128.0/255.0)*(float)normalization_value;

    cout.val[0]=(int)(1.164383562*((float)cin.val[0]-a)+
                      1.596026786*((float)cin.val[2]-b));

    cout.val[1]=(int)(1.164383562*((float)cin.val[0]-a)+
                      -0.39176229*((float)cin.val[1]-b)+
                      -0.812967647*((float)cin.val[2]-b));

    cout.val[2]=(int)(1.164383562*((float)cin.val[0]-a)+
                      2.017232143*((float)cin.val[1]-b));

    for(int i=0; i < 3; i++) {
        if (cout.val[i] < 0) cout.val[i] = 0;
        if (cout.val[i] > normalization_value) cout.val[i] = normalization_value;
    }

    return(cout);
}

iftColor iftYCbCrBT2020toRGB(iftColor cin, const int yCbCrBitDepth, const int rgbBitDepth)
{
    int minLum, minChr;
    double quantLum, quantChr;
    iftColor cout;

    switch (yCbCrBitDepth) {
        case 8:
            minLum = 16; // 16 * 2^(bitDepth-8)
            minChr = 128; // 128 * 2^(bitDepth-8)
            quantLum = 219.0; // 219 * 2^(bitDepth-8)
            quantChr = 224.0;  // 224 * 2^(bitDepth-8)
            break;
        case 10:
            minLum = 64; // 16 * 2^(bitDepth-8)
            minChr = 512; // 128 * 2^(bitDepth-8)
            quantLum = 876.0; // 219 * 2^(bitDepth-8)
            quantChr = 896.0;  // 224 * 2^(bitDepth-8)
            break;
        case 12:
            minLum = 256; // 16 * 2^(bitDepth-8)
            minChr = 2048; // 128 * 2^(bitDepth-8)
            quantLum = 3504.0; // 219 * 2^(bitDepth-8)
            quantChr = 3584.0;  // 224 * 2^(bitDepth-8)
            break;
        case 16:
            minLum = 4096; // 16 * 2^(bitDepth-8)
            minChr = 32768; // 128 * 2^(bitDepth-8)
            quantLum = 56064.0; // 219 * 2^(bitDepth-8)
            quantChr = 57344.0;  // 224 * 2^(bitDepth-8)
            break;
        default:
            iftError("Bit depth not specified in BT.2020", "iftYCbCrBT2020toRGB");
            cout.val[0] = cout.val[1] = cout.val[2] = 0;
            return cout;
    }

    double y = (cin.val[0] - minLum) / quantLum;
    double cb = (cin.val[1] - minChr) / quantChr;
    double cr = (cin.val[2] - minChr) / quantChr;

    double r = cr * 1.4746 + y;
    double b = cb * 1.8814 + y;
    double g = (y - 0.2627 * r - 0.0593 * b) / 0.6780;

    // clip rgb values to [0..1]
    if (r < 0.0) r = 0.0;
    else if (r > 1.0) r = 1.0;
    if (g < 0.0) g = 0.0;
    else if (g > 1.0) g = 1.0;
    if (b < 0.0) b = 0.0;
    else if (b > 1.0) b = 1.0;

    // perform quantization
    double maxRgbValue = (double) ((1 << rgbBitDepth) - 1);
    cout.val[0] = (int) (r * maxRgbValue);
    cout.val[1] = (int) (g * maxRgbValue);
    cout.val[2] = (int) (b * maxRgbValue);

    return cout;
}

iftFColor iftRGBtoLabNorm(iftColor rgb, int normalization_value)
{
    //RGB to XYZ

    float R = rgb.val[0]/(float)normalization_value;
    float G = rgb.val[1]/(float)normalization_value;
    float B = rgb.val[2]/(float)normalization_value;

    if(R <= 0.04045)	R = R/12.92;
    else	        R = pow((R+0.055)/1.055,2.4);

    if(G <= 0.04045)	G = G/12.92;
    else		G = pow((G+0.055)/1.055,2.4);

    if(B <= 0.04045)	B = B/12.92;
    else		B = pow((B+0.055)/1.055,2.4);

    float X = (0.4123955889674142161*R + 0.3575834307637148171*G + 0.1804926473817015735*B);
    float Y = (0.2125862307855955516*R + 0.7151703037034108499*G + 0.07220049864333622685*B);
    float Z = (0.01929721549174694484*R + 0.1191838645808485318*G + 0.9504971251315797660*B);

    //XYZ to lab
    X /= WHITEPOINT_X;
    Y /= WHITEPOINT_Y;
    Z /= WHITEPOINT_Z;
    X = LABF(X);
    Y = LABF(Y);
    Z = LABF(Z);
    float L = 116*Y - 16;
    float a = 500*(X - Y);
    float b = 200*(Y - Z);

    iftFColor lab;
    lab.val[0] = L;
    lab.val[1] = a;
    lab.val[2] = b;

    return lab;
}

iftColor iftRGBtoHSV(iftColor cin, int normalization_value) 
{
    float r = ((float)cin.val[0]/normalization_value),
            g = ((float)cin.val[1]/normalization_value),
            b = ((float)cin.val[2]/normalization_value), v, x, f;
    float a[3];
    int   i;
    iftColor cout;

    // RGB are each on [0, 1]. S and V are returned on [0, 1] and H is
    // returned on [0, 6].

    x = iftMin(iftMin(r, g), b);
    v = iftMax(iftMax(r, g), b);
    if (v == x) {
        a[0]=0.0;
        a[1]=0.0;
        a[2]=v;
    } else {
        f = (r == x) ? g - b : ((g == x) ? b - r : r - g);
        i = (r == x) ? 3 : ((g == x) ? 5 : 1);
        a[0]=((float)i)-f/(v-x);
        a[1]=(v-x)/v;
        a[2]=0.299*r+0.587*g+0.114*b;
    }

    // (un)normalize

    cout.val[0] = (int)(a[0]*60.0);
    cout.val[1] = (int)(a[1]*normalization_value);
    cout.val[2] = (int)(a[2]*normalization_value);

    return(cout);
}

iftColor iftHSVtoRGB(iftColor cin, int normalization_value) 
{
    // H is given on [0, 6]. S and V are given on [0, 1].
    // RGB are each returned on [0, 1].
    float h = ((float)cin.val[0]/60.0),
            s = ((float)cin.val[1]/normalization_value),
            v = ((float)cin.val[2]/normalization_value), m, n, f;
    float a[3]={0,0,0};
    int i;
    iftColor cout;

    if (s==0.0) {
        a[0]=a[1]=a[2]=v;
    } else {
        i = (int) floor(h);
        f = h - (float)i;
        if(!(i & 1)) f = 1 - f; // if i is even
        m = v * (1 - s);
        n = v * (1 - s * f);
        switch (i) {
            case 6:
            case 0: a[0]=v; a[1]=n; a[2]=m; break;
            case 1: a[0]=n; a[1]=v; a[2]=m; break;
            case 2: a[0]=m; a[1]=v; a[2]=n; break;
            case 3: a[0]=m; a[1]=n; a[2]=v; break;
            case 4: a[0]=n; a[1]=m; a[2]=v; break;
            case 5: a[0]=v; a[1]=m; a[2]=n; break;
        }
    }

    // (un)normalize
    for(i=0;i<3;i++)
        cout.val[i]=a[i]*normalization_value;

    return(cout);
}

// ---------- iftColor.c end
// ---------- iftDHeap.c start

iftDHeap *iftCreateDHeap(int n, double *value) 
{
    iftDHeap *H = NULL;
    int i;
    
    if (value == NULL) {
        iftError("Cannot create heap without priority value map", "iftCreateDHeap");
    }
    
    H = (iftDHeap *) iftAlloc(1, sizeof(iftDHeap));
    if (H != NULL) {
        H->n       = n;
        H->value   = value;
        H->color   = (char *) iftAlloc(sizeof(char), n);
        H->node    = (int *) iftAlloc(sizeof(int), n);
        H->pos     = (int *) iftAlloc(sizeof(int), n);
        H->last    = -1;
        H->removal_policy = MINVALUE;
        if (H->color == NULL || H->pos == NULL || H->node == NULL)
            iftError(MSG_MEMORY_ALLOC_ERROR, "iftCreateDHeap");
        for (i = 0; i < H->n; i++) {
            H->color[i] = IFT_WHITE;
            H->pos[i]   = -1;
            H->node[i] = -1;
        }
    }
    else
        iftError(MSG_MEMORY_ALLOC_ERROR, "iftCreateDHeap");
    
    return H;
}

void iftDestroyDHeap(iftDHeap **H) 
{
    iftDHeap *aux = *H;
    if (aux != NULL) {
        if (aux->node != NULL) iftFree(aux->node);
        if (aux->color != NULL) iftFree(aux->color);
        if (aux->pos != NULL)   iftFree(aux->pos);
        iftFree(aux);
        *H = NULL;
    }
}

char iftFullDHeap(iftDHeap *H) 
{
    if (H->last == (H->n - 1))
        return 1;
    else
        return 0;
}

char iftEmptyDHeap(iftDHeap *H) 
{
    if (H->last == -1){
        return 1;
    }else{
        return 0;
    }
}

char iftInsertDHeap(iftDHeap *H, int node) 
{
    
    if (!iftFullDHeap(H)) {
        H->last++;
        H->node[H->last] = node;
        H->color[node]   = IFT_GRAY;
        H->pos[node]     = H->last;
        iftGoUpDHeap(H, H->last);
        return 1;
    } else {
        iftWarning("DHeap is full","iftInsertDHeap");
        return 0;
    }
    
}

int iftRemoveDHeap(iftDHeap *H) 
{
    int node= IFT_NIL;
    
    if (!iftEmptyDHeap(H)) {
        node = H->node[0];
        H->pos[node]   = -1;
        H->color[node] = IFT_BLACK;
        H->node[0]     = H->node[H->last];
        H->pos[H->node[0]] = 0;
        H->node[H->last] = -1;
        H->last--;
        iftGoDownDHeap(H, 0);
    }else{
        iftWarning("DHeap is empty","iftRemoveDHeap");
    }
    
    return node;
    
}

void    iftRemoveDHeapElem(iftDHeap *H, int pixel)
{
    
    if(H->pos[pixel] == -1)
        iftError("Element is not in the Heap", "iftRemoveDHeapElem");
    
    double aux = H->value[pixel];
    
    if(H->removal_policy == MINVALUE)
        H->value[pixel] = IFT_INFINITY_DBL_NEG;
    else
        H->value[pixel] = IFT_INFINITY_DBL;
    
    iftGoUpDHeap(H, H->pos[pixel]);
    iftRemoveDHeap(H);
    
    H->value[pixel] = aux;
    H->color[pixel] = IFT_WHITE;
    
}

void  iftGoUpDHeap(iftDHeap *H, int i) 
{
    int j = iftDad(i);
    
    if(H->removal_policy == MINVALUE){
        
        while ((j >= 0) && (H->value[H->node[j]] > H->value[H->node[i]])) {
            iftSwap(H->node[j], H->node[i]);
            H->pos[H->node[i]] = i;
            H->pos[H->node[j]] = j;
            i = j;
            j = iftDad(i);
        }
    }
    else{ /* removal_policy == MAXVALUE */
        
        while ((j >= 0) && (H->value[H->node[j]] < H->value[H->node[i]])) {
            iftSwap(H->node[j], H->node[i]);
            H->pos[H->node[i]] = i;
            H->pos[H->node[j]] = j;
            i = j;
            j = iftDad(i);
        }
    }
}

void iftGoDownDHeap(iftDHeap *H, int i) 
{
    int j, left = iftLeftSon(i), right = iftRightSon(i);
    
    j = i;
    if(H->removal_policy == MINVALUE){
        
        if ((left <= H->last) &&
            (H->value[H->node[left]] < H->value[H->node[i]]))
            j = left;
        if ((right <= H->last) &&
            (H->value[H->node[right]] < H->value[H->node[j]]))
            j = right;
    }
    else{ /* removal_policy == MAXVALUE */
        
        if ((left <= H->last) &&
            (H->value[H->node[left]] > H->value[H->node[i]]))
            j = left;
        if ((right <= H->last) &&
            (H->value[H->node[right]] > H->value[H->node[j]]))
            j = right;
    }
    
    if(j != i) {
        iftSwap(H->node[j], H->node[i]);
        H->pos[H->node[i]] = i;
        H->pos[H->node[j]] = j;
        iftGoDownDHeap(H, j);
    }
}

void iftResetDHeap(iftDHeap *H)
{
    int i;
    
    for (i=0; i < H->n; i++) {
        H->color[i] = IFT_WHITE;
        H->pos[i]   = -1;
        H->node[i] = -1;
    }
    H->last = -1;
}

// ---------- iftDHeap.c end
// ---------- iftFIFO.c start 

iftFIFO *iftCreateFIFO(int n)
{
    iftFIFO *F=(iftFIFO *)iftAlloc(1,sizeof(iftFIFO));
    
    F->FIFO  = iftAllocIntArray(n);
    F->color = iftAllocCharArray(n);
    F->n     = n;
    F->first=F->last=0;
    
    return(F);
}

void iftDestroyFIFO(iftFIFO **F) 
{
    if (F != NULL) {
        iftFIFO *aux=*F;
        
        if (aux != NULL) {
            iftFree(aux->FIFO);
            iftFree(aux->color);
            iftFree(aux);
            *F = NULL;
        }
    }
}

char iftInsertFIFO(iftFIFO *F, int elem)
{
    if (iftFullFIFO(F)){
        iftWarning("FIFO is full","iftInsertFIFO");
        return 0;
    }
    F->color[elem]= IFT_GRAY;
    F->FIFO[F->last]=elem;  F->last++;
    
    return 1;
}

int      iftRemoveFIFO(iftFIFO *F)
{
    int node= IFT_NIL;
    
    if (!iftEmptyFIFO(F)){
        node = F->FIFO[F->first];  F->first++;
        F->color[node]= IFT_BLACK;
    }else{
        iftWarning("FIFO is empty","iftRemoveFIFO");
    }
    
    return node;
}

bool iftFullFIFO(iftFIFO *F)
{
    if (F->last==F->n) {
        return(1);
    }else
        return(0);
}

bool iftEmptyFIFO(iftFIFO *F)
{
    if (F->first == F->last) {
        // Changed by Falcao and Nikolas
        // iftResetFIFO(F);
        return(1);
    }else
        return(0);
}

void     iftResetFIFO(iftFIFO *F)
{
    int p;
    for (p=0; p < F->n; p++)
        F->color[p] = IFT_WHITE;
    F->first=F->last=0;
}

int    iftColorFIFO(iftFIFO *F, int pos)
{
    return F->color[pos];
}

// ---------- iftFIFO.c end
// ---------- iftFile.c start

bool iftFileExists(const char *pathname) 
{
    return (iftPathnameExists(pathname) && !iftDirExists(pathname));
}

const char *iftFileExt(const char *pathname) 
{
    if (pathname == NULL)
        iftError("Pathname is NULL", "iftFileExt");
    
    const char *dot = strrchr(pathname, '.'); // returns a pointer to the last occurrence of '.'
    
    if ( (!dot) || (dot == pathname)) {
        return ("");
    } else {
        if (iftRegexMatch(pathname, "^.*\\.tar\\.(gz|bz|bz2)$") || iftRegexMatch(pathname, "^.*\\.(scn|nii)\\.gz$")) {
            dot -= 4; // points to the penultimate dot '.'
        }
        
        return dot; // returns the extension with '.'
    }
}

char *iftJoinPathnames(long n, ...) 
{
    if (n <= 0)
        iftError("Number of pathnames to be concatenated is <= 0", "iftJoinPathnames");
    
    long out_str_size = 1; // '\0'
    
    // Counts the size of the concatenated string
    va_list path_list;
    va_start(path_list, n);
    for (int i = 0; i < n; i++)
        out_str_size += strlen(va_arg(path_list, char*)) + 1; // one char for '/' (separation char)
    va_end(path_list);
    
    char *joined_path = iftAllocCharArray(out_str_size);
    char *aux = iftAllocCharArray(out_str_size);
    
    va_start(path_list, n);
    strcpy(joined_path, va_arg(path_list, char*));
    
    for (int i = 1; i < n; i++) {
        char *path = va_arg(path_list, char*);
        if (iftStartsWith(path, IFT_SEP_C))
            path++; // skip the first char, which is the directory separator
        
        if (iftEndsWith(joined_path, IFT_SEP_C))
            sprintf(aux, "%s%s", joined_path, path);
        else
            sprintf(aux, "%s%s%s", joined_path, IFT_SEP_C, path);

        iftFree(joined_path);
        joined_path = iftCopyString(aux);
    }
    iftFree(aux);
    
    return joined_path;
}

char *iftFilename(const char *pathname, const char *suffix) 
{
    if (pathname == NULL)
        iftError("Pathname is NULL", "iftFilename");
    
    char *base = iftSplitStringAt(pathname, IFT_SEP_C, -1);
    
    if ((suffix != NULL) && (!iftCompareStrings(suffix, ""))) {
        char *out_base = iftRemoveSuffix(base, suffix);
        iftFree(base);
        base = out_base;
    }
    
    return base;
}

void iftDestroyFile(iftFile **f) {
    if (f != NULL) {
        iftFile *f_aux = *f;
        
        if (f_aux != NULL) {
            if (f_aux->path != NULL) {
                iftFree(f_aux->path);
                f_aux->path = NULL;
            }
            if(f_aux->suffix != NULL) {
                iftFree(f_aux->suffix);
                f_aux->suffix = NULL;
            }
            iftFree(f_aux);
            *f = NULL;
        }
    }
}

// ---------- iftFile.c end
// ---------- iftBoundingBox.h start 
// ---------- iftBoundingBox.h end
// ---------- iftImage.c start
#include "iftPng.h"
#include "jpeglib.h"

iftImage *iftReadImageByExt(const char *format, ...) 
{
    va_list args;
    char filename[IFT_STR_DEFAULT_SIZE];

    va_start(args, format);
    vsprintf(filename, format, args);
    va_end(args);

    if (!iftFileExists(filename))
        iftError("Image %s does not exist", "iftReadImageByExt", filename);

    iftImage *img = NULL;
    char *ext = iftLowerString(iftFileExt(filename));

    if(iftCompareStrings(ext, ".png")) {
        img = iftReadImagePNG(filename);
    }
    else if (iftCompareStrings(ext, ".pgm")){
        FILE *fp = fopen(filename,"r");
        char type[10];
        if(fscanf(fp,"%s",type)!=1) iftError("Reading Error", "iftReadImageByExt");
        if (iftCompareStrings(type,"P5")){
            fclose(fp);
            img   = iftReadImageP5(filename);
        } else {
            fclose(fp);
            img   = iftReadImageP2(filename);
        }
    } else if (iftCompareStrings(ext, ".ppm")){
        img   = iftReadImageP6(filename);
    } else if (iftCompareStrings(ext, ".scn")){
        img   = iftReadImage(filename);
    } else if (iftCompareStrings(ext, ".zscn") || iftCompareStrings(ext, ".scn.gz")) {
        img   = iftReadImageGZip(filename);
    } else if (iftCompareStrings(ext, ".jpg") || iftCompareStrings(ext, ".jpeg")){
        img = iftReadImageJPEG(filename);
    } else {
        iftError("Invalid image format: \"%s\" - Try .scn, .zscn, .scn.gz, .ppm, .pgm, .jpg, .png",
                 "iftReadImageByExt", ext);
    }

    iftFree(ext);
    return(img);
}

iftImage  *iftCreateImage(int xsize,int ysize,int zsize) 
{
    int *val = iftAllocIntArray(xsize*ysize*zsize);

    return iftCreateImageFromBuffer(xsize, ysize, zsize, val);
}

iftImage *iftCopyImage(const iftImage *img) 
{
    if (img == NULL)
        return NULL;

    iftImage *imgc=iftCreateImage(img->xsize,img->ysize,img->zsize);
    iftCopyImageInplace(img, imgc);

    return(imgc);
}

void iftDestroyImage(iftImage **img) 
{
    if(img != NULL) {
        iftImage *aux = *img;

        if (aux != NULL) {
            if (aux->val != NULL) iftFree(aux->val);
            if (aux->Cb != NULL) iftFree(aux->Cb);
            if (aux->Cr != NULL) iftFree(aux->Cr);
            if (aux->alpha != NULL) iftFree(aux->alpha);
            if (aux->tby != NULL) iftFree(aux->tby);
            if (aux->tbz != NULL) iftFree(aux->tbz);
            iftFree(aux);
            *img = NULL;
        }
    }
}

void iftWriteImageByExt(const iftImage *img, const char *format, ...) 
{
    if (img == NULL)
        iftWarning("Image is NULL... Nothing to write", "iftWriteImageByExt");
    else {
        char command[400];

        va_list args;
        char filename[IFT_STR_DEFAULT_SIZE];

        va_start(args, format);
        vsprintf(filename, format, args);
        va_end(args);

        char *parent_dir = iftParentDir(filename);
        if (!iftDirExists(parent_dir))
            iftMakeDir(parent_dir);
        iftFree(parent_dir);

        char *ext = iftLowerString(iftFileExt(filename));

        if(iftCompareStrings(ext, ".png")) {
            iftWriteImagePNG(img,filename);
        } else if (iftCompareStrings(ext, ".scn")) {
            iftWriteImage(img, filename);
        } else if (iftCompareStrings(ext, ".scn.gz") || iftCompareStrings(ext, ".zscn")) {
            iftWriteImageGZip(img, filename);
        }else if (iftCompareStrings(ext, ".pgm")) {
            if (iftMaximumValue(img)>255)
                iftWriteImageP2(img,filename);
            else
                iftWriteImageP5(img,filename);
        } else if (iftCompareStrings(ext, ".ppm")){
            iftWriteImageP6(img,filename);
        } else if (iftIsColorImage(img)){
            iftWriteImageP6(img,"temp.ppm");
            sprintf(command,"convert temp.ppm %s",filename);
            if (system(command)==-1)
                iftError("Program convert failed or is not installed", "iftWriteImageByExt");
            if (system("rm -f temp.ppm")==-1)
                iftError("Cannot remore temp.ppm", "iftWriteImageByExt");
        } else if(iftCompareStrings(ext, ".jpg") || iftCompareStrings(ext, ".jpeg")) {
            iftWriteImageJPEG(img,filename);
        } else {
            printf("Invalid image format: %s. Please select among the accepted ones: .scn, .zscn, .scn.gz, .ppm, .pgm, .png, .img\n",ext);
            exit(-1);
        }


        iftFree(ext);
    }
}

int iftMaximumValue(const iftImage *img) 
{
    if (img == NULL)
        iftError("Image is NULL", "iftMaximumValue");

    iftBoundingBox bb;
    bb.begin.x = bb.begin.y = bb.begin.z = 0;
    bb.end.x   = img->xsize-1;
    bb.end.y   = img->ysize-1;
    bb.end.z   = img->zsize-1;

    return iftMaximumValueInRegion(img, bb);
}

int iftMinimumValue(const iftImage *img) 
{
    int img_min_val = IFT_INFINITY_INT;

    for (int p = 0; p < img->n; p++)
        if (img_min_val > img->val[p])
            img_min_val = img->val[p];

    return img_min_val;
}

inline iftVoxel iftGetVoxelCoord(const iftImage *img, int p)
{
    /* old
     * u.x = (((p) % (((img)->xsize)*((img)->ysize))) % (img)->xsize)
     * u.y = (((p) % (((img)->xsize)*((img)->ysize))) / (img)->xsize)
     * u.z = ((p) / (((img)->xsize)*((img)->ysize)))
     */
    iftVoxel u;
    div_t res1 = div(p, img->xsize * img->ysize);
    div_t res2 = div(res1.rem, img->xsize);

    u.x = res2.rem;
    u.y = res2.quot;
    u.z = res1.quot;

    return u;
}

iftImage *iftSelectImageDomain(int xsize, int ysize, int zsize)
{
    iftImage *mask=iftCreateImage(xsize,ysize,zsize);

    iftSetImage(mask,1);
    return(mask);
}

iftBoundingBox iftMinBoundingBox(const iftImage *img, iftVoxel *gc_out) 
{
    if (img == NULL)
        iftError("Image is NULL", "iftMinBoundingBox");

    long n = 0; // number of spels non-background (non-zero)
    iftVoxel gc = {0.0, 0.0, 0.0};
    iftBoundingBox mbb;
    mbb.begin.x = mbb.begin.y = mbb.begin.z = IFT_INFINITY_INT;
    mbb.end.x = mbb.end.y = mbb.end.z = IFT_INFINITY_INT_NEG;

    for (long p = 0; p < img->n; p++) {
        if (img->val[p] != 0) {
            iftVoxel v = iftGetVoxelCoord(img, p);

            mbb.begin.x = iftMin(mbb.begin.x, v.x);
            mbb.begin.y = iftMin(mbb.begin.y, v.y);
            mbb.begin.z = iftMin(mbb.begin.z, v.z);

            mbb.end.x = iftMax(mbb.end.x, v.x);
            mbb.end.y = iftMax(mbb.end.y, v.y);
            mbb.end.z = iftMax(mbb.end.z, v.z);

            gc.x += v.x;
            gc.y += v.y;
            gc.z += v.z;
            n++;
        }
    }

    if (mbb.begin.x == IFT_INFINITY_INT) {
        mbb.begin.x = mbb.begin.y = mbb.begin.z = -1;
        mbb.end.x   = mbb.end.y   = mbb.end.z   = -1;
        gc.x        = gc.y        = gc.z        = -1.0;
    } else {
        gc.x /= n;
        gc.y /= n;
        gc.z /= n;
    }

    if (gc_out != NULL)
        *gc_out = gc;

    return mbb;
}

iftImage *iftReadImage(const char *format, ...) 
{
    iftImage *img    = NULL;
    FILE     *fp     = NULL;
    uchar    *data8  = NULL;
    ushort   *data16 = NULL;
    int      *data32 = NULL;
    char     type[10];
    int      p, v, xsize, ysize, zsize;

    va_list args;
    char    filename[IFT_STR_DEFAULT_SIZE];

    va_start(args, format);
    vsprintf(filename, format, args);
    va_end(args);

    fp = fopen(filename, "rb");

    if (fp == NULL) {
        iftError("Cannot open file: \"%s\"", "iftReadImage", filename);
    }
    if (fscanf(fp, "%s\n", type) != 1)
        iftError("Reading error: Image type", "iftReadImage");

    if (iftCompareStrings(type, "SCN")) {

        //iftSkipComments(fp);

        if (fscanf(fp, "%d %d %d\n", &xsize, &ysize, &zsize) != 3)
            iftError("Reading error: Image resolution/size", "iftReadImage");
        img = iftCreateImage(xsize, ysize, zsize);
        if (fscanf(fp, "%f %f %f\n", &img->dx, &img->dy, &img->dz) != 3) {
            iftError("Reading error: Pixel/Voxel size", "iftReadImage");
        }
        if (fscanf(fp, "%d", &v) != 1)
            iftError("Reading error", "iftReadImage");

        while (fgetc(fp) != '\n');

        if (v == 8) {
            data8 = iftAllocUCharArray(img->n);
            if (fread(data8, sizeof(uchar), img->n, fp) != img->n)
                iftError("Reading error", "iftReadImage");
            for (p = 0; p < img->n; p++)
                img->val[p] = (int) data8[p];
            iftFree(data8);
        } else if (v == 16) {
            data16 = iftAllocUShortArray(img->n);
            if (fread(data16, sizeof(ushort), img->n, fp) != img->n)
                iftError("Reading error 16 bits", "iftReadImage");
            for (p = 0; p < img->n; p++)
                img->val[p] = (int) data16[p];
            iftFree(data16);
        } else if (v == 32) {
            data32 = iftAllocIntArray(img->n);
            if (fread(data32, sizeof(int), img->n, fp) != img->n)
                iftError("Reading error", "iftReadImage");
            for (p = 0; p < img->n; p++)
                img->val[p] = data32[p];
            iftFree(data32);
        } else {
            iftError("Input scene must be 8, 16, or 32 bit", "iftReadImage");
        }
    } else {
        iftError("Invalid file type", "iftReadImage");
    }

    fclose(fp);
    return (img);
}

png_bytep* iftReadPngImageAux(const char *file_name, png_structp *png_ptr, png_infop *info_ptr)
{
    png_byte header[8];    // 8 is the maximum size that can be checked

    /* open file and test for it being a png */
    FILE *fp = fopen(file_name, "rb");
    if (!fp)
        iftError("File %s could not be opened for reading", "iftReadPngImageAux", file_name);
    if (fread(header, 1, 8, fp)!=8) iftError("Reading error", "iftReadPngImageAux");
    if (png_sig_cmp(header, 0, 8))
        iftError("File %s is not recognized as a PNG file", "iftReadPngImageAux", file_name);

    int height;
    png_bytep * row_pointers;

    /* initialize stuff */
    png_structp ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);

    if (!ptr)
        iftError("Internal error: png_create_read_struct failed", "iftReadImagePNG");

    *png_ptr = ptr;
    *info_ptr = png_create_info_struct(*png_ptr);
    int depth = png_get_bit_depth((*png_ptr), (*info_ptr));
    if(depth < 8){
        png_set_expand_gray_1_2_4_to_8(ptr);
    }


    if (!(*info_ptr))
        iftError("Internal error: png_create_info_struct failed", "iftReadImagePNG");

    if (setjmp(png_jmpbuf(*png_ptr)))
        iftError("Internal error: Error during init_io", "iftReadImagePNG");


    png_init_io(*png_ptr, fp);
    png_set_sig_bytes(*png_ptr, 8);

    png_read_info(*png_ptr, *info_ptr);
	// reduces the pixels back down to the original bit depth
	//png_color_8p sig_bit = NULL;
	//if (png_get_sBIT(*png_ptr, *info_ptr, &sig_bit)) {
	//	png_set_shift(*png_ptr, sig_bit);
	//}
    height = png_get_image_height(*png_ptr, *info_ptr);
    png_read_update_info(*png_ptr, *info_ptr);


    /* read file */
    if (setjmp(png_jmpbuf(*png_ptr)))
        iftError("Internal error: Error during read_image", "iftReadImagePNG");

    row_pointers = (png_bytep*) iftAlloc(height, sizeof(png_bytep));
    for (int y=0; y<height; y++)
        row_pointers[y] = (png_byte*) iftAlloc(png_get_rowbytes(*png_ptr, *info_ptr), 1);


    png_read_image(*png_ptr, row_pointers);

    fclose(fp);

    return row_pointers;
}

iftImage* iftReadImagePNG(const char* format, ...) 
{

    va_list args;
    char filename[IFT_STR_DEFAULT_SIZE];

    va_start(args, format);
    vsprintf(filename, format, args);
    va_end(args);

    png_infop info_ptr;
    png_structp png_ptr;
    png_bytep *row_pointers;

    row_pointers = iftReadPngImageAux(filename, &png_ptr, &info_ptr);

    int width, height, color_type, depth;

    width = png_get_image_width(png_ptr, info_ptr);
    height = png_get_image_height(png_ptr, info_ptr);
    color_type = png_get_color_type(png_ptr, info_ptr);
    depth = png_get_bit_depth(png_ptr, info_ptr);
    iftImage* img = iftCreateImage(width, height, 1);
    unsigned int numberChannels = png_get_channels(png_ptr, info_ptr);

    int byteshift = depth/8;

    int x, y;

    int p = 0;

    if(color_type==PNG_COLOR_TYPE_GRAY)//gray image
    {
        for (y=0; y<height; y++) {
            png_byte* row = row_pointers[y];
            for (x=0; x<width; x++) {
                png_byte* ptr = &(row[x*numberChannels*byteshift]);
                img->val[p] = ptr[0];
                if(depth==16) {
                    img->val[p] = (img->val[p]<<8)+ptr[1];
                }
                p++;
            }
        }
    }else if(color_type==PNG_COLOR_TYPE_GRAY_ALPHA ){
        if(img->alpha == NULL){
            iftSetAlpha(img,0);
        }
        for (y=0; y<height; y++) {
            png_byte* row = row_pointers[y];
            for (x=0; x<width; x++) {
                png_byte* ptr = &(row[x*numberChannels*byteshift]);
                if(depth == 8){
                    img->val[p] = ptr[0];
                    img->alpha[p] = ptr[1];
                }
                else if(depth==16) {
                    img->val[p] = ptr[0];
                    img->val[p] = (img->val[p]<<8)+ptr[1];
                    img->alpha[p] = ptr[2];
                    img->alpha[p] = (img->alpha[p]<<8)+ptr[3];
                }
                p++;
            }
        }
    }
    else if(color_type == PNG_COLOR_TYPE_RGB){//color image

        iftSetCbCr(img, 128);
        iftColor rgb, ycbcr;

        for (y=0; y<height; y++) {
            png_byte* row = row_pointers[y];

            for (x=0; x<width; x++) {
                png_byte* ptr = &(row[x*numberChannels*byteshift]);
                rgb.val[0] = ptr[0*byteshift];
                rgb.val[1] = ptr[1*byteshift];
                rgb.val[2] = ptr[2*byteshift];

                if(depth==16) { //read second byte in case of 16bit images
                    rgb.val[0] = (rgb.val[0]<<8) + ptr[1];
                    rgb.val[1] = (rgb.val[1]<<8) + ptr[3];
                    rgb.val[2] = (rgb.val[2]<<8) + ptr[5];
                }

                ycbcr = iftRGBtoYCbCrBT2020(rgb, depth, depth);

                img->val[p] = ycbcr.val[0];
                img->Cb[p]  = ycbcr.val[1];
                img->Cr[p]  = ycbcr.val[2];

                p++;
            }
        }
    }else if(color_type == PNG_COLOR_TYPE_RGB_ALPHA){
        iftSetCbCr(img, 128);
        iftColor rgb, ycbcr;
        if(img->alpha == NULL){
            iftSetAlpha(img,0);
        }

        for (y=0; y<height; y++) {
            png_byte* row = row_pointers[y];
            for (x=0; x<width; x++) {
                png_byte* ptr = &(row[x*numberChannels*byteshift]);
                rgb.val[0] = ptr[0*byteshift];
                rgb.val[1] = ptr[1*byteshift];
                rgb.val[2] = ptr[2*byteshift];
                ushort alpha = ptr[3*byteshift];

                if(depth==16) { //read second byte in case of 16bit images
                    rgb.val[0] = (rgb.val[0]<<8) + ptr[1];
                    rgb.val[1] = (rgb.val[1]<<8) + ptr[3];
                    rgb.val[2] = (rgb.val[2]<<8) + ptr[5];
                    alpha = (alpha<<8) +  ptr[7];
                }

                ycbcr = iftRGBtoYCbCr(rgb, depth==8?255:65535);

                img->val[p] = ycbcr.val[0];
                img->Cb[p] = ycbcr.val[1];
                img->Cr[p] = ycbcr.val[2];
                img->alpha[p] = alpha;

                p++;
            }
        }
    }

    for (y = 0; y < height; ++y) {
        iftFree(row_pointers[y]);
    }

    iftFree(row_pointers);

    png_destroy_read_struct(&png_ptr, &info_ptr, NULL);

    img->dz = 0.0;

    return img;
}

iftImage* iftReadImageJPEG(const char* format, ...) 
{
    va_list args;
    char filename[IFT_STR_DEFAULT_SIZE];

    va_start(args, format);
    vsprintf(filename, format, args);
    va_end(args);

    iftImage* image = NULL;
    //code based on externals/libjpeg/source/example.c
    /* This struct contains the JPEG decompression parameters and pointers to
* working space (which is allocated as needed by the JPEG library).
*/
    struct jpeg_decompress_struct cinfo;
    /* We use our private extension JPEG error handler.
* Note that this struct must live as long as the main JPEG parameter
* struct, to avoid dangling-pointer problems.
*/
    struct jpeg_error_mgr jerr;

    /* More stuff */
    FILE * infile;		/* source file */
    JSAMPARRAY buffer;		/* Output row buffer */
    int row_stride;		/* physical row width in output buffer */
    /* In this example we want to open the input file before doing anything else,
 * so that the setjmp() error recovery below can assume the file is open.
 * VERY IMPORTANT: use "b" option to fopen() if you are on a machine that
 * requires it in order to read binary files.
 */

    if ((infile = fopen(filename, "rb")) == NULL) {
        printf("[readImageJPEG] can't open %s\n",filename);
        return NULL;
    }

    /* Step 1: allocate and initialize JPEG decompression object */

    /* We set up the normal JPEG error routines, then override error_exit. */
    cinfo.err = jpeg_std_error(&jerr);
    //jerr.pub.error_exit = my_error_exit;

    /* Establish the setjmp return context for my_error_exit to use. */
    jmp_buf setjmp_buffer;
    if (setjmp(setjmp_buffer)) {
        /* If we get here, the JPEG code has signaled an error.
         * We need to clean up the JPEG object, close the input file, and return.
         */
        jpeg_destroy_decompress(&cinfo);
        printf("[readImageJPEG] code has signaled an error\n");
        fclose(infile);
        return NULL;
    }

    /* Now we can initialize the JPEG decompression object. */
    jpeg_create_decompress(&cinfo);
    /* Step 2: specify data source (eg, a file) */
    jpeg_stdio_src(&cinfo, infile);
    /* Step 3: read file parameters with jpeg_read_header() */
    (void) jpeg_read_header(&cinfo, TRUE);
    /* We can ignore the return value from jpeg_read_header since
     *   (a) suspension is not possible with the stdio data source, and
     *   (b) we passed TRUE to reject a tables-only JPEG file as an error.
     * See libjpeg.txt for more info.
     */

    /* Step 4: set parameters for decompression */

    /* In this example, we don't need to change any of the defaults set by
     * jpeg_read_header(), so we do nothing here.
     */

    /* Step 5: Start decompressor */
    (void) jpeg_start_decompress(&cinfo);
    /* We can ignore the return value since suspension is not possible
 * with the stdio data source.
 */

    /* We may need to do some setup of our own at this point before reading
 * the data.  After jpeg_start_decompress() we have the correct scaled
 * output image dimensions available, as well as the output colormap
 * if we asked for color quantization.
 * In this example, we need to make an output work buffer of the right size.
 */
    /* JSAMPLEs per row in output buffer */
    row_stride = cinfo.output_width * cinfo.output_components;
    /* Make a one-row-high sample array that will go away when done with image */
    buffer = (*cinfo.mem->alloc_sarray)
            ((j_common_ptr) &cinfo, JPOOL_IMAGE, row_stride, 1);


    /* Step 6: while (scan lines remain to be read) */
    /*           jpeg_read_scanlines(...); */

    /* Here we use the library's state variable cinfo.output_scanline as the
     * loop counter, so that we don't have to keep track ourselves.
     */
    image = iftCreateImage(cinfo.output_width,cinfo.output_height,1);

    //0 - JCS_GRAYSCALE
    //1 - JCS_RGB
    //2 - JCS_YCbCr
    //3 - JCS_CMYK
    //4 - JCS_YCCK
    //5 - JCS_BG_RGB
    //6 - JCS_BG_YCC
    unsigned  int imageRow = 0;
    unsigned int imageCol = 0;
    iftColor rgb;
    iftColor YCbCr;
    float scalingFactor = pow(2,cinfo.data_precision)-1;
    switch (cinfo.out_color_space){
        case JCS_GRAYSCALE:

            imageRow = 0;
            while (cinfo.output_scanline < cinfo.output_height){
                jpeg_read_scanlines(&cinfo, buffer, 1);
                imageCol = 0;
                for (unsigned int i = 0; i < (unsigned int)cinfo.output_width; i++) {
                    int shift = i*cinfo.num_components;
                    iftImgVal(image,imageCol,imageRow,0) = buffer[0][shift];
                    imageCol++;
                }
                imageRow++;
            }
            break;

        case JCS_RGB:

            iftSetCbCr(image,128);
            imageRow = 0;
            while (cinfo.output_scanline < cinfo.output_height){
                jpeg_read_scanlines(&cinfo, buffer, 1);
                imageCol = 0;
                for (unsigned int i = 0; i < (unsigned int)cinfo.output_width; i++) {
                    int shift = i*cinfo.num_components;
                    rgb.val[0] = buffer[0][shift+0];
                    rgb.val[1] = buffer[0][shift+1];
                    rgb.val[2] = buffer[0][shift+2];
                    YCbCr = iftRGBtoYCbCr(rgb,scalingFactor);
                    iftImgVal(image,imageCol,imageRow,0) = YCbCr.val[0];
                    iftImgCb(image,imageCol,imageRow,0) = YCbCr.val[1];
                    iftImgCr(image,imageCol,imageRow,0) = YCbCr.val[2];
                    imageCol++;
                }
                imageRow++;
            }

            break;

        case JCS_YCbCr:

            iftSetCbCr(image,128);
            imageRow = 0;
            while (cinfo.output_scanline < cinfo.output_height){
                jpeg_read_scanlines(&cinfo, buffer, 1);
                imageCol = 0;
                for (unsigned int i = 0; i < (unsigned int)cinfo.output_width; i++) {
                    int shift = i*cinfo.num_components;
                    iftImgVal(image,imageCol,imageRow,0) = buffer[0][shift+0];
                    iftImgCb(image,imageCol,imageRow,0) = buffer[0][shift+1];
                    iftImgCr(image,imageCol,imageRow,0) = buffer[0][shift+2];
                    imageCol++;
                }
                imageRow++;
            }

            break;
        case JCS_CMYK:

            iftSetCbCr(image,128);
            imageRow = 0;
            while (cinfo.output_scanline < cinfo.output_height){
                jpeg_read_scanlines(&cinfo, buffer, 1);
                imageCol = 0;
                for (unsigned int i = 0; i < (unsigned int)cinfo.output_width; i++) {
                    int shift = i*cinfo.num_components;
                    //convert CMYK to RGB (reference: http://www.rapidtables.com/convert/color/cmyk-to-rgb.htm)
                    rgb.val[0] = 255*(100-buffer[0][shift+0])*(100-buffer[0][shift+3]);
                    rgb.val[1] = 255*(100-buffer[0][shift+1])*(100-buffer[0][shift+3]);;
                    rgb.val[2] = 255*(100-buffer[0][shift+2])*(100-buffer[0][shift+3]);;
                    YCbCr = iftRGBtoYCbCr(rgb,scalingFactor);
                    iftImgVal(image,imageCol,imageRow,0) = YCbCr.val[0];
                    iftImgCb(image,imageCol,imageRow,0) = YCbCr.val[1];
                    iftImgCr(image,imageCol,imageRow,0) = YCbCr.val[2];
                    imageCol++;
                }
                imageRow++;
            }

            break;
        case JCS_YCCK:

            iftSetCbCr(image,128);
            imageRow = 0;
            iftWarning("Image is Y/Cb/Cr/K color space. The channel K is ignored", "iftReadImageJPEG");
            while (cinfo.output_scanline < cinfo.output_height){
                jpeg_read_scanlines(&cinfo, buffer, 1);
                imageCol = 0;
                for (unsigned int i = 0; i < (unsigned int)cinfo.output_width; i++) {
                    int shift = i*cinfo.num_components;
                    iftImgVal(image,imageCol,imageRow,0) = buffer[0][shift+0];
                    iftImgCb(image,imageCol,imageRow,0) = buffer[0][shift+1];
                    iftImgCr(image,imageCol,imageRow,0) = buffer[0][shift+2];
                    imageCol++;
                }
                imageRow++;
            }

            break;
        /*case JCS_BG_RGB:
    
            iftError("Big gamut red/green/blue color space not supported", "iftReadImageJPEG");

            break;
        case JCS_BG_YCC:
    
            iftError("Big gamut red/green/blue color space not supported", "iftReadImageJPEG");

            break;*/
        default:
    
            iftError("Unkwon color space", "iftReadImageJPEG");

            break;
    }

    /* Step 7: Finish decompression */
    (void) jpeg_finish_decompress(&cinfo);

    /* This is an important step since it will release a good deal of memory. */
    jpeg_destroy_decompress(&cinfo);

    /* After finish_decompress, we can close the input file.
     * Here we postpone it until after no more JPEG errors are possible,
     * so as to simplify the setjmp error logic above.  (Actually, I don't
     * think that jpeg_destroy can do an error exit, but why assume anything...)
     */

    fclose(infile);
    /* At this point you may want to check to see whether any corrupt-data
     * warnings occurred (test whether jerr.pub.num_warnings is nonzero).
     */

    //jerr.num_warnings; //useful to know about corrupted data
    //printf("%ld\n",jerr.num_warnings);

    return image;
}

iftImage *iftReadImageGZip(const char *format, ...) 
{
    iftImage *img    = NULL;
    iftGZipFile fp   = NULL;
    uchar    *data8  = NULL;
    ushort   *data16 = NULL;
    int      *data32 = NULL;
    char     type[10], header[IFT_STR_DEFAULT_SIZE];
    int      p, v, xsize, ysize, zsize;
    float    dx, dy, dz;

    va_list args;
    char    filename[IFT_STR_DEFAULT_SIZE];

    va_start(args, format);
    vsprintf(filename, format, args);
    va_end(args);

    fp = iftGZipOpen(filename, "rb", 1);

    if (fp == NULL) {
        iftError("Cannot open file: \"%s\"", "iftReadImageGZip", filename);
    }

    // Reading the header (everything until a \n is found)
    iftGZipGets(header, IFT_STR_DEFAULT_SIZE, fp);

    // Reading all info from the header
    if (sscanf(header, "%s %d %d %d %f %f %f %d", type, &xsize, &ysize, &zsize, &dx, &dy, &dz, &v) != 8)
        iftError("Reading error! Image header does not match what is expected: %s", "iftReadImageGZip", header);

    if (iftCompareStrings(type, "SCN")) {
        size_t nread;

        img = iftCreateImage(xsize, ysize, zsize);
        img->dx = dx;
        img->dy = dy;
        img->dz = dz;

        if (v == 8) {
            data8 = iftAllocUCharArray(img->n);
            if ((nread = iftGZipRead(data8, sizeof(uchar), img->n, fp)) != img->n)
                iftError("Reading error 8 bits. %lu voxels read when %d were expected", "iftReadImageGZip", nread,
                         img->n);
            for (p = 0; p < img->n; p++)
                img->val[p] = (int) data8[p];
            iftFree(data8);
        } else if (v == 16) {
            data16 = iftAllocUShortArray(img->n);
            if ((nread = iftGZipRead(data16, sizeof(ushort), img->n, fp)) != img->n)
                iftError("Reading error 16 bits. %lu voxels read when %d were expected", "iftReadImageGZip", nread,
                         img->n);
            for (p = 0; p < img->n; p++)
                img->val[p] = (int) data16[p];
            iftFree(data16);
        } else if (v == 32) {
            data32 = iftAllocIntArray(img->n);
            if ((nread = iftGZipRead(data32, sizeof(int), img->n, fp)) != img->n)
                iftError("Reading error 32 bits. %lu voxels read when %d were expected", "iftReadImageGZip", nread,
                         img->n);
            for (p = 0; p < img->n; p++)
                img->val[p] = data32[p];
            iftFree(data32);
        } else {
            iftError("Input scene must be 8, 16, or 32 bit", "iftReadImageGZip");
        }
    } else {
        iftError("Invalid file type", "iftReadImageGZip");
    }

    iftGZipClose(&fp);

    return (img);
}

iftImage *iftReadImageP5(const char *format, ...) 
{
    iftImage *img    = NULL;
    FILE     *fp     = NULL;
    uchar    *data8  = NULL;
    ushort   *data16 = NULL;
    char     type[10];
    int      p, v, xsize, ysize, zsize, hi, lo;

    va_list args;
    char    filename[IFT_STR_DEFAULT_SIZE];

    va_start(args, format);
    vsprintf(filename, format, args);
    va_end(args);

    fp = fopen(filename, "rb");
    if (fp == NULL) {
        iftError(MSG_FILE_OPEN_ERROR, "iftReadImageP5", filename);
    }

    if (fscanf(fp, "%s\n", type) != 1) {
        iftError("Reading error", "iftReadImageP5");
    }

    if (iftCompareStrings(type, "P5")) {

        iftSkipComments(fp);

        if (fscanf(fp, "%d %d\n", &xsize, &ysize) != 2)
            iftError("Reading error", "iftReadImageP5");
        zsize = 1;

        img = iftCreateImage(xsize, ysize, zsize);
        img->dz = 0.0;

        if (fscanf(fp, "%d", &v) != 1)
            iftError("Reading error", "iftReadImageP5");

        while (fgetc(fp) != '\n');

        if ((v <= 255) && (v > 0)) {
            data8 = iftAllocUCharArray(img->n);

            if (fread(data8, sizeof(uchar), img->n, fp) != img->n)
                iftError("Reading error", "iftReadImageP5");

            for (p = 0; p < img->n; p++)
                img->val[p] = (int) data8[p];

            iftFree(data8);

        } else if ((v <= 65535) && (v > 255)) {
            data16 = iftAllocUShortArray(img->n);

            for (p = 0; p < img->n; p++) {
                if ((hi = fgetc(fp)) == EOF)
                    iftError("Reading error", "iftReadImageP5");
                if ((lo = fgetc(fp)) == EOF)
                    iftError("Reading error", "iftReadImageP5");

                data16[p] = (hi << 8) + lo;
            }

            for (p = 0; p < img->n; p++)
                img->val[p] = (int) data16[p];

            iftFree(data16);

        } else {
            iftError("Invalid maximum value", "iftReadImageP5");
        }
    } else {
        iftError("Invalid image type", "iftReadImageP5");
    }

    fclose(fp);
    return (img);
}

iftImage *iftReadImageP6(const char *format, ...)
{
    iftImage  *img=NULL;
    FILE    *fp=NULL;
    char    type[10];
    int     p,v,xsize,ysize,zsize;
    ushort rgb16[3];
    iftColor RGB,YCbCr;

    va_list args;
    char filename[IFT_STR_DEFAULT_SIZE];

    va_start(args, format);
    vsprintf(filename, format, args);
    va_end(args);

    fp = fopen(filename,"r");
    if (fp == NULL){
        iftError(MSG_FILE_OPEN_ERROR, "iftReadImageP6", filename);
    }

    if(fscanf(fp,"%s\n",type)!=1)
        iftError("Reading error", "iftReadImageP6");
    if(iftCompareStrings(type,"P6")){

        iftSkipComments(fp);

        if(fscanf(fp,"%d %d\n",&xsize,&ysize)!=2)
            iftError("Reading error", "iftReadImageP6");

        zsize = 1;
        img = iftCreateImage(xsize,ysize,zsize);
        img->Cb = iftAllocUShortArray(img->n);
        img->Cr = iftAllocUShortArray(img->n);
        img->dz = 0.0;
        if (fscanf(fp,"%d",&v)!=1)
            iftError("Reading error", "iftReadImageP6");

        while(fgetc(fp) != '\n');

        if (v >= 0 && v < 256) {
            for (p=0; p < img->n; p++) {
                RGB.val[0] = fgetc(fp);
                RGB.val[1] = fgetc(fp);
                RGB.val[2] = fgetc(fp);
                YCbCr      = iftRGBtoYCbCr(RGB,255);
                img->val[p]=YCbCr.val[0];
                img->Cb[p] =(ushort)YCbCr.val[1];
                img->Cr[p] =(ushort)YCbCr.val[2];
            }
        } else if (v >= 256 && v <= 65536) {

            int rgbBitDepth = ceil(iftLog(v, 2));
            int ycbcrBitDepth = rgbBitDepth;

            if(ycbcrBitDepth<10)
                ycbcrBitDepth = 10;
            else if(ycbcrBitDepth < 12)
                ycbcrBitDepth = 12;
            else if(ycbcrBitDepth < 16)
                ycbcrBitDepth = 16;

            for (p=0; p < img->n; p++) {
                // read 6 bytes for each image pixel
                if (fread(rgb16, 2, 3, fp) == 3) {
                    // the PPM format specifies 2-byte integers as big endian,
                    // so we need to swap the bytes if the architecture is little endian
                    RGB.val[0]  = ((rgb16[0] & 0xff) << 8) | ((ushort) rgb16[0] >> 8);
                    RGB.val[1]  = ((rgb16[1] & 0xff) << 8) | ((ushort) rgb16[1] >> 8);
                    RGB.val[2]  = ((rgb16[2] & 0xff) << 8) | ((ushort) rgb16[2] >> 8);
                    YCbCr       = iftRGBtoYCbCrBT2020(RGB, rgbBitDepth, ycbcrBitDepth);
//                    YCbCr = iftRGBtoYCbCr(RGB, v);
                    img->val[p] = YCbCr.val[0];
                    img->Cb[p]  = (ushort)YCbCr.val[1];
                    img->Cr[p]  = (ushort)YCbCr.val[2];
                }
            }
        } else {
            iftError("Invalid maximum value", "iftReadImageP6");
        }
    }else{
        iftError("Invalid image type", "iftReadImageP6");
    }

    fclose(fp);
    return(img);
}

iftImage *iftReadImageP2(const char *format, ...) 
{
    iftImage *img = NULL;
    FILE     *fp  = NULL;
    char     type[10];
    int      p, v, xsize, ysize, zsize;

    va_list args;
    char    filename[IFT_STR_DEFAULT_SIZE];

    va_start(args, format);
    vsprintf(filename, format, args);
    va_end(args);

    fp = fopen(filename, "r");

    if (fp == NULL) {
        iftError(MSG_FILE_OPEN_ERROR, "iftReadImageP2", filename);
    }

    if (fscanf(fp, "%s\n", type) != 1)
        iftError("Reading error", "iftReadImageP2");

    if (iftCompareStrings(type, "P2")) {

        iftSkipComments(fp);

        if (fscanf(fp, "%d %d\n", &xsize, &ysize) != 2)
            iftError("Reading error", "iftReadImageP2");
        zsize = 1;
        img   = iftCreateImage(xsize, ysize, zsize);
        img->dz = 0.0;
        if (fscanf(fp, "%d", &v) != 1)
            iftError("Reading error", "iftReadImageP2");

        while (fgetc(fp) != '\n');

        for (p = 0; p < img->n; p++)
            if (fscanf(fp, "%d", &img->val[p]) != 1)
                iftError("Reading error", "iftReadImageP2");

    } else {
        iftError("Invalid image type", "iftReadImageP2");
    }

    fclose(fp);
    return (img);
}

iftImage *iftCreateImageFromBuffer(int xsize, int ysize, int zsize, int *val) 
{
    iftImage *img = NULL;
    int      y, z, xysize;

    img = (iftImage *) iftAlloc(1, sizeof(iftImage));
    if (img == NULL) {
        iftError(MSG_MEMORY_ALLOC_ERROR, "iftCreateImage");
    }

    img->val   = val;
    img->Cb    = img->Cr = NULL;
    img->alpha = NULL;
    img->xsize = xsize;
    img->ysize = ysize;
    img->zsize = zsize;
    img->dx    = 1.0;
    img->dy    = 1.0;
    img->dz    = 1.0;
    img->tby   = iftAllocIntArray(ysize);
    img->tbz   = iftAllocIntArray(zsize);
    img->n     = xsize * ysize * zsize;

    if (img->val == NULL || img->tbz == NULL || img->tby == NULL) {
        iftError(MSG_MEMORY_ALLOC_ERROR, "iftCreateImage");
    }

    img->tby[0] = 0;
    for (y = 1; y < ysize; y++)
        img->tby[y] = img->tby[y - 1] + xsize;

    img->tbz[0] = 0;
    xysize = xsize * ysize;
    for (z = 1; z < zsize; z++)
        img->tbz[z] = img->tbz[z - 1] + xysize;

    return (img);
}

void iftCopyImageInplace(const iftImage *src, iftImage *dest) 
{
    int p;

    iftVerifyImageDomains(src, dest, "iftCopyImageInplace");

    iftCopyVoxelSize(src, dest);

    for (p=0; p < src->n; p++)
        dest->val[p]= src->val[p];

    if (src->Cb != NULL) {
        if(dest->Cb == NULL)
            dest->Cb = iftAllocUShortArray(src->n);
        if(dest->Cr == NULL)
            dest->Cr = iftAllocUShortArray(src->n);
        for (p=0; p < src->n; p++) {
            dest->Cb[p]= src->Cb[p];
            dest->Cr[p]= src->Cr[p];
        }
    }
}

void iftWriteImage(const iftImage *img, const char *format, ...) 
{
    FILE   *fp     = NULL;
    int    p;
    uchar  *data8  = NULL;
    ushort *data16 = NULL;
    int    *data32 = NULL;

    va_list args;
    char    filename[IFT_STR_DEFAULT_SIZE];

    va_start(args, format);
    vsprintf(filename, format, args);
    va_end(args);

    int img_min_val = iftMinimumValue(img);
    int img_max_val = iftMaximumValue(img);

    if (img_min_val < 0) {
        char msg[200];
        sprintf(msg, "Shifting image values from [%d,%d] to [%d,%d] on the original image\n",
                img_min_val, img_max_val, 0, img_max_val - img_min_val);
        iftWarning(msg, "iftWriteImage");
        for (p = 0; p < img->n; p++)
            img->val[p] = img->val[p] - img_min_val;
        img_max_val = img_max_val - img_min_val;
    }

    fp = fopen(filename, "wb");
    if (fp == NULL)
        iftError("Cannot open file: \"%s\"", "iftWriteImage", filename);

    fprintf(fp, "SCN\n");
    fprintf(fp, "%d %d %d\n", img->xsize, img->ysize, img->zsize);
    fprintf(fp, "%f %f %f\n", img->dx, img->dy, img->dz);


    if (img_max_val < 256) {
        fprintf(fp, "%d\n", 8);
        data8 = iftAllocUCharArray(img->n);
        for (p = 0; p < img->n; p++)
            data8[p] = (uchar) img->val[p];
        fwrite(data8, sizeof(uchar), img->n, fp);
        iftFree(data8);
    } else if (img_max_val < 65536) {
        fprintf(fp, "%d\n", 16);
        data16 = iftAllocUShortArray(img->n);
        for (p = 0; p < img->n; p++)
            data16[p] = (ushort) img->val[p];
        fwrite(data16, sizeof(ushort), img->n, fp);

        iftFree(data16);
    } else if (img_max_val < IFT_INFINITY_INT) {
        fprintf(fp, "%d\n", 32);
        data32 = iftAllocIntArray(img->n);
        for (p = 0; p < img->n; p++)
            data32[p] = img->val[p];
        fwrite(data32, sizeof(int), img->n, fp);
        iftFree(data32);
    }

    fclose(fp);
}

void iftWriteImageGZip(const iftImage *img, const char *format, ...) 
{
    iftGZipFile fp = NULL;
    int    p;
    uchar  *data8  = NULL;
    ushort *data16 = NULL;
    int    *data32 = NULL;
    va_list args;
    char    filename[IFT_STR_DEFAULT_SIZE], header[IFT_STR_DEFAULT_SIZE];

    va_start(args, format);
    vsprintf(filename, format, args);
    va_end(args);

    int img_min_val = iftMinimumValue(img);
    int img_max_val = iftMaximumValue(img);

    if (img_min_val < 0) {
        char msg[200];
        sprintf(msg, "Shifting image values from [%d,%d] to [%d,%d] on the original image\n",
                img_min_val, img_max_val, 0, img_max_val - img_min_val);
        iftWarning(msg, "iftWriteImage");
        for (p = 0; p < img->n; p++)
            img->val[p] = img->val[p] - img_min_val;
        img_max_val = img_max_val - img_min_val;
    }

    fp = iftGZipOpen(filename, "wb", 1);
    if (fp == NULL)
        iftError("Cannot open file: \"%s\"", "iftWriteImageGZip", filename);


    if (img_max_val < 256) {
        // The header is saved with a single line because we use iftGZipGets to read it in iftReadImageGZip, and that
        // function reads everything until a '\n' is found. If the gzip file is uncompressed with an external tool,
        // iftReadImage may still be used to read it
        sprintf(header, "SCN %d %d %d %f %f %f %d\n", img->xsize, img->ysize, img->zsize,
                img->dx, img->dy, img->dz, 8);

        // Writes everything including the '\n', which is important for reading the header with iftGZipGets
        iftGZipPuts(header, fp);

        data8 = iftAllocUCharArray(img->n);
        for (p = 0; p < img->n; p++)
            data8[p] = (uchar) img->val[p];
        iftGZipWrite(data8, sizeof(uchar), img->n, fp);
        iftFree(data8);
    } else if (img_max_val < 65536) {
        // The header is saved with a single line because we use iftGZipGets to read it in iftReadImageGZip, and that
        // function reads everything until a '\n' is found. If the gzip file is uncompressed with an external tool,
        // iftReadImage may still be used to read it
        sprintf(header, "SCN %d %d %d %f %f %f %d\n", img->xsize, img->ysize, img->zsize,
                img->dx, img->dy, img->dz, 16);

        // Writes everything including the '\n', which is important for reading the header with iftGZipGets
        iftGZipPuts(header, fp);

        data16 = iftAllocUShortArray(img->n);
        for (p = 0; p < img->n; p++)
            data16[p] = (ushort) img->val[p];
        iftGZipWrite(data16, sizeof(ushort), img->n, fp);

        iftFree(data16);
    } else if (img_max_val < IFT_INFINITY_INT) {
        // The header is saved with a single line because we use iftGZipGets to read it in iftReadImageGZip, and that
        // function reads everything until a '\n' is found. If the gzip file is uncompressed with an external tool,
        // iftReadImage may still be used to read it
        sprintf(header, "SCN %d %d %d %f %f %f %d\n", img->xsize, img->ysize, img->zsize,
                img->dx, img->dy, img->dz, 32);
        // Writes everything including the '\n', which is important for reading the header with iftGZipGets
        iftGZipPuts(header, fp);

        data32 = iftAllocIntArray(img->n);
        for (p = 0; p < img->n; p++)
            data32[p] = img->val[p];
        iftGZipWrite(data32, sizeof(int), img->n, fp);
        iftFree(data32);
    }

    iftGZipClose(&fp);
}

void iftWriteImageP5(const iftImage *img, const char *format, ...) 
{
    FILE   *fp     = NULL;
    int    p, hi, lo;
    uchar  *data8  = NULL;
    ushort *data16 = NULL;

    va_list args;
    char    filename[IFT_STR_DEFAULT_SIZE];

    va_start(args, format);
    vsprintf(filename, format, args);
    va_end(args);

    fp = fopen(filename, "wb");
    if (fp == NULL)
        iftError(MSG_FILE_OPEN_ERROR, "iftWriteImageP5", filename);

    fprintf(fp, "P5\n");
    fprintf(fp, "%d %d\n", img->xsize, img->ysize);

    int img_max_val = iftMaximumValue(img);
    int img_min_val = iftMinimumValue(img);

    if ((img_max_val < 256) && (img_min_val >= 0)) {
        fprintf(fp, "%d\n", 255);
        data8 = iftAllocUCharArray(img->n);
        for (p = 0; p < img->n; p++)
            data8[p] = (uchar) img->val[p];
        fwrite(data8, sizeof(uchar), img->n, fp);
        iftFree(data8);
    } else if (img_max_val < 65536) {
        fprintf(fp, "%d\n", 65535);
        data16 = iftAllocUShortArray(img->n);
        for (p = 0; p < img->n; p++)
            data16[p] = (ushort) img->val[p];

        {
#define HI(num) (((num) & 0x0000FF00) >> 8)
#define LO(num) ((num) & 0x000000FF)
            for (p = 0; p < img->n; p++) {
                hi = HI(data16[p]);
                lo = LO(data16[p]);
                fputc(hi, fp);
                fputc(lo, fp);
            }
        }

        iftFree(data16);
    } else {
        char msg[200];
        sprintf(msg, "Cannot write image as P5 (%d/%d)", img_max_val, img_min_val);
        iftError(msg, "iftWriteImageP5");
    }
    fclose(fp);
}

void iftWriteImageP6(const iftImage *img, const char *format, ...) 
{
    FILE     *fp = NULL;
    int      p;
    ushort   rgb16[3];
    iftColor YCbCr, RGB;

    if (!iftIsColorImage(img))
        iftError("Image is not colored", "iftWriteImageP6");

    va_list args;
    char    filename[IFT_STR_DEFAULT_SIZE];

    va_start(args, format);
    vsprintf(filename, format, args);
    va_end(args);

    fp = fopen(filename, "w");
    if (fp == NULL)
        iftError(MSG_FILE_OPEN_ERROR, "iftWriteImageP6", filename);

    fprintf(fp, "P6\n");
    fprintf(fp, "%d %d\n", img->xsize, img->ysize);

    int img_max_val = iftMaximumValue(img);
    int img_min_val = iftMinimumValue(img);

    if (img_min_val < 0) {
        iftError("Cannot write image as P6", "iftWriteImageP6");
    }
    if (img_max_val < 256) {
        fprintf(fp, "%d\n", 255);
        for (p = 0; p < img->n; p++) {
            YCbCr.val[0] = img->val[p];
            YCbCr.val[1] = img->Cb[p];
            YCbCr.val[2] = img->Cr[p];

            RGB = iftYCbCrtoRGB(YCbCr, 255);

            fputc(((uchar) RGB.val[0]), fp);
            fputc(((uchar) RGB.val[1]), fp);
            fputc(((uchar) RGB.val[2]), fp);
        }
    } else if (img_max_val < 65536) {
//        int rgbBitDepth = 9;
//        // find the bit depth for the maximum value img_max_val
//        while ((1 << rgbBitDepth) <= img_max_val) {
//            rgbBitDepth++;
//        }

        int rgbBitDepth = ceil(iftLog(img_max_val, 2));

        fprintf(fp, "%d\n", (1 << rgbBitDepth) - 1);
        for (p = 0; p < img->n; p++) {
            YCbCr.val[0] = img->val[p];
            YCbCr.val[1] = img->Cb[p];
            YCbCr.val[2] = img->Cr[p];
            RGB = iftYCbCrBT2020toRGB(YCbCr, rgbBitDepth, rgbBitDepth);
//            RGB = iftYCbCrtoRGB(YCbCr, img_max_val);
            // the PPM format specifies 2-byte integers as big endian,
            // so we need to swap the bytes if the architecture is little endian
            rgb16[0] = ((RGB.val[0] & 0xff) << 8) | ((ushort) RGB.val[0] >> 8);
            rgb16[1] = ((RGB.val[1] & 0xff) << 8) | ((ushort) RGB.val[1] >> 8);
            rgb16[2] = ((RGB.val[2] & 0xff) << 8) | ((ushort) RGB.val[2] >> 8);
            // write 6 bytes for each image pixel
            if (fwrite(rgb16, 2, 3, fp) != 3) {
                iftError("Cannot write 16-bit image as P6", "iftWriteImageP6");
            }
        }
    } else {
        iftError("Cannot write image as P6", "iftWriteImageP6");
    }
    fclose(fp);
}

void iftWriteImageP2(const iftImage *img, const char *format, ...) 
{
    FILE *fp = NULL;
    int  p;

    va_list args;
    char    filename[IFT_STR_DEFAULT_SIZE];
    int     depth = iftImageDepth(img);
      
    va_start(args, format);
    vsprintf(filename, format, args);
    va_end(args);

    fp = fopen(filename, "w");
    if (fp == NULL)
        iftError(MSG_FILE_OPEN_ERROR, "iftWriteImageP2", filename);

    fprintf(fp, "P2\n");
    fprintf(fp, "%d %d\n", img->xsize, img->ysize);

    int img_max_val = iftMaxImageRange(depth);
    fprintf(fp, "%d\n", img_max_val);
    for (p = 0; p < img->n; p++) {
        fprintf(fp, "%d ", img->val[p]);
        if (iftGetXCoord(img, p) == (img->xsize - 1)) fprintf(fp, "\n");
    }

    fclose(fp);
}

void iftWritePngImageAux(const char *file_name, png_bytep *row_pointers, int width, int height, int bit_depth, int color_type) 
{

    /* create file */
    FILE *fp = fopen(file_name, "wb");
    if (!fp)
        iftError("Internal Error: File %s could not be opened for writing", "iftWritePngImageAux", file_name);


    /* initialize stuff */
    png_structp png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);

    if (!png_ptr)
        iftError("Internal Error: png_create_write_struct failed", "iftWriteImagePNG");

    png_infop info_ptr = png_create_info_struct(png_ptr);
    if (!info_ptr)
        iftError("Internal Error: png_create_info_struct failed", "iftWriteImagePNG");

    if (setjmp(png_jmpbuf(png_ptr)))
        iftError("Internal Error: Error during init_io", "iftWriteImagePNG");

    png_init_io(png_ptr, fp);


    /* write header */
    if (setjmp(png_jmpbuf(png_ptr)))
        iftError("Internal Error: Error during writing header", "iftWriteImagePNG");

    png_set_IHDR(png_ptr, info_ptr, width, height,
                 bit_depth, color_type, PNG_INTERLACE_NONE,
                 PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);

    png_write_info(png_ptr, info_ptr);


    /* write bytes */
    if (setjmp(png_jmpbuf(png_ptr)))
        iftError("Internal Error: Error during writing bytes", "iftWriteImagePNG");

    png_write_image(png_ptr, row_pointers);


    /* end write */
    if (setjmp(png_jmpbuf(png_ptr)))
        iftError("Internal Error: Error during end of write", "iftWriteImagePNG");

    png_write_end(png_ptr, NULL);

    png_destroy_write_struct(&png_ptr, &info_ptr);

    /* cleanup heap allocation */
    for (int y=0; y<height; y++)
        iftFree(row_pointers[y]);
    iftFree(row_pointers);

    fclose(fp);
}

void iftWriteImagePNG(const iftImage* img, const char* format, ...) 
{
    png_bytep *row_pointers;
    int width, height, depth, byteshift;

    va_list args;
    char filename[IFT_STR_DEFAULT_SIZE];

    va_start(args, format);
    vsprintf(filename, format, args);
    va_end(args);

    width = img->xsize;
    height = img->ysize;
    png_byte color_type;
    depth = iftImageDepth(img);

    if(depth<=8) {
        depth = 8;
    } else {
        depth = 16;
    }

    byteshift = depth/8;
    //int offset = depth==16?1:0;//to read the second byte first in cases of 16bit images

    size_t numberOfChannels=1;
    if(iftIsColorImage(img)){
        if(img->alpha == NULL){
            numberOfChannels = 3;//RGB
            color_type = PNG_COLOR_TYPE_RGB;
        }else{
            numberOfChannels = 4;//RGB_ALPHA
            color_type = PNG_COLOR_TYPE_RGB_ALPHA;
        }
    }else{
        if(img->alpha == NULL){
            numberOfChannels = 1;//GRAY
            color_type = PNG_COLOR_TYPE_GRAY;
        }else{
            numberOfChannels = 2;//GRAY_ALPHA
            color_type = PNG_COLOR_TYPE_GRAY_ALPHA;
        }
    }

    //size_t pixel_size = (iftIsColorImage(img)?3:1 ) * byteshift;
    row_pointers = (png_bytep*) iftAlloc(height, sizeof(png_bytep));
    for (int y=0; y<height; y++)
        row_pointers[y] = (png_byte*) iftAlloc(width, numberOfChannels*byteshift);

    if(color_type == PNG_COLOR_TYPE_GRAY){
        int p = 0;
        for (int y = 0; y < height; ++y) {
            png_byte* row = row_pointers[y];
            for (int x=0; x<width; x++) {
                png_byte* ptr = &(row[x*numberOfChannels*byteshift]);

                ptr[0] = img->val[p] & 0xFF;//get first byte

                if(depth==16) {//in 16bit image, we should store as big endian
                    ptr[1] = ptr[0];
                    ptr[0] = (img->val[p]>>8) & 0xFF;//get second byte
                }

                p++;
            }
        }
    }else if(color_type == PNG_COLOR_TYPE_GRAY_ALPHA){
        int p = 0;
        for (int y = 0; y < height; ++y) {
            png_byte* row = row_pointers[y];
            for (int x=0; x<width; x++) {
                png_byte* ptr = &(row[x*numberOfChannels*byteshift]);

                if(depth==8){
                    ptr[0] = img->val[p] & 0xFF;//get first byte
                    ptr[1] = img->alpha[p] & 0xFF;//get second byte
                }


                if(depth==16) {//in 16bit image, we should store as big endian
                    ptr[0] = img->val[p]>>8;//get first byte
                    ptr[1] = img->val[p] & 0xFF;//get second byte


                    ptr[2] = img->alpha[p]>>8;//get first byte;
                    ptr[3] = img->alpha[p] & 0xFF;;//get second byte
                }
                p++;
            }
        }
    }else if(color_type == PNG_COLOR_TYPE_RGB){
        iftColor rgb, ycbcr;
        int p = 0;
        for (int y = 0; y < height; ++y) {
            png_byte* row = row_pointers[y];
            for (int x=0; x<width; x++) {
                png_byte* ptr = &(row[x*numberOfChannels*byteshift]);

                ycbcr.val[0] = img->val[p];
                ycbcr.val[1] = img->Cb[p];
                ycbcr.val[2] = img->Cr[p];

                rgb = iftYCbCrBT2020toRGB(ycbcr, depth, depth);

                ptr[0*byteshift] = rgb.val[0] & 0xFF;//get first byte
                ptr[1*byteshift] = rgb.val[1] & 0xFF;
                ptr[2*byteshift] = rgb.val[2] & 0xFF;

                if(depth==16) {//in 16bit image, we should store as big endian
                    ptr[(0*byteshift)+1] = ptr[0*byteshift];
                    ptr[(1*byteshift)+1] = ptr[1*byteshift];
                    ptr[(2*byteshift)+1] = ptr[2*byteshift];

                    ptr[0*byteshift] = ((rgb.val[0]>>8) & 0xFF);//get second byte
                    ptr[1*byteshift] = ((rgb.val[1]>>8) & 0xFF);
                    ptr[2*byteshift] = ((rgb.val[2]>>8) & 0xFF);
                }

                p++;
            }
        }

    }else if(color_type == PNG_COLOR_TYPE_RGB_ALPHA){
        iftColor rgb, ycbcr;
        int p = 0;
        for (int y = 0; y < height; ++y) {
            png_byte* row = row_pointers[y];
            for (int x=0; x<width; x++) {
                png_byte* ptr = &(row[x*numberOfChannels*byteshift]);

                ycbcr.val[0] = img->val[p];
                ycbcr.val[1] = img->Cb[p];
                ycbcr.val[2] = img->Cr[p];
                ushort alpha = img->alpha[p];

                rgb = iftYCbCrBT2020toRGB(ycbcr, depth, depth);

                ptr[0*byteshift] = rgb.val[0] & 0xFF;//get first byte
                ptr[1*byteshift] = rgb.val[1] & 0xFF;
                ptr[2*byteshift] = rgb.val[2] & 0xFF;
                ptr[3*byteshift] = alpha & 0xFF;

                if(depth==16) {//in 16bit image, we should store as big endian
                    ptr[(0*byteshift)+1] = ptr[0*byteshift];
                    ptr[(1*byteshift)+1] = ptr[1*byteshift];
                    ptr[(2*byteshift)+1] = ptr[2*byteshift];
                    ptr[(3*byteshift)+1] = ptr[(3*byteshift)];

                    ptr[0*byteshift] = ((rgb.val[0]>>8) & 0xFF);//get second byte
                    ptr[1*byteshift] = ((rgb.val[1]>>8) & 0xFF);
                    ptr[2*byteshift] = ((rgb.val[2]>>8) & 0xFF);
                    ptr[(3*byteshift)] = ((alpha>>8) & 0xFF);
                }
                p++;
            }
        }

    }else{
        iftError("Unknwon color scape", "iftWriteImagePNG");
    };


    iftWritePngImageAux(filename, row_pointers, width, height, depth, color_type);
}

void iftWriteImageJPEG(const iftImage* img, const char* format, ...)
{
    va_list args;
    char filename[IFT_STR_DEFAULT_SIZE];

    va_start(args, format);
    vsprintf(filename, format, args);
    va_end(args);

    //code based on external/libjpeg/source/example.c
    /* This struct contains the JPEG compression parameters and pointers to
 * working space (which is allocated as needed by the JPEG library).
 * It is possible to have several such structures, representing multiple
 * compression/decompression processes, in existence at once.  We refer
 * to any one struct (and its associated working data) as a "JPEG object".
 */
    struct jpeg_compress_struct cinfo;

    /* This struct represents a JPEG error handler.  It is declared separately
 * because applications often want to supply a specialized error handler
 * (see the second half of this file for an example).  But here we just
 * take the easy way out and use the standard error handler, which will
 * print a message on stderr and call exit() if compression fails.
 * Note that this struct must live as long as the main JPEG parameter
 * struct, to avoid dangling-pointer problems.
 */
    struct jpeg_error_mgr jerr;

    /* More stuff */
    FILE * outfile;		/* target file */
    JSAMPARRAY buffer;
    /* Step 1: allocate and initialize JPEG compression object */

    /* We have to set up the error handler first, in case the initialization
     * step fails.  (Unlikely, but it could happen if you are out of memory.)
     * This routine fills in the contents of struct jerr, and returns jerr's
     * address which we place into the link field in cinfo.
     */
    cinfo.err = jpeg_std_error(&jerr);

    /* Now we can initialize the JPEG compression object. */
    jpeg_create_compress(&cinfo);
    /* Step 2: specify data destination (eg, a file) */
    /* Note: steps 2 and 3 can be done in either order. */

    /* Here we use the library-supplied code to send compressed data to a
     * stdio stream.  You can also write your own code to do something else.
     * VERY IMPORTANT: use "b" option to fopen() if you are on a machine that
     * requires it in order to write binary files.
     */
    if ((outfile = fopen(filename, "wb")) == NULL) {
        fprintf(stderr, "can't open %s\n", filename);
        exit(1);
    }
    jpeg_stdio_dest(&cinfo, outfile);

    /* First we supply a description of the input image.
* Four fields of the cinfo struct must be filled in:
*/


    cinfo.image_width = img->xsize; 	/* image width and height, in pixels */
    cinfo.image_height = img->ysize;
    cinfo.data_precision = iftImageDepth(img);

    cinfo.input_components = 3;
    cinfo.in_color_space = JCS_YCbCr;
    cinfo.jpeg_color_space = JCS_YCbCr;



    /* Now use the library's routine to set default compression parameters.
* (You must set at least cinfo.in_color_space before calling this,
* since the defaults depend on the source color space.)
*/
    jpeg_set_defaults(&cinfo);

    /* Now you can set any non-default parameters you wish to.
* Here we just illustrate the use of quality (quantization table) scaling:
*/
    int quality = 100;
    jpeg_set_quality(&cinfo, quality, TRUE /* limit to baseline-JPEG values */);

    /* Step 4: Start compressor */
    /* TRUE ensures that we will write a complete interchange-JPEG file.
     * Pass TRUE unless you are very sure of what you're doing.
     */
    jpeg_start_compress(&cinfo, TRUE);

    /* Step 5: while (scan lines remain to be written) */
    /*           jpeg_write_scanlines(...); */

    /* Here we use the library's state variable cinfo.next_scanline as the
     * loop counter, so that we don't have to keep track ourselves.
     * To keep things simple, we pass one scanline per call; you can pass
     * more if you wish, though.
     */
    int row_stride = cinfo.image_width * cinfo.num_components;
    buffer = (*cinfo.mem->alloc_sarray)
            ((j_common_ptr) &cinfo, JPOOL_IMAGE, row_stride, 1);
    unsigned int imageRow = 0;
    while (cinfo.next_scanline < cinfo.image_height) {
        /* jpeg_write_scanlines expects an array of pointers to scanlines.
         * Here the array is only one element long, but you could pass
         * more than one scanline at a time if that's more convenient.
         */
        unsigned int imageCol = 0;
        for (unsigned int i = 0; i < (unsigned int)cinfo.image_width; i++) {
            int shift = i*cinfo.num_components;
            buffer[0][(shift+0)] = (unsigned char)iftImgVal(img,imageCol,imageRow,0);
            buffer[0][(shift+1)] = (unsigned char)iftImgCb(img,imageCol,imageRow,0);
            buffer[0][(shift+2)] = (unsigned char)iftImgCr(img,imageCol,imageRow,0);

            imageCol++;
        }
        imageRow++;
        (void) jpeg_write_scanlines(&cinfo, buffer, 1);
    }

    /* Step 6: Finish compression */

    jpeg_finish_compress(&cinfo);
    /* After finish_compress, we can close the output file. */
    fclose(outfile);

    /* Step 7: release JPEG compression object */

    /* This is an important step since it will release a good deal of memory. */
    jpeg_destroy_compress(&cinfo);
}

int iftMaximumValueInRegion(const iftImage *img, iftBoundingBox bb) 
{
    // checkers
    if (img == NULL)
        iftError("Image is NULL", "iftMaximumValueInRegion");
    if (!iftValidVoxel(img, bb.begin))
        iftError("Beginning voxel (%d, %d, %d) from Region (Bound. Box) is not in the Image Domain\n" \
                 "Img (xsize, ysize, zsize): (%d, %d, %d)", "iftMaximumValueInRegion",
                 bb.begin.x, bb.begin.y, bb.begin.z, img->xsize, img->ysize, img->zsize);
    if (!iftValidVoxel(img, bb.end))
        iftError("Ending voxel (%d, %d, %d) from Region (Bound. Box) is not in the Image Domain\n" \
                 "Img (xsize, ysize, zsize): (%d, %d, %d)", "iftMaximumValueInRegion",
                 bb.end.x, bb.end.y, bb.end.z, img->xsize, img->ysize, img->zsize);

    int img_max_val = IFT_INFINITY_INT_NEG;

    iftVoxel v;
    for (v.z = bb.begin.z; v.z <= bb.end.z; v.z++) {
        for (v.y = bb.begin.y; v.y <= bb.end.y; v.y++) {
            for (v.x = bb.begin.x; v.x <= bb.end.x; v.x++) {
                int p = iftGetVoxelIndex(img, v);
                if (img_max_val < img->val[p]) {
                    img_max_val = img->val[p];
                }
            }
        }
    }

    return img_max_val;
}

void iftSetImage(iftImage *img, int value) 
{
    for (int p = 0; p < img->n; p++)
        img->val[p] = value;
}

void  iftSetAlpha(iftImage *img, ushort value)
{
    int p;
    if(img->alpha == NULL){
        img->alpha = iftAllocUShortArray(img->n);
    }

    for (p=0; p < img->n; p++) {
        img->alpha[p] = value;
    }
}

void    iftSetCbCr(iftImage *img, ushort value)
{
    int p;

    if (!iftIsColorImage(img)){
        img->Cb = iftAllocUShortArray(img->n);
        img->Cr = iftAllocUShortArray(img->n);
    }
    for (p=0; p < img->n; p++) {
        img->Cb[p] = value;
        img->Cr[p] = value;
    }
}

void iftVerifyImageDomains(const iftImage *img1, const iftImage *img2, const char *function) 
{
    if ((img1==NULL)||(img2==NULL)||(img1->xsize!=img2->xsize) || (img1->ysize!=img2->ysize) || (img1->zsize!=img2->zsize)) {
        iftError("Images with different domains:\n" \
                "img1 (xsize, ysize, zsize): (%d, %d, %d)\n" \
                "img2 (xsize, ysize, zsize): (%d, %d, %d)\n",
                 function,
                 img1->xsize, img1->ysize, img1->zsize, img2->xsize, img2->ysize, img2->zsize);
    }
}

uchar iftImageDepth(const iftImage *img) 
{
    int img_min, img_max;
    iftMinMaxValues(img, &img_min, &img_max);
    
    long max_range;

    if (img_min >= 0)
        max_range = iftNormalizationValue(img_max) + 1;
    else
        max_range = iftNormalizationValue(img_max - img_min) + 1;
    
    return (uchar) iftLog(max_range, 2);
}

void iftMinMaxValues(const iftImage *img, int *min, int *max) 
{
    *min = *max = img->val[0];

    for (int p = 1; p < img->n; p++) {
        if (img->val[p] < *min)
            *min = img->val[p];
        else if (img->val[p] > *max)
            *max = img->val[p];
    }
}

int iftNumberOfElements(const iftImage *mask)
{
    int p, nnodes = 0;

    for (p=0; p < mask->n; p++)
        if (mask->val[p]>0){
            nnodes++;
        }
    return(nnodes);
}

iftImage *iftImageGradientMagnitude(const iftImage *img, iftAdjRel *Ain)
{
    iftAdjRel *A = NULL;
    if (Ain == NULL) {
        if (iftIs3DImage(img))
            A = iftSpheric(1.0);
        else A = iftCircular(1.0);
    }
    else A = Ain;

    float   dist,gx,gy,gz, g, gmax;
    float   gxCb , gyCb , gzCb, gxCr , gyCr , gzCr;
    iftVoxel   u,v;
    int i;
    float     *mag =iftAllocFloatArray(A->n), *weight = iftAllocFloatArray(A->n);
    float     _2sigma2;
    int        dx, dy, dz;
    iftImage  *grad=iftCreateImage(img->xsize,img->ysize,img->zsize);

    iftCopyVoxelSize(img,grad);

    iftMaxAdjShifts(A, &dx, &dy, &dz);
    _2sigma2 = 2.0*(dx*dx+dy*dy+dz*dz)/9.0;
    for (i=0; i < A->n; i++){
        mag[i]=sqrtf(A->dx[i]*A->dx[i]+A->dy[i]*A->dy[i]+A->dz[i]*A->dz[i]);
        weight[i]=exp(-mag[i]/_2sigma2)/mag[i];
    }

    if ( (img->Cb == NULL) && (img->Cr == NULL) ){
        for (int z=0; z < img->zsize; z++)
            for (int y=0; y < img->ysize; y++)
#pragma omp parallel for shared(z,y,A,img,weight,grad) private(u,v,gx,gy,gz,i,dist)
                    for (int x=0; x < img->xsize; x++) {
                        u.x=x; u.y=y; u.z=z;
                        int p = iftGetVoxelIndex(img,u);
                        gx = gy = gz = 0.0;
                        for (i=1; i < A->n; i++) {
                            v.x = u.x + A->dx[i];
                            v.y = u.y + A->dy[i];
                            v.z = u.z + A->dz[i];
                            if (iftValidVoxel(img,v)){
                                int q = iftGetVoxelIndex(img,v);
                                dist = img->val[q]-img->val[p];
                                gx  += dist*A->dx[i]*weight[i];
                                gy  += dist*A->dy[i]*weight[i];
                                gz  += dist*A->dz[i]*weight[i];
                            }
                        }
                        grad->val[p]=(int)sqrtf(gx*gx + gy*gy + gz*gz);
                    }
    }else{ // colored image
        for (int z=0; z < img->zsize; z++)
            for (int y=0; y < img->ysize; y++)
#pragma omp parallel for shared(z,y,A,img,weight,grad) private(u,v,gx,gy,gz,i,dist,gxCb,gyCb,gzCb,gxCr,gyCr,gzCr,gmax,g)
                    for (int x=0; x < img->xsize; x++) {
                        u.x=x; u.y=y; u.z=z;
                        int p = iftGetVoxelIndex(img,u);
                        gx = gy = gz = 0.0;
                        gxCb = gyCb = gzCb = 0.0;
                        gxCr = gyCr = gzCr = 0.0;
                        for (i=1; i < A->n; i++) {
                            v.x = u.x + A->dx[i];
                            v.y = u.y + A->dy[i];
                            v.z = u.z + A->dz[i];
                            if (iftValidVoxel(img,v)){
                                int q = iftGetVoxelIndex(img,v);
                                dist = img->val[q]-img->val[p];
                                gx  += dist*A->dx[i]*weight[i];
                                gy  += dist*A->dy[i]*weight[i];
                                gz  += dist*A->dz[i]*weight[i];
                                dist = img->Cb[q]-img->Cb[p];
                                gxCb  += dist*A->dx[i]*weight[i];
                                gyCb  += dist*A->dy[i]*weight[i];
                                gzCb  += dist*A->dz[i]*weight[i];
                                dist = img->Cr[q]-img->Cr[p];
                                gxCr  += dist*A->dx[i]*weight[i];
                                gyCr  += dist*A->dy[i]*weight[i];
                                gzCr  += dist*A->dz[i]*weight[i];
                            }
                        }
                        gmax = sqrtf(gx*gx + gy*gy + gz*gz);
                        g    = sqrtf(gxCb*gxCb + gyCb*gyCb + gzCb*gzCb);
                        if (g > gmax)
                            gmax = g;
                        g    = sqrtf(gxCr*gxCr + gyCr*gyCr + gzCr*gzCr);
                        if (g > gmax)
                            gmax = g;
                        grad->val[p] = (int)gmax;
                    }
    }

    iftFree(mag);
    iftFree(weight);

    if (Ain == NULL) {
        iftDestroyAdjRel(&A);
    }

    return(grad);
}

iftImage *iftCreateImageFromImage(const iftImage *src) 
{
    iftImage *out = NULL;

    if (src != NULL) {
        if (iftIsColorImage(src)) {
            out = iftCreateColorImage(src->xsize, src->ysize, src->zsize, iftImageDepth(src));
        }
        else {
            out = iftCreateImage(src->xsize, src->ysize, src->zsize);
        }
        iftCopyVoxelSize(src, out);
    }

    return out;
}

iftImage  *iftCreateColorImage(int xsize,int ysize,int zsize, int depth)
{
    iftImage *img=NULL;
    img = iftCreateImage(xsize, ysize, zsize);

    iftSetCbCr(img, (iftMaxImageRange(depth)+1)/2);

    return(img);
}

// ---------- iftImage.c end 
// ---------- iftFImage.c start 

iftFImage *iftCreateFImage(int xsize, int ysize, int zsize) 
{
    iftFImage *img = NULL;
    int       y, z, xysize;


    img = (iftFImage *) iftAlloc(1, sizeof(iftFImage));
    if (img == NULL) {
        iftError(MSG_MEMORY_ALLOC_ERROR, "iftCreateFImage");
    }

    img->val   = iftAllocFloatArray(xsize * ysize * zsize);
    img->xsize = xsize;
    img->ysize = ysize;
    img->zsize = zsize;
    img->dx    = 1.0;
    img->dy    = 1.0;
    img->dz    = 1.0;
    img->tby   = iftAllocIntArray(ysize);
    img->tbz   = iftAllocIntArray(zsize);
    img->n     = xsize * ysize * zsize;

    if (img->val == NULL || img->tbz == NULL || img->tby == NULL) {
        iftError(MSG_MEMORY_ALLOC_ERROR, "iftCreateFImage");
    }

    img->tby[0] = 0;
    for (y = 1; y < ysize; y++)
        img->tby[y] = img->tby[y - 1] + xsize;

    img->tbz[0] = 0;
    xysize = xsize * ysize;
    for (z = 1; z < zsize; z++)
        img->tbz[z] = img->tbz[z - 1] + xysize;

    return (img);
}

void iftDestroyFImage(iftFImage **img) 
{
    iftFImage *aux;

    aux = *img;
    if (aux != NULL) {
        if (aux->val != NULL) iftFree(aux->val);
        if (aux->tby != NULL) iftFree(aux->tby);
        if (aux->tbz != NULL) iftFree(aux->tbz);
        iftFree(aux);
        *img = NULL;
    }
}

iftVoxel iftFGetVoxelCoord(const iftFImage *img, int p) 
{
    iftVoxel u;

    u.x = iftFGetXCoord(img, p);
    u.y = iftFGetYCoord(img, p);
    u.z = iftFGetZCoord(img, p);

    return (u);
}

char iftFValidVoxel(const iftFImage *img, iftVoxel v) 
{
    if ((v.x >= 0) && (v.x < img->xsize) &&
        (v.y >= 0) && (v.y < img->ysize) &&
        (v.z >= 0) && (v.z < img->zsize))
        return (1);
    else
        return (0);
}

// ---------- iftFImage.c end 
// ---------- iftGraphics.c start

void iftDrawPoint(iftImage *img, iftVoxel u, iftColor YCbCr, iftAdjRel *B, int normValue)
{
    int q,i;
    iftVoxel v;

    if (!iftValidVoxel(img,u))
        iftError("Point is outside the image domain", "iftDrawPoint");

    if (!iftIsColorImage(img)){
        img->Cb = iftAllocUShortArray(img->n);
        img->Cr = iftAllocUShortArray(img->n);
        for (q=0; q < img->n; q++) {
            img->Cb[q] = normValue/2;
            img->Cr[q] = normValue/2;
        }
    }

    for (i=0; i < B->n; i++) {
        v.x = u.x + B->dx[i];
        v.y = u.y + B->dy[i];
        v.z = u.z + B->dz[i];
        if (iftValidVoxel(img,v)){
            q = iftGetVoxelIndex(img,v);
            img->val[q]=YCbCr.val[0];
            img->Cb[q]=(ushort) YCbCr.val[1];
            img->Cr[q]=(ushort) YCbCr.val[2];
        }
    }

}

void iftDrawBorders(iftImage *img, iftImage *label, iftAdjRel *A, iftColor YCbCr, iftAdjRel *B)
{
    iftVoxel u,v;
    int i,p,q;
    int maxRangeValue = iftNormalizationValue(iftMaximumValue(img));

    if ((img->xsize != label->xsize)||
        (img->ysize != label->ysize)||
        (img->zsize != label->zsize))
        iftError("Images must have the same domain", "iftDrawBorders");


    if (!iftIsColorImage(img))
        iftSetCbCr(img,maxRangeValue/2);

    if (A->n > 1){
        for (p=0; p < img->n; p++) {
            u = iftGetVoxelCoord(label,p);
            for (i=0; i < A->n; i++) {
                v = iftGetAdjacentVoxel(A,u,i);
                if (iftValidVoxel(label,v)){
                    q = iftGetVoxelIndex(label,v);
                    if (label->val[p] != label->val[q]){
                        iftDrawPoint(img, u, YCbCr, B, maxRangeValue);
                        break;
                    }
                }
            }
        }
    } else {
        for (p=0; p < img->n; p++)
            if (label->val[p] != 0){
                u = iftGetVoxelCoord(label,p);
                iftDrawPoint(img, u, YCbCr, B, maxRangeValue);
            }
    }
}

void iftDrawBordersSingleLabel(iftImage *img, iftImage *labelMap, const int label, iftColor YCbCr)
{
  assert(img != NULL);
  assert(labelMap != NULL);
  assert(label > 0);
  assert(label <= iftMaximumValue(labelMap));
  assert(img->xsize == labelMap->xsize);
  assert(img->ysize == labelMap->ysize);
  assert(img->zsize == labelMap->zsize);

  int maxRangeValue = iftNormalizationValue(iftMaximumValue(img));
  iftAdjRel *A = iftCircular(1.0);
  iftAdjRel *brushShape = iftCircular(1.0);
  
  for (int p = 0; p < img->n; ++p) {
    if (labelMap->val[p] != label)
      continue;

    iftVoxel u = iftGetVoxelCoord(labelMap, p);

    for (int i = 1; i < A->n; ++i) {
      iftVoxel v = iftGetAdjacentVoxel(A, u, i);
      if (iftValidVoxel(labelMap, v)) {
        int q = iftGetVoxelIndex(labelMap, v);

        if (label != labelMap->val[q]) {
          iftDrawPoint(img, u, YCbCr, brushShape, maxRangeValue);
          break;
        }
      }
    }
  }
}

// ---------- iftGraphics.c end
// ---------- iftMatrix.c start 

iftMatrix *iftCreateMatrix(int ncols, int nrows) 
{
    iftMatrix *M = (iftMatrix *) iftAlloc(1, sizeof(iftMatrix));
    
    M->ncols = ncols;
    M->nrows = nrows;
    M->tbrow = (long*)iftAlloc(nrows,sizeof(long));
    //M->tbrow = iftAllocIntArray(nrows);
    for (long r = 0; r < (long)nrows; r++) {
        M->tbrow[r] = r * ncols;
    }
    M->n = (long) ncols * (long) nrows;
    M->allocated = true;
    
    M->val = iftAllocFloatArray(M->n);
    
    return (M);
}

iftMatrix *iftCopyMatrix(const iftMatrix *A) 
{
    iftMatrix *B;
    int i;

    B = iftCreateMatrix(A->ncols, A->nrows);

    for (i = 0; i < A->n; i++)
        B->val[i] = A->val[i];

    return (B);
}

void iftDestroyMatrix(iftMatrix **M) 
{
    if (M != NULL) {
        iftMatrix *aux = *M;
        
        if (aux != NULL) {
            if (aux->allocated && aux->val != NULL)
                iftFree(aux->val);
            iftFree(aux->tbrow);
            iftFree(aux);
        }
        *M = NULL;
    }
}

iftStrMatrix *iftCopyStrMatrix(char **val, int nrows, int ncols) 
{
    iftStrMatrix *mat = iftCreateStrMatrix(ncols, nrows);

    for (int i = 0; i < mat->n; i++)
        mat->val[i] = iftCopyString(val[i]);

    return mat;
}

void iftDestroyStrMatrix(iftStrMatrix **sM) 
{
    iftStrMatrix *aux = *sM;

    if (aux != NULL) {
        iftFree(aux->tbrow);
        if (aux->val != NULL) {
            for (int i = 0; i < aux->n; i++)
                iftFree(aux->val[i]);
            iftFree(aux->val);
        }
        iftFree(aux);
        *sM = NULL;
    }
}

iftStrMatrix *iftCreateStrMatrix(int ncols, int nrows) 
{
    iftStrMatrix *sM = (iftStrMatrix *) iftAlloc(1, sizeof(iftStrMatrix));

    sM->ncols = ncols;
    sM->nrows = nrows;

    sM->tbrow = iftAllocIntArray(nrows);
    for (int r = 0; r < nrows; r++) {
        sM->tbrow[r] = r * ncols;
    }

    sM->n = ncols * nrows;
    sM->val = (char **) iftAlloc(sM->n, sizeof(char *));
    for (int i = 0; i < sM->n; i++)
        sM->val[i] = iftAllocCharArray(1024);

    return sM;
}

// ---------- iftMatrix.c end
// ---------- iftMImage.c start 

iftMImage * iftCreateMImage(int xsize,int ysize,int zsize, int nbands)
{
  iftMImage *img=NULL;
  int        i,y,z,xysize;

  img = (iftMImage *) iftAlloc(1,sizeof(iftMImage));
  if (img == NULL){
      iftError(MSG_MEMORY_ALLOC_ERROR, "iftCreateMImage");
  }

  img->n       = xsize*ysize*zsize;
  img->m       = nbands;
  img->data    = iftCreateMatrix(img->m, img->n);
  img->val     = iftAlloc(img->n, sizeof *img->val);
  for (i = 0; i < img->n; i++)
      img->val[i] = iftMatrixRowPointer(img->data, i);
  img->xsize   = xsize;
  img->ysize   = ysize;
  img->zsize   = zsize;
  img->dx      = 1.0;
  img->dy      = 1.0;
  img->dz      = 1.0;
  img->tby     = iftAllocIntArray(ysize);
  img->tbz     = iftAllocIntArray(zsize);

  img->tby[0]=0;
  for (y=1; y < ysize; y++)
    img->tby[y]=img->tby[y-1] + xsize;

  img->tbz[0]=0; xysize = xsize*ysize;
  for (z=1; z < zsize; z++)
    img->tbz[z]=img->tbz[z-1] + xysize;

  return(img);
}

void iftDestroyMImage(iftMImage **img)
{
    iftMImage *aux;

    aux = *img;
    if (aux != NULL) {
        if (aux->val != NULL)
            iftFree(aux->val);
        if (aux->data != NULL)
            iftDestroyMatrix(&aux->data);
        if (aux->tby != NULL)
            iftFree(aux->tby);
        if (aux->tbz != NULL)
            iftFree(aux->tbz);
        iftFree(aux);
        *img = NULL;
    }
}

iftMImage * iftImageToMImage(const iftImage *img1, char color_space)
{
  iftMImage *img2=NULL;
  int normalization_value = iftNormalizationValue(iftMaximumValue(img1));

  switch (color_space) {

  case YCbCr_CSPACE:
    img2=iftCreateMImage(img1->xsize,img1->ysize,img1->zsize,3);
#pragma omp parallel for shared(img1, img2, normalization_value)
    for (int p=0; p < img2->n; p++) {
      img2->val[p][0]=((float)img1->val[p]);
      img2->val[p][1]=((float)img1->Cb[p]);
      img2->val[p][2]=((float)img1->Cr[p]);
    }
    break;

  case YCbCrNorm_CSPACE:

    img2=iftCreateMImage(img1->xsize,img1->ysize,img1->zsize,3);
#pragma omp parallel for shared(img1, img2, normalization_value)
    for (int p=0; p < img2->n; p++) {
      img2->val[p][0]=((float)img1->val[p])/(float)normalization_value;
      img2->val[p][1]=((float)img1->Cb[p])/(float)normalization_value;
      img2->val[p][2]=((float)img1->Cr[p])/(float)normalization_value;
    }
    break;

  case LAB_CSPACE:

    img2=iftCreateMImage(img1->xsize,img1->ysize,img1->zsize,3);
#pragma omp parallel for shared(img1, img2, normalization_value)
    for (int p=0; p < img2->n; p++) {
      iftColor  YCbCr,RGB;
      iftFColor Lab;
      YCbCr.val[0] = img1->val[p];
      YCbCr.val[1] = img1->Cb[p];
      YCbCr.val[2] = img1->Cr[p];
      RGB = iftYCbCrtoRGB(YCbCr,normalization_value);
      Lab = iftRGBtoLab(RGB,normalization_value);
      img2->val[p][0]=Lab.val[0];
      img2->val[p][1]=Lab.val[1];
      img2->val[p][2]=Lab.val[2];
    }
    break;

    case LABNorm_CSPACE:

      img2=iftCreateMImage(img1->xsize,img1->ysize,img1->zsize,3);
#pragma omp parallel for shared(img1, img2, normalization_value)
      for (int p=0; p < img2->n; p++) {
        iftColor  YCbCr,RGB;
        iftFColor Lab;
        YCbCr.val[0] = img1->val[p];
        YCbCr.val[1] = img1->Cb[p];
        YCbCr.val[2] = img1->Cr[p];
        RGB = iftYCbCrtoRGB(YCbCr,normalization_value);
        Lab = iftRGBtoLabNorm(RGB,normalization_value);
        img2->val[p][0]=Lab.val[0];
        img2->val[p][1]=Lab.val[1];
        img2->val[p][2]=Lab.val[2];
      }
      break;

  case LABNorm2_CSPACE:

    img2=iftCreateMImage(img1->xsize,img1->ysize,img1->zsize,3);
#pragma omp parallel for shared(img1, img2, normalization_value)
    for (int p=0; p < img2->n; p++) {
      iftColor  YCbCr,RGB;
      iftFColor LabNorm;
      YCbCr.val[0] = img1->val[p];
      YCbCr.val[1] = img1->Cb[p];
      YCbCr.val[2] = img1->Cr[p];
      RGB = iftYCbCrtoRGB(YCbCr,normalization_value);
      LabNorm = iftRGBtoLabNorm2(RGB,normalization_value);
      img2->val[p][0]=LabNorm.val[0];
      img2->val[p][1]=LabNorm.val[1];
      img2->val[p][2]=LabNorm.val[2];
    }
    break;

  case RGB_CSPACE:

    img2=iftCreateMImage(img1->xsize,img1->ysize,img1->zsize,3);
#pragma omp parallel for shared(img1, img2, normalization_value)
    for (int p=0; p < img2->n; p++) {
      iftColor YCbCr,RGB;
      YCbCr.val[0] = img1->val[p];
      YCbCr.val[1] = img1->Cb[p];
      YCbCr.val[2] = img1->Cr[p];
      RGB          = iftYCbCrtoRGB(YCbCr,normalization_value);
      img2->val[p][0]=((float)RGB.val[0]);
      img2->val[p][1]=((float)RGB.val[1]);
      img2->val[p][2]=((float)RGB.val[2]);
    }
    break;

  case RGBNorm_CSPACE:

    img2=iftCreateMImage(img1->xsize,img1->ysize,img1->zsize,3);
#pragma omp parallel for shared(img1, img2, normalization_value)
    for (int p=0; p < img2->n; p++) {
      iftColor YCbCr,RGB;
      YCbCr.val[0] = img1->val[p];
      YCbCr.val[1] = img1->Cb[p];
      YCbCr.val[2] = img1->Cr[p];
      RGB          = iftYCbCrtoRGB(YCbCr,normalization_value);
      img2->val[p][0]=((float)RGB.val[0])/(float)normalization_value;
      img2->val[p][1]=((float)RGB.val[1])/(float)normalization_value;
      img2->val[p][2]=((float)RGB.val[2])/(float)normalization_value;
    }
    break;

  case GRAY_CSPACE:

    img2=iftCreateMImage(img1->xsize,img1->ysize,img1->zsize,1);
#pragma omp parallel for shared(img1, img2, normalization_value)
    for (int p=0; p < img2->n; p++) {
      img2->val[p][0]=((float)img1->val[p]);
    }
    break;

  case GRAYNorm_CSPACE:

    img2=iftCreateMImage(img1->xsize,img1->ysize,img1->zsize,1);
#pragma omp parallel for shared(img1, img2, normalization_value)
    for (int p=0; p < img2->n; p++) {
      img2->val[p][0]=((float)img1->val[p])/(float)normalization_value;
    }

    break;

  case WEIGHTED_YCbCr_CSPACE:

    img2=iftCreateMImage(img1->xsize,img1->ysize,img1->zsize,3);

#pragma omp parallel for shared(img1, img2, normalization_value)
    for (int p=0; p < img2->n; p++) {
      img2->val[p][0]=(0.2/2.2)*((float)img1->val[p]/(float)normalization_value);
      img2->val[p][1]=(1.0/2.2)*((float)img1->Cb[p]/(float)normalization_value);
      img2->val[p][2]=(1.0/2.2)*((float)img1->Cr[p]/(float)normalization_value);
    }
    break;
  
  case HSV_CSPACE:

    img2=iftCreateMImage(img1->xsize,img1->ysize,img1->zsize,3);
#pragma omp parallel for shared(img1, img2, normalization_value)
    for (int p=0; p < img2->n; p++) {
      iftColor  YCbCr,RGB;
      iftColor HSV;
      YCbCr.val[0] = img1->val[p];
      YCbCr.val[1] = img1->Cb[p];
      YCbCr.val[2] = img1->Cr[p];
      RGB = iftYCbCrtoRGB(YCbCr,normalization_value);
      HSV = iftRGBtoHSV(RGB,normalization_value);
      img2->val[p][0]=HSV.val[0];
      img2->val[p][1]=HSV.val[1];
      img2->val[p][2]=HSV.val[2];
    }
    break;

  default:
      iftError("Invalid color space (see options in iftColor.h)", "iftImageToMImage");
  }

  img2->dx = img1->dx;
  img2->dy = img1->dy;
  img2->dz = img1->dz;

  return(img2);
}

iftImage * iftMImageToImage(const iftMImage *img1, int Imax, int band)
{
  iftImage *img2=iftCreateImage(img1->xsize,img1->ysize,img1->zsize);
  int p,b=band;
  double min = IFT_INFINITY_FLT, max = IFT_INFINITY_FLT_NEG;

  if ((band < 0)||(band >= img1->m))
      iftError("Invalid band", "iftMImageToImage");

  for (p=0; p < img1->n; p++) {
    if (img1->val[p][b] < min)
      min = img1->val[p][b];
    if (img1->val[p][b] > max)
      max = img1->val[p][b];
  }

  //printf("min %lf max %lf\n",min,max);

  if (max > min){
    for (p=0; p < img2->n; p++) {
      img2->val[p]=(int)(Imax*(img1->val[p][b]-min)/(max-min));
    }
  }else{
    char msg[100];
    sprintf(msg,"Image is empty: max = %f and min = %f\n",max,min);
    // iftWarning(msg,"iftMImageToImage");
  }

  img2->dx = img1->dx;
  img2->dy = img1->dy;
  img2->dz = img1->dz;

  return(img2);
}

iftImage *iftGridSampling(iftMImage *img, iftImage *mask1, int nsamples)
{
  iftImage  *mask2 = iftCreateImage(img->xsize,img->ysize,img->zsize);
  iftImage  *prob  = iftBorderProbImage(img);
  int        p, q, i, qmin, xsize,ysize,zsize;
  float      xspacing,yspacing,zspacing,deltax,deltay,deltaz;
  iftAdjRel *A;
  iftVoxel   u, v, m, uo, uf;

  /* Compute the extreme voxels that define the region of interest in
     mask1, and then compute the xsize, ysize, and zsize of that
     ROI. */
  iftBoundingBox mbb = iftMinBoundingBox(mask1, NULL);
  uo = mbb.begin;
  uf = mbb.end;
  xsize = uf.x - uo.x + 1;
  ysize = uf.y - uo.y + 1;
  zsize = uf.z - uo.z + 1;

  if (iftIs3DMImage(img)){

    A  = iftSpheric(sqrtf(3.0));
    /* finds displacements along each axis */
    /* uncomment the next 4 lines to use same number o superpixels per axis */
    /*
    float nsamples_per_axis = (float) pow((double)nsamples,1.0/3.0);
    xspacing          = (xsize/nsamples_per_axis);
    yspacing          = (ysize/nsamples_per_axis);
    zspacing          = (zsize/nsamples_per_axis);
    */
    /* uncomment the next 5 lines to obtain equally spaced seeds in every axis */
    float superpixelsize = 0.5+(float)(xsize*ysize*zsize)/(float)nsamples;
    float step = (float) pow((double)superpixelsize,1.0/3.0)+0.5;
    xspacing = step;
    yspacing = step;
    zspacing = step;

    deltax            = xspacing/2.0;
    deltay            = yspacing/2.0;
    deltaz            = zspacing/2.0;

    if ((deltax < 1.0)||(deltay < 1.0)||(deltaz < 1.0))
        iftError("The number of samples is too high", "iftGridSampling");

    for (m.z=(int)deltaz; m.z < (zsize-deltaz); m.z = (int)(m.z + zspacing)) {
      for (m.y=(int)deltay; m.y < (ysize-deltay); m.y = (int)(m.y + yspacing)) {
  for (m.x=(int)deltax; m.x < (xsize-deltax); m.x = (int)(m.x + xspacing)) {
    u.z = uo.z + m.z; u.y = uo.y + m.y; u.x = uo.x + m.x;
    p = iftGetVoxelIndex(mask1,u);
    if (mask1->val[p]!=0){
      for (i=1, qmin=p; i < A->n; i++) {
        v = iftGetAdjacentVoxel(A,u,i);
        q = iftGetVoxelIndex(mask1,v);
        if (iftValidVoxel(mask1,v) && (prob->val[q]<prob->val[qmin])&&
    (mask2->val[q]==0)&&
      (mask1->val[q]!=0))
    qmin=q;
      }
      mask2->val[qmin]=1;
    }
  }
      }
    }
  }else{
    A   = iftCircular(sqrtf(2.0));

    /* finds displacements along each axis  */
    /* uncomment the next 3 lines to use same number o superpixels per axis */
//    float nsamples_per_axis = (float) sqrt((double)nsamples);
//    xspacing          = (xsize/nsamples_per_axis);
//    yspacing          = (ysize/nsamples_per_axis);
    /* uncomment the next 4 lines to obtain equally spaced seeds in every axis */
    float superpixelsize = 0.5+(float)(xsize*ysize)/(float)nsamples;
    float step = sqrt(superpixelsize)+0.5;
    xspacing = step;
    yspacing = step;


    deltax            = xspacing/2.0;
    deltay            = yspacing/2.0;

    if ((deltax < 1.0)||(deltay < 1.0))
        iftError("The number of samples is too high", "iftGridSampling");

    u.z = 0;
    for (m.y=(int)deltay; m.y < ysize; m.y = (int)(m.y + yspacing)){
      for (m.x=(int)deltax; m.x < xsize; m.x = (int)(m.x + xspacing)) {
  u.y = uo.y + m.y; u.x = uo.x + m.x;
  p = iftGetVoxelIndex(mask1,u);
  if (mask1->val[p]!=0){
    for (i=1, qmin=p; i < A->n; i++) {
      v = iftGetAdjacentVoxel(A,u,i);
      if (iftValidVoxel(mask1,v)){
        q = iftGetVoxelIndex(mask1,v);
        if ((prob->val[q]<prob->val[qmin])&&
            (mask2->val[q]==0)&&
            (mask1->val[q]!=0))
          qmin=q;
      }
    }
    mask2->val[qmin]=1;
  }
      }
    }
  }

  //printf("nseeds: %d\n",iftNumberOfElements(mask2));

  iftDestroyImage(&prob);
  iftDestroyAdjRel(&A);

  return(mask2);
}

float iftNormalizedShannonEntropy(float *arr, int size ) 
{
  int i, nbins, b;
  int *histogram;
  float im, entropy, binsize;
  float minVal, maxVal, range, factor;
  minVal = IFT_INFINITY_FLT;
  maxVal = IFT_INFINITY_FLT_NEG;
  binsize = 5.0;
  entropy = 0.0;
  factor = 0.0;
  // Quantize
  for (i=0; i < size; i++) {
    if (arr[i] < minVal)
      minVal = arr[i];
    if (arr[i] > maxVal)
      maxVal = arr[i];
  }
  range = maxVal - minVal;
  nbins = iftMax((int)(range / binsize), 1);
  histogram = iftAllocIntArray(nbins);

  if (range > 0)
    factor = ((float)nbins)/range;

  for (i=0; i < size; i++) {
    b = (int)( ( arr[i]- minVal )*factor );
    b = iftMin(b, nbins - 1);
    histogram[b]++;
  }

  // Compute entropy
  im = (float)size;
  for (i=0; i < nbins; i++) {
    if (histogram[i] > 0)
      entropy += -((float)histogram[i]/im) * log((float)histogram[i]/im);
  }
  iftFree(histogram);

  entropy /= log(size);
  return entropy;

}

iftImage *iftAltMixedSamplingSecondStage(iftMImage *img, iftImage *mask, int nsamples)
{
  iftImage  *out, *quadMask, *quadGrid;
  iftMImage *quadImg;
  int i,j,k, t, p, initX, initY, initZ, endX, endY, endZ, x, y, z, origp, quadValuesSize, indexQV, indexQuad;
  int *quadNSamples;
  int nquad, nquadz, quadImgSize, quadImgSizeX, quadImgSizeY, quadImgSizeXY;
  float *quadEntropy, totalEntropy, *quadValues;
  nquad = 2;
  nquadz = 2;
  indexQuad = 0;
  totalEntropy = 0.0;
  if (img->zsize == 1)
    nquadz = 1;
  quadEntropy = iftAllocFloatArray(nquad*nquad*nquadz);
  quadNSamples = iftAllocIntArray(nquad*nquad*nquadz);

  out = iftCreateImage(img->xsize, img->ysize, img->zsize);

  for(k=0; k<nquadz; k++) {
    for (i = 0; i < nquad; i++) {
      for (j = 0; j < nquad; j++) {
        // Compute init and end coordinates
        initZ = k * (img->zsize / nquadz);
        initY = i * (img->ysize / nquad);
        initX = j * (img->xsize / nquad);
        endZ = (k + 1) * (img->zsize / nquadz);
        endY = (i + 1) * (img->ysize / nquad);
        endX = (j + 1) * (img->xsize / nquad);
        if (k == (nquadz - 1))
          endZ = img->zsize;
        if (img->zsize == 1)
          endZ = img->zsize;
        if (i == (nquad - 1))
          endY = img->ysize;
        if (j == (nquad - 1))
          endX = img->xsize;

        // Divide image in quadrants
        quadImgSizeY = (endY - initY);
        quadImgSizeX = (endX - initX);
        quadImgSizeXY = quadImgSizeX * quadImgSizeY;
        quadImgSize = quadImgSizeXY * (endZ - initZ);
        quadValues = iftAllocFloatArray(quadImgSize);
        indexQV = 0;
        for (p = 0; p < quadImgSize; p++) {
          z = p / quadImgSizeXY;
          y = (p % quadImgSizeXY) / quadImgSizeX;
          x = (p % quadImgSizeXY) % quadImgSizeX;
          origp = (x + initX) + (y + initY) * img->xsize + (z + initZ) * (img->xsize * img->ysize);
          quadValues[indexQV] = img->val[origp][0];
          indexQV++;
        }
        quadValuesSize = indexQV;
        // Compute entropy
        quadEntropy[indexQuad] = iftNormalizedShannonEntropy(quadValues, quadValuesSize);
        totalEntropy += quadEntropy[indexQuad];

        indexQuad++;
        iftFree(quadValues);
      }
    }
  }

  indexQuad = 0;
  for(k=0; k<nquadz; k++) {
    for (i = 0; i < nquad; i++) {
      for (j = 0; j < nquad; j++) {
        if (totalEntropy == 0)
          quadNSamples[indexQuad] = iftRound(((float) nsamples / (float) (nquad * nquad * nquadz) ));
        else
          quadNSamples[indexQuad] = iftRound((quadEntropy[indexQuad] / totalEntropy) * ((float) nsamples));

        if (quadNSamples[indexQuad] == 0)
          quadNSamples[indexQuad] = 1;

        // Compute init and end coordinates
        initZ = k * (img->zsize / nquadz);
        initY = i * (img->ysize / nquad);
        initX = j * (img->xsize / nquad);
        endZ = (k + 1) * (img->zsize / nquadz);
        endY = (i + 1) * (img->ysize / nquad);
        endX = (j + 1) * (img->xsize / nquad);
        if (k == (nquadz - 1))
          endZ = img->zsize;
        if (img->zsize == 1)
          endZ = img->zsize;
        if (i == (nquad - 1))
          endY = img->ysize;
        if (j == (nquad - 1))
          endX = img->xsize;

        // Divide image in quadrants
        quadImgSizeY = (endY - initY);
        quadImgSizeX = (endX - initX);
        quadImgSizeXY = quadImgSizeX * quadImgSizeY;
        quadImg = iftCreateMImage(quadImgSizeX, quadImgSizeY, (endZ - initZ), img->m);
        quadMask = iftCreateImage(quadImgSizeX, quadImgSizeY, (endZ - initZ));

        for (p = 0; p < quadImg->n; p++) {
          z = p / quadImgSizeXY;
          y = (p % quadImgSizeXY) / quadImgSizeX;
          x = (p % quadImgSizeXY) % quadImgSizeX;
          origp = (x + initX) + (y + initY) * img->xsize + (z + initZ) * (img->xsize * img->ysize);
          for (t = 0; t < quadImg->m; t++) {
            quadImg->val[p][t] = img->val[origp][t];
          }
          if (mask->val[origp] != 0)
            quadMask->val[p] = 1;
        }
        quadGrid = iftGridSampling(quadImg, quadMask, quadNSamples[indexQuad]);
        // Write seeds in the final out image
        for (p = 0; p < quadGrid->n; p++) {
          if (quadGrid->val[p] > 0) {
            z = p / (quadGrid->xsize*quadGrid->ysize);
            y = (p % (quadGrid->xsize*quadGrid->ysize)) / quadGrid->xsize;
            x = (p % (quadGrid->xsize*quadGrid->ysize)) % quadGrid->xsize;
            origp = (x + initX) + (y + initY) * img->xsize + (z + initZ) * (img->xsize * img->ysize);
            out->val[origp] = 1;
          }
        }

        indexQuad++;
        iftDestroyMImage(&quadImg);
        iftDestroyImage(&quadMask);
        iftDestroyImage(&quadGrid);

      }
    }
  }
  iftFree(quadEntropy);
  iftFree(quadNSamples);

  return out;
}

iftImage *iftAltMixedSampling(iftMImage *img, iftImage *mask, int nsamples)
{
  iftImage  *out, *quadMask, *quadGrid;
  iftMImage *quadImg;
  int i,j,k, t, p, initX, initY, initZ, endX, endY, endZ, x, y, z, origp, quadValuesSize, indexQV, indexQuad;
  int *quadNSamples;
  int nquad, nquadz, quadImgSize, quadImgSizeX, quadImgSizeY, quadImgSizeXY;
  float *quadEntropy, totalEntropy, *quadValues, meanEntropy, thrEntropy, sdEntropy;
  nquad = 2;
  nquadz =2;
  indexQuad = 0;
  totalEntropy = 0.0;
  if (img->zsize == 1)
    nquadz = 1;
  quadEntropy = iftAllocFloatArray(nquad*nquad*nquadz);
  quadNSamples = iftAllocIntArray(nquad*nquad*nquadz);

  out = iftCreateImage(img->xsize, img->ysize, img->zsize);

  for(k=0; k<nquadz; k++) {
    for(i=0; i<nquad; i++) {
      for(j=0; j<nquad; j++) {
        // Compute init and end coordinates
        initZ = k * (img->zsize / nquadz);
        initY = i * (img->ysize / nquad);
        initX = j * (img->xsize / nquad);
        endZ = (k + 1) * (img->zsize / nquadz);
        endY = (i + 1) * (img->ysize / nquad);
        endX = (j + 1) * (img->xsize / nquad);
        if (k == (nquadz - 1))
          endZ = img->zsize;
        if (img->zsize == 1)
          endZ = img->zsize;
        if (i == (nquad - 1))
          endY = img->ysize;
        if (j == (nquad - 1))
          endX = img->xsize;

        // Divide image in quadrants
        quadImgSizeY = (endY - initY);
        quadImgSizeX = (endX - initX);
        quadImgSizeXY = quadImgSizeX * quadImgSizeY;
        quadImgSize = quadImgSizeXY * (endZ - initZ);
        quadValues = iftAllocFloatArray(quadImgSize);
        indexQV = 0;
        for (p = 0; p < quadImgSize; p++) {
          z = p / quadImgSizeXY;
          y = (p % quadImgSizeXY) / quadImgSizeX;
          x = (p % quadImgSizeXY) % quadImgSizeX;
          origp = (x + initX) + (y + initY) * img->xsize + (z + initZ) * (img->xsize * img->ysize);
          quadValues[indexQV] = img->val[origp][0];
          indexQV++;
        }
        quadValuesSize = indexQV;
        // Compute entropy
        quadEntropy[indexQuad] = iftNormalizedShannonEntropy(quadValues, quadValuesSize);
        totalEntropy += quadEntropy[indexQuad];

        indexQuad++;
        iftFree(quadValues);
      }
    }
  }

  indexQuad = 0;
  meanEntropy = totalEntropy / (float)(nquad*nquad*nquadz);
  sdEntropy = 0.0;
  for(k=0; k<nquadz; k++) {
    for (i = 0; i < nquad; i++) {
      for (j = 0; j < nquad; j++) {
        sdEntropy = (quadEntropy[indexQuad] - meanEntropy) * (quadEntropy[indexQuad] - meanEntropy);
        indexQuad++;
      }
    }
  }
  sdEntropy /= (float)(nquad*nquad*nquadz-1);
  sdEntropy = sqrtf(sdEntropy);
  thrEntropy = meanEntropy + sdEntropy;

  indexQuad = 0;
  for(k=0; k<nquadz; k++) {
    for (i = 0; i < nquad; i++) {
      for (j = 0; j < nquad; j++) {
        if (totalEntropy == 0)
          quadNSamples[indexQuad] = iftRound(((float) nsamples / (float) (nquad * nquad * nquadz) ));
        else
          quadNSamples[indexQuad] = iftRound((quadEntropy[indexQuad] / totalEntropy) * ((float) nsamples));
        if (quadNSamples[indexQuad] == 0)
          quadNSamples[indexQuad] = 1;

        // Compute init and end coordinates
        initZ = k * (img->zsize / nquadz);
        initY = i * (img->ysize / nquad);
        initX = j * (img->xsize / nquad);
        endZ = (k + 1) * (img->zsize / nquadz);
        endY = (i + 1) * (img->ysize / nquad);
        endX = (j + 1) * (img->xsize / nquad);
        if (k == (nquadz - 1))
          endZ = img->zsize;
        if (img->zsize == 1)
          endZ = img->zsize;
        if (i == (nquad - 1))
          endY = img->ysize;
        if (j == (nquad - 1))
          endX = img->xsize;

        // Divide image in quadrants
        quadImgSizeY = (endY - initY);
        quadImgSizeX = (endX - initX);
        quadImgSizeXY = quadImgSizeX * quadImgSizeY;
        quadImg = iftCreateMImage(quadImgSizeX, quadImgSizeY, (endZ - initZ), img->m);
        quadMask = iftCreateImage(quadImgSizeX, quadImgSizeY, (endZ - initZ));
        for (p = 0; p < quadImg->n; p++) {
          z = p / quadImgSizeXY;
          y = (p % quadImgSizeXY) / quadImgSizeX;
          x = (p % quadImgSizeXY) % quadImgSizeX;
          origp = (x + initX) + (y + initY) * img->xsize + (z + initZ) * (img->xsize * img->ysize);
          for (t = 0; t < quadImg->m; t++) {
            quadImg->val[p][t] = img->val[origp][t];
          }
          if (mask->val[origp] != 0)
            quadMask->val[p] = 1;
        }

        if (quadEntropy[indexQuad] > thrEntropy) {
          // Execute Second Stage of mix sampling
          quadGrid = iftAltMixedSamplingSecondStage(quadImg, quadMask, quadNSamples[indexQuad]);
        } else {
          quadGrid = iftGridSampling(quadImg, quadMask, quadNSamples[indexQuad]);
        }

        // Write seeds in the final out image
        for (p = 0; p < quadGrid->n; p++) {
          if (quadGrid->val[p] > 0) {
            z = p / (quadGrid->xsize*quadGrid->ysize);
            y = (p % (quadGrid->xsize*quadGrid->ysize)) / quadGrid->xsize;
            x = (p % (quadGrid->xsize*quadGrid->ysize)) % quadGrid->xsize;
            origp = (x + initX) + (y + initY) * img->xsize + (z + initZ) * (img->xsize * img->ysize);
            out->val[origp] = 1;
          }
        }

        indexQuad++;
        iftDestroyMImage(&quadImg);
        iftDestroyImage(&quadMask);
        iftDestroyImage(&quadGrid);

      }
    }
  }
  iftFree(quadEntropy);
  iftFree(quadNSamples);

  return out;
}

iftImage *iftBorderProbImage(iftMImage *img) 
{
  iftAdjRel *A;
  iftImage  *prob;

  if (iftIs3DMImage(img))
    A = iftSpheric(sqrtf(3.0));
  else
    A = iftCircular(1.5);

  prob = iftMImageBasins(img, A);

  int      prob_max_val = iftMaximumValue(prob);
  for (int p            = 0; p < prob->n; p++)
    prob->val[p] = (int) ((float) prob->val[p] / prob_max_val * 100.0);

  iftDestroyAdjRel(&A);

  return (prob);
}

iftImage *iftMImageBasins(const iftMImage *img, iftAdjRel *A)
{
   iftImage   *basins=iftCreateImage(img->xsize,img->ysize,img->zsize);
   double     *grad=iftAllocDoubleArray(img->n);
   float      *w, wt=0.0, K;
   int         dx, dy, dz;

   iftMaxAdjShifts(A, &dx, &dy, &dz);
   K = sqrtf(dx*dx + dy*dy + dz*dz);

   w = iftAllocFloatArray(A->n);
   for (int i=1; i < A->n; i++) {
     w[i] = K / sqrtf(A->dx[i]*A->dx[i] + A->dy[i]*A->dy[i] + A->dz[i]*A->dz[i]);
     wt  += w[i];
   }
   for (int i=1; i < A->n; i++) {
     w[i] = w[i]/wt;
   }

   /* Compute basins image in the spatial domain */

#pragma omp parallel for shared(img,grad,A)
   for (int p=0; p < img->n; p++) {
     iftVoxel u   = iftMGetVoxelCoord(img,p);
     for (int i=1; i < A->n; i++) {
       iftVoxel v = iftGetAdjacentVoxel(A,u,i);
       double dist=0.0;
       if (iftMValidVoxel(img,v)){
   int q = iftMGetVoxelIndex(img,v);
   for (int b=0; b < img->m; b++) {
     dist += fabs(img->val[q][b]-img->val[p][b]);
   }
       }
       grad[p] += dist*w[i];
     }
   }
   
#pragma omp parallel for shared(grad,basins)
   for (int p=0; p < img->n; p++) {
     basins->val[p] = iftRound(grad[p]);
   }

   iftFree(grad);
   iftFree(w);
   
   return(basins);
 }

inline iftVoxel iftMGetVoxelCoord(const iftMImage *img, int p)
{
    /* old
     * u.x = (((p) % (((img)->xsize)*((img)->ysize))) % (img)->xsize)
     * u.y = (((p) % (((img)->xsize)*((img)->ysize))) / (img)->xsize)
     * u.z = ((p) / (((img)->xsize)*((img)->ysize)))
     */
    iftVoxel u;
    div_t res1 = div(p, img->xsize * img->ysize);
    div_t res2 = div(res1.rem, img->xsize);

    u.x = res2.rem;
    u.y = res2.quot;
    u.z = res1.quot;

  return(u);
}

inline char iftMValidVoxel(const iftMImage *img, iftVoxel v) 
{
  if ((v.x >= 0)&&(v.x < img->xsize)&&
      (v.y >= 0)&&(v.y < img->ysize)&&
      (v.z >= 0)&&(v.z < img->zsize))
    return(1);
  else
    return(0);
}

// ---------- iftMImage.c end
// ---------- iftMemory.c start 

#ifdef __linux__
#include <sys/sysinfo.h>
#include <malloc.h>
#endif

#ifdef __APPLE__
#include <mach/task.h>
#include <mach/mach_init.h>
#endif

#ifdef _WINDOWS
#include <windows.h>
#else
#include <sys/resource.h>
#endif

int *iftAllocIntArray(long n) 
{
    int *v = NULL;

    v = (int *) iftAlloc(n, sizeof(int));
    if (v == NULL)
        iftError("Cannot allocate memory space", "iftAllocIntArray");
    return(v);
}

void iftCopyIntArray(int *array_dst, const int *array_src, int nelems) 
{
    #pragma omp parallel for
    for (int i = 0; i < nelems; i++) {
        array_dst[i] = array_src[i];
    }
}

#ifndef  __cplusplus
long long *iftAllocLongLongIntArray(long n) 
{
    long long *v = NULL;

    v = (long long *) iftAlloc(n, sizeof(long long));
    if (v == NULL)
        iftError("Cannot allocate memory space", "iftAllocLongLongIntArray");

    return v;
}

void iftCopyLongLongIntArray(long long *array_dst, const long long *array_src, int nelems) 
{
    #pragma omp parallel for
    for (int i = 0; i < nelems; ++i) {
        array_dst[i] = array_src[i];
    }
}
#endif

float *iftAllocFloatArray(long n) 
{
    float *v = NULL;
    v = (float *) iftAlloc(n, sizeof(float));
    if (v == NULL)
        iftError("Cannot allocate memory space", "iftAllocFloatArray");
    return(v);
}

void iftCopyFloatArray(float *array_dst, float *array_src, int nelems) 
{
//    int i;
//
//    for (i = 0; i < nelems; i++) {
//      array_dst[i] = array_src[i];
//    }

    memmove(array_dst, array_src, nelems*sizeof(float));
}

double *iftAllocDoubleArray(long n) 
{
    double *v = NULL;

    v = (double *) iftAlloc(n, sizeof(double));
    if (v == NULL)
        iftError("Cannot allocate memory space", "iftAllocDoubleArray");
    return (v);
}

void iftCopyDoubleArray(double *array_dst, double *array_src, int nelems) 
{
    #pragma omp parallel for
    for (int i = 0; i < nelems; i++)
        array_dst[i] = array_src[i];
}

ushort *iftAllocUShortArray(long n) 
{
    ushort *v = NULL;

    v = (ushort *) iftAlloc(n, sizeof(ushort));
    if (v == NULL)
        iftError("Cannot allocate memory space", "iftAllocUShortArray");
    return(v);
}

uchar *iftAllocUCharArray(long n) 
{
    uchar *v = NULL;

    v = (uchar *) iftAlloc(n, sizeof(uchar));
    if (v == NULL)
        iftError("Cannot allocate memory space", "iftAllocUCharArray");
    return (v);
}

char *iftAllocCharArray(long n) 
{
    char *v = NULL;

    v = (char *) iftAlloc(n, sizeof(char));
    if (v == NULL)
        iftError("Cannot allocate memory space", "iftAllocCharArray");
    return (v);
}

char *iftAllocString(long n) 
{
    return (iftAllocCharArray(n+1));
}

// ---------- iftMemory.c end
// ---------- iftMetrics.c start

float iftFScoreGivenErrors(iftErrorClassification* error) 
{
    float precision = iftPrecisionGivenErrors(error);
    float recall = iftRecallGivenErrors(error);

    return (2 * precision * recall) / (precision + recall);
}

float iftPrecisionGivenErrors(iftErrorClassification* error) 
{
    return ((float) error->tp) / (error->tp + error->fp);
}

float iftRecallGivenErrors(iftErrorClassification* error)
{
    return ((float) error->tp) / (error->tp + error->fn);
}

// ---------- iftMetrics.c end 
// ---------- iftSegmentation.c start 

float iftBoundaryRecall(iftImage *gt, iftImage *border, float tolerance_dist)
{
  iftVerifyImageDomains(gt, border,"iftBoundaryRecall");
  float number_of_gt_voxels=0,number_of_matchings=0; 
  int   p, q, i;
  iftVoxel u, v;
  iftAdjRel *A;

  if (iftIs3DImage(gt))
    A = iftSpheric(tolerance_dist);
  else
    A = iftCircular(tolerance_dist);

  for(p = 0; p < gt->n; p++)
    if (gt->val[p]!=0){      
      number_of_gt_voxels++;
      u = iftGetVoxelCoord(gt,p);
        for (i=0; i < A->n; i++) {
        v = iftGetAdjacentVoxel(A,u,i);
        if (iftValidVoxel(gt,v)){
          q = iftGetVoxelIndex(gt,v);
          if (border->val[q]!=0){
            number_of_matchings++;
            break;
          }
        }
      }
    }

  if (number_of_gt_voxels == 0.0)
      iftError("Empty ground-truth image", "iftBoundaryRecall");

  iftDestroyAdjRel(&A);

  return (number_of_matchings/number_of_gt_voxels);
}

float iftUnderSegmentation(iftImage *gt_image, iftImage *label_image)
{
  iftVerifyImageDomains(gt_image, label_image, "iftUnderSegmentationMin");
  int i, j, p, num_obj, num_regions;
  float area, total_err;
  num_obj = iftMaximumValue(gt_image) + 1;
  num_regions = iftMaximumValue(label_image);
  
  int *num_pix_total = iftAllocIntArray(num_regions);
  for(p = 0; p < label_image->n; p++)
    num_pix_total[label_image->val[p] - 1]++;
  
  int *num_pix_obj = iftAllocIntArray(num_regions);
  
  int sum_err = 0;
  for(j = 0; j < num_obj; j++){
    for(p = 0; p < label_image->n; p++){
      if(gt_image->val[p] == j){
        num_pix_obj[label_image->val[p] - 1]++;
      }
    }
    for(i = 0; i < num_regions; i++){
      area =  num_pix_obj[i];
      if((num_pix_total[i] - num_pix_obj[i]) < area){
        area = (num_pix_total[i] - num_pix_obj[i]);
      }
      sum_err = sum_err + area;
      num_pix_obj[i] = 0;
    }
  }
  
  total_err = sum_err/(float)label_image->n;
  
  free(num_pix_total);
  free(num_pix_obj);
  
  return total_err;
}

iftImage *iftBorderImage(const iftImage *label, bool get_margins)
{
 iftAdjRel *A;
 iftImage  *border = iftCreateImage(label->xsize,label->ysize,label->zsize);
 int        p,q,i; 
 iftVoxel   u, v;
    
  if (iftIs3DImage(label))
    A = iftSpheric(1.0);
  else
    A = iftCircular(1.0);

  if (get_margins){
    for(p=0; p < label->n; p++){
      u = iftGetVoxelCoord(label, p);
      for(i=1; i < A->n; i++){
        v = iftGetAdjacentVoxel(A,u,i);
        if (iftValidVoxel(label, v)){
          q = iftGetVoxelIndex(label, v);
          if (label->val[p] != label->val[q]){
            border->val[p] = label->val[p];
            break;
          }
        } else {
          border->val[p] = label->val[p];
        }
      }
    }
  }
  else{
    for(p=0; p < label->n; p++) {
      u = iftGetVoxelCoord(label, p);
      for (i = 1; i < A->n; i++) {
        v = iftGetAdjacentVoxel(A, u, i);
        if (iftValidVoxel(label, v)) {
          q = iftGetVoxelIndex(label, v);
          if (label->val[p] != label->val[q]) {
            border->val[p] = label->val[p];
            break;
          }
        }
      }
    }
  }

    iftDestroyAdjRel(&A);
    return(border);
}

iftImage* iftBorderImageToLabelImage(iftImage* border)
{
  iftAdjRel *A; 
  iftImage  *marker = iftAddValue(border,1);
  iftImage  *label;
  
  if (iftIs3DImage(border))
    A = iftSpheric(sqrtf(3.0));
  else
    A = iftCircular(1);

  label = iftWaterGray(border, marker, A);

  iftDestroyAdjRel(&A);
  iftDestroyImage(&marker);

  return(label);
}

iftImage *iftWaterGray(iftImage *basins, iftImage *pathval, iftAdjRel  *A)
{
  iftImage   *label=NULL;
  iftGQueue  *Q=NULL;
  int         i,p,q,l=1,tmp;
  iftVoxel    u,v;
 
  // Initialization 
  
  label   = iftCreateImage(basins->xsize,basins->ysize,basins->zsize);
  Q       = iftCreateGQueue(iftMaximumValue(pathval)+2,pathval->n,pathval->val);

  for (p=0; p < basins->n; p++) {
    pathval->val[p] += 1;
    label->val[p]= IFT_NIL;
    iftInsertGQueue(&Q,p);
  }

  // Image Foresting Transform

  while(!iftEmptyGQueue(Q)) {
    p=iftRemoveGQueue(Q);

    if (label->val[p] == IFT_NIL) { // root voxel
      pathval->val[p]  -= 1;
      label->val[p]=l; l++;
    }

    basins->val[p] = pathval->val[p]; // set the reconstruction value

    u = iftGetVoxelCoord(basins,p);

    for (i=1; i < A->n; i++){
      v = iftGetAdjacentVoxel(A,u,i);
      if (iftValidVoxel(basins,v)){	
	q = iftGetVoxelIndex(basins,v);
	if (Q->L.elem[q].color != IFT_BLACK){
	  tmp = iftMax(pathval->val[p], basins->val[q]);
	  if (tmp < pathval->val[q]){ 
	    iftRemoveGQueueElem(Q,q);
	    label->val[q]      = label->val[p];
	    pathval->val[q]    = tmp;
	    iftInsertGQueue(&Q, q);
	  }
	}
      }
    }
  }
  
  iftDestroyGQueue(&Q);
  iftCopyVoxelSize(basins,label);

  return(label);
}

iftFImage *iftSmoothWeightImage(const iftImage *basins, float beta) 
{
    iftFImage *weight = iftCreateFImage(basins->xsize,basins->ysize,basins->zsize);

    #pragma omp parallel for
    for (int p = 0; p < basins->n; p++)
        weight->val[p] = 1.0 / (1.0 + (beta * basins->val[p])); // pseudo inverse
    

    return weight;
}

iftFImage *iftWeightNormFactor(const iftFImage *weight, iftAdjRel *A)
{
  iftFImage *norm_factor=iftCreateFImage(weight->xsize,weight->ysize,weight->zsize);

#pragma omp parallel for shared(weight,norm_factor)
  for (int p=0; p < weight->n; p++) {
    iftVoxel u = iftFGetVoxelCoord(weight,p);
    for (int i=1; i < A->n; i++) {
      iftVoxel v = iftGetAdjacentVoxel(A,u,i);
      if (iftFValidVoxel(weight,v)){
	     int q = iftFGetVoxelIndex(weight,v);
	     norm_factor->val[p] += weight->val[q];
      }
    }
  }

  return(norm_factor);
}

float *iftFScoreMultiLabel (iftImage *mlabel, iftImage *mlgt, int number_of_objects)
{
    int p = 0, k = 0;
    iftImage *tmp_label, *tmp_gt;

    float *output = iftAllocFloatArray(number_of_objects + 1);
    float accu_sum = 0.0;
    iftVerifyImageDomains(mlabel, mlgt, "iftFScoreMultiLabel");

    for (k = 1; k <= number_of_objects; k++)
    {
        tmp_label = iftCreateImage(mlabel->xsize, mlabel->ysize, mlabel->zsize);
        tmp_gt    = iftCreateImage(mlgt->xsize, mlgt->ysize, mlgt->zsize);
        for (p = 0; p < mlabel->n; p++)
        {
            if (mlabel->val[p] == k)
                tmp_label->val[p] = 1;
            if (mlgt->val[p] == k)
                tmp_gt->val[p] = 1;
        }
        output[k] = iftFScoreError(tmp_gt, tmp_label);
        accu_sum += output[k];
        iftDestroyImage(&tmp_label);
        iftDestroyImage(&tmp_gt);
    }
    output[0] = (float) accu_sum / (number_of_objects * 1.0);
    return output;
}

float iftCompactness2D(iftImage *label) 
{
  iftImage  *border = iftBorderImage(label,1);
  iftAdjRel *A=NULL;
  int      i, p , q, *area, number_of_regions;
  iftVoxel u, v; 
  float    *perim, co;
  
  number_of_regions = iftMaximumValue(label);
  perim = iftAllocFloatArray(number_of_regions);
  area  = iftAllocIntArray(number_of_regions);

  if (iftIs3DImage(label))
      iftError("Label image must be 2D", "iftCompactness2D");
  else
    A    = iftCircular(sqrtf(2.0));

  for (p=0; p < label->n; p++) {
    area[label->val[p] - 1]++; 
    if (border->val[p]) {
      u = iftGetVoxelCoord(label,p);
      for (i=1; i < A->n; i++) {
	v = iftGetAdjacentVoxel(A,u,i);
	if (iftValidVoxel(label,v)){
	  q = iftGetVoxelIndex(label,v);
	  if (border->val[p]==border->val[q]){
	    perim[border->val[p] - 1] += iftSmoothEuclideanDistance((float) iftSquaredVoxelDistance(u, v));
	  }
	}
      }
    }
  }
  for (i=0; i < number_of_regions; i++) 
    perim[i] /= 2.0;

  iftDestroyImage(&border);

  co = 0; 
  for (i=0; i < number_of_regions; i++) 
    if (perim[i] > 0.0) 
      co += (4.0 * IFT_PI * area[i] / (perim[i] * perim[i])) * ((float)area[i] / label->n);

  iftFree(area);
  iftFree(perim);
  iftDestroyAdjRel(&A);

  return(co);
}

float iftTopologyMeasure(iftImage *label) 
{
  iftMatrix *adjMatrix;
  int        number_of_regions, number_of_valid_regions;
  float     *average_size, topology, average_number_of_adjacent_segments;
  char      *valid_region;
  int        p, q, i, index, r, c, nadjs;
  iftAdjRel *A;
  iftVoxel   u, v;

  if (iftIs3DImage(label))
    A    = iftSpheric(1.0);
  else
    A    = iftCircular(1.0);

  number_of_regions = iftMaximumValue(label);
  adjMatrix         = iftCreateMatrix(number_of_regions, number_of_regions);
  valid_region      = iftAllocCharArray(number_of_regions);
  average_size      = iftAllocFloatArray(number_of_regions);

  /* Find the valid regions as those with no adjacent voxels outside
     the image domain */

  for (r=0; r < number_of_regions; r++)
    valid_region[r]=1;

  number_of_valid_regions = number_of_regions;

  for (p=0; p < label->n; p++) {
    if (valid_region[label->val[p] - 1]){
      u = iftGetVoxelCoord(label,p);
      for (i=1; i < A->n; i++) {
	v = iftGetAdjacentVoxel(A,u,i);
	if (!iftValidVoxel(label,v)){
	  valid_region[label->val[p] - 1]=0;
	  number_of_valid_regions--;
	  break;
	}
      }
    }
  }

  if (number_of_valid_regions==0){
    iftWarning("Infinity topology","iftTopologyMeasure");
    return(IFT_INFINITY_FLT);
  }

  /* For each valid region, compute in the adjacency matrix the size
     of its border segments with each adjacent region. */

  for (p=0; p < label->n; p++) {
    if (valid_region[label->val[p] - 1]){
      u = iftGetVoxelCoord(label,p);
      r = label->val[p] - 1;
      for (i=1; i < A->n; i++) {
	v = iftGetAdjacentVoxel(A,u,i);
	q = iftGetVoxelIndex(label,v);
	if (label->val[p] != label->val[q]) {
	  c     = label->val[q] - 1;
	  index = iftGetMatrixIndex(adjMatrix, r, c);
	  adjMatrix->val[index]++;
	}
      }
    }
  }

  /* For each valid region, compute the average size of its border
     segments between regions. Compute also the average number of adjacent segments. */

  average_number_of_adjacent_segments=0.0;
 
  for (r=0; r < adjMatrix->nrows; r++) {
    if (valid_region[r]){
      nadjs = 0;
      for (c=0; c < adjMatrix->ncols; c++) {
	index = iftGetMatrixIndex(adjMatrix, r, c);
	if (adjMatrix->val[index]>0.0){
	  average_size[r] += adjMatrix->val[index];
	  nadjs++;
	}
      }
      average_number_of_adjacent_segments += nadjs;
      average_size[r] /= nadjs;
    }
  }

  average_number_of_adjacent_segments /= number_of_valid_regions;

  /* compute the topology measure */

  topology = 0;
  for (r=0; r < adjMatrix->nrows; r++) {
    if (valid_region[r]){
      for (c=0; c < adjMatrix->ncols; c++) {
	index = iftGetMatrixIndex(adjMatrix, r, c);
	if (adjMatrix->val[index]>0.0){
	  topology += fabs(adjMatrix->val[index]-average_size[r]);
	}
      }
    }
  }

  if (average_number_of_adjacent_segments > IFT_EPSILON)
    topology /= (number_of_valid_regions * average_number_of_adjacent_segments);
  else{
    iftWarning("Infinity topology","iftTopologyMeasure");
    topology = IFT_INFINITY_FLT;
  }

  iftFree(average_size);
  iftDestroyMatrix(&adjMatrix);
  iftFree(valid_region);

  return(topology);
}

float iftFScoreError(iftImage *bin, iftImage *gt)
{
  float fscore = 0.0;
  iftErrorClassification errors;

  errors = iftSegmentationErrors(gt, bin);
  fscore = iftFScoreGivenErrors(&errors);

  return fscore;
}

iftErrorClassification iftSegmentationErrors(iftImage* gt_image, iftImage* cl_image)
{
  iftErrorClassification errors;
  errors.tp = 0; errors.fp = 0; errors.tn = 0; errors.fn = 0;

  int i;
  for(i = 0; i < cl_image->n; i++){
    //Object Voxel
    if(gt_image->val[i] != 0){
      if(cl_image->val[i] == gt_image->val[i])
        errors.tp++;
      else
        errors.fn++;
    }
    //Background Voxel
    else{
      if(cl_image->val[i] == 0)
        errors.tn++;
      else
        errors.fp++;
    }
  }
  return errors;
}

// ---------- iftSegmentation.c end
// ---------- iftSet.c start 

void iftInsertSet(iftSet **S, int elem)
{
    iftSet *p=NULL;
    
    p = (iftSet *) iftAlloc(1,sizeof(iftSet));
    if (p == NULL) iftError(MSG_MEMORY_ALLOC_ERROR, "iftInsertSet");
    if (*S == NULL){
        p->elem  = elem;
        p->next  = NULL;
    }else{
        p->elem  = elem;
        p->next  = *S;
    }
    *S = p;
}

int iftRemoveSet(iftSet **S)
{
    iftSet *p;
    int elem = IFT_NIL;
    
    if (*S != NULL){
        p    =  *S;
        elem = p->elem;
        *S   = p->next;
        iftFree(p);
    }
    
    return(elem);
}

void iftRemoveSetElem(iftSet **S, int elem)
{
    if (S == NULL || *S == NULL)
        return;
    
    iftSet *tmp = *S;
    
    if (tmp->elem == elem) {
        *S = tmp->next;
        iftFree(tmp);
    } else {
        while (tmp->next != NULL && tmp->next->elem != elem)
            tmp = tmp->next;
        if (tmp->next == NULL)
            return;
        iftSet *remove = tmp->next;
        tmp->next = remove->next;
        iftFree(remove);
    }
}

void iftDestroySet(iftSet **S)
{
    iftSet *p;
    while(*S != NULL){
        p = *S;
        *S = p->next;
        iftFree(p);
    }
    *S = NULL;
}

iftSet* 	iftSetUnion(iftSet *S1,iftSet *S2)
{
    iftSet *S = 0;
    
    iftSet *s = S1;
    while(s){
        iftInsertSet(&S,s->elem);
        s = s->next;
    }
    
    s = S2;
    while(s){
        iftUnionSetElem(&S,s->elem);
        s = s->next;
    }
    
    return S;
}

iftSet* 	iftSetConcat(iftSet *S1,iftSet *S2)
{
    iftSet *S = 0;
    
    iftSet *s = S1;
    while(s){
        iftInsertSet(&S,s->elem);
        s = s->next;
    }
    
    s = S2;
    while(s){
        iftInsertSet(&S,s->elem);
        s = s->next;
    }
    
    return S;
}

char iftUnionSetElem(iftSet **S, int elem)
{
    iftSet *aux=*S;
    
    while (aux != NULL) {
        if (aux->elem == elem)
            return(0);
        aux = aux->next;
    }
    iftInsertSet(S,elem);
    return(1);
}

void iftInvertSet(iftSet **S)
{
    iftSet *set=NULL;
    
    while (*S != NULL)
        iftInsertSet(&set,iftRemoveSet(S));
    *S = set;
}

int iftSetSize(const iftSet* S)
{
    const iftSet *s = S;
    
    int i = 0;
    while (s != NULL){
        i++;
        s = s->next;
    }
    
    return i;
    
}

iftSet* iftSetCopy(iftSet* S)
{
    return iftSetUnion(S,0);
}

int iftSetHasElement(iftSet *S, int elem)
{
    iftSet *s = S;
    while(s){
        if(s->elem == elem)
            return 1;
        
        s = s->next;
    }
    
    return 0;
}

// ---------- iftSet.c end
// ---------- iftSList.c start

iftSList *iftCreateSList() 
{
    iftSList *SL = (iftSList *) iftAlloc(1, sizeof(iftSList));

    // just to really force
    SL->n    = 0;
    SL->head = NULL;
    SL->tail = NULL;

    return SL;
}

void iftDestroySList(iftSList **SL) 
{
    if (SL != NULL) {
        iftSList *SL_aux = *SL;

        if (SL_aux != NULL) {
            iftSNode *snode = SL_aux->head;
            iftSNode *tnode = NULL;

            while (snode != NULL) {
                tnode = snode;
                snode = snode->next;

                if (tnode->elem != NULL)
                    iftFree(tnode->elem);
                iftFree(tnode);
            }
            iftFree(SL_aux);
            *SL = NULL;
        }
    }
}

void iftInsertSListIntoHead(iftSList *SL, const char *elem) 
{
    if (SL == NULL)
        iftError("The String Linked List SL is NULL. Allocated it first", "iftInsertSListIntoHead");

    iftSNode *snode = (iftSNode*) iftAlloc(1, sizeof(iftSNode));
    snode->elem     = iftAllocCharArray(512);
    snode->prev = NULL;
    snode->next     = NULL;
    strcpy(snode->elem, elem);


    // The String Linked List is empty
    if (SL->head == NULL) {
        SL->head = snode;
        SL->tail = snode;
        SL->n    = 1;
    }
    else {
        snode->next        = SL->head;
        SL->head->prev = snode;
        SL->head           = snode;
        SL->n++;
    }
}

void iftInsertSListIntoTail(iftSList *SL, const char *elem) 
{
    if (SL == NULL)
        iftError("The String Linked List SL is NULL. Allocated it first", "iftInsertSListIntoTail");
    if (elem == NULL)
        iftError("The Element to be Inserted is NULL", "iftInsertSListIntoTail");

    iftSNode *snode = (iftSNode*) iftAlloc(1, sizeof(iftSNode));
    snode->prev     = NULL;
    snode->next     = NULL;
    snode->elem     = iftAllocCharArray(strlen(elem) + 1);
    strcpy(snode->elem, elem);


    // The String Linked List is empty
    if (SL->head == NULL) {
        SL->head = snode;
        SL->tail = snode;
        SL->n    = 1;
    }
    else {
        SL->tail->next = snode;
        snode->prev    = SL->tail;
        SL->tail       = snode;
        SL->n++;
    }
}

char *iftRemoveSListHead(iftSList *SL) 
{
    if (SL == NULL)
        iftError("The String Linked List SL is NULL. Allocated it first", "iftRemoveSListHead");
    
    char *elem      = NULL;
    iftSNode *snode = NULL;

    // if there are elements
    if (SL->head != NULL) {
        snode    = SL->head;
        SL->head = SL->head->next;

        // checks if the list is empty now
        if (SL->head == NULL)
            SL->tail = NULL;
        else
            SL->head->prev = NULL;

        SL->n--;
        elem = snode->elem;
        snode->elem = NULL;
        iftFree(snode); // it does not deallocates snode->free
    }

    return elem;
}

char *iftRemoveSListTail(iftSList *SL) 
{
    if (SL == NULL)
        iftError("The String Linked List SL is NULL. Allocated it first", "iftRemoveSListTail");

    char *elem      = NULL;
    iftSNode *snode = NULL;

    // if there are elements
    if (SL->head != NULL) {
        snode    = SL->tail;
        SL->tail = SL->tail->prev;

        // checks if the list is empty now
        if (SL->tail == NULL)
            SL->head = NULL;
        else
            SL->tail->next = NULL;

        SL->n--;
        elem = snode->elem;
        snode->elem = NULL;
        iftFree(snode); // it does not deallocates snode->free
    }

    return elem;
}

// ---------- iftSList.c end
// ---------- iftDialog.c start 

void iftError(const char *msg, const char *func, ...) 
{
    va_list args;
    char    final_msg[4096];
    
    va_start(args, func);
    vsprintf(final_msg, msg, args);
    va_end(args);
    
    fprintf(stderr, "\nError in %s: \n%s\n", func, final_msg);
    fflush(stdout);
    exit(-1);
}

void iftWarning(const char *msg, const char *func, ...) 
{
    va_list args;
    char    final_msg[4096];
    
    va_start(args, func);
    vsprintf(final_msg, msg, args);
    va_end(args);
    
    fprintf(stdout, "\nWarning in %s: \n%s\n", func, final_msg);
}

void iftDeprecated(const char *old_function, const char *new_function, const char *message) 
{
    if(message != NULL) {
        fprintf(stderr,
                IFT_COLOR_YELLOW "\n-----\nThe function \"%s\" is deprecated.\nUse the function \"%s\".\n%s\n\n-----\n\n" IFT_COLOR_RESET,
                old_function, new_function, message);
    } else {
        fprintf(stderr,
                IFT_COLOR_YELLOW "\n-----\nThe function \"%s\" is deprecated.\nUse the function \"%s\".\n-----\n\n" IFT_COLOR_RESET,
                old_function, new_function);
    }
}

// ---------- iftDialog.c end
// ---------- iftStream.c start 

char *iftGetLine(FILE *stream) 
{
    if (stream == NULL || feof(stream))
        return NULL;
    
    size_t buffer_size = 0;
    char *line = NULL;
    
    // getline reads an entire line from stream, storing the text (including the newline and a terminating null character)
    // in a buffer and storing the buffer address in *line
    // If you set *line to a null pointer, and *buffer_size to zero, before the call, then getline
    // allocates the initial buffer for you by calling malloc
    // If an error occurs or end of file is reached without any bytes read, getline returns -1.
    if (getline(&line, &buffer_size, stream) == -1) {
        if (feof(stream))
            return NULL;
        else
            iftError("Error with getline command", "iftGetLine");
    }
    
    size_t str_len = strlen(line);
    if (line[str_len-1] == '\n')
        line[str_len-1] = '\0'; // without this, it reads the \n and does not put the \0 at the end
    
    return line;
}

// ---------- iftStream.c end
// ---------- iftString.c start 

void iftRightTrim(char* s, char c) 
{

    int idx = strlen(s) - 1;

    while(s[idx] == c) {
        idx--;
    }

    s[idx+1] = '\0';

}

iftSList *iftSplitString(const char *phrase, const char *delimiter) 
{
    if (phrase == NULL)
        iftError("String to be splitted is NULL", "iftSplitString");
    if (delimiter == NULL)
        iftError("Delimiter is NULL", "iftSplitString");

    char *buf       = iftAllocString(strlen(phrase)+1);
    const char *pr  = phrase;
    const char *loc = strstr(pr, delimiter); // pointer to first delimiter occurrence in the string

    size_t length = strlen(delimiter);
    size_t bytes;

    iftSList *SL = iftCreateSList();

    // build a list of sub-strings
    while (loc != NULL) {
        bytes = loc - pr;
        strncpy(buf, pr, bytes);
        buf[bytes] = '\0'; // \0 character must be added manually because strncpy does not do that, as opposed to other functions such as strcpy and sprintf

        iftInsertSListIntoTail(SL, buf);

        pr = loc + length;
        loc = strstr(pr, delimiter);
    }

    // Copies the last substring to the left of the last delimiter found OR
    // Copies the whole string if it doesn't have the delimiter 
    strcpy(buf, pr);
    iftInsertSListIntoTail(SL, buf);

    iftFree(buf);

    return SL;
}

char *iftLowerString(const char *str) 
{
    if (str == NULL)
        iftError("Input string is NULL", "iftLowerString");

    char *out_str = iftAllocCharArray(strlen(str)+1);

    for (size_t c = 0; c < strlen(str); c++)
        out_str[c] = tolower(str[c]);

    return out_str;
}

bool iftCompareStrings(const char *str1, const char *str2) 
{
    if (str1 == NULL)
        iftError("First String is NULL", "iftCompareStrings");
    if (str2 == NULL)
        iftError("Second String is NULL", "iftCompareStrings");

    return (strcmp(str1, str2) == 0);
}

char *iftSplitStringAt(const char *phrase, const char *delimiter, long position) 
{
    if (phrase == NULL)
        iftError("String to be splitted is NULL", "iftSplitStringAt");
    if (delimiter == NULL)
        iftError("Delimiter is NULL", "iftSplitStringAt");

    iftSList *SL = iftSplitString(phrase, delimiter);
    iftSNode *snode = NULL;

    // Copies the split sub-string of the position
    if (position >= 0) {
        if ((position+1) <= SL->n) {

            snode = SL->head;
            for (size_t i = 0; i < position; i++)
                snode = snode->next;
        }
        else {
            iftError("Invalid Position %ld\n-> Position Index must be < %ld\n",
                     "iftSplitStringAt", position, SL->n);
        }
    } else {
        if (labs(position) <= SL->n) {
            long real_pos = SL->n + position;
            snode = SL->tail;
            for (size_t i = SL->n-1; i > real_pos; i--) {
                snode = snode->prev;
            }
        }
        else {
            iftError("Invalid Negative Position %ld\n-> Negative Position Index must be >= %ld\n",
                     "iftSplitStringAt", position, -1 * SL->n);
        }
    }

    char *str = iftCopyString(snode->elem);

    iftDestroySList(&SL);

    return str;
}

char *iftCopyString(const char *format, ...) 
{
    va_list args;
    char str[IFT_STR_DEFAULT_SIZE];
    
    va_start(args, format);
    vsprintf(str, format, args);
    va_end(args);
    
    char *copy = iftAllocCharArray(strlen(str) + 1);
    strcpy(copy, str);

    return copy;
}

char *iftRemoveSuffix(const char *str, const char *suffix) 
{
    if (str == NULL)
        iftError("String is NULL", "iftRemoveSuffix");
    if (suffix == NULL)
        iftError("Suffix is NULL", "iftRemoveSuffix");

    if (!iftCompareStrings(suffix, "") && iftEndsWith(str, suffix)) {
        size_t shift = strlen(str) - strlen(suffix);
        char *out_str = iftCopyString(str);
        out_str[shift] = '\0';

        return (out_str);   
    }
    else {
        return iftCopyString(str);
    }
}

bool iftEndsWith(const char *str, const char *suffix) 
{
    if (str == NULL)
        iftError("String is NULL", "iftEndsWith");
    if (suffix == NULL)
        iftError("Suffix is NULL", "iftEndsWith");

    size_t len_suffix = strlen(suffix);
    size_t len_str    = strlen(str);

    if (len_suffix <= len_str) {
        size_t shift = len_str - len_suffix;
        return (strncmp(str+shift, suffix, len_suffix) == 0);
    }
    else
        return false;
}

bool iftStartsWith(const char *str, const char *prefix) 
{
    if (str == NULL)
        iftError("String is NULL", "iftStartsWith");
    if (prefix == NULL)
        iftError("Prefix is NULL", "iftStartsWith");

    size_t len_prefix = strlen(prefix);
    size_t len_str    = strlen(str);

    if (len_prefix <= len_str)
        return (strncmp(str, prefix, len_prefix) == 0);        
    else
        return false;
}

char *iftConcatStrings(int n, ...) 
{
    if (n <= 0)
        iftError("Number of Strings to be concatenated is <= 0", "iftConcatStrings");

    size_t out_str_size = 1; // '\0'

    // Counts the size of the concatenated string
    va_list strings;
    va_start(strings, n);
    for (int i = 0; i < n; i++)
        out_str_size += strlen(va_arg(strings, char*));
    va_end(strings);

    char *concat_str = iftAllocCharArray(out_str_size);

    va_start(strings, n);
    for (int i = 0; i < n; i++)
        strcat(concat_str, va_arg(strings, char*));
    va_end(strings);

    return concat_str;
}

char *iftRemovePrefix(const char *str, const char *prefix) 
{
    if (str == NULL)
        iftError("String is NULL", "iftRemovePrefix");
    if (prefix == NULL)
        iftError("Prefix is NULL", "iftRemovePrefix");

    if (!iftCompareStrings(prefix, "") && iftStartsWith(str, prefix)) {
        size_t shift = strlen(prefix);
        return (iftCopyString(str + shift));   
    }
    else 
        return (iftCopyString(str));
}

// ---------- iftString.c end
// ---------- iftIGraph.c start 

iftIGraph *iftMImageToIGraph(const iftMImage *img, const iftImage *mask)
{
    iftIGraph *igraph = (iftIGraph *)iftAlloc(1,sizeof(iftIGraph));
    int        p, i;

    igraph->nnodes  = iftNumberOfElements(mask);
    igraph->node    = (iftINode *)iftAlloc(igraph->nnodes,sizeof(iftINode));
    igraph->index   = iftCreateImage(img->xsize, img->ysize, img->zsize);
    igraph->nfeats  = img->m;
    igraph->feat    = (float **)iftAlloc(img->n,sizeof(float *));

    iftCopyVoxelSize(img, igraph->index);
    for (p=0; p < img->n; p++) {
        igraph->feat[p] = iftAllocFloatArray(igraph->nfeats);
        for (i=0; i < img->m; i++)
            igraph->feat[p][i] = img->val[p][i];
    }
    igraph->label   = iftAllocIntArray(img->n);
    igraph->marker  = iftAllocIntArray(img->n);
    igraph->root    = iftAllocIntArray(img->n);
    igraph->pred    = iftAllocIntArray(img->n);
    igraph->pvalue  = iftAllocDoubleArray(img->n);

    for (p=0, i=0; p < mask->n; p++) {
        igraph->index->val[p]     = IFT_NIL;
        if (mask->val[p]>0){
            igraph->node[i].adj     = NULL;
            igraph->node[i].voxel   = p;
            igraph->node[i].weight  = 0.0;
            igraph->index->val[p]   = i;
            i++;
        }
    }

    return(igraph);
}

iftIGraph *iftImplicitIGraph(iftMImage *img, const iftImage *mask, iftAdjRel *A)
{
    iftIGraph *igraph = iftMImageToIGraph(img,mask);

    igraph->A       = iftCopyAdjacency(A);
    igraph->type    = IMPLICIT;

    return(igraph);
}

void iftIGraphSetWeightForRegionSmoothing(iftIGraph *igraph, const iftImage *img)
{
    iftAdjRel *A;

    if (iftIs3DImage(img))
        A = iftSpheric(sqrtf(3.0));
    else
        A = iftCircular(sqrtf(2.0));

    iftImage  *grad      = iftImageGradientMagnitude(img,A);
    iftDestroyAdjRel(&A);

    iftFImage *weight    = iftSmoothWeightImage(grad,0.5);

    iftIGraphSetFWeight(igraph, weight);

    iftDestroyImage(&grad);
    iftDestroyFImage(&weight);
}

void iftIGraphSmoothRegions(iftIGraph *igraph, int num_smooth_iterations)
{
    iftImage  *prev_label,  *next_label;
    iftFImage *prev_weight, *next_weight, *norm_factor, *weight;
    float     *sum, max_membership;
    int        l, i, p, q, r, max_label, iter;
    iftVoxel   u, v;
    iftAdjRel *A = igraph->A;
    iftSet    *prev_frontier = NULL, *next_frontier = NULL, *S = NULL;
    iftBMap   *inFrontier;

    if (igraph->type != IMPLICIT)
        iftError("For implicit graphs only", "iftIGraphSmoothRegions");

    /* Initialization */

    prev_label  = iftIGraphLabel(igraph);
    next_label  = iftCopyImage(prev_label);
    weight      = iftIGraphWeight(igraph);
    norm_factor = iftWeightNormFactor(weight,A);
    inFrontier  = iftCreateBMap(prev_label->n);

    int prev_label_max_val = iftMaximumValue(prev_label);
    sum         = iftAllocFloatArray(prev_label_max_val + 1);
    prev_weight = iftCreateFImage(prev_label->xsize, prev_label->ysize, prev_label->zsize);
    next_weight = iftCreateFImage(next_label->xsize, next_label->ysize, next_label->zsize);
    prev_frontier = iftObjectBorderSet(prev_label, A);

    for (p = 0; p < prev_label->n; p++){
        prev_weight->val[p] = next_weight->val[p] = 1.0;
    }

    S = prev_frontier;
    while (S != NULL) {
        p = S->elem;
        iftBMapSet1(inFrontier,p);
        S = S->next;
    }

    /* Smooth frontier and reset its path values */

    for (iter = 0; iter < num_smooth_iterations; iter++)
    {
        while (prev_frontier != NULL)
        {
            p = iftRemoveSet(&prev_frontier);
            iftInsertSet(&next_frontier, p);
            u   = iftGetVoxelCoord(prev_label, p);

            for (l = 0; l <= prev_label_max_val; l++)
            {
                sum[l] = 0.0;
            }

            for (i = 1; i < A->n; i++)
            {
                v = iftGetAdjacentVoxel(A, u, i);
                if (iftValidVoxel(prev_label, v))
                {
                    q = iftGetVoxelIndex(prev_label, v);
                    sum[prev_label->val[q]] += prev_weight->val[q] * weight->val[q];
                    if (iftBMapValue(inFrontier, q) == 0) /* expand frontier */
                    {
                        if (igraph->pred[q] != IFT_NIL)
                        {
                            iftInsertSet(&next_frontier, q);
                            iftBMapSet1(inFrontier, q);
                        }
                    }
                }
            }

            for (l = 0; l <= prev_label_max_val; l++)
                sum[l]  = sum[l] / norm_factor->val[p];

            max_membership = IFT_INFINITY_FLT_NEG; max_label = IFT_NIL;
            for (l = 0; l <= prev_label_max_val; l++)
            {
                if (sum[l] > max_membership)
                {
                    max_membership = sum[l];
                    max_label      = l;
                }
            }
            next_label->val[p]  = max_label;
            next_weight->val[p] = sum[max_label];
        }

        prev_frontier = next_frontier;
        next_frontier = NULL;

        for (r = 0; r < prev_label->n; r++)
        {
            prev_weight->val[r] = next_weight->val[r];
            prev_label->val[r]  = next_label->val[r];
        }
    }

    iftFree(sum);
    iftDestroyFImage(&prev_weight);
    iftDestroyFImage(&next_weight);
    iftDestroyImage(&next_label);
    iftDestroyFImage(&weight);
    iftDestroyFImage(&norm_factor);
    iftDestroyBMap(&inFrontier);
    iftDestroySet(&prev_frontier);

    /* It fixes the label map, by eliminating the smallest regions and
       relabel them with the adjaceny labels */

    prev_label_max_val = iftMaximumValue(prev_label);
    next_label = iftSelectKLargestRegionsAndPropagateTheirLabels(prev_label, A, prev_label_max_val);
    for (p=0; p < next_label->n; p++)
        igraph->label[p]=next_label->val[p];

    iftDestroyImage(&next_label);
    iftDestroyImage(&prev_label);

}

void iftDestroyIGraph(iftIGraph **igraph)
{
    if(igraph != NULL && *igraph != NULL) {
        iftIGraph *aux = *igraph;
        int p, i;

        for (i = 0; i < aux->nnodes; i++) {
            if (aux->node[i].adj != NULL)
                iftDestroySet(&aux->node[i].adj);
        }
        for (p = 0; p < aux->index->n; p++) {
            iftFree(aux->feat[p]);
        }
        iftFree(aux->feat);
        iftFree(aux->label);
        iftFree(aux->marker);
        iftFree(aux->root);
        iftFree(aux->pred);
        iftFree(aux->pvalue);

        if (aux->type == IMPLICIT)
            iftDestroyAdjRel(&aux->A);

        iftFree(aux->node);
        iftDestroyImage(&aux->index);
        iftFree(aux);
        (*igraph) = NULL;
    }
}

iftImage *iftIGraphLabel(iftIGraph *igraph)
{
    int p;
    iftImage *label = iftCreateImage(igraph->index->xsize,igraph->index->ysize,igraph->index->zsize);

    for (p=0; p < label->n; p++) {
        label->val[p] = igraph->label[p];
    }

    iftCopyVoxelSize(igraph->index, label);

    return(label);
}

void iftIGraphSetFWeight(iftIGraph *igraph, iftFImage *weight)
{
    for (int s=0; s < igraph->nnodes; s++) {
        int p = igraph->node[s].voxel;
        igraph->node[s].weight = weight->val[p];
    }
}

iftFImage *iftIGraphWeight(iftIGraph *igraph)
{
    int p,s;
    iftFImage *weight = iftCreateFImage(igraph->index->xsize,igraph->index->ysize,igraph->index->zsize);

    for (s=0; s < igraph->nnodes; s++) {
        p = igraph->node[s].voxel;
        weight->val[p] = igraph->node[s].weight;
    }
    iftCopyVoxelSize(igraph->index, weight);

    return(weight);
}

// ---------- iftIGraph.c end
// ---------- iftDataSet.c start

float iftDistance1(float *f1, float *f2, float *alpha, int n)
{
    float dist=0.0f;
    for (int i=0; i < n; i++){
        dist += (f1[i]-f2[i])*(f1[i]-f2[i])*alpha[i];
    }

    return(sqrtf(dist));
}

float iftDistance10(float *f1, float *f2, float *alpha, int n) 
{
    int i;
    float dist=0.0f, sf1 = 0.0f, sf2 = 0.0f, num;

    for (i = 0; i < n; i++){
        sf1+=f1[i];
        sf2+=f2[i];
    }

    for (i=0; i < n; i++){
        num = f1[i] - f2[i];
        dist += (num*num)/(f1[i]+f2[i]+IFT_EPSILON);
    }
    return dist;
//    return(sqrtf(dist));
}

// ---------- iftDataSet.c end
// ---------- iftSeeds.c start

iftAdjRel *_iftFindBestAdjRelForGridSampling(float radius, bool is_3D_adj) 
{
    iftAdjRel *A = NULL;

    // Adjacency to compute the boundaries of the best adjacency relation for the grid sampling
    // We considered 8-neighborhood and 26-neighborhood to avoid possible leaks during sampling
    iftAdjRel *B = (is_3D_adj) ? iftSpheric(sqrtf(3.0)) : iftCircular(sqrtf(2.0));

    int it = 0;
    float epsilon = 0.1;
    float min_dist = 0.0;

    do {
        iftDestroyAdjRel(&A);

        A = (is_3D_adj) ? iftSpheric(radius + (it * epsilon)) : iftCircular(radius + (it * epsilon));
        it++;

        iftAdjRel *Abound = iftAdjacencyBoundaries(A, B);

        min_dist = IFT_INFINITY_FLT;
        for (int i = 0; i < Abound->n; i++) {
            float dist = sqrtf(Abound->dx[i]*Abound->dx[i] + Abound->dy[i]*Abound->dy[i] + Abound->dz[i]*Abound->dz[i]);
            if (dist < min_dist)
                min_dist = dist;
        }
        iftDestroyAdjRel(&Abound);
    } while (min_dist < radius);

    iftDestroyAdjRel(&B);    

    return A;
}

iftIntArray *iftGridSamplingOnMask(const iftImage *bin_mask, float radius, int initial_obj_voxel_idx, long n_samples) 
{
    int first_obj_voxel = initial_obj_voxel_idx;

    if (initial_obj_voxel_idx >= 0) {
      if (bin_mask->val[initial_obj_voxel_idx] == 0) {
        iftError("Initial Voxel Index %d is not an object voxel",
                 "iftGridSamplingOnMask", initial_obj_voxel_idx);
      }
    }
    else {
      // finds the first object voxel from the binary mask
      int p = 0;
      for (p = 0; p < bin_mask->n && bin_mask->val[p] == 0; p++) {}
      first_obj_voxel = p;
    }
                                     

    if (iftAlmostZero(radius)) {
      iftIntArray *grid_chosen = iftCreateIntArray(1);
      grid_chosen->val[0] = first_obj_voxel;
      return grid_chosen;
    }

    iftAdjRel *A = _iftFindBestAdjRelForGridSampling(radius, iftIs3DImage(bin_mask));
    iftFloatArray *dist = iftCreateFloatArray(A->n);
    for (int i = 0; i < A->n; i++)
        dist->val[i] = sqrtf(A->dx[i]*A->dx[i] + A->dy[i]*A->dy[i] + A->dz[i]*A->dz[i]);
    
    iftImage *prior = iftCreateImageFromImage(bin_mask);
    iftImage *label_img = iftCreateImageFromImage(bin_mask);
    iftGQueue *Q = iftCreateGQueue(IFT_QSIZE, prior->n, prior->val);
    iftSetRemovalPolicy(Q, MAXVALUE);
    
    prior->val[first_obj_voxel] = 1;
    iftInsertGQueue(&Q, first_obj_voxel);

    iftList *grid = iftCreateList();

    int label = 1;

    while (!iftEmptyGQueue(Q)) {
        int p = iftRemoveGQueue(Q);
        iftVoxel u = iftGetVoxelCoord(bin_mask, p);
        iftInsertListIntoTail(grid, p);

        for (int i = 0; i < A->n; i++) {
            iftVoxel v = iftGetAdjacentVoxel(A, u, i);

            if (iftValidVoxel(bin_mask, v) && (iftImgVoxelVal(bin_mask, v) != 0) && (iftImgVoxelVal(label_img, v) == 0)) {
                int q = iftGetVoxelIndex(bin_mask, v);

                // q is inside the sphere
                if (dist->val[i] < radius) {
                    label_img->val[q] = label;
                    
                    // voxel was border of another ball
                    if (Q->L.elem[q].color == IFT_GRAY) {
                        iftRemoveGQueueElem(Q, q);
                        prior->val[q] = 0;
                    }
                }
                // q is on the border of the ball with center p
                else {
                    if (Q->L.elem[q].color == IFT_WHITE) {
                        prior->val[q]++;
                        iftInsertGQueue(&Q, q);
                    }
                    // q is on the border of the intersection of two balls 
                    else if (Q->L.elem[q].color == IFT_GRAY) {
                        iftRemoveGQueueElem(Q, q);
                        prior->val[q]++;
                        iftInsertGQueue(&Q, q);
                    }
                }
            }
        }
        label++;
    }

    iftIntArray *grid_all = iftListToIntArray(grid);
    iftDestroyList(&grid);    

    iftIntArray *grid_chosen = NULL;
    if (n_samples <= 0)
        grid_chosen = grid_all;
    else if (n_samples >= grid_all->n) {
        printf("Warning: Number of required samples %ld is >= total number of sampling voxels %ld\n" \
               "All sampling points will be considered\n", n_samples, grid_all->n);
        grid_chosen = grid_all;
    }
    else {
        grid_chosen = iftCreateIntArray(n_samples);
        iftShuffleIntArray(grid_all->val, grid_all->n);

        #pragma omp parallel for
        for (long i = 0; i < n_samples; i++)
            grid_chosen->val[i] = grid_all->val[i];
        iftDestroyIntArray(&grid_all);
    }


    // cleaning up
    iftDestroyAdjRel(&A);
    iftDestroyImage(&prior);
    iftDestroyImage(&label_img);
    iftDestroyGQueue(&Q);
    iftDestroyFloatArray(&dist);

    return grid_chosen;
}

float iftEstimateGridOnMaskSamplingRadius(const iftImage *binMask, int initialObjVoxelIdx, int nSamples)
{
  if (nSamples == 1) { return 0.0; }

  bool is3D = iftIs3DImage(binMask);

  // Compute radius if each seed actually covered its entire radius 
  int totalArea = 0;

  for (int p = 0; p < binMask->n; ++p)
    if (binMask->val[p] != 0)
      totalArea += 1;
  double baseR;
  if (is3D)
    baseR = pow(((double)totalArea * 3.0) / (4.0 * IFT_PI * (double)nSamples), 1.0/3.0);
  else
    baseR = sqrt((double)totalArea/(IFT_PI * (double)nSamples));

  // Optimization method based on binary search
  double lowerBound = 0.0;
  double upperBound = IFT_INFINITY_FLT;
  int bestError = IFT_INFINITY_INT;
  int maxError = iftMax(nSamples/20, 5);
  // Arbitrary initial estimate
  double rad = baseR * (IFT_PI / 2.0);
  while (bestError > maxError && upperBound - lowerBound > IFT_EPSILON) {
    iftIntArray *seeds = iftGridSamplingOnMask(binMask, rad, initialObjVoxelIdx, 0);

    long error = labs(nSamples - seeds->n);
    if (error < bestError) {
      bestError = error;
    }

    if (seeds->n < nSamples && rad < upperBound)
      upperBound = rad;
    if (seeds->n > nSamples && rad > lowerBound)
      lowerBound = rad;

    if (upperBound != IFT_INFINITY_FLT)
      rad = (lowerBound + upperBound) / 2.0;
    else
      rad = rad * 2;

    iftDestroyIntArray(&seeds);
  }

  return rad;
}

iftSet *iftObjectBorderSet(const iftImage *label_img, iftAdjRel *Ain) 
{
    iftAdjRel *A = Ain;
    if (A == NULL) {
        if (iftIs3DImage(label_img))
            A = iftSpheric(1.0);
        else A = iftCircular(1.0);
    }

    iftSet *borders = NULL;

    for (int p = 0; p < label_img->n; p++) {
        if (label_img->val[p] != 0) {
            iftVoxel u = iftGetVoxelCoord(label_img, p);
            
            for (int i = 1; i < A->n; i++) {
                iftVoxel v = iftGetAdjacentVoxel(A, u, i);

                if (iftValidVoxel(label_img, v)) {
                    int q = iftGetVoxelIndex(label_img, v);
                    
                    if (label_img->val[q] != label_img->val[p]) {
                        iftInsertSet(&borders, p);
                        break;
                    }
                }
                else {
                    // p is an object pixel which is on the image's border
                    iftInsertSet(&borders, p);
                    break;
                }
            }
        }
    }


    if (Ain == NULL)
        iftDestroyAdjRel(&A);

    return borders;
}

iftImage *iftSelectKLargestRegionsAndPropagateTheirLabels(iftImage *label, iftAdjRel *A, int K)
{
  iftImage *nlabel[2];
  int       ncomps,p,i,*index,*size;
  iftSet *S = NULL;

  /* relabel components: it accounts for disconnected labels, which
     should be considered multiple components */

  nlabel[0] = iftRelabelRegions(label,A);

  ncomps = iftMaximumValue(nlabel[0]);
  size   = iftAllocIntArray(ncomps+1);
  index  = iftAllocIntArray(ncomps+1);

  if (K <= 0)
      iftError("Invalid number of components", "iftSelectKLargestRegions");

  if (K > ncomps)
    K = ncomps;

  for (i=0; i <= ncomps; i++) {
    index[i]=i;
    size[i] =0.0;
  }

  /* Do not consider label 0 --- background */
  for (p=0; p < nlabel[0]->n; p++)
    size[nlabel[0]->val[p]]++;
  size[0]=0;

  iftBucketSort(size, index, ncomps+1, IFT_DECREASING);
  iftFree(size);

  nlabel[1] = iftCreateImage(nlabel[0]->xsize,nlabel[0]->ysize,nlabel[0]->zsize);
  for (p=0; p < nlabel[0]->n; p++) {
    for (i=0; i < K; i++) {
      if (nlabel[0]->val[p]==index[i]){
	nlabel[1]->val[p]=i+1;
	iftInsertSet(&S,p);
      }
    }
  }

  iftCopyVoxelSize(label,nlabel[1]);
  iftFree(index);

  /* fill regions with label zero from the labeled voxels in the
     set */

  while(S != NULL) {
    int p      = iftRemoveSet(&S);
    iftVoxel u = iftGetVoxelCoord(nlabel[1],p);
    for (i=1; i < A->n; i++) {
      iftVoxel v = iftGetAdjacentVoxel(A,u,i);
      if (iftValidVoxel(nlabel[1],v)){
	int q = iftGetVoxelIndex(nlabel[1],v);
	if (nlabel[1]->val[q]==0){
	  nlabel[1]->val[q]=nlabel[1]->val[p];
	  iftInsertSet(&S,q);
	}
      }
    }
  }
	
  iftDestroySet(&S);
  iftDestroyImage(&nlabel[0]);

  return(nlabel[1]);
}

// ---------- iftSeeds.c end
// ---------- iftDir.c start 

int _iftCmpFiles(const void *a, const void *b) 
{
    iftFile **f1 = (iftFile**) a;
    iftFile **f2 = (iftFile**) b;
    
    return strcmp((*f1)->path, (*f2)->path);
}

int _iftCmpDirs(const void *a, const void *b) 
{
    iftDir **dir1 = (iftDir**) a;
    iftDir **dir2 = (iftDir**) b;
    
    return strcmp((*dir1)->path, (*dir2)->path);
}

void _iftCountFilesInDirectory(const char *dir_pathname, long *nfiles, long *nsubdirs) 
{
    //http://pubs.opengroup.org/onlinepubs/007908799/xsh/dirent.h.html
    //http://www.delorie.com/gnu/docs/glibc/libc_270.html
    DIR *system_dir;
    struct dirent *entry;
    char msg[512];
    char *pathname = NULL;
    *nfiles = 0;
    *nsubdirs = 0;
    
    system_dir = opendir(dir_pathname);
    if (system_dir) {
        while ((entry = readdir(system_dir)) != NULL)
            // it excludes the system_dir . and ..
            if ((strcmp(entry->d_name, ".") != 0) && (strcmp(entry->d_name, "..") != 0)) {
                pathname = iftJoinPathnames(2, dir_pathname, entry->d_name);
                
                if (iftDirExists(pathname)) {
                    (*nsubdirs)++;
                }
                else {
                    (*nfiles)++;
                }
                
                iftFree(pathname);
                pathname = NULL;
            }
        closedir(system_dir);
    }
    else {
        sprintf(msg, "Error opening directory path: \"%s\"", dir_pathname);
        iftError(msg, "_iftCountFilesInDirectory");
    }
}

void _iftListDirectoryRec(iftDir *dir, long hier_levels, long curr_level) 
{
    DIR *system_dir;
    struct dirent *entry;
    char *pathname = NULL;
    
    dir->files   = NULL;
    dir->subdirs = NULL;
    
    if (curr_level <= hier_levels) {
        system_dir = opendir(dir->path);
        if (system_dir) {
            _iftCountFilesInDirectory(dir->path, &dir->nfiles, &dir->nsubdirs);
            if (dir->nfiles != 0)
                dir->files   = (iftFile**) iftAlloc(dir->nfiles, sizeof(iftFile*));
            if (dir->nsubdirs != 0)
                dir->subdirs = (iftDir**) iftAlloc(dir->nsubdirs, sizeof(iftDir*));
            
            long i = 0, j = 0;
            while ((entry = readdir(system_dir)) != NULL) {
                // it excludes the dir . and ..
                if ((strcmp(entry->d_name, ".") != 0) && (strcmp(entry->d_name, "..") != 0)) {
                    pathname = iftJoinPathnames(2, dir->path, entry->d_name);
                    
                    if (iftDirExists(pathname)) { // it is a directory
                        iftDir *subdir = (iftDir*) iftAlloc(1, sizeof(iftDir));
                        subdir->path = pathname;
                        
                        subdir->nfiles   = 0;
                        subdir->nsubdirs = 0;
                        subdir->files    = NULL;
                        subdir->subdirs  = NULL;
                        
                        _iftListDirectoryRec(subdir, hier_levels, curr_level+1);
                        dir->subdirs[j++] = subdir;
                    }
                    else { // it is a File
                        iftFile *f = (iftFile*) iftAlloc(1, sizeof(iftFile));
                        f->path = pathname;
                        dir->files[i++] = f;
                        f->suffix = NULL;
                    }
                }
            }
            closedir(system_dir);
            
            /* sorts the pathnames using qsort functions */
            qsort(dir->files, dir->nfiles, sizeof(iftFile*), _iftCmpFiles);
            qsort(dir->subdirs, dir->nsubdirs, sizeof(iftDir*), _iftCmpDirs);
        }
        else {
            char msg[512];
            sprintf(msg, "Error opening directory path: \"%s\"", dir->path);
            iftError(msg, "_iftListDirectoryRec");
        }
    }
}

void _iftListDirectory(iftDir *root_dir, long hier_levels) 
{
    if (root_dir == NULL) {
        iftError("Directory is NULL", "_iftListDirectory");
    }
    else {
        if (hier_levels == 0)
            hier_levels = IFT_INFINITY_INT; // trick to set the hier_levels as the possible maximum
        
        root_dir->nfiles   = 0;
        root_dir->nsubdirs = 0;
        root_dir->files    = NULL;
        root_dir->subdirs  = NULL;
        
        long curr_level = 1;
        _iftListDirectoryRec(root_dir, hier_levels, curr_level);
    }
}

bool iftDirExists(const char *format, ...) 
{
    va_list args;
    char pathname[IFT_STR_DEFAULT_SIZE];
    
    va_start(args, format);
    vsprintf(pathname, format, args);
    va_end(args);
    
    struct stat st;
    
    if (stat(pathname, &st) == 0) {
        if (S_ISDIR(st.st_mode))
            return true; //it's a directory
    }
    return false ;
}

char *iftParentDir(const char *pathname) 
{
    char *filename   = NULL;
    char *parent_dir = NULL;
    
    filename   = iftSplitStringAt(pathname, IFT_SEP_C, -1);
    parent_dir = iftRemoveSuffix(pathname, filename);
    iftFree(filename);
    
    // if the parent_dir is empty
    if (strcmp(parent_dir, "") == 0) {
        strcpy(parent_dir, ".");
    }
    else { // eliminates the slash (dir. separator char) at the end
        parent_dir[strlen(parent_dir)-1] = '\0';
    }
    
    return (parent_dir);
}

void iftMakeDir(const char *dir_path) 
{
    if (!iftDirExists(dir_path)) {
        char *parent_dir = iftAllocCharArray(IFT_STR_DEFAULT_SIZE);
        strcpy(parent_dir, "");
        
        iftSList *SL    = iftSplitString(dir_path, IFT_SEP_C);
        char *inter_dir = iftRemoveSListHead(SL);
        
        while (inter_dir != NULL) {
            strcat(parent_dir, inter_dir);
            strcat(parent_dir, IFT_SEP_C);
            
            if (!iftDirExists(parent_dir)) {
            #if defined(__linux) || defined(__APPLE__)
                if (mkdir(parent_dir, 0777) == -1) // Create the directory
                    iftError("Problem to create the directory: %s", "iftMakeDir", dir_path);
            #else
                if (!CreateDirectory(parent_dir, NULL))
                        iftError("Problem to create the directory", "iftMakeDir");
            #endif
            }
            
            iftFree(inter_dir);
            inter_dir = iftRemoveSListHead(SL);
        }
        iftDestroySList(&SL);
        iftFree(parent_dir);
    }
}

iftDir *iftLoadFilesFromDirByRegex(const char *dir_pathname, const char *regex) 
{
    if (dir_pathname == NULL)
        iftError("Dir's Pathname is NULL", "iftLoadFilesFromDirByRegex");
    if (regex == NULL)
        iftError("Regex is NULL", "iftLoadFilesFromDirByRegex");
    
    iftDir *dir          = iftLoadDir(dir_pathname, 1);
    long n_all_files   = dir->nfiles;
    long n_final_files = 0;
    
    for (long f = 0; f < n_all_files; f++) {
        char *filename = iftFilename(dir->files[f]->path, NULL);
        
        if (iftRegexMatch(filename, regex))
            n_final_files++;
        
        iftFree(filename);
    }
    
    iftFile **file_array = dir->files;
    dir->files           = NULL;
    
    dir->files  = (iftFile**) iftAlloc(n_final_files, sizeof(iftFile*));
    dir->nfiles = n_final_files;
    
    long i = 0;
    for (long f = 0; f < n_all_files; f++) {
        char *filename = iftFilename(file_array[f]->path, NULL);
        
        if (iftRegexMatch(filename, regex)) {
            dir->files[i++] = file_array[f];
            file_array[f]   = NULL;
        }
        else iftDestroyFile(&file_array[f]);
        
        iftFree(filename);
    }
    iftFree(file_array);
    
    return dir;
}

iftDir *iftLoadDir(const char *dir_pathname, long hier_levels) 
{
    char msg[512];
    iftDir *dir = NULL;
    
    
    if (iftPathnameExists(dir_pathname)) {
        // it is really a directory and it exists
        if (iftDirExists(dir_pathname)) {
            dir = (iftDir*) iftAlloc(1, sizeof(iftDir));
            dir->path = iftAllocCharArray(strlen(dir_pathname) + 2); // one more char to put the separation '/'
            strcpy(dir->path, dir_pathname);
            
            // puts the '/' at the end of the pathname
            if (dir->path[strlen(dir->path) - 1] != IFT_SEP_C[0])
                strcat(dir->path, IFT_SEP_C);
            
            _iftListDirectory(dir, hier_levels);
        }
            // it is a File instead of a Directory
        else {
            sprintf(msg, "Pathname \"%s\" is a File", dir_pathname);
            iftError(msg, "iftLoadDir");
        }
    }
    else {
        sprintf(msg, "Pathname \"%s\" does not exist!", dir_pathname);
        iftError(msg, "iftLoadDir");
    }
    
    return dir;
}

void iftDestroyDir(iftDir **dir) 
{
    if (dir != NULL) {
        iftDir *dir_aux = *dir;
        
        if (dir_aux != NULL) {
            if (dir_aux->path != NULL)
                iftFree(dir_aux->path);
            
            // deallocates the files
            if (dir_aux->files != NULL) {
                for (long i = 0; i < dir_aux->nfiles; i++)
                    iftDestroyFile(&dir_aux->files[i]);
                
                iftFree(dir_aux->files);
                dir_aux->files = NULL;
            }
            
            if (dir_aux->subdirs != NULL) {
                // deallocates the subdirs
                for (long j = 0; j < dir_aux->nsubdirs; j++)
                    iftDestroyDir(&dir_aux->subdirs[j]);
                
                iftFree(dir_aux->subdirs);
                dir_aux->subdirs = NULL;
            }
            iftFree(dir_aux);
            *dir = NULL;
        }
    }
}

// ---------- iftDir.c end
// ---------- iftCompression.c start 

iftGZipFile iftGZipOpen(const char *path, const char *mode, bool use_compression)
{
    iftGZipFile file;
    file = (iftGZipFile) iftAlloc(1, sizeof(struct iftGZipPtr));
    if( file == NULL ){
        fprintf(stderr,"** IFT_ERROR: iftGZipOpen failed to alloc znzptr\n");
        return NULL;
    }

    file->nzfptr = NULL;

#ifdef HAVE_ZLIB
    file->zfptr = NULL;

    if (use_compression) {
        file->withz = 1;
        if((file->zfptr = gzopen(path,mode)) == NULL) {
            iftFree(file);
            file = NULL;
        }
    } else {
#endif

        file->withz = 0;
        if((file->nzfptr = fopen(path,mode)) == NULL) {
            iftFree(file);
            file = NULL;
        }

#ifdef HAVE_ZLIB
    }
#endif

    return file;
}

char *iftGZipGets(char *str, int size, iftGZipFile file)
{
    if (file==NULL) { return NULL; }
#ifdef HAVE_ZLIB
    if (file->zfptr!=NULL) return gzgets(file->zfptr,str,size);
#endif
    return fgets(str,size,file->nzfptr);
}

int iftGZipClose(iftGZipFile *file)
{
    int retval = 0;
    if (*file!=NULL) {
#ifdef HAVE_ZLIB
        if ((*file)->zfptr!=NULL)  { retval = gzclose((*file)->zfptr); }
#endif
        if ((*file)->nzfptr!=NULL) { retval = fclose((*file)->nzfptr); }

        iftFree(*file);
        *file = NULL;
    }
    return retval;
}

size_t iftGZipRead(void *buf, size_t size, size_t nmemb, iftGZipFile file)
{
    if (file==NULL) { return 0; }
#ifdef HAVE_ZLIB
    size_t     remain = size*nmemb;
    char     * cbuf = (char *)buf;
    unsigned   n2read;
    int        nread;

    if (file->zfptr!=NULL) {
        /* gzread/write take unsigned int length, so maybe read in int pieces
           (noted by M Hanke, example given by M Adler)   6 July 2010 [rickr] */
        while( remain > 0 ) {
            n2read = (remain < GZip_MAX_BLOCK_SIZE) ? remain : GZip_MAX_BLOCK_SIZE;
            nread = gzread(file->zfptr, (void *)cbuf, n2read);
            if( nread < 0 ) return nread; /* returns -1 on error */

            remain -= nread;
            cbuf += nread;

            /* require reading n2read bytes, so we don't get stuck */
            if( nread < (int)n2read ) break;  /* return will be short */
        }

        /* warn of a short read that will seem complete */
        if( remain > 0 && remain < size )
            fprintf(stderr,"** iftGZipRead: read short by %u bytes\n",(unsigned)remain);

        return nmemb - remain/size;   /* return number of members processed */
    }
#endif
    return fread(buf,size,nmemb,file->nzfptr);
}

int iftGZipPuts(const char *str, iftGZipFile file)
{
    if (file==NULL) { return 0; }
#ifdef HAVE_ZLIB
    if (file->zfptr!=NULL) return gzputs(file->zfptr,str);
#endif
    return fputs(str,file->nzfptr);
}

size_t iftGZipWrite(const void *buf, size_t size, size_t nmemb, iftGZipFile file)
{
    if (file==NULL) { return 0; }
#ifdef HAVE_ZLIB
    size_t     remain = size*nmemb;
    char     * cbuf = (char *)buf;
    unsigned   n2write;
    int        nwritten;

    if (file->zfptr!=NULL) {
        while( remain > 0 ) {
            n2write = (remain < GZip_MAX_BLOCK_SIZE) ? remain : GZip_MAX_BLOCK_SIZE;
            nwritten = gzwrite(file->zfptr, (void *)cbuf, n2write);

            /* gzread returns 0 on error, but in case that ever changes... */
            if( nwritten < 0 ) return nwritten;

            remain -= nwritten;
            cbuf += nwritten;

            /* require writing n2write bytes, so we don't get stuck */
            if( nwritten < (int)n2write ) break;
        }

        /* warn of a short write that will seem complete */
        if( remain > 0 && remain < size )
            fprintf(stderr,"** iftGZipWrite: write short by %u bytes\n",(unsigned)remain);

        return nmemb - remain/size;   /* return number of members processed */
    }
#endif
    return fwrite(buf,size,nmemb,file->nzfptr);
}

// ---------- iftCompression.c end
// ---------- iftRegex.c start

bool iftRegexMatch(const char *str, const char *regex_pattern, ...) {
    if (str == NULL)
        iftError("String is NULL", "iftRegexMatch");
    if (regex_pattern == NULL)
        iftError("Regular Expression is NULL", "iftRegexMatch");
    
    char error[IFT_STR_DEFAULT_SIZE];
    regex_t regex;
    int reti;
    
    va_list args;
    char final_regex_pattern[IFT_STR_DEFAULT_SIZE];
    
    va_start(args, regex_pattern);
    vsprintf(final_regex_pattern, regex_pattern, args);
    va_end(args);
    
    // Compile Regex
    if ((reti = regcomp(&regex, final_regex_pattern, REG_EXTENDED|REG_NOSUB)) != 0) {
        regerror(reti, &regex, error, sizeof(error));
        iftError("Regex Compilation Failed: \"%s\"\n" \
                 "IFT_ERROR: %s", "iftRegexMatch", final_regex_pattern, error);
    }
    
    // Execute Regex
    reti = regexec(&regex, str, (size_t) 0, NULL, 0);
    regfree(&regex);
    
    return (reti == 0);
}

// ---------- iftRegex.c end 
// ---------- iftImageMath.c start

iftImage *iftAddValue(const iftImage *img, int val)
{
  iftImage *nimg=NULL;
  int p;

  if (val < 0)
    iftWarning("Resulting image may have negative values","iftAddValue");

  if (iftIsColorImage(img)){
    nimg   = iftCreateColorImage(img->xsize,img->ysize,img->zsize, iftImageDepth(img));
    for (p=0; p < img->n; p++){
      nimg->val[p]=img->val[p] + val;
      nimg->Cb[p] =img->Cb[p];
      nimg->Cr[p] =img->Cr[p];
    }
  } else {
    nimg = iftCreateImage(img->xsize,img->ysize,img->zsize);
    for (p=0; p < img->n; p++){
      nimg->val[p]=img->val[p] + val;
    }
  }
  iftCopyVoxelSize(img,nimg);

  return(nimg);
}

// ---------- iftImageMath.c end 
// ---------- iftNumerical.c start

iftIntArray *iftIntRange(int begin, int end, int inc) 
{
    int n = ((end - begin) / inc) + 1;
    
    iftIntArray *space = iftCreateIntArray(n);
    
    for (int i = 0; i < n; ++i) {
        space->val[i] = begin + (inc*i);
    }
    
    return space;
}

// ---------- iftNumerical.c end
// ---------- iftSort.c start

void iftFQuickSort( float *value, int *index, int i0, int i1, uchar order ) 
{
  int m, d;


	if( i0 < i1 ) {
		/* random index to avoid bad pivots.*/
		// d = iftRandomInteger( i0, i1 );
    // to guarantee the same behavior on ties 
    d = (i0 + i1) / 2;
		iftSwap( value[ d ], value[ i0 ] );
		iftSwap( index[ d ], index[ i0 ] );
		m = i0;

		if(order == IFT_INCREASING ) {
			for( d = i0 + 1; d <= i1; d++ ) {
				if( value[ d ] < value[ i0 ] ) {
					m++;
					iftSwap( value[ d ], value[ m ] );
					iftSwap( index[ d ], index[ m ] );
				}
			}
		}
		else {
			for( d = i0 + 1; d <= i1; d++ ) {
				if( value[ d ] > value[ i0 ] ) {
					m++;
					iftSwap( value[ d ], value[ m ] );
					iftSwap( index[ d ], index[ m ] );
				}
			}
		}
		iftSwap( value[ m ], value[ i0 ] );
		iftSwap( index[ m ], index[ i0 ] );
		iftFQuickSort( value, index, i0, m - 1, order );
		iftFQuickSort( value, index, m + 1, i1, order );
	}
}


void iftBucketSort(int *value, int *index, int nelems, uchar order)
{
  int i,j,maxval= IFT_INFINITY_INT_NEG, *sval=NULL, *sind=NULL;
  iftGQueue *Q=NULL;

  for(i=0; i < nelems; i++) 
    if (value[i] > maxval)
      maxval = value[i]; 

 
  Q    = iftCreateGQueue(maxval+1,nelems,value);
  sval = iftAllocIntArray(nelems);
  sind = iftAllocIntArray(nelems);

  if (order == IFT_DECREASING){
    iftSetRemovalPolicy(Q,IFT_MAXVALUE);
  }

  for(i=0; i < nelems; i++)
    iftInsertGQueue(&Q,i);

  j = 0;
  while(!iftEmptyGQueue(Q)) {
    i = iftRemoveGQueue(Q);
    sval[j]  = value[i];
    sind[j]  = index[i];
    j++;
  }

  iftDestroyGQueue(&Q);

  for(i=0; i < nelems; i++){
    value[i]=sval[i];
    index[i]=sind[i];
  }
  iftFree(sval);
  iftFree(sind);
}

// ---------- iftSort.c end
// ---------- iftBMap.c start

iftBMap *iftCreateBMap(int n) 
{
    iftBMap *b;
    b= (iftBMap *) iftAlloc(1,sizeof(iftBMap));
    b->n        = n;
    b->nbytes   = n/8;
    if (n%8) b->nbytes++;
    b->val = (char *) iftAlloc(b->nbytes,sizeof(char));
    if (b->val==NULL){
        iftError(MSG_MEMORY_ALLOC_ERROR, "iftCreateBMap");
    }
    return b;
}

void iftDestroyBMap(iftBMap **bmap) 
{
    iftBMap *aux=*bmap;
    
    if (aux != NULL) {
        iftFree(aux->val);
        iftFree(aux);
        *bmap=NULL;
    }
    
}

// ---------- iftBMap.c end
// ---------- iftList.c start 

iftList *iftCreateList() 
{
    iftList *L = (iftList *) iftAlloc(1, sizeof(iftList));
    
    // just to really force
    L->n    = 0;
    L->head = NULL;
    L->tail = NULL;
    
    return L;
}

void iftDestroyList(iftList **L) 
{
    iftList *L_aux = *L;
    
    if (L_aux != NULL) {
        iftNode *snode = L_aux->head;
        iftNode *tnode = NULL;
        
        while (snode != NULL) {
            tnode = snode;
            snode = snode->next;
            iftFree(tnode);
        }
        iftFree(L_aux);
        *L = NULL;
    }
}

bool iftIsEmptyList(const iftList *L) 
{
    return (L->n == 0);
}

void iftInsertListIntoTail(iftList *L, int elem) 
{
    if (L == NULL)
        iftError("The Integer Linked List L is NULL. Allocated it firstly", "iftInsertListIntoTail");
    
    iftNode *node  = (iftNode*) iftAlloc(1, sizeof(iftNode));
    node->elem     = elem;
    node->previous = NULL;
    node->next     = NULL;
    
    
    // The String Linked List is empty
    if (L->head == NULL) {
        L->head = node;
        L->tail = node;
        L->n    = 1;
    }
    else {
        L->tail->next  = node;
        node->previous = L->tail;
        L->tail        = node;
        L->n++;
    }
}

int iftRemoveListTail(iftList *L) 
{
    if (L == NULL)
        iftError("The Integer Doubly Linked List L is NULL. Allocated it firstly", "iftRemoveListTail");
    
    int elem      = IFT_NIL;
    iftNode *node = NULL;
    
    // if there are elements
    if (L->head != NULL) {
        node    = L->tail;
        L->tail = L->tail->previous;
        
        // checks if the list is empty now
        if (L->tail == NULL)
            L->head = NULL;
        else
            L->tail->next = NULL;
        
        L->n--;
        elem = node->elem;
        iftFree(node);
    }
    
    return elem;
}

iftIntArray *iftListToIntArray(const iftList *L) 
{
    iftIntArray *arr = iftCreateIntArray(L->n);
    
    iftNode *node = L->head;
    
    for (int i = 0; i < L->n; i++) {
        arr->val[i] = node->elem;
        node = node->next;
    }
    
    return arr;
}

// ---------- iftList.c end
// ---------- iftRegion.c start 

iftImage* iftRelabelRegions(iftImage* labelled, iftAdjRel* adj_rel)
{

	iftImage *relabelled = iftCreateImage(labelled->xsize,labelled->ysize,labelled->zsize);
	iftCopyVoxelSize(labelled,relabelled);
	iftFIFO *F = iftCreateFIFO(labelled->n);

	int nlabels = 1;

	int i;
	for(i = 0; i < labelled->n; i++){
		if( (labelled->val[i] != 0) && (relabelled->val[i] == 0)){
			relabelled->val[i] = nlabels;
			iftInsertFIFO(F,i);

			while(!iftEmptyFIFO(F)){
				int p = iftRemoveFIFO(F);
				iftVoxel u = iftGetVoxelCoord(labelled,p);

				int j;
				for(j = 1; j < adj_rel->n; j++){
					iftVoxel v = iftGetAdjacentVoxel(adj_rel,u,j);

					if(iftValidVoxel(labelled,v)){
						int q = iftGetVoxelIndex(labelled,v);

//                        if((relabelled->val[q] == 0) && (labelled->val[p] == labelled->val[q]) ){
//                            relabelled->val[q] = nlabels;
//                            iftInsertFIFO(F,q);
//                        }

						if(relabelled->val[q] == 0){
                            if(labelled->val[p] == labelled->val[q])
                                iftInsertFIFO(F,q);
                            relabelled->val[q] = nlabels;
						}
					}
				}
			}

			nlabels++;
		}
	}

	iftDestroyFIFO(&F);

	return relabelled;
}

// ---------- iftRegion.c end 
// ---------- iftDict.c start

ulong iftPowerOf2GtN(size_t n, ulong *expo) 
{
    ulong pow2 = 2;
    
    (*expo) = 1;
    while (pow2 < n) {
        pow2 = pow2 << 1; // (pow2 << 1) = pow2 * 2^(1)
        (*expo)++;
    }
    
    return pow2;
}

ulong iftFastestPrimeFromPowerOf2(ulong lower_pow_of_2, ulong higher_power_of_2) 
{
    // Exceptional case
    if ((lower_pow_of_2 == 2) && (higher_power_of_2 == 4))
        return 4;
    
    ulong mean = (lower_pow_of_2 + higher_power_of_2) /2;;
    ulong candidates[2];
    
    candidates[0] = mean - 1;
    candidates[1] = mean + 1;
    
    while((candidates[0] > lower_pow_of_2) && (candidates[1] < higher_power_of_2)) {
        if (iftIsPrime(candidates[0]))
            return(candidates[0]);
        if (iftIsPrime(candidates[1]))
            return(candidates[1]);
        candidates[0]--;
        candidates[1]++;
    }
    
    iftError("There is no prime number between %ld and %ld", "iftFastestPrimeFromPowerOf2",
             lower_pow_of_2, higher_power_of_2);
    
    return -1; // just to avoid warnings... Never will be reached
}

long iftStrToLong(const char *str) 
{
    if (str == NULL)
        iftError("String is NULL", "iftStrToLong");
    
    long lval = 0;
    
    for(size_t i = 0; i < strlen(str); i++) {
        lval += str[i];
        lval += (lval << 10);
        lval ^= (lval >> 6);
    }
    lval += (lval << 3);
    lval ^= (lval >> 11);
    lval += (lval << 15);
    
    return lval;
}

size_t iftGetHashKey(iftGVal key) 
{
    long hash_key = 0;
    
    switch(key.type) {
        case IFT_CHAR_TYPE:
            hash_key = (size_t) key.char_val;
            break;
        case IFT_UCHAR_TYPE:
            hash_key = (size_t) key.uchar_val;
            break;
        case IFT_STR_TYPE:
            hash_key = iftStrToLong(key.str_val);
            break;
        case IFT_LONG_TYPE:
            hash_key = key.long_val;
            break;
        case IFT_ULONG_TYPE:
            hash_key = (long) key.ulong_val; // it could have loss of data
            break;
        case IFT_DBL_TYPE:
            hash_key = (key.dbl_val * 10000); // truncates the abs double value to int by increasing four decimal numbers
            break;
        default:
            iftError("Invalid key datatype %s for hashing", "iftGetHashKey",
                     iftCDataTypeToString(key.type));
    }
    
    return labs(hash_key);
}

size_t iftHashing(iftGVal key, iftDict *dict) 
{
    if (dict == NULL)
        iftError("Dict is NULL", "iftHashing");
    
    size_t hash_key = iftGetHashKey(key);
    
    size_t delta = 1;
    size_t h0, h;
    h0 = h = (hash_key % dict->capacity);
    
    while ((dict->table[h].key.type != IFT_UNTYPED) && (delta <= dict->capacity)) {
        // found a table position with the same hash_key... the new value will replace this bucket
        if (iftCompareGVal(dict->table[h].key, key)) {
            return h;
        }
        else {
            delta++;
            h = (h0 + (delta - 1) * (1 + (hash_key % (dict->capacity - 2)))) % dict->capacity;
        }
    }
    
    if (delta > dict->capacity)
        return IFT_HASH_FAIL;
    
    return h; // successful
}

void iftDictResizing(iftDict **dict) 
{
    if (dict == NULL)
        iftError("Dict Reference is NULL", "iftDictResizing");
    if (*dict == NULL)
        iftError("Dict is NULL", "iftDictResizing");
    
    iftDict *aux = *dict;
    size_t n = 2*aux->capacity + 1; // new number of elements to be inserted
    
    // Rehashing the elements
    iftDict *new_dict = iftCreateDictWithApproxSize(n); // results in a dict with a better size to avoid colisions
    for (iftKeyVal k = iftDictFirst(aux); iftDictValidKey(k); k = iftDictNext(aux, k))
        iftInsertKeyValIntoDict(k.key, k.val, new_dict);
    
    // just deallocates the table of the original dict in order to keep the reference to it
    // ps: the strings (char*) does not deallocated, because the new_dict points to them
    iftFree(aux->table);
    
    aux->capacity = new_dict->capacity;
    aux->table    = new_dict->table;
    aux->firstKey = new_dict->firstKey;
    aux->lastKey = new_dict->lastKey;
    new_dict->table = NULL; // avoid that the destroyer deallocates the table now pointed by the original dict
    iftDestroyDict(&new_dict);
}

bool iftDictContainGKey(iftGVal key, const iftDict *dict, size_t *key_index) 
{
    if (dict == NULL)
        iftError("Dict is NULL", "iftDictContainGKey");
    
    bool contain_gkey = false;
    
    if (!iftIsDictEmpty(dict)) {
        // It is basically the same hashing function of iftHashing(), however,
        // it ignores empty buckets and continue seeking the key (by double hashing) until
        // all buckets have been visited.
        // This enables that a key B rehashed due to a colision with a key A is able to reach
        // its real position even if the key A have been removed.
        size_t hash_key = iftGetHashKey(key);
        size_t delta    = 1;
        size_t h0, h;
        h0 = h = (hash_key % dict->capacity);
        
        while ((delta <= dict->capacity) && !iftCompareGVal(dict->table[h].key, key)) {
            delta++;
            h = (h0 + (delta - 1) * (1 + (hash_key % (dict->capacity - 2)))) % dict->capacity;
        }
        
        if ((delta > dict->capacity) || !iftCompareGVal(dict->table[h].key, key))
            contain_gkey = false;
        else {
            if (key_index != NULL)
                *key_index = h;
            contain_gkey = true;
        }
    }
    
    // deallocates if the input key is a string
    if (key.type == IFT_STR_TYPE) {
        iftFree(key.str_val);
        key.str_val = NULL;
    }
    
    return contain_gkey;
}

iftGVal iftGetGValFromDict(iftGVal key, const iftDict *dict, iftCDataType val_type) 
{
    size_t h;
    iftGVal val = {.type=IFT_UNTYPED};
    
    // since we can have compound keys (eg: "sub-dict:sub-sub-dict:test", we need to access the final/last dict to the element
    const iftDict *final_dict = dict;
    
    // if the key is a string, it could have intermediate dicts/sub-keys (e.g: abc:def - abc is a dict)
    if (key.type == IFT_STR_TYPE) {
        iftSList *SL = iftSplitString(key.str_val, ":");
        
        int n = SL->n;
        
        // finds the final/last dict to insert the element
        
        if (n > 1) {
            // scan up to the penultimate node
            for (int i = 0; i <= (n-2); i++) {
                char *sub_key = iftRemoveSListHead(SL);
                
                if (iftDictContainKey(sub_key, final_dict, &h)) {
                    if (final_dict->table[h].val.type != IFT_DICT_TYPE)
                        iftError("Sub-key \"%s\" is not a dict", "iftGetGValFromDict", sub_key);
                    final_dict = final_dict->table[h].val.dict_val;
                }
                else
                    iftError("Sub-key \"%s\" does not have a dict element", "iftGetGValFromDict", sub_key);
                
                iftFree(sub_key);
            }
        }
        char *final_key = iftRemoveSListHead(SL);
        iftDestroySList(&SL);
        
        iftFree(key.str_val);
        key = iftCreateGVal(final_key);
        iftFree(final_key);
    }
    
    if (iftDictContainKey(key, final_dict, &h))
        val = final_dict->table[h].val;
    else
        iftError("KeyError: %s (%s) does not exist", "iftGetGValFromDict", iftGValToString(key),
                 iftCDataTypeToString(key.type));
    
    if (val.type != val_type)
        iftError("Value %s from key %s is actually %s and not %s",
                 "iftGetValFromDict", iftGValToString(val), iftGValToString(key),
                 iftCDataTypeToString(val.type), iftCDataTypeToString(val_type));
    
    // deallocates if the input key is a string
    if (key.type == IFT_STR_TYPE) {
        iftFree(key.str_val);
        key.str_val = NULL;
    }
    
    return val;
}

iftDict *iftCreateDict() 
{
    return iftCreateDictWithInitialSize(IFT_INITIAL_DICT_SIZE);
}

iftDict *iftCreateDictWithInitialSize(size_t n) 
{
    iftDict *dict = (iftDict*) iftAlloc(1, sizeof(iftDict));
    
    dict->capacity             = n;
    dict->n_busy_buckets   = 0;
    dict->table            = (iftKeyVal*) iftAlloc(dict->capacity, sizeof(iftKeyVal));
    
    for (size_t i = 0; i < dict->capacity; i++) {
        dict->table[i].key.type = IFT_UNTYPED;
        dict->table[i].val.type = IFT_UNTYPED;
    }
    
    dict->erase_elements = false;
    dict->firstKey       = dict->lastKey = IFT_NIL;
    
    return dict;
}

bool iftIsDictEmpty(const iftDict *dict) 
{
    if (dict == NULL)
        iftError("Dict is NULL", "iftIsDictEmpty");
    
    return (dict->n_busy_buckets == 0);
}

iftDict *iftCopyDict(const iftDict *dict) 
{
    iftDict *copy = NULL;
    
    if (dict != NULL) {
        copy                 = (iftDict*) iftAlloc(1, sizeof(iftDict));
        copy->capacity       = dict->capacity;
        copy->n_busy_buckets = dict->n_busy_buckets;
        copy->table          = (iftKeyVal*) iftAlloc(copy->capacity, sizeof(iftKeyVal));
        copy->firstKey       = dict->firstKey;
        copy->lastKey        = dict->lastKey;
        copy->erase_elements = dict->erase_elements;
        
        int i = dict->firstKey;
        iftKeyVal keyval = iftDictFirst(dict);
        while (iftDictValidKey(keyval)) {
            copy->table[i].key  = iftCopyGVal(keyval.key);
            copy->table[i].val  = iftCopyGVal(keyval.val);
            copy->table[i].prev = keyval.prev;
            copy->table[i].next = keyval.next;
            
            i = keyval.next;
            keyval = iftDictNext(dict, keyval);
        }
    }
    
    return copy;
}

iftKeyVal iftDictFirst(const iftDict* dict) 
{
    return dict->table[dict->firstKey];
}

iftKeyVal iftDictNext(const iftDict* dict, const iftKeyVal keyval) 
{
    if(keyval.next!=IFT_NIL)
        return dict->table[keyval.next];
    else {
        iftKeyVal keyval;
        keyval.key.type = IFT_UNTYPED;
        keyval.val.type = IFT_UNTYPED;
        keyval.prev = IFT_NIL;
        keyval.next = IFT_NIL;
        return keyval;
    }
}

bool iftDictValidKey(const iftKeyVal keyval) 
{
    return keyval.key.type != IFT_UNTYPED && keyval.val.type != IFT_UNTYPED;
}

bool iftInsertKeyValIntoDict(iftGVal key, iftGVal val, iftDict *dict) 
{
    // INSERTIONS
    if (iftIsDictFull(dict)) {
        iftDictResizing(&dict);
    }
    
    // if the key is a string, it could have intermediate dicts.sub-keys (e.g: abc:def - abc is a dict)
    if (key.type == IFT_STR_TYPE) {
        iftSList *SL = iftSplitString(key.str_val, ":");
        
        
        // string key has a sub-key, which must be an intermediate dict
        if (SL->n > 1) {
            char *sub_key = iftRemoveSListHead(SL);
            iftDestroySList(&SL);
            
            iftDict *sub_dict = NULL;
            
            // gets the sub-dict
            size_t h;
            if (iftDictContainKey(sub_key, dict, &h)) {
                if (dict->table[h].val.type == IFT_DICT_TYPE) {
                    sub_dict = dict->table[h].val.dict_val;
                }
                else {
                    sub_dict = iftCreateDict();
                    iftInsertIntoDict(sub_key, sub_dict, dict);
                }
            }
            else {
                sub_dict = iftCreateDict();
                iftInsertIntoDict(sub_key, sub_dict, dict);
            }
            
            
            char *prefix = iftConcatStrings(2, sub_key, ":");
            iftFree(sub_key);
            
            sub_key = iftRemovePrefix(key.str_val, prefix);
            iftFree(prefix);
            
            // insert the remaining sub keys
            iftFree(key.str_val);
            key = iftCreateGVal(sub_key);
            iftFree(sub_key);
            return iftInsertKeyValIntoDict(key, val, sub_dict);
        }
        else iftDestroySList(&SL);
    }
    
    size_t h = iftHashing(key, dict);
    
    if (h == IFT_HASH_FAIL)
        return false;
    
    // if the chosen bucket is empty, it increments the n_busy_buckets
    dict->n_busy_buckets += (dict->table[h].key.type == IFT_UNTYPED);
    
    // bucket is empty, update the iterator
    if (dict->table[h].key.type == IFT_UNTYPED) {
        if (dict->lastKey!=IFT_NIL)
            dict->table[dict->lastKey].next = h;
        
        dict->table[h].next = IFT_NIL;
        dict->table[h].prev = dict->lastKey;
        
        dict->lastKey = h;
        
        if (dict->firstKey == IFT_NIL)
            dict->firstKey = h;
    }
        // busy bucket - deallocate its element
    else {
        // just to avoid memory leak when overwriting existent key elements
        if (dict->table[h].key.type == IFT_STR_TYPE) {
            iftFree(dict->table[h].key.str_val);
            dict->table[h].key.str_val = NULL;
        }
        if (dict->table[h].val.type == IFT_STR_TYPE)
            iftFree(dict->table[h].val.str_val);
        else if (dict->erase_elements)
            iftFreeGVal(dict->table[h].val);
    }
    
    dict->table[h].key = key;
    dict->table[h].val = val;
    
    return true;
}

bool iftIsDictFull(const iftDict *dict) 
{
    if (dict == NULL)
        iftError("Dict is NULL", "iftIsDictFull");
    
    return (dict->n_busy_buckets == dict->capacity);
}

iftDict *iftCreateDictWithApproxSize(size_t n) 
{
    if ((long) n <= 0)
        iftError("The Approximate Dictionary Size %ld is <= 0", "iftCreateDict", (long) n);
    
    ulong expo;
    ulong lower_pow_of_2  = iftPowerOf2GtN(n, &expo);
    ulong higher_pow_of_2 = 1 << (expo+1); // 1 * 2^(expo+1)
    size_t size           = iftFastestPrimeFromPowerOf2(lower_pow_of_2, higher_pow_of_2);
    
    return iftCreateDictWithInitialSize(size);
}

void iftDestroyDict(iftDict **dict) 
{
    if (dict != NULL) {
        iftDict *aux = *dict;
        
        if (aux != NULL) {
            if (aux->table != NULL) {
                if (!iftIsDictEmpty(aux) && (aux->erase_elements)) {
                    for (iftKeyVal keyval = iftDictFirst(aux); iftDictValidKey(keyval); keyval = iftDictNext(aux, keyval)) {
                        iftFreeGVal(keyval.key);
                        if (keyval.val.type == IFT_DICT_TYPE)
                            iftDestroyDict(&keyval.val.dict_val);
                        else iftFreeGVal(keyval.val);
                    }
                }
                iftFree(aux->table);
            }
            iftFree(aux);
            *dict = NULL;
        }
    }
}

// ---------- iftDict.c end 
// ---------- iftCommandLineParser.c start

void _iftValidateCmdLineOpts(int n_opts, const iftCmdLineOpt cmd_line_opts[]) 
{
    // EXCEPTION CHECKERS
    if (n_opts < 0)
        iftError("Number of command line options %d is < 0", "_iftValidateCmdLineOpts", n_opts);
    // checks each option
    for (int i = 0; i < n_opts; i++) {
        // short option name
        if (iftCompareStrings(cmd_line_opts[i].short_name, "")) {
            if (iftCompareStrings(cmd_line_opts[i].long_name, ""))
                iftError("Short and Long names from the option %d are empty \"\"",
                         "_iftValidateCmdLineOpts", i + 1);
        }
        else {
            if (iftRegexMatch(cmd_line_opts[i].short_name, "^-[a-zA-Z]([a-zA-Z]|[0-9])?$")) {
                if (iftCompareStrings(cmd_line_opts[i].short_name, "-h"))
                    iftError("The short option name \"-h\" is not allowed.\n" \
                             "It is only used to call the Usage Printing function",
                             "_iftValidateCmdLineOpts");
            }
            else
                iftError("Invalid command line short name: %s\n" \
                          "Try: -[a-zA-Z]([a-zA-Z]|[0-9])+\n\n" \
                          "Exs: -f, -a, -D, -Y, -x1", "_iftValidateCmdLineOpts", cmd_line_opts[i].short_name);
        }
        // short option name
        if (!iftCompareStrings(cmd_line_opts[i].long_name, "")) {
            if (iftRegexMatch(cmd_line_opts[i].long_name, "^--[a-zA-Z]([a-zA-Z]|[0-9]|-|_)+$")) {
                if (iftCompareStrings(cmd_line_opts[i].long_name, "--h") ||
                    iftCompareStrings(cmd_line_opts[i].long_name, "--help"))
                    iftError("The long option names \"--h\" and \"--help\" are not allowed.\n" \
                             "The areonly used to call the Usage Printing function.",
                             "iftCreateCmdLineParser");
            }
            else
                iftError("Invalid command line long name: %s\n" \
                          "Try: --[a-zA-Z]([a-zA-Z]|-|_)+\n\n" \
                          "Exs: --input-image, --outImage, --normalize, --score_file, --BOW",
                         "iftCreateCmdLineParser", cmd_line_opts[i].long_name);
        }
        // Check the Argument Type
        if (cmd_line_opts[i].has_arg) {
            if ((cmd_line_opts[i].arg_type != IFT_LONG_TYPE) &&
                (cmd_line_opts[i].arg_type != IFT_DBL_TYPE) &&
                (cmd_line_opts[i].arg_type != IFT_STR_TYPE)) {
                iftError("Invalid Argument Type. Try:\nIFT_LONG_TYPE or IFT_DBL_TYPE or IFT_STR_TYPE",
                         "iftCreateCmdLineParser");
            }
        }
    }
}

bool _iftHasRequiredOpts(const iftCmdLineParser *parser) 
{
    for (int i = 0; i < parser->n_opts; i++)
        if (parser->opts[i].required)
            return true;
    
    return false;
}

bool _iftHasOptionalOpts(const iftCmdLineParser *parser) 
{
    for (int i = 0; i < parser->n_opts; i++)
        if (!parser->opts[i].required)
            return true;
    
    return false;
}

void _iftCheckRequiredCmdLineOpt(const iftCmdLineParser *parser, const iftDict *args) 
{
    for (int i = 0; i < parser->n_opts; i++) {
        if (parser->opts[i].required) {
            
            char *key = NULL;
            
            // it sets the key which will be used to find out the option in the dict
            if (!iftCompareStrings(parser->opts[i].short_name, ""))
                key = parser->opts[i].short_name;
            else
                key = parser->opts[i].long_name;
            
            
            // the required option is missing
            if (!iftDictContainKey(key, args, NULL)) {
                char missing_opt_names[1024]; // just used to print the missing option names
                
                if (!iftCompareStrings(parser->opts[i].short_name, "") &&
                    !iftCompareStrings(parser->opts[i].long_name, "")) {
                    sprintf(missing_opt_names, "\'%s\' or \'%s\'",
                            parser->opts[i].short_name, parser->opts[i].long_name);
                }
                else if (!iftCompareStrings(parser->opts[i].short_name, ""))
                    strcpy(missing_opt_names, parser->opts[i].short_name);
                else
                    strcpy(missing_opt_names, parser->opts[i].long_name);
                
                fprintf(stderr, "Required option %s is missing.\nRun \'%s -h\' or \'%s --help\' to see a full list of " \
                "available command line options", missing_opt_names, parser->program_name, parser->program_name);
                exit(EXIT_FAILURE);
            }
        }
    }
}

bool _iftGetCmdLineOpt(const char *cmd_line_name, iftCmdLineParser *parser, iftCmdLineOpt *out_opt) 
{
    for (int i = 0; i < parser->n_opts; i++)
        if (iftCompareStrings(cmd_line_name, parser->opts[i].short_name) ||
            iftCompareStrings(cmd_line_name, parser->opts[i].long_name)) {
            *out_opt = parser->opts[i];
            return true;
        }
    
    return false;
}

iftCmdLineParser * iftCreateCmdLineParser(const char *description, int n_opts, iftCmdLineOpt cmd_line_opts[])
{
    _iftValidateCmdLineOpts(n_opts, cmd_line_opts);
    
    iftCmdLineParser *parser = (iftCmdLineParser*) iftAlloc(1, sizeof(iftCmdLineParser));
    parser->n_opts           = n_opts;
    parser->opts             = cmd_line_opts;
    parser->description      = iftCopyString(description);
    
    return parser;

}

iftDict *iftParseCmdLine(int argc, const char *argv[], iftCmdLineParser *parser) 
{
    // EXCEPTION CHECKERS
    if (argv == NULL)
        iftError("argv is NULL", "iftParseCmdLine");
    if (parser == NULL)
        iftError("Command Line Parser is NULL", "iftParseCmdLine");
    
    iftDict *arg_dict = NULL;
    strcpy(parser->program_name, argv[0]); // copies the program name to the parser
    
    char default_error_msg[IFT_STR_DEFAULT_SIZE];
    sprintf(default_error_msg, "- Run \'%s -h\' or \'%s --help\' to see a full list of available " \
                               "command line options\n", parser->program_name, parser->program_name);
    
    
    // checks if some help option was passed
    for (int i = 1; i < argc; i++)
        if (iftCompareStrings(argv[i], "-h") || iftCompareStrings(argv[i], "--h") ||
            iftCompareStrings(argv[i], "--help")) {
            iftPrintUsage(parser);
        }
    
    
    // gets the number of required options
    // this will enables to use a program that contains only optional options
    int n_req_opts = 0;
    for (int i = 0; i < parser->n_opts; i++)
        n_req_opts += (parser->opts[i].required);
    
    
    // Parse the Command Line
    if ((argc == 1) && (n_req_opts != 0))
        iftPrintUsage(parser);
    else {
        // the num of expected elements to be inserted into dict is this due to the possibility of two option names
        arg_dict = iftCreateDict();
        
        // for each cmd line option and argument (it ignores the program name in argv[0])
        for (int i = 1; i < argc; i++) {
            iftCmdLineOpt opt;
            if (_iftGetCmdLineOpt(argv[i], parser, &opt)) {
                iftGVal arg = {.type=IFT_UNTYPED, .ptr_val=NULL};
                
                if (opt.has_arg) {
                    i++;
                    if (i >= argc) {
                        fprintf(stderr, "- Argument for the option \'%s\' is missing\n%s", argv[i-1], default_error_msg);
                        fflush(stdout);
                        fflush(stderr);
                        exit(EXIT_FAILURE);
                    }
                    else if (iftRegexMatch(argv[i], "^-[a-zA-Z]([a-zA-Z]|[0-9])?$") || iftRegexMatch(argv[i], "^--[a-zA-Z]([a-zA-Z]|[0-9]|-|_)+$")) {
                        fprintf(stderr, "- Invalid argument \'%s\' for the option \'%s\'\n%s", argv[i], argv[i-1], default_error_msg);
                        fflush(stdout);
                        fflush(stderr);
                        exit(EXIT_FAILURE);
                    }
                    
                    if (opt.arg_type == IFT_STR_TYPE) {
                        arg = iftCreateGVal(argv[i]);
                    }
                        // checks if it's a valid number
                    else if (iftRegexMatch(argv[i], "^-?[0-9]+(\\.[0-9]*)?$")) {
                        if (opt.arg_type == IFT_LONG_TYPE) {
                            arg = iftCreateGVal(atol(argv[i]));
                        }
                        else if (opt.arg_type == IFT_DBL_TYPE)
                            arg = iftCreateGVal(atof(argv[i]));
                        else
                            iftError("Invalid Datatype for the Argument of the option: %s",
                                     "iftParseCmdLine", argv[i - 1]);
                    }
                    else
                        iftError("Invalid Argument for the of the option: %s\nArgument: %s",
                                 "iftParseCmdLine", argv[i - 1], argv[i]);
                }
                
                if (!iftCompareStrings(opt.short_name, ""))
                    iftInsertKeyValIntoDict(iftCreateGVal(opt.short_name), iftCopyGVal(arg), arg_dict);
                if (!iftCompareStrings(opt.long_name, ""))
                    iftInsertKeyValIntoDict(iftCreateGVal(opt.long_name), iftCopyGVal(arg), arg_dict);
                
                if (arg.type == IFT_STR_TYPE)
                    iftFree(arg.str_val);
            }
            else {
                fprintf(stderr, "- Unknown option: %s\n%s", argv[i], default_error_msg);
                fflush(stdout);
                fflush(stderr);
                exit(EXIT_FAILURE);
            }
        }
        
        // Checks if all required/mandatory option were passed
        _iftCheckRequiredCmdLineOpt(parser, arg_dict);
    }
    
    return arg_dict;
}

void iftDestroyCmdLineParser(iftCmdLineParser **parser) 
{
    iftCmdLineParser *aux = *parser;
    
    if (aux != NULL) {
        iftFree(aux->description);
        iftFree(aux);
        *parser = NULL;
    }
}

void iftPrintUsage(const iftCmdLineParser *parser) 
{
    if (parser == NULL)
        iftError("Parser is NULL", "iftPrintUsage");
    
    int length_for_printing_opts = 0; // used to align the option descriptions (helps)
    for (int i = 0; i < parser->n_opts; i++) {
        int length = strlen(parser->opts[i].short_name) + strlen(parser->opts[i].long_name);
        
        if (length_for_printing_opts < length)
            length_for_printing_opts = length;
    }
    length_for_printing_opts += 8; // first 2 spaces, 1 comma between short and long names, 4 spaces after long name
    
    ////////////////////// PRINTING //////////////////////
    char opt_names[1024];
    int n_spaces = 0;
    fprintf(stdout, "--------------\n%s\n\n", parser->description);
    fprintf(stdout, "Usage\n  %s [OPTIONS]\n\n", parser->program_name);
    
    // Prints the Required Options
    if (_iftHasRequiredOpts(parser)) {
        fprintf(stdout, "Required Options\n");
        
        for (int i = 0; i < parser->n_opts; i++) {
            if (parser->opts[i].required) {
                if (!iftCompareStrings(parser->opts[i].short_name, "")) {
                    sprintf(opt_names, "  %s", parser->opts[i].short_name);
                    
                    if (!iftCompareStrings(parser->opts[i].long_name, ""))
                        sprintf(opt_names, "%s, %s", opt_names, parser->opts[i].long_name);
                }
                else sprintf(opt_names, "  %s", parser->opts[i].long_name);
                
                n_spaces = length_for_printing_opts-strlen(opt_names);
                sprintf(opt_names, "%s%*s", opt_names, n_spaces, ""); // aligning the option descriptions (helps)
                fprintf(stdout, "%s", opt_names);
                if (parser->opts[i].has_arg)
                    fprintf(stdout, "[HAS ARG] %s\n", parser->opts[i].help);
                else
                    fprintf(stdout, "%s\n", parser->opts[i].help);
            }
        }
        puts("");
    }
    
    // Prints the Optional Options
    if (_iftHasOptionalOpts(parser)) {
        fprintf(stdout, "Optional\n");
        
        for (int i = 0; i < parser->n_opts; i++) {
            if (!parser->opts[i].required) {
                if (!iftCompareStrings(parser->opts[i].short_name, "")) {
                    sprintf(opt_names, "  %s", parser->opts[i].short_name);
                    
                    if (!iftCompareStrings(parser->opts[i].long_name, ""))
                        sprintf(opt_names, "%s, %s", opt_names, parser->opts[i].long_name);
                }
                else sprintf(opt_names, "  %s", parser->opts[i].long_name);
                
                n_spaces = length_for_printing_opts-strlen(opt_names);
                sprintf(opt_names, "%s%*s", opt_names, n_spaces, ""); // aligning the option descriptions (helps)
                fprintf(stdout, "%s", opt_names);
                if (parser->opts[i].has_arg)
                    fprintf(stdout, "[HAS ARG] %s\n", parser->opts[i].help);
                else
                    fprintf(stdout, "%s\n", parser->opts[i].help);
            }
        }
        puts("");
    }
    
    // help options
    fprintf(stdout, "Help Options\n");
    sprintf(opt_names, "  -h, --help");
    n_spaces = length_for_printing_opts-strlen(opt_names);
    sprintf(opt_names, "%s%*s", opt_names, n_spaces, ""); // aligning the option descriptions (helps)
    fprintf(stdout, "%s", opt_names);
    fprintf(stdout, "%s", "Show help options\n\n");
    
    exit(EXIT_FAILURE);
    ///////////////////////////////////////////////////////
}

// ---------- iftCommandLineParser.c end 

