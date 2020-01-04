/*------------------------------------------------------------------------
 *  Copyright 2007-2009 (c) Jeff Brown <spadix@users.sourceforge.net>
 *
 *  This file is part of the ZBar Bar Code Reader.
 *
 *  The ZBar Bar Code Reader is free software; you can redistribute it
 *  and/or modify it under the terms of the GNU Lesser Public License as
 *  published by the Free Software Foundation; either version 2.1 of
 *  the License, or (at your option) any later version.
 *
 *  The ZBar Bar Code Reader is distributed in the hope that it will be
 *  useful, but WITHOUT ANY WARRANTY; without even the implied warranty
 *  of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser Public License
 *  along with the ZBar Bar Code Reader; if not, write to the Free
 *  Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 *  Boston, MA  02110-1301  USA
 *
 *  http://sourceforge.net/projects/zbar
 *------------------------------------------------------------------------*/
#ifndef _CYCLIC_H_
#define _CYCLIC_H_

#include <string.h>

#define CodeElementLength 12

typedef struct {
    const char* name;
    int16_t elementSequence[CodeElementLength];
} CyclicCode;

extern CyclicCode CyclicCodes[];

//#define CyclicCodesCount  (sizeof(CyclicCodes) / sizeof(CyclicCodes[0]))
#define CyclicCodesCount  55

typedef struct CyclicCharacterTreeNode_s {
    struct CyclicCharacterTreeNode_s* children[3];
    int16_t leafValue;
} CyclicCharacterTreeNode;

static inline void CyclicCharacterTreeNodeReset(CyclicCharacterTreeNode* node) {
    memset(node, 0, sizeof(CyclicCharacterTreeNode));
    node->leafValue = -1;
}

CyclicCharacterTreeNode* CyclicCharacterTreeNodeCreate();

typedef enum {
    CyclicTrackerPossible = 1,
    CyclicTrackerConfirmed = 2,
    CyclicTrackerUncertain = 0,
    CyclicTrackerFailed = -1,
} CyclicTrackerResult;

typedef struct CodeTracker_s {
    int16_t candidate;
    int16_t fedElementsCount;
    struct CodeTracker_s* next;
//    int16_t window[CodeElementLength];
//    int16_t startIndex;
    float probabilities[CyclicCodesCount];
    float e[CodeElementLength - 1];///!!!For Debug
} CodeTracker;

/* Cyclic specific decode state */
typedef struct cyclic_decoder_s {
    CyclicCharacterTreeNode** codeTreeRoots;
//    CyclicCharacterTreeNode*** charSeekers;//One group for each elements-of-character number
    int16_t maxCodeLength;
//    int16_t characterPhase;// This means sum of 2 elements - 2
    int16_t* s12OfChars;
    int16_t minS12OfChar;
    int16_t maxS12OfChar;
    
    unsigned s12;                /* character width */
    
//    int16_t** candidates;
//    int16_t** repeatingCounts;
    CodeTracker codeTracker;

    unsigned config;
    int configs[NUM_CFGS];      /* int valued configurations */
} cyclic_decoder_t;

/* reset Cyclic specific state */
void cyclic_reset (cyclic_decoder_t *dcodeCyclic);

void cyclic_destroy (cyclic_decoder_t *dcodeCyclic);

/* decode Code 128 symbols */
zbar_symbol_type_t _zbar_decode_cyclic(zbar_decoder_t *dcode);

#endif
