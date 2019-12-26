/*------------------------------------------------------------------------
 *  Copyright 2007-2010 (c) Jeff Brown <spadix@users.sourceforge.net>
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

#include <config.h>
#include <string.h>     /* memmove */

#include <zbar.h>

#ifdef DEBUG_CYCLIC
# define DEBUG_LEVEL (DEBUG_CYCLIC)
#endif
#include "debug.h"
#include "decoder.h"

#define MinRepeatingRequired 4

//#define USE_SINGLE_TREE
//#define USE_SINGLE_ELEMENT_WIDTH

//#define TestCyclic

//static int g_mallocedNodesCount = 0;

typedef struct {
    const char* name;
    int16_t elementSequence[12];
} CyclicCode;

static CyclicCode Codes[] = {
    {"DK", {1,1,1,1,2,1,2,1,1,1,1,2}},//DK-D6=2
    {"DQ", {1,1,1,1,1,2,1,2,1,1,2,2}},//H4-DQ=1
    {"DJ", {1,1,1,1,1,2,2,1,1,2,1,2}},//D7-DJ=5
    {"DX", {1,1,1,1,2,1,1,2,2,1,1,2}},//D4-DX=1
    {"D9", {1,1,1,1,2,1,2,2,1,1,1,2}},//H5-D9=1
    {"D8", {1,1,1,1,2,2,2,1,1,1,1,2}},//D8-S3=1
    {"D7", {1,1,1,1,2,1,1,2,2,2,1,2}},//D7-DJ=5
    {"D6", {1,1,1,1,2,1,2,2,1,1,2,2}},//DK-D6=2
    {"D5", {1,1,1,1,2,2,1,1,2,1,2,2}},//HA-D5=3
    {"D4", {1,1,1,1,2,2,1,2,2,1,1,2}},//D4-DX=1
    {"D3", {1,1,1,1,2,2,2,2,1,1,1,2}},//H2-D3=1
    {"D2", {1,1,1,1,2,1,1,1,2,1,2,2}},//D2-S7=1
    {"DA", {1,1,1,1,1,2,2,1,1,1,1,2}},//DA-S2=1

    {"SK", {1,1,1,1,2,1,1,1,1,2,2,1}},//S4-SK=3
    {"SQ", {1,1,1,1,1,2,1,1,1,1,1,2}},//SA-SQ=1
    {"SJ", {1,1,1,1,1,2,1,2,1,2,1,2}},//JS-SJ=1
    {"SX", {1,1,1,1,2,1,1,1,2,2,1,2}},//HK-SX=1
    {"S9", {1,1,1,1,2,1,2,1,1,1,2,2}},//H6-S9=1
    {"S8", {1,1,1,1,2,2,1,1,1,2,1,2}},//H8-S8=1
    {"S7", {1,1,1,1,2,1,1,1,2,2,2,2}},//D2-S7=1
    {"S6", {1,1,1,1,2,1,2,1,1,2,2,2}},//H7-S6=2
    {"S5", {1,1,1,1,2,1,2,1,2,2,1,2}},//C6-S5=1
    {"S4", {1,1,1,1,2,2,2,2,1,2,2,1}},//S4-SK=3
    {"S3", {1,1,1,1,2,2,2,1,1,1,2,2}},//D8-S3=1
    {"S2", {1,1,1,1,1,2,2,1,2,1,1,2}},//DA-S2=1
    {"SA", {1,1,1,1,1,2,1,1,1,2,1,2}},//SA-SQ=1

    {"HK", {1,1,1,1,2,1,1,1,2,1,1,2}},//HK-SX=1
    {"HQ", {1,1,1,1,1,2,1,1,2,2,1,1}},//CQ-HQ=1
    {"HJ", {1,1,1,1,1,2,1,2,2,1,1,2}},//JB-HJ=1
    {"HX", {1,1,1,1,2,1,1,2,1,1,2,2}},//C7-HX=1
    {"H9", {1,1,1,1,2,1,2,1,1,2,1,2}},//H3-H9=1
    {"H8", {1,1,1,1,2,1,1,1,1,2,1,2}},//H8-S8=1
    {"H7", {1,1,1,1,2,1,1,2,1,2,2,2}},//H7-S6=2
    {"H6", {1,1,1,1,2,1,2,1,2,1,2,2}},//H6-S9=1
    {"H5", {1,1,1,1,2,1,2,2,2,1,1,2}},//H5-D9=1
    {"H4", {1,1,1,1,2,2,1,2,1,1,2,2}},//H4-DQ=1
    {"H3", {1,1,1,1,2,2,2,1,1,2,1,2}},//H3-H9=1
    {"H2", {1,1,1,1,1,2,2,2,1,1,1,1}},//H2-D3=1
    {"HA", {1,1,1,1,1,2,1,1,2,1,1,1}},//HA-D5=3

    {"CK", {1,1,1,1,2,1,1,2,1,1,1,2}},//CJ-CK=1
    {"CQ", {1,1,1,1,1,2,1,1,2,2,1,2}},//CQ-HQ=1
    {"CJ", {1,1,1,1,1,1,1,2,1,1,1,2}},//CJ-CK=1
    {"CX", {1,1,1,1,2,1,1,2,1,2,1,2}},//C4-CX=1
    {"C9", {1,1,1,1,2,1,2,1,2,1,1,2}},//C3-C9=1
    {"C8", {1,1,1,1,2,2,1,2,1,1,1,2}},//CA-C8=1
    {"C7", {1,1,1,1,2,1,1,2,2,1,2,2}},//C7-HX=1
    {"C6", {1,1,1,1,2,1,2,1,2,1,1,1}},//C6-S5=1
    {"C5", {1,1,1,1,2,2,1,1,1,2,2,2}},//C2-C5=1
    {"C4", {1,1,1,1,2,2,1,2,1,2,1,2}},//C4-CX=1
    {"C3", {1,1,1,1,2,2,2,1,2,1,1,2}},//C3-C9=1
    {"C2", {1,1,1,1,2,1,1,1,1,2,2,2}},//C2-C5=1
    {"CA", {1,1,1,1,1,2,1,2,1,1,1,2}},//CA-C8=1

    {"JB", {1,1,1,1,1,2,2,2,2,1,1,2}},//JB-HJ=1
    {"JS", {1,1,1,1,1,2,2,2,1,2,1,2}},//JS-SJ=1
    {"AD", {1,1,1,1,2,1,2,1,1,2,2,1}},
};
#ifdef TestCyclic
static int16_t TestPairWidths[] = {
//    0,0,0,0,1,1,0,0,1,1,1,
//    0,0,0,0,1,2,1,1,1,0,1,
//    0,0,0,1,1,1,2,1,1,1,1,

    0,0,0,0,1,1,0,0,0,0,1,2,1,1,1,0,0,0,0,1,1,0,0,1,1,1,
    0,0,0,1,1,1,2,1,1,1,1,
    
//    0,1,1,1,
//    1,1,0,1,
//    1,1,1,1,
};
#endif //#ifdef TestCyclic
void minHammingDistance() {
    int codesCount = sizeof(Codes) / sizeof(Codes[0]);
    int codeLength = sizeof(Codes[0].elementSequence) / sizeof(Codes[0].elementSequence[0]);
    int minHD = codeLength, minPairI = 0, minPairJ = 0;
    for (int i = codesCount - 1; i > 0; --i)
    {
        for (int j = i - 1; j >= 0; --j)
        {
            int HD = 0;
            for (int k = codeLength - 1; k >= 0; --k)
            {
                if (Codes[i].elementSequence[k] != Codes[j].elementSequence[k])
                    HD++;
            }
            if (HD < minHD)
            {
                minHD = HD;
                minPairI = i;
                minPairJ = j;
            }
        }
    }
    printf("#Barcodes# Min hamming distance is %d of '%s' and '%s'\n", minHD, Codes[minPairI].name, Codes[minPairJ].name);
}

void CyclicCharacterTreeAdd(CyclicCharacterTreeNode* root, int16_t leafValue, int16_t* path, int length) {
    if (!root) return;
    
    if (0 == length)
    {//printf("#Barcodes# leafValue=%d\n", leafValue);
        if (-1 != root->leafValue)
        {
            printf("#Barcodes# Collides between '%s' and '%s'\n", Codes[root->leafValue].name, Codes[leafValue].name);
        }
        root->leafValue = leafValue;
        return;
    }
#ifdef USE_SINGLE_ELEMENT_WIDTH
    int16_t c = path[0] - 1;
    if (c < 0 || c > 1) return;
#else
    int16_t c = path[0] + path[1] - 2;
    if (c < 0 || c > 2) return;
#endif
    CyclicCharacterTreeNode* child = root->children[c];
    if (!child)
    {
        child = CyclicCharacterTreeNodeCreate();
        //printf("#Cyclic# Add child as #%d\n", c);
        root->children[c] = child;
    }
    
    CyclicCharacterTreeAdd(child, leafValue, path + 1, length - 1);
}

CyclicCharacterTreeNode* CyclicCharacterTreeNodeCreate() {
    CyclicCharacterTreeNode* ret = (CyclicCharacterTreeNode*) malloc(sizeof(CyclicCharacterTreeNode));
    //printf("#Cyclic# New node: %d\n", ++g_mallocedNodesCount);
    CyclicCharacterTreeNodeReset(ret);
    return ret;
}

CyclicCharacterTreeNode* CyclicCharacterTreeNodeNext(const CyclicCharacterTreeNode* current, int16_t c) {
    if (!current) return NULL;
    
    if (c < 0 || c > 2) return NULL;
    
    return current->children[c];
}

void cyclic_destroy (cyclic_decoder_t *decoder)
{//TODO:
//    dcode128->direction = 0;
//    dcode128->element = 0;
//    dcode128->character = -1;
//    dcode128->s6 = 0;
    if (decoder->charTrees)
    {
        for (int i = decoder->maxS12OfChar - decoder->minS12OfChar; i >= 0; --i)
        {
            CyclicCharacterTreeNode* head = CyclicCharacterTreeNodeCreate();
            CyclicCharacterTreeNode* tail = head;
            head->children[0] = decoder->charTrees[i]; // children[0] as value, children[1] as next
            while (head)
            {
                for (int i = 0; i < 3; ++i)
                {
                    CyclicCharacterTreeNode* parent = head->children[0];
                    if (!parent->children[i]) continue;
                    
                    tail->children[1] = CyclicCharacterTreeNodeCreate();
                    tail->children[1]->children[0] = parent->children[i];
                    tail = tail->children[1];
                }
                
                free(head->children[0]);
                //printf("#Cyclic# Delete node: %d\n", --g_mallocedNodesCount);
                CyclicCharacterTreeNode* next = head->children[1];
                free(head);
                //printf("#Cyclic# Delete node: %d\n", --g_mallocedNodesCount);
                head = next;
            }
        }
        free(decoder->charTrees);
        decoder->charTrees = NULL;
    }
    
    if (decoder->charSeekers)
    {
        for (int i = decoder->charSeekersCount - 1; i >= 0; --i)
        {
            free(decoder->charSeekers[i]);
        }
        free(decoder->charSeekers);
        free(decoder->candidates);
        free(decoder->repeatingCounts);
        free(decoder->nonRepeatingSpans);
        decoder->charSeekers = NULL;
    }
    if (decoder->s12OfChars)
    {
        free(decoder->s12OfChars);
        decoder->s12OfChars = NULL;
    }
}

void cyclic_reset (cyclic_decoder_t *decoder)
{//TODO:
//    dcode128->direction = 0;
//    dcode128->element = 0;
//    dcode128->character = -1;
//    dcode128->s6 = 0;
    cyclic_destroy(decoder);
//    minHammingDistance();///!!!For Debug
    const int CodesCount = sizeof(Codes) / sizeof(Codes[0]);
    decoder->s12 = 0;
//    CyclicCharacterTreeNode*** charSeekers;//One group for each elements-of-character number
//    int16_t maxCharacterLength;
//    int16_t characterPhase;// This means sum of 2 elements - 2
//    int16_t* s12OfChar;
//    int16_t minS12OfChar;
//    int16_t maxS12OfChar;
    decoder->charSeekersCount = 0;
    decoder->s12OfChars = (int16_t*) malloc(sizeof(int16_t) * CodesCount);
    decoder->minS12OfChar = 127;
    decoder->maxS12OfChar = 0;
    for (int i = CodesCount - 1; i >= 0; --i)
    {
        int16_t* seq = Codes[i].elementSequence;
        int length = sizeof(Codes[i].elementSequence) / sizeof(Codes[i].elementSequence[0]);
        if (length > decoder->charSeekersCount)
        {
            decoder->charSeekersCount = length;
        }
        
        int16_t s12OfChar = 0;
        for (int j = length - 1; j >= 0; --j)
        {
            s12OfChar += seq[j];
        }
        decoder->s12OfChars[i] = s12OfChar;
        if (s12OfChar > decoder->maxS12OfChar)
        {
            decoder->maxS12OfChar = s12OfChar;
        }
        if (s12OfChar < decoder->minS12OfChar)
        {
            decoder->minS12OfChar = s12OfChar;
        }
    }
#ifndef USE_SINGLE_ELEMENT_WIDTH
    decoder->charSeekersCount--;
#endif
    const int uniqueS12Count = decoder->maxS12OfChar - decoder->minS12OfChar + 1;
//    printf("#Barcodes# minS12OfChar:%d, maxS12OfChar:%d, charSeekersCount:%d\n", decoder->minS12OfChar, decoder->maxS12OfChar, decoder->charSeekersCount);
    decoder->charSeekers = (CyclicCharacterTreeNode***) malloc(sizeof(CyclicCharacterTreeNode**) * decoder->charSeekersCount);
    decoder->candidates = (int16_t*) malloc(sizeof(int16_t) * decoder->charSeekersCount);
    decoder->repeatingCounts = (int16_t*) malloc(sizeof(int16_t) * decoder->charSeekersCount);
    decoder->nonRepeatingSpans = (int16_t*) malloc(sizeof(int16_t) * decoder->charSeekersCount);
    for (int i = decoder->charSeekersCount - 1; i >= 0; --i)
    {
        decoder->charSeekers[i] = (CyclicCharacterTreeNode**) malloc(sizeof(CyclicCharacterTreeNode*) * uniqueS12Count);
        memset(decoder->charSeekers[i], 0, sizeof(CyclicCharacterTreeNode*) * uniqueS12Count);
        decoder->candidates[i] = -1;
        decoder->repeatingCounts[i] = 0;
        decoder->nonRepeatingSpans[i] = 0;
    }
    decoder->characterPhase = 0;

    decoder->charTrees = (CyclicCharacterTreeNode**) malloc(sizeof(CyclicCharacterTreeNode*) * uniqueS12Count);
    for (int i = uniqueS12Count - 1; i >= 0; --i)
    {
        decoder->charTrees[i] = CyclicCharacterTreeNodeCreate();
    }
    for (int i = CodesCount - 1; i >= 0; --i)
    {
        int16_t* seq = Codes[i].elementSequence;
        int length = sizeof(Codes[i].elementSequence) / sizeof(Codes[i].elementSequence[0]);
#ifdef USE_SINGLE_ELEMENT_WIDTH

#ifdef USE_SINGLE_TREE
        CyclicCharacterTreeAdd(decoder->charTrees[0], i, seq, length);
#else //#ifdef USE_SINGLE_TREE
        CyclicCharacterTreeAdd(decoder->charTrees[decoder->s12OfChars[i] - decoder->minS12OfChar], i, seq, length);
#endif //#ifdef USE_SINGLE_TREE

#else //#ifdef USE_SINGLE_ELEMENT_WIDTH

#ifdef USE_SINGLE_TREE
        CyclicCharacterTreeAdd(decoder->charTrees[0], i, seq, length - 1);
#else //#ifdef USE_SINGLE_TREE
        CyclicCharacterTreeAdd(decoder->charTrees[decoder->s12OfChars[i] - decoder->minS12OfChar], i, seq, length - 1);
#endif //#ifdef USE_SINGLE_TREE

#endif //#ifdef USE_SINGLE_ELEMENT_WIDTH
    }
#ifdef TestCyclic
    ///!!!For Test:
    for (int i = 0; i < sizeof(TestPairWidths) / sizeof(TestPairWidths[0]); ++i)
    {
        cyclic_feed_element(decoder, TestPairWidths[i]);
    }
    ///!!!:For Test
#endif
}

zbar_symbol_type_t _zbar_decode_cyclic (zbar_decoder_t *dcode)
{
    cyclic_decoder_t* decoder = &dcode->cyclic;
    decoder->s12 -= get_width(dcode, 12);
    decoder->s12 += get_width(dcode, 0);
#ifdef USE_SINGLE_ELEMENT_WIDTH
    int16_t pairWidth = get_width(dcode, 0);
#else
    int16_t pairWidth = pair_width(dcode, 0);
#endif
    int16_t c = -1;
    for (int iPhase = decoder->charSeekersCount - 1; iPhase >= 0; --iPhase)
    {
        CyclicCharacterTreeNode** charSeekers = decoder->charSeekers[iPhase];
        for (int iS12OfChar = decoder->maxS12OfChar - decoder->minS12OfChar;
                 iS12OfChar >= 0; --iS12OfChar)
        {
    #ifdef TestCyclic
            int16_t e = pairWidth;///!!!For Test
    #else //#ifdef TestCyclic
            int16_t s12OfChar = decoder->minS12OfChar + iS12OfChar;
    #ifdef USE_SINGLE_ELEMENT_WIDTH
            int e = decode_e(pairWidth, decoder->s12, s12OfChar) + 1;
    #else //#ifdef USE_SINGLE_ELEMENT_WIDTH
            int e = decode_e(pairWidth, decoder->s12, s12OfChar);
    #endif //#ifdef USE_SINGLE_ELEMENT_WIDTH
    //        printf("#Barcodes# e=%d. pairWidth=%d, s12=%d, n=%d\n", e, pairWidth, decoder->s12, s12OfChar);
    #endif //#ifdef TestCyclic
    #ifdef USE_SINGLE_ELEMENT_WIDTH
            if (e < 0 || e > 1)
    #else
            if (e < 0 || e > 2)
    #endif
            {
                charSeekers[iS12OfChar] = NULL;
            }
            else if (!charSeekers[iS12OfChar])
            {
                if (decoder->characterPhase == iPhase)
                {
    #ifdef USE_SINGLE_TREE
                    charSeekers[iS12OfChar] = decoder->charTrees[0]->children[e];
    #else
                    charSeekers[iS12OfChar] = decoder->charTrees[iS12OfChar]->children[e];
    #endif
                }
            }
            else
            {
                charSeekers[iS12OfChar] = charSeekers[iS12OfChar]->children[e];
                if (charSeekers[iS12OfChar] && charSeekers[iS12OfChar]->leafValue > -1)
                {
                    c = charSeekers[iS12OfChar]->leafValue;//TODO: One c for each S12
                    printf("#Barcodes# A character found: %s, s12=%d; n=%d,i=%d\n", Codes[charSeekers[iS12OfChar]->leafValue].name, decoder->s12OfChars[charSeekers[iS12OfChar]->leafValue], s12OfChar, iPhase);
                    charSeekers[iS12OfChar] = NULL;
                }
            }
        }
        
        if (c > -1)
        {
            if (c == decoder->candidates[iPhase])
            {
                printf("#Barcodes# %d th(nd) time found '%s'\n", decoder->repeatingCounts[iPhase] + 1, Codes[c].name);
                decoder->nonRepeatingSpans[iPhase] = 0;
                if (++decoder->repeatingCounts[iPhase] == MinRepeatingRequired)
                {
                    release_lock(dcode, ZBAR_CYCLIC);
                    
                    decoder->candidates[iPhase] = -1;
                    decoder->repeatingCounts[iPhase] = 0;

                    int length = (int)(strlen(Codes[c].name) + 1);
                    size_buf(dcode, length);
                    memcpy(dcode->buf, Codes[c].name, length);
                    dcode->buflen = length;
                    printf("#Barcodes# Confirm '%s'\n", Codes[c].name);
                    return(ZBAR_CYCLIC);
                }
            }
            else
            {
                decoder->repeatingCounts[iPhase] = 1;
                decoder->candidates[iPhase] = c;
                decoder->nonRepeatingSpans[iPhase] = 0;
                printf("#Barcodes# First time found '%s'\n", Codes[c].name);
                acquire_lock(dcode, ZBAR_CYCLIC);
                
                return(ZBAR_PARTIAL);
            }
        }
        else if (decoder->candidates[iPhase] > -1)
        {
            printf("#Barcodes# Break '%s', nonRepeatingSpan=%d\n", Codes[decoder->candidates[iPhase]].name, decoder->nonRepeatingSpans[iPhase] + 1);
            if (++decoder->nonRepeatingSpans[iPhase] > sizeof(Codes[decoder->candidates[iPhase]].elementSequence) / sizeof(Codes[decoder->candidates[iPhase]].elementSequence[0]))
            {
                printf("#Barcodes# Abort '%s'\n", Codes[decoder->candidates[iPhase]].name);
                decoder->candidates[iPhase] = -1;
                decoder->repeatingCounts[iPhase] = 0;
                decoder->nonRepeatingSpans[iPhase] = 0;
                release_lock(dcode, ZBAR_CYCLIC);
            }
        }
    }
    
    if (++decoder->characterPhase == decoder->charSeekersCount)
    {
        decoder->characterPhase = 0;
    }

    return(ZBAR_NONE);
}
