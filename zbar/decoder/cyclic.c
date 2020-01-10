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
#include <stdbool.h>

#include <zbar.h>

#ifdef DEBUG_CYCLIC
# define DEBUG_LEVEL (DEBUG_CYCLIC)
#endif
#include "debug.h"
#include "decoder.h"

#define MinRepeatingRequired 1

#define MaxSpaceBetweenPeriods 2

//#define USE_SINGLE_TREE
//#define USE_SINGLE_ELEMENT_WIDTH

//#define TestCyclic

//static int g_mallocedNodesCount = 0;

CyclicCode CyclicCodes[] = {
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

    {"SK", {1,1,1,1,2,1,1,1,1,2,1,2}},//S4-SK=3
    {"SQ", {1,1,1,1,1,2,1,1,1,2,2,2}},//SA-SQ=1
    {"SJ", {1,1,1,1,1,2,1,2,1,2,1,2}},//JS-SJ=1
    {"SX", {1,1,1,1,2,1,1,1,2,2,1,2}},//HK-SX=1
    {"S9", {1,1,1,1,2,1,2,1,1,1,2,2}},//H6-S9=1
    {"S8", {1,1,1,1,2,2,1,1,1,2,1,2}},//H8-S8=1
    {"S7", {1,1,1,1,2,1,1,1,2,2,2,2}},//D2-S7=1
    {"S6", {1,1,1,1,2,1,2,1,1,2,2,2}},//H7-S6=2
    {"S5", {1,1,1,1,2,1,2,2,1,2,1,2}},//C6-S5=1
    {"S4", {1,1,1,1,2,2,1,1,2,2,1,2}},//S4-SK=3
    {"S3", {1,1,1,1,2,2,2,1,1,1,2,2}},//D8-S3=1
    {"S2", {1,1,1,1,1,2,2,1,2,1,1,2}},//DA-S2=1
    {"SA", {1,1,1,1,1,2,1,1,1,2,1,2}},//SA-SQ=1

    {"HK", {1,1,1,1,2,1,1,1,2,1,1,2}},//HK-SX=1
    {"HQ", {1,1,1,1,1,2,1,1,2,1,2,2}},//CQ-HQ=1
    {"HJ", {1,1,1,1,1,2,1,2,2,1,1,2}},//JB-HJ=1
    {"HX", {1,1,1,1,2,1,1,2,1,1,2,2}},//C7-HX=1
    {"H9", {1,1,1,1,2,1,2,1,1,2,1,2}},//H3-H9=1
    {"H8", {1,1,1,1,2,2,1,1,2,1,1,2}},//H8-S8=1
    {"H7", {1,1,1,1,2,1,1,2,1,2,2,2}},//H7-S6=2
    {"H6", {1,1,1,1,2,1,2,1,2,1,2,2}},//H6-S9=1
    {"H5", {1,1,1,1,2,1,2,2,2,1,1,2}},//H5-D9=1
    {"H4", {1,1,1,1,2,2,1,2,1,1,2,2}},//H4-DQ=1
    {"H3", {1,1,1,1,2,2,2,1,1,2,1,2}},//H3-H9=1
    {"H2", {1,1,1,1,1,2,2,2,1,1,1,2}},//H2-D3=1
    {"HA", {1,1,1,1,1,2,1,1,2,1,1,2}},//HA-D5=3

    {"CK", {1,1,1,1,2,1,1,2,1,1,1,2}},//CJ-CK=1
    {"CQ", {1,1,1,1,1,2,1,1,2,2,1,2}},//CQ-HQ=1
    {"CJ", {1,1,1,1,1,2,2,1,1,1,2,2}},//CJ-CK=1
    {"CX", {1,1,1,1,2,1,1,2,1,2,1,2}},//C4-CX=1
    {"C9", {1,1,1,1,2,1,2,1,2,1,1,2}},//C3-C9=1
    {"C8", {1,1,1,1,2,2,1,2,1,1,1,2}},//CA-C8=1
    {"C7", {1,1,1,1,2,1,1,2,2,1,2,2}},//C7-HX=1
    {"C6", {1,1,1,1,2,1,2,1,2,2,1,2}},//C6-S5=1
    {"C5", {1,1,1,1,2,2,1,1,1,2,2,2}},//C2-C5=1
    {"C4", {1,1,1,1,2,2,1,2,1,2,1,2}},//C4-CX=1
    {"C3", {1,1,1,1,2,2,2,1,2,1,1,2}},//C3-C9=1
    {"C2", {1,1,1,1,2,1,1,1,1,2,2,2}},//C2-C5=1
    {"CA", {1,1,1,1,1,2,1,2,1,1,1,2}},//CA-C8=1

    {"JB", {1,1,1,1,1,2,2,2,2,1,1,2}},//JB-HJ=1
    {"JS", {1,1,1,1,1,2,2,2,1,2,1,2}},//JS-SJ=1
    {"AD", {1,1,1,1,1,2,2,2,1,1,2,2}},
};

//#ifdef TestCyclic
//static int16_t TestPairWidths[] = {
////    0,0,0,0,1,1,0,0,1,1,1,
////    0,0,0,0,1,2,1,1,1,0,1,
////    0,0,0,1,1,1,2,1,1,1,1,
//
//    0,0,0,0,1,1,0,0,0,0,1,2,1,1,1,0,0,0,0,1,1,0,0,1,1,1,
//    0,0,0,1,1,1,2,1,1,1,1,
//
////    0,1,1,1,
////    1,1,0,1,
////    1,1,1,1,
//};
//#endif //#ifdef TestCyclic
void minHammingDistance() {
    int codesCount = sizeof(CyclicCodes) / sizeof(CyclicCodes[0]);
    int codeLength = sizeof(CyclicCodes[0].elementSequence) / sizeof(CyclicCodes[0].elementSequence[0]);
    int minHD = codeLength, minPairI = 0, minPairJ = 0;
    for (int i = codesCount - 1; i > 0; --i)
    {
        for (int j = i - 1; j >= 0; --j)
        {
            int HD = 0;
            for (int k = codeLength - 1; k >= 0; --k)
            {
                if (CyclicCodes[i].elementSequence[k] != CyclicCodes[j].elementSequence[k])
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
    dbprintf(0, "#Barcodes# Min hamming distance is %d of '%s' and '%s'\n", minHD, CyclicCodes[minPairI].name, CyclicCodes[minPairJ].name);
}

void CyclicCharacterTree0Add(CyclicCharacterTreeNode* root, int16_t leafValue, int16_t* path, int length) {
    if (!root) return;
    
    if (0 == length)
    {//dbprintf(0, "#Barcodes# leafValue=%d\n", leafValue);
        if (-1 != root->leafValue)
        {
            dbprintf(DEBUG_CYCLIC, "#Barcodes# Collides between '%s' and '%s'\n", CyclicCodes[root->leafValue].name, CyclicCodes[leafValue].name);
        }
        root->leafValue = leafValue;
        return;
    }
//#ifdef USE_SINGLE_ELEMENT_WIDTH
    int16_t c = path[0] - 1;
    if (c < 0 || c > 1) return;
//#else
//    int16_t c = path[0] + path[1] - 2;
//    if (c < 0 || c > 2) return;
//#endif
    CyclicCharacterTreeNode* child = root->children[c];
    if (!child)
    {
        child = CyclicCharacterTreeNodeCreate();
        //dbprintf(0, "#Cyclic# Add child as #%d\n", c);
        root->children[c] = child;
    }
    
    CyclicCharacterTree0Add(child, leafValue, path + 1, length - 1);
}

void CyclicCharacterTree1Add(CyclicCharacterTreeNode* root, int16_t leafValue, int16_t* path, int length) {
    if (!root) return;
    
    if (0 == length)
    {//dbprintf(0, "#Barcodes# leafValue=%d\n", leafValue);
        if (-1 != root->leafValue)
        {
            dbprintf(DEBUG_CYCLIC, "#Barcodes# Collides between '%s' and '%s'\n", CyclicCodes[root->leafValue].name, CyclicCodes[leafValue].name);
        }
        root->leafValue = leafValue;
        return;
    }
//#ifdef USE_SINGLE_ELEMENT_WIDTH
//    int16_t c = path[0] - 1;
//    if (c < 0 || c > 1) return;
//#else
    int16_t c = path[0] + path[1] - 2;
    if (c < 0 || c > 2) return;
//#endif
    CyclicCharacterTreeNode* child = root->children[c];
    if (!child)
    {
        child = CyclicCharacterTreeNodeCreate();
        //dbprintf(0, "#Cyclic# Add child as #%d\n", c);
        root->children[c] = child;
    }
    
    CyclicCharacterTree1Add(child, leafValue, path + 1, length - 1);
}

CyclicCharacterTreeNode* CyclicCharacterTreeNodeCreate() {
    CyclicCharacterTreeNode* ret = (CyclicCharacterTreeNode*) malloc(sizeof(CyclicCharacterTreeNode));
    //dbprintf(0, "#Cyclic# New node: %d\n", ++g_mallocedNodesCount);
    CyclicCharacterTreeNodeReset(ret);
    return ret;
}

CyclicCharacterTreeNode* CyclicCharacterTreeNodeNext(const CyclicCharacterTreeNode* current, int16_t c) {
    if (!current) return NULL;
    
    if (c < 0 || c > 2) return NULL;
    
    return current->children[c];
}

void CodeTrackerReset(CodeTracker* ct) {
    ct->candidate = -1;
    ct->fedElementsCount = 0;
    //    ct->startIndex = idx;
    //    for (int i=0; i<CodeElementLength; ++i) ct->window[i] = -1;
    for (int i=0; i<CyclicCodesCount; ++i) ct->probabilities[i] = 0.f;
}

CodeTracker* CodeTrackerCreate(/*int idx*/) {
    CodeTracker* ret = (CodeTracker*) malloc(sizeof(CodeTracker));
    CodeTrackerReset(ret);
    ret->next = NULL;
    return ret;
}

CodeTracker* CodeTrackerClone(const CodeTracker* ct) {
    CodeTracker* ret = (CodeTracker*) malloc(sizeof(CodeTracker));
    memcpy(ret, ct, sizeof(CodeTracker));
    return ret;
}

CyclicTrackerResult CodeTrackerFeedElement(CodeTracker* tracker, zbar_decoder_t* dcode) {
    if (tracker->candidate > -1 &&
        tracker->fedElementsCount - CodeElementLength > MaxSpaceBetweenPeriods)
    {// Failed:
        return CyclicTrackerFailed;
    }
//    int idx = tracker->fedElementsCount % CodeElementLength;
//    tracker->window[idx] = element;
    cyclic_decoder_t* decoder = &dcode->cyclic;

    int16_t c = -1;
    if (++tracker->fedElementsCount >= CodeElementLength
        && ZBAR_BAR == ((dcode->idx + CodeElementLength - 1) & 0x1))
    {
        for (int iS12 = decoder->maxS12OfChar - decoder->minS12OfChar;
             iS12 >= 0; --iS12)
        {
            int16_t s12 = decoder->minS12OfChar + iS12;
            CyclicCharacterTreeNode* nodes[2];
            nodes[0] = decoder->codeTreeRoots[0][iS12];
            nodes[1] = decoder->codeTreeRoots[1][iS12];
            int i = CodeElementLength - 1;
            for (; i > 0; --i)
            {
                int16_t singleWidth = get_width(dcode, i);
                int e0 = decode_e(singleWidth, decoder->s12, s12) + 1;
                if (e0 < 0 || e0 > 1) break;
                nodes[0] = nodes[0]->children[e0];
                if (!nodes[0]) break;
                
                int16_t pairWidth = pair_width(dcode, i - 1);
                int e1 = decode_e(pairWidth, decoder->s12, s12);
                tracker->e[CodeElementLength - 1 - i] = e1;///!!!For Debug
                if (e1 < 0 || e1 > 2) break;
                nodes[1] = nodes[1]->children[e1];
                if (!nodes[1]) break;
            }
            if (0 == i)
            {
                int16_t singleWidth = get_width(dcode, i);
                int e0 = decode_e(singleWidth, decoder->s12, s12) + 1;
                if (e0 < 0 || e0 > 1) continue;
                nodes[0] = nodes[0]->children[e0];
                if (!nodes[0]) continue;
                
                if (nodes[0]->leafValue > -1 && nodes[1]->leafValue > -1
                    && nodes[0]->leafValue == nodes[1]->leafValue &&
                    (-1 == tracker->candidate
                     || tracker->fedElementsCount - CodeElementLength <= MaxSpaceBetweenPeriods))
                {
                    c = nodes[0]->leafValue;
                    tracker->probabilities[c] += 0.7;
                }
            }
        }
    }

    if (c > -1)
    {
        if (tracker->candidate > -1)
        {
            tracker->fedElementsCount -= CodeElementLength;
        }
        else
        {
            tracker->fedElementsCount = 0;
        }
        
        float maxP = tracker->probabilities[c];
        int iMaxP = c;
        for (int i=0; i<CyclicCodesCount; ++i)
        {
            if (tracker->probabilities[i] > maxP)
            {
                maxP = tracker->probabilities[i];
                iMaxP = i;
            }
        }
        tracker->candidate = iMaxP;
        const float epsilon = 0.0001;
        int minRepeatingRequired = decoder->configs[ZBAR_CFG_MIN_REPEATING_REQUIRED - ZBAR_CFG_MIN_LEN];
        if (0 == minRepeatingRequired) minRepeatingRequired = MinRepeatingRequired;
        if (maxP >= 0.7 * minRepeatingRequired - epsilon)
        {// Confirmed:
            return CyclicTrackerConfirmed;
        }
        else
        {// Possible:
            return CyclicTrackerPossible;
        }
    }
    else if (tracker->candidate > -1 &&
             tracker->fedElementsCount - CodeElementLength > MaxSpaceBetweenPeriods)
    {// Failed:
        return CyclicTrackerFailed;
    }
    else
    {// Uncertain:
        return CyclicTrackerUncertain;
    }
}

void cyclic_destroy (cyclic_decoder_t *decoder)
{
//    dbprintf(DEBUG_CYCLIC, "#Barcodes# cyclic_destroy()\n");
    for (int j=0; j<2; ++j)
    {
        if (decoder->codeTreeRoots[j])
        {
            for (int i = decoder->maxS12OfChar - decoder->minS12OfChar; i >= 0; --i)
            {
                CyclicCharacterTreeNode* head = CyclicCharacterTreeNodeCreate();
                CyclicCharacterTreeNode* tail = head;
                head->children[0] = decoder->codeTreeRoots[j][i]; // children[0] as value, children[1] as next
                while (head)
                {
                    for (int k = 1 + j; k >= 0; --k)
                    {
                        CyclicCharacterTreeNode* parent = head->children[0];
                        if (!parent->children[k]) continue;
                        
                        tail->children[1] = CyclicCharacterTreeNodeCreate();
                        tail->children[1]->children[0] = parent->children[k];
                        tail = tail->children[1];
                    }
                    
                    free(head->children[0]);
                    //dbprintf(0, "#Cyclic# Delete node: %d\n", --g_mallocedNodesCount);
                    CyclicCharacterTreeNode* next = head->children[1];
                    free(head);
                    //dbprintf(0, "#Cyclic# Delete node: %d\n", --g_mallocedNodesCount);
                    head = next;
                }
            }
            free(decoder->codeTreeRoots[j]);
            decoder->codeTreeRoots[j] = NULL;
        }
    }

    if (decoder->s12OfChars)
    {
        free(decoder->s12OfChars);
        decoder->s12OfChars = NULL;
    }
    
    CodeTracker* ct = decoder->codeTracker.next;
    decoder->codeTracker.next = NULL;
    while (ct)
    {
        CodeTracker* tmp = ct;
        ct = ct->next;
        tmp->next = NULL;
        free(tmp);
    }
}

void cyclic_reset (cyclic_decoder_t *decoder)
{
//    dbprintf(DEBUG_CYCLIC, "#Barcodes# cyclic_reset()\n");
    cyclic_destroy(decoder);
    decoder->s12 = 0;
//    CyclicCharacterTreeNode*** charSeekers;//One group for each elements-of-character number
//    int16_t maxCharacterLength;
//    int16_t characterPhase;// This means sum of 2 elements - 2
//    int16_t* s12OfChar;
//    int16_t minS12OfChar;
//    int16_t maxS12OfChar;
    decoder->s12OfChars = (int16_t*) malloc(sizeof(int16_t) * CyclicCodesCount);
    decoder->minS12OfChar = 127;
    decoder->maxS12OfChar = 0;
    decoder->maxCodeLength = 0;
    for (int i = CyclicCodesCount - 1; i >= 0; --i)
    {
        int16_t* seq = CyclicCodes[i].elementSequence;
        int length = sizeof(CyclicCodes[i].elementSequence) / sizeof(CyclicCodes[i].elementSequence[0]);
        if (length > decoder->maxCodeLength)
        {
            decoder->maxCodeLength = length;
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
    const int uniqueS12Count = decoder->maxS12OfChar - decoder->minS12OfChar + 1;
//    dbprintf(0, "#Barcodes# minS12OfChar:%d, maxS12OfChar:%d, charSeekersCount:%d\n", decoder->minS12OfChar, decoder->maxS12OfChar, decoder->charSeekersCount);
//    decoder->charSeekers = (CyclicCharacterTreeNode***) malloc(sizeof(CyclicCharacterTreeNode**) * decoder->maxCodeLength);
//    decoder->candidates = (int16_t**) malloc(sizeof(int16_t*) * decoder->maxCodeLength);
//    decoder->repeatingCounts = (int16_t**) malloc(sizeof(int16_t*) * decoder->maxCodeLength);
//    for (int i = decoder->maxCodeLength - 1; i >= 0; --i)
//    {
//        decoder->charSeekers[i] = (CyclicCharacterTreeNode**) malloc(sizeof(CyclicCharacterTreeNode*) * uniqueS12Count);
//        memset(decoder->charSeekers[i], 0, sizeof(CyclicCharacterTreeNode*) * uniqueS12Count);
//        decoder->candidates[i] = (int16_t*) malloc(sizeof(int16_t) * uniqueS12Count);
//        decoder->repeatingCounts[i] = (int16_t*) malloc(sizeof(int16_t) * uniqueS12Count);
//        for (int j = uniqueS12Count - 1; j >= 0; --j)
//        {
//            decoder->candidates[i][j] = -1;
//            decoder->repeatingCounts[i][j] = 0;
//        }
//    }
//    decoder->characterPhase = 0;

    for (int j=0; j<2; ++j)
    {
        decoder->codeTreeRoots[j] = (CyclicCharacterTreeNode**) malloc(sizeof(CyclicCharacterTreeNode*) * uniqueS12Count);
        for (int i = uniqueS12Count - 1; i >= 0; --i)
        {
            decoder->codeTreeRoots[j][i] = CyclicCharacterTreeNodeCreate();
        }
        for (int i = CyclicCodesCount - 1; i >= 0; --i)
        {
            int16_t* seq = CyclicCodes[i].elementSequence;
            int length = sizeof(CyclicCodes[i].elementSequence) / sizeof(CyclicCodes[i].elementSequence[0]);
            if (0 == j)
            {
#ifdef USE_SINGLE_TREE
                CyclicCharacterTree0Add(decoder->codeTreeRoots[0][0], i, seq, length);
#else //#ifdef USE_SINGLE_TREE
                CyclicCharacterTree0Add(decoder->codeTreeRoots[0][decoder->s12OfChars[i] - decoder->minS12OfChar], i, seq, length);
#endif //#ifdef USE_SINGLE_TREE
            }
            else
            {
#ifdef USE_SINGLE_TREE
                CyclicCharacterTree1Add(decoder->codeTreeRoots[1][0], i, seq, length - 1);
#else //#ifdef USE_SINGLE_TREE
                CyclicCharacterTree1Add(decoder->codeTreeRoots[1][decoder->s12OfChars[i] - decoder->minS12OfChar], i, seq, length - 1);
#endif //#ifdef USE_SINGLE_TREE
            }

        }
    }
    
    CodeTrackerReset(&decoder->codeTracker);
    CodeTracker* ct = decoder->codeTracker.next;
    decoder->codeTracker.next = NULL;
    while (ct)
    {
        CodeTracker* tmp = ct;
        ct = ct->next;
        tmp->next = NULL;
        free(tmp);
    }
}

zbar_symbol_type_t _zbar_decode_cyclic (zbar_decoder_t *dcode)
{
//    if (!dcode) return(ZBAR_NONE);
//    if (dcode->scanDX != 1 || dcode->scanDY != 0) return(ZBAR_NONE);
    switch (dcode->rotationZ) {
        case 0:
            if (dcode->scanDX != 1 || dcode->scanDY != 0)
                return(ZBAR_NONE);
            break;
        case 1:
            if (dcode->scanDY != -1 || dcode->scanDX != 0)
                return(ZBAR_NONE);
            break;
        case 2:
            if (dcode->scanDX != -1 || dcode->scanDY != 0)
                return(ZBAR_NONE);
            break;
        case 3:
            if (dcode->scanDY != 1 || dcode->scanDX != 0)
                return(ZBAR_NONE);
            break;
        default:
            break;
    }

    zbar_symbol_type_t ret = ZBAR_NONE;
    cyclic_decoder_t* decoder = &dcode->cyclic;
    decoder->s12 -= get_width(dcode, CodeElementLength);
    decoder->s12 += get_width(dcode, 0);
    
    CodeTracker* tracker = &decoder->codeTracker;
    while (tracker)
    {
        CyclicTrackerResult result = CodeTrackerFeedElement(tracker, dcode);
        if (CyclicTrackerConfirmed == result)
        {
            release_lock(dcode, ZBAR_CYCLIC);
            int length = (int)(strlen(CyclicCodes[tracker->candidate].name) + 1);
            size_buf(dcode, length);
            memcpy(dcode->buf, CyclicCodes[tracker->candidate].name, length);
            dcode->buflen = length;
            dbprintf(DEBUG_CYCLIC, "#Barcodes# Confirm '%s', dx=%d, dy=%d\n", CyclicCodes[tracker->candidate].name, dcode->scanDX, dcode->scanDY);
            ret = ZBAR_CYCLIC;
            return ret;
        }
        else if (&decoder->codeTracker == tracker)
        {
            if (CyclicTrackerPossible == result)
            {
                acquire_lock(dcode, ZBAR_CYCLIC);
                ret = ZBAR_PARTIAL;
                dbprintf(DEBUG_CYCLIC, "#Barcodes# Possible '%s', dx=%d, dy=%d\n", CyclicCodes[tracker->candidate].name, dcode->scanDX, dcode->scanDY);
                bool isNewFound = true;
                for (CodeTracker* node = tracker->next;
                     node != NULL;
                     node = node->next)
                {
                    if (node->candidate == tracker->candidate)
                    {
                        isNewFound = false;//TODO: Merge the results
                        break;
                    }
                }
                CodeTracker* newNode;
                if (isNewFound)
                {
                    newNode = CodeTrackerClone(tracker);
                    newNode->next = tracker->next;
                    tracker->next = newNode;
                }
                else
                {
                    newNode = tracker;
                }
                
                CodeTrackerReset(tracker);
                tracker->fedElementsCount = CodeElementLength - 1;
                
                tracker = newNode->next;
                continue;
            }
        }
        else if (CyclicTrackerFailed == result)
        {
            //TODO:
        }
        
        tracker = tracker->next;
    }

    return(ret);
}
