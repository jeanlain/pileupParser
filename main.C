//
//  main.c
//  piluepParser
//
//  Created by Jean Peccoud on 26/04/2014.
//  Copyright (c) 2014 Jean Peccoud. All rights reserved.
//

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

int main(int argc, const char * argv[]) {
    if (argc < 0) {
        printf("pileupParser by Jean Peccoud\nUsage: pileupParser min_read_count(integer) min_number_of_alleles(integer>=1, default = 2) min_phred_value(integer, default 0) max_individuals_to_consider(integer, default all, max 100)\n");
        exit(0);
    }
    short MAC = 1;
    if (argc > 1) MAC = atoi(argv[1]);
    short minAllele = 2;
    if (argc > 2) minAllele = atoi(argv[2]);
    char minPhred = 0+33;
    if (argc > 3) minPhred = atoi(argv[3])+33;
    short nInds = -1;
    if (argc > 4) nInds = atoi(argv[4]);

    size_t const max = 2000000;     //max size of a line to read
    char line[max];
    char string[200000];           //will hold the base string (CIGAR) for a given individual (tolerates coverage up to 10000), bases will be encoded as numbers corresponding to the indices of counts arrays below


    int i, j;    //incrementers
    int p = 0; short field = 1; //position in the pileup line, indice of next field (string between tabs) to scan

    
    short indelSize = 0;        //size of a new described indel (used to know which character to skip in the cigar)
    char bases[] = "xACGT";     // correspondance between bases numbers and letters, base 0 is the base that is in the reference fasta, whatever its nature
    int counts[100][6] = {0};    //allele counts for all bases in individuals at a position, tolerates up to 100 individuals
    int allCounts[6] = {0};      // total allele counts for each base
    int readCount, ind;         //total number of read in of individual (its coverage) at a position, individual number
    int nIndels;                //number of reads having indels at that positions (all inds)
    enum {period = 46, coma = 44, a = 97, A = 65, c = 99, C = 67, g = 103, G = 71, t = 116, T = 84, n = 110, N = 78, asterisk = 42, plus = 43, minus = 45, caret = 94, tab = 9, newLine = 10};  //ascii codes for different caracters that we need. I think it is faster than using a built-in function to convert ascii to their character representations, while still making it easy to understand the code.
        // this doesn't consider ">" or "<" in the CIGAR. I don't know what they're supposed to mean as the doc for mplileup is unclear. It just says "reference skip" (???)
        
    while (fgets(line, max, stdin) != NULL) {
        if (nInds == -1) {   //counts the number of individuals the first time
            while (line[p] != newLine) {
                if (line[p++] == tab) field++;  //increments the number of field after each tab
            }
            nInds = (field-3)/3;// IMPORTANT: assumes that pileup as a coverage field for each individual
            if (nInds < 1) {
                fprintf(stderr, "%d fields, insufficient number of fields in input file. Exiting...", field);
                exit(1);
            }
            if (nInds > 100) nInds = 100;  //limits the number of individuals to 100
        }
        p = 0;
        nIndels = 0,
        readCount = 0;
        ind = 0;
        field = 1;
        short d = 0;   //line position just after the reference base string
        do {
            while (field < 4) {     //we skip the first 3 fields for now, we're only interested in them if there is a snp
                while(line[p++] != tab)
                    ;  // reads the field (until tab character)
                field++;
                if (field == 3) {   // the reference base field
                    d = p-1;
                    bases[0] = line[p];
                }
            }
            if ((field-4) % 3 == 0) {     //we aslo skip the coverage field as we compute our own count
                while(line[p++] != tab)
                    ;
                field++;
            }
            if ((field-5) % 3 == 0) {   //the base string (cigar)
                i = 0;   //number of reads with actual bases
                ind++;   // we increment the number of individuals scanned here
                while(line[p++] != tab)  {
                    char x = line[p-1];
                    switch (x) {
                        case period:
                        case coma:          //both meaning the reference base
                            string[i++]=0;  //it will alway be counted in the first position of the count arrays
                            break;
                        case a:
                        case A:
                            string[i++]=1;
                            break;
                        case c:
                        case C:
                            string[i++]=2;
                            break;
                        case g:
                        case G:
                            string[i++]=3;
                            break;
                        case t:
                        case T:
                            string[i++]=4;
                            break;
                        case n:
                        case N:
                            string[i++]=5;  //this is lumped into base number 5, which we won't record in the final counts.
                            break;
                        case asterisk:    //"*" indicates a deletion
                            string[i++]=5;  //we record it since the read still has a phred value at that position (not sure what it corresponds to)
                            nIndels++;
                            break;
                        case plus:       // "+" when a new insertion is described
                        case minus:        // "-" for a deletion
                            j = 0; indelSize = 0;
                            nIndels++;
                            while (line[p+j] >= 48 && line[p+j] <= 57) {  //the following char(s) encode the size of indel in decimal, so are comprised between "0" and "9"
                                indelSize = indelSize*10 + line[p+j]-48;
                                j++;
                            }
                            p += j+indelSize;   //we skip the indel description as it doesn't correspond to reads
                            break;
                        case caret:       //"^", begining of a read
                            p++;        // we skip the next char, which is not a base (the mapping quality)
                            break;
                        default:
                            break;
                    }
                }
                field++;
                readCount = i;  //the readcount is used to read the quality string faster, without controling for the tab caracter as we know the field is at least this long
            }
            if ((field-6) % 3 == 0) {   //the base quality field. We use it to counts the alleles for each ind, since we have a quality filter
                for (i = 0; i < readCount; i++) {
                    if (line[p++]> minPhred) {
                        counts[ind-1][string[i]]++;
                        allCounts[string[i]]++;
                    }
                }
                while(line[p++] != tab);
                field++;
            }
        } while (ind < nInds);
        
        short basesToCount[5] = {0};  //will indicate which bases are present at that position (referred to by their numbers)
        short nAlleles = 0;     //numbers of alleles (different bases found)
        line[d++] = tab;
        j = 0;
        for (i = 0; i < 5; i++) {  //for each alternative base (remember that base 0 is the reference base)
            if (allCounts[i] >= MAC || (i == 0 & bases[0] != N)) {  //if we have enough reads for a base
                line[d++] = bases[i]; line[d++] = coma;  //we add the base to the line after the ref base field (so we reuse the first fields of the line for the output)
                nAlleles++;
                basesToCount[j++] = i;  //the base will be counted (note that basesToCount[0] == 0, for the reference base)
            }
        }
        if (nAlleles >= minAllele) {  // we only print out positions with variation
            line[d-1] = tab; //d++;   // adds a tab between fields, replaces the trailing coma after the base field
            line[d++] = nAlleles + 48;  //prints the number of alleles, +48 converts the number (from 0 to 9) to the ascii cacacter
            line[d++] = tab;  
            line[d] = 0;   //this turns the line into a string so we can concatenate new strings to it
            char countString[9];                    //temporary string that will hold count values
            for (i = 0; i < nAlleles; i++) {        //we write the total counts for each base
                d+= sprintf(countString, "%u,", allCounts[basesToCount[i]]);
                strcat(line, countString);
            }
            line[d-1] = tab;              //replaces last coma with tab
            for (j=0; j < nInds; j++) {   //we write the individual counts for each base
                for (i = 0; i < nAlleles; i++) {
                    d+= sprintf(countString, "%u,", counts[j][basesToCount[i]]);
                    strcat(line, countString);
                }
                line[d-1] = tab;
            }
            d+= sprintf(countString, "%u\n", nIndels);  //we finish by the number of indels
            strcat(line, countString);
          //  printf("%s", line);
            fputs(line,stdout);
            
        }
        for (i = 0; i < nAlleles; i++) {   //we reinitialize the counts
            allCounts[basesToCount[i]] = 0;
            for  (j=0; j < nInds; j++){
                counts[j][basesToCount[i]] = 0;
            }
        }
    }
    return 0;
}

