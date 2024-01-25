#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_LINE_LENGTH 2048

typedef struct {
    char* id;
    char* seq;
    char* qual;
} FastqRecord;

typedef struct {
    FastqRecord* records;
    size_t size;
    size_t capacity;
} FastqRecordArray;

void initFastqRecordArray(FastqRecordArray* array, size_t initialCapacity) {
    array->records = (FastqRecord*)malloc(initialCapacity * sizeof(FastqRecord));
    array->size = 0;
    array->capacity = initialCapacity;
}

void pushFastqRecord(FastqRecordArray* array, const char* id, const char* seq, const char* qual) {
    if (array->size == array->capacity) {
        array->capacity *= 2;
        array->records = (FastqRecord*)realloc(array->records, array->capacity * sizeof(FastqRecord));
    }

    array->records[array->size].id = strdup(id);
    array->records[array->size].seq = strdup(seq);
    array->records[array->size].qual = strdup(qual);

    array->size++;
}

void freeFastqRecordArray(FastqRecordArray* array) {
    for (size_t i = 0; i < array->size; ++i) {
        free(array->records[i].id);
        free(array->records[i].seq);
        free(array->records[i].qual);
    }
    free(array->records);
}

void shuffleIndices(size_t* indices, size_t size) {
    for (size_t i = size - 1; i > 0; --i) {
        size_t j = rand() % (i + 1);
        size_t temp = indices[i];
        indices[i] = indices[j];
        indices[j] = temp;
    }
}

void splitFastq(const char* inputFile1, const char* inputFile2, const char* outputPrefix, size_t numParts) {
    FILE* fp1 = fopen(inputFile1, "r");
    FILE* fp2 = fopen(inputFile2, "r");

    if (!fp1 || !fp2) {
        perror("Error opening input files");
        return;
    }

    FastqRecordArray records1, records2;
    initFastqRecordArray(&records1, 1000);
    initFastqRecordArray(&records2, 1000);

    char line[MAX_LINE_LENGTH];

    while (fgets(line, MAX_LINE_LENGTH, fp1) != NULL) {
        char id[MAX_LINE_LENGTH], seq[MAX_LINE_LENGTH], qual[MAX_LINE_LENGTH];
        sscanf(line, "%s", id);
        fgets(seq, MAX_LINE_LENGTH, fp1);  // Read the sequence line
        fgets(line, MAX_LINE_LENGTH, fp1);  // Skip '+'
        fgets(qual, MAX_LINE_LENGTH, fp1);  // Read the quality line

        pushFastqRecord(&records1, id, seq, qual);
    }

    while (fgets(line, MAX_LINE_LENGTH, fp2) != NULL) {
        char id[MAX_LINE_LENGTH], seq[MAX_LINE_LENGTH], qual[MAX_LINE_LENGTH];
        sscanf(line, "%s", id);
        fgets(seq, MAX_LINE_LENGTH, fp2);  // Read the sequence line
        fgets(line, MAX_LINE_LENGTH, fp2);  // Skip '+'
        fgets(qual, MAX_LINE_LENGTH, fp2);  // Read the quality line

        pushFastqRecord(&records2, id, seq, qual);
    }

    fclose(fp1);
    fclose(fp2);

    size_t numRecords = records1.size;
    if (numParts <= 0 || numParts > numRecords) {
        fprintf(stderr, "Invalid number of parts\n");
        freeFastqRecordArray(&records1);
        freeFastqRecordArray(&records2);
        return;
    }

    size_t* indices = (size_t*)malloc(numRecords * sizeof(size_t));
    for (size_t i = 0; i < numRecords; ++i) {
        indices[i] = i;
    }

    shuffleIndices(indices, numRecords);

    size_t recordsPerPart = numRecords / numParts;

    for (size_t partNumber = 0; partNumber < numParts; ++partNumber) {
        char partOutputFile1[MAX_LINE_LENGTH];
        char partOutputFile2[MAX_LINE_LENGTH];
        sprintf(partOutputFile1, "%s_part%zu_1.fastq", outputPrefix, partNumber + 1);
        sprintf(partOutputFile2, "%s_part%zu_2.fastq", outputPrefix, partNumber + 1);

        FILE* outFp1 = fopen(partOutputFile1, "w");
        FILE* outFp2 = fopen(partOutputFile2, "w");

        size_t startIdx = partNumber * recordsPerPart;
        size_t endIdx = (partNumber == numParts - 1) ? numRecords : (partNumber + 1) * recordsPerPart;

        for (size_t i = startIdx; i < endIdx; ++i) {
            fprintf(outFp1, "%s\n%s+\n%s", records1.records[indices[i]].id,
                    records1.records[indices[i]].seq, records1.records[indices[i]].qual);

            fprintf(outFp2, "%s\n%s+\n%s", records2.records[indices[i]].id,
                    records2.records[indices[i]].seq, records2.records[indices[i]].qual);
        }

        fclose(outFp1);
        fclose(outFp2);
    }

    free(indices);
    freeFastqRecordArray(&records1);
    freeFastqRecordArray(&records2);
}

int main(int argc, char* argv[]) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s <inputFile1> <inputFile2> <outputPrefix> <numParts>\n", argv[0]);
        return 1;
    }

    const char* inputFile1 = argv[1];
    const char* inputFile2 = argv[2];
    const char* outputPrefix = argv[3];
    size_t numParts = atoi(argv[4]);

    if (numParts <= 0) {
        fprintf(stderr, "Invalid number of parts\n");
        return 1;
    }

    splitFastq(inputFile1, inputFile2, outputPrefix, numParts);

    return 0;
}

