#include "FileGenerators.h"
#include <memory>



kmcdb::WriterSortedPlain<uint64_t>* BinFileGenerator::kmcDBWriter = nullptr;
bool BinFileGenerator::writerWasOpened = false;
DumpWriter* MatrixFileGenerator::dumpWriter = nullptr;
bool MatrixFileGenerator::writerWasOpened = false;
DumpWriter* FASTAFileGenerator::dumpWriter = nullptr;
bool FASTAFileGenerator::writerWasOpened = false;
