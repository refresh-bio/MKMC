#include "FileGenerators.h"
#include <memory>



kmcdb::WriterSortedPlain<uint64_t>* BinFileGenerator::kmcDBWriter = nullptr;
bool BinFileGenerator::writerWasOpened = false;
TextFileWriter* MatrixFileGenerator::dumpWriter = nullptr;
bool MatrixFileGenerator::writerWasOpened = false;
TextFileWriter* FASTAFileGenerator::dumpWriter = nullptr;
bool FASTAFileGenerator::writerWasOpened = false;
