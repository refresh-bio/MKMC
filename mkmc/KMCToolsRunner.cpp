#include <cstring>
#include <iostream>
#include <filesystem>
#include <cstdint>
#include <algorithm>
#if defined(WIN32) || defined(_WIN32)
#define NOMINMAX
#include <direct.h>
#include <shlwapi.h>
#else
#include <unistd.h>
#include <libgen.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <iterator>
#include <algorithm>
#include <sstream>
#endif
#include "KMCToolsRunner.h"



void KMCToolsRunner::operator()(TasksPool& tasksPool)
{
    std::string inputFile, outputFile;
    while (tasksPool.getTask(inputFile, outputFile))
    {
        KMC::Runner runner;

        std::ostringstream sstream;

        sstream << "-t" << params.mkmcParams.nThreads;
        sstream << " -hp";

        sstream << " transform " << inputFile;
        sstream << " sort " << outputFile;

        unsigned long result;
        if (!runCommand(KMC_TOOLS_EXECUTABLE_NAME, sstream.str(), result) || result != 0) {
            std::cerr << "ERROR: cannot run kmc_tools." << std::endl;
        }
    }
}

bool KMCToolsRunner::checkToolsRequired(const std::string& kmcOutputFile)
{
    CKMCFile file;
    if (file.OpenForListing(kmcOutputFile))
    {
        return file.IsKMC2();
    }
    else
    {
        std::cerr << "ERROR: cannot open temporary file " << kmcOutputFile << std::endl;
        std::exit(1);
        return false;
    }
}

void KMCToolsRunner::runKMCToolsParallel()
{
    std::vector<std::string> inputFiles, outputFiles;
    for (size_t i = 0; i < params.mkmcParams.kmcOutputFiles.size(); ++i)
    {
        const std::string& kmcOutputFile = params.mkmcParams.kmcOutputFiles[i];
        const std::string& toolsOutputFile = params.mkmcParams.toolsOutputFiles[i];
        if (checkToolsRequired(kmcOutputFile))
        {
            inputFiles.push_back(kmcOutputFile);
            outputFiles.push_back(toolsOutputFile);
        }
        else
        {
            std::filesystem::rename(kmcOutputFile + ".kmc_pre", toolsOutputFile + ".kmc_pre");
            std::filesystem::rename(kmcOutputFile + ".kmc_suf", toolsOutputFile + ".kmc_suf");
        }
    }
    TasksPool tasksPool(inputFiles, outputFiles);

    (*this)(tasksPool);

    //std::vector<std::thread> threads(std::min(static_cast<size_t>(params.mkmcParams.nKMCWorkers), inputFiles.size()));
    //for (uint32_t i_thred = 0; i_thred < std::min(static_cast<size_t>(params.mkmcParams.nKMCWorkers), inputFiles.size()); ++i_thred)
    //{
    //    threads[i_thred] = std::thread([this, &tasksPool] { (*this)(tasksPool); });
    //}

    //for (std::thread& thread : threads)
    //{
    //    thread.join();
    //}
}

#if defined(WIN32) || defined(_WIN32) //Windows

#pragma warning(push)
#pragma warning(disable: 4996)
bool KMCToolsRunner::runCommand(std::string command, const std::string& args, unsigned long& processResult, char* buffer /*= NULL*/, const int bufferSize /*= 0*/) {
    // Get directory path
    const DWORD fileNameBufferSize = 261;
    char fileNameBuffer[fileNameBufferSize];
    if (GetModuleFileName(NULL, fileNameBuffer, fileNameBufferSize) > 0) {
        PathRemoveFileSpec(fileNameBuffer);
        command = std::string(fileNameBuffer) + "\\" + command;
    }

    PROCESS_INFORMATION processInfo;
    STARTUPINFO processStartupInfo;

    ZeroMemory(&processInfo, sizeof(PROCESS_INFORMATION));

    ZeroMemory(&processStartupInfo, sizeof(STARTUPINFO));
    processStartupInfo.cb = sizeof(STARTUPINFO);
    processStartupInfo.hStdError = GetStdHandle(STD_OUTPUT_HANDLE);

    HANDLE pipeReadHandle = NULL;
    HANDLE pipeWriteHandle = NULL;

    if (buffer != NULL) {
        if (!createPipe(pipeReadHandle, pipeWriteHandle)) {
            return false;
        }
        processStartupInfo.hStdOutput = pipeWriteHandle;
    }
    else {
        processStartupInfo.hStdOutput = GetStdHandle(STD_OUTPUT_HANDLE);
    }

    processStartupInfo.hStdInput = NULL;
    processStartupInfo.dwFlags |= STARTF_USESTDHANDLES;

    // Fullfill the ANSI specifications by passing command name as argv[0].
    // Quote command to prevent from path spaces problems.
    std::string argsStr = "\"";
    argsStr += command;
    argsStr += "\"";

    char* cArgs = new char[argsStr.length() + 1 + args.length() + 1];
    std::size_t copied = argsStr.copy(cArgs, argsStr.length());
    cArgs[copied] = ' ';
    ++copied;
    copied += args.copy(cArgs + copied, args.length());
    cArgs[copied] = '\0';

    BOOL createProcessSuccess = CreateProcess(command.c_str(), cArgs, NULL, NULL, TRUE, 0, NULL, NULL, &processStartupInfo, &processInfo);

    delete[] cArgs;

    if (!createProcessSuccess) {
        return false;
    }
    else {
        WaitForSingleObject(processInfo.hProcess, INFINITE);

        if (buffer != NULL) {
            DWORD bytesRead;

            BOOL pipeReadSuccess = PeekNamedPipe(pipeReadHandle, buffer, bufferSize - 1, &bytesRead, NULL, NULL);

            if (!pipeReadSuccess || bytesRead == 0) {
                return false;
            }
            buffer[bytesRead] = '\0';

            destroyPipe(pipeReadHandle, pipeWriteHandle);
        }

        if (!GetExitCodeProcess(processInfo.hProcess, &processResult)) {
            return false;
        }

        CloseHandle(processInfo.hProcess);
        CloseHandle(processInfo.hThread);

        return true;
    }
}
#pragma warning(pop)

bool KMCToolsRunner::createPipe(HANDLE& pipeReadHandle, HANDLE& pipeWriteHandle) {
    SECURITY_ATTRIBUTES securityAttributes;
    securityAttributes.nLength = sizeof(SECURITY_ATTRIBUTES);
    securityAttributes.bInheritHandle = TRUE;
    securityAttributes.lpSecurityDescriptor = NULL;

    if (!CreatePipe(&pipeReadHandle, &pipeWriteHandle, &securityAttributes, 0)) {
        return false;
    }

    if (!SetHandleInformation(pipeReadHandle, HANDLE_FLAG_INHERIT, 0)) {
        return false;
    }

    return true;
}

void KMCToolsRunner::destroyPipe(HANDLE& pipeReadHandle, HANDLE& pipeWriteHandle) {
    CloseHandle(pipeReadHandle);
    CloseHandle(pipeWriteHandle);

    pipeReadHandle = NULL;
    pipeWriteHandle = NULL;
}

#else // Linux

bool KMCToolsRunner::runCommand(std::string command, const std::string& args, unsigned long& processResult, char* buffer /*= NULL*/, const int bufferSize /*= 0*/) {
    // Get directory path
    const int linkNameBufferSize = 32;
    const int fileNameBufferSize = 261;
    char linkName[linkNameBufferSize];
    char fileNameBuffer[fileNameBufferSize];
    sprintf(linkName, "/proc/%d/exe", getpid());
    int bytesInserted = readlink(linkName, fileNameBuffer, fileNameBufferSize - 1);
    if (bytesInserted >= 0) {
        fileNameBuffer[bytesInserted] = '\0';
        dirname(fileNameBuffer);
        command = std::string(fileNameBuffer) + std::filesystem::path::preferred_separator + command;
    }

    if (buffer == NULL) {
        int pid = fork();
        if (pid == -1) {
            return false;
        }
        else if (pid == 0) {
            // Fullfill the ANSI specifications by passing command name as argv[0].
            std::istringstream sstream(args);
            std::vector<std::string> splittedArgs;
            splittedArgs.push_back(command);
            copy(std::istream_iterator<std::string>(sstream), std::istream_iterator<std::string>(), back_inserter(splittedArgs));
            char** cSplittedArgs = new char* [splittedArgs.size() + 1];
            for (unsigned i = 0; i < splittedArgs.size(); ++i) {
                cSplittedArgs[i] = new char[splittedArgs[i].length() + 1];
                std::size_t copied = splittedArgs[i].copy(cSplittedArgs[i], splittedArgs[i].length());
                cSplittedArgs[i][copied] = '\0';
            }
            cSplittedArgs[splittedArgs.size()] = NULL;

            execv(cSplittedArgs[0], cSplittedArgs);
            for (unsigned i = 0; i < splittedArgs.size(); ++i) {
                delete[] cSplittedArgs[i];
            }
            delete[] cSplittedArgs;
            exit(EXIT_FAILURE);
        }
        else {
            int status;
            if (waitpid(pid, &status, 0) < 0) {
                return false;
            }

            processResult = WEXITSTATUS(status);
            return true;
        }
    }
    else {
        // Replace space to \space
        std::size_t commandPos = 0;
        while ((commandPos = command.find(" ", commandPos)) != std::string::npos) {
            command.replace(commandPos, 1, "\\ ");
            commandPos += 2;
        }

        std::string childCommand = command + " " + args;

        FILE* pipe = popen(childCommand.c_str(), "r");
        if (pipe == NULL) {
            return false;
        }

        std::size_t bytesRead = fread(buffer, sizeof(char), bufferSize - 1, pipe);
        if (bytesRead == 0) {
            return false;
        }
        buffer[bytesRead] = '\0';

        int result = pclose(pipe);
        processResult = WEXITSTATUS(result);

        return true;
    }
}

#endif