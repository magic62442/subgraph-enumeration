//
// Created by anonymous authors on 2024/3/4.
//

#include "command.h"
#include <iostream>

Command::Command(int argc, char **argv) : CommandParser(argc, argv){
    // 检查是否请求帮助
    if (commandOptionExists("--help") || commandOptionExists("-h")) {
        printHelp();
        exit(0);
    }
    
    optionsKey[OptionKeyword::QueryGraphPath] = "-q";
    optionsKey[OptionKeyword::DataGraphPath] = "-d";
    optionsKey[OptionKeyword::ResultPath] = "-r";
    optionsKey[OptionKeyword::IntersectType] = "-i";
    optionsKey[OptionKeyword::UseIEP] = "-u";
    optionsKey[OptionKeyword::ThreadNumber] = "-t";
    optionsKey[OptionKeyword::Help] = "--help";
    intOptionValue[OptionKeyword::ThreadNumber] = 1;
    booleanOptionValue[OptionKeyword::UseIEP] = false;
    booleanOptionValue[OptionKeyword::Help] = false;
    booleanOptionValue[OptionKeyword::IntersectType] = false;

    processOptions();
}

void Command::processOptions() {
    optionsValue[OptionKeyword::QueryGraphPath] = getCommandOption(optionsKey[OptionKeyword::QueryGraphPath]);
    optionsValue[OptionKeyword::DataGraphPath] = getCommandOption(optionsKey[OptionKeyword::DataGraphPath]);
    optionsValue[OptionKeyword::ResultPath] = getCommandOption(optionsKey[OptionKeyword::ResultPath]);
    intOptionValue[OptionKeyword::ThreadNumber] = commandOptionType(optionsKey[OptionKeyword::ThreadNumber]);
    booleanOptionValue[OptionKeyword::UseIEP] = commandOptionExists(optionsKey[OptionKeyword::UseIEP]);
    booleanOptionValue[OptionKeyword::Help] = commandOptionExists(optionsKey[OptionKeyword::Help]);
    booleanOptionValue[OptionKeyword::IntersectType] = commandOptionExists(optionsKey[OptionKeyword::IntersectType]);
}

int Command::getThreadNumber() {
    return intOptionValue[OptionKeyword::ThreadNumber];
}


void Command::printHelp() const {
    std::cout << "Usage: graph_mining_system [options]\n"
              << "Options:\n"
              << "  -q <path>        Query graph file path\n"
              << "  -d <path>        Data graph file path\n"
              << "  -r <path>        Result file path (optional, print to stdout if not specified)\n"
              << "  -i               Input type (default: false)\n"
              << "  -u               Use IEP (default: false)\n"
              << "  -t <number>      Number of threads for parallel execution (default: 1)\n"
              << "  --help, -h       Show this help message\n"
              << std::endl;
}