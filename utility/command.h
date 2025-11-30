//
// Created by anonymous authors on 2024/3/4.
//

#ifndef ASDMATCH_COMMAND_H
#define ASDMATCH_COMMAND_H

#include "command_parser.h"
#include <iostream>
#include <map>

enum OptionKeyword {
    QueryGraphPath = 1,     // -q, the query graph file path
    DataGraphPath = 2,      // -d, the data graph file path
    ResultPath = 3,         // -r, the result file path
    IntersectType = 4,          // -i, intersect type
    UseIEP = 5,             // -u, flag of using iep
    ThreadNumber = 6,       // -t, the number of threads
    Help = 7               // --help, show help information
};

class Command : public CommandParser {
private:
    std::map<OptionKeyword, std::string> optionsKey;
    std::map<OptionKeyword, std::string> optionsValue;
    std::map<OptionKeyword, bool> booleanOptionValue;
    std::map<OptionKeyword, int> intOptionValue;

private:
    void processOptions();

public:
    Command(int argc, char **argv);

    std::string getQueryGraphPath() {
        return optionsValue[OptionKeyword::QueryGraphPath];
    }

    std::string getDataGraphPath() {
        return optionsValue[OptionKeyword::DataGraphPath];
    }

    std::string getResultPath() {
        return optionsValue[OptionKeyword::ResultPath];
    }

    bool getIntersectType() {
        return booleanOptionValue[OptionKeyword::IntersectType];
    }

    bool getUseIEP() {
        return booleanOptionValue[OptionKeyword::UseIEP];
    }

    int getThreadNumber();

    void printHelp() const;
};


#endif //ASDMATCH_COMMAND_H
