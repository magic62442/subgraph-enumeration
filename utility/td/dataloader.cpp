#ifndef DATALOADER_HZY_CPP
#define DATALOADER_HZY_CPP

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <set>
#include "hypergraph.h"


HyperG BuildHyperGraph(const char * file_name) {
    /*
        the data should be organized as edge_name(element 1, element 2, ...), line by line, the detail can refer Hyperbench
    */
    std::ifstream fin;
    // ofstream fout;

    fin.open(file_name, std::ios::in);

    std::string str;
    std::map<std::string, size_t> f;
    size_t cnt = 0;
    std::vector <std::set<std::string> > v;
    while(getline(fin, str)) {
        size_t p_L = str.find_first_of("(");
        size_t p_R = str.find_last_of(")");
        if(p_L == str.npos || p_R == str.npos || p_R  <= p_L + 1)
            continue; // skip this line
            
        str = str.substr(p_L + 1, p_R - p_L - 1);
        std::string curr = "";
        std::set <std::string> temp;
        for(size_t i = 0; i < str.size(); ++i) {
            if(str[i] == ' ' && curr == "")
                continue;
            if(str[i] == ',') {
                if(curr != "") {
                    if(!f.count(curr))
                        f[curr] = cnt++;
                    temp.insert(curr);
                    curr = "";
                }
            }
            else
                curr += str[i];
        }
        if(curr != "") {
            temp.insert(curr);
            if(!f.count(curr))
                f[curr] = cnt++;
        }
        v.push_back(temp);
    }
    
    return HyperG(cnt, v, f);
}

#endif