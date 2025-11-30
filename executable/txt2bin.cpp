#include <iostream>
#include <string>
#include "../relation/graph.h"

int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cout << "Usage: " << argv[0] << " <input_text_file> <output_binary_file>" << std::endl;
        std::cout << "Convert a text graph file to binary format." << std::endl;
        return 1;
    }

    std::string inputFile = argv[1];
    std::string outputFile = argv[2];

    try {
        // Load graph from text file
        DataGraph graph;
        graph.loadDataGraph(inputFile);

        std::cout << "Loaded graph with " << graph.getNumVertices()
                  << " vertices and " << graph.getNumEdges() / 2 << " edges." << std::endl;

        // Save graph to binary file
        graph.saveBinaryGraph(outputFile);

        std::cout << "Successfully converted " << inputFile << " to " << outputFile << std::endl;

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    } catch (...) {
        std::cerr << "Unknown error occurred during conversion." << std::endl;
        return 1;
    }

    return 0;
}