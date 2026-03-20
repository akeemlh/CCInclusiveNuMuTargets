#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include "TFile.h"
#include "TTree.h"

//run with root -l -b -q 'convertTxtToRoot.cpp("path/to/your_data.txt", "path/to/output.root")'

void convertTxtToRoot(const char* inputPath, const char* outputPath) {
    // 1. Open the input file
    std::ifstream infile(inputPath);
    if (!infile.is_open()) {
        std::cerr << "Error: Could not open input file: " << inputPath << std::endl;
        return;
    }

    // 2. Create output ROOT file and Tree
    TFile *outFile = new TFile(outputPath, "RECREATE");
    TTree *tree = new TTree("FriendTree", "Tree converted from CSV");

    // 3. Define variables for branches
    Int_t run, subrun, gate_num, evt_num, truth_plane, pred_plane;
    
    // Columns 6 to 220 (inclusive) = 215 elements
    const int nProbs = 215; 
    Double_t probs[nProbs];

    // 4. Create Branches
    tree->Branch("run", &run, "run/I");
    tree->Branch("subrun", &subrun, "subrun/I");
    tree->Branch("gate_num", &gate_num, "gate_num/I");
    tree->Branch("evt_num", &evt_num, "evt_num/I");
    tree->Branch("truth_plane", &truth_plane, "truth_plane/I");
    tree->Branch("pred_plane", &pred_plane, "pred_plane/I");
    
    // Array branch: formatting is "name[length]/Type"
    tree->Branch("probs", probs, TString::Format("probs[%d]/D", nProbs));

    std::cout << "Processing file: " << inputPath << "..." << std::endl;

    // 5. Loop over lines
    std::string line;
    std::string cell;
    long long lineCount = 0;

    while (std::getline(infile, line)) {
        std::stringstream lineStream(line);
        std::vector<std::string> row;
        
        // Fast parsing of comma separated line
        while(std::getline(lineStream, cell, ',')) {
            row.push_back(cell);
        }

        // Safety check for empty lines or malformed rows
        // We expect at least 221 columns (0-220)
        if (row.size() < 221) {
            continue; 
        }

        try {
            // Fill Integers (Cols 0-5)
            run         = std::stoi(row[0]);
            subrun      = std::stoi(row[1]);
            gate_num    = std::stoi(row[2]);
            evt_num     = std::stoi(row[3]);
            truth_plane = std::stoi(row[4]);
            pred_plane  = std::stoi(row[5]);

            // Fill Array (Cols 6-220)
            for (int i = 0; i < nProbs; ++i) {
                // row index is i + 6
                probs[i] = std::stod(row[i + 6]);
            }

            tree->Fill();
            lineCount++;

            // Progress printout every 10k events
            if (lineCount % 10000 == 0) {
                std::cout << "\rProcessed " << lineCount << " events..." << std::flush;
            }

        } catch (const std::exception& e) {
            std::cerr << "\nError parsing line " << lineCount + 1 << ": " << e.what() << std::endl;
        }
    }

    // 6. Write and Close
    std::cout << "\nFinalizing... Total events: " << lineCount << std::endl;
    tree->Write();
    outFile->Close();
    
    delete outFile; // Clean up memory (tree is owned by file)
    std::cout << "Done. Output saved to: " << outputPath << std::endl;
}