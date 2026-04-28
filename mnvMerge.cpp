#include <iostream>
#include <string>
#include <vector>
#include <set>
#include <filesystem>
#include <stdexcept>
#include <cmath>

#include <TFile.h>
#include <TKey.h>
#include <TClass.h>
#include <TDirectory.h>

#include "PlotUtils/MnvH1D.h"
#include "PlotUtils/MnvH2D.h"

namespace fs = std::filesystem;
using namespace PlotUtils;

// Helper to parse the comma-separated exempt list
std::set<std::string> ParseExemptList(const std::string& str) {
    std::set<std::string> result;
    size_t start = 0;
    size_t end = str.find(',');
    while (end != std::string::npos) {
        result.insert(str.substr(start, end - start));
        start = end + 1;
        end = str.find(',', start);
    }
    result.insert(str.substr(start));
    return result;
}

// Helper to copy a TDirectory structure when missing from the other file
void CopyDirectory(TDirectory* source, TDirectory* dest) {
    TIter next(source->GetListOfKeys());
    TKey* key;
    while ((key = (TKey*)next())) {
        TObject* obj = key->ReadObj();
        if (obj->IsA()->InheritsFrom(TDirectory::Class())) {
            dest->mkdir(obj->GetName());
            dest->cd(obj->GetName());
            CopyDirectory((TDirectory*)obj, gDirectory);
        } else {
            dest->cd();
            obj->Write();
        }
        delete obj;
    }
}

// Templated merge function to handle both MnvH1D and MnvH2D
template <typename T>
T* MergeMnvHistograms(T* h1, T* h2, bool overwrite, const std::set<std::string>& exempt) {
    // 1. Check Central Value consistency using a tolerance
    double tolerance = 1e-9;
    if (h1->GetNcells() != h2->GetNcells()) {
        throw std::runtime_error(std::string("Bin counts differ for ") + h1->GetName());
    }
    for (int i = 0; i < h1->GetNcells(); ++i) {
        if (std::abs(h1->GetBinContent(i) - h2->GetBinContent(i)) > tolerance) {
            throw std::runtime_error(std::string("Central values differ outside tolerance for ") + h1->GetName() + " at bin " + std::to_string(i));
        }
    }

    // Clone h1 as the base for the output
    T* hOut = (T*)h1->Clone();

    // 2. Vertical Error Bands
    for (const auto& name : h2->GetVertErrorBandNames()) {
        bool hasBand = hOut->HasVertErrorBand(name);
        bool isExempt = exempt.count(name) > 0;
        
        if (hasBand && overwrite && !isExempt) {
            delete hOut->PopVertErrorBand(name);
            hasBand = false;
        } else if (hasBand) {
            continue; // Skip, keep file 1's band
        }
        if (!hasBand) {
            hOut->PushErrorBand(name, new PlotUtils::MnvVertErrorBand(*h2->GetVertErrorBand(name)));
        }
    }

    // 3. Lateral Error Bands
    for (const auto& name : h2->GetLatErrorBandNames()) {
        bool hasBand = hOut->HasLatErrorBand(name);
        bool isExempt = exempt.count(name) > 0;
        
        if (hasBand && overwrite && !isExempt) {
            delete hOut->PopLatErrorBand(name);
            hasBand = false;
        } else if (hasBand) {
            continue;
        }
        if (!hasBand) {
            hOut->PushErrorBand(name, new PlotUtils::MnvLatErrorBand(*h2->GetLatErrorBand(name)));
        }
    }

    // 4. Uncorrelated Errors
    for (const auto& name : h2->GetUncorrErrorNames()) {
        bool hasBand = hOut->HasUncorrError(name);
        bool isExempt = exempt.count(name) > 0;
        
        if (hasBand && overwrite && !isExempt) {
            delete hOut->PopUncorrError(name);
            hasBand = false;
        } else if (hasBand) {
            continue;
        }
        if (!hasBand) {
            hOut->PushUncorrError(name, new TH1D(*h2->GetUncorrError(name)));
        }
    }

    // 5. Covariance Matrices
    for (const auto& name : h2->GetCovMatricesNames()) {
        bool hasBand = hOut->HasErrorMatrix(name);
        bool isExempt = exempt.count(name) > 0;
        
        if (hasBand && overwrite && !isExempt) {
            delete hOut->PopSysErrorMatrix(name);
            hasBand = false;
        } else if (hasBand) {
            continue;
        }
        if (!hasBand) {
            hOut->PushCovMatrix(name, h2->GetSysErrorMatrix(name));
        }
    }

    return hOut;
}

// Recursively traverse TDirectories and merge matched objects
void MergeRootDirectories(TDirectory* dir1, TDirectory* dir2, TDirectory* dirOut, bool overwrite, const std::set<std::string>& exempt) {
    TIter next1(dir1->GetListOfKeys());
    TKey* key1;
    std::set<std::string> processedKeys;

    // Process file 1 elements (and merge if present in file 2)
    while ((key1 = (TKey*)next1())) {
        std::string name = key1->GetName();
        processedKeys.insert(name);

        TObject* obj1 = key1->ReadObj();
        TObject* obj2 = dir2->Get(name.c_str());

        if (obj1->IsA()->InheritsFrom(TDirectory::Class())) {
            dirOut->mkdir(name.c_str());
            dirOut->cd(name.c_str());
            TDirectory* newDirOut = gDirectory;
            if (obj2 && obj2->IsA()->InheritsFrom(TDirectory::Class())) {
                MergeRootDirectories((TDirectory*)obj1, (TDirectory*)obj2, newDirOut, overwrite, exempt);
            } else {
                CopyDirectory((TDirectory*)obj1, newDirOut);
            }
        } 
        else if (obj1->IsA()->InheritsFrom(PlotUtils::MnvH1D::Class())) {
            dirOut->cd();
            if (obj2 && obj2->IsA()->InheritsFrom(PlotUtils::MnvH1D::Class())) {
                MnvH1D* hOut = MergeMnvHistograms((MnvH1D*)obj1, (MnvH1D*)obj2, overwrite, exempt);
                hOut->Write(name.c_str(), TObject::kOverwrite);
                delete hOut;
            } else {
                obj1->Write();
            }
        } 
        else if (obj1->IsA()->InheritsFrom(PlotUtils::MnvH2D::Class())) {
            dirOut->cd();
            if (obj2 && obj2->IsA()->InheritsFrom(PlotUtils::MnvH2D::Class())) {
                MnvH2D* hOut = MergeMnvHistograms((MnvH2D*)obj1, (MnvH2D*)obj2, overwrite, exempt);
                hOut->Write(name.c_str(), TObject::kOverwrite);
                delete hOut;
            } else {
                obj1->Write();
            }
        } 
        else {
            dirOut->cd();
            obj1->Write();
        }
        delete obj1;
        if (obj2) delete obj2;
    }

    // Process file 2 elements entirely missing from file 1
    TIter next2(dir2->GetListOfKeys());
    TKey* key2;
    while ((key2 = (TKey*)next2())) {
        std::string name = key2->GetName();
        if (processedKeys.count(name) > 0) continue; 

        TObject* obj2 = key2->ReadObj();
        if (obj2->IsA()->InheritsFrom(TDirectory::Class())) {
            dirOut->mkdir(name.c_str());
            dirOut->cd(name.c_str());
            CopyDirectory((TDirectory*)obj2, gDirectory);
        } else {
            dirOut->cd();
            obj2->Write();
        }
        delete obj2;
    }
}

// Orchestrates merging of two discrete ROOT files
void ProcessRootFile(const fs::path& p1, const fs::path& p2, const fs::path& pOut, bool overwrite, const std::set<std::string>& exempt) {
    std::cout << "Merging: " << p1 << " and " << p2 << " -> " << pOut << std::endl;
    
    TFile* f1 = TFile::Open(p1.c_str(), "READ");
    TFile* f2 = TFile::Open(p2.c_str(), "READ");
    if (!f1 || f1->IsZombie() || !f2 || f2->IsZombie()) {
        throw std::runtime_error("Failed to open input files.");
    }

    fs::create_directories(pOut.parent_path());
    TFile* fOut = TFile::Open(pOut.c_str(), "RECREATE");
    
    MergeRootDirectories(f1, f2, fOut, overwrite, exempt);

    fOut->Close();
    f1->Close();
    f2->Close();
    delete fOut; delete f1; delete f2;
}

int main(int argc, char** argv) {
    bool overwrite = false;
    std::set<std::string> exemptBands;
    std::string path1, path2, outPath;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--overwrite") {
            overwrite = true;
        } else if (arg == "--exempt") {
            if (i + 1 < argc) {
                exemptBands = ParseExemptList(argv[++i]);
            } else {
                std::cerr << "Error: --exempt requires a comma-separated list." << std::endl;
                return 1;
            }
        } else if (path1.empty()) {
            path1 = arg;
        } else if (path2.empty()) {
            path2 = arg;
        } else if (outPath.empty()) {
            outPath = arg;
        } else {
            std::cerr << "Unrecognized trailing argument: " << arg << std::endl;
            return 1;
        }
    }

    if (path1.empty() || path2.empty() || outPath.empty()) {
        std::cerr << "Usage: " << argv[0] << " [options] <input1> <input2> <output>\n"
                  << "Options:\n"
                  << "  --overwrite              Overwrite error bands with the same name (File 2 overrides File 1).\n"
                  << "  --exempt <band1,band2>   Comma-separated list of bands exempt from overwriting.\n\n"
                  << "Inputs can be two ROOT files or two directories." << std::endl;
        return 1;
    }

    fs::path p1(path1), p2(path2), pOut(outPath);

    try {
        if (fs::is_regular_file(p1) && fs::is_regular_file(p2)) {
            // File to File merging
            ProcessRootFile(p1, p2, pOut, overwrite, exemptBands);
        } 
        else if (fs::is_directory(p1) && fs::is_directory(p2)) {
            // Directory to Directory merging
            fs::create_directories(pOut);
            
            // Loop through dir 1
            for (const auto& entry : fs::recursive_directory_iterator(p1)) {
                if (entry.is_regular_file() && entry.path().extension() == ".root") {
                    fs::path relPath = fs::relative(entry.path(), p1);
                    fs::path p2Target = p2 / relPath;
                    fs::path pOutTarget = pOut / relPath;

                    if (fs::exists(p2Target)) {
                        ProcessRootFile(entry.path(), p2Target, pOutTarget, overwrite, exemptBands);
                    } else {
                        std::cout << "Copying (only in dir1): " << entry.path() << " -> " << pOutTarget << std::endl;
                        fs::create_directories(pOutTarget.parent_path());
                        fs::copy(entry.path(), pOutTarget, fs::copy_options::overwrite_existing);
                    }
                }
            }

            // Loop through dir 2 to catch files that weren't in dir 1
            for (const auto& entry : fs::recursive_directory_iterator(p2)) {
                if (entry.is_regular_file() && entry.path().extension() == ".root") {
                    fs::path relPath = fs::relative(entry.path(), p2);
                    fs::path p1Target = p1 / relPath;
                    if (!fs::exists(p1Target)) { // Handled matched ones in the previous loop
                        fs::path pOutTarget = pOut / relPath;
                        std::cout << "Copying (only in dir2): " << entry.path() << " -> " << pOutTarget << std::endl;
                        fs::create_directories(pOutTarget.parent_path());
                        fs::copy(entry.path(), pOutTarget, fs::copy_options::overwrite_existing);
                    }
                }
            }
        } 
        else {
            std::cerr << "Error: Inputs must both be regular files or both be directories." << std::endl;
            return 1;
        }
    } catch (const std::exception& e) {
        std::cerr << "Fatal error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}