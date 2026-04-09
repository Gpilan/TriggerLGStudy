#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <iostream>

void check_electron() {
    TFile* file = TFile::Open("electron_test_0.root");
    if (!file || file->IsZombie()) {
        std::cout << "Error: Cannot open ROOT file" << std::endl;
        return;
    }
    
    std::cout << "=== ROOT File Contents ===" << std::endl;
    file->ls();
    
    TTree* tree = (TTree*)file->Get("CBDsim");
    if (tree) {
        std::cout << "\n=== Tree Information ===" << std::endl;
        std::cout << "Number of entries: " << tree->GetEntries() << std::endl;
        
        if (tree->GetEntries() > 0) {
            std::cout << "\n=== First Entry Data ===" << std::endl;
            tree->Show(0);
            
            // 히스토그램 생성
            TH1F* h_energy = new TH1F("h_energy", "Energy Distribution;Energy (MeV);Counts", 50, 0, 100);
            tree->Draw("Edeps.Edep>>h_energy", "", "goff");
            
            TCanvas* c = new TCanvas("c", "Electron Test Results", 800, 600);
            h_energy->Draw();
            c->SaveAs("electron_energy_plot.png");
            std::cout << "Plot saved as electron_energy_plot.png" << std::endl;
        }
    } else {
        std::cout << "Error: Cannot find CBDsim tree" << std::endl;
    }
    
    file->Close();
}
