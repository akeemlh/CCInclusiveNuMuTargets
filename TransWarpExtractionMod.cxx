#include "TransWarpExtractionMod.h" 
#include "MinervaUnfold/MnvUnfold.h"
#include "RooUnfold/RooUnfold.h"
#ifndef MNVROOT6
#ifndef NCINTEX
#include "Cintex/Cintex.h"
#endif
#endif

#include "PlotUtils/TargetUtils.h"
#include "PlotUtils/MnvPlotter.h"
#include "PlotUtils/MnvH1D.h"
#include "PlotUtils/MnvH2D.h"
#include "PlotUtils/MnvH3D.h"

#include "TFile.h"
#include "TDirectory.h"
#include "TRandom3.h"
#include "TParameter.h"
#include "TProfile.h"

#include <iostream>
#include <numeric>
#include <getopt.h>
#include <string>
#include <cassert>
#include <algorithm>
#include <iterator>
#include <stdlib.h>// or cstdlib is c++

#include "util/GetIngredient.h"
#include "util/NukeUtils.h"

#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/Minimizer.h"

//////////////////////////////////////////////////////////////////////////
// TransWarpExtraction.cxx
// Dan Ruterbories wrote this
// Aaron Bercellie editted it and came up with the Star Trek reference
//////////////////////////////////////////////////////////////////////////
void split(const std::string& s,std::vector<std::string>& v, char delim) {
    int i = 0;
    std::size_t pos = s.find(delim);
    if(pos == std::string::npos) v.push_back(s.substr(i, s.length()));//For the case there is no delim
    while (pos != std::string::npos) {
      v.push_back(s.substr(i, pos-i));
      i = ++pos;
      pos = s.find(delim, pos);

      if (pos == std::string::npos)
         v.push_back(s.substr(i, s.length()));
    }
}

template <class MnvH>
struct Chi2Calculator {
    MnvH* dataHist;
    MnvH* mcHist;

    Chi2Calculator(MnvH* data, MnvH* mc) : dataHist(data), mcHist(mc) {}

    double operator()(const double* val) const {
      double scale  = val[0];
      int ndof;
      PlotUtils::MnvPlotter plotter;
      plotter.ApplyStyle(PlotUtils::kCCQENuInclusiveStyle);
      auto tmpHist = dataHist->Clone("tmp_minimization_hist");
      double chiSq = plotter.Chi2DataMC(tmpHist, mcHist, ndof, scale, true, true);
      delete tmpHist; 
      return chiSq;
    }
};

namespace PlotUtils
{
  template <class MnvH> TransWarpExtraction<MnvH>::TransWarpExtraction( std::string data_file, std::string data_name, 
                                                                        std::string data_truth_file, std::string data_truth_name, 
                                                                        std::string reco_file, std::string reco_name, 
                                                                        std::string truth_file, std::string truth_name, 
                                                                        std::string migration_file, std::string migration_name,  
                                                                        std::vector<int> iterations, std::vector<int> exclude_chi2_bins,
                                                                        int random_seed, int stat_universes, bool bIterLogScale, std::string targetStr, std::string variable, int PDG)
  {
    bookHistos( m_data       , data_file, data_name, targetStr, variable, PDG);
    bookHistos( m_data_truth , data_truth_file, data_truth_name, targetStr, variable, PDG);
    bookHistos( m_reco       , reco_file, reco_name, targetStr, variable, PDG);
    bookHistos( m_truth      , truth_file, truth_name, targetStr, variable, PDG);
    bookHistos( m_migration  , migration_file, migration_name, targetStr, variable, PDG);
    m_migration_original = m_migration->Clone("Original_Migration");
    //    m_reco_original = m_reco->Clone("Original_Reco");
    //    m_truth_original = m_truth->Clone("Original_Truth");
    for(uint it = 0; it<iterations.size(); ++it)
    {
       m_avg_unfoldingCovMatrices[it] = NULL;
       m_avg_data_unfolded[it]        = NULL;
    }

    m_chi2_md_md_iter_chi2 = NULL;
    m_chi2_md_td_iter_chi2 = NULL;
    m_chi2_md_tmc_iter_chi2 = NULL;
    m_chi2_td_tmc_iter_chi2 = NULL;

    m_chi2_md_md_iter_stat = NULL;
    m_chi2_md_td_iter_stat = NULL;
    m_chi2_md_tmc_iter_stat = NULL;
    m_chi2_td_tmc_iter_stat = NULL;

    plotter = new MnvPlotter();
    myrandom = new TRandom3(42);
    if(random_seed >= 0 ) myrandom->SetSeed(random_seed);

    doTailAvgThrows = false;
    m_iterations = iterations;
    m_exclude_chi2_bins = exclude_chi2_bins;
    m_nStatUniverses = stat_universes+1;
    m_ChiSquareGuess = 5000;
    m_ChiSquareStep  = 50;

/*    for(uint bin=0; bin<m_exclude_chi2_bins.size(); ++bin)
    {
      m_data_truth->SetBinContent(m_exclude_chi2_bins[bin],0);
      m_truth->SetBinContent(m_exclude_chi2_bins[bin],0);
      m_data_truth->SetBinError(m_exclude_chi2_bins[bin],0);
      m_truth->SetBinError(m_exclude_chi2_bins[bin],0);

    } */ //commented out so that closure test would pass   
    
    for(uint it = 0; it < m_iterations.size(); ++it) m_iterations_bins.push_back((double)m_iterations[it]);
    if( bIterLogScale ) m_iterations_bins.push_back(m_iterations[m_iterations.size()-1]*10);
    else
    {
      if(m_iterations.size()>1) m_iterations_bins.push_back( 2*m_iterations[m_iterations.size()-1]-m_iterations[m_iterations.size()-2] );
      else m_iterations_bins.push_back( m_iterations[m_iterations.size()-1]+1 );
    }
  }
  template <class MnvH> TransWarpExtraction<MnvH>::~TransWarpExtraction()
  {
    delete plotter;
    delete myrandom;

    delete m_data_truth;
    delete m_data;
    delete m_reco;
    delete m_truth;
    delete m_migration;
    
    m_data_pull.clear();

    m_avg_data_error.clear();      
    m_avg_data_pull.clear();
    m_bin_data_errors.clear();
    m_bin_data_pulls.clear();
    m_bin_avg_data_errors.clear();
    m_bin_avg_data_pulls.clear();

    m_ratio_md_td.clear();
    m_ratio_md_input.clear();
    m_bin_iter_ratio_md_td.clear();
    m_bin_iter_ratio_md_input.clear();
    m_bin_iter_avg_ratio_md_td.clear();
    m_bin_iter_avg_ratio_md_input.clear();

    m_avg_unfoldingCovMatrices.clear();
    m_avg_data_unfolded.clear();

    m_unfoldingCovMatrices.clear();
    m_outputresults.clear();
    m_inputresults.clear();
    m_valmat_vec.clear();

    m_bin_chi2_iter.clear();
    m_avg_bin_chi2_iter.clear();
    
    //delete m_avg_chi2_md_md_iter_chi2;
    //delete m_avg_chi2_md_td_iter_chi2;
    //delete m_avg_chi2_md_tmc_iter_chi2;
    //delete m_avg_chi2_td_tmc_iter_chi2;

    //delete m_chi2_md_md_iter_chi2;
    //delete m_chi2_md_td_iter_chi2;
    //delete m_chi2_md_tmc_iter_chi2;
    //delete m_chi2_td_tmc_iter_chi2;

    //delete m_chi2_md_md_iter_stat;
    //delete m_chi2_md_td_iter_stat;
    //delete m_chi2_md_tmc_iter_stat;
    //delete m_chi2_td_tmc_iter_stat;

    //delete m_avg_chi2_md_td_iter_chi2_trunc;
    //delete m_chi2_md_td_iter_chi2_trunc;
  }
  template <class MnvH> void TransWarpExtraction<MnvH>::Setup2DChi2Hists( )
  {
    uint nIterBins = m_iterations_bins.size();
    Double_t iter_bins[nIterBins];
    for(uint iter = 0; iter< nIterBins; ++iter){
      iter_bins[iter]=m_iterations_bins[iter];
    }

    uint nChiBins = (int)m_ChiSquareGuess/m_ChiSquareStep;
    Double_t chi_bins[nChiBins+1];
    for(uint chi = 0; chi <= nChiBins; ++chi) chi_bins[chi]=(double)chi*m_ChiSquareStep;

    Double_t stat_bins[m_nStatUniverses+1];
    for(uint su = 0; su <= m_nStatUniverses; ++su) stat_bins[su]=su;

    int ndf = m_data->GetNbinsX()*m_data->GetNbinsY()*m_data->GetNbinsZ();
 
    m_avg_chi2_md_md_iter_chi2  = new TProfile("m_avg_chi2_modelData_modelData_iter_chi2",Form(";#bf{Averaged (Unfolded Data:Unfolded Data)} # of Iterations;#Chi^{2}(ndf=%i)", ndf), nIterBins-1, iter_bins);
    m_avg_chi2_md_td_iter_chi2  = new TProfile("m_avg_chi2_modelData_trueData_iter_chi2" ,Form(";#bf{Averaged (Unfolded Data:True Data)} # of Iterations;#Chi^{2}(ndf=%i)"    , ndf), nIterBins-1, iter_bins);
    m_avg_chi2_md_tmc_iter_chi2 = new TProfile("m_avg_chi2_modelData_trueMC_iter_chi2"   ,Form(";#bf{Averaged (Unfolded Data:True MC)} # of Iterations;#Chi^{2}(ndf=%i)"      , ndf), nIterBins-1, iter_bins);
    m_avg_chi2_td_tmc_iter_chi2 = new TProfile("m_avg_chi2_trueData_trueMC_iter_chi2"    ,Form(";#bf{Averaged (True Data:True MC)} # of Iterations;#Chi^{2}(ndf=%i)"          , ndf), nIterBins-1, iter_bins);

    //Model Data - Model Data
    m_chi2_md_md_iter_chi2  = new MnvH2D("h_chi2_modelData_modelData_iter_chi2",Form(";#bf{(Unfolded Data:Unfolded Data)} # of Iterations;#Chi^{2}(ndf=%i)", ndf), nIterBins-1, iter_bins, nChiBins, chi_bins);
    m_chi2_md_md_iter_stat  = new MnvH2D("h_chi2_modelData_modelData_iter_stat",";#bf{(Unfolded Data:Unfolded Data)} # of Iterations;Stat Universe #", nIterBins-1, iter_bins, m_nStatUniverses, stat_bins);
    //Model Data - True Data
    m_chi2_md_td_iter_chi2  = new MnvH2D("h_chi2_modelData_trueData_iter_chi2",Form(";#bf{(Unfolded Data:True Data)} # of Iterations;#Chi^{2}(ndf=%i)", ndf),  nIterBins-1, iter_bins, nChiBins, chi_bins);
    m_chi2_md_td_iter_stat  = new MnvH2D("h_chi2_modelData_trueData_iter_stat",";#bf{(Unfolded Data:True Data)} # of Iterations;Stat Universe #",  nIterBins-1, iter_bins, m_nStatUniverses, stat_bins);
    //Model Data - True MC 
    m_chi2_md_tmc_iter_chi2  = new MnvH2D("h_chi2_modelData_trueMC_iter_chi2",Form(";#bf{(Unfolded Data:True MC)} # of Iterations;#Chi^{2}(ndf=%i)", ndf),  nIterBins-1, iter_bins, nChiBins, chi_bins);
    m_chi2_md_tmc_iter_stat  = new MnvH2D("h_chi2_modelData_trueMC_iter_stat",";#bf{(Unfolded Data:True MC)} # of Iterations;Stat Universe #",  nIterBins-1, iter_bins, m_nStatUniverses, stat_bins);
    //True Data - True MC
    m_chi2_td_tmc_iter_chi2  = new MnvH2D("h_chi2_trueData_trueMC_iter_chi2",Form(";#bf{(True Data:True MC)} # of Iterations;#Chi^{2}(ndf=%i)", ndf), nIterBins-1, iter_bins, nChiBins, chi_bins); 
    m_chi2_td_tmc_iter_stat  = new MnvH2D("h_chi2_trueData_trueMC_iter_stat",";#bf{(True Data:True MC)} # of Iterations;Stat Universe #", nIterBins-1, iter_bins, m_nStatUniverses, stat_bins); 
  }
  
  template <class MnvH> void TransWarpExtraction<MnvH>::ClearInputErrorBands()
  {
    m_data->ClearAllErrorBands();
    m_migration->ClearAllErrorBands();
    m_reco->ClearAllErrorBands();
    m_truth->ClearAllErrorBands();
  }
  template <class MnvH> void TransWarpExtraction<MnvH>::ScaleDataPOT( double data_pot_norm ) 
  {
    m_data->Scale(data_pot_norm);      
    m_data_truth->Scale(data_pot_norm);
    m_migration->Scale(data_pot_norm);
    m_reco->Scale(data_pot_norm);
    m_truth->Scale(data_pot_norm);
  }
  template <class MnvH> void TransWarpExtraction<MnvH>::ScalePOT( double pot_scale ) 
  {
    m_migration->Scale(pot_scale);
    m_reco->Scale(pot_scale);
    m_truth->Scale(pot_scale);
  }
  template <class MnvH> void TransWarpExtraction<MnvH>::ScaleMigStatUnc( double stat_scale )
  {
    int nbins_x = m_migration->GetNbinsX()+2;
    int nbins_y = m_migration->GetNbinsY()+2;
    
    for(int i=0;i<nbins_x;i++){
      for(int j=0;j<nbins_y;j++){
	double orig_unc = m_migration->GetBinError(i,j);
	m_migration->SetBinError(i,j,orig_unc*stat_scale);
      }
    }

  }
  template <class MnvH> void TransWarpExtraction<MnvH>::WriteOutput( std::string output_filename, bool b_WriteChi2Verbose, bool b_WriteLinearHists ) 
  {
  //==================================================================
  // Create ROOT file to store histograms
  // Write to file all the created histograms 
  //==================================================================
    TFile *f_output = TFile::Open( output_filename.c_str(), "RECREATE");
    f_output->cd();
    std::cout<< "--------------------------------------------"  << std::endl;
    std::cout<< "Writing histograms to file"  << std::endl;
    std::cout<< "--------------------------------------------"  << std::endl;

    TDirectory *input_dir = f_output->mkdir("Input_Hists");
    input_dir->cd();

    m_data       -> Write("h_data");
    m_data_truth -> Write("h_data_truth");
    m_reco       -> Write("h_mc_reco");
    m_truth      -> Write("h_mc_truth");
    m_migration  -> Write("h_migration_matrix");

    if(b_WriteLinearHists){
      MnvH1D* tmp_data_Linear       = LinearizeHist( m_data       );  
      MnvH1D* tmp_data_truth_Linear = LinearizeHist( m_data_truth );
      MnvH1D* tmp_reco_Linear       = LinearizeHist( m_reco       );
      MnvH1D* tmp_truth_Linear      = LinearizeHist( m_truth      );

      tmp_data_Linear       -> SetDirectory(input_dir);
      tmp_data_truth_Linear -> SetDirectory(input_dir);
      tmp_reco_Linear       -> SetDirectory(input_dir);
      tmp_truth_Linear      -> SetDirectory(input_dir);

      tmp_data_Linear       -> Write("h_data_Linear");
      tmp_data_truth_Linear -> Write("h_data_truth_Linear");
      tmp_reco_Linear       -> Write("h_mc_reco_Linear");
      tmp_truth_Linear      -> Write("h_mc_truth_Linear");
    }
    f_output->cd();

    TDirectory *out_dir = f_output->mkdir("Unfolded_Data");
    out_dir->cd();
    for(std::vector<int>::iterator num_iter = m_iterations.begin(); num_iter != m_iterations.end(); ++num_iter)
    {
      for(uint i=0;i<m_outputresults[*num_iter].size();i++){
        if(m_outputresults[*num_iter][i]){
          m_outputresults[*num_iter][i]->SetDirectory(out_dir);
          m_outputresults[*num_iter][i]->Write();
          if(b_WriteLinearHists){
            MnvH1D* tmp_outputresult = LinearizeHist( m_outputresults[*num_iter][i] );
            tmp_outputresult->SetDirectory(out_dir);
            tmp_outputresult->Write();
          }
        }
      }
      m_avg_data_unfolded[*num_iter]->SetDirectory(out_dir);
      m_avg_data_unfolded[*num_iter]->Write();
      if(b_WriteLinearHists){
        MnvH1D* tmp_avg_data_unfolded = LinearizeHist(m_avg_data_unfolded[*num_iter]);
        tmp_avg_data_unfolded->SetDirectory(out_dir);
        tmp_avg_data_unfolded->Write();
      }
    }
    f_output->cd();

    TDirectory *in_dir = f_output->mkdir("Stat_Varied_Smeared_Data");
    in_dir->cd();
    for(std::vector<int>::iterator num_iter = m_iterations.begin(); num_iter != m_iterations.end(); ++num_iter)
    {
      for(uint i=0;i<m_inputresults[*num_iter].size();i++){
        if(m_inputresults[*num_iter][i]){
          m_inputresults[*num_iter][i]->SetDirectory(in_dir);
          m_inputresults[*num_iter][i]->Write();
          if(b_WriteLinearHists){
            MnvH1D* tmp_inputresult = LinearizeHist( m_inputresults[*num_iter][i] );
            tmp_inputresult->SetDirectory(in_dir);
            tmp_inputresult->Write();
          }
        }
      }
    }
    f_output->cd();

    TDirectory *chi_dir = f_output->mkdir("Chi2Maps");
    chi_dir->cd();
    for(std::vector<int>::iterator num_iter = m_iterations.begin(); num_iter != m_iterations.end(); ++num_iter)
    {
     //m_avg_valmat_vec[*num_iter]->Write(Form("AvgChi2Map_Iter_%d",*num_iter));
      for(uint i=0;i<m_valmat_vec[*num_iter].size();i++){
        m_valmat_vec[*num_iter][i]->Write(Form("Chi2Map_Stat_%d_Iter_%d",i,*num_iter));
      }
    }
    f_output->cd();

    TDirectory *cov_dir = f_output->mkdir("Cov_Matrix");
    cov_dir->cd();
    for(std::vector<int>::iterator num_iter = m_iterations.begin(); num_iter != m_iterations.end(); ++num_iter)
    {
      m_avg_unfoldingCovMatrices[*num_iter]->Write(Form("AvgCovMatrix_Iter_%d",*num_iter));
      for(uint i=0;i<m_unfoldingCovMatrices[*num_iter].size();i++){
        m_unfoldingCovMatrices[*num_iter][i].Write(Form("CovMatrix_Stat_%d_Iter_%d",i,*num_iter));
      }
    }
    f_output->cd();
    
    TDirectory *pull_dir_bin = f_output->mkdir("Bin_Pull_Histograms");
    pull_dir_bin->cd();
    for(std::map<int,MnvH2D*>::iterator it=m_bin_data_pulls.begin(); it!=m_bin_data_pulls.end(); ++it)
    {
      if(it->second){
        it->second->SetDirectory(pull_dir_bin);
        it->second->Write();
      }
    }
    for(std::map<int,MnvH2D*>::iterator it=m_bin_data_errors.begin(); it!=m_bin_data_errors.end(); ++it)
    {
      if(it->second){
        it->second->SetDirectory(pull_dir_bin);
        it->second->Write();
      }
    }
    for(std::map<int,TProfile*>::iterator it=m_bin_avg_data_pulls.begin(); it!=m_bin_avg_data_pulls.end(); ++it)
    {
      if(it->second){
        it->second->SetDirectory(pull_dir_bin);
        it->second->Write();
      }
    }
    for(std::map<int,TProfile*>::iterator it=m_bin_avg_data_errors.begin(); it!=m_bin_avg_data_errors.end(); ++it)
    {
      if(it->second){
        it->second->SetDirectory(pull_dir_bin);
        it->second->Write();
      }
    }
    f_output->cd();

    TDirectory *pull_dir = f_output->mkdir("Pull_Histograms");
    pull_dir->cd();
    for(std::vector<int>::iterator num_iter = m_iterations.begin(); num_iter != m_iterations.end(); ++num_iter)
    {
      for(uint i=0;i<m_data_pull[*num_iter].size();i++){
        if(m_data_pull[*num_iter][i]){
          m_data_pull[*num_iter][i]->SetDirectory(pull_dir);
          m_data_pull[*num_iter][i]->Write();
          if(b_WriteLinearHists){
            MnvH1D* tmp_data_pull = LinearizeHist( m_data_pull[*num_iter][i] );
            tmp_data_pull->SetDirectory(pull_dir);
            tmp_data_pull->Write();
          }
        }
      }
    }
    f_output->cd();

    TDirectory *pull_dir_avg = f_output->mkdir("Average_Pull_Histograms");
    pull_dir_avg->cd();
    for(std::vector<int>::iterator num_iter = m_iterations.begin(); num_iter != m_iterations.end(); ++num_iter)
    {
      if(m_avg_data_error[*num_iter]){
        m_avg_data_error[*num_iter]->SetDirectory(pull_dir_avg);
        m_avg_data_pull[*num_iter] ->SetDirectory(pull_dir_avg);
        m_avg_data_error[*num_iter]->Write();
        m_avg_data_pull[*num_iter] ->Write();
        if(b_WriteLinearHists)
        {
          MnvH1D* tmp_avg_data_error = LinearizeHist( m_avg_data_error[*num_iter] ); 
          MnvH1D* tmp_avg_data_pull  = LinearizeHist( m_avg_data_pull[*num_iter]  );

          tmp_avg_data_error->SetDirectory(pull_dir_avg); 
          tmp_avg_data_pull ->SetDirectory(pull_dir_avg); 
          tmp_avg_data_error->Write(); 
          tmp_avg_data_pull ->Write(); 
        }
      }
    }
    f_output->cd();

    TDirectory *ratio_ut_dir = f_output->mkdir("Ratio_Unfolded_True_Histograms");
    ratio_ut_dir->cd();
    for(std::vector<int>::iterator num_iter = m_iterations.begin(); num_iter != m_iterations.end(); ++num_iter)
    {
      for(uint i=0;i<m_ratio_md_td[*num_iter].size();i++){
        if(m_ratio_md_td[*num_iter][i]){
          m_ratio_md_td[*num_iter][i]->SetDirectory(ratio_ut_dir);
          m_ratio_md_td[*num_iter][i]->Write();
          if(b_WriteLinearHists){
            MnvH1D* tmp_ratio_md_td = LinearizeHist( m_ratio_md_td[*num_iter][i] );
            tmp_ratio_md_td->SetDirectory(ratio_ut_dir);
            tmp_ratio_md_td->Write();
          }
        }
      }
    }
    f_output->cd();

    TDirectory *ratio_ui_dir = f_output->mkdir("Ratio_Unfolded_Input_Histograms");
    ratio_ui_dir->cd();
    for(std::vector<int>::iterator num_iter = m_iterations.begin(); num_iter != m_iterations.end(); ++num_iter)
    {
      for(uint i=0;i<m_ratio_md_input[*num_iter].size();i++){
        if(m_ratio_md_input[*num_iter][i]){
          m_ratio_md_input[*num_iter][i]->SetDirectory(ratio_ui_dir);
          m_ratio_md_input[*num_iter][i]->Write();
          if(b_WriteLinearHists){
            MnvH1D* tmp_ratio_md_input = LinearizeHist( m_ratio_md_input[*num_iter][i] );
            tmp_ratio_md_input->SetDirectory(ratio_ui_dir);
            tmp_ratio_md_input->Write();
          }
        }
      }
    }
    f_output->cd();

    TDirectory *ratio_ut_bin_dir = f_output->mkdir("Bin_Ratio_Unfolded_True_Histograms");
    ratio_ut_bin_dir->cd();
    for(std::map<int,MnvH2D*>::iterator it=m_bin_iter_ratio_md_td.begin(); it!=m_bin_iter_ratio_md_td.end(); ++it)
    {
      if(it->second){
        it->second->SetDirectory(ratio_ut_bin_dir);
        it->second->Write();
      }
    }
    for(std::map<int,TProfile*>::iterator it=m_bin_iter_avg_ratio_md_td.begin(); it!=m_bin_iter_avg_ratio_md_td.end(); ++it)
    {
      if(it->second){
        it->second->SetDirectory(ratio_ut_bin_dir);
        it->second->Write();
      }
    }
    f_output->cd();

    TDirectory *ratio_ui_bin_dir = f_output->mkdir("Bin_Ratio_Unfolded_Input_Histograms");
    ratio_ui_bin_dir->cd();
    for(std::map<int,MnvH2D*>::iterator it=m_bin_iter_ratio_md_input.begin(); it!=m_bin_iter_ratio_md_input.end(); ++it)
    {
      if(it->second){
        it->second->SetDirectory(ratio_ui_bin_dir);
        it->second->Write();
      }
    }
    for(std::map<int,TProfile*>::iterator it=m_bin_iter_avg_ratio_md_input.begin(); it!=m_bin_iter_avg_ratio_md_input.end(); ++it)
    {
      if(it->second){
        it->second->SetDirectory(ratio_ui_bin_dir);
        it->second->Write();
      }
    }
    f_output->cd();

    TDirectory *bin_chi2_dists = f_output->mkdir("Bin_Chi2_Iteration_Dists");
    bin_chi2_dists->cd();
    for(std::map<int,MnvH2D*>::iterator it=m_bin_chi2_iter.begin(); it!=m_bin_chi2_iter.end(); ++it)
    {
      if(it->second){
        it->second->SetDirectory(bin_chi2_dists);
        it->second->Write();
      }
    }
    for(std::map<int,TProfile*>::iterator it=m_avg_bin_chi2_iter.begin(); it!=m_avg_bin_chi2_iter.end(); ++it)
    {
      if(it->second){
        it->second->SetDirectory(bin_chi2_dists);
        it->second->Write();
      }
    }
    f_output->cd();

    //Chi^2 metrics
    TDirectory *chi2_dists = f_output->mkdir("Chi2_Iteration_Dists");
    chi2_dists->cd();
    if(m_chi2_md_md_iter_chi2)
    {
      m_chi2_md_td_iter_chi2 ->SetDirectory(chi2_dists);
      m_chi2_md_tmc_iter_chi2->SetDirectory(chi2_dists);

      m_chi2_md_td_iter_stat ->SetDirectory(chi2_dists);
      m_chi2_md_tmc_iter_stat->SetDirectory(chi2_dists);

      m_avg_chi2_md_td_iter_chi2 ->SetDirectory(chi2_dists);
      m_avg_chi2_md_tmc_iter_chi2->SetDirectory(chi2_dists);
    
      m_chi2_md_td_iter_chi2 ->Write();
      m_chi2_md_tmc_iter_chi2->Write();

      m_chi2_md_td_iter_stat ->Write();
      m_chi2_md_tmc_iter_stat->Write();

      m_avg_chi2_md_td_iter_chi2 ->Write();
      m_avg_chi2_md_tmc_iter_chi2->Write();
    
      if(b_WriteChi2Verbose)
      {
        m_chi2_md_md_iter_chi2 ->SetDirectory(chi2_dists);
        m_chi2_td_tmc_iter_chi2->SetDirectory(chi2_dists);

        m_chi2_md_md_iter_stat ->SetDirectory(chi2_dists);
        m_chi2_td_tmc_iter_stat->SetDirectory(chi2_dists);
      
        m_avg_chi2_md_md_iter_chi2 ->SetDirectory(chi2_dists);
        m_avg_chi2_td_tmc_iter_chi2->SetDirectory(chi2_dists);

        m_chi2_md_md_iter_chi2 ->Write();
        m_chi2_td_tmc_iter_chi2->Write();

        m_chi2_md_md_iter_stat ->Write();
        m_chi2_td_tmc_iter_stat->Write();
      
        m_avg_chi2_md_md_iter_chi2 ->Write();
        m_avg_chi2_td_tmc_iter_chi2->Write();
      }
    }

    if(m_chi2_md_td_iter_chi2_trunc)
    {
      m_avg_chi2_md_td_iter_chi2_trunc ->SetDirectory(chi2_dists);
      m_chi2_md_td_iter_chi2_trunc     ->SetDirectory(chi2_dists);

      m_avg_chi2_md_td_iter_chi2_trunc ->Write();
      m_chi2_md_td_iter_chi2_trunc     ->Write();
    }

    if( m_chi2_md_td_iter_chi2_percentile )
    {
      m_median_chi2_md_td_iter_chi2    ->SetDirectory(chi2_dists);
      m_chi2_md_td_iter_chi2_percentile->SetDirectory(chi2_dists);

      m_median_chi2_md_td_iter_chi2    ->Write();
      m_chi2_md_td_iter_chi2_percentile->Write();
    }
    f_output->cd();
  
    f_output->Close();
    delete f_output;
  }
  template <class MnvH> void TransWarpExtraction<MnvH>::UnfoldStatUniverses( )
  {
    for(uint su = 0;su<m_nStatUniverses;su++){
      MnvH *stat_varied = su==0 ? (MnvH*)m_data->Clone("Stat-Input-0"):throwStat(m_data, myrandom ,su, "Stat");
      for(std::vector<int>::iterator num_iter = m_iterations.begin(); num_iter != m_iterations.end(); ++num_iter){
        MnvH* h_data_unfolded = NULL;
        if (*num_iter == 0) {
          //------------------------------------------------------------------
          // Special case code where there are no requested iterations
          // Copies the hist then changes the name to end in _unfold
          //------------------------------------------------------------------
          h_data_unfolded = new MnvH(*stat_varied);
          h_data_unfolded->SetName(Form( "%s_unfold", h_data_unfolded->GetName() ));
          if( !m_avg_data_unfolded[*num_iter] ){
             m_avg_data_unfolded[*num_iter] = new MnvH( h_data_unfolded->GetCVHistoWithStatError() );
             m_avg_data_unfolded[*num_iter]->SetName( Form("Avg_Iter%d", *num_iter) );
          }
          else
          {
            MnvH* h_tmp_unfold = new MnvH(h_data_unfolded->GetCVHistoWithStatError());
            m_avg_data_unfolded[*num_iter]->Add( h_tmp_unfold ); 
          }
          //!!!FIX THIS
          TMatrixD unfoldingCovMatrixOrig = h_data_unfolded->GetTotalErrorMatrix();

          int nCovCols = unfoldingCovMatrixOrig.GetNcols();
          int nCovRows = unfoldingCovMatrixOrig.GetNrows();
          if( !m_avg_unfoldingCovMatrices[*num_iter] ) m_avg_unfoldingCovMatrices[*num_iter] = new TMatrixD( nCovRows, nCovCols );
          TMatrixD* avg_unfoldingCovMatrix = m_avg_unfoldingCovMatrices[*num_iter]; 

          for( int row = 0; row < nCovRows; ++row)
          {
            for( int col = 0; col < nCovCols; ++col)
            {
              if(col==row) continue;
              avg_unfoldingCovMatrix[0][row][col] += (double)unfoldingCovMatrixOrig(row,col);
            }
          }
          m_avg_unfoldingCovMatrices[*num_iter] = avg_unfoldingCovMatrix; 
        } 
        else 
        {
          //Get a stat varied universe. Each bin is given a Poisson throw
          std::cout<< "---------------------------------------"  << std::endl;
          std::cout<< "Now unfolding statistical universe " << su << "/" << m_nStatUniverses-1 << " with " << *num_iter << " iterations"<< std::endl;
          std::cout<< "---------------------------------------"  << std::endl;
          bool data_unfolded2 = false;
          data_unfolded2 = m_fakes?
            UnfoldDataWithFakes( h_data_unfolded, stat_varied, m_reco, m_truth, m_migration, *num_iter ):
            UnfoldData( h_data_unfolded, stat_varied, m_reco, m_truth, m_migration, *num_iter );
          if( !data_unfolded2 )
          {
            std::cout << "Unfolding failed for either data or MC. Please check " << std::endl;
            return;
          }
        }
        m_outputresults[*num_iter].push_back((MnvH*)h_data_unfolded->Clone(Form("Stat_%d_Iter_%d",su,*num_iter)));
        m_inputresults[*num_iter].push_back((MnvH*)stat_varied->Clone(Form("Input_Stat_%d_Iter_%d",su,*num_iter)));
      }
    }
  }
  template <class MnvH> bool TransWarpExtraction<MnvH>::UnfoldDataWithFakes( MnvH* &h_data_unfolded, MnvH* stat_varied, MnvH* h_reco, MnvH* h_truth, MnvH2D* h_migration, int num_iter) {
    bool data_unfolded = false;
    TMatrixD CovMatrix;
    //unfold.setUseBetterStatErrorCalc(true);
    data_unfolded = unfold.UnfoldHistoWithFakes( h_data_unfolded, CovMatrix, h_migration, stat_varied, h_reco, h_truth, nullptr, num_iter, false,false);
    if( !data_unfolded)
      {
        std::cout << "Unfolding failed for either data or MC. Please check " << std::endl;
        return false;  
      }
    for (int i =0;i<CovMatrix.GetNrows();++i) {
      CovMatrix[i][i]=0;
      for(uint bin=0; bin<m_exclude_chi2_bins.size(); ++bin)
        {
          CovMatrix[i][m_exclude_chi2_bins[bin]] = 0;
          CovMatrix[m_exclude_chi2_bins[bin]][i] = 0;
        }
    }
    m_unfoldingCovMatrices[num_iter].push_back(CovMatrix);
    //h_data_unfolded->PushCovMatrix("unfoldingCov",CovMatrix);
    if( !m_avg_data_unfolded[num_iter] ) {
      m_avg_data_unfolded[num_iter] = new MnvH( h_data_unfolded->GetCVHistoWithStatError() );
      m_avg_data_unfolded[num_iter]->SetName( Form("Avg_Iter%d", num_iter) );
    }
    else {
      MnvH* h_tmp_unfold = new MnvH(h_data_unfolded->GetCVHistoWithStatError());
      m_avg_data_unfolded[num_iter]->Add( h_tmp_unfold );
    }
    int nCovCols = CovMatrix.GetNcols();
    int nCovRows = CovMatrix.GetNrows();
    if( !m_avg_unfoldingCovMatrices[num_iter] ) m_avg_unfoldingCovMatrices[num_iter] = new TMatrixD( nCovRows, nCovCols );
    TMatrixD* avg_unfoldingCovMatrix = m_avg_unfoldingCovMatrices[num_iter];
    for( int row = 0; row < nCovRows; ++row)
      {
        for( int col = 0; col < nCovCols; ++col)
          {
            avg_unfoldingCovMatrix[0][row][col] += (double)CovMatrix(row,col);
          }
      }
    m_avg_unfoldingCovMatrices[num_iter] = avg_unfoldingCovMatrix; 

    return true;

  } 

  template <class MnvH> bool TransWarpExtraction<MnvH>::UnfoldData( MnvH1D* &h_data_unfolded, MnvH1D* stat_varied, MnvH1D* h_reco, MnvH1D* h_truth, MnvH2D* h_migration, int num_iter ) 
  {
    bool data_unfolded = false;
    TMatrixD dummyCovMatrix;
    data_unfolded = unfold.UnfoldHisto( h_data_unfolded, dummyCovMatrix, h_migration, stat_varied, RooUnfold::kBayes, num_iter, false, false );
    
    if( !data_unfolded )
    {
      std::cout << "Unfolding failed for either data or MC. Please check " << std::endl;
      return false;  
    }

    TMatrixD unfoldingCovMatrixOrig = UnfoldDummy( stat_varied, h_data_unfolded, h_reco, h_truth, h_migration, num_iter ); 
    const int highBinX = h_data_unfolded->GetNbinsX() + 1;
    const int highBin  = h_data_unfolded->GetBin( highBinX );
    const int lowBin = 0;
    for( int i = lowBin; i <= highBin; ++i )
    {
      for(uint bin=0; bin<m_exclude_chi2_bins.size(); ++bin)
      {
        unfoldingCovMatrixOrig[i][m_exclude_chi2_bins[bin]] = 0;
        unfoldingCovMatrixOrig[m_exclude_chi2_bins[bin]][i] = 0;
      }
    }

    if( m_corr_factor > 0 ) {
      double uncmod  = (1. + 1./m_corr_factor); // To apply to the covariance matrix
      double sqrtmod = sqrt(uncmod);            // To apply to the diagonal error

      // Scale covariance matrix
      unfoldingCovMatrixOrig *= uncmod;

      // Scale diagonal errors
      for ( int i = 0; i < h_data_unfolded->GetNbinsX()+2; ++i )
	h_data_unfolded -> SetBinError(i, h_data_unfolded->GetBinError(i) * sqrtmod);
    }

    m_unfoldingCovMatrices[num_iter].push_back(unfoldingCovMatrixOrig);
    h_data_unfolded->PushCovMatrix("unfoldingCov",unfoldingCovMatrixOrig);
    
    //Getting hists and covariance matrix for average chi2 calculations
    if( !m_avg_data_unfolded[num_iter] ){
       m_avg_data_unfolded[num_iter] = new MnvH( h_data_unfolded->GetCVHistoWithStatError() );
       m_avg_data_unfolded[num_iter]->SetName( Form("Avg_Iter%d", num_iter) );
    }
    else
    {
      MnvH1D* h_tmp_unfold = new MnvH1D(h_data_unfolded->GetCVHistoWithStatError());
      m_avg_data_unfolded[num_iter]->Add( h_tmp_unfold ); 
    }

    int nCovCols = unfoldingCovMatrixOrig.GetNcols();
    int nCovRows = unfoldingCovMatrixOrig.GetNrows();
    if( !m_avg_unfoldingCovMatrices[num_iter] ) m_avg_unfoldingCovMatrices[num_iter] = new TMatrixD( nCovRows, nCovCols );
    TMatrixD* avg_unfoldingCovMatrix = m_avg_unfoldingCovMatrices[num_iter]; 

    for( int row = 0; row < nCovRows; ++row)
    {
      for( int col = 0; col < nCovCols; ++col)
      {
        avg_unfoldingCovMatrix[0][row][col] += (double)unfoldingCovMatrixOrig(row,col);
      }
    }
    m_avg_unfoldingCovMatrices[num_iter] = avg_unfoldingCovMatrix; 

    return true;
  }
  template <class MnvH> bool TransWarpExtraction<MnvH>::UnfoldData( MnvH2D* &h_data_unfolded, MnvH2D* stat_varied, MnvH2D* h_reco, MnvH2D* h_truth, MnvH2D* h_migration, int num_iter ) 
  {
    bool data_unfolded = false;
    data_unfolded = unfold.UnfoldHisto2D(h_data_unfolded, h_migration, h_reco, h_truth, stat_varied, num_iter, false, false );
    
    if( !data_unfolded )
    {
      std::cout << "Unfolding failed for either data or MC. Please check " << std::endl;
      return false;  
    }

    TMatrixD unfoldingCovMatrixOrig = UnfoldDummy( stat_varied, h_data_unfolded, h_reco, h_truth, h_migration, num_iter ); 
    const int highBinX = h_data_unfolded->GetNbinsX() + 1;
    const int highBinY = h_data_unfolded->GetNbinsY() + 1;
    const int highBin  = h_data_unfolded->GetBin( highBinX, highBinY );
    const int lowBin = 0;
    for( int i = lowBin; i <= highBin; ++i )
    {
      for(uint bin=0; bin<m_exclude_chi2_bins.size(); ++bin)
      {
        unfoldingCovMatrixOrig[i][m_exclude_chi2_bins[bin]] = 0;
        unfoldingCovMatrixOrig[m_exclude_chi2_bins[bin]][i] = 0;
      }
    }

    if( m_corr_factor > 0 ) {
      double uncmod  = (1. + 1./m_corr_factor); // To apply to the covariance matrix
      double sqrtmod = sqrt(uncmod);            // To apply to the diagonal error

      // Scale covariance matrix
      unfoldingCovMatrixOrig *= uncmod;

      // Scale diagonal errors
      for( int i = 0; i < h_data_unfolded->GetNbinsX() + 2; i++ ){
        for( int j = 0; j < h_data_unfolded->GetNbinsY() + 2; j++ ){
          h_data_unfolded->SetBinError( i, j, h_data_unfolded->GetBinError(i, j)*sqrtmod );
        }
      }
    }

    m_unfoldingCovMatrices[num_iter].push_back(unfoldingCovMatrixOrig);
    h_data_unfolded->PushCovMatrix("unfoldingCov",unfoldingCovMatrixOrig);

    //Getting hists and covariance matrix for average chi2 calculations
    if( !m_avg_data_unfolded[num_iter] ){
       m_avg_data_unfolded[num_iter] = new MnvH( h_data_unfolded->GetCVHistoWithStatError() );
       m_avg_data_unfolded[num_iter]->SetName( Form("Avg_Iter%d", num_iter) );
    }
    else
    {
      MnvH2D* h_tmp_unfold = new MnvH2D(h_data_unfolded->GetCVHistoWithStatError());
      m_avg_data_unfolded[num_iter]->Add( h_tmp_unfold ); 
    }

    int nCovCols = unfoldingCovMatrixOrig.GetNcols();
    int nCovRows = unfoldingCovMatrixOrig.GetNrows();
    if( !m_avg_unfoldingCovMatrices[num_iter] ) m_avg_unfoldingCovMatrices[num_iter] = new TMatrixD( nCovRows, nCovCols );
    TMatrixD* avg_unfoldingCovMatrix = m_avg_unfoldingCovMatrices[num_iter]; 

    for( int row = 0; row < nCovRows; ++row)
    {
      for( int col = 0; col < nCovCols; ++col)
      {
        avg_unfoldingCovMatrix[0][row][col] += (double)unfoldingCovMatrixOrig(row,col);
      }
    }
    m_avg_unfoldingCovMatrices[num_iter] = avg_unfoldingCovMatrix; 

    return true;
  }
  //All 3d methods are suspect.  Putting this here, but THINK CAREFULLY before using this
  template <class MnvH> bool TransWarpExtraction<MnvH>::UnfoldData( MnvH3D* &h_data_unfolded, MnvH3D* stat_varied, MnvH3D* h_reco, MnvH3D* h_truth, MnvH2D* h_migration, int num_iter ) 
  {
    bool data_unfolded = false;
    data_unfolded = unfold.UnfoldHisto3D(h_data_unfolded, h_migration, h_reco, h_truth, stat_varied, num_iter, false, false);
    
    if( !data_unfolded )
    {
      std::cout << "Unfolding failed for either data or MC. Please check " << std::endl;
      return false;  
    }

    //MnvH3D doesn't have a method to push covariance matrix.  Removing this for now
    //TMatrixD unfoldingCovMatrixOrig = UnfoldDummy( stat_varied, h_data_unfolded, h_reco, h_truth, h_migration, num_iter ); 

    //h_data_unfolded->PushCovMatrix("unfoldingCov",unfoldingCovMatrixOrig);
    
    ////Getting hists and covariance matrix for average chi2 calculations
    //if( !m_avg_data_unfolded[num_iter] ) m_avg_data_unfolded[num_iter] = new MnvH( h_data_unfolded->GetCVHistoWithStatError() );
    //else
    //{
    //  TH3D* h_tmp_unfold = new TH2D(h_data_unfolded->GetCVHistoWithStatError());
    //  m_avg_data_unfolded[num_iter]->Add( h_tmp_unfold ); 
    //}

    //int nCovCols = unfoldingCovMatrixOrig.GetNcols();
    //int nCovRows = unfoldingCovMatrixOrig.GetNrows();
    //if( !m_avg_unfoldingCovMatrices[num_iter] ) m_avg_unfoldingCovMatrices[num_iter] = new TMatrixD( nCovRows, nCovCols );
    //TMatrixD* avg_unfoldingCovMatrix = m_avg_unfoldingCovMatrices[num_iter]; 

    //for( int row = 0; row < nCovRows; ++row)
    //{
    //  for( int col = 0; col < nCovCols; ++col)
    //  {
    //    avg_unfoldingCovMatrix[0][row][col] += unfoldingCovMatrixOrig(row,col)/m_nStatUniverses;
    //  }
    //}
    //m_avg_unfoldingCovMatrices[num_iter] = avg_unfoldingCovMatrix; 
    
    return true;
  }
  template <class MnvH> TMatrixD TransWarpExtraction<MnvH>::UnfoldDummy( MnvH1D* stat_varied, MnvH1D* h_data_unfolded, MnvH1D* h_reco, MnvH1D* h_truth, MnvH2D* h_migration, int num_iter ) 
  {
    std::cout << "Getting the covariance of the unfolding" << std::endl;
    TMatrixD unfoldingCovMatrixOrig;
    int correctNbins;
    int matrixRows;  

    TH1D* hUnfoldedDummy  = new TH1D(h_data_unfolded->GetCVHistoWithStatError());
    TH1D* hRecoDummy      = new TH1D(h_reco->GetCVHistoWithStatError());
    TH1D* hTruthDummy     = new TH1D(h_truth->GetCVHistoWithStatError());
    TH1D* hBGSubDataDummy = new TH1D(stat_varied->GetCVHistoWithStatError());
    TH2D* hMigrationDummy = new TH2D(h_migration->GetCVHistoWithStatError());
    unfold.UnfoldHisto(hUnfoldedDummy, unfoldingCovMatrixOrig, hMigrationDummy, hRecoDummy, hTruthDummy, hBGSubDataDummy,RooUnfold::kBayes, num_iter);//Stupid RooUnfold.  This is dummy, we don't need iterations

    correctNbins=hUnfoldedDummy->fN;
    matrixRows=unfoldingCovMatrixOrig.GetNrows();

    if(correctNbins!=matrixRows){
      std::cout << "****************************************************************************" << std::endl;
      std::cout << "*  Fixing unfolding matrix size because of RooUnfold bug. From " << matrixRows << " to " << correctNbins << std::endl;
      std::cout << "****************************************************************************" << std::endl;
      // It looks like this, since the extra last two bins don't have any content
      unfoldingCovMatrixOrig.ResizeTo(correctNbins, correctNbins);
    }

    for(int i=0; i<unfoldingCovMatrixOrig.GetNrows(); ++i) unfoldingCovMatrixOrig(i,i)=0;
    delete hUnfoldedDummy;
    delete hMigrationDummy;
    delete hRecoDummy;
    delete hTruthDummy;
    delete hBGSubDataDummy;

    return unfoldingCovMatrixOrig;
  }
  template <class MnvH> TMatrixD TransWarpExtraction<MnvH>::UnfoldDummy( MnvH2D* stat_varied, MnvH2D* h_data_unfolded, MnvH2D* h_reco, MnvH2D* h_truth, MnvH2D* h_migration, int num_iter ) 
  {
    std::cout << "Getting the covariance of the unfolding" << std::endl;
    TMatrixD unfoldingCovMatrixOrig;
    int correctNbins;
    int matrixRows;  

    TH2D* hUnfoldedDummy  = new TH2D(h_data_unfolded->GetCVHistoWithStatError());
    TH2D* hRecoDummy      = new TH2D(h_reco->GetCVHistoWithStatError());
    TH2D* hTruthDummy     = new TH2D(h_truth->GetCVHistoWithStatError());
    TH2D* hBGSubDataDummy = new TH2D(stat_varied->GetCVHistoWithStatError());
    TH2D* hMigrationDummy = new TH2D(h_migration->GetCVHistoWithStatError());
    unfold.UnfoldHisto2D(hUnfoldedDummy, unfoldingCovMatrixOrig, hMigrationDummy, hRecoDummy, hTruthDummy, hBGSubDataDummy, num_iter);

    correctNbins=hUnfoldedDummy->fN;
    matrixRows=unfoldingCovMatrixOrig.GetNrows();

    if(correctNbins!=matrixRows){
      std::cout << "****************************************************************************" << std::endl;
      std::cout << "*  Fixing unfolding matrix size because of RooUnfold bug. From " << matrixRows << " to " << correctNbins << std::endl;
      std::cout << "****************************************************************************" << std::endl;
      // It looks like this, since the extra last two bins don't have any content
      unfoldingCovMatrixOrig.ResizeTo(correctNbins, correctNbins);
    }

    for(int i=0; i<unfoldingCovMatrixOrig.GetNrows(); ++i) unfoldingCovMatrixOrig(i,i)=0;
    delete hUnfoldedDummy;
    delete hMigrationDummy;
    delete hRecoDummy;
    delete hTruthDummy;
    delete hBGSubDataDummy;

    return unfoldingCovMatrixOrig;
  }
  template <class MnvH> TMatrixD TransWarpExtraction<MnvH>::UnfoldDummy( MnvH3D* stat_varied, MnvH3D* h_data_unfolded, MnvH3D* h_reco, MnvH3D* h_truth, MnvH2D* h_migration, int num_iter ) 
  {
    std::cout << "Getting the covariance of the unfolding" << std::endl;
    TMatrixD unfoldingCovMatrixOrig;
    int correctNbins;
    int matrixRows;  

    TH3D* hUnfoldedDummy  = new TH3D(h_data_unfolded->GetCVHistoWithStatError());
    TH3D* hRecoDummy      = new TH3D(h_reco->GetCVHistoWithStatError());
    TH3D* hTruthDummy     = new TH3D(h_truth->GetCVHistoWithStatError());
    TH3D* hBGSubDataDummy = new TH3D(stat_varied->GetCVHistoWithStatError());
    TH2D* hMigrationDummy = new TH2D(h_migration->GetCVHistoWithStatError());
    unfold.UnfoldHisto3D(hUnfoldedDummy, unfoldingCovMatrixOrig, hMigrationDummy, hRecoDummy, hTruthDummy, hBGSubDataDummy, num_iter);

    correctNbins=hUnfoldedDummy->fN;
    matrixRows=unfoldingCovMatrixOrig.GetNrows();

    if(correctNbins!=matrixRows){
      std::cout << "****************************************************************************" << std::endl;
      std::cout << "*  Fixing unfolding matrix size because of RooUnfold bug. From " << matrixRows << " to " << correctNbins << std::endl;
      std::cout << "****************************************************************************" << std::endl;
      // It looks like this, since the extra last two bins don't have any content
      unfoldingCovMatrixOrig.ResizeTo(correctNbins, correctNbins);
    }

    for(int i=0; i<unfoldingCovMatrixOrig.GetNrows(); ++i) unfoldingCovMatrixOrig(i,i)=0;
    delete hUnfoldedDummy;
    delete hMigrationDummy;
    delete hRecoDummy;
    delete hTruthDummy;
    delete hBGSubDataDummy;

    return unfoldingCovMatrixOrig;
  }
  template <class MnvH> void TransWarpExtraction<MnvH>::CalcChi2( )
  {
    for(std::vector<int>::iterator num_iter = m_iterations.begin(); num_iter != m_iterations.end(); ++num_iter)
    {
      MnvH *input_truth = (MnvH*)m_data_truth->Clone("input_truth");
      MnvH *unfoldtruth = (MnvH*)m_truth->Clone("unfoldtruth");
      
      input_truth->ClearAllErrorBands();//Do this so we don't pick up all the systematics on the result!
      //Need to set inputreco stat error to be "data like" this will tell us bias scale!

      for(int bin=0;bin<input_truth->fN;bin++)
      {
        double bincont = input_truth->GetBinContent(bin);
        input_truth->SetBinError(bin,sqrt(bincont));
      }

      for(uint bin=0; bin<m_exclude_chi2_bins.size(); ++bin)
      {
        input_truth->SetBinContent(m_exclude_chi2_bins[bin],0);  
        unfoldtruth->SetBinContent(m_exclude_chi2_bins[bin],0);  
        input_truth->SetBinError(m_exclude_chi2_bins[bin],0);  
        unfoldtruth->SetBinError(m_exclude_chi2_bins[bin],0);  
      }

      for(uint su = 0; su<m_nStatUniverses; ++su) 
      {
        if(!m_outputresults[*num_iter][su])
        {
          std::cout<<"Trying to call CalcChi2 without an unfolded data histogram!"<<std::endl;
          std::cout<<"Try to call UnfoldStatUniverses first"<<std::endl;
          return;
        }
        for(uint bin=0; bin<m_exclude_chi2_bins.size(); ++bin)
        {
          m_outputresults[*num_iter][su]->SetBinContent(m_exclude_chi2_bins[bin],0);  
          m_outputresults[*num_iter][su]->SetBinError(m_exclude_chi2_bins[bin],0);  
        }
        FillChi2Dists( su, m_outputresults[*num_iter][su], input_truth, unfoldtruth, *num_iter );
      }

      for(uint bin=0; bin<m_exclude_chi2_bins.size(); ++bin)
      {
        m_avg_data_unfolded[*num_iter]->SetBinContent(m_exclude_chi2_bins[bin],0);  
        m_avg_data_unfolded[*num_iter]->SetBinError(m_exclude_chi2_bins[bin],0);  
      }
      m_avg_data_unfolded[*num_iter] = ScaleAverageUnfold(m_avg_data_unfolded[*num_iter],*num_iter);

      /*This is what was used to get the chi^2 of the "average" unfolded... not particularly interesting, rife with problems */
      //int ndf = 0;
      //TMatrixD *valmat = new TMatrixD(500,500);
      //double avg_chi2_md_md   = plotter->Chi2DataMC(m_avg_data_unfolded[*num_iter],m_avg_data_unfolded[*num_iter],ndf,1,true,false,false); //(1)
      //double avg_chi2_md_td   = plotter->Chi2DataMC(m_avg_data_unfolded[*num_iter],input_truth,ndf,1,true,false,false,valmat); //(3)
      //double avg_chi2_md_tmc  = plotter->Chi2DataMC(m_avg_data_unfolded[*num_iter],unfoldtruth,ndf,1,true,false,false); // (4)
      //double avg_chi2_td_tmc  = plotter->Chi2DataMC(input_truth,unfoldtruth,ndf,1,true,false); // (6)

      //m_avg_valmat_vec[*num_iter]=(TMatrixD*)valmat->Clone();
      
      //m_avg_chi2_md_md_iter_chi2 ->Fill(*num_iter,avg_chi2_md_md );
      //m_avg_chi2_md_td_iter_chi2 ->Fill(*num_iter,avg_chi2_md_td );
      //m_avg_chi2_md_tmc_iter_chi2->Fill(*num_iter,avg_chi2_md_tmc);
      //m_avg_chi2_td_tmc_iter_chi2->Fill(*num_iter,avg_chi2_td_tmc);

      //m_avg_chi2_md_md_iter_chi2 ->SetBinError( m_avg_chi2_md_md_iter_chi2 ->FindFixBin( *num_iter), 0 ); 
      //m_avg_chi2_md_td_iter_chi2 ->SetBinError( m_avg_chi2_md_td_iter_chi2 ->FindFixBin( *num_iter), 0 );
      //m_avg_chi2_md_tmc_iter_chi2->SetBinError( m_avg_chi2_md_tmc_iter_chi2->FindFixBin( *num_iter), 0 );
      //m_avg_chi2_td_tmc_iter_chi2->SetBinError( m_avg_chi2_td_tmc_iter_chi2->FindFixBin( *num_iter), 0 );

    }
  }
  template <class MnvH> void TransWarpExtraction<MnvH>::MakeBinChi2Dists( )
  {
    //Make sure we have something to work with
    if( m_valmat_vec[*m_iterations.begin()].size()==0 ) return;

    //Get the rows and columns 
    int nCols = m_valmat_vec[*m_iterations.begin()][0][0].GetNcols();
    int nRows = m_valmat_vec[*m_iterations.begin()][0][0].GetNrows();

    //Setting up the binning for the chi^2.  Assuming it's going to be proportionally divided by bin 
    uint nChiBins = (int)m_ChiSquareGuess/m_ChiSquareStep;
    Double_t chi_bins[nChiBins+1];
    for(uint chi = 0; chi <= nChiBins; ++chi) chi_bins[chi]=(double)chi*(m_ChiSquareStep/(double)nRows);

    uint nIterBins = m_iterations_bins.size();
    Double_t iter_bins[nIterBins];
    for(uint iter = 0; iter< nIterBins; ++iter) iter_bins[iter]=m_iterations_bins[iter];

    for( int iRow = 0; iRow < nRows; ++iRow)
    {
      m_bin_chi2_iter[iRow+1]     = new MnvH2D(Form("Bin_%d_Iter_Chi2",iRow+1), Form(";#bf{ Bin %d } Iterations; #chi^{2}",iRow+1),(int)nIterBins-1,iter_bins,(int)nChiBins,chi_bins);        
      m_avg_bin_chi2_iter[iRow+1] = new TProfile(Form("Avg_Bin_%d_Iter_Chi2",iRow+1),Form(";#bf{ Bin %d } Iterations; #chi^{2}",iRow+1),(int)nIterBins-1,iter_bins);        
      for(std::vector<int>::iterator num_iter = m_iterations.begin(); num_iter != m_iterations.end(); ++num_iter)
      {
        for( uint su = 0; su <m_nStatUniverses; ++su )
        { 
          double tmp_chi2 = 0;
          for( int iCol = 0; iCol < nCols; ++iCol) tmp_chi2 += m_valmat_vec[*num_iter][su][0][iRow][iCol];
          m_bin_chi2_iter[iRow+1]->Fill(*num_iter,tmp_chi2,1);
          m_avg_bin_chi2_iter[iRow+1]->Fill(*num_iter,tmp_chi2);
        }
      }
    }
  }
  template <class MnvH> void TransWarpExtraction<MnvH>::MakeMedianChi2Dists( double percentile )
  {
    //Percentile is the percentage around the median of universes we're including.  0.5 mean 25% on either side
    if( !m_chi2_md_td_iter_stat ) return;

    uint nIterBins = m_iterations_bins.size();
    Double_t iter_bins[nIterBins];
    for(uint iter = 0; iter< nIterBins; ++iter) iter_bins[iter]=m_iterations_bins[iter];

    m_chi2_md_td_iter_chi2_percentile = (MnvH2D*)m_chi2_md_td_iter_chi2->Clone(Form("%s_percentile",m_chi2_md_td_iter_chi2->GetName()));
    m_median_chi2_md_td_iter_chi2     = new MnvH1D("h_median_chi2_modelData_trueData_iter_chi2",
                                                   Form(";%s;Median %s",m_chi2_md_td_iter_chi2->GetXaxis()->GetTitle(),m_chi2_md_td_iter_chi2->GetYaxis()->GetTitle()),
                                                   nIterBins-1, iter_bins );

    //Find the beginning, median, and end of the percentile
    uint lower_per = (uint)( m_nStatUniverses*( 0.5*(1 - percentile) ) );
    uint median = (uint)( m_nStatUniverses*0.5 );
    uint upper_per = (uint)( m_nStatUniverses*( 0.5*(1 + percentile) ) );

    std::vector<double> chi2_iter_values;
    for(uint iter = 0; iter <  m_iterations.size(); ++iter)
    {
      chi2_iter_values.clear();
      for(uint su = 1; su <= m_nStatUniverses; ++su) chi2_iter_values.push_back( m_chi2_md_td_iter_stat->GetBinContent(iter+1,su) );
      //Sort the chi2 
      std::sort( chi2_iter_values.begin(), chi2_iter_values.end() );     

      m_median_chi2_md_td_iter_chi2->Fill( m_iterations[iter], chi2_iter_values.at(median) ); 
      for(uint su = lower_per; su <= upper_per; ++su)
      {
        m_chi2_md_td_iter_chi2_percentile->Fill( m_iterations[iter], chi2_iter_values.at(su) ); 
      }
    }      
  }
  template <class MnvH> void TransWarpExtraction<MnvH>::MakeTruncatedChi2Dists( double nStdDev )
  {
    if( !m_chi2_md_td_iter_stat ) return;
    m_chi2_md_td_iter_chi2_trunc = (MnvH2D*)m_chi2_md_td_iter_chi2->Clone(Form("%s_truncated",m_chi2_md_td_iter_chi2->GetName()));
    m_avg_chi2_md_td_iter_chi2_trunc = (TProfile*)m_avg_chi2_md_td_iter_chi2->Clone(Form("%s_truncated",m_avg_chi2_md_td_iter_chi2->GetName()));

    for(uint iter = 0; iter <  m_iterations.size(); ++iter)
    {
      std::vector<double> chi2_iter_values;
      for(uint su = 1; su <= m_nStatUniverses; ++su) chi2_iter_values.push_back( m_chi2_md_td_iter_stat->GetBinContent(iter+1,su) );
      
      double mean = std::accumulate( chi2_iter_values.begin(), chi2_iter_values.end(), 0.0)/chi2_iter_values.size();
      double rms  = std::inner_product( chi2_iter_values.begin(), chi2_iter_values.end(), chi2_iter_values.begin(), 0.0 )/chi2_iter_values.size();
      double std_dev = sqrt( rms - mean*mean );

      for(uint su = 1; su <= m_nStatUniverses; ++su)
      {
        if( fabs(m_chi2_md_td_iter_stat->GetBinContent(iter,su)-mean) > nStdDev*std_dev ) continue;
        m_chi2_md_td_iter_chi2_trunc->Fill( m_iterations[iter], m_chi2_md_td_iter_stat->GetBinContent(iter,su) ); 
        m_avg_chi2_md_td_iter_chi2_trunc->Fill( m_iterations[iter], m_chi2_md_td_iter_stat->GetBinContent(iter,su) ); 
      }
    }      
  }
  template <class MnvH> void TransWarpExtraction<MnvH>::MakeErrorHists( )
  {
    int nUnfoldedBins = m_outputresults[*m_iterations.begin()][0]->fN; 

    uint nIterBins = m_iterations_bins.size();
    Double_t iter_bins[nIterBins];
    for(uint iter = 0; iter< nIterBins; ++iter) iter_bins[iter]=m_iterations_bins[iter];

    int nPullBins = 100;
    Double_t pull_bins[nPullBins+1];
    for(int pull = 0; pull < nPullBins+1; ++pull) pull_bins[pull]=-2.0+0.04*pull; 

    int nErrorBins = (int)1.5*sqrt(m_outputresults[*m_iterations.begin()][0]->GetMaximum()); //Guesstimating.  Assuming the max error found will be about ~ 1.5 max error (sqrt of content)
    Double_t error_bins[nErrorBins+1];
    for(int error = 0; error < nErrorBins+1; ++error) error_bins[error]=error; 

    for( int iBin = 0; iBin < nUnfoldedBins; ++iBin )
    {
      m_bin_data_pulls[iBin] =  new MnvH2D(Form("Bin_%i_Iter_pull", iBin),Form(";Iterations (#bf{Bin} %i);#frac{reco-true}{unc}",iBin),nIterBins-1, iter_bins,nPullBins ,pull_bins); 
      m_bin_data_errors[iBin] = new MnvH2D(Form("Bin_%i_Iter_error",iBin),Form(";Iterations (#bf{Bin} %i);error"                ,iBin),nIterBins-1, iter_bins,nErrorBins,error_bins); 
      m_bin_avg_data_pulls[iBin] =  new TProfile(Form("Avg_Bin_%i_Iter_pull", iBin),Form(";Iterations (#bf{Bin} %i);#frac{reco-true}{unc}",iBin),nIterBins-1, iter_bins,"s"); 
      m_bin_avg_data_errors[iBin] = new TProfile(Form("Avg_Bin_%i_Iter_error",iBin),Form(";Iterations (#bf{Bin} %i);error"                ,iBin),nIterBins-1, iter_bins,"s"); 
    }

    for(std::vector<int>::iterator num_iter = m_iterations.begin(); num_iter != m_iterations.end(); ++num_iter)
    {
      MnvH *input_truth = (MnvH*)m_data_truth->Clone("input_truth2");
      
      input_truth->ClearAllErrorBands();//Do this so we don't pick up all the systematics on the result!

      for(int bin=0;bin<input_truth->fN;bin++){
        double bincont = input_truth->GetBinContent(bin);
        input_truth->SetBinError(bin,sqrt(bincont));
      }

      MnvH* avg_error = (MnvH*)m_outputresults[*num_iter][0]->Clone(Form("avg_Iter_%i_error",*num_iter));
      MnvH* avg_pull  = (MnvH*)m_outputresults[*num_iter][0]->Clone(Form("avg_Iter_%i_pull",*num_iter));
      MnvH* avg_rms   = (MnvH*)m_outputresults[*num_iter][0]->Clone(Form("avg_Iter_%i_pull_rms",*num_iter));
      avg_error->Reset(); 
      avg_pull ->Reset(); 
      avg_rms ->Reset(); 

      avg_error->GetXaxis()->SetTitle(m_outputresults[*num_iter][0]->GetXaxis()->GetTitle());
      avg_pull->GetXaxis()->SetTitle(m_outputresults[*num_iter][0]->GetXaxis()->GetTitle());
      avg_error->GetYaxis()->SetTitle("error");
      avg_pull->GetYaxis()->SetTitle("#frac{reco-true}{unc}");

      for(uint su = 0; su<m_nStatUniverses; ++su) 
      {
        for(int bin = 0; bin<nUnfoldedBins; ++bin)
        {
          avg_error->SetBinContent(bin,avg_error->GetBinContent(bin) + m_outputresults[*num_iter][su]->GetBinError(bin));
          m_bin_data_errors[bin]->Fill(*num_iter,m_outputresults[*num_iter][su]->GetBinError(bin));
          m_bin_avg_data_errors[bin]->Fill(*num_iter,m_outputresults[*num_iter][su]->GetBinError(bin));
        }

        MnvH* tmp_pull  = new MnvH(m_outputresults[*num_iter][su]->GetCVHistoWithStatError());
        tmp_pull->SetName(Form("%s_pull",m_outputresults[*num_iter][su]->GetName()));
        tmp_pull->GetYaxis()->SetTitle("#frac{reco-true}{unc}");
        tmp_pull->Add(input_truth,-1);
        for(int bin = 0; bin<nUnfoldedBins; ++bin){
          if( fabs(m_outputresults[*num_iter][su]->GetBinError(bin)) > 0 ){
            tmp_pull->SetBinContent(bin, tmp_pull->GetBinContent(bin)/m_outputresults[*num_iter][su]->GetBinError(bin));
            tmp_pull->SetBinError(bin, 0);
            m_bin_data_pulls[bin]->Fill(*num_iter,tmp_pull->GetBinContent(bin));
            m_bin_avg_data_pulls[bin]->Fill(*num_iter,tmp_pull->GetBinContent(bin));
            avg_rms ->SetBinContent(bin, avg_rms->GetBinContent(bin)+(tmp_pull->GetBinContent(bin)*tmp_pull->GetBinContent(bin)) );
          }
          else
          {
            tmp_pull->SetBinContent(bin, 0);
            tmp_pull->SetBinError(bin, 0);
            m_bin_data_pulls[bin]->Fill(*num_iter,0);
            m_bin_avg_data_pulls[bin]->Fill(*num_iter,0);
            avg_rms ->SetBinContent(bin, avg_rms->GetBinContent(bin)+0);
          }
        }

        avg_pull->Add(tmp_pull);
        m_data_pull[*num_iter].push_back(tmp_pull);
      }

      avg_error->Scale(1.0/m_nStatUniverses);
      avg_pull->Scale(1.0/m_nStatUniverses);
      avg_rms ->Scale(1.0/m_nStatUniverses);
      
      for( int iBin = 0; iBin <= avg_pull->GetNbinsX()+1; iBin++)
      {
        avg_pull->SetBinError(iBin,sqrt(avg_rms->GetBinContent(iBin)));
        avg_error->SetBinError(iBin,0);
      }

      m_avg_data_error[*num_iter] = avg_error;
      m_avg_data_pull[*num_iter]  = avg_pull;
    }
  }
  template <class MnvH> void TransWarpExtraction<MnvH>::MakeRatioHists( )
  {
    int nUnfoldedBins = m_outputresults[*m_iterations.begin()][0]->fN; 

    uint nIterBins = m_iterations_bins.size();
    Double_t iter_bins[nIterBins];
    for(uint iter = 0; iter< nIterBins; ++iter) iter_bins[iter]=m_iterations_bins[iter];

    int nRatioBins = 200;
    Double_t ratio_bins[nRatioBins+1];
    for(int ratio = 0; ratio < nRatioBins+1; ++ratio) ratio_bins[ratio]=0.01*ratio; 

    for( int iBin = 0; iBin < nUnfoldedBins; ++iBin )
    {
      m_bin_iter_ratio_md_td[iBin]    = new MnvH2D(Form("Bin_%i_Iter_Ratio_modelData_trueData",iBin),   Form(";Iterations (#bf{Bin} %i);Unfolded Data/True Data",iBin) ,nIterBins-1, iter_bins, nRatioBins, ratio_bins);
      m_bin_iter_ratio_md_input[iBin] = new MnvH2D(Form("Bin_%i_Iter_Ratio_modelData_input",iBin),Form(";Iterations (#bf{Bin} %i);Unfolded Data/Input Data",iBin),nIterBins-1, iter_bins, nRatioBins, ratio_bins);
      m_bin_iter_avg_ratio_md_td[iBin]    = new TProfile(Form("Avg_Bin_%i_Iter_Ratio_modelData_trueData",iBin),   Form(";Iterations (#bf{Bin} %i);Unfolded Data/True Data",iBin) ,nIterBins-1, iter_bins);
      m_bin_iter_avg_ratio_md_input[iBin] = new TProfile(Form("Avg_Bin_%i_Iter_Ratio_modelData_input",iBin),Form(";Iterations (#bf{Bin} %i);Unfolded Data/Input Data",iBin),nIterBins-1, iter_bins);
    }

    MnvH *input_truth = (MnvH*)m_data_truth->Clone("input_truth3");
    input_truth->ClearAllErrorBands();//Do this so we don't pick up all the systematics on the result!

    for(int bin=0;bin<input_truth->fN;bin++){
      double bincont = input_truth->GetBinContent(bin);
      input_truth->SetBinError(bin,sqrt(bincont));
    }

    for(std::vector<int>::iterator num_iter = m_iterations.begin(); num_iter != m_iterations.end(); ++num_iter)
    {
      for(uint su = 0; su<m_nStatUniverses; ++su) 
      {
        MnvH* tmp_ratio_md_td    = new MnvH(m_outputresults[*num_iter][su]->GetCVHistoWithStatError());
        MnvH* tmp_ratio_md_input = new MnvH(m_outputresults[*num_iter][su]->GetCVHistoWithStatError());
        tmp_ratio_md_td    ->SetName(Form("Ratio_modelData_trueData_Stat_%i_Iter_%i",su,*num_iter));
        tmp_ratio_md_input ->SetName(Form("Ratio_modelData_input_Stat_%i_Iter_%i",su,*num_iter));
        tmp_ratio_md_td    ->Divide(tmp_ratio_md_td   ,input_truth);
        tmp_ratio_md_input ->Divide(tmp_ratio_md_input,m_inputresults[*num_iter][su]);

        for(int bin = 0; bin<nUnfoldedBins; ++bin)
        {
          tmp_ratio_md_td    ->SetBinError(bin,0);
          tmp_ratio_md_input ->SetBinError(bin,0);
          
          m_bin_iter_ratio_md_td[bin]   ->Fill(*num_iter, tmp_ratio_md_td    ->GetBinContent(bin));
          m_bin_iter_ratio_md_input[bin]->Fill(*num_iter, tmp_ratio_md_input ->GetBinContent(bin));
          m_bin_iter_avg_ratio_md_td[bin]   ->Fill(*num_iter, tmp_ratio_md_td    ->GetBinContent(bin)); 
          m_bin_iter_avg_ratio_md_input[bin]->Fill(*num_iter, tmp_ratio_md_input ->GetBinContent(bin)); 
        }
        m_ratio_md_td[*num_iter].push_back(tmp_ratio_md_td);
        m_ratio_md_input[*num_iter].push_back(tmp_ratio_md_input);
      }
    }
  }
  template <class MnvH> void TransWarpExtraction<MnvH>::FillChi2Dists( int stat_uni, MnvH* h_data_unfolded, MnvH* input_truth, MnvH* unfoldtruth, int num_iter )
  {
    //I expect (1) to be 0
    //I expect (2) and (4) to be identical (maybe there is an error matrix I don't know about???)
    //      * I have verified that the outputs of (2) and (4) are iteration independent and identical
    //I expect (3) to minimize after some number of iterations.
    //do chi2
    int ndf = 0;

    TMatrixD *valmat = new TMatrixD(500,500);
    //Model Data - Model Data
    double chi2_md_md   = plotter->Chi2DataMC(h_data_unfolded,h_data_unfolded,ndf,1,true,false,false); //(1)
    //Model Data - True Data
    double chi2_md_td   = plotter->Chi2DataMC(h_data_unfolded,input_truth,ndf,1,true,false,false,valmat); //(3)
    //Model Data - True MC
    double chi2_md_tmc  = plotter->Chi2DataMC(h_data_unfolded,unfoldtruth,ndf,1,true,false,false); // (4)
    //True Data - True MC
    double chi2_td_tmc  = plotter->Chi2DataMC(input_truth,unfoldtruth,ndf,1,true,false); // (6)

    if(m_chi2_md_md_iter_chi2)
    {
      //printf("A Iter: %3i %3i %5.3e\n",num_iter,stat_uni,chi2_md_td);
      m_chi2_md_md_iter_chi2 ->Fill(num_iter,chi2_md_md , 1);
      m_chi2_md_td_iter_chi2 ->Fill(num_iter,chi2_md_td , 1);
      m_chi2_md_tmc_iter_chi2->Fill(num_iter,chi2_md_tmc, 1);
      m_chi2_td_tmc_iter_chi2->Fill(num_iter,chi2_td_tmc, 1);

      m_chi2_md_md_iter_stat ->Fill(num_iter,stat_uni,chi2_md_md);
      m_chi2_md_td_iter_stat ->Fill(num_iter,stat_uni,chi2_md_td);
      m_chi2_md_tmc_iter_stat->Fill(num_iter,stat_uni,chi2_md_tmc);
      m_chi2_td_tmc_iter_stat->Fill(num_iter,stat_uni,chi2_td_tmc);

      //These if statements are here are to stop stat universes that explode
      if( chi2_md_md > 0) m_avg_chi2_md_md_iter_chi2 ->Fill(num_iter,chi2_md_md );
      if( chi2_md_td > 0) m_avg_chi2_md_td_iter_chi2 ->Fill(num_iter,chi2_md_td );
      if( chi2_md_tmc> 0) m_avg_chi2_md_tmc_iter_chi2->Fill(num_iter,chi2_md_tmc);
      if( chi2_td_tmc> 0) m_avg_chi2_td_tmc_iter_chi2->Fill(num_iter,chi2_td_tmc);
    }

    m_valmat_vec[num_iter].push_back((TMatrixD*)valmat->Clone());
    delete valmat;
  }                  
  template <class MnvH> MnvH* TransWarpExtraction<MnvH>::ScaleAverageUnfold( MnvH1D* avg_unfold, int num_iter )
  {
    //Scale the bin content by number of universes
    MnvH* h_unfolded = (MnvH*)avg_unfold->Clone(Form("%s_Scaled",avg_unfold->GetName()));
    h_unfolded->ClearAllErrorBands();
    for( int bin = 0; bin < h_unfolded->fN; bin++){
      h_unfolded->SetBinContent(bin,h_unfolded->GetBinContent(bin)/m_nStatUniverses);
    }
    //Scaling the covariance matrix
    int nCovCols = m_avg_unfoldingCovMatrices[num_iter]->GetNcols();
    int nCovRows = m_avg_unfoldingCovMatrices[num_iter]->GetNrows();

    for( int row = 0; row < nCovRows; ++row)
    {
      for( int col = 0; col < nCovCols; ++col)
      {
        m_avg_unfoldingCovMatrices[num_iter][0][row][col] = m_avg_unfoldingCovMatrices[num_iter][0][row][col]/m_nStatUniverses;
      }
    }

    //Scale bin error,  I'm assuming i'm doing this right...
    for(uint su = 0; su<m_nStatUniverses; ++su) 
    {
      TH1D tmp_stat_hist = m_outputresults[num_iter][su]->GetCVHistoWithStatError();
      for( int bin = 0; bin < h_unfolded->fN; bin++) 
      {
        if( su == 0 ) h_unfolded->SetBinError(bin,tmp_stat_hist.GetBinError(bin)/m_nStatUniverses);
        else h_unfolded->SetBinError(bin,h_unfolded->GetBinError(bin)+tmp_stat_hist.GetBinError(bin)/m_nStatUniverses);
      }
    }   
    h_unfolded->PushCovMatrix("unfoldingCov",*m_avg_unfoldingCovMatrices[num_iter]);
    return h_unfolded;
  }
  template <class MnvH> MnvH* TransWarpExtraction<MnvH>::ScaleAverageUnfold( MnvH2D* avg_unfold, int num_iter )
  {
    //Scale the bin content by number of universes
    MnvH* h_unfolded = (MnvH*)avg_unfold->Clone(Form("%s_Scaled",avg_unfold->GetName()));
    h_unfolded->ClearAllErrorBands();
    for( int bin = 0; bin < h_unfolded->fN; bin++){
      h_unfolded->SetBinContent(bin,h_unfolded->GetBinContent(bin)/m_nStatUniverses);
    }
    //Scaling the covariance matrix
    int nCovCols = m_avg_unfoldingCovMatrices[num_iter]->GetNcols();
    int nCovRows = m_avg_unfoldingCovMatrices[num_iter]->GetNrows();

    for( int row = 0; row < nCovRows; ++row)
    {
      for( int col = 0; col < nCovCols; ++col)
      {
        m_avg_unfoldingCovMatrices[num_iter][0][row][col] = m_avg_unfoldingCovMatrices[num_iter][0][row][col]/m_nStatUniverses;
      }
    }

    //Scale bin error,  I'm assuming i'm doing this right...
    for(uint su = 0; su<m_nStatUniverses; ++su) 
    {
      TH2D tmp_stat_hist = m_outputresults[num_iter][su]->GetCVHistoWithStatError();
      for( int bin = 0; bin < h_unfolded->fN; bin++) 
      {
        if( su == 0 ) h_unfolded->SetBinError(bin,tmp_stat_hist.GetBinError(bin)/m_nStatUniverses);
        else h_unfolded->SetBinError(bin,h_unfolded->GetBinError(bin)+tmp_stat_hist.GetBinError(bin)/m_nStatUniverses);
      }
    }   
    h_unfolded->PushCovMatrix("unfoldingCov",*m_avg_unfoldingCovMatrices[num_iter]);
    return h_unfolded;
  }
  template <class MnvH> MnvH* TransWarpExtraction<MnvH>::ScaleAverageUnfold( MnvH3D* avg_unfold, int num_iter )
  {
    //Scale the bin content by number of universes
    MnvH* h_unfolded = (MnvH*)avg_unfold->Clone(Form("%s_Scaled",avg_unfold->GetName()));
    h_unfolded->ClearAllErrorBands();
    for( int bin = 0; bin < h_unfolded->fN; bin++){
      h_unfolded->SetBinContent(bin,h_unfolded->GetBinContent(bin)/m_nStatUniverses);
    }
    //Scaling the covariance matrix
    int nCovCols = m_avg_unfoldingCovMatrices[num_iter]->GetNcols();
    int nCovRows = m_avg_unfoldingCovMatrices[num_iter]->GetNrows();

    for( int row = 0; row < nCovRows; ++row)
    {
      for( int col = 0; col < nCovCols; ++col)
      {
        m_avg_unfoldingCovMatrices[num_iter][0][row][col] = m_avg_unfoldingCovMatrices[num_iter][0][row][col]/m_nStatUniverses;
      }
    }

    //Scale bin error,  I'm assuming i'm doing this right...
    for(uint su = 0; su<m_nStatUniverses; ++su) 
    {
      TH3D tmp_stat_hist = m_outputresults[num_iter][su]->GetCVHistoWithStatError();
      for( int bin = 0; bin < h_unfolded->fN; bin++) 
      {
        if( su == 0 ) h_unfolded->SetBinError(bin,tmp_stat_hist.GetBinError(bin)/m_nStatUniverses);
        else h_unfolded->SetBinError(bin,h_unfolded->GetBinError(bin)+tmp_stat_hist.GetBinError(bin)/m_nStatUniverses);
      }
    }   
    h_unfolded->PushCovMatrix("unfoldingCov",*m_avg_unfoldingCovMatrices[num_iter]);
    return h_unfolded;
  }
  template <class MnvH>
  template <class Hist> 
  void TransWarpExtraction<MnvH>::bookHistos( Hist*& ret_hist, std::string filename, std::string hist_name, std::string targetStr, std::string variable, int PDG)
  {
    std::vector <std::string> targetsInTgt;
    int waterFilledEmpty = -1;
    bool istracker = false;
    if (targetStr=="Iron") targetsInTgt = {"2026", "3026", "5026"};
    else if (targetStr=="Carbon") targetsInTgt = {"3006"};
    else if (targetStr=="Lead") targetsInTgt = {"2082", "3082", "4082", "5082"};
    else if (targetStr=="WaterEmpty")
    {
      targetsInTgt = {"6000"};
      waterFilledEmpty = 2;
    }
    else if (targetStr=="WaterFull")
    {
      targetsInTgt = {"6000"};
      waterFilledEmpty = 1;
    }
    else if (targetStr=="Water")
    {
      targetsInTgt = {"6000"};
      waterFilledEmpty = 0;
    }
    else if (targetStr=="Tracker")
    {
      targetsInTgt = {"_petal_0", "_petal_1", "_petal_2", "_petal_3", "_petal_4", "_petal_5", "_petal_6", "_petal_7", "_petal_8", "_petal_9", "_petal_10", "_petal_11"};
      istracker=true;
    }
    else targetsInTgt = {targetStr};
    std::vector<std::string> dirs;
    if (!istracker) dirs = util::findContainingDirectories(filename, "Targets", true, true);
    else dirs = util::findContainingDirectories(filename, "Tracker", true, true);
    if (dirs.size()==0) dirs = {"./"};
    double totalPOT = 0;
    Hist* out_Hist = new Hist();
    for (int c = 0; c < dirs.size(); c++)
    {
      for (auto tgt : targetsInTgt)
      {
        std::string eventLoopBase;
        if (!istracker) eventLoopBase = "/runEventLoopTargetsMC";
        else eventLoopBase = "/runEventLoopTrackerMC";
        if (variable == "pTmu_pZmu" && (hist_name=="migration" || hist_name=="reco" || hist_name=="truth"))
        {
          if (!istracker) eventLoopBase = "/runEventLoopTargets2DMigration";
          else eventLoopBase = "/runEventLoopTracker2DMigration";
          std::cout<<"Getting 2D hists: " << hist_name << " from file " << eventLoopBase << std::endl;
        }
          
        
        std::cout<<"Investigating target "<< tgt << " within directory " << dirs[c] << " which is " << c+1 <<"/"<<dirs.size() <<" playlists identified" <<std::endl;
        std::string migpath = dirs[c] + eventLoopBase + tgt + ".root";

        auto migFile = TFile::Open(migpath.c_str(), "READ");
        if (!migFile)
        {
          std::cerr << "Failed to open MC file " << migpath.c_str() << ".\n";
          Error("CrossSectionPlots",(std::string("Could not get ") + migpath + " or it was empty.").c_str());
          continue;
        }

        std::string mcpath;
        if (!istracker) mcpath = dirs[c] + "/runEventLoopTargetsMC" + tgt + ".root";
        else  mcpath = dirs[c] + "/runEventLoopTrackerMC" + tgt + ".root";
        auto mcFile = TFile::Open(mcpath.c_str(), "READ");
        if (!mcFile)
        {
          std::cerr << "Failed to open MC file " << mcpath.c_str() << ".\n";
          Error("CrossSectionPlots",(std::string("Could not get ") + mcpath + " or it was empty.").c_str());
          migFile->Close();
          continue;
        }
        std::string playlistUsed = util::GetIngredientSmart<TNamed>(*mcFile, "PlaylistUsed")->GetTitle();
        int filledorempty = util::filledOrEmptyMEPlaylist(playlistUsed);

        int nuoranu = util::nuOrAntiNuMode(playlistUsed);
        int nupdg;
        std::string FHCorRHC = "";
        if (nuoranu == 1) nupdg = 14;
        else if (nuoranu == 2) nupdg = -14;
        if (nupdg != PDG){
          migFile->Close();
          mcFile->Close();
          continue;
        }

        double mc_pot = util::GetIngredientSmart<TParameter<double>>(*mcFile, "POTUsed")->GetVal();
        if (waterFilledEmpty>0 && (waterFilledEmpty != filledorempty)) continue;
        util::AddHist(*out_Hist,util::GetIngredientSmart<Hist>(*migFile, (variable+"_"+hist_name)).get(), mc_pot);
        totalPOT+=mc_pot;
        migFile->Close();
        mcFile->Close();
      }
    }
    //out_Hist->Scale(1/totalPOT);
    ret_hist = out_Hist;
  }

  template <class MnvH>
  template <class Hist> 
  void TransWarpExtraction<MnvH>::bookDataHisto( Hist*& ret_hist, std::string filename, std::string hist_name, std::string targetStr, std::string variable, int PDG)
  {
    if (targetStr=="Tracker") return bookHistos(ret_hist, filename,  hist_name, targetStr, variable, PDG);
    TCanvas can("Test");
    std::vector<std::string> dirs;
    dirs = util::findContainingDirectories(filename, "Targets", true, true);
    if (dirs.size()==0) dirs = {"./"};
    int waterFilledEmpty = -1; //Only used with water target

    double waterFilledPOTData = 0;
    double waterEmptyPOTData = 0;
    double waterFilledPOTMC = 0;
    double waterEmptyPOTMC = 0;
    std::cout<<"Investigating directory\n";
    std::vector <std::string> targetsInTgt;
    if (targetStr=="Iron") targetsInTgt = {"2026", "3026", "5026"};
    else if (targetStr=="Carbon") targetsInTgt = {"3006"};
    else if (targetStr=="Lead") targetsInTgt = {"2082", "3082", "4082", "5082"};
    else if (targetStr=="WaterEmpty")
    {
      targetsInTgt = {"6000"};
      waterFilledEmpty = 2;
    }
    else if (targetStr=="WaterFull")
    {
      targetsInTgt = {"6000"};
      waterFilledEmpty = 1;
    }
    else if (targetStr=="Water")
    {
      targetsInTgt = {"6000"};
      waterFilledEmpty = 0; //Do subtractive analysis
      for (int c = 0; c<dirs.size(); c++)
      {
        std::cout<<"Investigating directory " << dirs[c] << " which is " << c+1 <<"/"<<dirs.size() <<" playlists identified" <<std::endl;
        std::string datapath = dirs[c] + "/runEventLoopTargetsData" + targetsInTgt[0] + ".root";
        std::string mcpath = dirs[c] + "/runEventLoopTargetsMC" + targetsInTgt[0] + ".root";

        auto dataFile = TFile::Open(datapath.c_str(), "READ");
        if (!dataFile)
        {
          Error("CrossSectionPlots",(std::string("Could not get ") + datapath + " or it was empty.").c_str());
        }
        double data_pot = util::GetIngredientSmart<TParameter<double>>(*dataFile, "POTUsed")->GetVal();
        auto mcFile = TFile::Open(mcpath.c_str(), "READ");
        if (!mcFile)
        {
          Error("CrossSectionPlots",(std::string("Could not get ") + mcpath + " or it was empty.").c_str());
        }
        double mc_pot = util::GetIngredientSmart<TParameter<double>>(*mcFile, "POTUsed")->GetVal();
        std::string playlistUsed = util::GetIngredientSmart<TNamed>(*mcFile, "PlaylistUsed")->GetTitle();
        int filledorempty = util::filledOrEmptyMEPlaylist(playlistUsed);

        int nuoranu = util::nuOrAntiNuMode(playlistUsed);
        int nupdg;
        if (nuoranu==1) nupdg = 14;
        else if (nuoranu==2) nupdg = -14;
        if (PDG!=nupdg)
        {
          std::cout<<"Skipping this set of files because this playlist pdg is " << nupdg << " but you specified " << PDG << std::endl;
          continue;
        }

        if (filledorempty == 2) //Empty
        {
          waterEmptyPOTData += data_pot;
          waterEmptyPOTMC += mc_pot;
        }
        if (filledorempty == 1) //Filled
        {
          waterFilledPOTData += data_pot;
          waterFilledPOTMC += mc_pot;
        }
      }
      std::cout<<"Total MC POT Filled: "<< waterFilledPOTMC << std::endl;
      std::cout<<"Total MC POT Empty: "<< waterEmptyPOTMC << std::endl;
      std::cout<<"Total Data POT Filled: "<< waterFilledPOTData << std::endl;
      std::cout<<"Total Data POT Empty: "<< waterEmptyPOTData << std::endl;
    }
    else targetsInTgt = {targetStr};
    std::cout<<"Investigating directory1\n";
    //Flux parameters
    int n_flux_universes = 100; // Is this right
    const bool use_nue_constraint = true;
    const std::string project_dir = "targets_2345_jointNueIMD";
    double min_energy = 0;
    double max_energy = 120;

    double mcPOT = 0;
    double dataPOT = 0;
    std::cout<<"Investigating directory2\n";

    double totalMCPOTFull = 0; //Only needed for water subtractive analysis
    double totalDataPOTFull = 0;
    double totalMCPOTEmpty = 0; //Only needed for water subtractive analysis
    double totalDataPOTEMpt = 0;
    double culmulativeMCPOT = 0; //Only needed for water subtractive analysis
    double culmulativeDataPOT = 0;

    PlotUtils::MnvH1D* folded  = new PlotUtils::MnvH1D();

    PlotUtils::MnvH1D* USSidebandSignal  = new PlotUtils::MnvH1D();
    PlotUtils::MnvH1D* USSidebandDS  = new PlotUtils::MnvH1D();
    PlotUtils::MnvH1D* USSidebandUS  = new PlotUtils::MnvH1D();
    PlotUtils::MnvH1D* USSidebandOther  = new PlotUtils::MnvH1D();
    PlotUtils::MnvH1D* DSSidebandSignal  = new PlotUtils::MnvH1D();
    PlotUtils::MnvH1D* DSSidebandDS  = new PlotUtils::MnvH1D();
    PlotUtils::MnvH1D* DSSidebandUS  = new PlotUtils::MnvH1D();
    PlotUtils::MnvH1D* DSSidebandOther  = new PlotUtils::MnvH1D();
    PlotUtils::MnvH1D* DataUSSideband  = new PlotUtils::MnvH1D();
    PlotUtils::MnvH1D* DataDSSideband  = new PlotUtils::MnvH1D();
    PlotUtils::MnvH1D* DataSignal  = new PlotUtils::MnvH1D();

    PlotUtils::MnvH1D* BackgroundWrongSign  = new PlotUtils::MnvH1D();
    PlotUtils::MnvH1D* BackgroundNC  = new PlotUtils::MnvH1D();
    PlotUtils::MnvH1D* BackgroundWaterTank  = new PlotUtils::MnvH1D();
    PlotUtils::MnvH1D* BackgroundDSPlastic  = new PlotUtils::MnvH1D();
    PlotUtils::MnvH1D* BackgroundUSPlastic  = new PlotUtils::MnvH1D();
    PlotUtils::MnvH1D* BackgroundOtherTarget  = new PlotUtils::MnvH1D();
    PlotUtils::MnvH1D* BackgroundVtxElsewhere  = new PlotUtils::MnvH1D();
    PlotUtils::MnvH1D* BackgroundOther  = new PlotUtils::MnvH1D();

    int sidebandSubtraction = 1;
    
    std::cout<<"dirs " << dirs.size()<<"\n";
    for (int c = 0; c<dirs.size(); c++)
    {
      std::cout<<"Investigating directory " << dirs[c] << " which is " << c+1 <<"/"<<dirs.size() <<" playlists identified" <<std::endl;
      for (int d = 0; d<targetsInTgt.size(); d++)
      {
        std::cout<<"Investigating target " << targetsInTgt[d] << " which is " << d+1 <<"/"<<targetsInTgt.size() <<" targets selected" <<std::endl;
        std::string datapath = dirs[c] + "/runEventLoopTargetsData" + targetsInTgt[d] + ".root";
        std::string mcpath = dirs[c] + "/runEventLoopTargetsMC" + targetsInTgt[d] + ".root";

        auto dataFile = TFile::Open(datapath.c_str(), "READ");
        if (!dataFile)
        {
          Error("CrossSectionPlots",(std::string("Could not get ") + datapath + " or it was empty.").c_str());

        }

        auto mcFile = TFile::Open(mcpath.c_str(), "READ");
        if (!mcFile)
        {
          Error("CrossSectionPlots",(std::string("Could not get ") + mcpath + " or it was empty.").c_str());
        }
        
        std::string playlistUsed = util::GetIngredientSmart<TNamed>(*mcFile, "PlaylistUsed")->GetTitle();
        int filledorempty = util::filledOrEmptyMEPlaylist(playlistUsed);
        std::cout<<"playlistUsed: " << playlistUsed << " filledorempty " << filledorempty <<std::endl;

        int nuoranu = util::nuOrAntiNuMode(playlistUsed);
        int nupdg;
        if (nuoranu==1) nupdg = 14;
        else if (nuoranu==2) nupdg = -14;

        if (waterFilledEmpty>0 && (waterFilledEmpty != filledorempty)) continue;
        if (PDG!=nupdg)
        {
          std::cout<<"Skipping this set of files because this playlist pdg is " << nupdg << " but you specified " << PDG << std::endl;
          continue;
        }

        double mcscale = 1; //Only needed for water subtractive analysis
        double datascale = 1;

          //We scale to a standard POT, perform the necessary addition of subtraction for each playlist, keeping track of the culmulative POT (Standard x number of playlists), then scale back down to the actual total POT
        if (waterFilledEmpty == 0 && filledorempty==2)
        {
          mcscale = -(waterFilledPOTMC/waterEmptyPOTMC);
          datascale = -(waterFilledPOTData/waterEmptyPOTData);
          std::cout<<"datascale: " << datascale << std::endl;
          std::cout<<"mcscale: " << mcscale << std::endl;
        }

        std::cout<<"1: " << std::endl;
        //util::AddHist(*flux,util::GetIngredientSmart<PlotUtils::MnvH1D>(*mcFile, (std::string("reweightedflux_integrated")), variable).get(), mcscale);
        util::AddHist(*folded,util::GetIngredientSmart<PlotUtils::MnvH1D>(*mcFile, (std::string("data")), variable).get(), mcscale);
        std::cout<<"2: " << std::endl;
        //If doing scalar sideband tuning
        if (sidebandSubtraction == 1)
        {
          std::cout<<"3: " << std::endl;
          util::AddHist(*USSidebandSignal,util::GetIngredientSmart<PlotUtils::MnvH1D>(*mcFile, (std::string("segment_US_sideband_Signal"))).get(), mcscale);
          util::AddHist(*USSidebandDS,util::GetIngredientSmart<PlotUtils::MnvH1D>(*mcFile, (std::string("segment_US_sideband_DS"))).get(), mcscale);
          util::AddHist(*USSidebandUS,util::GetIngredientSmart<PlotUtils::MnvH1D>(*mcFile, (std::string("segment_US_sideband_US"))).get(), mcscale);
          util::AddHist(*USSidebandOther,util::GetIngredientSmart<PlotUtils::MnvH1D>(*mcFile, (std::string("segment_US_sideband_Other"))).get(), mcscale);
          util::AddHist(*DSSidebandSignal,util::GetIngredientSmart<PlotUtils::MnvH1D>(*mcFile, (std::string("segment_DS_sideband_Signal"))).get(), mcscale);
          std::cout<<"4: " << std::endl;
          util::AddHist(*DSSidebandDS,util::GetIngredientSmart<PlotUtils::MnvH1D>(*mcFile, (std::string("segment_DS_sideband_DS"))).get(), mcscale);
          util::AddHist(*DSSidebandUS,util::GetIngredientSmart<PlotUtils::MnvH1D>(*mcFile, (std::string("segment_DS_sideband_US"))).get(), mcscale);
          util::AddHist(*DSSidebandOther,util::GetIngredientSmart<PlotUtils::MnvH1D>(*mcFile, (std::string("segment_DS_sideband_Other"))).get(), mcscale);
          util::AddHist(*DataUSSideband,util::GetIngredientSmart<PlotUtils::MnvH1D>(*dataFile, (std::string("segment_US_Sideband"))).get(), datascale);
          util::AddHist(*DataDSSideband,util::GetIngredientSmart<PlotUtils::MnvH1D>(*dataFile, (std::string("segment_DS_Sideband"))).get(), datascale);
          util::AddHist(*DataSignal,util::GetIngredientSmart<PlotUtils::MnvH1D>(*dataFile, (std::string("segment_data"))).get(), datascale);
        }

        std::cout<<"5: " << std::endl;
        util::AddHist(*BackgroundWrongSign,util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, (std::string("background_Wrong_Sign_Bkg")), variable), mcscale);
        util::AddHist(*BackgroundNC,util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, (std::string("background_NC_Bkg")), variable), mcscale);
        util::AddHist(*BackgroundWaterTank,util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, (std::string("background_Water_Tank_Bkg")), variable), mcscale);
        std::cout<<"6: " << std::endl;
        util::AddHist(*BackgroundDSPlastic,util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, (std::string("background_Downstream_Plastic_Bkg")), variable), mcscale);
        util::AddHist(*BackgroundUSPlastic,util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, (std::string("background_Upstream_Plastic_Bkg")), variable), mcscale);
        util::AddHist(*BackgroundOtherTarget,util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, (std::string("background_True_In_Other_Target_Bkg")), variable), mcscale);
        util::AddHist(*BackgroundVtxElsewhere,util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, (std::string("background_True_Vtx_Elsewhere_Bkg")), variable), mcscale);
        util::AddHist(*BackgroundOther,util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, (std::string("background_Other")), variable), mcscale);
        std::cout<<"6: " << std::endl;
        if (d==0) //Only get POT once per playlist. E.g we don't wanna double count the POT for the same runs just because we're looking at target 2026 and 3026
        {
          mcPOT += util::GetIngredientSmart<TParameter<double>>(*mcFile, "POTUsed")->GetVal();
          std::cout<<"mcPOT: " << mcPOT << std::endl;
          dataPOT+= util::GetIngredientSmart<TParameter<double>>(*dataFile, "POTUsed")->GetVal();
          std::cout<<"dataPOT: " << dataPOT << std::endl;
        }

        dataFile->Close();
        mcFile->Close();
      }
    }

    std::cout<<"mcPOT: " << mcPOT << std::endl;
    std::cout<<"dataPOT: " << dataPOT << std::endl;

    folded->Draw("HIST");
    can.Print((variable+hist_name+"folded.png").c_str());


    //************************************************
    //Sideband chisq minimisation - single scale factor across a given sideband region
    //************************************************

    bool isNum = true;
    for (char c : targetStr) {
      if (!std::isdigit(c)) {
        isNum = false;
        break;
      }
    }

    if (isNum)
    {
      int tgtCode = std::stoi(targetStr);
      if (tgtCode>=7 && tgtCode <= 19) sidebandSubtraction = 0; // No plastic sideband tuning possible for the pseudotargets
    } 
    double USScaleFactor = 1;
    double DSScaleFactor = 1;
    if (sidebandSubtraction==1)
    {
      PlotUtils::MnvH1D* combinedUS = USSidebandSignal->Clone("CombinedUSMC");
      combinedUS->Add(USSidebandUS);
      combinedUS->Add(USSidebandDS);
      combinedUS->Add(USSidebandOther);
      //GetCVHistoWithStatError();
      combinedUS->Scale(dataPOT/mcPOT);

      PlotUtils::MnvH1D* combinedDS = DSSidebandSignal->Clone("CombinedUSMC");
      combinedDS->Add(DSSidebandUS);
      combinedDS->Add(DSSidebandDS);
      combinedDS->Add(DSSidebandOther);
      combinedDS->Scale(dataPOT/mcPOT);

      //Plot(*combinedUS, "combinedUS", prefix, tgtname);
      //Plot(*DataUSSideband, "DataUSSideband", prefix, tgtname);
      double minValUS, minValDS;
      {          //US Sideband
        Chi2Calculator<MnvH1D> chi2Calc(DataUSSideband, combinedUS);
        ROOT::Math::Minimizer* minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2");
        minimizer->SetMaxFunctionCalls(1000000);
        minimizer->SetMaxIterations(100000);
        minimizer->SetTolerance(0.01);
        minimizer->SetPrintLevel(1);
        ROOT::Math::Functor f(chi2Calc, 1);
        minimizer->SetFunction(f);
        minimizer->SetVariable(0, "scale", 1, 0.001);
        minimizer->Minimize();
        const double *bestfit = minimizer->X();
        double bestChiSq = minimizer->MinValue();
        USScaleFactor = bestfit[0];
        std::cout<<"Minised chisq for USSideband. Scale factor: " << USScaleFactor << " best chisq: " << bestChiSq << std::endl;
      }

      {          //DS Sideband
        Chi2Calculator<MnvH1D> chi2Calc(DataDSSideband, combinedDS);
        ROOT::Math::Minimizer* minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2");
        minimizer->SetMaxFunctionCalls(1000000);
        minimizer->SetMaxIterations(100000);
        minimizer->SetTolerance(0.01);
        minimizer->SetPrintLevel(1);
        ROOT::Math::Functor f(chi2Calc, 1);
        minimizer->SetFunction(f);
        minimizer->SetVariable(0, "scale", 1, 0.001);
        minimizer->Minimize();
        const double *bestfit = minimizer->X();
        double bestChiSq = minimizer->MinValue();
        DSScaleFactor = bestfit[0];
        std::cout<<"Minised chisq for DSSideband. Scale factor: " << DSScaleFactor << " best chisq: " << bestChiSq << std::endl;
      }
      BackgroundUSPlastic->Scale(USScaleFactor);
      BackgroundDSPlastic->Scale(DSScaleFactor);
    }
    //************************************************
    //End - Sideband chisq minimisation - single scale factor across a given sideband region
    //************************************************

    // There are no error bands in the data, but I need somewhere to put error bands on the results I derive from it.
    //folded->AddMissingErrorBandsAndFillWithCV(*migration);

    std::vector<PlotUtils::MnvH1D*> sidebandScaled = {BackgroundUSPlastic, BackgroundDSPlastic};
    std::vector<PlotUtils::MnvH1D*> backgrounds = { BackgroundWrongSign, BackgroundNC, BackgroundOther, BackgroundWaterTank, BackgroundOtherTarget, BackgroundVtxElsewhere};

    //sum up backgrounds without sideband scaling
    auto bkgSubtracted = std::accumulate(backgrounds.begin(), backgrounds.end(), folded->Clone(),
                                          [mcPOT, dataPOT](auto sum, const auto hist)
                                          {
                                            std::cout << "Subtracting " << hist->GetName() << " scaled by " << -dataPOT / mcPOT << " from " << sum->GetName() << "\n";
                                            sum->Add(hist, -1);
                                            return sum;
                                          });

    //Both background subtracted and sideband scaled
    auto bkgScaledSubtracted = std::accumulate(sidebandScaled.begin(), sidebandScaled.end(), bkgSubtracted->Clone(),
                                          [mcPOT, dataPOT](auto sum, const auto hist)
                                          {
                                            std::cout << "Subtracting " << hist->GetName() << " scaled by " << -dataPOT / mcPOT << " from " << sum->GetName() << "\n";
                                            sum->Add(hist, -1);
                                            return sum;
                                          });

    //bkgScaledSubtracted->Scale(1/mcPOT);
    bkgScaledSubtracted->Draw("HIST");
    can.Print((variable+hist_name+"bkgScaledSubtracted.png").c_str());
    ret_hist = bkgScaledSubtracted;
  }

  /* template <class MnvH> void TransWarpExtraction<MnvH>::bookHistos( MnvH2D*& ret_hist, std::string filename, std::string hist_name, std::string targetStr)
  {
    TFile *f = TFile::Open( filename.c_str(), "READ" );
    if (f->IsZombie() || f->GetListOfKeys()->IsEmpty()){
      Error("CrossSectionPlots","Could not get histogram ROOT file or it was empty.");
      return;
    }
    ret_hist = (MnvH2D*)f->Get(hist_name.c_str());

    if ( !ret_hist ){
      Error("CrossSectionPlots",TString::Format("No histogram found with name %s",hist_name.c_str()));
      return;
    }
  }
  template <class MnvH> void TransWarpExtraction<MnvH>::bookHistos( MnvH3D*& ret_hist, std::string filename, std::string hist_name, std::string targetStr)
  {
    TFile *f = TFile::Open( filename.c_str(), "READ" );
    if (f->IsZombie() || f->GetListOfKeys()->IsEmpty()){
      Error("CrossSectionPlots","Could not get histogram ROOT file or it was empty.");
      return;
    }
    ret_hist = (MnvH3D*)f->Get(hist_name.c_str());

    if ( !ret_hist ){
      Error("CrossSectionPlots",TString::Format("No histogram found with name %s",hist_name.c_str()));
      return;
    }
  } */
  template< class MnvH > MnvH1D* TransWarpExtraction<MnvH>::LinearizeHist( MnvH* hist )
  {
    MnvH1D* linear_hist = new MnvH1D(Form("%s_Linear",hist->GetName()),Form(";%s;%s",hist->GetXaxis()->GetTitle(),hist->GetYaxis()->GetTitle()),hist->fN, 0.0, (double)hist->fN);
    for( int i=0; i<hist->fN; ++i )
    {
      linear_hist->SetBinContent(i+1,hist->GetBinContent(i));
      linear_hist->SetBinError(i+1,hist->GetBinError(i));
    }
    linear_hist->GetXaxis()->SetTitle(Form("Bin Number (%s)",linear_hist->GetXaxis()->GetTitle()));
    return linear_hist;      
  }

  template< class MnvH > MnvH* TransWarpExtraction<MnvH>::throwStat(MnvH* hist, TRandom3 *gen, int myuni, std::string prefix)
  {
    //Do the original Throwing method (fake data ONLY)
    //This is the default behavior
    if(m_stat_scale<0){
      std::cout << "Staring to throw stat (fakedata only!!) " << myuni << std::endl;
      MnvH* myhist = (MnvH*)hist->Clone(Form("%s-Input-%d",prefix.c_str(),myuni));
      myhist->Reset();
      
      for(int i=0;i<hist->fN;i++){
        double randomNumber = gen->PoissonD(hist->GetBinContent(i));
        myhist->SetBinContent(i,randomNumber);
        myhist->SetBinError(i,sqrt(randomNumber));
      }
      return myhist;
    }

    std::string classname = m_reco->ClassName();
    std::cout << "I'm working with classname\t" << classname << std::endl;
    //Make containers for the new throws
    MnvH* myhist = (MnvH*)hist->Clone(Form("%s-Input-%d",prefix.c_str(),myuni));
    MnvH* myhist_truth = (MnvH*)m_truth->Clone(Form("%s-Input-%d_truth",prefix.c_str(),myuni));
    MnvH* myhist_bkg = (MnvH*)m_reco->Clone(Form("%s-Input-%d_bkg",prefix.c_str(),myuni));
    myhist->Reset();
    myhist_truth->Reset();
    MnvH1D * my_hist_signal = m_migration_original->ProjectionX("_px",0,-1);
    for (int i=0;i<myhist_bkg->fN;++i) {
      myhist_bkg->AddBinContent(i,-1.0*my_hist_signal->GetBinContent(i));
    }
    delete my_hist_signal;
    my_hist_signal=nullptr;


    //Populate the migration matrix with the original migration matrix
    for(int i=0;i<m_migration->GetNbinsX()+2;i++){
      for(int j=0;j<m_migration->GetNbinsY()+2;j++){
	m_migration->SetBinContent(i,j,m_migration_original->GetBinContent(i,j));
	m_migration->SetBinError(i,j,m_migration_original->GetBinError(i,j));
      }
    }

    //Now we do the throwing
    //This agnostic to 1D or 2D results.
    for(int i=0;i<m_migration->GetNbinsX()+2;i++){
      for(int j=0;j<m_migration->GetNbinsY()+2;j++){
	double randomNumber = 0;
	double content = m_migration->GetBinContent(i,j);
	//This is the more intelligent throwing?
	//default configuration doesn't run this mode
	if(doTailAvgThrows){
	 if(content>40)randomNumber = gen->PoissonD(content*m_stat_scale)/m_stat_scale;
	 else{
	   double content_ineg_j = i-1>=0?m_migration->GetBinContent(i-1,j):0;
	   double content_ineg_jpos = i-1>=0&&j+1<m_migration->GetNbinsY()+1?m_migration->GetBinContent(i-1,j+1):0;
	   double content_ineg_jneg = i-1>=0&&j-1>=0?m_migration->GetBinContent(i-1,j-1):0;
	   double content_i_j = m_migration->GetBinContent(i,j);
	   double content_i_jpos = j+1<m_migration->GetNbinsY()+1?m_migration->GetBinContent(i,j+1):0;
	   double content_i_jneg = j-1>=0?m_migration->GetBinContent(i,j-1):0;
	   double content_ipos_j = i+1<m_migration->GetNbinsX()+1?m_migration->GetBinContent(i+1,j):0;
	   double content_ipos_jpos = i+1<m_migration->GetNbinsX()+1?m_migration->GetBinContent(i+1,j+1):0;
	   double content_ipos_jneg = i+1<m_migration->GetNbinsX()+1&&j-1>=0?m_migration->GetBinContent(i+1,j-1):0;

	   double avg = (content_ineg_j+content_ineg_jpos+content_ineg_jneg+content_i_j+content_i_jpos+content_i_jneg+content_ipos_j+content_ipos_jpos+content_ipos_jneg)/9;
	   randomNumber = gen->PoissonD(avg*m_stat_scale)/m_stat_scale;
	 }
	}
	else if (m_stat_scale>0){
	  randomNumber=gen->PoissonD(content*m_stat_scale)/m_stat_scale;
	} else {
    randomNumber = content;
  }
	m_migration->SetBinContent(i,j,randomNumber);
	m_migration->SetBinError(i,j,sqrt(randomNumber));
      }
    }
    double n_model_truth=0,n_data_truth=0;
    for (int i=0;i<m_truth->fN;++i) {
      n_model_truth += m_truth->GetBinContent(i);
      n_data_truth += m_data_truth->GetBinContent(i);
    }
    double scale = n_data_truth/n_model_truth;
    if(classname=="PlotUtils::MnvH2D"){
      std::cout << "2D analysis stat throw (migration plus fakedata)" << std::endl;
    } else if(classname=="PlotUtils::MnvH1D"){
        //1D      
        //1D analyses don't use an encoded space so they just use ufof and standard bins
      std::cout << "1D analysis stat throw (migration plus fakedata)" << std::endl;
    } else {
      std::cout << "You are trying to run in 3D mode which is not supported" << std::endl;
      exit(1);
    }
    for(int i=0;i<m_migration->GetNbinsX()+2;i++){//reco axis
      for(int j=0;j<m_migration->GetNbinsY()+2;j++){//truth axis
        myhist_truth->AddBinContent(j,m_migration->GetBinContent(i,j));
      }
    }
    myhist_truth->Divide(myhist_truth,m_data_truth);
    for(int i=0;i<m_migration->GetNbinsX()+2;i++){//reco
      for(int j=0;j<m_migration->GetNbinsY()+2;j++){//truth
        double correction = myhist_truth->GetBinContent(j)>0?1/myhist_truth->GetBinContent(j):1;
        myhist->AddBinContent(i,m_migration->GetBinContent(i,j)*correction);
      }
    }
    for(int i=0;i<hist->fN;i++){
      myhist->SetBinError(i,sqrt(myhist->GetBinContent(i)));
    }
    myhist->SetEntries(1);
    myhist_bkg->Scale(scale);
    //Now throw the fake data again
    for(int i=0;i<hist->fN;i++){
      //Add bkg to myhist.
      double randomNumber = gen->PoissonD(myhist->GetBinContent(i)+myhist_bkg->GetBinContent(i));
      myhist->SetBinContent(i,randomNumber);
      myhist->SetBinError(i,sqrt(randomNumber));
    }

    //Now put the migration matrix back in the original configuration for unfolding the thrown universe.
    for(int i=0;i<m_migration->GetNbinsX()+2;i++){
      for(int j=0;j<m_migration->GetNbinsY()+2;j++){
	double spec_err = sqrt(m_migration->GetBinContent(i,j)*m_stat_scale);
	m_migration->SetBinContent(i,j,m_migration_original->GetBinContent(i,j));
	m_migration->SetBinError(i,j,spec_err);
      }
    }

    return myhist;
  }

  /*This is the original throwing mechanics.
  template< class MnvH > MnvH* TransWarpExtraction<MnvH>::throwStat(MnvH* hist, TRandom3 *gen, int myuni, std::string prefix)
  {
    std::cout << "Staring to throw stat " << myuni << std::endl;
    MnvH* myhist = (MnvH*)hist->Clone(Form("%s-Input-%d",prefix.c_str(),myuni));
    myhist->Reset();
    
    for(int i=0;i<hist->fN;i++){
      double randomNumber = gen->PoissonD(hist->GetBinContent(i));
      myhist->SetBinContent(i,randomNumber);
      myhist->SetBinError(i,sqrt(randomNumber));
    }
    return myhist;
  }
  */
  
  
}

using namespace PlotUtils;
int main( int argc, char **argv)
{
  #ifndef MNVROOT6
  #ifndef NCINTEX
  ROOT::Cintex::Cintex::Enable();
  #endif
  #endif
  TH1::AddDirectory(false);

  const char* const short_options = "o:D:d:I:i:M:m:R:r:T:t:n:lLu:z:c:C:p:P:f:x:s:vVhb:g:a";
  static struct option long_options[]=
  {
    {"output_file",     required_argument, nullptr ,'o'},
    {"data_file",       required_argument, nullptr, 'D'},
    {"data",            required_argument, nullptr, 'd'},
    {"data_truth_file", required_argument, nullptr, 'I'},
    {"data_truth",      required_argument, nullptr, 'i'},
    {"migration_file",  required_argument, nullptr, 'M'},
    {"migration",       required_argument, nullptr, 'm'},
    {"reco_file",       required_argument, nullptr, 'R'},
    {"reco",            required_argument, nullptr, 'r'},
    {"truth_file",      required_argument, nullptr, 'T'},
    {"truth",           required_argument, nullptr, 't'},
    {"num_iter",        required_argument, nullptr, 'n'},
    {"log_scale",       no_argument,       nullptr, 'l'},
    {"linearize",       no_argument,       nullptr, 'L'},
    {"num_uni" ,        required_argument, nullptr, 'u'},
    {"num_dim" ,        required_argument, nullptr, 'z'},
    {"max_chi2",        required_argument, nullptr, 'c'},
    {"step_chi2",       required_argument, nullptr, 'C'},
    {"pot_scale",       required_argument, nullptr, 'p'},
    {"data_pot_norm",   required_argument, nullptr, 'P'},
    {"stat_scale",      required_argument, nullptr, 'f'},
    {"corr_factor",     required_argument, nullptr, 'F'},
    {"exclude_bins",    required_argument, nullptr, 'x'},
    {"random_seed",     required_argument, nullptr, 's'},
    {"verbhist",        no_argument,       nullptr, 'V'},
    {"help",            no_argument,       nullptr, 'h'},
    {"target_string",   required_argument,       nullptr, 'g'},
    {"variable",        required_argument,       nullptr, 'a'},
    {"PDG",             required_argument,       nullptr, 'y'},
    {nullptr,           0,                 nullptr, 0}
  };

  std::string output_file     ; bool opt_output_file     = false; 
  std::string data            ; bool opt_data            = false; 
  std::string data_file       ; bool opt_data_file       = false; 
  std::string data_truth      ; bool opt_data_truth      = false; 
  std::string data_truth_file ; bool opt_data_truth_file = false; 
  std::string migration       ; bool opt_migration       = false; 
  std::string migration_file  ; bool opt_migration_file  = false; 
  std::string reco            ; bool opt_reco            = false; 
  std::string reco_file       ; bool opt_reco_file       = false; 
  std::string truth           ; bool opt_truth           = false; 
  std::string truth_file      ; bool opt_truth_file      = false; 
  std::string num_iter        ;
  std::string exclude_bins    ;
  std::string target_string   ;
  int PDG = 14 ;
  std::string variable        ;
  bool bIterLogScale = false;
  bool fakes         = false;
  bool bLinearize    = false;
  int num_uni  = 1;
  int num_dim  = 1;
  double max_chi2      = 5000;
  double step_chi2     = 50;
  double pot_scale     = 1.0;
  double data_pot_norm = 1.0;
  double stat_scale    = -999; //Negative means standard throwing
  double corr_factor   = -999; //Negative means no correction in covariance matrix and unfolded histogram
  int random_seed      = -1;
  bool verbhist        = false;
  std::cout<<"Tmp1\n";
  int cc;
  while ((cc = getopt_long(argc, argv, short_options, long_options, nullptr)) != -1) {
    std::cout << "cc: " << (char)cc << "\t optarg: " << optarg <<std::endl;
    switch (cc)
    {
      case 'o': output_file     = std::string(optarg) ; opt_output_file     = true; break;
      case 'd': data            = std::string(optarg) ; opt_data            = true; break;
      case 'D': data_file       = std::string(optarg) ; opt_data_file       = true; break;
      case 'i': data_truth      = std::string(optarg) ; opt_data_truth      = true; break;
      case 'I': data_truth_file = std::string(optarg) ; opt_data_truth_file = true; break;
      case 'm': migration       = std::string(optarg) ; opt_migration       = true; break;
      case 'M': migration_file  = std::string(optarg) ; opt_migration_file  = true; break;
      case 'r': reco            = std::string(optarg) ; opt_reco            = true; break;
      case 'R': reco_file       = std::string(optarg) ; opt_reco_file       = true; break;
      case 't': truth           = std::string(optarg) ; opt_truth           = true; break;
      case 'T': truth_file      = std::string(optarg) ; opt_truth_file      = true; break;
      case 'n': num_iter        = std::string(optarg) ; break;
      case 'l': bIterLogScale   = true ; break;
      case 'L': bLinearize      = true ; break;
      case 'u': num_uni         = atoi(optarg); break;
      case 'z': num_dim         = atoi(optarg); break;
      case 'c': max_chi2        = atof(optarg); break;
      case 'C': step_chi2       = atof(optarg); break;
      case 'p': pot_scale       = atof(optarg); break;
      case 'P': data_pot_norm   = atof(optarg); break;
      case 'f': stat_scale      = atof(optarg); break;
      case 'F': corr_factor     = atof(optarg); break;
      case 'V': verbhist        = true; break;
      case 'x': exclude_bins    = std::string(optarg); break;
      case 's': random_seed     = atoi(optarg); break;
      case 'b' : fakes          = true; break;
      case '?':
      case 'g' : target_string  = std::string(optarg); break;
      case 'a' : variable       = std::string(optarg); break;
      case 'y' : PDG            = atoi(optarg); break;
      case 'h':
      default:
        std::cout << "|********************************* Options *****************************************|" << std::endl
             << "| This program outputs the unfolded data distribution, as well as various           |" << std::endl
             << "| chi squared metrics.  See Dan Ruterbories talk: docdb 16130                       |" << std::endl
             << "|                                                                                   |" << std::endl
             << "| This can unfold 1,2, and 3 dimensional analyses                                   |" << std::endl
             << "| Recommend using MnvResponse for dimensions > 1                                    |" << std::endl
             << "| The program needs four histograms: data, mc reco, mc truth, and migration matrix  |" << std::endl
             << "| Further instructions here:                                                        |" << std::endl
             << "| https://minerva-docdb.fnal.gov/cgi-bin/private/ShowDocument?docid=18488           |" << std::endl
             << "|                                                                                   |" << std::endl
             << "| 'Data' is a warped MC distributions                                               |" << std::endl
             << "| All these histograms are signal only                                              |" << std::endl
             << "|                                                                                   |" << std::endl
             << "| Options:                                                                          |" << std::endl
             << "|   -o, --output_file     : (Required) Name of output file                          |" << std::endl
             << "|   -d, --data            : (Required) MnvH*D Name of data histogram                |" << std::endl
             << "|   -D, --data_file       : (Required) File with data histogram                     |" << std::endl
             << "|   -i, --data_truth      : (Required) MnvH*D Name of data histogram                |" << std::endl
             << "|   -I, --data_truth_file : (Required) File with data histogram                     |" << std::endl
             << "|   -m, --migration       : (Required) MnvH*D Name of migration histogram           |" << std::endl
             << "|                         : between MC reco and MC truth                            |" << std::endl
             << "|   -M, --migration_file  : (Required) File with migration histogram                |" << std::endl
             << "|   -r, --reco            : (Required) MnvH*D Name of MC reco histogram             |" << std::endl
             << "|   -R, --reco_file       : (Required) File with MC reco histogram                  |" << std::endl
             << "|   -t, --truth           : (Required) MnvH*D Name of MC truth histogram            |" << std::endl
             << "|   -T, --truth_file      : (Required) File with C truth histogram                  |" << std::endl
             << "|   -n, --num_iter        : Number of times to unfold data histogram                |" << std::endl
             << "|                         : Can input multiple iterations, comma separated          |" << std::endl
             << "|   -l, --log_scale       : Logarithmic binning for iterations (Default: linear)    |" << std::endl
             << "|   -L, --linearize       : Output folders of 'flattened' histograms  where         |" << std::endl
             << "|                         : each bin is the one bin on a 1D hist (under/overflow    |" << std::endl
             << "|                         : included).  Bin width is equal                          |" << std::endl
             << "|   -u, --num_uni         : Number of statistical universes       (Default: 1)      |" << std::endl
             << "|   -z, --num_dim         : Number of dimensions in this analysis (Default: 1)      |" << std::endl
             << "|   -c, --max_chi2        : Maximum Chi2 on 2D plots (Default: 5000)                |" << std::endl
             << "|   -C, --step_chi2       : Size of binning of Chi2  (Default: 50)                  |" << std::endl
             << "|   -p, --pot_scale       : POT scaling between data and MC (Default: 1)            |" << std::endl
             << "|                         : You don't need to touch this unless you're looking at   |" << std::endl
             << "|                         : chi2 comparisons between data and MC                    |" << std::endl
             << "|   -P, --data_pot_norm   : POT norm of data (playlist norm) (Default: 1)           |" << std::endl
             << "|                         : You'll want to use this to scale your data              |" << std::endl
             << "|                         : to the entire ME playlist, for example                  |" << std::endl
	     << "|   -f, --stat_scale      : Used to rescale unctertainty in migration matrix bins   |" << std::endl
	     << "|   -F, --corr_factor     : Used to apply a correction to the covariance matrix     |" << std::endl
	     << "|                         : and the diagonal uncertainty of the unfolded histogram  |" << std::endl
             << "|   -x, --exclude_bins    : List of global bins to exclude from chi2 calculation    |" << std::endl
             << "|                         : Can input bins, comma separated                         |" << std::endl
             << "|                         : Use FindBin(x,y,z) to get global bin number             |" << std::endl
             << "|   -s, --random_seed     : Can set a consistent random seed                        |" << std::endl
             << "|   -V, --verbHist        : Writes out some debug histograms                        |" << std::endl
             << "|   -g, --target_string   : Which target over which to perform warping study        |" << std::endl
             << "|   -a, --variable        : Which variable over which to perform warping study      |" << std::endl
             << "|   -y, --PDG             : Which PDG over which to perform warping study (Default: 14)     |" << std::endl
             << "|   -h, --help            : Print this.                                             |" << std::endl
             << "|***********************************************************************************|" << std::endl;
         return 0;
    }
  }
  std::cout<<"Tmp2\n";
  if( !opt_output_file )     { std::cout<<"Need an output file (-o, --output_file )"<<std::endl           ; return 1 ; }
  if( !opt_data        )     { std::cout<<"Need a data histogram name (-d, --data )"<<std::endl           ; return 1 ; }
  if( !opt_data_file  )      { std::cout<<"Need an data file (-D, --data_file )"<<std::endl               ; return 1 ; }
  if( !opt_data_truth )      { std::cout<<"Need a data truth histogram name (-i, --data_truth )"<<std::endl           ; return 1 ; }
  if( !opt_data_truth_file ) { std::cout<<"Need an data truth file (-I, --data_truth_file )"<<std::endl               ; return 1 ; }
  if( !opt_migration   )     { std::cout<<"Need a migration histogram name (-m, --migration )"<<std::endl ; return 1 ; }
  if( !opt_migration_file  ) { std::cout<<"Need an migration file (-M, --migration_file )"<<std::endl     ; return 1 ; }
  if( !opt_reco        )     { std::cout<<"Need a mc reco histogram name (-r, --reco )"<<std::endl        ; return 1 ; }
  if( !opt_reco_file  )      { std::cout<<"Need an reco file (-R, --reco_file )"<<std::endl               ; return 1 ; }
  if( !opt_truth       )     { std::cout<<"Need a mc truth histogram name (-t, --truth )"<<std::endl      ; return 1 ; }
  if( !opt_truth_file  )     { std::cout<<"Need an truth file (-T, --truth_file )"<<std::endl             ; return 1 ; }
  std::cout<<"Tmp3\n";
  std::cout << " Output File:               " << output_file    << std::endl
            << " Data File:                 " << data_file       << std::endl
            << " Data Histogram Name:       " << data            << std::endl
            << " Data Truth File:           " << data_truth_file << std::endl
            << " Data Truth Histogram Name: " << data_truth      << std::endl
            << " Migration File:            " << migration_file  << std::endl
            << " Migration Histogram Name:  " << migration       << std::endl
            << " MC Reco File:              " << reco_file       << std::endl
            << " MC Reco Histogram Name:    " << reco            << std::endl
            << " MC Truth File:             " << truth_file      << std::endl
            << " MC Truth Histogram Name:   " << truth           << std::endl
            << " Number of iterations:      " << num_iter        << std::endl
            << " Iteration Log Scale:       " << (int)bIterLogScale << std::endl
            << " Linearized Hists:          " << (int)bLinearize    << std::endl
            << " Number of Stat Universes:  " << num_uni        << std::endl
            << " Number of Dimensions:      " << num_dim        << std::endl
            << " Max Chi2:                  " << max_chi2       << std::endl
            << " Step Chi2:                 " << step_chi2      << std::endl
            << " Data POT/MC POT:           " << pot_scale      << std::endl
            << " Data POT Scale:            " << data_pot_norm  << std::endl
            << " Stat. Unc. Scale:          " << stat_scale     << std::endl
            << " Correction Factor:         " << corr_factor    << std::endl
            << " Excluded Bins:             " << exclude_bins   << std::endl
            << " Random Seed (if -1 random) " << random_seed    << std::endl
            << " With Fakes?                " << fakes          << std::endl
            << " Target                     " << target_string  << std::endl
            << " PDG                        " << PDG            << std::endl
            << " Variable                   " << variable       << std::endl;
  
  //Splicing up the number of iterations
  std::vector<std::string> str_iterations;
  std::vector<int> iterations;
  split(num_iter,str_iterations,',');
  for(uint i = 0; i<str_iterations.size(); ++i) iterations.push_back( atoi(str_iterations[i].c_str()) );

  //Splicing up the number of bins to exclude
  std::vector<std::string> str_exclude_bins;
  std::vector<int> exclude_chi2_bins;
  split(exclude_bins,str_exclude_bins,',');
  for(uint i = 0; i<str_exclude_bins.size(); ++i) exclude_chi2_bins.push_back( atoi(str_exclude_bins[i].c_str()) );
  TH1::AddDirectory(kFALSE);
  switch(num_dim){
    case 1: 
    {
      TransWarpExtraction< MnvH1D > excelsior(data_file, data, data_truth_file, data_truth, reco_file, reco, truth_file, truth, migration_file, migration, iterations, exclude_chi2_bins, random_seed, num_uni, bIterLogScale, target_string, variable, PDG);
      excelsior.m_fakes=fakes;
      excelsior.ClearInputErrorBands();
      excelsior.ScalePOT( pot_scale ); 
      excelsior.ScaleDataPOT( data_pot_norm ); 
      //      excelsior.ScaleMigStatUnc( stat_scale );
      excelsior.SetChi2MaxAndStep( max_chi2, step_chi2 );
      excelsior.SetStatScale(stat_scale);
      excelsior.SetCorrFactor(corr_factor);
      excelsior.Setup2DChi2Hists();
      excelsior.UnfoldStatUniverses( );
      excelsior.CalcChi2( );
      excelsior.MakeBinChi2Dists( );
      excelsior.MakeMedianChi2Dists( );
      excelsior.MakeTruncatedChi2Dists( );
      excelsior.MakeErrorHists( );
      excelsior.MakeRatioHists( );
      excelsior.WriteOutput( output_file, verbhist, bLinearize ); 
      return 0;
    }
    case 2: 
    {
      /* TransWarpExtraction< MnvH2D > excelsior(data_file, data, data_truth_file, data_truth, reco_file, reco, truth_file, truth, migration_file, migration, iterations, exclude_chi2_bins, random_seed, num_uni, bIterLogScale, target_string, variable, PDG);
      excelsior.m_fakes=fakes;
      excelsior.ClearInputErrorBands();
      excelsior.ScalePOT( pot_scale ); 
      excelsior.ScaleDataPOT( data_pot_norm );
      //      excelsior.ScaleMigStatUnc( stat_scale ); 
      excelsior.SetChi2MaxAndStep( max_chi2, step_chi2 );
      excelsior.SetStatScale(stat_scale);
      excelsior.SetCorrFactor(corr_factor);
      excelsior.Setup2DChi2Hists();
      excelsior.UnfoldStatUniverses( );
      excelsior.CalcChi2( );
      excelsior.MakeBinChi2Dists( );
      excelsior.MakeMedianChi2Dists( );
      excelsior.MakeTruncatedChi2Dists( );
      excelsior.MakeErrorHists( );
      excelsior.MakeRatioHists( );
      excelsior.WriteOutput( output_file, verbhist, bLinearize );  */
      return 0;
    }
    //TODO: This block doesn't compile because MnvH3 has no GetBin(int, int) method
    case 3: 
    {
      std::cerr << "MnvH3D not support by TransWarpExtractor because Andrew doesn't know how to fix the GetBin(int, int) signatures.  Compile it yourself!\n";
      return 1;

      /*TransWarpExtraction< MnvH3D > excelsior(data_file, data, data_truth_file, data_truth, reco_file, reco, truth_file, truth, migration_file, migration, iterations, exclude_chi2_bins, random_seed, num_uni, bIterLogScale );
      excelsior.ClearInputErrorBands();
      excelsior.ScalePOT( pot_scale ); 
      excelsior.ScaleDataPOT( data_pot_norm ); 
      excelsior.SetChi2MaxAndStep( max_chi2, step_chi2 );
      excelsior.SetStatScale(stat_scale);
      excelsior.SetCorrFactor(corr_factor);
      excelsior.Setup2DChi2Hists();
      excelsior.UnfoldStatUniverses( );
      //excelsior.CalcChi2( );//This doesn't exist for 3d yet.  So until then...
      //excelsior.MakeBinChi2Dists( );
      //excelsior.MakeMedianChi2Dists( );
      //excelsior.MakeTruncatedChi2Dists( );
      excelsior.MakeErrorHists( );
      excelsior.MakeRatioHists( );
      excelsior.WriteOutput( output_file, verbhist, bLinearize ); 
      return 0;*/
    }
  }
}
