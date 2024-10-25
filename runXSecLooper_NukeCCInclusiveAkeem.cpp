#include "GENIEXSecExtract/XSecLooper.h"

#include <cstdlib>
#include "TFile.h"
#include "TH1.h"
#include <PlotUtils/MnvH1D.h>
#include <PlotUtils/TargetUtils.h>

#include "NukeCCInclusive_bins.cc"

using namespace std;
using namespace PlotUtils;

namespace {
    const double TRACKER_ZMIN = 6117;
    const double TRACKER_ZMAX = 8193;
}

double GetNormValue( int targetID, int targetZ, int targetNucleon = 0 )
{
    
   double trackerAtomsC = TargetUtils::Get().GetTrackerElementNAtoms( 6, 108, true );
    
    cout<<trackerAtomsC<<"  "<<2.22311e+27 * 108.<<endl;
    
    
    //    const double trackerAtomsC = 2.22311e+27 * 92.;
    //     const double trackerAtomsC = 2.22311e+27 * 108.;
    double passiveNucleons = -999.;
    
    if ( targetNucleon == 0 ){
        
        if(targetID < 10 ){
            
            passiveNucleons = TargetUtils::Get().GetPassiveTargetNNucleons(targetID,targetZ, true);
            
        }
        if(targetID > 10){
            passiveNucleons = TargetUtils::Get().GetTrackerNNucleons( 12, true, 850. );
        }
        
    }

    else
    {
        assert( false && "Target nucleons can be all, proton or neutron only" );
    }
    
    if( passiveNucleons < 0 )
        assert( false && "Normalizations only known for Z = 1,6,26,82 and 0(scint)" );
    
    cout<<"The normalization factor for this analysis is "<<trackerAtomsC / passiveNucleons<<endl;
    
    return trackerAtomsC / passiveNucleons;
}



//======================================
//! NukeCC XSec
//======================================
class NukeCCXSec : public XSec
{
    
public:
    
    NukeCCXSec( const char* name, int target = 0, int targetZ = 0, bool isCCInclusive=0, bool isDIS = 0, bool InHex = 0, bool PassTrueDistToDivisionCut = 0 ) :
    XSec( name ),
    m_target( target ),
    m_nucleus( targetZ ),
    m_isCCInclusive( isCCInclusive),
    m_isDIS( isDIS),
    m_isinHex( InHex),
    m_truedisttoDiv(PassTrueDistToDivisionCut ){};
    
    double Qsquared(ChainWrapper& chw, int entry){
        const int PROTON_PDG = 2212;
        const int NEUTRON_PDG = 2112;
        const double M_neutron = 939.56;
        const double M_proton  = 938.27;
        const double M_nucleon  = ( 1.5*M_neutron + M_proton ) / 2.5;
        
        const double PI = atan(1.0)*4;
        const double deg_to_rad = PI / 180.0;
        const double rad_to_deg = 1. / deg_to_rad;
        
        bool signal = true;
        
        //good muon
        double muon_E     = (double)chw.GetValue("mc_primFSLepton",entry,3);
        double muon_theta = (double)chw.GetValue("truth_muon_theta",entry);
        
        double neutrino_E = (double)chw.GetValue("mc_incomingE",entry);
        double Q2 = 4*neutrino_E*muon_E*pow( sin( muon_theta /2 ), 2. )/(1000*1000);
        //        cout<<"The Q2 is "<<Q2<<endl;
        return Q2;
    }
    
    double W_had(ChainWrapper& chw, int entry){
        const int PROTON_PDG = 2212;
        const int NEUTRON_PDG = 2112;
        const double M_neutron = 939.56;
        const double M_proton  = 938.27;
        const double M_nucleon  = ( 1.5*M_neutron + M_proton ) / 2.5;
        
        const double PI = atan(1.0)*4;
        const double deg_to_rad = PI / 180.0;
        const double rad_to_deg = 1. / deg_to_rad;
        
        int target_nucleon = (int)chw.GetValue("mc_targetNucleon",entry);
        
        //good muon
        double muon_E     = (double)chw.GetValue("mc_primFSLepton",entry,3);
        double muon_theta = (double)chw.GetValue("truth_muon_theta",entry);
        double neutrino_E = (double)chw.GetValue("mc_incomingE",entry);
        double Q2 = Qsquared(chw,entry)*1000*1000;
        double W = -1.0;
        double nuclMass = M_nucleon;
        if( NEUTRON_PDG  == target_nucleon )
            nuclMass = M_neutron;
        else if( PROTON_PDG == target_nucleon )
            nuclMass = M_proton;
        double W2 = ( pow(nuclMass, 2) +  2. * ( neutrino_E -  muon_E ) * nuclMass - Q2);
        
        if(W2 > 0){
            W = sqrt(W2)/1000.;
            //            cout<<"The W value is "<<W<<endl;
        }
        else{
            return -999.; //if W2 < 0, the event is not DIS...
        }
        return W;
    }
    
    bool InHex( ChainWrapper& chw, int entry, double apothem )
    {
        if( apothem == 0. )
            return false;
        
        double x = (double)chw.GetValue("mc_vtx",entry,0);
        double y = (double)chw.GetValue("mc_vtx",entry,1);
        
        //        cout<<x<<"  "<<y<<endl;
        
        if( pow(x,2) + pow(y,2) < pow(apothem,2) )
            return true;
        
        //Hexagon is symmetric about its x and y
        x = fabs(x);
        y = fabs(y);
        
        double lenOfSide = apothem * ( 2 / sqrt(3) );
        
        if( x > apothem )
            return false;
        
        if( y < lenOfSide/2.0 )
            return true;
        
        double slope = (lenOfSide / 2.0) / apothem;
        if( y < lenOfSide - x*slope )
            return true;
        
        return false;
    }
    
    bool PassTrueDistToDivisionCut( ChainWrapper& chw, int entry, int targetID, int num_targetZ, double xySep /* = 25. */ )
    {
        
        double true_target_dist_to_division = (double)chw.GetValue("truth_target_dist_to_division",entry);
        
        
        int true_module   = (int)chw.GetValue("truth_vtx_module",entry);
        int true_targetID = (int)chw.GetValue("truth_targetID",entry);
        int true_targetZ  = (int)chw.GetValue("truth_targetZ",entry);
        //cout << true_targetID <<" "<< true_target_dist_to_division <<" "<< true_module << endl;
        if( targetID < 10 && true_targetID != targetID) return false;
                
        if( targetID < 10 && true_targetZ != num_targetZ) return false;
        if( targetID > 10 && ( true_module  < 27 || true_module > 32)) return false;
        if( targetID == 24 && ( true_module < 33 || true_module > 38)) return false;
        if( targetID == 34 && ( true_module < 39 || true_module > 44)) return false;
        if( targetID == 44 && ( true_module < 45 || true_module > 50)) return false;
        if( targetID == 54 && ( true_module < 51 || true_module > 56)) return false;
        if( targetID == 64 && ( true_module < 57 || true_module > 62)) return false;
        if( targetID == 74 && ( true_module < 63 || true_module > 68)) return false;
        if( targetID == 84 && ( true_module < 69 || true_module > 74)) return false;
        if( targetID == 94 && ( true_module < 75 || true_module > 80)) return false;
        
        
        //            only relevant for passive targets 1235
        if( 0 < true_targetID && true_targetID < 10 && 4 != true_targetID)
            
            return ( xySep < true_target_dist_to_division );
               //cout<<"fiducial event "<<targetID<<endl;
       
        return true;
    }
    
    
    
    bool isCCInclusive(ChainWrapper& chw, int entry){
    
        const double PI = atan(1.0)*4;
        const double deg_to_rad = PI / 180.0;
        const double rad_to_deg = 1. / deg_to_rad;
        
        //bool signal = true;
        
        //right particle
        int is_neutrino = (int)chw.GetValue("mc_incoming",entry);
        //        cout<<"am I a neutrino? "<<is_neutrino<<endl;
        if(is_neutrino!=14) return false;
        
        //right interaction type
        int is_CC     = (int)chw.GetValue("mc_current",entry);
        //        cout<<"am I CC? "<<is_CC<<endl;
        if(is_CC!=1) return false;
        
        //right genie kindmaitcs
        //int genie_dis = (int)chw.GetValue("mc_intType",entry);
        //        cout<<"am I a GENIE DIS Event? "<<genie_dis<<endl;
        //if(genie_dis==1 || genie_dis==2) return false;
        
        
        double muon_E     = (double)chw.GetValue("mc_primFSLepton",entry,3);
        double muon_theta = (double)chw.GetValue("truth_muon_theta",entry);
        //        cout<<"My muons are "<<muon_E<<"  "<<muon_theta<<endl;
        if( muon_E <= 2000.0 || muon_E > 50000.0 ) return false;
        //         if( muon_E < 2000.0) return false;
        if( muon_theta*rad_to_deg >= 17.0 || muon_theta*rad_to_deg < 0.0 ) return false;
    
        return true;
    }
    
    //! Return QE-like signal
    bool isDIS( ChainWrapper& chw, int entry ) {
        
       bool signal = true;
        const double PI = atan(1.0)*4;
        const double deg_to_rad = PI / 180.0;
        const double rad_to_deg = 1. / deg_to_rad;
        
       /* 
        //right particle
        int is_neutrino = (int)chw.GetValue("mc_incoming",entry);
        //        cout<<"am I a neutrino? "<<is_neutrino<<endl;
        if(is_neutrino!=14) return false;
        
        //right interaction type
        int is_CC     = (int)chw.GetValue("mc_current",entry);
        //        cout<<"am I CC? "<<is_CC<<endl;
        if(is_CC!=1) return false;
        
        //right genie kindmaitcs
        int genie_dis = (int)chw.GetValue("mc_intType",entry);
        //        cout<<"am I a GENIE DIS Event? "<<genie_dis<<endl;
        if(genie_dis==1 || genie_dis==2) return false;
        
        
        double muon_E     = (double)chw.GetValue("mc_primFSLepton",entry,3);
        double muon_theta = (double)chw.GetValue("truth_muon_theta",entry);
        //        cout<<"My muons are "<<muon_E<<"  "<<muon_theta<<endl;
        if( muon_E <= 2000.0 || muon_E > 50000.0 ) return false;
        //         if( muon_E < 2000.0) return false;
        if( muon_theta*rad_to_deg >= 17.0 || muon_theta*rad_to_deg < 0.0 ) return false;
        
       */ 
        
        double Q2 = Qsquared(chw, entry);
        double W  = W_had(chw, entry);
        
        if(Q2 < 1.0) return false;
        if(W < 2.0)  return false;
        
        return signal;
    }
    
    //! Override this method from the base class to decide what events to include in this selection
    virtual bool passesCuts( ChainWrapper& chw, int entry ) {
        
        //        cout<<chw.GetValue("mc_primFSLepton",entry,0)/1000.<<"  "<<chw.GetValue("mc_primFSLepton",entry,1)/1000.<<"  "<<chw.GetValue("mc_primFSLepton",entry,2)/1000.<<"  "<<chw.GetValue("mc_primFSLepton",entry,3)/1000.<<endl;
        
        //! pass signal
        if(!isCCInclusive(chw, entry)) return false;
       // if(!isDIS(chw,entry)) return false;
        if(!InHex(chw,entry,850.0)) return false;
        if(!PassTrueDistToDivisionCut( chw, entry, m_target, m_nucleus, 25.0)) return false;
        
        //this stuff is nice to have for event by event checking... 
        int run = (int)chw.GetValue("mc_run",entry);
        int subrun = (int)chw.GetValue("mc_subrun",entry);
        int gate = (int)chw.GetValue("mc_nthEvtInFile",entry);
        
        const double PI = atan(1.0)*4;
        const double deg_to_rad = PI / 180.0;
        const double rad_to_deg = 1. / deg_to_rad;
        
        
//                cout<<run<<" , "<<subrun<<" , "<<gate<<" , ";
//                cout.precision(5);
//                cout<<Qsquared(chw,entry)<<" , ";
//                cout<<W_had(chw,entry)<<" , ";
//                cout<<chw.GetValue("mc_primFSLepton",entry,3)/1000.<<" , ";
//                cout<<chw.GetValue("truth_muon_theta",entry)*rad_to_deg<<" , ";
//        //        cout.precision(10);
//                cout<<chw.GetValue("mc_vtx",entry,0)<<" , ";
//                cout<<chw.GetValue("mc_vtx",entry,1)<<" , ";
//                cout<<chw.GetValue("mc_vtx",entry,2)<<" , ";
//        //        cout<<chw.GetValue("mc_vtx",entry,3)<<" , ";
//        //        cout.precision(2);
//                cout<<chw.GetValue("mc_incomingE",entry)/1000.<<endl;
        
        
        
        return true;
        
    } //! end function passesCuts
    
    //! input parameters
    int m_target;
    int  m_nucleus;
    bool m_isCCInclusive;
    bool m_isDIS;
    bool m_isinHex;
    bool m_truedisttoDiv;
    
};

//======================
//! main:
//======================
//void runXSecLooper_NukeCC( string& target_type )
void runXSecLooper_NukeCCInclusive(std::string playlistFile)
{
    TH1::AddDirectory(0);
    
    //! get the ntuple directory names
    string directory   = "genie_xsection_ntuples";
    string output_name = "genie_xsections";
    
//    output_name += "_" + target_type;
    
    
    output_name += "_histos.root";
    
    cout << "    creating a genie file with name = " << output_name << endl;
    
    //! Create the XSecLooper and tell it the input files
    
//    XSecLooper loop("/minerva/data/users/norrick/NukeTuples/v1_WithTargetCode_MLNX/NukeTruth_AnaTuple_Tgt*minervame1D*.root");
    XSecLooper loop(playlistFile.c_str());
    
    //! Tell the XSecLooper which neutrino type we're considering (mandatory)
    loop.setNuPDG(14);
    
    //! Setting the number of Universes in the GENIE error band (default 100, put 0 if you do not want to include the universes)
    loop.setNumUniv(0); //Tammy used 100
    
    //! store cross section names
    const char* channel = "CC";
//    std::vector< string > targets;
//    if( target_type == "nuclei" ) {
//        targets.push_back( "carbon" );
//        targets.push_back( "iron" );
//        targets.push_back( "lead" );
//    } else { targets.push_back( "plastic" ); }
    
//    string processes[1] = { "dis"};
    
    vector<XSec::EVariable> vars;
    vars.push_back( XSec::kENu );
//vars.push_back( XSec::kEHad );
//    vars.push_back( XSec::kELep );
    vars.push_back( XSec::kxExp );
//    vars.push_back( XSec::ky );
////    //    //For Debugging
//    vars.push_back( XSec::kExpW );
//    vars.push_back( XSec::kExpQ2 );
//    vars.push_back( XSec::kThetaLep );
    
    vector<int> targetZs;
    targetZs.push_back(6);
    targetZs.push_back(26);
    targetZs.push_back(82);
    
    vector<int> targetIDs;
      targetIDs.push_back(1);
      targetIDs.push_back(2);
      targetIDs.push_back(3);
      targetIDs.push_back(4);
      targetIDs.push_back(5);
      targetIDs.push_back(14);
      targetIDs.push_back(24);
      targetIDs.push_back(34);
      targetIDs.push_back(44);
      targetIDs.push_back(54);
      targetIDs.push_back(64);
      targetIDs.push_back(74);
      targetIDs.push_back(84);
      targetIDs.push_back(94);
    
    for(int i=0; i <targetIDs.size(); ++i){
        
        for(int j=0; j <targetZs.size(); ++j){
            
            int num_target = targetIDs[i];
            int num_targetZ = targetZs[j];
            
            if(targetIDs[i]!=3 && targetZs[j]==6)
                continue;
            if(targetIDs[i]==4 && targetZs[j]!=82)
                continue;
            if(targetIDs[i]>10 && targetZs[j]!=82)
                continue;
            
            //! loop over the cross section extraction
            for( vector<XSec::EVariable>::iterator var = vars.begin(); var != vars.end(); ++var ){
                
                string varName = "";
                switch(*var)
                {
                    case XSec::kx:
                        varName = "xGen";
                        break;
                    case XSec::kxExp:
                        varName = "x";
                        break;
                    case XSec::kENu:
                        varName = "Enu";
                        break;
                    case XSec::kELep:
                        varName = "Emu";
                        break;
                    case XSec::kThetaLep:
                        varName = "ThetaMu";
                        break;
                    case XSec::kExpW:
                        varName = "W";
                        break;
                    case XSec::kExpQ2:
                        varName = "Q2";
                        break;
                    case XSec::kEHad:
                        varName = "Ehad";
                        break;
                    case XSec::ky:
                        varName = "y";
                        break;
                    default:
                        varName = "unknown";
                        break;
                }
                string nuc_name = "";
                bool makeDISHists = false;
                string sample = (makeDISHists) ? "dis" : "inc";
                if(num_targetZ==6) nuc_name = "carbon";
                if(num_targetZ==26) nuc_name = "iron";
                if(num_targetZ==82 && num_target< 10) nuc_name = "lead";
                if(num_targetZ==82 && num_target> 10) nuc_name = "tracker";
                
                
                char* name = Form( "%s_%d_%s_std_%s", nuc_name.c_str(), num_target, varName.c_str(), sample.c_str());
                cout<<"making cross section for "<<num_target<<"  "<<num_targetZ<<endl;
                NukeCCXSec* xsec = new NukeCCXSec( name, num_target, num_targetZ, false, makeDISHists, true, true );
                
                //! container for bins
                std::vector<double> bins;
                GetBins( *var, bins, false, false );
                int num_bins = (int)bins.size()-1;
                //! set cross section data
                xsec->setBinEdges(num_bins,&bins.front());
                xsec->setNormalizationValue( GetNormValue( num_target, num_targetZ ) );
                
                xsec->setVariable(*var);
                if(*var !=XSec::kENu) xsec->setIsFluxIntegrated(true);
                else  xsec->setIsFluxIntegrated(false);
                
                xsec->setFluxIntLimits(0.,120.);
                xsec->setUniverses(100);
                //                xsec->setNormalizationType(XSec::kPerNucleon);
                xsec->setNormalizationType(XSec::kSelfNorm);
                //! add cross section
                loop.addXSec(xsec);
                
            } //! end loop over cross sections
        }
    }
    
    //! run
     //! get the output histograms and save them to file
    
    loop.runLoop();
    
   string fname = Form("/pnfs/minerva/persistent/users/${USER}/NukeFiles/NukeCC_xsec_%s_MuonKludged.root" , playlist.c_str());
//     string fname = "/minerva/data/users/${USER}/NukeFiles/NukeCC_xsec_CRAP.root";
    TFile fout(fname.c_str(), "RECREATE");

    for(unsigned int i = 0; i < loop.getXSecs().size(); i++) {
       // loop.getXSecs().at(i)->getXSecHist()->GetDirectory()->cd();
        loop.getXSecs().at(i)->getXSecHist()->Write();
       // loop.getXSecs().at(i)->getEvRateHist()->GetDirectory()->cd();
        loop.getXSecs().at(i)->getEvRateHist()->Write();
        
    }
    //loop.getFluxHist()->GetDirectory()->cd();
    loop.getFluxHist()->Write();
    fout.Write();
    return;
}

int main( int argc, char *argv[] )
{
    cout << "Enter running the GENIE Xsection Extraction for the NukeCC analysis" << endl;
    
//    if( argc == 1 ) {
//        cout << "============================================================================================================" << endl;
//        cout << "    unable to run the xsection extraction, must specify correct arguments" << endl;
//        cout << "      nuclei, plastic ( option: default = plastic )" << endl;
//        cout << "      1track_style, minos_match, all_muons ( option: default = all_muons )" << endl;
//        cout << "      no_cut_proton, exactly_one_proton_450MeV, atleast_one_proton_450MeV ( option: default = no_cut_proton )" << endl;
//        cout << "      run_xsection, run_no_fsi_xsection, run_hN_xsection ( option: run_xsection )" << endl;
//        cout << "    " << endl;
//        cout << "=============================================================================================================" << endl;
//    } else {
//
//        string target_type      = "plastic";
//        string muon_type        = "all_muons";
//        string proton_threshold = "no_cut_proton";
//        string run_type         = "";
    
        
        cout << "running genie xsection extraction" << endl;
        
//        runXSecLooper_NukeCC( target_type );
            runXSecLooper_NukeCCInclusive(argv[1]);
        
        cout << "completed genie xsection extraction" << endl;
        
//    }
    
    cout << "Exit running the GENIE XSection Extraction for the NukeCC analysis" << endl;
    return 0;
}
