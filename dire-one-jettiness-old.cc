//********************************
// Will compute the 1-Jettiness Distribution in DIS
// Right now it computes the xB and Q^2 distributions for DIS
// Usage:
//        make dire-one-jettiness
//        ./dire-one-jettiness dire-one-jettiness.cmnd
//************************


#include "Dire/Dire.h"  //inlcude package that allows us to run DIS

#include "Pythia8/Pythia.h" //inlcude the Pythia package
#include "fastjet/PseudoJet.hh"    //include the PseudoJet class of the fastjet package
#include "fastjet/ClusterSequence.hh" //include the Cluster sequence class of the fastjet package

using namespace std;
using namespace Pythia8;
using namespace fastjet;


// Include analysis. You can replace this include with another file
// to include your own analysis.
#include "DirePlugins/analyses/AnalysisDummy.h"
//#include "DirePlugins/analyses/AnalysisDIS.h"
//==========================================================================

int main( int argc, char* argv[] ){

  // Check that correct number of command-line arguments
  if (argc < 2) {
    cerr << " Unexpected number of command-line arguments ("<<argc-1<<"). \n"
         << " You are expected to provide the arguments." << endl
         << argc-1 << " arguments provided:";
         for ( int i=1; i<argc; ++i) cerr << " " << argv[i];
         cerr << "\n Program stopped. " << endl;
    return 1;
  }
// *************************************************
    
  //Output file name
  
    const char *path="../../../../Pythia-Output-Data-Files/one-jettiness/";
    string filename="";
    char numstr[21];
    string outname = "dire20_output_eCM_90";
    ofstream outfile;
    
    vector<PseudoJet> fjInputs;     //define vector pseudojet object that will, for each event, store all final particle momenta
    
    // create a jet definition:
    // a jet algorithm with a given radius parameter
    double R = 1.0;
    fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, R);
    
    cout << jet_def.description() << endl;
    
    
  Pythia pythia;

  // Create and initialize DIRE shower plugin.
  Dire dire;
  dire.init(pythia, argv[1]);

  // Histogram the weight.
  //Hist histWT("weight",100000,-5000.,5000.);

  // Initialize Pythia analysis.
//  MyAnalysis analysis;
 // analysis.init();


  cout << endl << endl << endl;
  cout << "Start generating events" << endl;

  int nEvent = pythia.settings.mode("Main:numberOfEvents");

/*  double wmax =-1e15;
  double wmin = 1e15;
  double sumwt = 0.;
  double sumwtsq = 0.; */

    
    double eCM=pythia.info.eCM();
    cout << "The center of mass energy is: " << eCM << endl;
    
    filename=path+outname+".dat";
    outfile.open(filename);
    
  // Start event generation loop
  for( int iEvent=0; iEvent<nEvent; ++iEvent ){

      pythia.next();
      int eIid=-1;
      
      double electronIp[4]={eCM/2,0,0,eCM/2};
      double protonIp[4]={eCM/2,0,0,-eCM/2};
      
      double electronFp[4]={0,0,0,0};
      double q[4]={0,0,0,0};
      double Q2=0;
      double xB=0;   //Bjorken-x
      double Pdotq=0;
      
      //for each event clear Pseudojet of previous event data before reading particles for current event
      fjInputs.resize(0);
      
      for(int i=0; i < pythia.event.size(); ++i){        // loop over current event record
        
            if (pythia.event[i].id()==11 && pythia.event[i].isFinal()){   //final electron momentum
                electronFp[0]=pythia.event[i].e();
                electronFp[1]=pythia.event[i].px();
                electronFp[2]=pythia.event[i].py();
                electronFp[3]=pythia.event[i].pz();
            }
          
           //store momenta of all final particles (except the final electron) and push to fjinputs
            if (pythia.event[i].id()!=11 && pythia.event[i].isFinal()){
                PseudoJet particle(
                             pythia.event[i].px(),
                             pythia.event[i].py(),
                             pythia.event[i].pz(),
                             pythia.event[i].e()
                             );
          //Push back particle with its four momentum and other identifying info to fjinputs, the Pseudojet object which stores info about all particles in the current event
                fjInputs.push_back(particle);
            }
      } //end particle loop

      //fastjet analysis for current event
      fastjet::ClusterSequence clust_seq(fjInputs, jet_def); //cluster particles of current event,contained in fjinputs, into jets
      double ptmin = 1.0;  //jet veto GeV
      vector<fastjet::PseudoJet> inclusive_jets = sorted_by_pt(clust_seq.inclusive_jets(ptmin));
       if (inclusive_jets.size() == 0) continue; //if no jets pass the veto, continue to the next event
      
  //    cout << " " << endl;
  //    cout << "Leading Jet pT: " << inclusive_jets[0].pt() << endl;
  //    cout << "Leading Jet pseudorapidity: " << inclusive_jets[0].eta() << endl;
  //    cout << " " << endl;
      
      
      for(int k=0; k < 4; k++){
          q[k]=electronFp[k]-electronIp[k];
      }
      
      Pdotq=protonIp[0]*electronIp[0]
           -protonIp[1]*electronIp[1]
           -protonIp[2]*electronIp[2]
           -protonIp[3]*electronIp[3];
      
      Q2=pow(electronFp[1],2)
         +pow(electronFp[2],2)
         +pow(electronFp[3]-electronIp[3],2)
        -pow(electronFp[0]-electronIp[0],2);
      
      xB=Q2/2/Pdotq;
      
      if (Q2 > 1.0 && Q2 < 8.0){
         outfile << Q2 <<  "   " << xB << endl;
      };
      
  }
    outfile.close();
  return 0;


}
