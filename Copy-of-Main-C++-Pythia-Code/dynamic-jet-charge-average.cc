//----------------------------------------------------------------------
//   Program that calculates the "Standard" and "Dynamic" Jet Charges
//
//    Designed to help produce the plots in 1509.05190 (ATLAS Jet Charge Paper)
//    Computes and writes out Standard and Dynamic Jet Charges the for Leading and Subleading jets
//    Classifies the leading and subleading jets as "more forward" or "more central" jets
//    Writes out the jet charges for the "more forward" and "more central" jets
//
// Usage:
//       make dynamic-jet-charge-average
//       ./dynamic-jet-charge-average joblists-dynamic-jet-charge-average/*.dat *.cmnd
//
// Input Files:
//            *.dat: feeds in a list of jobs where each job is characterized by various parameters and settings
//            *.cmnd: Pythia command file that determines process settings
//
//
// Hard process channels:
//           channel 1: dg --> W^- u
//           channel 2: ug --> W^+ d
//
//           channel 3: gg2gg
//           channel 4: qq2qq (no restriction on quark and anti-quark flavors)
//
//           channel 5: HardQCD:all = on
//
// Types of jet charge:
//
//        Standard Jet Charge
//           jetchargedef=0:
//                             QJ = Sum_i [(pTi/pTJ)^kappa Q_i]
//
//        Dynamic Jet Charge
//           jetchargedef=8:
//                             QJ = Sum_i xi^{kappa(xi,xcut,klt,kgt)}*Q_i
//                             where,
//                               xi=pTi/pTJ
//                               kappa(xi,xcut,klt,kgt)= klt for xi < xcut
//                               kappa(xi,xcut,klt,kgt)= kgt for xi >= xcut
//
//----------------------------------------------------------------------


#include "Pythia8/Pythia.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/contrib/SoftDrop.hh"
#include <iostream>
#include <cstdio>
#include <istream>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

using namespace std;
using namespace Pythia8;
using namespace fastjet;


//----------------Declarartions of Classes and Functions------------------------------------------------


//***********
/* Class MyInfo is derived (derived class) from PseudoJet::UserInfo
It passes additional information such as id, charge, etc. of each particle to Pseudojet */

class MyInfo: public PseudoJet::UserInfoBase {
    public:
    
     /*Define constructor MyInfo of MyInfo class and initialize
      the member protected variables to the arguments to read in by an object of the MyInfo class */
    
    MyInfo(int event_id_in, double event_charge_in, int event_chargeType_in, string event_name_in) :
    _event_id(event_id_in),_event_charge(event_charge_in),_event_chargeType(event_chargeType_in),_event_name(event_name_in){} //initialize protected member variables to the input values
    
    //member functions that return values of the protected variables, which are initialized above
    int event_id() const {return _event_id;}
    double event_charge() const {return _event_charge;}
    int event_chargeType() const {return _event_chargeType;}
    string event_name() const {return _event_name;}
    
    
    //protected member variables that are initialized above to the input values taken in by MyInfo
    protected:
    int _event_id;
    double _event_charge;
    int _event_chargeType;
    string _event_name;
};

//**********
/*function that reads in pTbin and HeavyIonswitch and sets the pT bin range [ptmin,ptmax] according to the
ATLAS paper 1509.05190 and set Pythia command file pTHatmin to be at 80% of pTmin*/

void pTbinRange(int pTbin, int HeavyIonswitch, double &ptmin, double &ptmax, double &pTHatMin);

//**********
/*function that reads in parameters from joblist input file and determines the output filename*/
string output_filename(int jobno, int channel, int jetchargedef,double R, int pTbin, int SDswitch,
                       double &z_cut, double &beta, double &kappa, double &xcut, double &klt,
                       double &kgt, int Hadswitch, int MPIswitch, int HeavyIonswitch);

//**********
/*function that reads in parameters from joblist input file and determines the output filename for data on average jet charge vs pT*/
string output_filename2(int jobno, int channel, int jetchargedef,double R, int SDswitch,
                       double &z_cut, double &beta, double &kappa, double &xcut, double &klt,
                       double &kgt, int Hadswitch, int MPIswitch, int HeavyIonswitch);

//**********
/* Define derived class PythiaSettings (derived from Pythia class) to set pTHatmin,channel, Had, MPI in Pythia */
class PythiaSettings: public Pythia{
    public:
        PythiaSettings(Pythia &pythia, double pTHatMin, int channel, int Hadswitch, int MPIswitch);  //pass Pythia object by reference into constructor function
};



//**********
/* Define derived class JetChargeAlg (derived from the PseudoJet class) to calculate the jet charge*/
class JetChargeAlg: public PseudoJet{
    public:
        JetChargeAlg(vector<PseudoJet>  &constituentsIN, double pJTIN, double kappaIN, double xcutIN, double kltIN, double kgtIN); //contructor
        double StandardJetCharge(); //returns standard jet charge
        double DynamicJetCharge(); //returns dynmaic jet charge
    
    private:
        vector<PseudoJet> constituents;  //constiuents  in the Jet
        double pJT,kappa,xcut,klt,kgt;
    
};

//************
/* Define function that computes the "more central (QJC)" and "more forward(QJF)" jet charges from the leading (QJ) and SUBleading (QJ1) jet chares*/
void CentralForwardJetCharges(double &QJC, double &QJF, double QJ, double QJ1, double etaJ, double etaJ1);

//--------------------------------------------------------------------------------------------------------------------











//-------------------- Begin Main Program ----------------------------------------------------------
//-------------------- Begin Main Program ----------------------------------------------------------
//-------------------- Begin Main Program ----------------------------------------------------------
//-------------------- Begin Main Program ----------------------------------------------------------
//-------------------- Begin Main Program ----------------------------------------------------------
//-------------------- Begin Main Program ----------------------------------------------------------
//-------------------- Begin Main Program ----------------------------------------------------------
//-------------------- Begin Main Program ----------------------------------------------------------
//-------------------- Begin Main Program ----------------------------------------------------------


int main(int argc, char* argv[]){
    
    string filename_in_1=argv[1]; //read in joblist *.dat filename from command line
    string filename_in_2=argv[2]; //read in *.cmnd command filename from command line
 
    //Parameters to be read in from filename_in_1
    int jobno; //index number for each job (each line of jet-charge.dat)
    int channel; //hard process channel (see code for channels at the top of the file)
    int jetchargedef; //0 = standard jet charge,  8= dynamic new jet charge
    
    double R; // jet radius
    int pTbin;  //integer value to denotes the pTbin=[pTmin (GeV),pTmax (GeV)];
                // 0=default : [220,270] for pp, [80,150] for Pb-Pb
                //pp-collisions: 0=[220,270], 1=[50,100), 2=[100,200), 3=[200,300), 4=[300,400), 5=[400,500),
                //               6=[500,600), 7=[600,800), 8=[800,1000), 9=[1000,1200), 10=[1200,1500]
    
    int SDswitch;  //SoftDrop Switch: 0=off, 1=on
    double z_cut; //SD grooming parameter
    double beta;  //SD grooming parameter
    
    double kappa; //kappa parameter in standard jet charge definition
    
    //parameters in the "dynamic kappa" profiles
    double xcut;
    double klt;
    double kgt;

    int Hadswitch; //Hadronization Switch: 0=off, 1=on
    int MPIswitch; //MPI Switch: 0=off, 1=on
    int HeavyIonswitch; //Heavy Ion Switch: 0=off (pp), 1=on (Pb-Pb)
    
    
    // set up an ifstream object and open filename_in_1 (joblist file)
    ifstream infile;
    infile.open(filename_in_1);
    
    // set path of output directories and initialize output filename variables
    const char *path="../../Pythia-Output-Data-Files/jet-charge-paper-plots-13-TeV/";
    const char *path2="../../Pythia-Output-Data-Files/jet-charge-average-paper-plots-13-TeV/";
    string outname="";
    string outname2="";
    string filename_out="";    //file will contain list of QJF and QJC values for a given pTbin
    string filename_out2="";   //file will contain columns of the average values QJCavg and QJFavg over each pT bin and the central value of the pTbin
    ofstream outfile, outfile2; // Set up ofstream object to output values

    double ptmin,ptmax,etacut; //parameters for kinematic cuts for leading jet
    
    
    Pythia pythia;    //define pythia object
    pythia.readFile(filename_in_2);   //read in Pythia command file
    int nEvent = pythia.mode("Main:numberOfEvents");   //read in number of events to be generated by Pythia from the *.cmnd file
    
    
    int NoOfpTbins=10; //default value is 10. Check the "pTbinRange" function for the total number of pTbins
    //read in parameters for each job from the joblist *.dat file and generate events and write data to ouput file
    
    bool isOpen=false;
    
    while(infile >> jobno >> channel >> jetchargedef  >> R >> pTbin >> SDswitch >> z_cut >> beta >>  kappa >> xcut >>  klt >> kgt >> Hadswitch >> MPIswitch >> HeavyIonswitch){
        
        
        int dijet=0;  //for channels 1 and 2 (qg-->Wq') we are only interested in the "leading" jet
        if (channel==3 || channel==4 || channel==5){
            dijet=1;            //for channels 3,4,5 we are interested in the "leading" and "SUBleading" jets
        }
        
        
        //call function that returns output filename based on input joblist parameters
        outname=output_filename(jobno, channel, jetchargedef, R, pTbin, SDswitch, z_cut, beta, kappa, xcut, klt, kgt, Hadswitch, MPIswitch, HeavyIonswitch);
        filename_out=path+outname;   //add directory path to output filename
        outfile.open(filename_out); //open output file
        
        
        if (isOpen==false){   //open output file 2 only once as we loop through pTbins 1 through 10
                outname2=output_filename2(jobno, channel, jetchargedef, R, SDswitch, z_cut, beta, kappa, xcut, klt, kgt, Hadswitch, MPIswitch, HeavyIonswitch);
            
                filename_out2=path2+outname2;   //add directory path to output filename
                outfile2.open(filename_out2); //open output file
        }
        isOpen=true; //set to false so that filename_out2 is only opened once
        
        
        //set rapidity cuts for jets
        if (HeavyIonswitch==0){
            etacut= 2.1;   // pp collisions
        }
        else{
            etacut= 0.9;   // Pb-Pb collisions
        }
        
        
        //set pT bin ranges [ptmin,ptmax] according to the ATLAS paper 1509.05190 and set Pythia phasespace pTHatmin to be at 80% of pTmin
        double pTHatMin;
        pTbinRange(pTbin, HeavyIonswitch, ptmin, ptmax, pTHatMin);  //read in pTbin and HeavyIonswitch and pass ptmin,ptmax,pTHatMin by reference to set their values
      
        PythiaSettings pythiasettings(pythia, pTHatMin, channel, Hadswitch, MPIswitch); //creat class object to activate constructor to set pTHatmin,channel, Had, MPI in Pythia
        pythia.init(); // Initialize pythia
        
        double QJsum=0.0;  //reset (initialize) sum over "leading" jet charge in current pTbin
        double QJFsum=0.0; //reset (initialize) sum over "more forward" jet charge in current pTbin
        double QJCsum=0.0; //reset (initialize) sum over "more central" jet charge in current pTbin
        
        double QJavg=0.0; //reset (initialize) average value of "leading" jet charges in current pT bin
        double QJFavg=0.0; //reset (initialize) average value of "more forward" jet charges in current pT bin
        double QJCavg=0.0; //reset (initialize) average value of "more central" jet charges in current pT bin
        
        int NoOfEvents=0; // reset (initialize) the counter for the number of events that pass phase cuts and will be recorded to output file
    
        
       /*define vector pseudojet object that will, for each event, store all final particle momenta
        and other information through the MyInfo class and be inputed to fastjet which will run a jet algorithm*/
        vector<PseudoJet> fjInputs;
        
        fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, R); //choose the jet algorithm
        contrib::SoftDrop sd(beta,z_cut,R); // create softdrop object sd. Soft drop condition: z > z_cut (theta_ij/R)^beta
    
        
        int id; //temporary variable that will store particle id inside the particle loop
        
        // begin loop over pythia events
        //********************************** begin event loop *******************************
        for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
            if (!pythia.next()) continue;  //generate next event or skip event if false

            // select dg --> W^- u channel-1 from the qg2Wq setting
            if (channel==1){
                  if ((pythia.info.id1()!=1 && pythia.info.id1()!=21  ) || (pythia.info.id2()!=1 && pythia.info.id2()!=21) ){
                    continue;
                }
            }
            
            // select ug --> W^+ d channel-2 from the qg2Wq setting
            if (channel==2){
                if ((pythia.info.id1()!=2 && pythia.info.id1()!=21  ) || (pythia.info.id2()!=2 && pythia.info.id2()!=21) ){
                    continue;
                }
            }

        
            //for each event clear Pseudojet object fjInputs of previous event data before reading particles for current event
            fjInputs.resize(0);
        
            //begin loop over particles in each event
            //********************************** begin particle loop *******************************
            for (int i = 0; i < pythia.event.size(); ++i){
                
               if (pythia.event[i].isFinal()){
                 
                   id = pythia.event[i].id();  //store final state particle id
                   
                   //do not send leptons, electroweak bosons, and higgses into jet algorithm (fastjet)
                   if (id==11||id==-11||id==12||id==-12||id==13||id==-13||id==14||id==-14||id==16||id==-16||id==23||id==24||id==-24||id==25){
                       continue;
                   }
                   
                   
                   PseudoJet particle(                              //store final particle momentum in Pseudojet object "particle"
                                   pythia.event[i].px(),
                                   pythia.event[i].py(),
                                   pythia.event[i].pz(),
                                   pythia.event[i].e()
                                   );
                
                   
                   particle.set_user_info(new MyInfo(               //Pass other info about the particle to the "particle" object through MyInfo, the derived class of PseudoJet
                                                  pythia.event[i].id(),
                                                  pythia.event[i].charge(),
                                                  pythia.event[i].chargeType(),
                                                  pythia.event[i].name()
                                                  )
                                          );
                   
                   particle.set_user_index(i);   //store particle loop index for the current particle
                   
                   fjInputs.push_back(particle);  //Push back "particle" with its four momentum and other identifying info to fjinputs
                   
              }                                   //end if statement for "pythia.event[i].isFinal()"
            }                                     //end particle loop
        //********************************** end particle loop *******************************
            
            
        if (fjInputs.size() == 0) continue; //if no "particles" found in event, continue to the next event
            
        
        //fastjet analysis for current event
        fastjet::ClusterSequence clust_seq(fjInputs, jet_def); //cluster particles of current event,contained in fjinputs, into jets
        
        //sort list of jets by pT and veto jets for which pT < pTveto and store in vector PseudoJet object
        double pTveto=ptmin;
        if (dijet==0){
            pTveto=0.0;  //if only interested in leading get, remove pTveto on
        }
        vector<fastjet::PseudoJet> inclusive_jets = sorted_by_pt(clust_seq.inclusive_jets(pTveto));
            
            
      
        
        if (inclusive_jets.size() == 0) continue; //if no jets pass the veto, continue to the next event
        if (inclusive_jets[0].perp() > ptmax) continue; //set leading maximum allowed pT for leading jet
        if (abs(inclusive_jets[0].eta()) > etacut ) continue; //set rapidity cut on leading jet
        if (inclusive_jets.size() < 2) continue; //require at least two jets (with pTveto=0 for dijet=0 qg->Wq' channels)

        
        //for dijet channels, set cuts on leading and SUBleading jets, see 1509.05190 (ATLAS Jet Charge Paper)
        if(dijet==1){
                if (inclusive_jets[1].perp() > ptmax) continue; //set SUBleading maximum allowed pT for leading jet
                if (abs(inclusive_jets[1].eta()) > etacut ) continue; //set rapidity cut on SUBleading jet
                if(inclusive_jets[0].perp()/inclusive_jets[1].perp() > 1.5) continue;
        }
        

            
        //apply soft drop grooming to leading and SUBleading jets
        
        PseudoJet sd_jet=sd(inclusive_jets[0]);
        assert(sd_jet != 0); //soft drop is a groomer (not a tagger), it should always return a soft-dropped jet
        PseudoJet sd_jet1=sd(inclusive_jets[1]);
        assert(sd_jet1 != 0); //soft drop is a groomer (not a tagger), it should always return a soft-dropped jet
        
            
        
   
        //Pseudojet objects that can store the constituents of a jet
        vector<fastjet::PseudoJet> constituents, constituents1;
        double pJT,pJT1;//jet transverse momentum
            
        if (SDswitch==0){  //soft drop off
            //leading jet
            pJT=inclusive_jets[0].perp();
            constituents = inclusive_jets[0].constituents();

            //SUBleading jet
            pJT1=inclusive_jets[1].perp();
            constituents1 = inclusive_jets[1].constituents();
        }
        else{               //soft drop on
            //leading jet
            pJT=sd_jet.perp();
            constituents = sd_jet.constituents();
 
            //SUBleading jet
            pJT1=sd_jet1.perp();
            constituents1 = sd_jet1.constituents();
        }
            
        //Initialize Jet Charges
        double QJ=0.0; //leading jet
        double QJ1=0.0; //SUBleading jet
            
        //Define objects for the JetChargeAlg class for the leaing and SUBleading jets
        JetChargeAlg QJL(constituents, pJT, kappa, xcut, klt, kgt), QJSL(constituents1, pJT1, kappa, xcut, klt, kgt);
     
            
        //calculate standard jet charge or dynamic groomed standard jet charge
        if (jetchargedef==0){
            //Standard Jet Charge
            QJ=QJL.StandardJetCharge();          //leading jet
            QJ1=QJSL.StandardJetCharge();        //SUBLEADING jet
         }
        //jetchargedef=8 (xi-dependent dynamic standard jet charge )
        else if (jetchargedef==8){
            //Dynamic Jet Charge
            QJ=QJL.DynamicJetCharge();            //leading jet
            QJ1=QJSL.DynamicJetCharge();          //SUBleading jet
        }
            
           
        //Call function that computes the more central (QJC) and more forward (QJF) jet charges
        double QJC=0.0,QJF=0.0; //initialize
        CentralForwardJetCharges(QJC, QJF, QJ, QJ1, inclusive_jets[0].eta(), inclusive_jets[1].eta());

        QJsum = QJsum + QJ;             //update QJsum by adding the charge of the current leading jet
        QJFsum = QJFsum + QJF;          //update QJFsum by adding the charge of the current forward jet
        QJCsum = QJCsum + QJC;           //update QJCsum by adding the charge of the current central jet
            
        NoOfEvents = NoOfEvents+1;      //update number of jets that pass cuts and will be recorded to output file
            
        //write out jet charges to output file
        if (dijet==1){
            outfile << setprecision(10) << QJC << " " << QJF <<  endl;
        }
        else{
             outfile << setprecision(10) << QJ << endl;
        }
            

    }//end event loop
        
  
        //write to a file the QJCavg, QJFavg, and central pTbin values
        if (dijet==1){
            QJFavg=QJFsum/(double)NoOfEvents;   //compute average forward jet charge average
            QJCavg=QJCsum/(double)NoOfEvents;   //compute average central jet charge average
            outfile2 << setprecision(10) << QJCavg << " " << QJFavg <<  " " << (ptmin+ptmax)/2.0 << endl;
        }
        else{
            QJavg= QJsum/(double)NoOfEvents;    //compute average leading jet charge average
            outfile2 << setprecision(10) << QJavg << " " << (ptmin+ptmax)/2.0 << endl;
        }
        
        
        outfile.close();    //close output file

        if (pTbin==NoOfpTbins){
            outfile2.close();   //close output file2 after going over all pTbins 1 through 10
        }
        
  } // end the jobs loop
    

  return 0;
}



//-------------------- End Main Program ----------------------------------------------------------
//-------------------- End Main Program ----------------------------------------------------------
//-------------------- End Main Program ----------------------------------------------------------
//-------------------- End Main Program ----------------------------------------------------------
//-------------------- End Main Program ----------------------------------------------------------
//-------------------- End Main Program ----------------------------------------------------------
//-------------------- End Main Program ----------------------------------------------------------
//-------------------- End Main Program ----------------------------------------------------------
//-------------------- End Main Program ----------------------------------------------------------














//*************   Function Definitions *******************************




//**********************Function 1******************

void pTbinRange(int pTbin, int HeavyIonswitch, double &ptmin, double &ptmax, double &pTHatMin){
    
    
    if (HeavyIonswitch==0){    //pp collisons

        
        if(pTbin==0){  //default pTbin
            ptmin = 220.0;  //jet pT-veto in GeV
            ptmax = 270.0;  //max pT of leading jet in GeV
        }
        if(pTbin==1){  //default pTbin
            ptmin = 50.0;  //jet pT-veto in GeV
            ptmax = 99.9;  //max pT of leading jet in GeV
        }
        if(pTbin==2){  //default pTbin
            ptmin = 100.0;  //jet pT-veto in GeV
            ptmax = 199.9;  //max pT of leading jet in GeV
        }
        if(pTbin==3){  //default pTbin
            ptmin = 200.0;  //jet pT-veto in GeV
            ptmax = 299.9;  //max pT of leading jet in GeV
        }
        if(pTbin==4){  //default pTbin
            ptmin = 300.0;  //jet pT-veto in GeV
            ptmax = 399.9;  //max pT of leading jet in GeV
        }
        if(pTbin==5){  //default pTbin
            ptmin = 400.0;  //jet pT-veto in GeV
            ptmax = 499.9;  //max pT of leading jet in GeV
        }
        if(pTbin==6){  //default pTbin
            ptmin = 500.0;  //jet pT-veto in GeV
            ptmax = 599.9;  //max pT of leading jet in GeV
        }
        if(pTbin==7){  //default pTbin
            ptmin = 600.0;  //jet pT-veto in GeV
            ptmax = 799.9;  //max pT of leading jet in GeV
        }
        if(pTbin==8){  //default pTbin
            ptmin = 800.0;  //jet pT-veto in GeV
            ptmax = 999.9;  //max pT of leading jet in GeV
        }
        if(pTbin==9){  //default pTbin
            ptmin = 1000.0;  //jet pT-veto in GeV
            ptmax = 1199.9;  //max pT of leading jet in GeV
        }
        if(pTbin==10){  //default pTbin
            ptmin = 1200.0;  //jet pT-veto in GeV
            ptmax = 1500.0;  //max pT of leading jet in GeV
        }
        
    }
    else{        //Ion Pb-Pb collisions
        if(pTbin==0){       //default pTbin
            ptmin = 80.0;  //jet pT-veto in GeV
            ptmax = 150.0;  //max pT of leading jet in GeV
        }
    }
    
    // minimal pT scale in hard process in Pythia (80% of ptmin)
    pTHatMin=80.0/100.0*ptmin;
    
}



//*******************************Function 2*********************

string output_filename(int jobno, int channel, int jetchargedef,double R, int pTbin, int SDswitch,
                       double &z_cut, double &beta, double &kappa, double &xcut, double &klt,
                       double &kgt, int Hadswitch, int MPIswitch, int HeavyIonswitch){
    
    
   
    string outname=""; //string that will store the filename that will be returned by this function
   // string filename_out=""; //string that will store the path+filename that will be returned by this function
    
    //Depending on jetchardedef option, set unused parameters to zero to clean up output file name
    if (jetchargedef ==0){   //if jet charge is not "dynamic" set kappa-profile parameters to zero, which will be reflected in the output file name
        xcut=0;
        kgt=0;
        klt=0;
    }
    
    if (jetchargedef ==8){ //set kappa=0 in file name for "dynamic" jet charge and instead use the "dynamic kappa"
        kappa=0;
    }
    
    if (SDswitch ==0){   //if soft drop is off, set gromming parameters to zero, which will be reflected in the file name
        z_cut=0;
        beta=0;
    }
    
    
    
    
    
    
    //convert input parameter values to strings to construct output file name
    string jobnostr=to_string(jobno);
    string channelstr=to_string(channel);
    string jetchargedefstr=to_string(jetchargedef);
    
    string Rx10str=to_string((int)(R*10));
    string pTbinstr=to_string((int)(pTbin));
    
    
    string SDswitchstr=to_string((int)SDswitch);
    string z_cutx100str=to_string((int)(z_cut*100));
    string betax10str=to_string((int)(beta*10));
    
    string kappax100str=to_string((int)(kappa*100));
    
    string xcutx100str=to_string((int)(xcut*100));
    string kgtx100str=to_string((int)(kgt*100));
    string kltx100str=to_string((int)(klt*100));
    
    string Hadswitchstr=to_string((int)Hadswitch);
    string MPIswitchstr=to_string((int)MPIswitch);
    
    string HeavyIonswitchstr=to_string((int)HeavyIonswitch);
    
    
    //construct output filename and path
    outname="job-"+jobnostr+"-channel-"+channelstr+"-charge-"+jetchargedefstr
    
    +"-Rx10-"+Rx10str
    +"-pTbin-"+pTbinstr
    
    +"-SD-"+SDswitchstr
    +"-zcutx100-"+z_cutx100str
    +"-betax10-"+betax10str
    
    +"-kappax100-"+kappax100str
    
    +"-xcutx100-"+xcutx100str
    +"-kgtx100-"+kgtx100str
    +"-kltx100-"+kltx100str
    
    +"-Had-"+Hadswitchstr
    +"-MPI-"+MPIswitchstr
    
    +"-HeavyIon-"+HeavyIonswitchstr;
    
    outname = outname+".dat";

    
    return outname;
    
}




//*******************************Function 2.5 *********************

string output_filename2(int jobno, int channel, int jetchargedef,double R, int SDswitch,
                       double &z_cut, double &beta, double &kappa, double &xcut, double &klt,
                       double &kgt, int Hadswitch, int MPIswitch, int HeavyIonswitch){
    
    
    
    string outname2=""; //string that will store the filename that will be returned by this function
    // string filename_out=""; //string that will store the path+filename that will be returned by this function
    
    //Depending on jetchardedef option, set unused parameters to zero to clean up output file name
    if (jetchargedef ==0){   //if jet charge is not "dynamic" set kappa-profile parameters to zero, which will be reflected in the output file name
        xcut=0;
        kgt=0;
        klt=0;
    }
    
    if (jetchargedef ==8){ //set kappa=0 in file name for "dynamic" jet charge and instead use the "dynamic kappa"
        kappa=0;
    }
    
    if (SDswitch ==0){   //if soft drop is off, set gromming parameters to zero, which will be reflected in the file name
        z_cut=0;
        beta=0;
    }
    
    
    
    
    
    
    //convert input parameter values to strings to construct output file name
    string jobnostr=to_string(jobno);
    string channelstr=to_string(channel);
    string jetchargedefstr=to_string(jetchargedef);
    
    string Rx10str=to_string((int)(R*10));

    
    
    string SDswitchstr=to_string((int)SDswitch);
    string z_cutx100str=to_string((int)(z_cut*100));
    string betax10str=to_string((int)(beta*10));
    
    string kappax100str=to_string((int)(kappa*100));
    
    string xcutx100str=to_string((int)(xcut*100));
    string kgtx100str=to_string((int)(kgt*100));
    string kltx100str=to_string((int)(klt*100));
    
    string Hadswitchstr=to_string((int)Hadswitch);
    string MPIswitchstr=to_string((int)MPIswitch);
    
    string HeavyIonswitchstr=to_string((int)HeavyIonswitch);
    
    
    //construct output filename and path
    outname2="job-"+jobnostr+"-channel-"+channelstr+"-charge-"+jetchargedefstr
    
    +"-Rx10-"+Rx10str
    
    +"-SD-"+SDswitchstr
    +"-zcutx100-"+z_cutx100str
    +"-betax10-"+betax10str
    
    +"-kappax100-"+kappax100str
    
    +"-xcutx100-"+xcutx100str
    +"-kgtx100-"+kgtx100str
    +"-kltx100-"+kltx100str
    
    +"-Had-"+Hadswitchstr
    +"-MPI-"+MPIswitchstr
    
    +"-HeavyIon-"+HeavyIonswitchstr;
    
    outname2 = outname2+".dat";
    
    
    return outname2;
    
}




//*******************************Function 3*********************

//Constructor function of the derived class PythiaSettings (derived from Pythia class)

PythiaSettings::PythiaSettings(Pythia &pythia, double pTHatMin, int channel, int Hadswitch, int MPIswitch){
    

    
    //set minimal pT scale in hard process in Pythia (80% of ptmin)
    string pTHatMinstr=to_string((double)pTHatMin);
    string pythiapTHatMinstr="PhaseSpace:pTHatMin ="+pTHatMinstr;
    pythia.readString(pythiapTHatMinstr);
    
    //hard process selection
    if (channel==1 || channel==2){
        pythia.readString("WeakBosonAndParton:qg2Wq = on");
    }
    
    if (channel==3){
        pythia.readString("HardQCD:gg2gg = on");
    }
    
    if (channel==4){
        pythia.readString("HardQCD:qq2qq = on");
    }
    
    if (channel==5){
        pythia.readString("HardQCD:all = on");
        // pythia.readString("HardQCD:gg2gg = on");
        // pythia.readString("HardQCD:qq2qq = on");
    }
    
    
    //Turn hadronization on or off based on settings specified in filename_in
    if (Hadswitch==1) {
        pythia.readString("HadronLevel:Hadronize=on");
    }
    else{
        pythia.readString("HadronLevel:Hadronize=off");
    }
    
    //Turn MPI on or off based on settings specified in filename_in
    if (MPIswitch==1) {
        pythia.readString("PartonLevel:MPI=on");
    }
    else{
        pythia.readString("PartonLevel:MPI=off");
    }
    
}


//*******************************Function 4*********************


//Constructor for the JetChargeAlg class
JetChargeAlg::JetChargeAlg(vector<PseudoJet>  &constituentsIN, double pJTIN, double kappaIN, double xcutIN, double kltIN, double kgtIN){
    
    constituents=constituentsIN;
    pJT=pJTIN;
    kappa=kappaIN;
    xcut=xcutIN;
    klt=kltIN;
    kgt=kgtIN;
}


//Standard Jet Charge function member of the JetChargeAlg class
double JetChargeAlg::StandardJetCharge(){
    //calculate standard jet charge or dynamic groomed standard jet charge

    double QJet=0.0; //initialize jet charge to zero (clear any previously stored value)

        //jet charge calculation
        for (unsigned int k=0; k < constituents.size(); k++){
            QJet=QJet+constituents[k].user_info<MyInfo>().event_charge()*pow(constituents[k].pt(),kappa);
        }
        QJet=QJet/pow(pJT,kappa); //divide by the leading jet pT^kappa
       return QJet;
}


//Dynamic Jet Charge function member of the JetChargeAlg class
double JetChargeAlg::DynamicJetCharge(){
    //calculate standard jet charge or dynamic groomed standard jet charge
    
    double QJet=0.0; //initialize jet charge to zero (clear any previously stored value)
    
    double xk=0.0;  //initialize pT fraction for leading and subleading jet
    double kappaeff=0.0; //initialize variables  for leading and subleading jet
    
    //dynamic jet charge calculation
    for (unsigned int k=0; k < constituents.size(); k++){
        xk=constituents[k].pt()/pJT;  //define the xk=pTk/pJT variable
        
        if (xk < xcut){
            kappaeff=klt;
        }
        else{
            kappaeff=kgt;
        }
        
        QJet=QJet+constituents[k].user_info<MyInfo>().event_charge()*pow(xk,kappaeff);
    
    }
    return QJet;
}



//*****************************************Function ****************************
void CentralForwardJetCharges(double &QJC, double &QJF, double QJ, double QJ1, double etaJ, double etaJ1){
    
    double abs_etaJ=abs(etaJ);
    double abs_etaJ1=abs(etaJ1);
    
    if(abs_etaJ >=abs_etaJ1){
        QJF=QJ;
        QJC=QJ1;
    }
    else{
        QJF=QJ1;
        QJC=QJ;
    }
}
