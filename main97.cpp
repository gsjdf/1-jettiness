//////////////////////////////////////////////////////////////////////////////
//
// File to simulate DIS process e+p->e'+jet+X 				  ////
//
/////////////////////////////////////////////////////////////////////////////
#include "Pythia8/Pythia.h"
#include <cmath>
#include "fastjet/ClusterSequence.hh"
using namespace Pythia8;
using namespace std;
using namespace fastjet;

///////////some short hands///////////
double dot(Vec4 a,Vec4 b){
    return a.e()*b.e()-a.px()*b.px()-a.py()*b.py()-a.pz()*b.pz();
}

///////////end some short hands///////////

///////////constants///////////
// range of tau1
double maxx = 1.005;
double minx = 0.005;
int bin=100;

double pTjetMin=0.0; // lower limit of ptjet
double radius = 1.0; // jet cone radius
//range of Q2 and y
double Q2max=80.0;
double Q2min=60.0;
double ymax=0.6;
double ymin=0.2;
///////////end constants///////////
int main() {
    //Storage Path
    ofstream coutx("/home/matata/CLionProjects/jettiness/tau1ah.txt");
    ofstream couty("/home/matata/CLionProjects/jettiness/sigma1ah.txt");
    ofstream coutyerr("/home/matata/CLionProjects/jettiness/sigma1aerrh.txt");
    ofstream couttot("/home/matata/CLionProjects/jettiness/tot1ah.txt");

    Pythia pythia;
    Event &event = pythia.event;
    Settings &settings = pythia.settings;
    ParticleData &particleData = pythia.particleData;

    // Common settings for all the subruns
    pythia.readFile("/home/matata/CLionProjects/jettiness/main97.cmnd");

    // Number of events
    int nEvent = 1.e8;//1.e7;
    settings.parm("PhaseSpace:Q2Min", Q2min);
    vector <PseudoJet> fjInputs; /* Fastjet input */
    // create a jet definition:
    // a jet algorithm with a given radius parameter
    //----------------------------------------------------------
    JetDefinition jet_def(antikt_algorithm, radius);
    pythia.init();
    // Begin event loop. Generate event, Skip if error. list first few.


    //data for x and y, here x is tau1a and y is sigma
    double xbin[bin], ybin[bin];
    for (int i = 0; i < bin; ++i) {
        xbin[i] = (maxx - minx) / bin * i + minx;
        ybin[i] = 0;
    }

    for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
        if (!pythia.next()) continue;

        // fastjet analysis
        // create a fastjet::PseudoJet with these components and put it onto
        // back of the input_particles vector
        fjInputs.resize(0);
        for (int i = 0; i < event.size(); ++i) if (event[i].isFinal() && event[i].id()!=11 && event[i].mother1()!=6) {
                PseudoJet particle( event[i].px(), event[i].py(),
                                    event[i].pz(), event[i].e() );
                fjInputs.push_back(particle);
            }
        // run the jet clustering with the above jet definition
        //----------------------------------------------------------
        ClusterSequence clust_seq(fjInputs, jet_def);
        vector <PseudoJet> inclusiveJets, sortedJets;
        // get the resulting jets ordered in pt
        //----------------------------------------------------------
        inclusiveJets = clust_seq.inclusive_jets(pTjetMin);
        sortedJets = sorted_by_pt(inclusiveJets);

        ///////////////////////////////calculate tau
        Vec4 pProton = event[1].p();
        Vec4 peIn = event[4].p();
        Vec4 peOut = event[6].p();
        Vec4 pPhoton = peIn - peOut;
        double Q2=-dot(pPhoton, pPhoton);
        double y = dot(pProton, pPhoton)/ dot(pProton, peIn);
        if (Q2<=Q2min||Q2>=Q2max||y<=ymin||y>=ymax){continue;}
        double xB = -dot(pPhoton, pPhoton) / 2.0 / dot(pProton, pPhoton);
        double qB[4]={xB*pProton.e(),0.0,0.0,xB*pProton.pz()};
        double qJ[4]={sortedJets[0].Et()*cosh(sortedJets[0].rap()),sortedJets[0].px(),sortedJets[0].py(),sortedJets[0].Et()*sinh(sortedJets[0].rap())};
        double tau1a=0;
        for (int i = 0; i < event.size(); ++i) if (event[i].isFinal() && event[i].id()!=11 && event[i].mother1()!=6) {
            tau1a+=min(event[i].e()*qB[0]-event[i].px()*qB[1]-event[i].py()*qB[2]-event[i].pz()*qB[3],event[i].e()*qJ[0]-event[i].px()*qJ[1]-event[i].py()*qJ[2]-event[i].pz()*qJ[3]);
        }
        tau1a=tau1a*2/Q2;
        if (minx < tau1a && tau1a < maxx) {
            int j = int((tau1a - minx) / (maxx - minx) * bin);
            ybin[j] = ybin[j] + 1;
        }
    }
    // end of event loop. statistics
    pythia.stat();
    float sigma = pythia.info.sigmaGen() * 1.0e12; // total cross section in fb
    float norm = 1.0 * bin/(maxx - minx)*sigma/pythia.info.weightSum();

    //////save data
    cout<<pythia.info.weightSum()<<endl;
    for (int i = 0; i < bin; ++i) {
        if (i == 0) {
            coutx << xbin[i];
            couty << ybin[i]*norm;
            coutyerr << 1/sqrt(ybin[i]);
        } else { coutx << "," << xbin[i];
            couty << "," << ybin[i]*norm;
            coutyerr << "," << 1/sqrt(ybin[i]);}
    }
    couttot<<sigma;
    return 0;
}
