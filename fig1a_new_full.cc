//------ DON'T TOUCH THIS PART! ------
#include <bits/dis-phasespace.h>
#include <bits/dis-process.h>
#include <bits/dis-jetfunc.h>
#include <nlo++-module_add.h>

//----- used namespaces -----
using namespace nlo;
using namespace std;


//----- declaration of the user defined functons -----
void inputfunc(unsigned int&, unsigned int&, unsigned int&);
void psinput(phasespace_dis *, double&, double&, double&, double&, double&, double&, double&, double&);
user_base_dis * userfunc();

typedef unsigned long int (*module_add_type)(bool, const list<basic_string<char> >&, const basic_string<char>&);
extern  module_add_type module_add;


//----- array of the symbols symbols -----
extern "C"{

struct {
    const char *name;
    void *address;
} user_defined_functions[] =
        {
                //   process index: hhc for e + p --> jets (DIS)
                {"procindex", (void *) "dis"},

                //   input function
                {"inputfunc", (void *) inputfunc},

                //   phase space input function
                {"psinput", (void *) psinput},

                //   user defined functions
                {"userfunc",  (void *) userfunc},

                //   module to generate the readable result
                {"main_module_add", (void *) module_add},

                //  end of the list
                {0, 0}
        };
}
//------ END OF THE DO-NOT-TOUCH-PART ------


//------ USER defined part starts here ------
//#include "pdf-cteq6.h"
#include "lha_pdf_full.h"

class UserDIS : public user1h_dis
{
public:
    //   init and user function
    void initfunc(unsigned int);
    void userfunc(const event_dis&, const amplitude_dis&);

private:
    //   pdf
    //pdf_cteq6 pdf;
    lha_pdf pdf;
};

//----- defines the module to sum up the results of the different runs -----
module_add_type module_add =  main_module_add<basic_user_result<user1h_dis::distbook_type> >;


user_base_dis * userfunc() {
    return new UserDIS;
}


void inputfunc(unsigned int& nj, unsigned int& nu, unsigned int& nd)
{
    //  number of jets
    nj = 2U;

    //  number of the up and down type flavours
    nu = 2U;
    nd = 3U;
}
//constants
double dtau1amax=1.005;
double dtau1amin=0.005;
int nbin=100;
double pTjetMin=0.;
double radius = 1.; // jet cone radius


void psinput(phasespace_dis *ps, double& el, double& eh, double& q2min, double& q2max,
             double& xmin, double& xmax, double& ymin, double& ymax)
{
    //   we use the default phase space generator
    ps = 0;

    el = 159.5;    // GeV
    eh = 159.5;     // GeV

    //  Q^2 cuts
    q2min = 60.0;     // GeV^2
    q2max = 80.0;  // GeV^2

    //   xB cuts
    xmin = 0.0;
    xmax = 1.0;

    //   y cuts
    ymin = 0.2;
    ymax = 0.6;
}

void UserDIS::initfunc(unsigned int njet)
{
    phys(1, "tau1ac", nbin, dtau1amin, dtau1amax);
    phys(2, "tau1au", nbin, dtau1amin, dtau1amax);
    phys(3, "tau1al", nbin, dtau1amin, dtau1amax);
    phys(4, "number", nbin, dtau1amin, dtau1amax);
}

extern"C" double xalfaem_(double *);

double xalpha_em(double mq2) {
    return xalfaem_(&mq2);
}

#include "fastjet/ClusterSequence.hh"
using namespace fastjet;
void UserDIS::userfunc(const event_dis& p, const amplitude_dis& amp)
{


    //----- weight of the event -----

    // fastjet analysis
    // create a fastjet::PseudoJet with these components and put it onto
    // back of the input_particles vector
    vector <PseudoJet> fjInputs; /* Fastjet input */
    fjInputs.resize(0);
    for (unsigned int i=1;i<=p.upper();i++)  {
        PseudoJet particle( p[i].X(), p[i].Y(),
                            p[i].Z(), p[i].T() );
        fjInputs.push_back(particle);
    }
    // create a jet definition:
    // a jet algorithm with a given radius parameter
    //----------------------------------------------------------
    JetDefinition jet_def(antikt_algorithm, radius);
    // run the jet clustering with the above jet definition
    //----------------------------------------------------------
    ClusterSequence clust_seq(fjInputs, jet_def);
    vector <PseudoJet> inclusiveJets, sortedJets;
    // get the resulting jets ordered in pt
    //----------------------------------------------------------
    inclusiveJets = clust_seq.inclusive_jets(pTjetMin);
    sortedJets = sorted_by_pt(inclusiveJets);

    double Q2 = -((p[-1] - p[-2]).mag2());
    double xB=Q2/2.0/dot(p[hadron(0)],(p[-1] - p[-2]));
    
    double hardscale=Q2;
    double hardscale1=hardscale/4;
    double hardscale2=hardscale*4;
    
    double qB[4]={xB*p[hadron(0)].T(),0.0,0.0,xB*p[hadron(0)].Z()};
    double qJ[4]={sortedJets[0].Et()*cosh(sortedJets[0].rap()),sortedJets[0].px(),sortedJets[0].py(),sortedJets[0].Et()*sinh(sortedJets[0].rap())};
    double tau1a=0;
    for (int i=1;i<=p.upper();i++){
        tau1a+=min(p[i].T()*qB[0]-p[i].X()*qB[1]-p[i].Y()*qB[2]-p[i].Z()*qB[3],p[i].T()*qJ[0]-p[i].X()*qJ[1]-p[i].Y()*qJ[2]-p[i].Z()*qJ[3]);
    }
    tau1a=tau1a*2/Q2;
    double alem = 1.0/128.0;//xalpha_em(Q2);
    double coef = 389379323000*alem*alem;
    weight_dis wtc = amp(&pdf, hardscale, hardscale, coef);
    weight_dis wtu = amp(&pdf, hardscale1, hardscale1, coef);
    weight_dis wtl = amp(&pdf, hardscale2, hardscale2, coef);
    physfilld(1, tau1a, wtc);
    physfilld(2, tau1a, wtu);
    physfilld(3, tau1a, wtl);
    physfilld(4, tau1a, wtc/_S_conv(wtc)*(dtau1amax-dtau1amin)/nbin);

}


