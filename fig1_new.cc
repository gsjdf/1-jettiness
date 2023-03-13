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
#include "lha_pdf.h"

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
double dtau1max=13.1;
double dtau1min=0.1;
int nbin=130;
double pTjetMin=0.0;
double radius = 1.0; // jet cone radius
double Pjtmax=30.0;
double Pjtmin=20.0;
double yjmax=2.5;
double yjmin=-2.5;

void psinput(phasespace_dis *ps, double& el, double& eh, double& q2min, double& q2max,
             double& xmin, double& xmax, double& ymin, double& ymax)
{
    //   we use the default phase space generator
    ps = 0;

    el = 45.0;    // GeV
    eh = 45.0;     // GeV

    //  Q^2 cuts
    q2min = 0.0001;     // GeV^2
    q2max = 8100.0;  // GeV^2

    //   xB cuts
    xmin = 0.0;
    xmax = 1.0;

    //   y cuts
    ymin = 0.0;
    ymax = 1.0;
}

void UserDIS::initfunc(unsigned int njet)
{
    phys(1, "tau1c", nbin, dtau1min, dtau1max);
    phys(2, "tau1u", nbin, dtau1min, dtau1max);
    phys(3, "tau1l", nbin, dtau1min, dtau1max);
    phys(4, "number", nbin, dtau1min, dtau1max);
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
    //----------1-jettiness
    double Q2 = -((p[-1] - p[-2]).mag2());
    double xB=Q2/2.0/dot(p[hadron(0)],(p[-1] - p[-2]));
    double s = 4.0*p[-1].T()*p[hadron(0)].T();
    double qB[4]={xB*p[hadron(0)].T(),0.0,0.0,xB*p[hadron(0)].Z()};
    double QB = xB*sqrt(s);
    double qJ[4],QJ,Pj[4],Pjt,yj;
    for (int jetnum=0;jetnum<sortedJets.size();jetnum++){
        qJ[0]=sortedJets[jetnum].Et()*cosh(sortedJets[jetnum].rap());
        qJ[1]=sortedJets[jetnum].px();
        qJ[2]=sortedJets[jetnum].py();
        qJ[3]=sortedJets[jetnum].Et()*sinh(sortedJets[jetnum].rap());
        QJ=2*sortedJets[jetnum].Et()*cosh(sortedJets[jetnum].rap());
        Pj[0] =0.0;
        Pj[1] =0.0;
        Pj[2] =0.0;
        Pj[3] =0.0;
        for (int i=1;i<=p.upper();i++){
            if((p[i].T()*qB[0]-p[i].X()*qB[1]-p[i].Y()*qB[2]-p[i].Z()*qB[3]/QB)>(p[i].T()*qJ[0]-p[i].X()*qJ[1]-p[i].Y()*qJ[2]-p[i].Z()*qJ[3])/QJ){
                Pj[0]+=p[i].T();
                Pj[1]+=p[i].X();
                Pj[2]+=p[i].Y();
                Pj[3]+=p[i].Z();
            }
        }
        Pjt=sqrt(pow(Pj[1],2)+pow(Pj[2],2));
        yj=0.5*log((Pj[0]+Pj[3])/(Pj[0]-Pj[3]));
        if(Pjt>Pjtmax|| Pjt<Pjtmin||yj>yjmax||yj<yjmin)
        {if(jetnum<sortedJets.size()-1)continue;
            if(jetnum==sortedJets.size()-1)return;}
        break
    }
    
    double hardscale=pow(2*Pjt,2);
    //double hardscale=pow(2*Pjtmin,2);
    double hardscale1=hardscale/4;
    double hardscale2=hardscale*4;
    
    double tau1=0.0;
    for (int i=1;i<=p.upper();i++){
        tau1+=min((p[i].T()*qB[0]-p[i].X()*qB[1]-p[i].Y()*qB[2]-p[i].Z()*qB[3])/QB,(p[i].T()*qJ[0]-p[i].X()*qJ[1]-p[i].Y()*qJ[2]-p[i].Z()*qJ[3])/QJ);
    }
    double alem = 1.0/128.0;//xalpha_em(Q2);
    double coef = 389379323000*alem*alem;
    weight_dis wtc = amp(&pdf, hardscale, hardscale, coef);
    weight_dis wtu = amp(&pdf, hardscale1, hardscale1, coef);
    weight_dis wtl = amp(&pdf, hardscale2, hardscale2, coef);
    physfilld(1, tau1, wtc);
    physfilld(2, tau1, wtu);
    physfilld(3, tau1, wtl);
    physfilld(4, tau1, wtc/_S_conv(wtc)*(dtau1max-dtau1min)/nbin);
}

