#include "constants.h"


/* ********************************************************* */
// RANDOM NUMBER GENERATORS
/////////////////////////////////////////////////////
long double UniformRand(){
  return ( (long double)(rand()) + 1. )/( (long double)(RAND_MAX) + 1. );
}

///////////////////////////////////////////////////////
long double NormalRand(long double mean, long double stddev)
{//Box muller method
    static long double n2 = 0.0;
    static int n2_cached = 0;
    if (!n2_cached)
    {
        long double x, y, r;
        do
        {
            x = 2.0*rand()/RAND_MAX - 1;
            y = 2.0*rand()/RAND_MAX - 1;

            r = x*x + y*y;
        }
        while (r == 0.0 || r > 1.0);
        {
            long double d = sqrt(-2.0*log(r)/r);
            long double n1 = x*d;
            n2 = y*d;
            long double result = n1*stddev + mean;
            n2_cached = 1;
            return result;
        }
    }
    else
    {
        n2_cached = 0;
        return n2*stddev + mean;
    }
}


/* ********************************************************* */
// HELPER FUNCTIONS
int vector3d_contains_nans(const std::vector<std::vector<std::vector<long double>>> &P){
    for (int i = 0; i < P.size(); i++)
    {
        for (int j = 0; j < P[i].size(); j++)
        {
            for (int k = 0; k < P[i][j].size(); k++)
            {           
                if (P[i][j][k]!=P[i][j][k])
                {
                    return 1;
                }
            }
        }
    }
    return 0;
}
int vector_of_vectors_contains_nans(const std::vector<std::vector<long double>> &P){
    for (int i = 0; i < P.size(); i++)
    {
        for (int j = 0; j < P[i].size(); j++)
        {
            if (P[i][j]!=P[i][j])
            {
                return 1;
            }
        }
    }
    return 0;
}
int vector_contains_nans(const std::vector<long double> &P){
    for (int i = 0; i < P.size(); i++)
    {
        if (P[i]!=P[i])
        {
            return 1;
        }
    }
    return 0;
}
int is_nan(long double x){
  if(x!=x)
    {
      return 1;
    } 
  return 0;
}

void pretty_print(std::string CHANNEL){  
  std::cout<<"Trident channel: "<<CHANNEL<<"\n\n";
}


std::vector<std::string> list_of_trident_channels(std::vector<std::string> &sList){
  sList.push_back("eeee");
  sList.push_back("mmmm");
  sList.push_back("emme");
  sList.push_back("mmee");
  sList.push_back("eemm");
  sList.push_back("meem");
  sList.push_back("eett");
  sList.push_back("mmtt");
  sList.push_back("tttt");
  sList.push_back("ette");
  sList.push_back("mttm");
  sList.push_back("teet");
  sList.push_back("tmmt");
  sList.push_back("ttee");
  sList.push_back("ttmm");

  return sList;
}

///////////////////////////////////////////////////////
std::string getLastLine(std::ifstream& in)
{
    std::string line;
    while (in >> std::ws && std::getline(in, line)) // skip empty lines
        ;

    return line;
}

///////////////////////////////////////////////////////
std::vector<std::vector<std::vector<long double>>> make_3d_vector(int z, int y, int x, long double C = 0)
{
    return std::vector<std::vector<std::vector<long double>>>(z, std::vector<std::vector<long double>>(y, std::vector<long double>(x, C)));
}

/////////////////////////////////////
trident_channel::trident_channel(int C){
  switch (abs(C))
    {
      case eeee_flag:
        PDG_nu_out = nue_flag;
        PDG_nu_inc = nue_flag;
        PDG_lm = e_flag;
        PDG_lp = -e_flag;
        nu_alpha = e_flag;
        l1 = e_flag;
        l2 = e_flag;
        nu1_name="nue";
        nu2_name="nue";
        l1_name="e+";
        l2_name="e-";
        break;
      case mmmm_flag:
        PDG_nu_out = numu_flag;
        PDG_nu_inc = numu_flag;
        PDG_lm = mu_flag;
        PDG_lp = -mu_flag;
        nu_alpha = mu_flag;
        l1 = mu_flag;
        l2 = mu_flag;
        nu1_name="numu";
        nu2_name="numu";
        l1_name="mu+";
        l2_name="mu-";
        break;
      case emme_flag:
        PDG_nu_out = numu_flag;
        PDG_nu_inc = nue_flag;
        PDG_lm = e_flag;
        PDG_lp = -mu_flag;
        nu_alpha = e_flag;
        l1 = e_flag;
        l2 = mu_flag;
        nu1_name="nue";
        nu2_name="numu";
        l1_name="mu+";
        l2_name="e-";
        break;
      case mmee_flag:
        PDG_nu_out = numu_flag;
        PDG_nu_inc = numu_flag;
        PDG_lm = e_flag;
        PDG_lp = -e_flag;
        nu_alpha = mu_flag;
        l1 = e_flag;
        l2 = e_flag;
        nu1_name="numu";
        nu2_name="numu";
        l1_name="e+";
        l2_name="e-";
        break;
      case eemm_flag:
        PDG_nu_out = nue_flag;
        PDG_nu_inc = nue_flag;
        PDG_lm = mu_flag;
        PDG_lp = -mu_flag;
        nu_alpha = e_flag;
        l1 = mu_flag;
        l2 = mu_flag;
        nu1_name="nue";
        nu2_name="nue";
        l1_name="mu+";
        l2_name="mu-";
        break;
      case meem_flag:
        PDG_nu_out = nue_flag;
        PDG_nu_inc = numu_flag;
        PDG_lm = mu_flag;
        PDG_lp = -e_flag;
        nu_alpha = mu_flag;
        l1 = mu_flag;
        l2 = e_flag;
        nu1_name="numu";
        nu2_name="nue";
        l1_name="e+";
        l2_name="mu-";
        break;
      
      case eett_flag:
        PDG_nu_out = nue_flag;
        PDG_nu_inc = nue_flag;
        PDG_lm = tau_flag;
        PDG_lp = -tau_flag;
        nu_alpha = e_flag;
        l1 = tau_flag;
        l2 = tau_flag;
        nu1_name="nue";
        nu2_name="nue";
        l1_name="tau+";
        l2_name="tau-";
        break;
      case mmtt_flag:
        PDG_nu_out = nue_flag;
        PDG_nu_inc = numu_flag;
        PDG_lm = tau_flag;
        PDG_lp = -tau_flag;
        nu_alpha = mu_flag;
        l1 = tau_flag;
        l2 = tau_flag;
        nu1_name="nue";
        nu2_name="nue";
        l1_name="tau+";
        l2_name="tau-";
        break;
      case ette_flag:
        PDG_nu_out = nutau_flag;
        PDG_nu_inc = nue_flag;
        PDG_lm = e_flag;
        PDG_lp = -tau_flag;
        nu_alpha = e_flag;
        l1 = e_flag;
        l2 = tau_flag;
        nu1_name="nue";
        nu2_name="nutau";
        l1_name="tau+";
        l2_name="e-";
        break;
      case ttee_flag:
        PDG_nu_out = nutau_flag;
        PDG_nu_inc = nutau_flag;
        PDG_lm = e_flag;
        PDG_lp = -e_flag;
        nu_alpha = tau_flag;
        l1 = e_flag;
        l2 = e_flag;
        nu1_name="nutau";
        nu2_name="nutau";
        l1_name="e+";
        l2_name="e-";
        break;
      case ttmm_flag:
        PDG_nu_out = nutau_flag;
        PDG_nu_inc = nutau_flag;
        PDG_lm = mu_flag;
        PDG_lp = -mu_flag;
        nu_alpha = tau_flag;
        l1 = mu_flag;
        l2 = mu_flag;
        nu1_name="nutau";
        nu2_name="nutau";
        l1_name="mu+";
        l2_name="mu-";
        break;
      case teet_flag:
        PDG_nu_out = nue_flag;
        PDG_nu_inc = nutau_flag;
        PDG_lm = tau_flag;
        PDG_lp = -e_flag;
        nu_alpha = tau_flag;
        l1 = tau_flag;
        l2 = e_flag;
        nu1_name="nutau";
        nu2_name="nue";
        l1_name="e+";
        l2_name="tau-";        
        break;

        case tmmt_flag:
        PDG_nu_out = numu_flag;
        PDG_nu_inc = nutau_flag;
        PDG_lm = tau_flag;
        PDG_lp = -mu_flag;
        nu_alpha = tau_flag;
        l1 = tau_flag;
        l2 = mu_flag;
        nu1_name="nutau";
        nu2_name="numu";
        l1_name="mu+";
        l2_name="tau-";
        break;
      case tttt_flag:
        PDG_nu_out = nutau_flag;
        PDG_nu_inc = nutau_flag;
        PDG_lm = tau_flag;
        PDG_lp = -tau_flag;
        nu_alpha = tau_flag;
        l1 = tau_flag;
        l2 = tau_flag;
        nu1_name="nutau";
        nu2_name="nutau";
        l1_name="tau+";
        l2_name="tau-";
        break;
      case mttm_flag:
        PDG_nu_out = nutau_flag;
        PDG_nu_inc = numu_flag;
        PDG_lm = mu_flag;
        PDG_lp = -tau_flag;
        nu_alpha = mu_flag;
        l1 = mu_flag;
        l2 = tau_flag;
        nu1_name="numu";
        nu2_name="nutau";
        l1_name="tau+";
        l2_name="mu-";
        break;
    }

    if (C<0)
    {
      channel_name=nu1_name+"bar_to_"+nu2_name+"bar_"+l1_name+"_"+l2_name;
      PDG_nu_inc = -PDG_nu_inc;
      PDG_nu_out = -PDG_nu_out;    
    }
    else{
      channel_name=nu1_name+"_to_"+nu2_name+"_"+l1_name+"_"+l2_name;
    }

    ////////////////
    // !!!!! FIX ME !!!!!!!!
    PDG_had = 1;

}