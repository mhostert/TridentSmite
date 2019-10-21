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
int Momentum_contains_nans(const std::vector<long double> &P){
 
  if(P[0]!=P[0]||P[1]!=P[1]||P[2]!=P[2]||P[3]!=P[3])
    {
      return 1;
    } 
  return 0;
}

void pretty_print(std::string CHANNEL){  
  std::cout<<"TRIDENT CHANNEL:\n"<<CHANNEL<<"\n\n";
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




trident_channel::trident_channel(int C){
  switch (C)
    {
      case eeee_flag:
        nu_alpha = e_flag;
        l1 = e_flag;
        l2 = e_flag;
        channel_name = "nue -> nue e+ e-";
        break;
      case mmmm_flag:
        nu_alpha = mu_flag;
        l1 = mu_flag;
        l2 = mu_flag;
        channel_name = "numu -> numu mu+ mu-";
        break;
      case emme_flag:
        nu_alpha = e_flag;
        l1 = e_flag;
        l2 = mu_flag;
        channel_name = "nue -> numu mu+ e-";
        break;
      case mmee_flag:
        nu_alpha = mu_flag;
        l1 = e_flag;
        l2 = e_flag;
        channel_name = "numu -> numu e+ e-";
        break;
      case eemm_flag:
        nu_alpha = e_flag;
        l1 = mu_flag;
        l2 = mu_flag;
        channel_name = "nue -> nue mu+ mu-";
        break;
      case meem_flag:
        nu_alpha = mu_flag;
        l1 = mu_flag;
        l2 = e_flag;
        channel_name = "numu -> nue e+ mu-";
        break;
      
      case eett_flag:
        nu_alpha = e_flag;
        l1 = tau_flag;
        l2 = tau_flag;
        channel_name = "nue -> nue tau+ tau-";
        break;
      case mmtt_flag:
        nu_alpha = mu_flag;
        l1 = tau_flag;
        l2 = tau_flag;
        channel_name = "nue -> nue tau+ tau-";
        break;
      case ette_flag:
        nu_alpha = e_flag;
        l1 = e_flag;
        l2 = tau_flag;
        channel_name = "nue -> nutau tau+ e-";
        break;
      case ttee_flag:
        nu_alpha = tau_flag;
        l1 = e_flag;
        l2 = e_flag;
        channel_name = "nutau -> nutau e+ e-";
        break;
      case ttmm_flag:
        nu_alpha = tau_flag;
        l1 = mu_flag;
        l2 = mu_flag;
        channel_name = "nutau -> nutau mu+ mu-";
        break;
      case teet_flag:
        nu_alpha = tau_flag;
        l1 = tau_flag;
        l2 = e_flag;
        channel_name = "nutau -> nue e+ tau-";        
        break;

        case tmmt_flag:
        nu_alpha = tau_flag;
        l1 = tau_flag;
        l2 = mu_flag;
        channel_name = "nutau -> numu mu+ tau-";
        break;
      case tttt_flag:
        nu_alpha = tau_flag;
        l1 = tau_flag;
        l2 = tau_flag;
        channel_name = "nutau -> nutau tau+ tau-";
        break;
      case mttm_flag:
        nu_alpha = mu_flag;
        l1 = mu_flag;
        l2 = tau_flag;
        channel_name = "numu -> nutau tau+ mu-";
        break;
    }
}