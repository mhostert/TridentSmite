#include <boost/program_options.hpp>
namespace po = boost::program_options;
#include "integrator.h"

///////////////////////////////////////////////////////////////////////////
// MAIN 
///////////////////////////////////////////////////////////////////////////


int main(int argc, char* argv[]) {


  int C,A,Z,Block;
  std::string run_number;
  long double Mz,gl,gnu;
  long double Emin,Emax;

  po::options_description desc("Allowed options");
  desc.add_options()
      ("help,h", 
          "produce help message")
      ("channel,c", po::value<int>(& C)->default_value(1),
          "Trident channel")
      ("znumber,z", po::value<int>(& Z)->default_value(1),
          "Proton number")
      ("anumber,a", po::value<int>(& A)->default_value(1),
          "Atomic number")
      ("pb,b",
          "Pauli Blocking")
      ("emin,l", po::value<long double>(& Emin)->default_value(0.1),
          "Minimum Enu")
      ("emax,u", po::value<long double>(& Emax)->default_value(40.0),
          "Maximum Enu")
      ;

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);
  
  // Print help if needed
  if (vm.count("help")) {
      std::cout << desc << "\n";
      return 1;
  }

  /**********************/
  long double myparams[] = {ZPRIME,Mz,gl,gnu};
  std::vector<long double> my_BSM (myparams, myparams + sizeof(myparams) / sizeof(long double) );

  //////////////////////////////////
  // Trident channel of interest
  trident_channel channel(C);
  pretty_print(channel.channel_name);
  std::string samplesfile = "samples/temp_"+channel.channel_name+"_"+std::to_string(Z)+"_"+std::to_string(A);


  ///////////////////////////////////
  // FLUX FILE TO USE
  std::string fluxfile;
  fluxfile   = "fluxes/uniform_0.1_200_GeV.dat";

  ///////////////////////////////
  // Set up the MC 
  tridentMC MC(C, Z, A, (long double) (A*m_AVG));
  
  MC.params = my_BSM;

  // FLUX INTEGRATED
  MC.E_FLAG = WFLUX;
  // if 4, includes the sum of diagrams as well. Necessary because of numerical precision issues.
  MC.TOTAL_DIAGRAMS = 4;


  MC.open_flux_file(fluxfile, (long double) Emin, (long double) Emax);

//////////////// COHERENT ////////////////////////
  // List of names
  std::string target;
  switch (A){
    case 40: {
      MC.integrate_wflux_wsamples((void *)coh_BSM_flux, 9, samplesfile);
      target = "coh_40Ar";
      break;
    }
    case 16: {
      MC.integrate_wflux_wsamples((void *)coh_BSM_flux, 9, samplesfile);
      target = "coh_16O";
      break;
    }
    case 12: {
      MC.integrate_wflux_wsamples((void *)coh_BSM_flux, 9, samplesfile);
      target = "coh_12C";
      break;
    }

//////////////// ELASTIC PROTON ////////////////////////
    case 1: {

      if ( Z == 1 )
      {
        if (vm.count("pb"))
        {
          MC.PAULI_BLOCKING = W_BLOCKING;
          MC.integrate_wflux_wsamples((void *)dif_BSM_flux, 9, samplesfile);
          target  = "proton_PB";
        }
        else
        {
          MC.PAULI_BLOCKING = NO_BLOCKING;
          MC.integrate_wflux_wsamples((void *)dif_BSM_flux, 9, samplesfile);
          target  = "proton_noPB";
        }
      }


//////////////// ELASTIC NEUTRON ////////////////////////
      if ( Z == 0 )
      {
        if (vm.count("pb"))
        {
          MC.PAULI_BLOCKING = W_BLOCKING;
          MC.integrate_wflux_wsamples((void *)dif_BSM_flux, 9, samplesfile);
          target  = "neutron_PB";
        }
        else
        {
          MC.PAULI_BLOCKING = NO_BLOCKING;
          MC.integrate_wflux_wsamples((void *)dif_BSM_flux, 9, samplesfile);
          target  = "neutron_noPB";
        }
      }
      break;

    }
  }

  // SAVE EVENTS TO FILE
  MC.generate_events(samplesfile, "events/MC_events_"+channel.channel_name+"_"+target+"_"+std::to_string((int)Emin)+"_"+std::to_string((int)Emax)+"_GeV.dat");

  return 0;
}
