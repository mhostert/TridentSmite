#include <boost/program_options.hpp>
namespace po = boost::program_options;
#include "integrator.h"

///////////////////////////////////////////////////////////////////////////
// MAIN 
///////////////////////////////////////////////////////////////////////////


int main(int argc, char* argv[]) {


  int C,A,Z,Block,NUMBER_OF_EVENTS;
  std::string run_number;
  long double Mz,gl,gnu;
  long double Emin,Emax;
  std::string fluxfile;

  po::options_description desc("Allowed options");
  desc.add_options()
      ("help,h", 
          "produce help message")
      ("nevents,N", po::value<int>(& NUMBER_OF_EVENTS)->default_value(1e4),
          "Number of HEPevt events to generate")
     ("channel,c", po::value<int>(& C)->default_value(1),
          "Trident channel to use (see README for definition)")
      ("znumber,z", po::value<int>(& Z)->default_value(1),
          "Target proton number")
      ("anumber,a", po::value<int>(& A)->default_value(1),
          "Target mass number (e.g. A = 12 for Carbon)")
      ("pb,b",
          "Pauli blocking (only include if computing scattering on bound protons)")
      ("emin,l", po::value<long double>(& Emin)->default_value(0.10),
          "Minimum Enu to sample from")
      ("emax,u", po::value<long double>(& Emax)->default_value(40.0),
          "Maximum Enu to sample from")
      ("fluxfile,f", po::value<std::string>(& fluxfile)->default_value("fluxes/uniform_0.1_200_GeV.dat"),
          "Neutrino flux file to use")
      ("mzprime,m", po::value<long double>(& Mz)->default_value(1.0),
          "Zprime mass")      
      ("gprime,g", po::value<long double>(& gl)->default_value(0.0),
          "Zprime coupling")
      ;

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);
  
  // Print help if needed
  if (vm.count("help")) {
      std::cout << desc << "\n";
      return 1;
  }


  //////////////////////////////////
  // Trident channel of interest
  trident_channel channel(C);
  pretty_print(channel.channel_name);
  std::string samplesfile = "samples/temp_"+channel.channel_name+"_"+std::to_string(Z)+"_"+std::to_string(A);


  //////////////////////////////////////////////
  // Set up the MC 
  tridentMC MC(C, Z, A, (long double) (A*m_AVG));
  /**********************/
  // FIX ME -- NOT FUNCTIONAL YET!!
  long double myparams[] = {ZPRIME,Mz,gl,gnu};
  std::vector<long double> my_BSM (myparams, myparams + sizeof(myparams) / sizeof(long double) );
  MC.params = my_BSM;
  // if 4, does the sum of diagrams
  MC.TOTAL_DIAGRAMS = 4;

  // FLUX INTEGRATED
  MC.E_FLAG = WFLUX;

  MC.open_flux_file(fluxfile, (long double) Emin, (long double) Emax);
//////////////////////////////////////////////


  //////////////// COHERENT ////////////////////////
  // List of names
  std::string target;
  MC.integrate_wflux_wsamples((void *)coh_BSM_flux, 9, samplesfile);
  target="coh_"+std::to_string(Z)+"_"+std::to_string(A);

  //////////////// ELASTIC PROTON ////////////////////////
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


  // SAVE EVENTS TO FILE
  std::string eventsfile="events/MC_events_"+channel.channel_name+"_"+target+"_"+std::to_string((int)Emin)+"_"+std::to_string((int)Emax)+".dat";
  std::string hepevtfile="HEPevt/MC_events_"+channel.channel_name+"_"+target+".dat";

  MC.generate_events(samplesfile,eventsfile);
  MC.HEPevt_format(eventsfile,hepevtfile,NUMBER_OF_EVENTS);

  return 0;
}
