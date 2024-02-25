#include <boost/program_options.hpp>
#include <filesystem>
namespace po = boost::program_options;
#include "integrator.h"
#include <sys/stat.h>

using namespace std::chrono;

///////////////////////////////////////////////////////////////////////////
// MAIN 
int main(int argc, char* argv[]) {


  int C, A, Z, n_vegas, n_hepevt_events;
  long double mzprime, gprimeV, gprimeA;
  long double Emin, Emax;
  std::string fluxfile, dir_path;

  po::options_description desc("Allowed options");
  desc.add_options()
    ("help,h", 
        "produce help message")
    // number of VEGAS events
    ("nevents,N", po::value<int>(& n_vegas)->default_value(1e5),
        "Number of VEGAS events to compute (not number of vegas events)")
    ("text", "Print all vegas events (weighted) in text (.dat) format")
    ("nonpy", "Do not print all vegas events (weighted) in numpy (.npy) format")

    // number of HEPevt events to produce
    ("n_hep_events,Nhep", po::value<int>(& n_hepevt_events)->default_value(1e2),
        "Number of HEPevt events to generate (not number of vegas events)")
    ("hepevt", "Print unweighted events in hepevt format")
    
    // trident channel
    ("channel,c", po::value<int>(& C)->default_value(1),
        "Trident channel to use (see README for definition)")

    // Hadronic target
    ("znumber,z", po::value<int>(& Z)->default_value(1),
        "Target proton number")
    ("anumber,a", po::value<int>(& A)->default_value(1),
        "Target mass number (e.g. A = 12 for Carbon)")
    ("pb,b",
        "Pauli blocking (only include if computing scattering on bound protons)")
    
    // neutrino energies
    ("emin,l", po::value<long double>(& Emin)->default_value(1.0),
        "Minimum Enu to sample from")
    ("emax,u", po::value<long double>(& Emax)->default_value(2.0),
        "Maximum Enu to sample from")
    ("fluxfile,f", po::value<std::string>(& fluxfile)->default_value("./fluxes/uniform_0.1_200_GeV.dat"),
        "Neutrino flux file to use")
    
    ("path,dir", po::value<std::string>(& dir_path)->default_value("./data/default/"),
        "Neutrino flux file to use")

    // BSM params
    ("mzprime,m", po::value<long double>(& mzprime)->default_value(1.0),
        "Zprime mass")      
    ("gprimeV,gV", po::value<long double>(& gprimeV)->default_value(0.0),
        "Zprime vector coupling")
    ("gprimeA,gA", po::value<long double>(& gprimeA)->default_value(0.0),
        "Zprime axial coupling");

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
  // pass BSM params
  long double myparams[] = {ZPRIME,mzprime,gprimeV,gprimeA};
  std::vector<long double> my_BSM (myparams, myparams + sizeof(myparams) / sizeof(long double));

  // create object
  tridentMC MC(C, Z, A, (long double) (A*m_AVG), my_BSM, n_vegas);
  
  // if 4, does the sum of diagrams
  MC.TOTAL_DIAGRAMS = 4;

  // FLUX INTEGRATED
  MC.E_FLAG = WFLUX;

  MC.open_flux_file(fluxfile, (long double) Emin, (long double) Emax);
  //////////////////////////////////////////////


  // List of names
  std::string target;
 
  if (Z>1){
      //////////////// COHERENT ////////////////////////      
      MC.integrate_wflux_wsamples((void *)coh_BSM_flux, 9, samplesfile);
      target="coh_"+std::to_string(Z)+"_"+std::to_string(A);
  }
  else if (Z==1)
  {
      //////////////// ELASTIC PROTON ////////////////////////
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
  else if (Z==0){
      //////////////// ELASTIC NEUTRON ////////////////////////
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
  else
  {
    std::cout << "Could not identify hadronic target, A="<< A << " Z="<<Z << std::endl;
    return 1;
  }

  ////////////////////////////////////
  // reconstruct kinematics in LAB frame
  MC.generate_kinematics();
 

  std::string gprimeVstring = std::to_string((float) gprimeV); // s == "1e+05"
  std::string gprimeAstring = std::to_string((float) gprimeA); // s == "1e+05"
  std::string mzprimestring = std::to_string((float) mzprime); // s == "1e+05"

  // auto start = std::chrono::system_clock::now();
  // SAVE EVENTS TO FILE
  if (vm.count("hepevt"))
  {
    std::filesystem::create_directories("hepevt/");
    std::string hepevtfile="HEPevt/MC_events_"+channel.channel_name+"_"+target+"_gv_"+gprimeVstring+"_ga_"+gprimeAstring+"_mz_"+mzprimestring+".dat";
    std::cout<<"Saving events to file: "<<hepevtfile<<std::endl;
    MC.save_to_HEPevt(hepevtfile,n_hepevt_events);
  }
  
  if (vm.count("text"))
  {
    std::filesystem::create_directories(dir_path.c_str());
    std::string eventsfile=dir_path+"MC_events_"+channel.channel_name+"_"+target+"_gv_"+gprimeVstring+"_ga_"+gprimeAstring+"_mz_"+mzprimestring+".dat";
    std::cout<<"Saving events to file: "<<eventsfile<<std::endl;
    MC.save_to_text(eventsfile);
  }

  if (!(vm.count("nonpy")))
  {
    std::filesystem::create_directories(dir_path.c_str());
    std::string eventsfile=dir_path+"MC_events_"+channel.channel_name+"_"+target+"_gv_"+gprimeVstring+"_ga_"+gprimeAstring+"_mz_"+mzprimestring+".npy";
    std::cout<<"Saving events to file: "<<eventsfile<<std::endl;
    MC.save_to_npy(eventsfile);
  }

  // record end time
  // auto end = std::chrono::system_clock::now();
  // std::chrono::duration<double> diff = end-start;
  // std::cout << "Time to write to file " << " ints : " << diff.count() << " s\n";
  return 0;
}
