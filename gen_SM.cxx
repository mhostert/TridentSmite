#include "integrator.h"


///////////////////////////////////////////////////////////////////////////
// MAIN 
///////////////////////////////////////////////////////////////////////////


int main(int argc, char* argv[]) {

  // string with trident channel lists
  std::vector<std::string> sList;
  list_of_trident_channels(sList);

  int nu_alpha, l1, l2;
  int C,A,Z,Block;
  std::string run_number;
  long double Mz,gl,gnu;
  long double Emin,Emax;

  /**********************/
  if(argc < 6) {
    printf("You must provide at least 6 arguments!\\nUsage: \n\n./gen_SM CHANNEL Znumber Anumber Block Enumin Enumax\n\n");
    exit(0);
  }
  else 
  {  
    C = std::stoi(argv[1]);
    Z = std::stoi(argv[2]);
    A = std::atoi(argv[3]);
    Block = std::atoi(argv[4]);
    Emin = std::stof(argv[5]);
    Emax = std::atof(argv[6]);
  }
  /**********************/

  std::string fluxfile;
  std::string eventsfile_coh, eventsfile_p, eventsfile_n;

  long double myparams[] = {ZPRIME,Mz,gl,gnu};
  std::vector<long double> my_BSM (myparams, myparams + sizeof(myparams) / sizeof(long double) );


///////////////////////////////////////////////////////////////////////////////////

  trident_channel channel(C);
  pretty_print(channel.channel_name);
  std::string samplesfile = "samples/temp_"+sList[C]+"_"+std::to_string(Z)+"_"+std::to_string(A)+"_"+std::to_string(Block)+"_"+std::to_string(Emax);


  // fluxfile   = "fluxes/DUNE/DUNE_flux_ND_numode_459m_3horn_62.4GeV_1.83POT_2m.dat";
  // fluxfile   = "fluxes/DUNE/DUNE_flux_ND_nubarmode_459m_3horn_62.4GeV_1.83POT_2m.dat";
  fluxfile   = "fluxes/uniform_0.1_200_GeV.dat";


  ///////////////////////////////
  // Set up the MC 
  tridentMC MC(C, Z, A, (long double) (A*m_AVG));
  
  MC.params = my_BSM;
  MC.IS_NUBAR = 0;

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
        if (Block == 0 )
        {
          MC.PAULI_BLOCKING = NO_BLOCKING;
          MC.integrate_wflux_wsamples((void *)dif_BSM_flux, 9, samplesfile);
          target  = "proton_noPB";
        }
        if (Block == 1 )
        {
          MC.PAULI_BLOCKING = W_BLOCKING;
          MC.integrate_wflux_wsamples((void *)dif_BSM_flux, 9, samplesfile);
          target  = "proton_PB";
        }
      }

//////////////// ELASTIC NEUTRON ////////////////////////
      if ( Z == 0 )
      {
        if ( Block == 0 )
        {
          MC.PAULI_BLOCKING = NO_BLOCKING;
          MC.integrate_wflux_wsamples((void *)dif_BSM_flux, 9, samplesfile);
          target  = "neutron_noPB";
        }
        if ( Block ==1 )
        {
          MC.PAULI_BLOCKING = W_BLOCKING;
          MC.integrate_wflux_wsamples((void *)dif_BSM_flux, 9, samplesfile);
          target  = "neutron_PB";
        }
      }
      break;

    }
  }

  // SAVE EVENTS TO FILE
  MC.generate_events(samplesfile, "events/MC_events_"+sList[C]+"_"+target+"_"+std::to_string((int)Emin)+"_"+std::to_string((int)Emax)+"_GeV.dat");


  return 0;
}
