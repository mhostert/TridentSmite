#include "integrator.h"



// Constructor defines nu flavours, trident channel and target nucleus
tridentMC::tridentMC(int C, int Z_arg, int A_arg, long double Mn_arg, std::vector<long double> params){
  


  // Is it an antineutrino initial state?
  if (C < 0)
  {
    C = -C;
    IS_NUBAR = 1; 
  }else{IS_NUBAR=0;}


  // Determine the trident channel
  trident_channel channel(C);
  nu_alpha = channel.nu_alpha;
  l1 = channel.l1;
  l2 = channel.l2;
  PDG_lp = channel.PDG_lp; 
  PDG_lm = channel.PDG_lm; 
  PDG_nu_inc = channel.PDG_nu_inc;
  PDG_nu_out = channel.PDG_nu_out;
  PDG_had = channel.PDG_had;

  // Now assign the charged lepton masses
  switch (abs(l1))
    {
      case e_flag:
        ml1 = m_e;
        break;
      case mu_flag:
        ml1 = m_mu;
        break;
      case tau_flag:
        ml1 = m_tau;
        break;
    }

  switch (abs(l2))
    {
      case e_flag:
        ml2 = m_e;
        break;
      case mu_flag:
        ml2 = m_mu;
        break;
      case tau_flag:
        ml2 = m_tau;
        break;
    }


  // Initialize the vector of integrals
  for (int i = 0; i < 36; ++i)
  {
    dsigma.push_back(0.0);
  }

  A = A_arg;
  Z = Z_arg;
  Mn = Mn_arg;

  for (int i = 0; i < 7; ++i)
  {
    terms.push_back(1.0);
  }

  //   Define the proper axial and vector coefficients for nu(nu_alpha), l-(l1) and l+(l2)
  if (nu_alpha == abs(l1))
  {
    if (abs(l1) == abs(l2)){Aijk = 0.5; Vijk = 0.5 + 2*sw2;}
    else if (abs(l1) != abs(l2)) {Aijk = 1.0; Vijk = 1.0;}
    else {printf("Error! Flags for the leptons not well defined or not listed.");}
  }
  else if (nu_alpha != abs(l1)){Aijk = -0.5; Vijk = -0.5 + 2*sw2;}
  else {printf("Error! Flags for the leptons not well defined or not listed.");}


  // ANTINEUTRINO CROSS SECTION!
  if (IS_NUBAR == 1)
  {
    Aijk = - Aijk;
  }


  //////////
  // BSM parameters
  CHARGE = 1.0; // always 1 -- can redefine coupling if needed.
  mzprime = params[1];
  gprimeV = params[2];
  gprimeA = params[3];

  // fixed allocation for physicar vars vector 
  xphys.resize(12, 0);

  // start with no events
  total_vegas_events = 0;
}

void tridentMC::open_flux_file(std::string flux_file, long double elow, long double eup){
  
  Evec.clear();
  dPHIdE.clear();
  
  my_flux.open(flux_file.c_str());
  my_flux.precision(17);

  ////////////////////////////////////////
  // WARNING: Flux file needs to have 7 columns
  std::vector<long double> line(7);
  std::string dummyline;

  long double f_max = 0.0;

  if (my_flux.is_open())
  {
    my_flux.clear();
    my_flux.seekg(0);

    std::getline(my_flux,dummyline);  // Skip header!
    my_flux.clear();
    my_flux.seekg(0);
    
    if (dummyline.front() == '#')
    {
      std::getline(my_flux,dummyline);  // Skip header!
    }


    while(!my_flux.eof())
    {
      for (std::vector<int>::size_type i = 0; i < line.size(); i++)
      {
        my_flux >> line[i];
      }
 
      Evec.push_back(line[0]);
      dPHIdE.push_back(line[ ( abs(nu_alpha)-11)/2 +1+3*IS_NUBAR]);
      
    }
  }
  else
  {
    std::cout << "Could not open the flux file!" << std::endl;
  }

  dE = Evec[1] - Evec[0];
  Ei = 0;

  Emax = eup;

  if (elow > (SQR(ml1 + ml2) + 2*(ml1 + ml2)*Mn)/2.0/Mn)
  {
     Emin = elow;
  }
  else{ Emin = (SQR(ml1 + ml2) + 2*(ml1 + ml2)*Mn)/2.0/Mn;}
}

////////////////////////////////////////
long double tridentMC::integrate_wflux(void * integrando, int ndim){

  SAMPLES_FLAG = NO_SAMPLES;

  for (std::vector<int>::size_type i = 0; i < terms.size(); i++)
  {
    terms[i] = 1.0;
  }
  terms[6] = SMandBSM;
  terms[3]     = Vijk*Vijk;
  terms[4]     = Aijk*Aijk;
  terms[5]     = Vijk*Aijk;
  if (Emax <=  2*(SQR(ml1 + ml2) + 2*(ml1 + ml2)*Mn)/2.0/Mn  )
  {
    integral[0] = 0.0;
    error[0] = 0.0;
    chi2prob[0] = 0.0;
    comp = 0;
  }else{
      // Evaluate the integral
    // Vegas(ndim, NCOMP, (integrand_t)integrando, (void *)this, NVEC,
    // EPSREL, EPSABS, VERBOSE, SEED,
    // MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
    // GRIDNO, STATEFILE, SPIN,
    // &neval, &fail, integral, error, chi2prob);

    Suave(ndim, NCOMP, (integrand_t)integrando, (void *)this, NVEC,
    EPSREL, EPSABS, 0, SEED,
    MINEVAL, MAXEVAL, NNEW, NMIN,
    FLATNESS, STATEFILE, SPIN, &regions,
    &neval, &fail, integral, error, chi2prob);

  }

  // Close the input and output files
  my_flux.close();

  return integral[0];

}

////////////////////////////////////////
long double tridentMC::integrate_wflux_wsamples(void * integrando, int ndim, std::string samplesfile){


  // turn on all contributions
  for (std::vector<int>::size_type i = 0; i < terms.size(); i++)
  {
    terms[i] = 1.0;
  }
  terms[6]     = SMandBSM;
  terms[3]     = Vijk*Vijk;
  terms[4]     = Aijk*Aijk;
  terms[5]     = Vijk*Aijk;

  if (Emax <=  2*(SQR(ml1 + ml2) + 2*(ml1 + ml2)*Mn)/2.0/Mn  )
  {
    integral[0] = 0.0;
    error[0] = 0.0;
    chi2prob[0] = 0.0;
    comp = 0;
  }else{
    counter = 0;


    // // Evaluate the integral
    Vegas(ndim, NCOMP, (integrand_t)integrando, (void *)this, NVEC,
    EPSREL, EPSABS, VERBOSE, SEED,
    MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
    GRIDNO, STATEFILE, SPIN,
    &neval, &fail, integral, error, chi2prob);

    // Suave(ndim, NCOMP, (integrand_t)integrando, (void *)this, NVEC,
    // EPSREL, EPSABS, VERBOSE, SEED,
    // MINEVAL, MAXEVAL, NNEW, NMIN,
    // FLATNESS, STATEFILE, SPIN, &regions,
    // &neval, &fail, integral, error, chi2prob);
    // std::cout<<"\n\nNeeded "<<regions<<" regions and "<<neval<<" evaluations in Suave (my counter ="<<counter<<")"<<std::endl;
    
    std::cout<<"Integral: "<<integral[0]<<"+-"<<error[0]<<", Chi2 prob: "<<chi2prob[0]<<std::endl;

  }

  // Close the input and output files
  my_flux.close();

  return integral[0];
}

////////////////////////////////////////
long double tridentMC::integrate_energy(void * integrando, long double Enu, int ndim){

  SAMPLES_FLAG = NO_SAMPLES;

  nu_energy = Enu;
  // What component of the vector integrand to use (not applicable to use yet)
  if (nu_energy <=  2*(SQR(ml1 + ml2) + 2*(ml1 + ml2)*Mn)/2.0/Mn  )
  {

    integral[0] = 0.0;
    error[0] = 0.0;
    chi2prob[0] = 0.0;
    comp = 0;

  }else{

    // Evaluate the integral
    Vegas(ndim, NCOMP, (integrand_t)integrando, (void *)this, NVEC,
    EPSREL, EPSABS, VERBOSE, SEED,
    MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
    GRIDNO, STATEFILE, SPIN,
    &neval, &fail, integral, error, chi2prob);
  }
  // Close the input and output files
  return integral[0];
}

////////////////////////////////////////
long double tridentMC::integrate_energy_wsamples(void * integrando, long double Enu, int ndim, std::string samplesfile){

  nu_energy = Enu;
  // What component of the vector integrand to use (not applicable to use yet)
  if (nu_energy <=  2*(SQR(ml1 + ml2) + 2*(ml1 + ml2)*Mn)/2.0/Mn  )
  {

    integral[0] = 0.0;
    error[0] = 0.0;
    chi2prob[0] = 0.0;
    comp = 0;

  }else{

    // Evaluate the integral
    Vegas(ndim, NCOMP, (integrand_t)integrando, (void *)this, NVEC,
    EPSREL, EPSABS, VERBOSE, SEED,
    MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
    GRIDNO, STATEFILE, SPIN,
    &neval, &fail, integral, error, chi2prob);
  }

  return integral[0];
}

////////////////////////////////////////
long double tridentMC::compute_total_xsec_energy(void * integrando, long double Enu, int ndim){

  SAMPLES_FLAG = NO_SAMPLES;

  nu_energy = Enu;

  // turn on all contributions
  for (std::vector<int>::size_type i = 0; i < terms.size(); i++)
  {
    terms[i] = 1.0;
  }
  terms[6] = SMandBSM;
  terms[3]     = Vijk*Vijk;
  terms[4]     = Aijk*Aijk;
  terms[5]     = Vijk*Aijk;
  // What component of the vector integrand to use (not applicable to use yet)
  if (nu_energy <=  2*(SQR(ml1 + ml2) + 2*(ml1 + ml2)*Mn)/2.0/Mn  )
  {

    integral[0] = 0.0;
    error[0] = 0.0;
    chi2prob[0] = 0.0;
    comp = 0;

  }else{

    // Evaluate the integral
    Vegas(ndim, NCOMP, (integrand_t)integrando, (void *)this, NVEC,
    EPSREL, EPSABS, VERBOSE, SEED,
    MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
    GRIDNO, STATEFILE, SPIN,
    &neval, &fail, integral, error, chi2prob);
  }

  // reset terms contributions
  for (std::vector<int>::size_type i = 0; i < terms.size(); i++)
  {
    terms[i] = 0.0;
  }

  return integral[0];
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////
int tridentMC::compute_P_LAB()
{
  long double Enu,Q2,w,f;
  int iter,maxiter;

  // last vegas iteration
  maxiter = (int) cuba_samples.back().back();

  // total number of events in vegas
  total_vegas_events = cuba_samples.size();

  vector3d_of_P.resize(total_vegas_events);
  // total_vegas_events = make_3d_vector(int z, int y, int x, long double C = 0)
  //////////////////////////////////////////////////////////
  // WRITING EVENT FILES WITH WEIGHT
  std::vector<long double> sums(maxiter,0), sums2(maxiter,0), norm(maxiter,0);
  std::vector<int> evals(maxiter,0);
  std::vector<long double> avg(maxiter,0),std(maxiter,0);

  for(int i =0; i< total_vegas_events; i++)
  {
    iter = cuba_samples[i].back();
    f = cuba_samples[i][9];
    w = cuba_samples[i][10];

    if (iter <= maxiter)
    {
      sums[iter-1] += f*w; // real weight
      sums2[iter-1] += (f*f*w*w); 
      evals[iter-1]++;      // norm[iter-1]+= 1.0/w;
    }
    else std::cout << "ERROR! Iter is greater than the maximum found." <<std::endl;
  }

  long double avg_f = 0.0;
  long double std_f = 0.0;

  for (int i = 0; i < maxiter; i++)
  {
    avg[i] = sums[i];
    std[i] = ((evals[i])*sums2[i]-(sums[i])*(sums[i]))/(evals[i]-1);

    std_f += 1.0/std[i];
    avg_f += avg[i]/std[i];

  }

  // computing the correction factor
  std_f = 1.0/(std_f);
  avg_f *= std_f;
  weight_correction = std_f / (std[iter-1]);
  std::cout<<"correction: "<< std::endl;
  

  int count_nans = 0;

  for (int i = 0; i < total_vegas_events; i++)
  {
      vector3d_of_P[i] = P_LAB(ml1/*l-*/, ml2/*l-*/, Mn, cuba_samples[i]);

      f = cuba_samples[i][9];
      w = cuba_samples[i][10];
      Enu = cuba_samples[i][8];
      Q2 = cuba_samples[i][0];

      // destroy NaN entries and record 
      if ( (vector3d_contains_nans(vector3d_of_P)) || (is_nan(Q2)) || (is_nan(Enu)) || (is_nan(f*w * std_f / (std[iter-1]))) )
      {
        vector3d_of_P[i].clear();
        count_nans++;
      }
    }
  
  std::cout<<"Total number of events: "<<total_vegas_events<<std::endl;
  std::cout<<"Number of NaNed events: "<<count_nans<<std::endl;

  return 0;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////
int tridentMC::generate_events(std::string eventsfile)
{

  //////////////////////////////////////////////////////////
  // WRITING EVENT FILES WITH WEIGHTS  
  std::ofstream outfile(eventsfile.c_str());
  std::cout<<"opening new eventsfile..."<<eventsfile.c_str()<<std::endl;
  outfile.precision(17);

  outfile <<"# Units of GeV for momenta and (zb = 1e-45 cm^2)*(neutrino flux units) for the weights."<<std::endl;
  outfile <<"# Enu Q^2 P+(0) P+(1) P+(2) P+(3) P-(0) P-(1) P-(2) P-(3) weight"<<std::endl;


  long double Enu,Q2,w,f;
  for (int i = 0; i < total_vegas_events; i++)
  {
   
    // still want to print some low level information to file
    f = cuba_samples[i][9];
    w = cuba_samples[i][10];
    Enu = cuba_samples[i][8];
    Q2 = cuba_samples[i][0];

    outfile << Enu  << " "
    << Q2 << ""
    << vector3d_of_P[i][0][0]<< " " // Pnu
    << vector3d_of_P[i][0][1]<< " " // Pnu
    << vector3d_of_P[i][0][2]<< " " // Pnu
    << vector3d_of_P[i][0][3]<< " " // Pnu
    << vector3d_of_P[i][1][0]<< " " // Pplus
    << vector3d_of_P[i][1][1]<< " " // Pplus
    << vector3d_of_P[i][1][2]<< " " // Pplus
    << vector3d_of_P[i][1][3]<< " " // Pplus
    << vector3d_of_P[i][2][0]<< " " // Pminus
    << vector3d_of_P[i][2][1]<< " " // Pminus
    << vector3d_of_P[i][2][2]<< " " // Pminus
    << vector3d_of_P[i][2][3]<< " " // Pminus 
    << f*w*weight_correction << std::endl;

    }

  outfile.close();
  std::cout<<"Written vegas events to file."<< std::endl;

  return 0;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////
int tridentMC::HEPevt_format(std::string HEPevtfile, int NUMBER_OF_EVENTS)
{
  ///////////////////////////////////////////////////////
  // READING SAMPLES FILE A
  long double Enu,Q2,wf;
  
  int line_count=0;
  long double sum_wf = 0;

  for (int i = 0; i < total_vegas_events; i++)
  {
    sum_wf += cuba_samples[i].back();
  }
  std::cout<<"Total integral from events = "<<sum_wf<<std::endl;
  //////////////////////////////////////////////////////////

  // WRITING EVENTS TO FILE AFTER ACCEPT/REJECT
  std::ofstream outfile(HEPevtfile.c_str());
  std::cout<<"opening HEPevtfile..."<<HEPevtfile.c_str()<<std::endl;
  outfile.precision(17);

  int i = 0;
  while (i < NUMBER_OF_EVENTS)
  {
    wf = cuba_samples[i].back();
    // Accept/reject to  
    if ( wf/sum_wf > UniformRand() )
    {
      
      Enu = cuba_samples[i][8];
      
      // event header
      outfile <<i<<" "<<5<<std::endl;
      // incoming neutrino
      outfile <<"0 "<<PDG_nu_inc<<" 0 0 0 0 "<<Enu<<" 0.0 0.0 "<<Enu
                                            <<" 0.0 0.0 0.0 0.0 0.0"<<std::endl;
      /////////////////////////////////////////////////
      // !!!! FIX ME -- FIX THE TARGET HADRON INFO !!!!!
      outfile <<"0 "<<PDG_lm<<" 0 0 0 0 "<<vector3d_of_P[i][2][0]<<" "
                                          <<vector3d_of_P[i][2][1]<<" "
                                          <<vector3d_of_P[i][2][2]<<" "
                                          <<vector3d_of_P[i][2][3]<<" "
                                          <<ml1<<" 0.0 0.0 0.0 0.0"<<std::endl;
      // outgoing neutrino
      outfile <<"1 "<<PDG_nu_out<<" 0 0 0 0 "<<vector3d_of_P[i][0][0]<<" "
                                          <<vector3d_of_P[i][0][1]<<" "
                                          <<vector3d_of_P[i][0][2]<<" "
                                          <<vector3d_of_P[i][0][3]<<
                                          " 0.0 0.0 0.0 0.0 0.0"<<std::endl;
      // positive charged lepton
      outfile <<"1 "<<PDG_lp<<" 0 0 0 0 "<<vector3d_of_P[i][1][0]<<" "
                                          <<vector3d_of_P[i][1][1]<<" "
                                          <<vector3d_of_P[i][1][2]<<" "
                                          <<vector3d_of_P[i][1][3]<<" "
                                          <<ml2<<" 0.0 0.0 0.0 0.0"<<std::endl;
      // negative charged lepton
      outfile <<"1 "<<PDG_lm<<" 0 0 0 0 "<<vector3d_of_P[i][2][0]<<" "
                                          <<vector3d_of_P[i][2][1]<<" "
                                          <<vector3d_of_P[i][2][2]<<" "
                                          <<vector3d_of_P[i][2][3]<<" "
                                          <<ml1<<" 0.0 0.0 0.0 0.0"<<std::endl;
      /////////////////////////////////////////////////
      // !!!! FIX ME -- FIX THE STRUCK HADRON INFO !!!!!
      outfile <<"1 "<<PDG_lm<<" 0 0 0 0 "<<vector3d_of_P[i][2][0]<<" "
                                          <<vector3d_of_P[i][2][1]<<" "
                                          <<vector3d_of_P[i][2][2]<<" "
                                          <<vector3d_of_P[i][2][3]<<" "
                                          <<ml1<<" 0.0 0.0 0.0 0.0"<<std::endl;
      i++;
    }
      
  }
  
  outfile.close();

  return 0;
}

