#include "integrator.h"



// Constructor defines nu flavours, trident channel and target nucleus
tridentMC::tridentMC(int C, int Z_arg, int A_arg, long double Mn_arg){
  

  if (C < 0)
  {
    C = -C;
    IS_NUBAR = 1; 
  }else{IS_NUBAR=0;}

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
  else if (nu_alpha != abs(l1))
  {
    Aijk = -0.5; Vijk = -0.5 + 2*sw2;
  }
  else {printf("Error! Flags for the leptons not well defined or not listed.");}


  // ANTINEUTRINO CROSS SECTION!!!
  if (IS_NUBAR == 1)
  {
    Aijk = - Aijk;
  }

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
    // std::getline(my_flux,dummyline);  // Skip header!

    while(!my_flux.eof())
    {
      for (std::vector<int>::size_type i = 0; i < line.size(); i++)
      {
        my_flux >> line[i];
      }

      Evec.push_back(line[0]);
      dPHIdE.push_back(line[ ( abs(nu_alpha)-12)/2 +1+3*IS_NUBAR]);
    
    }
  }
  else
  {
    std::cout << "Could not open the flux file!" << std::endl;
  }

  dE = Evec[1] - Evec[0];
  Ei = 0;

  Emax = eup;
  Emin = elow;
}

////////////////////////////////////////
long double tridentMC::integrate_wflux(void * integrando, int ndim){

  SAMPLES_FLAG = NO_SAMPLES;

  // for (std::vector<int>::size_type i = 0; i < terms.size(); i++)
  // {
    // terms[i] = 1.0;
  // }
  // terms[6] = SMonly;
  // terms[3]     = Vijk*Vijk;
  // terms[4]     = Aijk*Aijk;
  // terms[5]     = Vijk*Aijk;
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


  SAMPLES_FLAG = PRINT_SAMPLES;
  my_samples.open(samplesfile.c_str());
  my_samples.precision(17);
  // turn on all contributions
  for (std::vector<int>::size_type i = 0; i < terms.size(); i++)
  {
    terms[i] = 1.0;
  }
  terms[6] = SMonly;
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

  my_samples.close();
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

  SAMPLES_FLAG = PRINT_SAMPLES;
  my_samples.open(samplesfile.c_str());
  my_samples.precision(17);

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
  my_samples.close();
  std::cout<<"DONE WITH THE FIRST INTEGRAL"<<std::endl;
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
  terms[6] = SMonly;
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
void tridentMC::generate_events(std::string samplesfile, std::string eventsfile)
{

  ///////////////////////////////////////////////////////
  // READING SAMPLES FILE A
  std::ifstream ifile(samplesfile.c_str());
  ifile.precision(17);
  std::cout<<"reading integration file..."<<samplesfile.c_str()<<std::endl;

  long double x1,x2,x3,x4,x5,x6,x7,x8,x9,w,f,prob;
  int iter,maxiter;

  // Get the last number of iterations
  std::string line;
  if (ifile.is_open())
  {
    line = getLastLine(ifile);
    std::istringstream iss(line);
    iss >> x1 >> x2 >> x3 >> x4 >> x5 >> x6 >> x7 >> x8 >> x9 >> f >> w >> iter;
  }
  else{
    std::cout << "ERROR! Not able to read integration file." << std::endl;
  }

  maxiter = iter;

  //////////////////////////////////////////////////////////
  // WRITING EVENT FILES WITH WEIGHT

  std::vector<long double> sums(maxiter,0), sums2(maxiter,0), norm(maxiter,0);
  std::vector<int> evals(maxiter,0);
  std::vector<long double> avg(maxiter,0),std(maxiter,0);


  // Start reading file again
  ifile.clear();
  ifile.seekg(0);

  while(std::getline(ifile,line))
  {
    
    std::istringstream iss(line);
    if(!(iss >> x1 >> x2 >> x3 >> x4 >> x5 >> x6 >> x7 >> x8 >> x9 >> f >> w >> iter))
    {
      // NaN's in the sample file.
      continue;
      // FIX ME -- add a measure of NaNs here too.

    }

    if (iter <= maxiter)
    {
      sums[iter-1] += (f*w);
      sums2[iter-1] += (f*w)*(f*w);
      evals[iter-1]++;
      // norm[iter-1]+= 1.0/w;
    }else std::cout << "ERROR! Iter is greater than the maximum found." <<std::endl;
  }

  long double avg_f = 0.0;
  long double std_f = 0.0;

  for (int i = 0; i < maxiter; ++i)
  {
    avg[i] = sums[i];
    std[i] = ((evals[i])*sums2[i]-(sums[i])*(sums[i]))/(evals[i]-1);

    std_f += 1.0/std[i];
    avg_f += avg[i]/std[i];

  }

  std_f = 1.0/(std_f);
  avg_f *= std_f;
  
  // Start reading file again
  ifile.clear();
  ifile.seekg(0);


  //////////////////////////////////////////////////////////
  // WRITING EVENT FILES WITH WEIGHTS  

  std::ofstream outfile(eventsfile.c_str());
  std::cout<<"opening new eventsfile..."<<eventsfile.c_str()<<std::endl;
  outfile.precision(17);


  outfile <<"# MC events with uniform energy profile. All units are in GeV or GeV^2. "<<std::endl;
  outfile <<"# Enu Q^2 P+(0) P+(1) P+(2) P+(3) P-(0) P-(1) P-(2) P-(3) weight"<<std::endl;

  // Must have the same number of entries as columns in the sampes#.txt file!!
  std::vector<long double> row(9);
  long double mj = ml1; // l-
  long double mk = ml2; // l+
  std::vector<long double> Pm(4), Pp(4);


  int count_nans = 0;
  int counter2 = 0;

  while(std::getline(ifile,line))
  {
    std::istringstream iss(line);
    iss >> x1 >> x2 >> x3 >> x4 >> x5 >> x6 >> x7 >> x8 >> x9 >> f >> w >> iter;

    // if (itmaxitermaxitermaxitermaxiterer == maxiter)
    {
      row[0] = x1;
      row[1] = x2;
      row[2] = x3;
      row[3] = x4;
      row[4] = x5;
      row[5] = x6;
      row[6] = x7;
      row[7] = x8;
      row[8] = x9;

      /////////////////////////////////////////////////
      // get the 4-momenta
      Pm = P_minus_LAB(mj, mk, Mn, row);
      Pp = P_plus_LAB(mj, mk, Mn, row);
      // Pm = P_minus_Sframe(mj, mk, Mn, row);
      // Pp = P_plus_Sframe(mj, mk, Mn, row);


      if (!Momentum_contains_nans(Pp) && !Momentum_contains_nans(Pm) && !(is_nan(x1)) && !(is_nan(x9)) && !(is_nan(f*w * std_f / (std[iter-1]))) )
      {

        // if (counter2 >= NSTART)
        {
        //   /* code */
            outfile << x9 /*Enu*/  << " "
            << x1  /*Q^2*/ << " "
            << Pp[0]<< " "
            << Pp[1]<< " "
            << Pp[2]<< " "
            << Pp[3]<< " "
            << Pm[0]<< " "
            << Pm[1]<< " "
            << Pm[2]<< " "
            << Pm[3]<< " "
            << f*w * std_f / (std[iter-1])  << std::endl;

        }
        counter2+=1;
      }
      else{count_nans++;}

    }
  }

  ifile.close();
  outfile.close();
  
  std::cout<<"Number of printed events: "<<counter2<<std::endl;
  std::cout<<"Number of NaNed events: "<<count_nans<<std::endl;
  std::cout<<"Written to observables file."<< std::endl;

  ///////////////////////////////////
  // REMOVE integration samples file
  remove(samplesfile.c_str());
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////
void tridentMC::HEPevt_format(std::string eventsfile, std::string HEPevtfile, int NUMBER_OF_EVENTS)
{
  ///////////////////////////////////////////////////////
  // READING SAMPLES FILE A
  std::ifstream ifile(eventsfile.c_str());
  ifile.precision(17);
  std::cout<<"reading events file..."<<eventsfile.c_str()<<std::endl;
  long double E,Q2,pp0,pp1,pp2,pp3,pm0,pm1,pm2,pm3,wf;
  
  int line_count=0;
  long double sum_wf = 0;

  // Get the last number of iterations
  std::string line,dummyline;
  
  if (ifile.is_open())
  {
    std::getline(my_flux, dummyline);  // Skip header!
    std::getline(my_flux, dummyline);  // Skip header!

    while(std::getline(ifile,line))
    {
      std::istringstream iss(line);
  
      iss >> E >> Q2 >> pp0 >> pp1 >> pp2 >> pp3 >> pm0 >> pm1 >> pm2 >> pm3 >> wf;
      if (sum_wf<wf)
      {
        sum_wf = wf;
      }

      line_count++;
    }
  }
  else{
    std::cout << "ERROR! Not able to open events file." << std::endl;
  }

  std::cout<<"Max weight = "<<sum_wf<<std::endl;
  //////////////////////////////////////////////////////////

  // WRITING EVENT FILES AFTER ACCEPT/REJECT

  std::ofstream outfile(HEPevtfile.c_str());
  std::cout<<"opening HEPevtfile..."<<HEPevtfile.c_str()<<std::endl;
  outfile.precision(17);

  // Start reading file again
  ifile.clear();
  ifile.seekg(0);
  std::getline(my_flux, dummyline);  // Skip header!
  std::getline(my_flux, dummyline);  // Skip header!

  int i = 0;
  while (i < NUMBER_OF_EVENTS)
  {
      if(std::getline(ifile,line))
      {
        std::istringstream iss(line);
        if(!(iss >> E >> Q2 >> pp0 >> pp1 >> pp2 >> pp3 >> pm0 >> pm1 >> pm2 >> pm3 >> wf)){
          // std::cout<<E <<" "<< Q2 <<" "<< pp0 <<" "<< pp1 <<" "<< pp2 <<" "<< pp3 <<" "<< pm0 <<" "<< pm1 <<" "<< pm2 <<" "<< pm3 <<" "<< wf<<std::endl;
          continue;
        }

        if ( wf/sum_wf > UniformRand() )
        {
          // FIX ME -- INCLUDE THE STRUCK HADRON!
          outfile <<i<<" "<<3<<std::endl;
          outfile <<"1 "<<PDG_nu_inc<<" 0 0 0 0 "<<E<<" 0.0 0.0 "<<E<<" 0.0 0.0 0.0 0.0 0.0"<<std::endl;
          outfile <<"2 "<<PDG_lp<<" 0 0 0 0 "<<pp0<<" "<<pp1<<" "<<pp2<<" "<<pp3<<" "<<ml2<<" 0.0 0.0 0.0 0.0"<<std::endl;
          outfile <<"2 "<<PDG_lm<<" 0 0 0 0 "<<pm0<<" "<<pm1<<" "<<pm2<<" "<<pm3<<" "<<ml1<<" 0.0 0.0 0.0 0.0"<<std::endl;
          i++;
        }
      }
      else
      { 
        ifile.clear();
        ifile.seekg(0);
        std::getline(my_flux, dummyline);  // Skip header!
        std::getline(my_flux, dummyline);  // Skip header!

      }
  }
  
  ifile.close();
  outfile.close();
  
  ///////////////////////////////////
  // REMOVE original events file
  remove(eventsfile.c_str());
}

