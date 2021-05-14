functions {
  vector SEIURD(
  vector state, //(S, E, IU, RU, ID, UD, RD, DD)
  real[] theta // Parameters: (betaU, betaD, alpha, gammaU, eta, nu, gammaD, deltaD)
  ) {
   vector[8] dstate_dt;
   real Nt;
   
   Nt = state[1] + state[2] + state[3] + state[4] + state[5] + state[6] + state[7];
   
   dstate_dt[1] = -(theta[1]*state[3] + theta[2]*state[5])*state[1]/Nt;
   dstate_dt[2] = (theta[1]*state[3] + theta[2]*state[5])*state[1]/Nt - theta[3]*state[2];
   dstate_dt[3] = theta[3]*state[2] - (theta[4] + theta[5])*state[3];
   dstate_dt[4] = theta[4]*state[3];
   dstate_dt[5] = theta[5]*state[3] - theta[6]*state[5];
   dstate_dt[6] = theta[6]*state[5] - (theta[7] + theta[8])*state[6];
   dstate_dt[7] = theta[7]*state[6];
   dstate_dt[8] = theta[8]*state[6];
   return dstate_dt;
  }


  // computing number of new cases from updated state
  real compute_new_cases(
    vector prev_state,
    vector new_state
  ) {
    real new_cases;
    new_cases = (new_state[5] - prev_state[5]) + (new_state[6] - prev_state[6]) + (new_state[7] - prev_state[7]) + (new_state[8] - prev_state[8]);
    return new_cases;
  }
  // computing number of new deaths from updated state
  real compute_new_deaths(
    vector prev_state,
    vector new_state
  ){
    real new_deaths;
    new_deaths = (new_state[8] - prev_state[8]);
    return new_deaths;
  }
}