Louis = function (t, Cb2, Cg2, Tb2, Tg2, bdervativeCi, bCi, s, d, Xp)
{
  # variance estimation
  Cg2[Cg2 < 0.00001] = 0.00001
  Tg2[Tg2 < 0.00001] = 0.00001
  k = length(Cb2)
  n = length(d)
  nods_k = length(Cg2)
  bi = t + (bCi %*% as.vector(Tg2)) * exp(Xp %*% Tb2) + (bCi %*% as.vector(Cg2)) * exp(Xp %*% Cb2)
  ci = t + (bCi %*% as.vector(Cg2)) * exp(Xp %*% Cb2)
  Eetai = (1-d)*(t + s)/bi + d * ((t + s) / ci) * ((1 - (ci/bi)^(t + s + 1 ))/(1 - (ci/bi)^(t + s)))
  Elogetai = (1-d)*(digamma(t + s) - log(bi)) +
    d * ( digamma(t + s) - (log(ci) - (ci/bi)^(t + s)*log(bi)) / (1 - (ci/bi)^(t + s)) )
  EVil = t(t(bdervativeCi)*as.vector(Cg2))/as.vector(bdervativeCi %*% as.vector(Cg2))
  EZi  = d * (bCi %*% Tg2) * exp(Xp %*% Tb2) * (t + s)/ci * bi^(t + s )/(bi^(t + s) - ci^(t + s))
  EZil = t(t(bCi) * Tg2)/as.vector(bCi %*% Tg2) * as.vector(EZi)

  Q2dtau2 =  - n * trigamma(t) + n / t
  Q2dbc2  =  - t(as.vector( bCi %*% as.vector(Cg2) * exp(Xp %*% Cb2) * Eetai) * Xp) %*% Xp    #2 * 2
  Q2dbt2  =  - t(as.vector( bCi %*% as.vector(Tg2) * exp(Xp %*% Tb2) * Eetai) * Xp) %*% Xp    #2 * 2

  Q2dbcrc =  - t(as.vector(exp(Xp %*% Cb2) * Eetai) * bCi ) %*% Xp
  Q2dbtrt =  - t(as.vector(exp(Xp %*% Tb2) * Eetai) * bCi ) %*% Xp

  Q2drcrc =  - diag(as.vector(colSums(s * EVil)/ Cg2 / Cg2))                                 #7 * 7
  Q2drtrt =  - diag(as.vector(colSums(EZil)/ Tg2 / Tg2))                                     #7 * 7
  #                       tau      bc      bt       rc       rt
  # Qdtheta2 = matrix(
  #           tau        Q2dtau2,  0,      0,       0,       0,
  #           bc           0,      Q2dbc2, 0,       Q2dbcrc, 0,
  #           bt           0,      0,      Q2dbt2,  0,       Q2dbtrt,
  #           rc           0,      Q2dbcrc,0,       Q2drcrc, 0,
  #           rt           0,      0,      Q2dbtrt, 0,       Q2drtrt )

  Qdtheta2 =
    rbind(           c(Q2dtau2,       matrix(0,1,2*nods_k+2*k)                              ),
                     cbind(matrix(0,k,1), Q2dbc2, matrix(0,k,k), t(Q2dbcrc),    matrix(0,k,nods_k)),
                     cbind(matrix(0,k,1+k),         Q2dbt2,        matrix(0,k,nods_k), t(Q2dbtrt)   ),
                     cbind(matrix(0,nods_k,1), Q2dbcrc,matrix(0,nods_k,k), Q2drcrc,       matrix(0,nods_k,nods_k)),
                     cbind(matrix(0,nods_k,1+k),         Q2dbtrt,       matrix(0,nods_k,nods_k), Q2drtrt      ))



  Var_Vil = EVil * (1 - EVil)
  Var_eta = (1-d) * (t + s)/(bi^2) +
    d * ( (t + s)*(t + s + 1)/(ci^2) * (1 - (ci/bi)^(t + s + 2))/(1 - (ci/bi)^(t + s)) -
            ((t + s) / ci * (1 - (ci/bi)^(t + s + 1 ))/(1 - (ci/bi)^(t + s)))^2 )
  Var_logeta = (1-d) * trigamma(t + s) +
    d * ( ( trigamma(t + s) + (digamma(t + s) - log(ci))^2 ) -
            (ci/bi)^(t + s) * ( trigamma(t + s) + (digamma(t + s) - log(bi))^2 ))/(1 - (ci/bi)^(t + s)) -
    d * ( digamma(t + s) - ( log(ci) - (ci/bi)^(t + s)* log(bi))/(1 - (ci/bi)^(t + s)))^2

  Cov_eta_logeta = (1-d)/bi + d * ((t + s) / ci * ( digamma(t + s + 1) - log(ci) -
                                                      (ci/bi)^(t + s + 1) *(digamma(t + s + 1) - log(bi))) / (1 - (ci/bi)^(t + s)) -
                                     (t + s) / ci * (1 - (ci / bi)^(t + s + 1)) / (1 - (ci / bi)^(t + s)) *
                                     (digamma(t + s) - (log(ci) - (ci/bi)^(t + s) * log(bi)) / (1 - (ci/bi)^(t + s) ) ))

  Var_Zi = d * ((bCi %*% as.vector(Tg2)) * exp(Xp %*% Tb2))^2 *
    (t + s + 1) * (t + s) / (ci^2) / (1 - (ci/bi)^(t + s)) + EZi - EZi^2

  Cov_Zij_Zi = t(t(bCi) * as.vector(Tg2))/as.vector(bCi %*% as.vector(Tg2)) * as.vector(Var_Zi)

  #Cov_Zi_etai = E(Zi etai) - E(Zi)E(etai)
  Cov_Zi_etai = d * (bCi %*% as.vector(Tg2)) * exp(Xp %*% Tb2) *
    (t + s + 1) * (t + s) / (ci^2)  / (1 - (ci/bi)^(t + s)) - EZi * Eetai
  Cov_Zij_etai = t(t(bCi) * Tg2)/as.vector(bCi %*% Tg2) * as.vector(Cov_Zi_etai)

  #Cov_Zi_logetai
  Cov_Zi_logetai = d * (bCi %*% as.vector(Tg2)) * exp(Xp %*% Tb2) *
    (t + s) * (digamma(t + s + 1) - log(ci)) / ci  / (1 - (ci/bi)^(t + s)) - EZi * Elogetai
  Cov_Zij_logetai = t(t(bCi) * Tg2)/as.vector(bCi %*% Tg2) * as.vector(Cov_Zi_logetai)

  ##########sum Cov_Zij_Zij'
  pij = t(t(bCi) * as.vector(Tg2))/as.vector(bCi %*% as.vector(Tg2))
  CovZZ = t(pij * as.vector(Var_Zi - EZi) ) %*% pij + diag(colSums(pij * as.vector(EZi)))
  p1ij = t(t(bdervativeCi) * as.vector(Cg2))/as.vector(bdervativeCi %*% as.vector(Cg2))
  CovVV = - t(s*p1ij) %*% p1ij + diag(colSums(s *p1ij))

  tau2 = sum(Var_eta + Var_logeta - 2 * Cov_eta_logeta)
  tau_bc = - t(bCi %*% as.vector(Cg2) * exp(Xp %*% Cb2) * (Cov_eta_logeta - Var_eta)) %*% Xp
  tau_bt = - t(bCi %*% as.vector(Tg2) * exp(Xp %*% Tb2) * (Cov_eta_logeta - Var_eta)) %*% Xp +
    as.vector(Cov_Zi_logetai - Cov_Zi_etai) %*% Xp
  tau_rc = - colSums(bCi * as.vector(exp(Xp %*% Cb2) * (Cov_eta_logeta - Var_eta)))
  tau_rt = - colSums(bCi * as.vector(exp(Xp %*% Tb2) * (Cov_eta_logeta - Var_eta))) +
    colSums(Cov_Zij_logetai - Cov_Zij_etai) / Tg2

  bc2 = t(as.vector((bCi %*% as.vector(Cg2) * exp(Xp %*% Cb2))^2 * Var_eta) * Xp) %*% Xp
  bt2 = t(as.vector((bCi %*% as.vector(Tg2) * exp(Xp %*% Tb2))^2 * Var_eta) * Xp) %*% Xp +
    t(Xp * as.vector(Var_Zi)) %*% Xp -
    2 * t(as.vector(bCi %*% as.vector(Tg2) * exp(Xp %*% Tb2) * Cov_Zi_etai) * Xp) %*% Xp
  bc_bt = t(as.vector(bCi %*% as.vector(Cg2) * exp(Xp %*% Cb2) * Var_eta) * Xp) %*%
    (as.vector(bCi %*% as.vector(Tg2) * exp(Xp %*% Tb2)) * Xp)-
    t(as.vector(bCi %*% as.vector(Cg2) * exp(Xp %*% Cb2) * Cov_Zi_etai) * Xp) %*% Xp

  ##########
  bc_rc = t(as.vector(bCi %*% as.vector(Cg2) * (exp(Xp %*% Cb2))^2 * Var_eta) * Xp) %*% bCi
  bc_rt = t(as.vector(bCi %*% as.vector(Cg2) * exp(Xp %*% Cb2) * exp(Xp %*% Tb2) * Var_eta) * Xp) %*% bCi -
    t(t(t(as.vector(bCi %*% as.vector(Cg2) * exp(Xp %*% Cb2)) * Xp) %*% Cov_Zij_etai) / Tg2)
  bt_rc = t(as.vector(bCi %*% as.vector(Tg2) * exp(Xp %*% Tb2) * exp(Xp %*% Cb2) * Var_eta) * Xp) %*% bCi -
    t(as.vector(exp(Xp %*% Cb2) * Cov_Zi_etai) * Xp) %*% bCi
  bt_rt = t((t(Cov_Zij_Zi)/as.vector(Tg2)) %*% Xp) +
    t(as.vector(bCi %*% as.vector(Tg2) * (exp(Xp %*% Tb2))^2 * Var_eta) * Xp) %*% bCi -
    t(t(t(as.vector(bCi %*% as.vector(Tg2) * exp(Xp %*% Tb2)) * Xp) %*% Cov_Zij_etai)/Tg2) -
    t(as.vector(exp(Xp %*% Tb2) * Cov_Zi_etai) * Xp) %*% bCi

  rc2 = t(CovVV / as.vector(Cg2)) / as.vector(Cg2) + t(bCi * as.vector((exp(Xp %*% Cb2))^2 * Var_eta)) %*% bCi
  rt2 = t(CovZZ / as.vector(Tg2)) / as.vector(Tg2) + t(bCi * as.vector((exp(Xp %*% Tb2))^2 * Var_eta)) %*% bCi -
    t(bCi * as.vector(exp(Xp %*% Tb2))) %*% t(t(Cov_Zij_etai)/Tg2) -
    t(t(bCi * as.vector(exp(Xp %*% Tb2))) %*% t(t(Cov_Zij_etai)/Tg2))
  rc_rt = t(bCi * as.vector(exp(Xp %*% Cb2))) %*% (bCi * as.vector(exp(Xp %*% Tb2) * Var_eta)) -
    t(bCi * as.vector(exp(Xp %*% Cb2))) %*% t(t(Cov_Zij_etai)/Tg2)

  #                             tau      bc           bt           rc            rt
  # Var_logL_theta = matrix(
  #               tau         c(tau2,    tau_bc  ,    tau_bt,      tau_rc,       tau_rt,
  #               bc            tau_bc,  bc2(2*2),    bc_bt(2*2),  bc_rc(2*2),   bc_rt(2*7),
  #               bt            tau_bt,  bc_bt(7*2),  bt2(7*7),    rc_bt(7*2),   bt_rt(7*7),
  #               rc            tau_rc,  bc_rc(2*2),  rc_bt(2*7),  rc2(2*2),     rc_rt(2*7),
  #               rt            tau_rt,  rt_bc(7*2),  rt_bt(7*7),  rt_rc(7*2),   rt2(7*7))    )


  Var_logL_theta = rbind(   c(  tau2    ,   tau_bc,   tau_bt,   tau_rc,  tau_rt),
                            cbind(t(tau_bc),      bc2,    bc_bt,    bc_rc,   bc_rt),
                            cbind(t(tau_bt),    bc_bt,      bt2,    bt_rc,   bt_rt),
                            cbind(matrix(tau_rc,nods_k,1), t(bc_rc),   t(bt_rc),   rc2,   rc_rt),
                            cbind(matrix(tau_rt,nods_k,1), t(bc_rt), t(bt_rt),    rc_rt,     rt2))

  return(- Qdtheta2 - Var_logL_theta )
}
