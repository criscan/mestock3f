DATA_SECTION
 init_int ymin
 init_int ymax
 int ntime
 !!ntime=ymax-ymin+1;
 
 init_int nedades
 init_int minedad

 init_matrix mdatos(1,ntime,1,20)
 init_int ntallas
 init_vector Tallas(1,ntallas)

 init_int N_farr
 init_ivector nanos_farr(1,N_farr)
 init_matrix Carr(1,N_farr,1,ntallas)

 init_int N_fpal
 init_ivector nanos_fpal(1,N_fpal)
 init_matrix Cpal(1,N_fpal,1,ntallas)

 init_int N_fesp
 init_ivector nanos_fesp(1,N_fesp)
 init_matrix Cesp(1,N_fesp,1,ntallas)

 init_int N_fcru
 init_ivector nanos_fcru(1,N_fcru)
 init_matrix Ccru(1,N_fcru,1,ntallas)

 //!!ad_comm::change_datafile_name("mestock3F.ctl");
 init_number sigmaR
 init_vector dt(1,4)

 init_number L50ms
 init_number L95ms
 init_number lnaw
 init_number bw

 init_vector Par_bio(1,8)
 init_number cv_Loo
 init_number cv_k
 init_number cv_L0
 init_number cv_alfa
 init_number cv_beta
 init_number cv_M
 init_number cv_h
 init_number cv_R0
 init_int fase_Loo
 init_int fase_k
 init_int fase_L0
 init_int fase_alfa
 init_int fase_beta
 init_int fase_M
 init_int fase_h
 init_int fase_R0

  number log_Loo_prior
  number log_k_prior
  number log_L0_prior
  number log_alfa_prior
  number log_beta_prior
  number log_M_prior
  number log_h_prior
  number log_R0_prior
  
  !! log_Loo_prior = log(Par_bio(1));
  !! log_k_prior = log(Par_bio(2));
  !! log_L0_prior= log(Par_bio(3));
  !! log_alfa_prior = log(Par_bio(4)+1e-5);
  !! log_beta_prior = log(Par_bio(5)+1e-5);
  !! log_M_prior= log(Par_bio(6));
  !! log_h_prior = log(Par_bio(7));
  !! log_R0_prior = log(Par_bio(8));

 //Arrastre
 init_vector parsel1(1,3)
 init_vector cv_parsel1(1,3)
 init_int fase_parsel1_1
 init_int fase_parsel1_2
 init_int fase_parsel1_3

  number log_L50a_prior
  number log_sigma1a_prior
  number log_sigma2a_prior
 !! log_L50a_prior=log(parsel1(1));
 !! log_sigma1a_prior=log(parsel1(2));
 !! log_sigma2a_prior=log(parsel1(3));

 //Palangre
 init_vector parsel2(1,3)
 init_vector cv_parsel2(1,3)
 init_int fase_parsel2_1
 init_int fase_parsel2_2
 init_int fase_parsel2_3

  number log_L50p_prior
  number log_sigma1p_prior
  number log_sigma2p_prior
 !! log_L50p_prior=log(parsel2(1));
 !! log_sigma1p_prior=log(parsel2(2));
 !! log_sigma2p_prior=log(parsel2(3));

 //Espinel
 init_vector parsel3(1,3)
 init_vector cv_parsel3(1,3)
 init_int fase_parsel3_1
 init_int fase_parsel3_2
 init_int fase_parsel3_3

  number log_L50e_prior
  number log_sigma1e_prior
  number log_sigma2e_prior
 !! log_L50e_prior=log(parsel3(1));
 !! log_sigma1e_prior=log(parsel3(2));
 !! log_sigma2e_prior=log(parsel3(3));

 //Cruceros
 init_vector parsel4(1,3)
 init_vector cv_parsel4(1,3)
 init_int fase_parsel4_1
 init_int fase_parsel4_2
 init_int fase_parsel4_3

 number log_L50c_prior
 number log_sigma1c_prior
 number log_sigma2c_prior
 !! log_L50c_prior=log(parsel4(1));
 !! log_sigma1c_prior=log(parsel4(2));
 !! log_sigma2c_prior=log(parsel4(3));


 // bloques de selectividad
 init_int    nbloques_1
 init_vector ybloques1(1,nbloques_1)

 init_int    nbloques_2
 init_vector ybloques2(1,nbloques_2)

 init_int    nbloques_3
 init_vector ybloques3(1,nbloques_3)

 init_int    nbloques_4
 init_vector ybloques4(1,nbloques_4)

 // bloques de capturabilidad
 init_int    nqbloques1
 init_vector yqbloques1(1,nqbloques1)
 init_int    faseq_1

 init_int    nqbloques2
 init_vector yqbloques2(1,nqbloques2)
 init_int    faseq_2

 init_int    nqbloques3
 init_vector yqbloques3(1,nqbloques3)
 init_int    faseq_3

 init_int    nqbloques4
 init_vector yqbloques4(1,nqbloques4)
 init_int    faseq_4



 init_int    opt_Fa //_mortalidad_por_pesc
 init_int    opt_Fp //_mortalidad_por_pesc
 init_int    opt_Fe //_mortalidad_por_pesc

 init_int    opt_devRt

 init_int    npbr

 init_number  pbr_tar // PBR pivote

 init_vector pbr(1,npbr)
 init_int ntime_sim //Años_a_proyectar


INITIALIZATION_SECTION

  log_Loo        log_Loo_prior
  log_k          log_k_prior
  log_L0         log_L0_prior
  log_alfa       log_alfa_prior
  log_beta       log_beta_prior
  log_Rmed       log_R0_prior
  log_M          log_M_prior
  log_h          log_h_prior

  log_Fa         log_M_prior
  log_Fp         log_M_prior
  log_Fe         log_M_prior

  log_L50a        log_L50a_prior
  log_sigma1a     log_sigma1a_prior
  log_sigma2a     log_sigma2a_prior

  log_L50p        log_L50p_prior
  log_sigma1p     log_sigma1p_prior
  log_sigma2p     log_sigma2p_prior

  log_L50e        log_L50e_prior
  log_sigma1e     log_sigma1e_prior
  log_sigma2e     log_sigma2e_prior

  log_L50c        log_L50c_prior
  log_sigma1c     log_sigma1c_prior
  log_sigma2c     log_sigma2c_prior



PARAMETER_SECTION

// parametros reclutamientos y mortalidades)
 init_number log_Rmed(fase_R0)
 init_vector log_desv_Rt(1,ntime,opt_devRt)

 init_vector log_Fa(1,ntime,opt_Fa) // log  mortalidad por pesca por flota
 init_vector log_Fp(1,ntime,opt_Fp) // log  mortalidad por pesca por flota
 init_vector log_Fe(1,ntime,opt_Fe) // log  mortalidad por pesca por flota

// Parámetros biológicos
 init_number log_Loo(fase_Loo)
 init_number log_k(fase_k)
 init_number log_L0(fase_L0)
 init_number log_alfa(fase_alfa)
 init_number log_beta(fase_beta)
 init_number log_M(fase_M)
 init_number log_h(fase_h)

// selectividad paramétrica a la edad
 init_vector log_L50a(1,nbloques_1,fase_parsel1_1)
 init_vector log_sigma1a(1,nbloques_1,fase_parsel1_2)
 init_vector log_sigma2a(1,nbloques_1,fase_parsel1_3)

 init_vector log_L50p(1,nbloques_2,fase_parsel2_1)
 init_vector log_sigma1p(1,nbloques_2,fase_parsel2_2)
 init_vector log_sigma2p(1,nbloques_2,fase_parsel2_3)

 init_vector log_L50e(1,nbloques_3,fase_parsel3_1)
 init_vector log_sigma1e(1,nbloques_3,fase_parsel3_2)
 init_vector log_sigma2e(1,nbloques_3,fase_parsel3_3)

 init_vector log_L50c(1,nbloques_4,fase_parsel4_1)
 init_vector log_sigma1c(1,nbloques_4,fase_parsel4_2)
 init_vector log_sigma2c(1,nbloques_4,fase_parsel4_3)

// capturabilidades
 init_vector log_q1(1,nqbloques1,faseq_1)
 init_vector log_q2(1,nqbloques2,faseq_2)
 init_vector log_q3(1,nqbloques3,faseq_3)
 init_vector log_q4(1,nqbloques4,faseq_4)


//---------------------------------------------------------------------------------
//Defino las variables de estado

 matrix Sel_tot(1,ntime,1,nedades)
 matrix S1(1,nbloques_1,1,nedades)
 matrix S2(1,nbloques_2,1,nedades)
 matrix S3(1,nbloques_3,1,nedades)
 matrix S4(1,nbloques_4,1,nedades)
 matrix Sela(1,ntime,1,nedades)
 matrix Selp(1,ntime,1,nedades)
 matrix Sele(1,ntime,1,nedades)
 matrix Sel_cru(1,ntime,1,nedades)
 matrix Y(1,3,1,ntime)
 matrix cvY(1,3,1,ntime)
 matrix CPUE(1,3,1,ntime)
 matrix cvCPUE(1,3,1,ntime)
 matrix nm(1,4,1,ntime)
 vector Bacu(1,ntime)
 vector cvBacu(1,ntime)

 vector edades(1,nedades)
 vector yrs(1,ntime)
 vector dt_cru(1,ntime)
 vector Unos_edades(1,nedades) 
 vector Unos_year(1,ntime)
 
 vector mu_edad(1,nedades)
 vector sigma_edad(1,nedades)
 vector msex(1,ntallas)
 vector Wmed(1,ntallas)

 number M
 number h
 number alfa_sr
 number beta_sr
 number SSBo
 number Linf
 number k
 number rango
 number L0

 matrix Farr(1,ntime,1,nedades) 
 matrix Fpal(1,ntime,1,nedades) 
 matrix Fesp(1,ntime,1,nedades)
 matrix Z(1,ntime,1,nedades) 
 matrix S(1,ntime,1,nedades) 
 vector No(1,nedades)

 matrix N(1,ntime,1,nedades)
 vector Rpred(1,ntime)
 matrix ND(1,ntime,1,ntallas)

 matrix NVarr(1,ntime,1,ntallas)
 matrix NVpal(1,ntime,1,ntallas)
 matrix NVes(1,ntime,1,ntallas)
 matrix NVacu(1,ntime,1,nedades)
 vector BVarr(1,ntime)
 vector BVpal(1,ntime)
 vector BVes(1,ntime)
 vector BVcru(1,ntime)

 matrix pred_Carr(1,ntime,1,ntallas)
 matrix pred_Cpal(1,ntime,1,ntallas)
 matrix pred_Cesp(1,ntime,1,ntallas)

 vector pred_Yarr(1,ntime)
 vector pred_Ypal(1,ntime)
 vector pred_Yesp(1,ntime)

 matrix ppred_arr(1,N_farr,1,ntallas)
 matrix pobs_arr(1,N_farr,1,ntallas)
 matrix ppred_pal(1,N_fpal,1,ntallas)
 matrix pobs_pal(1,N_fpal,1,ntallas)
 matrix ppred_esp(1,N_fesp,1,ntallas)
 matrix pobs_esp(1,N_fesp,1,ntallas)
 matrix ppred_cru(1,N_fcru,1,ntallas)
 matrix pobs_cru(1,N_fcru,1,ntallas)

 matrix Prob_talla(1,nedades,1,ntallas)


 vector pred_Lmed_arr(1,N_farr)
 vector obs_Lmed_arr(1,N_farr)
 vector pred_Lmed_pal(1,N_fpal)
 vector obs_Lmed_pal(1,N_fpal)
 vector pred_Lmed_esp(1,N_fesp)
 vector obs_Lmed_esp(1,N_fesp)
 vector pred_Lmed_cru(1,N_fcru)
 vector obs_Lmed_cru(1,N_fcru)

 vector pred_CPUE_arr(1,ntime)
 vector pred_CPUE_pal(1,ntime)
 vector pred_CPUE_es(1,ntime)
 vector pred_Bcru(1,ntime)

 number suma1
 number suma2
 number suma3
 number suma4
 number suma5
 number suma6
 number suma7
 number suma8

 vector likeval(1,11)
 vector penalty(1,20)

 number eps
 number maxF


 vector Np(1,nedades)
 vector Zpbr(1,nedades)
 vector Fpbr(1,nedades)
 vector Sp(1,nedades)

 matrix Bp(1,npbr,1,ntime_sim)
 vector CTPp(1,nedades)
 matrix Yp(1,npbr,1,ntime_sim)
 matrix Fproy(1,npbr,1,ntime_sim)
 number BDp
 number Npplus
 number Bref
 number plus

 sdreport_vector Ftot(1,ntime) 
 sdreport_vector Reclutas(1,ntime)
 sdreport_vector BD(1,ntime)
 sdreport_vector SPR(1,ntime)
 sdreport_vector Redstock(1,npbr)

 objective_function_value f

PRELIMINARY_CALCS_SECTION

 eps=1.e-10;
 yrs=column(mdatos,1);

  Y(1)=column(mdatos,2);
  Y(2)=column(mdatos,4);
  Y(3)=column(mdatos,6);

  cvY(1)=column(mdatos,3);
  cvY(2)=column(mdatos,5);
  cvY(3)=column(mdatos,7);

  CPUE(1)=column(mdatos,8);
  CPUE(2)=column(mdatos,10);
  CPUE(3)=column(mdatos,12);
  Bacu=column(mdatos,14);
  
  cvCPUE(1)=column(mdatos,9);
  cvCPUE(2)=column(mdatos,11);
  cvCPUE(3)=column(mdatos,13);
  cvBacu=column(mdatos,15);

  nm(1)=column(mdatos,16);
  nm(2)=column(mdatos,17);
  nm(3)=column(mdatos,18);
  nm(4)=column(mdatos,19);

  dt_cru=column(mdatos,20);


  Unos_edades=1;// lo uso en operaciones matriciales con edades
  Unos_year=1;// lo uso en operaciones matriciales con el año

  Wmed=exp(lnaw)*pow(Tallas,bw);
  rango=L95ms-L50ms;
  msex=1/(1+exp(-log(19)*(Tallas-L50ms)/rango));

 
PROCEDURE_SECTION

  Eval_prob_talla_edad();
  Eval_selectividad();
  Eval_mortalidades();
  Eval_abundancia();
  Eval_biomasas();
  Eval_capturas_predichas();
  Eval_indices();
  Eval_logverosim();
  Eval_funcion_objetivo();
  if(last_phase()){Eval_Proyeccion();}

FUNCTION Eval_prob_talla_edad


 Linf=exp(log_Loo);
 k=exp(log_k);
 L0=exp(log_L0);

 int i, j;

// genero una clave edad-talla para otros calculos. Se modela desde L(1)
 mu_edad(1)=L0;
 for (i=2;i<=nedades;i++)
  {
  mu_edad(i)=Linf*(1-exp(-k))+exp(-k)*mu_edad(i-1);
  }

  sigma_edad=exp(log_alfa)+exp(log_beta)*mu_edad;

  Prob_talla = ALK( mu_edad, sigma_edad, Tallas);

FUNCTION Eval_selectividad

 int i,j;

 // Arrastre...................

 for (j=1;j<=nbloques_1;j++){

  S1(j)=exp(-0.5*square(mu_edad-exp(log_L50a(j)))/square(exp(log_sigma1a(j))));

    for (i=1;i<=nedades;i++){
      if(mu_edad(i)>=exp(log_L50a(j))){
          S1(j,i)= exp(-0.5*square(mu_edad(i)-exp(log_L50a(j)))/square(exp(log_sigma2a(j))));
      }
 }}

   for (i=1;i<=ntime;i++){
      for (j=1;j<=nbloques_1;j++){
              if (yrs(i)>=ybloques1(j)){
                Sela(i)=S1(j);}
       }
   }

 


 // Palangre...................

 for (j=1;j<=nbloques_2;j++){

 S2(j)=exp(-0.5*square(mu_edad-exp(log_L50p(j)))/square(exp(log_sigma1p(j))));
    for (i=1;i<=nedades;i++){
      if(mu_edad(i)>=exp(log_L50p(j))){
     S2(j,i)= exp(-0.5*square(mu_edad(i)-exp(log_L50p(j)))/square(exp(log_sigma2p(j))));

      }
 }}
   for (i=1;i<=ntime;i++){
      for (j=1;j<=nbloques_2;j++){
              if (yrs(i)>=ybloques2(j)){
                Selp(i)=S2(j);}
       }
   }

 // Espinel...................

 for (j=1;j<=nbloques_3;j++){

 S3(j)=exp(-0.5*square(mu_edad-exp(log_L50e(j)))/square(exp(log_sigma1e(j))));
    for (i=1;i<=nedades;i++){
      if(mu_edad(i)>=exp(log_L50e(j))){
      S3(j,i)= exp(-0.5*square(mu_edad(i)-exp(log_L50e(j)))/square(exp(log_sigma2e(j))));
      }
 }}
   for (i=1;i<=ntime;i++){
      for (j=1;j<=nbloques_3;j++){
              if (yrs(i)>=ybloques3(j)){
                Sele(i)=S3(j);}
       }
   }

 // Crucero ...................
 for (j=1;j<=nbloques_4;j++){

  S4(j)=exp(-0.5*square(mu_edad-exp(log_L50c(j)))/square(exp(log_sigma1c(j))));
    for (i=1;i<=nedades;i++){
      if(mu_edad(i)>=exp(log_L50c(j))){
      S4(j,i)= exp(-0.5*square(mu_edad(i)-exp(log_L50c(j)))/square(exp(log_sigma2c(j))));
      }
 }}


   for (i=1;i<=ntime;i++){
      for (j=1;j<=nbloques_4;j++){
              if (yrs(i)>=ybloques4(j)){
                Sel_cru(i)=S4(j);}
       }
   }



FUNCTION Eval_mortalidades

 M=exp(log_M);

 for (int i=1;i<=ntime;i++){
 Farr(i)=Sela(i)*mfexp(log_Fa(i));
 Fpal(i)=Selp(i)*mfexp(log_Fp(i));
 Fesp(i)=Sele(i)*mfexp(log_Fe(i));
 }

 Z=Farr+Fpal+Fesp+M;
 S=mfexp(-1.0*Z);

 for (int i=1;i<=ntime;i++){
 Ftot(i)=max(Farr(i)+Fpal(i)+Fesp(i));
 }
 
 // Calculo de la selectividad total
 for (int i=1;i<=ntime;i++){
 Sel_tot(i)=(Farr(i)+Fpal(i)+Fesp(i))/(max(Farr(i)+Fpal(i)+Fesp(i)));
 }


FUNCTION Eval_abundancia
 int i, j;

// genero la poblacion en equilibrio virginal de LP;
  No(1)=mfexp(log_Rmed);//-0.5*pow(sigmaR,2));
  for (int j=2;j<=nedades;j++){
  No(j)=No(j-1)*mfexp(-M); // 
  }
  No(nedades)=No(nedades)/(1-mfexp(-M)); // plus grup

  h=mfexp(log_h);
  SSBo=sum(elem_prod(No*mfexp(-dt(1)*M)*Prob_talla,elem_prod(Wmed,msex)));
  alfa_sr=4*h*mfexp(log_Rmed)/(5*h-1);//
  beta_sr=(1-h)*SSBo/(5*h-1);// Reclutamiento

 // -----------------primer año 
 // Condición inicial año 1 en equilibrio en torno a Z y R
    N(1,1)=mfexp(log_Rmed+log_desv_Rt(1));
    for (int j=2;j<=nedades;j++)
    {N(1,j)=N(1,j-1)*mfexp(-1.*Z(1,j-1));}
     N(1,nedades)=N(1,nedades)/(1-mfexp(-1.*Z(1,nedades)));


  ND(1)=elem_prod(elem_prod(N(1),mfexp(-dt(1)*Z(1)))*Prob_talla,msex);
  BD(1)=sum(elem_prod(Wmed,ND(1)));
  
  Rpred=mfexp(log_Rmed);
  Reclutas(1)=N(1,1);


// --------------------dinamica anual
 
  for (int y=1;y<ntime;y++){

        if(y>minedad){
        Rpred(y+1)=(alfa_sr*BD(y-minedad)/(beta_sr+BD(y-minedad)));
        }
        Reclutas(y+1)=Rpred(y+1)*mfexp(log_desv_Rt(y+1));

        N(y+1,1)=Reclutas(y+1);  // 
        N(y+1)(2,nedades)=++elem_prod(N(y)(1,nedades-1),S(y)(1,nedades-1));
        N(y+1,nedades)+=N(y,nedades)*S(y,nedades);// grupo plus

        ND(y+1)=elem_prod(elem_prod(N(y+1),mfexp(-dt(1)*Z(y+1)))*Prob_talla,msex);
        BD(y+1)=sum(elem_prod(Wmed,ND(y+1)));
   }
     



FUNCTION Eval_biomasas

 NVarr=elem_prod(elem_prod(N,mfexp(-dt(2)*(Z))),Sela)*Prob_talla;
 NVpal=elem_prod(elem_prod(N,mfexp(-dt(3)*(Z))),Selp)*Prob_talla;
 NVes=elem_prod(elem_prod(N,mfexp(-dt(4)*(Z))),Sele)*Prob_talla;


 for (int i=1;i<=ntime;i++){

  NVacu(i)=elem_prod(elem_prod(N(i),mfexp(-dt_cru(i)*Z(i))),Sel_cru(i));

 }

 BVarr=NVarr*Wmed;
 BVpal=NVpal*Wmed;
 BVes=NVes*Wmed;
 BVcru=(NVacu*Prob_talla)*Wmed;

 SPR=BD/SSBo;



FUNCTION Eval_capturas_predichas


// matrices de capturas predichas por talla y año
 pred_Carr=elem_prod(elem_div(Farr,Z),elem_prod(1-S,N))*Prob_talla;
 pred_Cpal=elem_prod(elem_div(Fpal,Z),elem_prod(1-S,N))*Prob_talla;
 pred_Cesp=elem_prod(elem_div(Fesp,Z),elem_prod(1-S,N))*Prob_talla;


// vectores de desembarques predichos por año
 pred_Yarr=pred_Carr*Wmed;
 pred_Ypal=pred_Cpal*Wmed;
 pred_Yesp=pred_Cesp*Wmed;


 for (int i=1;i<=N_farr;i++){
 ppred_arr(i)=pred_Carr(nanos_farr(i)-ymin+1)/sum(pred_Carr(nanos_farr(i)-ymin+1)+eps);
 pobs_arr(i)=Carr(i)/sum(Carr(i)+eps);
 }

 for (int i=1;i<=N_fpal;i++){
 ppred_pal(i)=pred_Cpal(nanos_fpal(i)-ymin+1)/sum(pred_Cpal(nanos_fpal(i)-ymin+1)+eps);
 pobs_pal(i)=Cpal(i)/sum(Cpal(i)+eps);
 }

 for (int i=1;i<=N_fesp;i++){
 ppred_esp(i)=pred_Cesp(nanos_fesp(i)-ymin+1)/sum(pred_Cesp(nanos_fesp(i)-ymin+1)+eps);
 pobs_esp(i)=Cesp(i)/sum(Cesp(i)+eps);
 }

 for (int i=1;i<=N_fcru;i++){
 ppred_cru(i)=NVacu(nanos_fcru(i)-ymin+1)*Prob_talla/sum(NVacu(nanos_fcru(i)-ymin+1)*Prob_talla+eps);
 pobs_cru(i)=Ccru(i)/sum(Ccru(i)+eps);
 }


 // Edades promedio por flota

 pred_Lmed_arr=Tallas*trans(ppred_arr);
 obs_Lmed_arr=Tallas*trans(pobs_arr);


 pred_Lmed_pal=Tallas*trans(ppred_pal);
 obs_Lmed_pal=Tallas*trans(pobs_pal);

 pred_Lmed_esp=Tallas*trans(ppred_esp);
 obs_Lmed_esp=Tallas*trans(pobs_esp);

 pred_Lmed_cru=Tallas*trans(ppred_cru);
 obs_Lmed_cru=Tallas*trans(pobs_cru);


FUNCTION Eval_indices

   for (int i=1;i<=ntime;i++){
      for (int j=1;j<=nqbloques1;j++){
              if (yrs(i)>=yqbloques1(j)){ //   -ymin+1){
                 pred_CPUE_arr(i)=exp(log_q1(j))*BVarr(i);}
       }}

   for (int i=1;i<=ntime;i++){
      for (int j=1;j<=nqbloques2;j++){
              if (yrs(i)>=yqbloques2(j)){ //-ymin+1){
                 pred_CPUE_pal(i)=exp(log_q2(j))*BVpal(i);}
       }}

   for (int i=1;i<=ntime;i++){
      for (int j=1;j<=nqbloques3;j++){
              if (yrs(i)>=yqbloques3(j)){ //-ymin+1){
                 pred_CPUE_es(i)=exp(log_q3(j))*BVes(i);}
       }}

   for (int i=1;i<=ntime;i++){
      for (int j=1;j<=nqbloques4;j++){
              if (yrs(i)>=yqbloques4(j)){ // -ymin+1){
                 pred_Bcru(i)=exp(log_q4(j))*BVcru(i);}
       }}




FUNCTION Eval_logverosim
 
 suma1=0; suma2=0; suma3=0; suma4=0; 

 for (int y=1;y<=ntime;y++)
 {
  if (CPUE(1,y)>0){
    suma1+=square(log(CPUE(1,y)/pred_CPUE_arr(y))*1/cvCPUE(1,y));}

  if (CPUE(2,y)>0){
    suma2+=square(log(CPUE(2,y)/pred_CPUE_pal(y))*1/cvCPUE(2,y));}

  if (CPUE(3,y)>0){
    suma3+=square(log(CPUE(3,y)/pred_CPUE_es(y))*1/cvCPUE(3,y));}

  if (Bacu(y)>0){
    suma4+=square(log(Bacu(y)/pred_Bcru(y))*1/cvBacu(y));}
 }



FUNCTION Eval_funcion_objetivo


 suma5=0;suma6=0;suma7=0;suma8=0; 

 likeval(1)=0.5*suma1;//CPUE_arr
 likeval(2)=0.5*suma2;//CPUE_pal
 likeval(3)=0.5*suma3;//CPUE_esp
 likeval(4)=0.5*suma4;//Crucero

 likeval(5)=0.5*norm2(elem_div(log(elem_div(Y(1)+eps,pred_Yarr)),cvY(1)));// Arrastre
 likeval(6)=0.5*norm2(elem_div(log(elem_div(Y(2)+eps,pred_Ypal)),cvY(2)));// Palangre
 likeval(7)=0.5*norm2(elem_div(log(elem_div(Y(3)+eps,pred_Yesp)),cvY(3)));// Palangre



 for (int i=1;i<=N_farr;i++){
 //suma5+=-nm(1,nanos_farr(i)-min(nanos_farr)+1)*sum(elem_prod(pobs_arr(i),log(ppred_arr(i)+eps)));}
 suma5+=-nm(1,nanos_farr(i)-ymin+1)*sum(elem_prod(pobs_arr(i),log(ppred_arr(i)+eps)));}

 for (int i=1;i<=N_fpal;i++){
  //suma6+=-nm(2,nanos_fpal(i)-min(nanos_fpal)+1)*sum(elem_prod(pobs_pal(i),log(ppred_pal(i)+eps)));}
 suma6+=-nm(2,nanos_fpal(i)-ymin+1)*sum(elem_prod(pobs_pal(i),log(ppred_pal(i)+eps)));}

 for (int i=1;i<=N_fesp;i++){
 //suma7+=-nm(3,nanos_fesp(i)-min(nanos_fesp)+1)*sum(elem_prod(pobs_esp(i),log(ppred_esp(i)+eps)));}
 suma7+=-nm(3,nanos_fesp(i)-ymin+1)*sum(elem_prod(pobs_esp(i),log(ppred_esp(i)+eps)));}

 for (int i=1;i<=N_fcru;i++){
 //suma8+=-nm(4,nanos_fcru(i)-min(nanos_fcru)+1)*sum(elem_prod(pobs_cru(i),log(ppred_cru(i)+eps)));}
  suma8+=-nm(4,nanos_fcru(i)-ymin+1)*sum(elem_prod(pobs_cru(i),log(ppred_cru(i)+eps)));}


 likeval(8)=suma5;//
 likeval(9)=suma6;//
 likeval(10)=suma7;//
 likeval(11)=suma8;//


 //Priors de parámetros y procesos
 
 penalty(1)=1./(2*square(sigmaR))*(norm2(log_desv_Rt));

 penalty(2)=1./(2*square(cv_Loo))*square(log_Loo-log_Loo_prior);
 penalty(3)=1./(2*square(cv_k))*square(log_k-log_k_prior);
 penalty(4)=1./(2*square(cv_L0))*square(log_L0-log_L0_prior);
 penalty(5)=1./(2*square(cv_alfa))*square(log_alfa-log_alfa_prior);
 penalty(6)=1./(2*square(cv_beta))*square(log_beta-log_beta_prior);
 penalty(7)=1./(2*square(cv_M))*square(log_M-log_M_prior);
 penalty(8)=1./(2*square(cv_h))*square(log_h-log_h_prior);


 penalty(9)=1./(2*square(cv_parsel1(1)))*(norm2(log_L50a_prior-log_L50a));
 penalty(10)=1./(2*square(cv_parsel1(2)))*(norm2(log_sigma1a_prior-log_sigma1a));
 penalty(11)=1./(2*square(cv_parsel1(3)))*(norm2(log_sigma2a_prior-log_sigma2a));

 penalty(12)=1./(2*square(cv_parsel2(1)))*(norm2(log_L50p_prior-log_L50p));
 penalty(13)=1./(2*square(cv_parsel2(2)))*(norm2(log_sigma1p_prior-log_sigma1p));
 penalty(14)=1./(2*square(cv_parsel2(3)))*(norm2(log_sigma2p_prior-log_sigma2p));

 penalty(15)=1./(2*square(cv_parsel3(1)))*(norm2(log_L50e_prior-log_L50e));
 penalty(16)=1./(2*square(cv_parsel3(2)))*(norm2(log_sigma1e_prior-log_sigma1e));
 penalty(17)=1./(2*square(cv_parsel3(3)))*(norm2(log_sigma2e_prior-log_sigma2e));

 penalty(18)=1./(2*square(cv_parsel4(1)))*(norm2(log_L50c_prior-log_L50c));
 penalty(19)=1./(2*square(cv_parsel4(2)))*(norm2(log_sigma1c_prior-log_sigma1c));
 penalty(20)=1./(2*square(cv_parsel4(3)))*(norm2(log_sigma2c_prior-log_sigma2c));
 
 f=sum(likeval)+sum(penalty);


FUNCTION  Eval_Proyeccion
//-----------------------------------------------------------------

  for (int i=1;i<=npbr;i++){ // ciclo de PBR

  Np=N(ntime);
  Sp=S(ntime);
  Bref=BD(ntime);

   for (int j=1;j<=ntime_sim;j++){ // ciclo de años

    plus=Np(nedades)*Sp(nedades);
    if(j<=minedad){
    Np(1)=alfa_sr*BD(ntime-minedad+j)/(beta_sr+BD(ntime-minedad+j));
    } // cuando j<minedad

    if(j>minedad){
    Np(1)=alfa_sr*Bp(i,j-minedad)/(beta_sr+Bp(i,j-minedad));
    } // cuando j>minedad

    Np(2,nedades)=++elem_prod(Np(1,nedades-1),Sp(1,nedades-1));
    Np(nedades)+=plus;

    Fpbr=(Farr(ntime)+Fpal(ntime)+Fesp(ntime))*pbr(i);
    


    if(Bref/SSBo<=pbr_tar){ // regla de control

     Fpbr=Fpbr*Bref/SSBo;
    }

    Zpbr=Fpbr+M;

    Bp(i,j)=sum(elem_prod(elem_prod(Np,mfexp(-dt(1)*Zpbr))*Prob_talla,elem_prod(msex,Wmed)));
    CTPp=elem_prod(elem_div(Fpbr,Zpbr),elem_prod(1.-exp(-1.*Zpbr),Np));
    Yp(i,j)=sum(elem_prod(CTPp*Prob_talla,Wmed));
    Fproy(i,j)=max(Fpbr);

    Sp=exp(-1.*Zpbr);
    Bref=Bp(i,j);

  }
   Redstock(i)=Bref/SSBo;
  }



REPORT_SECTION

  report << "Likelihood components" << endl;
 report << "_____________________" << endl;
 report << "CPUE_F1   :" <<" "<<likeval(1)<< endl;
 report << "CPUE_F2   :" <<" "<<likeval(2)<< endl;
 report << "CPUE_F3   :" <<" "<<likeval(3)<< endl;
 report << "Biom_camp :" <<" "<<likeval(4)<< endl;
 report << "Y_F1      :" <<" "<<likeval(5)<< endl;
 report << "Y_F2      :" <<" "<<likeval(6)<< endl;
 report << "Y_F3      :" <<" "<<likeval(7)<< endl;
 report << "frec_F1   :" <<" "<<likeval(8)<< endl;
 report << "frec_F2   :" <<" "<<likeval(9)<< endl;
 report << "frec_F3   :" <<" "<<likeval(10)<< endl;
 report << "frec_camp :" <<" "<<likeval(11)<< endl;
 report << "_____________________" << endl;
 report << "LL total  :" <<" "<<sum(likeval)<< endl;
 
 report <<" "<<endl;
 report << "Priors y penalties" << endl;
 report << "_____________________" << endl;
 report << "devR      :" <<" "<<penalty(1)<< endl;
 report << "Loo       :" <<" "<<penalty(2)<< endl;
 report << "k         :" <<" "<<penalty(3)<< endl;
 report << "Lo        :" <<" "<<penalty(4)<< endl;
 report << "alfa      :" <<" "<<penalty(5)<< endl;
 report << "beta      :" <<" "<<penalty(6)<< endl;
 report << "M         :" <<" "<<penalty(7)<< endl;
 report << "h         :" <<" "<<penalty(8)<< endl;
 report << "L50_F1    :" <<" "<<penalty(9)<< endl;
 report << "s1_F1     :" <<" "<<penalty(10)<< endl;
 report << "s2_F1     :" <<" "<<penalty(11)<< endl;
 report << "L50_F2    :" <<" "<<penalty(12)<< endl;
 report << "s1_F2     :" <<" "<<penalty(13)<< endl;
 report << "s2_F2     :" <<" "<<penalty(14)<< endl;
 report << "L50_F3    :" <<" "<<penalty(15)<< endl;
 report << "s1_F3     :" <<" "<<penalty(16)<< endl;
 report << "s2_F3     :" <<" "<<penalty(17)<< endl;
 report << "L50_camp  :" <<" "<<penalty(18)<< endl;
 report << "s1_camp   :" <<" "<<penalty(19)<< endl;
 report << "s2_camp   :" <<" "<<penalty(20)<< endl;
 report << "_____________________" << endl;
 report << "Pp total  :" <<" "<<sum(penalty)<< endl;
 report << "_____________________" << endl;
 report << "Obj Function  :" <<" "<<sum(penalty)+sum(likeval)<< endl;



 
RUNTIME_SECTION
  convergence_criteria 1.e-3,1.e-4,1e-4,1e-5
  maximum_function_evaluations 1000,2000,3000,3500


TOP_OF_MAIN_SECTION
 time(&start);
 arrmblsize = 90000000; 
 gradient_structure::set_GRADSTACK_BUFFER_SIZE(1.e7); 
 gradient_structure::set_CMPDIF_BUFFER_SIZE(1.e7); 
 gradient_structure::set_MAX_NVAR_OFFSET(5000); 
 gradient_structure::set_NUM_DEPENDENT_VARIABLES(5000); 



GLOBALS_SECTION
 #include <admodel.h>
 #include <stdio.h>
 #include <time.h>
 time_t start,finish;
 long hour,minute,second;
 double elapsed_time;
 ofstream mcmc_report("mcmc.csv");


FINAL_SECTION
 time(&finish);
 elapsed_time=difftime(finish,start);
 hour=long(elapsed_time)/3600;
 minute=long(elapsed_time)%3600/60;
 second=(long(elapsed_time)%3600)%60;
 cout<<endl<<endl<<"*********************************************"<<endl;
 cout<<"--Start time:  "<<ctime(&start)<<endl;
 cout<<"--Finish time: "<<ctime(&finish)<<endl;
 cout<<"--Runtime: ";
 cout<<hour<<" hours, "<<minute<<" minutes, "<<second<<" seconds"<<endl;
 cout<<"*********************************************"<<endl;



 ofstream print_R("For_R.rep");

 print_R << "Yrs" << endl;
 print_R << yrs << endl;
 print_R << "Y_obs_pred_F1" << endl;
 print_R << Y(1) << endl;
 print_R << pred_Yarr << endl;
 print_R << "Y_obs_pred_F2" << endl;
 print_R << Y(2) << endl;
 print_R << pred_Ypal << endl;
 print_R << "Y_obs_pred_F3" << endl;
 print_R << Y(3) << endl;
 print_R << pred_Yesp << endl;
 print_R << "CPUE_obs_pred_F1" << endl;
 print_R << CPUE(1) << endl;
 print_R << pred_CPUE_arr << endl;
 print_R << "CPUE_obs_pred_F2" << endl;
 print_R << CPUE(2) << endl;
 print_R << pred_CPUE_pal << endl;
 print_R << "CPUE_obs_pred_F3" << endl;
 print_R << CPUE(3) << endl;
 print_R << pred_CPUE_es << endl;
 print_R << "Bacu_obs_pred" << endl;
 print_R << Bacu << endl;
 print_R << pred_Bcru << endl;
 print_R << "Tallas" << endl;
 print_R << Tallas << endl;
 print_R << "Msex_edad" << endl;
 print_R << 1/(1+mfexp(-log(19)*(mu_edad-L50ms)/rango)) << endl;
 print_R << "Peso_edad" << endl;
 print_R << mfexp(lnaw)*pow(mu_edad,bw) << endl;
 print_R << "Lmed_edad" << endl;
 print_R << mu_edad << endl;
 print_R << "sigma_edad" << endl;
 print_R << sigma_edad << endl;

 print_R << "Prob_talla" << endl;
 print_R << Prob_talla << endl;

 print_R << "No" << endl;
 print_R << No << endl;



 print_R << "Sel_F1" << endl;
 print_R << Sela << endl;
 print_R << "Sel_F2" << endl;
 print_R << Selp << endl;
 print_R << "Sel_F3" << endl;
 print_R << Sele << endl;
 print_R << "Sel_cru" << endl;
 print_R << Sel_cru << endl;
 print_R << "Sel_tot" << endl;
 print_R << Sel_tot << endl;
 print_R << "Lm_F1" << endl;
 print_R<< nanos_farr<< endl;
 print_R << obs_Lmed_arr << endl;
 print_R << pred_Lmed_arr << endl;
 print_R << "Lm_F2" << endl;
 print_R<< nanos_fpal<< endl;
 print_R << obs_Lmed_pal << endl;
 print_R << pred_Lmed_pal << endl;
 print_R << "Lm_F3" << endl;
 print_R<< nanos_fesp<< endl;
 print_R << obs_Lmed_esp << endl;
 print_R << pred_Lmed_esp << endl;
 print_R << "Lm_cru" << endl;
 print_R<< nanos_fcru<< endl;
 print_R << obs_Lmed_cru << endl;
 print_R << pred_Lmed_cru << endl;
 print_R << "Frec_marg_F1" << endl;
 print_R << colsum(pobs_arr) << endl;
 print_R << colsum(ppred_arr) << endl;
 print_R << "Frec_marg_F2" << endl;
 print_R << colsum(pobs_pal) << endl;
 print_R << colsum(ppred_pal) << endl;
 print_R << "Frec_marg_F3" << endl;
 print_R << colsum(pobs_esp) << endl;
 print_R << colsum(ppred_esp) << endl;
 print_R << "Frec_marg_cru" << endl;
 print_R << colsum(pobs_cru) << endl;
 print_R << colsum(ppred_cru) << endl;
 print_R << "B_B0" << endl;
 print_R << SPR << endl;
 print_R << "SPR_sd" << endl;
 print_R << SPR.sd << endl;
 print_R << "Ftot" << endl;
 print_R << Ftot << endl;
 print_R << "F_sd" << endl;
 print_R << Ftot.sd << endl;
 print_R << "SSB" << endl;
 print_R << BD << endl;
 print_R << "BD_sd" << endl;
 print_R << BD.sd << endl;
 print_R << "BT" << endl;
 print_R << (N*Prob_talla)*Wmed<< endl;
 print_R << "Reclutas" << endl;
 print_R << Reclutas << endl;
 print_R << Rpred << endl;
 print_R << log_desv_Rt << endl;
 print_R << "R_sd" << endl;
 print_R << Reclutas.sd << endl;
 print_R << "Nanos_frec_F1" << endl;
 print_R<< nanos_farr<< endl;
 print_R << "pobs_F1" << endl;
 print_R << pobs_arr << endl;
 print_R << "ppred_F1" << endl;
 print_R << ppred_arr << endl;

 print_R << "Nanos_frec_F2" << endl;
 print_R<< nanos_fpal<< endl;
 print_R << "pobs_F2" << endl;
 print_R << pobs_pal << endl;
 print_R << "ppred_F2" << endl;
 print_R << ppred_pal << endl;
 
 print_R << "Nanos_frec_F3" << endl;
 print_R<< nanos_fesp<< endl;
 print_R << "pobs_F3" << endl;
 print_R << pobs_esp << endl;
 print_R << "ppred_F3" << endl;
 print_R << ppred_esp << endl;
 
 print_R << "Nanos_frec_cru" << endl;
 print_R<< nanos_fcru<< endl;
 print_R << "pobs_cru" << endl;
 print_R << pobs_cru << endl;
 print_R << "ppred_cru" << endl;
 print_R << ppred_cru << endl;
 print_R << "M" << endl;
 print_R << M << endl;
 print_R << "h" << endl;
 print_R << h << endl;
 print_R << "dts" << endl;
 print_R << dt(1) << endl;
 print_R << "min_edad" << endl;
 print_R << minedad << endl;
 print_R << "R0" << endl;
 print_R << mfexp(log_Rmed) << endl;
 print_R << "B0" << endl;
 print_R << SSBo << endl;
 print_R << "Linf_k_Lo_alfa_beta_M_h" << endl;
 print_R <<Linf<<" "<<k<<" "<<L0<<" "<<mfexp(log_alfa)<<" "<<mfexp(log_beta)<<" "<<M<<" "<<h<< endl;
 print_R << "L50ms" << endl;
 print_R << L50ms << endl;
 print_R << "Wtalla" << endl;
 print_R << Wmed << endl;
 print_R << "MStalla" << endl;
 print_R << msex << endl;
 print_R << "N_tot" << endl;
 print_R << N << endl;
 print_R << "Z_tot" << endl;
 print_R << Z << endl;

 print_R << "Fa_F1" << endl;
 print_R << Farr << endl;

 print_R << "Fa_F2" << endl;
 print_R << Fpal << endl;

 print_R << "Fa_F3" << endl;
 print_R << Fesp << endl;

 print_R << "F_F1" << endl;
 print_R << mfexp(log_Fa) << endl;
 print_R << "F_F2" << endl;
 print_R << mfexp(log_Fp) << endl;
 print_R << "F_F3" << endl;
 print_R << mfexp(log_Fe) << endl;

 print_R << "Mult_F" << endl;
 print_R << pbr << endl;
 print_R << "Bio_proy" << endl;
 print_R << trans(Bp) << endl;
 print_R << "Capt_proy" << endl;
 print_R << trans(Yp) << endl;
 print_R << "F_proy" << endl;
 print_R << trans(Fproy) << endl;
 print_R << "Red_stock" << endl;
 print_R << Redstock << endl;
 print_R << Redstock.sd << endl;

 print_R << "Likeval" << endl;
 print_R <<  likeval << endl;
 print_R << "MaxGrad" << endl;
 print_R <<objective_function_value::pobjfun->gmax<<endl;
 print_R << "FunObj" << endl;
 print_R <<  f << endl;
