#V3.30.18.00;_fast(opt);_compile_date:_Sep 30 2021;_Stock_Synthesis_by_Richard_Methot_(NOAA)_using_ADMB_12.3
#_Stock_Synthesis_is_a_work_of_the_U.S._Government_and_is_not_subject_to_copyright_protection_in_the_United_States.
#_Foreign_copyrights_may_apply._See_copyright.txt_for_more_information.
#_User_support_available_at:NMFS.Stock.Synthesis@noaa.gov
#_User_info_available_at:https://vlab.noaa.gov/group/stock-synthesis
#_Source_code_at:_https://github.com/nmfs-stock-synthesis/stock-synthesis
#C L budegassa (black angler) in areas 7,8abde
#_data_and_control_files: data.ss // control.ss
0  # 0 means do not read wtatage.ss; 1 means read and use wtatage.ss and also read and use growth parameters
1  #_N_Growth_Patterns (Growth Patterns, Morphs, Bio Patterns, GP are terms used interchangeably in SS)
1 #_N_platoons_Within_GrowthPattern 
#_Cond 1 #_Platoon_within/between_stdev_ratio (no read if N_platoons=1)
#_Cond  1 #vector_platoon_dist_(-1_in_first_val_gives_normal_approx)
#
4 # recr_dist_method for parameters:  2=main effects for GP, Area, Settle timing; 3=each Settle entity; 4=none (only when N_GP*Nsettle*pop==1)
1 # not yet implemented; Future usage: Spawner-Recruitment: 1=global; 2=by area
1 #  number of recruitment settlement assignments 
0 # unused option
#GPattern month  area  age (for each settlement assignment)
 1 1 1 0
#
#_Cond 0 # N_movement_definitions goes here if Nareas > 1
#_Cond 1.0 # first age that moves (real age at begin of season, not integer) also cond on do_migration>0
#_Cond 1 1 1 2 4 10 # example move definition for seas=1, morph=1, source=1 dest=2, age1=4, age2=10
#
0 #_Nblock_Patterns
#_Cond 0 #_blocks_per_pattern 
# begin and end years of blocks
#
# controls for all timevary parameters 
1 #_time-vary parm bound check (1=warn relative to base parm bounds; 3=no bound check); Also see env (3) and dev (5) options to constrain with base bounds
#
# AUTOGEN
 1 1 1 1 1 # autogen: 1st element for biology, 2nd for SR, 3rd for Q, 4th reserved, 5th for selex
# where: 0 = autogen time-varying parms of this category; 1 = read each time-varying parm line; 2 = read then autogen if parm min==-12345
#
#_Available timevary codes
#_Block types: 0: P_block=P_base*exp(TVP); 1: P_block=P_base+TVP; 2: P_block=TVP; 3: P_block=P_block(-1) + TVP
#_Block_trends: -1: trend bounded by base parm min-max and parms in transformed units (beware); -2: endtrend and infl_year direct values; -3: end and infl as fraction of base range
#_EnvLinks:  1: P(y)=P_base*exp(TVP*env(y));  2: P(y)=P_base+TVP*env(y);  3: P(y)=f(TVP,env_Zscore) w/ logit to stay in min-max;  4: P(y)=2.0/(1.0+exp(-TVP1*env(y) - TVP2))
#_DevLinks:  1: P(y)*=exp(dev(y)*dev_se;  2: P(y)+=dev(y)*dev_se;  3: random walk;  4: zero-reverting random walk with rho;  5: like 4 with logit transform to stay in base min-max
#_DevLinks(more):  21-25 keep last dev for rest of years
#
#_Prior_codes:  0=none; 6=normal; 1=symmetric beta; 2=CASAL's beta; 3=lognormal; 4=lognormal with biascorr; 5=gamma
#
# setup for M, growth, wt-len, maturity, fecundity, (hermaphro), recr_distr, cohort_grow, (movement), (age error), (catch_mult), sex ratio 
#_NATMORT
3 #_natM_type:_0=1Parm; 1=N_breakpoints;_2=Lorenzen;_3=agespecific;_4=agespec_withseasinterpolate;_5=BETA:_Maunder_link_to_maturity
 #_Age_natmort_by sex x growthpattern (nest GP in sex)
 1 0.72 0.44 0.34 0.32 0.32 0.32 0.32 0.32 0.32 0.32 0.32 0.32 0.32 0.32 0.32
 1 0.72 0.44 0.34 0.32 0.32 0.32 0.32 0.32 0.32 0.32 0.32 0.32 0.32 0.32 0.32
#
1 # GrowthModel: 1=vonBert with L1&L2; 2=Richards with L1&L2; 3=age_specific_K_incr; 4=age_specific_K_decr; 5=age_specific_K_each; 6=NA; 7=NA; 8=growth cessation
1 #_Age(post-settlement)_for_L1;linear growth below this
999 #_Growth_Age_for_L2 (999 to use as Linf)
-999 #_exponential decay for growth above maxage (value should approx initial Z; -999 replicates 3.24; -998 to not allow growth above maxage)
0  #_placeholder for future growth feature
#
0 #_SD_add_to_LAA (set to 0.1 for SS2 V1.x compatibility)
0 #_CV_Growth_Pattern:  0 CV=f(LAA); 1 CV=F(A); 2 SD=F(LAA); 3 SD=F(A); 4 logSD=F(A)
#
1 #_maturity_option:  1=length logistic; 2=age logistic; 3=read age-maturity matrix by growth_pattern; 4=read age-fecundity; 5=disabled; 6=read length-maturity
2 #_First_Mature_Age
1 #_fecundity option:(1)eggs=Wt*(a+b*Wt);(2)eggs=a*L^b;(3)eggs=a*Wt^b; (4)eggs=a+b*L; (5)eggs=a+b*W
0 #_hermaphroditism option:  0=none; 1=female-to-male age-specific fxn; -1=male-to-female age-specific fxn
1 #_parameter_offset_approach for M, G, CV_G:  1- direct, no offset**; 2- male=fem_parm*exp(male_parm); 3: male=female*exp(parm) then old=young*exp(parm)
#_** in option 1, any male parameter with value = 0.0 and phase <0 is set equal to female parameter
#
#_growth_parms
#_ LO HI INIT PRIOR PR_SD PR_type PHASE env_var&link dev_link dev_minyr dev_maxyr dev_PH Block Block_Fxn
# Sex: 1  BioPattern: 1  NatMort
# Sex: 1  BioPattern: 1  Growth
 3 30 12.5 12.5 0 0 -2 0 0 0 0 0 0 0 # L_at_Amin_Fem_GP_1
 50 150 129 100 1 0 -2 0 0 0 0 0 0 0 # L_at_Amax_Fem_GP_1
 0.01 0.5 0.101 0.1 0 0 -2 0 0 0 0 0 0 0 # VonBert_K_Fem_GP_1
 0.05 0.5 0.244 0.1 0 0 -2 0 0 0 0 0 0 0 # CV_young_Fem_GP_1
 0.05 0.5 0.1 0.1 0 0 -2 0 0 0 0 0 0 0 # CV_old_Fem_GP_1
# Sex: 1  BioPattern: 1  WtLen
 -1 1 1.77e-05 2.7e-05 0.05 0 -3 0 0 0 0 0 0 0 # Wtlen_1_Fem_GP_1
 2 4 2.95 2.839 0.05 0 -3 0 0 0 0 0 0 0 # Wtlen_2_Fem_GP_1
# Sex: 1  BioPattern: 1  Maturity&Fecundity
 40 70 65 65 5 0 -3 0 0 0 0 0 0 0 # Mat50%_Fem_GP_1
 -1 1 -0.15 -0.5 0.5 0 -3 0 0 0 0 0 0 0 # Mat_slope_Fem_GP_1
 -3 3 1 1 0.8 0 -3 0 0 0 0 0 0 0 # Eggs/kg_inter_Fem_GP_1
 -3 3 0 0 0.8 0 -3 0 0 0 0 0 0 0 # Eggs/kg_slope_wt_Fem_GP_1
# Sex: 2  BioPattern: 1  NatMort
# Sex: 2  BioPattern: 1  Growth
 3 30 12.5 12.5 0 0 -2 0 0 0 0 0 0 0 # L_at_Amin_Mal_GP_1
 50 150 78 60 1 0 -2 0 0 0 0 0 0 0 # L_at_Amax_Mal_GP_1
 0.01 0.5 0.197 0.305 0 0 -2 0 0 0 0 0 0 0 # VonBert_K_Mal_GP_1
 0.05 0.5 0.244 0.1 0 0 -2 0 0 0 0 0 0 0 # CV_young_Mal_GP_1
 0.05 0.5 0.1 0.1 0 0 -2 0 0 0 0 0 0 0 # CV_old_Mal_GP_1
# Sex: 2  BioPattern: 1  WtLen
 -1 1 1.77e-05 2.7e-05 0.05 0 -3 0 0 0 0 0 0 0 # Wtlen_1_Mal_GP_1
 2 4 2.95 2.839 0.05 0 -3 0 0 0 0 0 0 0 # Wtlen_2_Mal_GP_1
# Hermaphroditism
#  Recruitment Distribution  
#  Cohort growth dev base
 0.1 10 1 1 1 0 -1 0 0 0 0 0 0 0 # CohortGrowDev
#  Movement
#  Age Error from parameters
#  catch multiplier
#  fraction female, by GP
 1e-06 0.999999 0.5 0.5 0.5 0 -99 0 0 0 0 0 0 0 # FracFemale_GP_1
#  M2 parameter for each predator fleet
#
#_no timevary MG parameters
#
#_seasonal_effects_on_biology_parms
 0 0 0 0 0 0 0 0 0 0 #_femwtlen1,femwtlen2,mat1,mat2,fec1,fec2,Malewtlen1,malewtlen2,L1,K
#_ LO HI INIT PRIOR PR_SD PR_type PHASE
#_Cond -2 2 0 0 -1 99 -2 #_placeholder when no seasonal MG parameters
#
3 #_Spawner-Recruitment; Options: 1=NA; 2=Ricker; 3=std_B-H; 4=SCAA; 5=Hockey; 6=B-H_flattop; 7=survival_3Parm; 8=Shepherd_3Parm; 9=RickerPower_3parm
0  # 0/1 to use steepness in initial equ recruitment calculation
0  #  future feature:  0/1 to make realized sigmaR a function of SR curvature
#_          LO            HI          INIT         PRIOR         PR_SD       PR_type      PHASE    env-var    use_dev   dev_mnyr   dev_mxyr     dev_PH      Block    Blk_Fxn #  parm_name
           1.5            30       11.9215            12             0             0          1          0          0          0          0          0          0          0 # SR_LN(R0)
           0.2         0.999          0.93          0.93           0.2             0         -5          0          0          0          0          0          0          0 # SR_BH_steep
           0.1             2           0.5          0.62           0.2             0         -3          0          0          0          0          0          0          0 # SR_sigmaR
            -5             5             0             0             2             0         -2          0          0          0          0          0          0          0 # SR_regime
             0             0             0             0            -1             0        -99          0          0          0          0          0          0          0 # SR_autocorr
#_no timevary SR parameters
1 #do_recdev:  0=none; 1=devvector (R=F(SSB)+dev); 2=deviations (R=F(SSB)+dev); 3=deviations (R=R0*dev; dev2=R-f(SSB)); 4=like 3 with sum(dev2) adding penalty
1986 # first year of main recr_devs; early devs can preceed this era
2025 # last year of main recr_devs; forecast devs start in following year
2 #_recdev phase 
1 # (0/1) to read 13 advanced options
 0 #_recdev_early_start (0=none; neg value makes relative to recdev_start)
 -3 #_recdev_early_phase
 0 #_forecast_recruitment phase (incl. late recr) (0 value resets to maxphase+1)
 1 #_lambda for Fcast_recr_like occurring before endyr+1
1946.0   #_last_early_yr_nobias_adj_in_MPD 
1992.4   #_first_yr_fullbias_adj_in_MPD 
2025.9   #_last_yr_fullbias_adj_in_MPD 
2026.0   #_first_recent_yr_nobias_adj_in_MPD 
0.9397  #_max_bias_adj_in_MPD (1.0 to mimic pre-2009 models) 
  0 #_period of cycles in recruitment (N parms read below)
 -5 #min rec_dev
 5 #max rec_dev
 0 #_read_recdevs
#_end of advanced SR options
#
#_placeholder for full parameter lines for recruitment cycles
# read specified recr devs
#_Yr Input_value
#
# all recruitment deviations
#  1986R 1987R 1988R 1989R 1990R 1991R 1992R 1993R 1994R 1995R 1996R 1997R 1998R 1999R 2000R 2001R 2002R 2003R 2004R 2005R 2006R 2007R 2008R 2009R 2010R 2011R 2012R 2013R 2014R 2015R 2016R 2017R 2018R 2019R 2020R 2021R 2022R 2023F
#  0.0367464 -0.456123 -0.288607 -0.139215 -0.309589 -0.104274 -0.00899076 0.0836757 -0.527035 -0.245062 -0.340599 -0.958691 -0.849222 -0.455147 0.965419 -1.08284 -1.20495 -0.771847 0.790415 0.178464 0.0610966 0.524358 0.464641 -0.862543 -0.0861641 0.433608 0.402915 0.717749 0.748063 0.110517 0.323259 0.518906 0.247722 0.733804 1.45807 -0.184764 0.0762489 0
#
#Fishing Mortality info 
0.3 # F ballpark value in units of annual_F
-2001 # F ballpark year (neg value to disable)
4 # F_Method:  1=Pope midseason rate; 2=F as parameter; 3=F as hybrid; 4=fleet-specific parm/hybrid (#4 is superset of #2 and #3 and is recommended)
2.9 # max F (methods 2-4) or harvest fraction (method 1)
# read list of fleets that do F as parameter; unlisted fleets stay hybrid, bycatch fleets must be included with start_PH=1, high F fleets should switch early
# (A) fleet, (B) F_initial_value (used if start_PH=1), (C) start_PH for parms (99 to stay in hybrid)
# (A) (B) (C)  (terminate list with -9999 for fleet)
 1 0.3 1 # FL1_Trawls
 2 0.05 1 # FL2_Gillnets
-9999 1 1 # end of list
7 #_number of loops for hybrid tuning; 4 good; 3 faster; 2 enough if switching to parms is enabled
#
#_initial_F_parms; for each fleet x season that has init_catch; nest season in fleet; count = 1
#_for unconstrained init_F, use an arbitrary initial catch and set lambda=0 for its logL
#_ LO HI INIT PRIOR PR_SD  PR_type  PHASE
 1e-06 1 0.01 0.05 0.1 -6 -1 # InitF_seas_1_flt_1FL1_Trawls
#
# F rates by fleet x season
# Yr:  1950 1951 1952 1953 1954 1955 1956 1957 1958 1959 1960 1961 1962 1963 1964 1965 1966 1967 1968 1969 1970 1971 1972 1973 1974 1975 1976 1977 1978 1979 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 2004 2005 2006 2007 2008 2009 2010 2011 2012 2013 2014 2015 2016 2017 2018 2019 2020 2021 2022 2023
# seas:  1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
# FL1_Trawls 0.00988474 0.0101955 0.0108922 0.0117135 0.012548 0.0142634 0.0160165 0.0178145 0.0196638 0.0215706 0.0236489 0.0277258 0.0261368 0.0297285 0.0336764 0.0410418 0.0477212 0.0617752 0.0643162 0.0739124 0.0797065 0.0910897 0.0793111 0.0592874 0.0695613 0.0754863 0.0798845 0.085526 0.107104 0.109196 0.182202 0.211777 0.229814 0.271322 0.259762 0.275866 0.196377 0.173461 0.186899 0.202343 0.233881 0.22836 0.224783 0.181532 0.150495 0.170698 0.195769 0.202838 0.23427 0.213943 0.237109 0.209121 0.182878 0.168552 0.213412 0.153946 0.142994 0.10414 0.1459 0.155855 0.14992 0.151468 0.141741 0.176371 0.15977 0.112339 0.139817 0.13639 0.086883 0.100412 0.0843388 0.0761888 0.0906329 0.0869797
# FL2_Gillnets 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.000288845 0.0033162 0.00389625 0.00556266 0.000617341 0.00222155 0.00178664 0.0123944 0.014573 0.0174239 0.0173428 0.0117376 0.00804523 0.0106263 0.00950482 0.00822053 0.00554769 0.00881745 0.00885067 0.0128488 0.0133938 0.00733455 0.0119502 0.00806354 0.00938194 0.00839985 0.00435435 0.00458394 0.00737872
#
#_Q_setup for fleets with cpue or survey data
#_1:  fleet number
#_2:  link type: (1=simple q, 1 parm; 2=mirror simple q, 1 mirrored parm; 3=q and power, 2 parm; 4=mirror with offset, 2 parm)
#_3:  extra input for link, i.e. mirror fleet# or dev index number
#_4:  0/1 to select extra sd parameter
#_5:  0/1 for biasadj or not
#_6:  0/1 to float
#_   fleet      link link_info  extra_se   biasadj     float  #  fleetname
         3         1         0         1         0         1  #  FR_IE_IBTS
         4         1         0         1         0         1  #  IE_IAMS
-9999 0 0 0 0 0
#
#_Q_parms(if_any);Qunits_are_ln(q)
#_          LO            HI          INIT         PRIOR         PR_SD       PR_type      PHASE    env-var    use_dev   dev_mnyr   dev_mxyr     dev_PH      Block    Blk_Fxn  #  parm_name
           -15            -5      -10.0727           -13             1             0         -1          0          0          0          0          0          0          0  #  LnQ_base_FR_IE_IBTS(3)
             0             1             0           0.1           0.1             0         -4          0          0          0          0          0          0          0  #  Q_extraSD_FR_IE_IBTS(3)
           -15            -5       -7.3398           -13             1             0         -1          0          0          0          0          0          0          0  #  LnQ_base_IE_IAMS(4)
             0             1             0           0.1           0.1             0         -4          0          0          0          0          0          0          0  #  Q_extraSD_IE_IAMS(4)
#_no timevary Q parameters
#
#_size_selex_patterns
#Pattern:_0;  parm=0; selex=1.0 for all sizes
#Pattern:_1;  parm=2; logistic; with 95% width specification
#Pattern:_2;  parm=6; modification of pattern 24 with improved sex-specific offset
#Pattern:_5;  parm=2; mirror another size selex; PARMS pick the min-max bin to mirror
#Pattern:_11; parm=2; selex=1.0  for specified min-max population length bin range
#Pattern:_15; parm=0; mirror another age or length selex
#Pattern:_6;  parm=2+special; non-parm len selex
#Pattern:_43; parm=2+special+2;  like 6, with 2 additional param for scaling (average over bin range)
#Pattern:_8;  parm=8; double_logistic with smooth transitions and constant above Linf option
#Pattern:_9;  parm=6; simple 4-parm double logistic with starting length; parm 5 is first length; parm 6=1 does desc as offset
#Pattern:_21; parm=2+special; non-parm len selex, read as pairs of size, then selex
#Pattern:_22; parm=4; double_normal as in CASAL
#Pattern:_23; parm=6; double_normal where final value is directly equal to sp(6) so can be >1.0
#Pattern:_24; parm=6; double_normal with sel(minL) and sel(maxL), using joiners
#Pattern:_25; parm=3; exponential-logistic in length
#Pattern:_27; parm=special+3; cubic spline in length; parm1==1 resets knots; parm1==2 resets all 
#Pattern:_42; parm=special+3+2; cubic spline; like 27, with 2 additional param for scaling (average over bin range)
#_discard_options:_0=none;_1=define_retention;_2=retention&mortality;_3=all_discarded_dead;_4=define_dome-shaped_retention
#_Pattern Discard Male Special
 1 1 0 0 # 1 FL1_Trawls
 1 1 0 0 # 2 FL2_Gillnets
 24 0 0 0 # 3 FR_IE_IBTS
 1 0 0 0 # 4 IE_IAMS
#
#_age_selex_patterns
#Pattern:_0; parm=0; selex=1.0 for ages 0 to maxage
#Pattern:_10; parm=0; selex=1.0 for ages 1 to maxage
#Pattern:_11; parm=2; selex=1.0  for specified min-max age
#Pattern:_12; parm=2; age logistic
#Pattern:_13; parm=8; age double logistic
#Pattern:_14; parm=nages+1; age empirical
#Pattern:_15; parm=0; mirror another age or length selex
#Pattern:_16; parm=2; Coleraine - Gaussian
#Pattern:_17; parm=nages+1; empirical as random walk  N parameters to read can be overridden by setting special to non-zero
#Pattern:_41; parm=2+nages+1; // like 17, with 2 additional param for scaling (average over bin range)
#Pattern:_18; parm=8; double logistic - smooth transition
#Pattern:_19; parm=6; simple 4-parm double logistic with starting age
#Pattern:_20; parm=6; double_normal,using joiners
#Pattern:_26; parm=3; exponential-logistic in age
#Pattern:_27; parm=3+special; cubic spline in age; parm1==1 resets knots; parm1==2 resets all 
#Pattern:_42; parm=2+special+3; // cubic spline; with 2 additional param for scaling (average over bin range)
#Age patterns entered with value >100 create Min_selage from first digit and pattern from remainder
#_Pattern Discard Male Special
 0 0 0 0 # 1 FL1_Trawls
 0 0 0 0 # 2 FL2_Gillnets
 0 0 0 0 # 3 FR_IE_IBTS
 0 0 0 0 # 4 IE_IAMS
#
#_          LO            HI          INIT         PRIOR         PR_SD       PR_type      PHASE    env-var    use_dev   dev_mnyr   dev_mxyr     dev_PH      Block    Blk_Fxn  #  parm_name
# 1   FL1_Trawls LenSelex
             5            30        11.894            15            99             0          5          0          0          0          0          0          0          0  #  Size_inflection_FL1_Trawls(1)
          0.01             5       1.65529            15            99             0          5          0          0          0          0          0          0          0  #  Size_95%width_FL1_Trawls(1)
            15            40       25.2715            25            99            -6          5          0         23       2005       2025          7          0          0  #  Retain_L_infl_FL1_Trawls(1)
             1             4       2.10552             2            99            -6          6          0          0          0          0          0          0          0  #  Retain_L_width_FL1_Trawls(1)
          -999           999           999            10            99            -6         -3          0          0          0          0          0          0          0  #  Retain_L_asymptote_logit_FL1_Trawls(1)
             0             0             0             0            -1          0.01         -3          0          0          0          0          0          0          0  #  Retain_L_maleoffset_FL1_Trawls(1)
# 2   FL2_Gillnets LenSelex
             5            30       11.7471            15            99             0          5          0          0          0          0          0          0          0  #  Size_inflection_FL2_Gillnets(2)
          0.01             5       1.76301            15            99             0          5          0          0          0          0          0          0          0  #  Size_95%width_FL2_Gillnets(2)
            15            40       34.8807            25            99            -6          5          0          0          0          0          0          0          0  #  Retain_L_infl_FL2_Gillnets(2)
             1             4       2.61899             2            99            -6          6          0          0          0          0          0          0          0  #  Retain_L_width_FL2_Gillnets(2)
          -999           999           999            10            99            -6         -3          0          0          0          0          0          0          0  #  Retain_L_asymptote_logit_FL2_Gillnets(2)
             0             0             0             0            -1          0.01         -3          0          0          0          0          0          0          0  #  Retain_L_maleoffset_FL2_Gillnets(2)
# 3   FR_IE_IBTS LenSelex
            10            25       16.7654            15            99             0          5          0          0          0          0          0          0          0  #  Size_DblN_peak_FR_IE_IBTS(3)
            -5             2     -0.514901           1.5            99             0         -6          0          0          0          0          0          0          0  #  Size_DblN_top_logit_FR_IE_IBTS(3)
             1            15       3.96614             3            99             0          6          0          0          0          0          0          0          0  #  Size_DblN_ascend_se_FR_IE_IBTS(3)
           -25            15        6.1057             5            99             0          6          0          0          0          0          0          0          0  #  Size_DblN_descend_se_FR_IE_IBTS(3)
          -999          -999          -999          -999             0             0         -2          0          0          0          0          0          0          0  #  Size_DblN_start_logit_FR_IE_IBTS(3)
          -999          -999          -999          -999             0             0         -2          0          0          0          0          0          0          0  #  Size_DblN_end_logit_FR_IE_IBTS(3)
# 4   IE_IAMS LenSelex
             5            30       13.1515            15            99             0          5          0          0          0          0          0          0          0  #  Size_inflection_IE_IAMS(4)
          0.01            10       2.87246            15            99             0          5          0          0          0          0          0          0          0  #  Size_95%width_IE_IAMS(4)
# 1   FL1_Trawls AgeSelex
# 2   FL2_Gillnets AgeSelex
# 3   FR_IE_IBTS AgeSelex
# 4   IE_IAMS AgeSelex
#_Dirichlet parameters
#_multiple_fleets_can_refer_to_same_parm;_but_list_cannot_have_gaps
            -5             5       1.19083             0         1.813             6          8          0          0          0          0          0          0          0  #  ln(DM_theta)_1
            -5             5       4.51438             0         1.813             6          8          0          0          0          0          0          0          0  #  ln(DM_theta)_2
            -5             5       1.81934             0         1.813             6          8          0          0          0          0          0          0          0  #  ln(DM_theta)_3
            -5             5        4.5972             0         1.813             6          8          0          0          0          0          0          0          0  #  ln(DM_theta)_4
# timevary selex parameters 
#_          LO            HI          INIT         PRIOR         PR_SD       PR_type    PHASE  #  parm_name
        0.0001            12           2.5           0.5           0.5             6      -5  # Retain_L_infl_FL1_Trawls(1)_dev_se
         -0.99          0.99             0             0           0.5             6      -6  # Retain_L_infl_FL1_Trawls(1)_dev_autocorr
# info on dev vectors created for selex parms are reported with other devs after tag parameter section 
#
0   #  use 2D_AR1 selectivity(0/1)
#_no 2D_AR1 selex offset used
#
# Tag loss and Tag reporting parameters go next
0  # TG_custom:  0=no read and autogen if tag data exist; 1=read
#_Cond -6 6 1 1 2 0.01 -4 0 0 0 0 0 0 0  #_placeholder if no parameters
#
# deviation vectors for timevary parameters
#  base   base first block   block  env  env   dev   dev   dev   dev   dev
#  type  index  parm trend pattern link  var  vectr link _mnyr  mxyr phase  dev_vector
#      5     3     1     0     0     0     0     1     3  2005  2022     7 1.35886 -0.885866 1.21071 0.293259 0.183967 0.0780434 0.186033 0.410773 -0.261179 0.525751 0.594156 0.0281647 -0.357143 0.0381957 -0.348989 -0.845477 1.46327 0.239291
     #
# Input variance adjustments factors: 
 #_1=add_to_survey_CV
 #_2=add_to_discard_stddev
 #_3=add_to_bodywt_CV
 #_4=mult_by_lencomp_N
 #_5=mult_by_agecomp_N
 #_6=mult_by_size-at-age_N
 #_7=mult_by_generalized_sizecomp
#_Factor  Fleet  Value
      4      1         1
      4      2         1
      4      3         1
      4      4         1
 -9999   1    0  # terminator
#
1 #_maxlambdaphase
1 #_sd_offset; must be 1 if any growthCV, sigmaR, or survey extraSD is an estimated parameter
# read 0 changes to default Lambdas (default value is 1.0)
# Like_comp codes:  1=surv; 2=disc; 3=mnwt; 4=length; 5=age; 6=SizeFreq; 7=sizeage; 8=catch; 9=init_equ_catch; 
# 10=recrdev; 11=parm_prior; 12=parm_dev; 13=CrashPen; 14=Morphcomp; 15=Tag-comp; 16=Tag-negbin; 17=F_ballpark; 18=initEQregime
#like_comp fleet  phase  value  sizefreq_method
-9999  1  1  1  1  #  terminator
#
# lambdas (for info only; columns are phases)
#  0 #_CPUE/survey:_1
#  0 #_CPUE/survey:_2
#  1 #_CPUE/survey:_3
#  1 #_CPUE/survey:_4
#  1 #_discard:_1
#  1 #_discard:_2
#  0 #_discard:_3
#  0 #_discard:_4
#  1 #_lencomp:_1
#  1 #_lencomp:_2
#  1 #_lencomp:_3
#  1 #_lencomp:_4
#  1 #_init_equ_catch1
#  1 #_init_equ_catch2
#  1 #_init_equ_catch3
#  1 #_init_equ_catch4
#  1 #_recruitments
#  1 #_parameter-priors
#  1 #_parameter-dev-vectors
#  1 #_crashPenLambda
#  0 # F_ballpark_lambda
0 # (0/1/2) read specs for more stddev reporting: 0 = skip, 1 = read specs for reporting stdev for selectivity, size, and numbers, 2 = add options for M,Dyn. Bzero, SmryBio
 # 0 2 0 0 # Selectivity: (1) fleet, (2) 1=len/2=age/3=both, (3) year, (4) N selex bins
 # 0 0 # Growth: (1) growth pattern, (2) growth ages
 # 0 0 0 # Numbers-at-age: (1) area(-1 for all), (2) year, (3) N ages
 # -1 # list of bin #'s for selex std (-1 in first bin to self-generate)
 # -1 # list of ages for growth std (-1 in first bin to self-generate)
 # -1 # list of ages for NatAge std (-1 in first bin to self-generate)
999

