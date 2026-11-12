#V3.30.18.00;_safe;_compile_date:_Sep 30 2021;_Stock_Synthesis_by_Richard_Methot_(NOAA)_using_ADMB_12.3
#_Stock_Synthesis_is_a_work_of_the_U.S._Government_and_is_not_subject_to_copyright_protection_in_the_United_States.
#_Foreign_copyrights_may_apply._See_copyright.txt_for_more_information.
#_User_support_available_at:NMFS.Stock.Synthesis@noaa.gov
#_User_info_available_at:https://vlab.noaa.gov/group/stock-synthesis
#_Source_code_at:_https://github.com/nmfs-stock-synthesis/stock-synthesis

#_data_and_control_files: Pandalus.dat // Pandalus.ctl
0  # 0 means do not read wtatage.ss; 1 means read and use wtatage.ss and also read and use growth parameters
1  #_N_Growth_Patterns (Growth Patterns, Morphs, Bio Patterns, GP are terms used interchangeably in SS)
1 #_N_platoons_Within_GrowthPattern 
#_Cond 1 #_Platoon_within/between_stdev_ratio (no read if N_platoons=1)
#_Cond  1 #vector_platoon_dist_(-1_in_first_val_gives_normal_approx)
#
3 # recr_dist_method for parameters:  2=main effects for GP, Area, Settle timing; 3=each Settle entity; 4=none (only when N_GP*Nsettle*pop==1)
1 # not yet implemented; Future usage: Spawner-Recruitment: 1=global; 2=by area
2 #  number of recruitment settlement assignments 
0 # unused option
#GPattern month  area  age (for each settlement assignment)
 1 1 1 0
 1 1 2 0
#
0 #_N_movement_definitions
#_Cond 1.0 # first age that moves (real age at begin of season, not integer) if do_migration>0
#_Cond 1 1 1 2 4 10 # example move definition for seas=1, GP=1, source=1 dest=2, age1=4, age2=10
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
1 #_natM_type:_0=1Parm; 1=N_breakpoints;_2=Lorenzen;_3=agespecific;_4=agespec_withseasinterpolate;_5=BETA:_Maunder_link_to_maturity
4 #_N_breakpoints
 0.5 1.5 3.5 7.5 # age(real) at M breakpoints
#
1 # GrowthModel: 1=vonBert with L1&L2; 2=Richards with L1&L2; 3=age_specific_K_incr; 4=age_specific_K_decr; 5=age_specific_K_each; 6=NA; 7=NA; 8=growth cessation
0.5 #_Age(post-settlement)_for_L1;linear growth below this
999 #_Growth_Age_for_L2 (999 to use as Linf)
-999 #_exponential decay for growth above maxage (value should approx initial Z; -999 replicates 3.24; -998 to not allow growth above maxage)
0  #_placeholder for future growth feature
#
0 #_SD_add_to_LAA (set to 0.1 for SS2 V1.x compatibility)
1 #_CV_Growth_Pattern:  0 CV=f(LAA); 1 CV=F(A); 2 SD=F(LAA); 3 SD=F(A); 4 logSD=F(A)
#
1 #_maturity_option:  1=length logistic; 2=age logistic; 3=read age-maturity matrix by growth_pattern; 4=read age-fecundity; 5=disabled; 6=read length-maturity
2 #_First_Mature_Age
1 #_fecundity option:(1)eggs=Wt*(a+b*Wt);(2)eggs=a*L^b;(3)eggs=a*Wt^b; (4)eggs=a+b*L; (5)eggs=a+b*W
-1 #_hermaphroditism option:  0=none; 1=female-to-male age-specific fxn; -1=male-to-female age-specific fxn
-1 # Hermaphro_season.first_age (seas=-1 means all seasons; first_age must be 0 to 9)
0 # fraction_of_maleSSB_added_to_total_SSB 
2 #_parameter_offset_approach for M, G, CV_G:  1- direct, no offset**; 2- male=fem_parm*exp(male_parm); 3: male=female*exp(parm) then old=young*exp(parm)
#_** in option 1, any male parameter with value = 0.0 and phase <0 is set equal to female parameter
#
#_growth_parms
#_ LO HI INIT PRIOR PR_SD PR_type PHASE env_var&link dev_link dev_minyr dev_maxyr dev_PH Block Block_Fxn
# Sex: 1  BioPattern: 1  NatMort
 0.15 3 1.709 0.2 10 0 -3 0 0 0 0 0 0 0 # NatM_break_1_Fem_GP_1
 0.15 1.4 1.019 0.2 10 0 -3 0 0 0 0 0 0 0 # NatM_break_2_Fem_GP_1
 0.15 0.8 0.699 0.2 10 0 -3 0 0 0 0 0 0 0 # NatM_break_3_Fem_GP_1
 0.15 0.8 0.577 0.2 10 0 -3 0 0 0 0 0 0 0 # NatM_break_4_Fem_GP_1
# Sex: 1  BioPattern: 1  Growth
 0 4 0.695867 0.7 5 0 2 0 0 0 0 0 0 0 # L_at_Amin_Fem_GP_1
 2 4 2.66268 2.9 10 0 2 0 0 0 0 0 0 0 # L_at_Amax_Fem_GP_1
 0.2 0.8 0.470716 0.39 0.05 0 2 0 0 0 0 0 0 0 # VonBert_K_Fem_GP_1
 0.005 0.4 0.118852 0.2 0.1 0 2 0 0 0 0 0 0 0 # CV_young_Fem_GP_1
 0.005 0.4 0.0489836 0.2 0.1 0 2 0 0 0 0 0 0 0 # CV_old_Fem_GP_1
# Sex: 1  BioPattern: 1  WtLen
 -1 1 0.0016 0.0009 0.05 0 -2 0 0 0 0 0 0 0 # Wtlen_1_Fem_GP_1
 -3 4 2.7532 2.7532 0.05 0 -2 0 0 0 0 0 0 0 # Wtlen_2_Fem_GP_1
# Sex: 1  BioPattern: 1  Maturity&Fecundity
 0 3 1.974 2 5.1 0 -2 0 0 0 0 0 0 0 # Mat50%_Fem_GP_1
 -35 5 -10.253 -21.253 0.5 0 -2 0 0 0 0 0 0 0 # Mat_slope_Fem_GP_1
 -3 3 1 1 0.8 0 -2 0 0 0 0 0 0 0 # Eggs/kg_inter_Fem_GP_1
 -3 3 0 0 0.8 0 -2 0 0 0 0 0 0 0 # Eggs/kg_slope_wt_Fem_GP_1
# Sex: 2  BioPattern: 1  NatMort
 -1 1 0 0 0 0 -3 0 0 0 0 0 0 0 # NatM_break_1_Mal_GP_1
 -1 1 0 0 0 0 -3 0 0 0 0 0 0 0 # NatM_break_2_Mal_GP_1
 -1 1 0 0 0 0 -3 0 0 0 0 0 0 0 # NatM_break_3_Mal_GP_1
 -1 1 0 0 0 0 -3 0 0 0 0 0 0 0 # NatM_break_4_Mal_GP_1
# Sex: 2  BioPattern: 1  Growth
 -1 1 0 0 0 0 -3 0 0 0 0 0 0 0 # L_at_Amin_Mal_GP_1
 -1 1 0 0 0 0 -4 0 0 0 0 0 0 0 # L_at_Amax_Mal_GP_1
 -1 1 0 0 0 0 -4 0 0 0 0 0 0 0 # VonBert_K_Mal_GP_1
 -1 1 0 0 0 0 -3 0 0 0 0 0 0 0 # CV_young_Mal_GP_1
 -1 1 0 0 0 0 -3 0 0 0 0 0 0 0 # CV_old_Mal_GP_1
# Sex: 2  BioPattern: 1  WtLen
 -1 1 0.0016 0.0009 0.05 0 -2 0 0 0 0 0 0 0 # Wtlen_1_Mal_GP_1
 -3 4 2.7532 2.7532 0.05 0 -2 0 0 0 0 0 0 0 # Wtlen_2_Mal_GP_1
# Hermaphroditism
 0 8 2.1 1.2 0 0 -2 0 0 0 0 0 0 0 # Herm_Infl_age
 0.1 21.253 0.75 10.1407 0 0 -2 0 0 0 0 0 0 0 # Herm_stdev
 0 1 1 0.890681 0 0 -2 0 0 0 0 0 0 0 # Herm_asymptote
#  Recruitment Distribution  
 -35 25 0 0 0 0 -2 0 0 0 0 0 0 0 # RecrDist_GP_1_area_1_month_1
 -35 25 -0.252805 0 1 0 2 0 23 1984 2024 3 0 0 # RecrDist_GP_1_area_2_month_1
#  Cohort growth dev base
 0.1 10 1 1 1 0 -1 0 0 0 0 0 0 0 # CohortGrowDev
#  Movement
#  Age Error from parameters
#  catch multiplier
#  fraction female, by GP
 1e-06 1 1e-06 1 0 0 -99 0 0 0 0 0 0 0 # FracFemale_GP_1
#  M2 parameter for each predator fleet
#
# timevary MG parameters 
#_ LO HI INIT PRIOR PR_SD PR_type  PHASE
 0.5 2 1.5 1 0.5 6 -7 # RecrDist_GP_1_area_2_month_1_dev_se
 -0.99 0.99 0 0 0.5 6 -6 # RecrDist_GP_1_area_2_month_1_dev_autocorr
# info on dev vectors created for MGparms are reported with other devs after tag parameter section 
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
             3            30       18.1772             6             5             0          1          0          0          0          0          0          0          0 # SR_LN(R0)
           0.1             1          0.72         0.999          0.34             0         -2          0          0          0          0          0          0          0 # SR_BH_steep
             0             2           0.5           1.2           0.8             0         -2          0          0          0          0          0          0          0 # SR_sigmaR
            -5             5             0             0             1             0        -99          0          0          0          0          0          0          0 # SR_regime
             0             1          0.18         0.254         0.054             0         -2          0          0          0          0          0          0          0 # SR_autocorr
#_no timevary SR parameters
1 #do_recdev:  0=none; 1=devvector (R=F(SSB)+dev); 2=deviations (R=F(SSB)+dev); 3=deviations (R=R0*dev; dev2=R-f(SSB)); 4=like 3 with sum(dev2) adding penalty
1984 # first year of main recr_devs; early devs can preceed this era
2024 # last year of main recr_devs; forecast devs start in following year
2 #_recdev phase 
1 # (0/1) to read 13 advanced options
 -6 #_recdev_early_start (0=none; neg value makes relative to recdev_start)
 3 #_recdev_early_phase
 0 #_forecast_recruitment phase (incl. late recr) (0 value resets to maxphase+1)
 1 #_lambda for Fcast_recr_like occurring before endyr+1
 1977 #_last_yr_nobias_adj_in_MPD; begin of ramp
 1983 #_first_yr_fullbias_adj_in_MPD; begin of plateau
 2023 #_last_yr_fullbias_adj_in_MPD
 2024 #_end_yr_for_ramp_in_MPD (can be in forecast to shape ramp, but SS sets bias_adj to 0.0 for fcast yrs)
 0.938 #_max_bias_adj_in_MPD (typical ~0.8; -3 sets all years to 0.0; -2 sets all non-forecast yrs w/ estimated recdevs to 1.0; -1 sets biasadj=1.0 for all yrs w/ recdevs)
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
#  1978E 1979E 1980E 1981E 1982E 1983E 1984R 1985R 1986R 1987R 1988R 1989R 1990R 1991R 1992R 1993R 1994R 1995R 1996R 1997R 1998R 1999R 2000R 2001R 2002R 2003R 2004R 2005R 2006R 2007R 2008R 2009R 2010R 2011R 2012R 2013R 2014R 2015R 2016R 2017R 2018R 2019R 2020R 2021F 2022F 2023F
#  -0.0210723 -0.0512967 -0.176138 -0.382713 -0.336842 -0.0106347 0.241159 -0.360898 -0.0751685 -0.79452 0.370824 0.460548 0.6935 0.148852 0.654718 0.238996 0.334776 0.485059 0.347759 -0.391614 -0.0839029 0.013161 0.435525 0.605406 0.282763 0.367376 -0.0455556 0.139019 -0.333998 -0.251834 -0.955803 -0.985384 -0.669633 0.0312256 0.236621 0.906896 0.271369 -0.256522 -0.0845618 -0.386233 0.306604 -0.347223 -1.5493 -0.386144 -0.0962413 -0.0239869
#
#Fishing Mortality info 
0.537 # F ballpark value in units of annual_F
-2008 # F ballpark year (neg value to disable)
4 # F_Method:  1=Pope midseason rate; 2=F as parameter; 3=F as hybrid; 4=fleet-specific parm/hybrid (#4 is superset of #2 and #3 and is recommended)
10 # max F (methods 2-4) or harvest fraction (method 1)
# read list of fleets that do F as parameter; unlisted fleets stay hybrid, bycatch fleets must be included with start_PH=1, high F fleets should switch early
# (A) fleet, (B) F_initial_value (used if start_PH=1), (C) start_PH for parms (99 to stay in hybrid)
# (A) (B) (C)  (terminate list with -9999 for fleet)
 1 0.05 1 # DK3a
 2 0.05 99 # DK4a
 3 0.05 1 # NOR3a
 4 0.05 1 # NOR4a
 5 0.05 1 # SWE3a
 6 0.05 99 # SWE4a
-9999 1 1 # end of list
4 #_number of loops for hybrid tuning; 4 good; 3 faster; 2 enough if switching to parms is enabled
#
#_initial_F_parms; for each fleet x season that has init_catch; nest season in fleet; count = 16
#_for unconstrained init_F, use an arbitrary initial catch and set lambda=0 for its logL
#_ LO HI INIT PRIOR PR_SD  PR_type  PHASE
 0 8 0.025124 1 99 0 1 # InitF_seas_1_flt_1DK3a
 0 8 0.0398011 1 99 0 1 # InitF_seas_1_flt_3NOR3a
 0 8 0.0865064 1 99 0 1 # InitF_seas_1_flt_4NOR4a
 0 8 0.0616421 1 99 0 1 # InitF_seas_1_flt_5SWE3a
 0 8 0.0364625 1 99 0 1 # InitF_seas_2_flt_1DK3a
 0 8 0.0385716 1 99 0 1 # InitF_seas_2_flt_3NOR3a
 0 8 0.0871043 1 99 0 1 # InitF_seas_2_flt_4NOR4a
 0 8 0.0784148 1 99 0 1 # InitF_seas_2_flt_5SWE3a
 0 8 0.036935 1 99 0 1 # InitF_seas_3_flt_1DK3a
 0 8 0.0295664 1 99 0 1 # InitF_seas_3_flt_3NOR3a
 0 8 0.0407975 1 99 0 1 # InitF_seas_3_flt_4NOR4a
 0 8 0.0671539 1 99 0 1 # InitF_seas_3_flt_5SWE3a
 0 8 0.0232964 1 99 0 1 # InitF_seas_4_flt_1DK3a
 0 8 0.0251639 1 99 0 1 # InitF_seas_4_flt_3NOR3a
 0 8 0.0201022 1 99 0 1 # InitF_seas_4_flt_4NOR4a
 0 8 0.0456068 1 99 0 1 # InitF_seas_4_flt_5SWE3a
#
# F rates by fleet x season
# Yr:  1970 1970 1970 1970 1971 1971 1971 1971 1972 1972 1972 1972 1973 1973 1973 1973 1974 1974 1974 1974 1975 1975 1975 1975 1976 1976 1976 1976 1977 1977 1977 1977 1978 1978 1978 1978 1979 1979 1979 1979 1980 1980 1980 1980 1981 1981 1981 1981 1982 1982 1982 1982 1983 1983 1983 1983 1984 1984 1984 1984 1985 1985 1985 1985 1986 1986 1986 1986 1987 1987 1987 1987 1988 1988 1988 1988 1989 1989 1989 1989 1990 1990 1990 1990 1991 1991 1991 1991 1992 1992 1992 1992 1993 1993 1993 1993 1994 1994 1994 1994 1995 1995 1995 1995 1996 1996 1996 1996 1997 1997 1997 1997 1998 1998 1998 1998 1999 1999 1999 1999 2000 2000 2000 2000 2001 2001 2001 2001 2002 2002 2002 2002 2003 2003 2003 2003 2004 2004 2004 2004 2005 2005 2005 2005 2006 2006 2006 2006 2007 2007 2007 2007 2008 2008 2008 2008 2009 2009 2009 2009 2010 2010 2010 2010 2011 2011 2011 2011 2012 2012 2012 2012 2013 2013 2013 2013 2014 2014 2014 2014 2015 2015 2015 2015 2016 2016 2016 2016 2017 2017 2017 2017 2018 2018 2018 2018 2019 2019 2019 2019 2020 2020 2020 2020 2021 2021 2021 2021 2022 2022 2022 2022 2023 2023 2023 2023
# seas:  1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4
# DK3a 0.0105145 0.0153513 0.0155457 0.0098005 0.0114961 0.0168272 0.0171315 0.0108642 0.0101295 0.0147313 0.0149851 0.00946462 0.00763035 0.0110822 0.0112195 0.00713611 0.00536749 0.00789411 0.00794386 0.00499787 0.008367 0.012256 0.0123826 0.00785825 0.0125469 0.0185276 0.0188855 0.0119425 0.0121529 0.0176431 0.0177658 0.0111045 0.015586 0.0225951 0.0228221 0.0143498 0.0113001 0.0164328 0.0167102 0.0106077 0.0182617 0.0272849 0.0280997 0.0179276 0.0307117 0.0463943 0.0481651 0.0311365 0.0504635 0.0773354 0.0820875 0.0537848 0.0289678 0.0440888 0.0461922 0.0293414 0.031875 0.0458577 0.0459954 0.028585 0.074699 0.104855 0.10122 0.0616446 0.0709632 0.105639 0.112489 0.0745691 0.0792939 0.125794 0.137861 0.093838 0.15604 0.283579 0.179445 0.136539 0.1295 0.366127 0.350862 0.133696 0.110184 0.162215 0.212791 0.107919 0.147501 0.173521 0.183305 0.143697 0.189345 0.152709 0.174253 0.138792 0.0745751 0.108797 0.171803 0.072561 0.0555132 0.109723 0.06016 0.0608686 0.0709588 0.119442 0.0841893 0.0958057 0.136306 0.15712 0.0917037 0.11587 0.0977448 0.12126 0.107094 0.112157 0.0836692 0.115085 0.120388 0.0817846 0.0203703 0.0655954 0.0905022 0.0598396 0.0740655 0.0415662 0.119701 0.117142 0.0351195 0.0629132 0.0695016 0.0418054 0.0491494 0.0560811 0.0804993 0.0821171 0.0793735 0.078428 0.0856885 0.097952 0.0693291 0.127249 0.12277 0.119758 0.071698 0.108512 0.107542 0.0940441 0.124726 0.11623 0.10745 0.106585 0.0901217 0.0711804 0.132147 0.120303 0.0622341 0.0904951 0.17507 0.148314 0.122015 0.131855 0.124143 0.18268 0.136964 0.117254 0.139471 0.113943 0.11896 0.218211 0.279973 0.203005 0.131912 0.175688 0.224553 0.190255 0.243939 0.300782 0.267556 0.201906 0.159018 0.334867 0.293437 0.28641 0.238553 0.231138 0.278678 0.197354 0.203626 0.0569837 0.236028 0.189559 0.248827 0.256557 0.197964 0.281849 0.202996 0.252445 0.325024 0.195838 0.203263 0.323644 0.358081 0.242297 0.359224 0.218758 0.456509 0.65781 0.121227 0.152881 0.17158 0.144029 0.128744 0.179024 0.220849 0.195972 0.172819 0.21701 0.23737 0.194031
# DK4a 0.0272836 0.0305375 0.0105498 0.00943605 0.0293471 0.0330659 0.0114982 0.0104188 0.0254679 0.0288076 0.00999314 0.0090168 0.019348 0.0217102 0.0075851 0.00683382 0.0135855 0.0152108 0.00529378 0.00478745 0.0210042 0.0236665 0.00827224 0.00738683 0.0313667 0.0353461 0.0122751 0.0111201 0.0295678 0.0332146 0.0115952 0.0104161 0.0387119 0.043584 0.0150472 0.0135951 0.0282146 0.0316557 0.0109941 0.00982568 0.0450376 0.0511822 0.0179826 0.0161909 0.0730359 0.0844008 0.0300688 0.0274953 0.118441 0.140774 0.0513042 0.0475713 0.0675554 0.0792143 0.0283789 0.025721 0.074459 0.0837501 0.0285843 0.0253063 0.181373 0.210631 0.0750249 0.0682343 0.221281 0.269217 0.0991966 0.0928683 0.25978 0.302821 0.101132 0.0867327 0.214316 0.130623 0.0146004 0.00477043 0.0804983 0.0765486 0.0711984 0.0265245 0.0176422 0.0512329 0.0207451 0.0337862 0.0220686 0.0347405 0.00708709 0.00861129 0.0529019 0.0548554 0.00862091 0.00824359 0.0537468 0.0339762 0.00355894 0.00363976 0.0160791 0.0124448 0.00180245 0.000160837 0.0043382 0.00971821 0.00656029 0.00278721 0.0290686 0.0195437 0.00571212 0.0076686 0.0243672 0.032687 0.00270723 0.0102252 0.0495682 0.0154261 0.0170122 0.00099782 0.105231 0.0214597 0.0309686 0.0015913 0.0476675 0.052001 0.0241057 0.0130706 0.140996 0.0355989 0.00994343 0 0.037146 0.036459 0.0255314 0.00524783 0.0695097 0.0338191 0.00603584 0.00224652 0.110834 0.0210317 0.00316045 0.000296964 0.059056 0.0130616 0.00320921 0.000353103 0.0333922 0.00851339 0.00122796 0.00253444 0.0182675 0.00307185 0.000739491 0 0.0106383 0.01426 0.00431916 0.000324357 0.040665 0.0310374 0.0119144 0.00787777 0.0244546 0.0173376 0.00609972 0.0146515 0.0916134 0.0620766 0.0364404 0.00340889 0.102174 0.028466 0.0149153 0.0058805 0.0325436 0.00140548 0.000675122 0 0.00329909 0.00443439 0.00608944 0.0188408 0.0341389 0.00326338 0.0275543 0.0166033 0.0256746 0.0220676 0.00725295 0.00991156 0.0702595 0 0.0122046 0.0155775 0.0729969 0.00309279 0.00507413 0.00807517 0.0338651 0 0.0829837 0.0346767 0.0504912 0.001838 0.0433881 0.0188683 0.0781067 0.0609034 0.0268236 0.0196689 0.0812914 0.0674868 0.0315685 0.0239188 0.0979397 0.0766789 0.0334507 0.0241332
# NOR3a 0.041859 0.0407953 0.0312844 0.0265878 0.0586786 0.0574755 0.0443086 0.0378787 0.0545351 0.053241 0.0408763 0.0347953 0.0619694 0.0603835 0.0462589 0.0392874 0.0737552 0.0719976 0.0551045 0.0468949 0.0830254 0.0810033 0.0618926 0.0526726 0.111014 0.109495 0.0844204 0.0722057 0.110329 0.0760893 0.0502132 0.041975 0.066642 0.0645466 0.0646497 0.0373306 0.0707117 0.098805 0.0644134 0.0810327 0.116368 0.157358 0.110842 0.0815788 0.175217 0.136625 0.120021 0.118081 0.210531 0.176993 0.140852 0.126694 0.236164 0.244575 0.151841 0.111909 0.181347 0.150135 0.308453 0.095062 0.217336 0.267989 0.236778 0.181501 0.194473 0.208242 0.261431 0.231653 0.343424 0.347268 0.250206 0.289856 0.297945 0.280074 0.289085 0.187235 0.332745 0.253993 0.473024 0.246291 0.216565 0.212806 0.234784 0.203614 0.191166 0.137826 0.228107 0.171882 0.306441 0.190581 0.218887 0.166449 0.187059 0.125274 0.245567 0.19571 0.169737 0.201401 0.15554 0.134992 0.270973 0.181178 0.171978 0.169103 0.190369 0.15833 0.160691 0.224654 0.208041 0.17925 0.186553 0.122415 0.285823 0.174066 0.267729 0.195661 0.152145 0.126948 0.215444 0.202213 0.18874 0.126012 0.235587 0.175366 0.172568 0.167131 0.144522 0.104632 0.116048 0.117851 0.204357 0.126975 0.146818 0.104027 0.151904 0.133302 0.188982 0.141379 0.224158 0.17882 0.190224 0.156198 0.232219 0.166602 0.204512 0.229466 0.244229 0.239495 0.318642 0.215199 0.324162 0.324775 0.373022 0.339865 0.353118 0.334613 0.291706 0.253467 0.2736 0.198191 0.386928 0.356361 0.233413 0.139259 0.460328 0.643258 0.651489 0.569986 0.681285 0.232573 0.514805 0.755776 0.741219 0.601307 0.556145 0.525249 0.722918 0.683073 0.630482 0.540904 0.926475 0.467897 0.495716 0.453608 0.799232 0.268689 0.639434 0.644457 0.521476 0.432691 0.730311 0.642572 0.462815 0.491718 0.684612 0.602833 0.551858 0.48273 0.563112 0.38196 0.707654 0.586002 0.770349 0.674213 0.363893 0.296853 0.350582 0.302609 0.381754 0.340966 0.445208 0.410424 0.520052 0.426435 0.492048 0.413379
# NOR4a 0.0347781 0.034812 0.0163082 0.00801728 0.0504257 0.0510103 0.0239614 0.0118638 0.0634248 0.0643624 0.0303655 0.015058 0.0521307 0.0525408 0.0246389 0.0122512 0.0299578 0.0300859 0.0140907 0.00693164 0.036701 0.03712 0.0173623 0.00856552 0.047424 0.0480696 0.0225415 0.0112037 0.0444882 0.0392774 0.024338 0.00544248 0.0351947 0.0321513 0.00886325 0.0061299 0.0232937 0.0287378 0.0126182 0.00961041 0.0492813 0.0595407 0.0260531 0.0115011 0.0756934 0.0687687 0.0433559 0.0221859 0.106226 0.0871379 0.0379551 0.0268542 0.130565 0.124462 0.0420818 0.0257893 0.112205 0.0729039 0.0571072 0.0223187 0.119931 0.12584 0.0680676 0.0647607 0.129161 0.127353 0.0749183 0.0877972 0.313598 0.387902 0.0985016 0.150932 0.507153 0.472499 0.298331 0.219539 0.320972 0.296156 0.312898 0.133858 0.20543 0.260577 0.164283 0.114252 0.206 0.148169 0.126629 0.0822151 0.21089 0.154428 0.122361 0.0649464 0.195114 0.132119 0.19517 0.109191 0.107351 0.177605 0.0841837 0.0387055 0.117062 0.166777 0.0951377 0.0695798 0.14467 0.160881 0.0625812 0.0568631 0.152495 0.169478 0.105204 0.0624158 0.205962 0.169838 0.122329 0.067976 0.192119 0.146317 0.165818 0.108884 0.206389 0.177906 0.149948 0.0824302 0.280631 0.264132 0.244121 0.0619751 0.141372 0.220251 0.181343 0.112603 0.215785 0.192496 0.151789 0.101187 0.239813 0.212738 0.177197 0.10212 0.205326 0.248969 0.179404 0.0861139 0.200668 0.195616 0.196316 0.0684168 0.195401 0.19399 0.130879 0.0868344 0.140267 0.233053 0.219048 0.075966 0.212231 0.220992 0.200893 0.0828339 0.305876 0.277626 0.23837 0.186102 0.393788 0.485042 0.497653 0.108714 0.335725 0.254419 0.221221 0.137052 0.277658 0.227957 0.275916 0.109309 0.318467 0.223776 0.197624 0.122268 0.386091 0.234851 0.238032 0.170873 0.602655 0.489786 0.381623 0.354375 0.784139 0.403407 0.329765 0.239126 0.780857 0.286493 0.28986 0.224179 0.539848 0.183271 0.238876 0.114423 0.62991 0.188681 0.414675 0.290319 0.289005 0.236588 0.193883 0.115759 0.295614 0.254156 0.221808 0.138759 0.360423 0.299973 0.244527 0.143965
# SWE3a 0.0884416 0.113246 0.0979246 0.0662017 0.0945898 0.121847 0.10604 0.0720803 0.0845908 0.108522 0.0941468 0.0638168 0.0724847 0.092917 0.0802598 0.0543035 0.0688121 0.0883067 0.0763684 0.0516147 0.0604829 0.0775679 0.0670059 0.0452638 0.0839908 0.109051 0.0951818 0.0647656 0.06958 0.0893379 0.0769186 0.051608 0.0368061 0.0671211 0.0805515 0.0437431 0.0379366 0.0935148 0.0630611 0.0663244 0.0855286 0.105936 0.0913176 0.0493518 0.146245 0.0934795 0.0830145 0.0547656 0.0616942 0.0876392 0.0720072 0.0515738 0.0607563 0.0816287 0.0620474 0.0375819 0.040451 0.0467377 0.100929 0.048382 0.0406034 0.111118 0.113223 0.0876579 0.0491372 0.0989565 0.0860383 0.0641773 0.0571496 0.0866111 0.0912265 0.0882049 0.0731568 0.132447 0.106952 0.128636 0.136796 0.243822 0.268968 0.0876635 0.133143 0.139931 0.196727 0.122371 0.117387 0.138751 0.145471 0.0905094 0.156567 0.134151 0.118534 0.113726 0.0908288 0.154376 0.158568 0.12244 0.116748 0.154503 0.0994158 0.0982648 0.12023 0.171692 0.114431 0.114228 0.0962546 0.116517 0.084748 0.120577 0.067041 0.107488 0.0953064 0.0946101 0.0735143 0.104224 0.104322 0.101106 0.0523552 0.107036 0.142889 0.167293 0.0715823 0.123788 0.134539 0.146847 0.0673843 0.147243 0.116265 0.129226 0.0582483 0.118256 0.126904 0.0945516 0.0680799 0.101418 0.092146 0.108353 0.0919446 0.117079 0.113974 0.0921889 0.0767171 0.111048 0.0964635 0.121078 0.10038 0.132435 0.112401 0.130997 0.0728446 0.104925 0.150529 0.18046 0.147828 0.187577 0.205755 0.221071 0.12014 0.215745 0.206266 0.282222 0.195828 0.304391 0.306504 0.295341 0.200757 0.513032 0.402085 0.486462 0.316891 0.473212 0.451253 0.24725 0.137328 0.17983 0.197876 0.204243 0.244519 0.260344 0.220265 0.216929 0.126438 0.171752 0.182454 0.193429 0.123566 0.26136 0.100498 0.135385 0.187498 0.166898 0.18492 0.409251 0.0643645 0.293113 0.24301 0.310585 0.137599 0.282134 0.265858 0.178337 0.120715 0.411056 0.420404 0.325827 0.131782 0.203993 0.190176 0.176702 0.136481 0.228456 0.236418 0.237616 0.188307 0.297062 0.272418 0.246251
# SWE4a 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.0223803 0.00884743 0.00362609 0.00109751 0.0127914 0.00513233 0.00212672 0.000557726 0.0110056 0.00432599 0.00185975 0.000561843 0.00709436 0.00173906 0.00172352 0.000986758 0.00397996 0.000433966 0.000716955 0.000282091 0.00400882 0.00132837 0 0.000146254 0.00284623 0.001411 0.000638095 0 0.00500003 0.00626257 0.00148769 0 0.00913066 0.00698812 0.00260037 0.00302911 0.00446009 0.00992575 0.00249096 0.00369697 0.00661202 0.0116419 0.00494598 0.000231952 0.0132714 0.0121175 0 0.00343901 0.0238538 0.0227515 0.00578858 0.0040357 0.0278515 0.0251414 0.00587456 0.00484504 0.020634 0.00817877 0.00621802 0.000896468 0.0233677 0.0069985 0.000259652 0.000995437 0.0319871 0.00223007 0 0 0.0269113 0 0 0 0.0298662 0 0 0 0.0113299 0.00850117 0.00147145 0.00158531 0.0117 0.00886983 0.00157254 0.00169266 0 0.000174357 0.00720755 0.00868037 0.0383798 0.0100359 0 0.00623274 0.0314989 0.00191056 0.000200938 0 0.0408741 0.00319185 0 0 0.0412059 0.0036594 0.00349938 0.00140263 0.0337765 0 0.00159177 0 0.0276436 0.00300629 0 0 0.0262454 0.000173117 0 0 0.0179657 0.00820503 0 0.00114609 0.0117036 0.00799198 0.00399129 0.00299047 0.0222047 0.0141001 0 0 0.027205 0.0134311 0.000264629 0 0.028242 0.017993 0 0 0.0630477 0.00202549 0 0 0.0883238 0.00112884 0.000601565 0.00250924 0.133774 0 0 0 0.149381 0.00103177 0.00189316 0.00295397 0.0619188 0.00875495 0.00077572 0.00221322 0.0181783 0.0253941 0.014992 0.0230903 0.0445361 0.00966926 0.0133871 0.000612553 0.0808863 0.0113749 0.00161023 0.00249921 0.102909 0.004451 0.00880591 0.00834828 0.109582 0.00724209 0 0.00116949 0.0497864 0.0196211 0.00507904 0.0141938 0.0300611 0.00351545 0.0061464 0.00769915 0.0402493 0.00838015 0.00287584 0.00290995 0.0421883 0.00939509 0.00341794 0.00355498 0.0505215 0.010517 0.00357213 0.0035543
#
#_Q_setup for fleets with cpue or survey data
#_1:  fleet number
#_2:  link type: (1=simple q, 1 parm; 2=mirror simple q, 1 mirrored parm; 3=q and power, 2 parm; 4=mirror with offset, 2 parm)
#_3:  extra input for link, i.e. mirror fleet# or dev index number
#_4:  0/1 to select extra sd parameter
#_5:  0/1 for biasadj or not
#_6:  0/1 to float
#_   fleet      link link_info  extra_se   biasadj     float  #  fleetname
         7         1         0         1         0         1  #  NORSURVEY3a
         8         1         0         1         0         1  #  NORSURVEY4a
-9999 0 0 0 0 0
#
#_Q_parms(if_any);Qunits_are_ln(q)
#_          LO            HI          INIT         PRIOR         PR_SD       PR_type      PHASE    env-var    use_dev   dev_mnyr   dev_mxyr     dev_PH      Block    Blk_Fxn  #  parm_name
           -15            15      -2.39701             0             1             0         -1          0          0          0          0          0          0          0  #  LnQ_base_NORSURVEY3a(7)
             0             1             0        0.0001          0.01             0         -4          0          0          0          0          0          0          0  #  Q_extraSD_NORSURVEY3a(7)
           -15            15      -1.45302             0             1             0         -1          0          0          0          0          0          0          0  #  LnQ_base_NORSURVEY4a(8)
             0             1             0             0          0.01             0         -4          0          0          0          0          0          0          0  #  Q_extraSD_NORSURVEY4a(8)
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
 1 1 0 0 # 1 DK3a
 1 1 0 0 # 2 DK4a
 1 1 0 0 # 3 NOR3a
 1 1 0 0 # 4 NOR4a
 1 1 0 0 # 5 SWE3a
 1 1 0 0 # 6 SWE4a
 1 0 0 0 # 7 NORSURVEY3a
 1 0 0 0 # 8 NORSURVEY4a
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
 0 0 0 0 # 1 DK3a
 0 0 0 0 # 2 DK4a
 0 0 0 0 # 3 NOR3a
 0 0 0 0 # 4 NOR4a
 0 0 0 0 # 5 SWE3a
 0 0 0 0 # 6 SWE4a
 0 0 0 0 # 7 NORSURVEY3a
 0 0 0 0 # 8 NORSURVEY4a
#
#_          LO            HI          INIT         PRIOR         PR_SD       PR_type      PHASE    env-var    use_dev   dev_mnyr   dev_mxyr     dev_PH      Block    Blk_Fxn  #  parm_name
# 1   DK3a LenSelex
          0.35          3.45       1.61594           0.1            10             0          5          0          0          0          0          0          0          0  #  Size_inflection_DK3a(1)
           0.1          3.45      0.362517           0.1            10             0          6          0          0          0          0          0          0          0  #  Size_95%width_DK3a(1)
           0.2            30       1.11612            27          0.01             0          5          0         23       2019       2024          6          0          0  #  Retain_L_infl_DK3a(1)
            -5            30      0.239188           0.8          0.01             0          6          0          0          0          0          0          0          0  #  Retain_L_width_DK3a(1)
           -10           999           999            10          0.01             0         -6          0          0          0          0          0          0          0  #  Retain_L_asymptote_logit_DK3a(1)
             0             0             0             0          0.01             0         -6          0          0          0          0          0          0          0  #  Retain_L_maleoffset_DK3a(1)
# 2   DK4a LenSelex
          0.35          3.45       1.52137           0.1            10             0          5          0          0          0          0          0          0          0  #  Size_inflection_DK4a(2)
           0.1           4.5      0.356558           0.1            10             0          6          0          0          0          0          0          0          0  #  Size_95%width_DK4a(2)
           0.2            30       1.24319            27          0.01             0          5          0         23       2019       2024          6          0          0  #  Retain_L_infl_DK4a(2)
            -5            30      0.260519           0.8          0.01             0          6          0          0          0          0          0          0          0  #  Retain_L_width_DK4a(2)
           -10           999           999            10          0.01             0         -6          0          0          0          0          0          0          0  #  Retain_L_asymptote_logit_DK4a(2)
             0             0             0             0          0.01             0         -6          0          0          0          0          0          0          0  #  Retain_L_maleoffset_DK4a(2)
# 3   NOR3a LenSelex
          0.35          3.45       1.66415           0.1            10             0          5          0          0          0          0          0          0          0  #  Size_inflection_NOR3a(3)
           0.1           4.5       0.36529           0.1            10             0          6          0          0          0          0          0          0          0  #  Size_95%width_NOR3a(3)
           0.2            30       1.30302            27          0.01             0          5          0         23       2019       2024          6          0          0  #  Retain_L_infl_NOR3a(3)
          -5.1            30      0.143753           0.8          0.01             0          6          0          0          0          0          0          0          0  #  Retain_L_width_NOR3a(3)
           -10           999           999            10          0.01             0         -6          0          0          0          0          0          0          0  #  Retain_L_asymptote_logit_NOR3a(3)
             0             0             0             0          0.01             0         -6          0          0          0          0          0          0          0  #  Retain_L_maleoffset_NOR3a(3)
# 4   NOR4a LenSelex
          0.35          3.45        1.6294           0.1            10             0          5          0          0          0          0          0          0          0  #  Size_inflection_NOR4a(4)
           0.1           4.5      0.363968           0.1            10             0          6          0          0          0          0          0          0          0  #  Size_95%width_NOR4a(4)
           0.2            30      0.798612            27          0.01             0          5          0          0          0          0          0          0          0  #  Retain_L_infl_NOR4a(4)
          -5.1            30      0.240279           0.8          0.01             0          6          0          0          0          0          0          0          0  #  Retain_L_width_NOR4a(4)
           -10           999           999            10          0.01             0         -6          0          0          0          0          0          0          0  #  Retain_L_asymptote_logit_NOR4a(4)
             0             0             0             0          0.01             0         -6          0          0          0          0          0          0          0  #  Retain_L_maleoffset_NOR4a(4)
# 5   SWE3a LenSelex
          0.35          3.45       1.73129           0.1            10             0          5          0          0          0          0          0          0          0  #  Size_inflection_SWE3a(5)
           0.1           4.5      0.370335           0.1            10             0          6          0          0          0          0          0          0          0  #  Size_95%width_SWE3a(5)
           0.2            30       1.65231            27          0.01             0          5          0         23       2019       2024          6          0          0  #  Retain_L_infl_SWE3a(5)
          -5.1            30      0.144021           0.8          0.01             0          6          0          0          0          0          0          0          0  #  Retain_L_width_SWE3a(5)
           -10           999           999            10          0.01             0         -6          0          0          0          0          0          0          0  #  Retain_L_asymptote_logit_SWE3a(5)
             0             0             0             0          0.01             0         -6          0          0          0          0          0          0          0  #  Retain_L_maleoffset_SWE3a(5)
# 6   SWE4a LenSelex
          0.35          3.45       1.47946           0.1            10             0          5          0          0          0          0          0          0          0  #  Size_inflection_SWE4a(6)
           0.1           4.5       0.33725           0.1            10             0          6          0          0          0          0          0          0          0  #  Size_95%width_SWE4a(6)
           0.2            30       1.66808            27          0.01             0          5          0         23       2019       2024          6          0          0  #  Retain_L_infl_SWE4a(6)
          -5.1            30      0.109142           0.8          0.01             0          6          0          0          0          0          0          0          0  #  Retain_L_width_SWE4a(6)
           -10           999           999            10          0.01             0         -6          0          0          0          0          0          0          0  #  Retain_L_asymptote_logit_SWE4a(6)
             0             0             0             0          0.01             0         -6          0          0          0          0          0          0          0  #  Retain_L_maleoffset_SWE4a(6)
# 7   NORSURVEY3a LenSelex
          0.35          3.45       1.54924           0.1            10             0          5          0          0          0          0          0          0          0  #  Size_inflection_NORSURVEY3a(7)
            -5           4.5      0.535606           0.1            10             0          6          0          0          0          0          0          0          0  #  Size_95%width_NORSURVEY3a(7)
# 8   NORSURVEY4a LenSelex
          0.35          3.45       1.62084           0.1            10             0          5          0          0          0          0          0          0          0  #  Size_inflection_NORSURVEY4a(8)
           0.1           4.5       0.43794           0.1            10             0          6          0          0          0          0          0          0          0  #  Size_95%width_NORSURVEY4a(8)
# 1   DK3a AgeSelex
# 2   DK4a AgeSelex
# 3   NOR3a AgeSelex
# 4   NOR4a AgeSelex
# 5   SWE3a AgeSelex
# 6   SWE4a AgeSelex
# 7   NORSURVEY3a AgeSelex
# 8   NORSURVEY4a AgeSelex
#_No_Dirichlet parameters
# timevary selex parameters 
#_          LO            HI          INIT         PRIOR         PR_SD       PR_type    PHASE  #  parm_name
           0.5             2           1.5             1           0.5             6      -7  # Retain_L_infl_DK3a(1)_dev_se
         -0.99          0.99             0             0           0.5             6      -6  # Retain_L_infl_DK3a(1)_dev_autocorr
           0.5             2           1.5             1           0.5             6      -7  # Retain_L_infl_DK4a(2)_dev_se
         -0.99          0.99             0             0           0.5             6      -6  # Retain_L_infl_DK4a(2)_dev_autocorr
           0.5             2           1.5             1           0.5             6      -7  # Retain_L_infl_NOR3a(3)_dev_se
         -0.99          0.99             0             0           0.5             6      -6  # Retain_L_infl_NOR3a(3)_dev_autocorr
           0.5             2           1.5             1           0.5             6      -7  # Retain_L_infl_SWE3a(5)_dev_se
         -0.99          0.99             0             0           0.5             6      -6  # Retain_L_infl_SWE3a(5)_dev_autocorr
           0.5             2           1.5             1           0.5             6      -7  # Retain_L_infl_SWE4a(6)_dev_se
         -0.99          0.99             0             0           0.5             6      -6  # Retain_L_infl_SWE4a(6)_dev_autocorr
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
#      1    31     1     0     0     0     0     1     3  1984  2020     3 -0.349682 0.149717 0.60587 -0.155229 -0.0264751 -0.504616 0.331032 -0.0753733 -0.223275 0.250467 -0.0781108 -0.271319 0.070237 -0.0183626 0.158613 0.299994 -0.27664 0.094737 0.0711257 -0.0245331 -0.392618 -0.103155 0.082669 -0.194251 0.257721 0.157828 -0.45181 0.176877 -0.0174566 -0.096235 -0.131472 0.297178 -0.100503 0.02426 0.532521 -0.230665 -0.141003
#      5     3     3     0     0     0     0     2     3  2017  2020     6 -0.159553 -0.130948 0.25647 -0.0949946
#      5     9     5     0     0     0     0     3     3  2017  2020     6 0.130119 -0.0183394 -0.33782 -0.0138783
#      5    15     7     0     0     0     0     4     3  2017  2020     6 -0.0694148 0.0319803 0.0244834 -0.0976033
#      5    27     9     0     0     0     0     5     3  2017  2020     6 -0.179095 -0.0691265 0.152638 -0.0124351
#      5    33    11     0     0     0     0     6     3  2017  2020     6 0.0446424 0.0144834 -0.173222 0.0266183
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
 -9999   1    0  # terminator
#
4 #_maxlambdaphase
1 #_sd_offset; must be 1 if any growthCV, sigmaR, or survey extraSD is an estimated parameter
# read 0 changes to default Lambdas (default value is 1.0)
# Like_comp codes:  1=surv; 2=disc; 3=mnwt; 4=length; 5=age; 6=SizeFreq; 7=sizeage; 8=catch; 9=init_equ_catch; 
# 10=recrdev; 11=parm_prior; 12=parm_dev; 13=CrashPen; 14=Morphcomp; 15=Tag-comp; 16=Tag-negbin; 17=F_ballpark; 18=initEQregime
#like_comp fleet  phase  value  sizefreq_method
-9999  1  1  1  1  #  terminator
#
# lambdas (for info only; columns are phases)
#  0 0 0 0 #_CPUE/survey:_1
#  0 0 0 0 #_CPUE/survey:_2
#  0 0 0 0 #_CPUE/survey:_3
#  0 0 0 0 #_CPUE/survey:_4
#  0 0 0 0 #_CPUE/survey:_5
#  0 0 0 0 #_CPUE/survey:_6
#  1 1 1 1 #_CPUE/survey:_7
#  1 1 1 1 #_CPUE/survey:_8
#  1 1 1 1 #_discard:_1
#  1 1 1 1 #_discard:_2
#  1 1 1 1 #_discard:_3
#  1 1 1 1 #_discard:_4
#  1 1 1 1 #_discard:_5
#  1 1 1 1 #_discard:_6
#  0 0 0 0 #_discard:_7
#  0 0 0 0 #_discard:_8
#  1 1 1 1 #_lencomp:_1
#  1 1 1 1 #_lencomp:_2
#  1 1 1 1 #_lencomp:_3
#  1 1 1 1 #_lencomp:_4
#  1 1 1 1 #_lencomp:_5
#  1 1 1 1 #_lencomp:_6
#  1 1 1 1 #_lencomp:_7
#  1 1 1 1 #_lencomp:_8
#  1 1 1 1 #_init_equ_catch1
#  1 1 1 1 #_init_equ_catch2
#  1 1 1 1 #_init_equ_catch3
#  1 1 1 1 #_init_equ_catch4
#  1 1 1 1 #_init_equ_catch5
#  1 1 1 1 #_init_equ_catch6
#  1 1 1 1 #_init_equ_catch7
#  1 1 1 1 #_init_equ_catch8
#  1 1 1 1 #_recruitments
#  1 1 1 1 #_parameter-priors
#  1 1 1 1 #_parameter-dev-vectors
#  1 1 1 1 #_crashPenLambda
#  0 0 0 0 # F_ballpark_lambda
0 # (0/1/2) read specs for more stddev reporting: 0 = skip, 1 = read specs for reporting stdev for selectivity, size, and numbers, 2 = add options for M,Dyn. Bzero, SmryBio
 # 0 2 0 0 # Selectivity: (1) fleet, (2) 1=len/2=age/3=both, (3) year, (4) N selex bins
 # 0 0 # Growth: (1) growth pattern, (2) growth ages
 # 0 0 0 # Numbers-at-age: (1) area(-1 for all), (2) year, (3) N ages
 # -1 # list of bin #'s for selex std (-1 in first bin to self-generate)
 # -1 # list of ages for growth std (-1 in first bin to self-generate)
 # -1 # list of ages for NatAge std (-1 in first bin to self-generate)
999