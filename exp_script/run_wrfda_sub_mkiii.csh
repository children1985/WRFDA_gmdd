#!/bin/csh
#
### Project name
#PBS -A NMMM0015

### Job name
#PBS -N wrfda

### Wallclock time
#PBS -l walltime=01:00:00

### Queue
#PBS -q regular

### Merge output and error files
#PBS -j oe                    

### Select 2 nodes with 36 CPUs, for 72 MPI processes 
#PBS -l select=4:ncpus=36:mpiprocs=36:ompthreads=30  

set echo
############case dependent######################

set DATE         = AADATE #2013122306

set RUN_TYPE     = cycle   # cold or cycle
set DOMAIN_ID    = 01
set NPROC        = 36

set DA_SRCDIR    = AADA_SRCDIR #/glade/u/home/wangs/WRFDA
set BIN_DIR      = ${DA_SRCDIR}/var/build
set OB_DATDIR    = AAOB_DATDIR  #/glade/work/wangs/wrfda_test/cycling/ob
set BE_DATDIR    = AABE_DATDIR  #/glade/work/wangs/wrfda_test/cycling/be
set RC_DATDIR    = AARC_DATDIR  #/glade/work/wangs/wrfda_test/cycling/rc
set RUN_BASEDIR  = AARUN_BASEDIR  #/glade/work/wangs/wrfda_test/cycling/run
set DA_RUNDIR    = ${RUN_BASEDIR}/${DATE}/da/d${DOMAIN_ID}

set TIMEWINDOW1  =  -10m
set TIMEWINDOW2  =  10m
set CYCLE_PERIOD =  AACYCLE_PERIODm

set nouterloop = AAnouterloop
set ninnerloop = AAninnerloop
set cloud_cv_options = AAcloud_cv_options
set use_cv_w = AAuse_cv_w
set twc_opt = AAtwc_opt
set radar_rf_opt = AAradar_rf_opt
set use_radar_rv = AAuse_radar_rv
set use_radar_rf = AAuse_radar_rf
set use_radar_rqv = AAuse_radar_rqv
set use_radar_rhv = AAuse_radar_rhv
set rf_noice = AArf_noice
set preproc = AApreproc
set uvlenscl = AAuvlenscl
set tlenscl = AAtlenscl
set qvlenscl = AAqvlenscl
set hydrolenscl = AAhydrolenscl

set nx1 = 450
set ny1 = 450
set nz1 = 43
set nx2 = 450
set ny2 = 450 
set nz2 = 43

set dx1 = 3000
set dy1 = 3000
set dx2 = 1000
set dy2 = 1000
##############################################################
set mphyopt = 7
set cuopt = 0
##############################################################

set var_scaling1 = $uvlenscl
set var_scaling2 = $uvlenscl
set var_scaling3 = $tlenscl
set var_scaling4 = $qvlenscl
set var_scaling5 = 1.0
set var_scaling6 = $hydrolenscl
set var_scaling7 = $hydrolenscl
set var_scaling8 = $hydrolenscl
set var_scaling9 = $hydrolenscl
set var_scaling10 = $hydrolenscl
set var_scaling11 = 1.0
set len_scaling1 = 1.0
set len_scaling2 = 1.0
set len_scaling3 = 1.0
set len_scaling4 = 1.0
set len_scaling5 = 1.0
set len_scaling6 = 1.0
set len_scaling7 = 1.0
set len_scaling8 = 1.0
set len_scaling9 = 1.0
set len_scaling10 = 1.0
set len_scaling11 = 1.0

##############################################################

set ccyy = `echo $DATE | cut -c1-4`
set   mm = `echo $DATE | cut -c5-6`
set   dd = `echo $DATE | cut -c7-8`
set   hh = `echo $DATE | cut -c9-10`
set   nn = `echo $DATE | cut -c11-12`

set DATE_MIN = `${BIN_DIR}/da_advance_time.exe ${DATE} ${TIMEWINDOW1} -w`
set DATE_MAX = `${BIN_DIR}/da_advance_time.exe ${DATE} ${TIMEWINDOW2} -w`

set PREV_DATE = `${BIN_DIR}/da_advance_time.exe ${DATE} -${CYCLE_PERIOD} -f ccyymmddhhnn`
set VARBC_PREV_DATE = `${BIN_DIR}/da_advance_time.exe ${DATE} -${CYCLE_PERIOD} -f ccyymmddhhnn`

if ( ! -d ${DA_RUNDIR} ) mkdir -p ${DA_RUNDIR}
cd ${DA_RUNDIR}

if ( ${DOMAIN_ID} == '01' ) then
   set WEST_EAST_GRID_NUMBER   = $nx1
   set SOUTH_NORTH_GRID_NUMBER = $ny1
   set VERTICAL_GRID_NUMBER    = $nz1
   set GRID_DISTANCE           = $dx1
   set BE_FILE = be.dat
else if ( ${DOMAIN_ID} == '02' ) then
   set WEST_EAST_GRID_NUMBER   = $nx2
   set SOUTH_NORTH_GRID_NUMBER = $ny2
   set VERTICAL_GRID_NUMBER    = $nz2
   set GRID_DISTANCE           = $dx2
   set BE_FILE = be_d02.dat
endif

rm -f ${DA_RUNDIR}/rsl.*
#
# link some constant files
#
ln -fs ${DA_SRCDIR}/run/LANDUSE.TBL  .
#
# link first-guess, observations, background-error
#
# link either ob.bufr (for ob_format=1) or ob.ascii (for ob_format=2)
#ln -fs ${OB_DATDIR}/${DATE}/prepbufr.gdas.${ccyy}${mm}${dd}.t${hh}z.nr    ./ob.bufr
#ln -fs ${OB_DATDIR}/${DATE}/ob.ascii                           ./ob.ascii
# BE file for cv_options=5
ln -fs ${BE_DATDIR}/${BE_FILE}                                 ./be.dat
# BE file for cv_options=3
#ln -fs ${DA_SRCDIR}/var/run/be.dat.cv3                      ./be.dat
#
# link radiance-related files
#
#ln -fs ${OB_DATDIR}/${DATE}/gdas.1bamua.t${hh}z.${ccyy}${mm}${dd}.bufr   ./amsua.bufr
#ln -fs ${OB_DATDIR}/${DATE}/gdas.1bamub.t${hh}z.${ccyy}${mm}${dd}.bufr   ./amsub.bufr
ln -fs ${DA_SRCDIR}/var/run/radiance_info .
# link crtm_coeffs if using CRTM
ln -fs ${DA_SRCDIR}/var/run/crtm_coeffs   .
# link rttov_coeffs if using RTTOV
#ln -fs /kumquat/wrfhelp/external/rttov11/rtcoef_rttov11/rttov7pred54L ./rttov_coeffs
#
# link radar data
#ln -sf ${OB_DATDIR}/AAradardir/ob.radar.${ccyy}${mm}${dd}${hh}${nn}00 ./ob.radar
ln -sf AAradardir/ob.radar.${ccyy}${mm}${dd}${hh}${nn}00 ./ob.radar

set VARBC_DIR = ${RUN_BASEDIR}/${VARBC_PREV_DATE}/da/d${DOMAIN_ID}
if ( ! -e ${VARBC_DIR}/VARBC.out ) then
   ln -fs ${DA_SRCDIR}/var/run/VARBC.in ./VARBC.in
else
   ln -fs ${VARBC_DIR}/VARBC.out  ./VARBC.in
endif

if ( ${RUN_TYPE} == 'cold' ) then
   if ( ! -e ${RC_DATDIR}/${DATE}/wrfinput_d${DOMAIN_ID} ) then
      echo "ERROR in run_wrfda.csh : first guess ${RC_DATDIR}/${DATE}/wrfinput_d${DOMAIN_ID} not found..."
      exit 1
   endif
   ln -fs ${RC_DATDIR}/${DATE}/wrfinput_d${DOMAIN_ID} fg
else  # cycling
   if ( ! -e ${RUN_BASEDIR}/${PREV_DATE}/wrf/wrfout_d${DOMAIN_ID}_${ccyy}-${mm}-${dd}_${hh}:${nn}:00 ) then
      echo "ERROR in run_wrfda.csh : first guess ${RUN_BASEDIR}/${PREV_DATE}/wrf/wrfout_d${DOMAIN_ID}_${ccyy}-${mm}-${dd}_${hh}:${nn}:00 not found..."
      exit 1
   else
     if( AAidastep == 1 ) then
        ln -fs ${RUN_BASEDIR}/${PREV_DATE}/wrf/wrfout_d${DOMAIN_ID}_${ccyy}-${mm}-${dd}_${hh}:${nn}:00 fg_orig
        cp -p ${RUN_BASEDIR}/${PREV_DATE}/wrf/wrfout_d${DOMAIN_ID}_${ccyy}-${mm}-${dd}_${hh}:${nn}:00 fg
     else
        cp fg fg_AAidastep_bak
        rm fg
        cp -f wrfvar_output fg
     endif
   endif
endif

if ( ${RUN_TYPE} == 'cycle' ) then
   # when not cold-starting, update lower boundary first, before running wrfvar
   # fields to be updated are TSK over water, TMN, SST, VEGFRA, ALBBCK, SEAICE, LANDMASK, IVGTYP, ISLTYP
   cd ${DA_RUNDIR}
   cat >! ${DA_RUNDIR}/parame.in << EOF
&control_param
 da_file            = '${DA_RUNDIR}/fg'
! wrf_input          = '${RC_DATDIR}/${DATE}/wrfinput_d${DOMAIN_ID}'
 wrf_input          = '${RUN_BASEDIR}/${PREV_DATE}/wrf/wrfinput_d${DOMAIN_ID}'
 domain_id          = ${DOMAIN_ID}
 debug              = .false.
 update_lateral_bdy = .false.
 update_low_bdy     = .true.
 update_lsm         = .false.
 iswater            = 16 /
EOF
   ln -fs ${DA_SRCDIR}/var/build/da_update_bc.exe .
   time ./da_update_bc.exe >&! update_low_bc_${DATE}.log
   mv parame.in parame.in.lowbdy
   # check status
   grep "Update_bc completed successfully" update_low_bc_${DATE}.log
   if ( $status != 0 ) then
      echo "ERROR in run_wrfda.csh : update low bdy failed..."
      exit 1
   endif
endif
#
# create namelist for WRFDA run
#
cat >! namelist.input << EOF
 &wrfvar1
 PRINT_DETAIL_RAD =  F,
 PRINT_DETAIL_XA  =  F,
 PRINT_DETAIL_XB  =  F,
 PRINT_DETAIL_OBS =  F,
 PRINT_DETAIL_F_OBS =  F,
 PRINT_DETAIL_MAP   =  F,
 PRINT_DETAIL_GRAD  =  F,
 PRINT_DETAIL_TESTING =  F,
 CHECK_MAX_IV_PRINT   =  F,
 print_detail_radar=.true.,
 print_detail_grad=false
 /
 &wrfvar2
 /
 &wrfvar3
 ob_format=1
 /
 &wrfvar4
 thin_mesh_conv = 28*180.0
 USE_SYNOPOBS =  F,
 USE_SHIPSOBS =  F,
 USE_METAROBS =  F,
 USE_SOUNDOBS =  F,
 USE_MTGIRSOBS =  F,
 USE_TAMDAROBS =  F,
 USE_PILOTOBS =  F,
 USE_AIREPOBS =  F,
 USE_GEOAMVOBS =  F,
 USE_POLARAMVOBS =  F,
 USE_BOGUSOBS =  F,
 USE_BUOYOBS =  F,
 USE_PROFILEROBS =  F,
 USE_SATEMOBS =  F,
 USE_GPSZTDOBS =  F,
 USE_GPSPWOBS =  F,
 USE_GPSREFOBS =  F,
 USE_QSCATOBS =  F,
 USE_AIRSRETOBS =  F,
 use_ssmiretrievalobs=false
 use_amsuaobs = F,
 use_amsubobs = F,
 use_mhsobs   = F,
 use_airsobs  = F,
 use_eos_amsuaobs = F,
 USE_OBS_ERRFAC=F
 use_radarobs = .true.
 use_radar_rv = $use_radar_rv
 use_radar_rf = $use_radar_rf
 radar_rf_opt = $radar_rf_opt
 twc_opt=$twc_opt  ! total water conservation constraint
 radar_rf_rscl = 1.0
 radar_rv_rscl = 1.0
 rf_noice = $rf_noice
 use_radar_rqv = $use_radar_rqv
 use_radar_rhv = $use_radar_rhv
! use_3dvar_phy = .true.
 /
 &wrfvar5
 CHECK_MAX_IV    = T,
 MAX_ERROR_T =    5.0     ,
 MAX_ERROR_UV =    5.0     ,
 MAX_ERROR_PW =    5.0     ,
 MAX_ERROR_REF =    5.0     ,
 MAX_ERROR_Q =    5.0     ,
 MAX_ERROR_P =    5.0   ,
 max_error_rv= 5.0,
 max_error_rf= 50.0,
 /
 &wrfvar6
 max_ext_its=$nouterloop,
 ntmax=$ninnerloop,$ninnerloop,$ninnerloop,$ninnerloop,$ninnerloop,$ninnerloop,$ninnerloop,$ninnerloop,$ninnerloop,$ninnerloop,$ninnerloop
 eps=0.01,
 orthonorm_gradient=.true.,
 /
 &wrfvar7
 cv_options=7
 USE_RF  = T,
 RF_PASSES =  2,
 var_scaling1 = $var_scaling1,$var_scaling1,$var_scaling1,$var_scaling1,$var_scaling1,$var_scaling1,$var_scaling1,$var_scaling1,$var_scaling1
 var_scaling2 = $var_scaling2,$var_scaling2,$var_scaling2,$var_scaling2,$var_scaling2,$var_scaling2,$var_scaling2,$var_scaling2,$var_scaling2
 var_scaling3 = $var_scaling3,$var_scaling3,$var_scaling3,$var_scaling3,$var_scaling3,$var_scaling3,$var_scaling3,$var_scaling3,$var_scaling3
 var_scaling4 = $var_scaling4,$var_scaling4,$var_scaling4,$var_scaling4,$var_scaling4,$var_scaling4,$var_scaling4,$var_scaling4,$var_scaling4
 var_scaling5 = $var_scaling5,$var_scaling5,$var_scaling5,$var_scaling5,$var_scaling5,$var_scaling5,$var_scaling5,$var_scaling5,$var_scaling5
 var_scaling6 = $var_scaling6,$var_scaling6,$var_scaling6,$var_scaling6,$var_scaling6,$var_scaling6,$var_scaling6,$var_scaling6,$var_scaling6
 var_scaling7 = $var_scaling7,$var_scaling7,$var_scaling7,$var_scaling7,$var_scaling7,$var_scaling7,$var_scaling7,$var_scaling7,$var_scaling7
 var_scaling8 = $var_scaling8,$var_scaling8,$var_scaling8,$var_scaling8,$var_scaling8,$var_scaling8,$var_scaling8,$var_scaling8,$var_scaling8
 var_scaling9 = $var_scaling9,$var_scaling9,$var_scaling9,$var_scaling9,$var_scaling9,$var_scaling9,$var_scaling9,$var_scaling9,$var_scaling9
 var_scaling10 = $var_scaling10,$var_scaling10,$var_scaling10,$var_scaling10,$var_scaling10,$var_scaling10,$var_scaling10,$var_scaling10,$var_scaling10 
 var_scaling11 = $var_scaling11,$var_scaling11,$var_scaling11,$var_scaling11,$var_scaling11,$var_scaling11,$var_scaling11,$var_scaling11,$var_scaling11
 len_scaling1 = $len_scaling1,$len_scaling1,$len_scaling1,$len_scaling1,$len_scaling1,$len_scaling1,$len_scaling1,$len_scaling1,$len_scaling1
 len_scaling2 = $len_scaling2,$len_scaling2,$len_scaling2,$len_scaling2,$len_scaling2,$len_scaling2,$len_scaling2,$len_scaling2,$len_scaling2
 len_scaling3 = $len_scaling3,$len_scaling3,$len_scaling3,$len_scaling3,$len_scaling3,$len_scaling3,$len_scaling3,$len_scaling3,$len_scaling3
 len_scaling4 = $len_scaling4,$len_scaling4,$len_scaling4,$len_scaling4,$len_scaling4,$len_scaling4,$len_scaling4,$len_scaling4,$len_scaling4
 len_scaling5 = $len_scaling5,$len_scaling5,$len_scaling5,$len_scaling5,$len_scaling5,$len_scaling5,$len_scaling5,$len_scaling5,$len_scaling5
 len_scaling6 = $len_scaling6,$len_scaling6,$len_scaling6,$len_scaling6,$len_scaling6,$len_scaling6,$len_scaling6,$len_scaling6,$len_scaling6
 len_scaling7 = $len_scaling7,$len_scaling7,$len_scaling7,$len_scaling7,$len_scaling7,$len_scaling7,$len_scaling7,$len_scaling7,$len_scaling7
 len_scaling8 = $len_scaling8,$len_scaling8,$len_scaling8,$len_scaling8,$len_scaling8,$len_scaling8,$len_scaling8,$len_scaling8,$len_scaling8
 len_scaling9 = $len_scaling9,$len_scaling9,$len_scaling9,$len_scaling9,$len_scaling9,$len_scaling9,$len_scaling9,$len_scaling9,$len_scaling9
 len_scaling10 = $len_scaling10,$len_scaling10,$len_scaling10,$len_scaling10,$len_scaling10,$len_scaling10,$len_scaling10,$len_scaling10,$len_scaling10
 len_scaling11 = $len_scaling11,$len_scaling11,$len_scaling11,$len_scaling11,$len_scaling11,$len_scaling11,$len_scaling11,$len_scaling11,$len_scaling11
 cloud_cv_options=$cloud_cv_options,
 use_cv_w=$use_cv_w, ! used when cloud_cv_options=2
 /
 &wrfvar8
 /
 &wrfvar9
 /
 &wrfvar10
 test_transforms=false,
 test_gradient=false,
 /
 &wrfvar11
 calculate_cg_cost_fn = .true.
 SFC_ASSI_OPTIONS =      2,
 check_rh=1
 /
 &wrfvar12
 /
 &wrfvar13
 /
 &wrfvar14
 rtminit_nsensor=5,
 rtminit_platform=1,1,10,1,1,
 rtminit_satid=16,18,2,16,17,
 rtminit_sensor=3,3,3,4,4,
 thinning_mesh=30*200.0
 thinning=true,
 qc_rad=true,
 write_iv_rad_ascii=false,
 write_oa_rad_ascii=true,
 rtm_option=2,
 only_sea_rad=false,
 use_varbc=.true.
 use_crtm_kmatrix=.true.
 /
 &wrfvar15
 /
 &wrfvar16
 /
 &wrfvar17
! analysis_type="QC-OBS"
  analysis_type="3D-VAR",
 /
 &wrfvar18
 analysis_date="${ccyy}-${mm}-${dd}_${hh}:${nn}:00",
 /
 &wrfvar19
 /
 &wrfvar20
 /
 &wrfvar21
 time_window_min="${DATE_MIN}",
 /
 &wrfvar22
 time_window_max="${DATE_MAX}",
 /
 &wrfvar23
 /
 &time_control
 start_year                          = ${ccyy}
 start_month                         = ${mm}
 start_day                           = ${dd}
 start_hour                          = ${hh}
 start_minute                        = 00
 start_second                        = 00
 end_year                            = ${ccyy}
 end_month                           = ${mm}
 end_day                             = ${dd}
 end_hour                            = ${hh}
 end_minute                          = 00
 end_second                          = 00
 /
 &dfi_control
 /
 &domains
 s_we                                = 1
 e_we                                = ${WEST_EAST_GRID_NUMBER}
 s_sn                                = 1
 e_sn                                = ${SOUTH_NORTH_GRID_NUMBER}
 s_vert                              = 1
 e_vert                              = ${VERTICAL_GRID_NUMBER}
 dx                                  = ${GRID_DISTANCE}
 dy                                  = ${GRID_DISTANCE}
 /
 &physics
 mp_physics                          = ${mphyopt},
 mp_physics_ad                       =  99,
 sf_sfclay_physics                   = 1,
 sf_surface_physics                  = 2,
 bl_pbl_physics                      = 1,
 cu_physics                          = ${cuopt},
 num_soil_layers                     = 4,
 num_land_cat                        = 21,
 num_soil_cat                        = 16,
 /
 &fdda
 /
 &dynamics
 /
 &bdy_control
 /
 &grib2
 /
 &namelist_quilt
 nio_tasks_per_group = 0,
 nio_groups = 1,
 /
EOF

ln -sf ${DA_SRCDIR}/var/build/da_wrfvar.exe .
if ( $NPROC > 1 ) then
#  cp /glade/u/home/wangs/test/WRFDA/tools/testdbz/bak_be/be_dvar_ascii.dat ./
  if( $preproc == 1 ) then
    cp /glade/u/home/wangs/test/WRFDA/tools/testdbz/rdr2wrf.exe ./
    setenv OMP_NUM_THREADS 30
    ./rdr2wrf.exe
  endif
#   mpirun -np ${NPROC} ./da_wrfvar.exe  >&! wrfda_${DATE}_d${DOMAIN_ID}.log
  mpiexec_mpt ./da_wrfvar.exe  >&! wrfda_${DATE}_d${DOMAIN_ID}.log
else
   time ./da_wrfvar.exe  >&! wrfda_${DATE}_d${DOMAIN_ID}.log
endif

# check status
if ( $NPROC > 1 ) then
   grep "WRF-Var completed successfully" rsl.out.0000
else
   grep "WRF-Var completed successfully" wrfda_${DATE}_d${DOMAIN_ID}.log
endif
if ( $status != 0 ) then
   echo "ERROR in run_wrfda.csh : da_wrfvar.exe failed..."
   exit 1
endif

rm gts_omb_oma_01.* unpert_obs.*
foreach file_rej ( `find . -maxdepth 1 -name "rej_obs_conv_01.*" | xargs ls` )
   cat ${file_rej} >> rej_obs_conv_01
end
rm rej_obs_conv_01.* filtered_obs.0*

# update lateral bdy for coarse domain
if ( ${DOMAIN_ID} == '01' ) then
   cd ${DA_RUNDIR}
#   cp -p ${RC_DATDIR}/${DATE}/wrfbdy_d${DOMAIN_ID} ${DA_RUNDIR}/wrfbdy_d${DOMAIN_ID}
   cp -p ${RUN_BASEDIR}/${PREV_DATE}/wrf/wrfbdy_d${DOMAIN_ID} ${DA_RUNDIR}/wrfbdy_d${DOMAIN_ID}  
   cat >! ${DA_RUNDIR}/parame.in << EOF
&control_param
 da_file            = '${DA_RUNDIR}/wrfvar_output'
 wrf_bdy_file       = '${DA_RUNDIR}/wrfbdy_d${DOMAIN_ID}'
! wrf_input          = '${RC_DATDIR}/${DATE}/wrfinput_d${DOMAIN_ID}'
 wrf_input          = '${RUN_BASEDIR}/${PREV_DATE}/wrf/wrfinput_d${DOMAIN_ID}'
 domain_id          = ${DOMAIN_ID}
 debug              = .false.
 update_lateral_bdy = .true.
 update_low_bdy     = .false
 update_lsm         = .false.
 iswater            = 16 /
EOF
   ln -sf ${DA_SRCDIR}/var/build/da_update_bc.exe .
   time ./da_update_bc.exe >&! update_lat_bc_${DATE}.log
   mv parame.in parame.in.latbdy
   # check status
   grep "Update_bc completed successfully" update_lat_bc_${DATE}.log
   if ( $status != 0 ) then
      echo "ERROR in run_wrfda.csh : update lateral bdy failed..."
      exit 1
   endif
endif

echo "Done run_wrfda.csh for domain ${DOMAIN_ID} `date`"
