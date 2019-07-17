#!/bin/csh

set echo
set fcscript = run_wrf_from_real_sub.csh   # first cycle

set scscript = run_wrf_from_wrfda_sub.csh  # subsequent cycles

set dascript0 = run_wrfda_sub_mkiii.csh # data assimilation

set stime = 201806020000
set etime = 201806020100

#set stime = 201806020000
#set etime = $stime
set tint = 60   # minute

set coldstart = 0
set analdomain = 01
set IC_dir = /glade/scratch/wangs/02jun2018/noDA00/

set nouterloop = 6
set ninnerloop = 150
set cloud_cv_options = 2
set use_cv_w = '.true.'
set twc_opt = 0
set radar_rf_opt = 2
set use_radar_rv = '.false.'
set use_radar_rf = '.false.'
set use_radar_rqv = '.false.'
set use_radar_rhv = '.true.'
set rf_noice = 0
set preproc = 1
set uvlenscl = 1.0
set tlenscl = 1.0
set qvlenscl = 1.0
set hydrolenscl = 1.0

set expid = $1
if( $expid == ""  ) then
  echo "need to enter expid looking like exp001"
  exit
endif
set dascript = run_wrfda_sub_${expid}.csh
if( ! -e $dascript ) then
  cp $dascript0 $dascript
endif 
# for wrf run

set BIN_DIR     = /glade/u/home/wangs/test/WRFDA/var/build
set RC_DATDIR   = /glade/scratch/wangs/02jun2018/wrfda_cycling/rc/
set WRF_SRCDIR  = /glade/u/home/wangs/WRFV3_test/
set RUN_BASEDIR = /glade/scratch/wangs/02jun2018/wrfda_cycling/$expid

# for da run

set DA_SRCDIR    = /glade/u/home/wangs/test/WRFDA
#set radardir     = /glade/work/wangs/02jun2018/radar/dbzerr1.0_ijskip2_altlim4km/
set radardir     = /glade/work/wangs/02jun2018/radar/dbzerr2.0_ijskip2/
set OB_DATDIR    = /glade/work/wangs/02jun2018/wrfda_cycling/ob
set BE_DATDIR    = /glade/work/wangs/jun20/genbe/test_cv7_cvcld2/
set RC_DATDIR    = /glade/scratch/wangs/02jun2018/wrfda_cycling/rc
set RUN_BASEDIR  = /glade/scratch/wangs/02jun2018/wrfda_cycling/$expid

set da_type = ( only_rv only_rf  )
set da_type = ( only_rf  )
set ndatype = $#da_type

# set a fake stime to creat a psudo directory of previous forecast cycle
set stime = `${BIN_DIR}/da_advance_time.exe $stime -${tint}m -f ccyymmddhhnn`

set done = 0
@ itime = $stime
while( $done == 0 )

    if( $coldstart  == 1 ) then
      echo "Cold start"
      set START_DATE = `${BIN_DIR}/da_advance_time.exe $stime 0 -w`
      set END_DATE = `${BIN_DIR}/da_advance_time.exe $stime $tint -w`
      set ccyy_s = `echo $START_DATE | cut -c1-4`
      set mm_s   = `echo $START_DATE | cut -c6-7`
      set dd_s   = `echo $START_DATE | cut -c9-10`
      set hh_s   = `echo $START_DATE | cut -c12-13`
      set ccyy_e = `echo $END_DATE | cut -c1-4`
      set mm_e   = `echo $END_DATE | cut -c6-7`
      set dd_e   = `echo $END_DATE | cut -c9-10`
      set hh_e   = `echo $END_DATE | cut -c12-13`
      set runscript = fcscript_$ccyy_s$mm_s$dd_s$hh_s.csh 
      cp $fcscript $runscript

#      sed -e "s@@@" -i $runscript
      sed -e "s@AADATE@$ccyy_s$mm_s$dd_s$hh_s@" -i $runscript
      sed -e "s@AABIN_DIR@$BIN_DIR@" -i $runscript 
      sed -e "s@AARC_DATDIR@$RC_DATDIR@" -i $runscript
      sed -e "s@AAWRF_SRCDIR@$WRF_SRCDIR@" -i $runscript
      sed -e "s@AARUN_BASEDIR@$RUN_BASEDIR@" -i $runscript

      set chkwrd = 'Done run_wrf_from_real.csh'
      set jobid = `qsub $runscript`
      set jobid = `echo $jobid | awk -F "." '{printf "%s ",$1}'`
      echo "jobid is $jobid"
      set chkdone = 0
      while( $chkdone == 0 )

         set chkfile = wrf.o$jobid
         if( -e $chkfile ) then
           set tem1 = `grep "$chkwrd" $chkfile `
           if( "$tem1" != "" ) then
             echo $tem1
             set chkdone = 1
           endif
         endif
         sleep 10

      end      

    endif

    set START_DATE = `${BIN_DIR}/da_advance_time.exe $itime 0 -w`
    set END_DATE = `${BIN_DIR}/da_advance_time.exe $itime ${tint}m -w`
    set ccyy_s = `echo $START_DATE | cut -c1-4`
    set mm_s   = `echo $START_DATE | cut -c6-7`
    set dd_s   = `echo $START_DATE | cut -c9-10`
    set hh_s   = `echo $START_DATE | cut -c12-13`
    set nn_s   = `echo $START_DATE | cut -c15-16`

    set ccyy_e = `echo $END_DATE | cut -c1-4`
    set mm_e   = `echo $END_DATE | cut -c6-7`
    set dd_e   = `echo $END_DATE | cut -c9-10`
    set hh_e   = `echo $END_DATE | cut -c12-13` 
    set nn_e   = `echo $END_DATE | cut -c15-16`  

    echo "Do DA at $ccyy_e$mm_e$dd_e$hh_e$nn_e"

    if( $itime == $stime ) then
      set prevdir = $RUN_BASEDIR/$ccyy_s$mm_s$dd_s$hh_s$nn_s/wrf/
      mkdir -p $prevdir
      cp $IC_dir/wrfout_d${analdomain}_${ccyy_e}-${mm_e}-${dd_e}_${hh_e}:${nn_e}:00 \
         $prevdir
      cp $IC_dir/wrfbdy_d${analdomain} $prevdir
      cp $IC_dir/wrfinput_d${analdomain} $prevdir
    endif

    @ idastep = 1
    while( $idastep <= $ndatype )

      set runscript = dascript_${expid}_$ccyy_e$mm_e$dd_e$hh_e${nn_e}_$idastep.csh
      cp $dascript $runscript

      if( $idastep == 1 ) then
        sed -e "s@AApreproc@$preproc@" -i $runscript        
      else
        sed -e "s@AApreproc@0@" -i $runscript
      endif

#     sed -e "s@@@" -i $runscript
      sed -e "s@AADATE@$ccyy_e$mm_e$dd_e$hh_e$nn_e@" -i $runscript
      sed -e "s@AADA_SRCDIR@$DA_SRCDIR@" -i $runscript
      sed -e "s@AAOB_DATDIR@$OB_DATDIR@" -i $runscript
      sed -e "s@AABE_DATDIR@$BE_DATDIR@" -i $runscript
      sed -e "s@AARC_DATDIR@$RC_DATDIR@" -i $runscript
      sed -e "s@AARUN_BASEDIR@$RUN_BASEDIR@" -i $runscript
      sed -e "s@AAradardir@$radardir@" -i $runscript
      sed -e "s@AACYCLE_PERIOD@$tint@" -i $runscript
      sed -e "s@AAidastep@$idastep@g" -i $runscript
      #sed -e "s@AAnouterloop@$nouterloop@" -i $runscript
      sed -e "s@AAninnerloop@$ninnerloop@" -i $runscript
      sed -e "s@AAcloud_cv_options@$cloud_cv_options@" -i $runscript
      sed -e "s@AAuse_cv_w@$use_cv_w@" -i $runscript
      sed -e "s@AAtwc_opt@$twc_opt@" -i $runscript
      sed -e "s@AAradar_rf_opt@$radar_rf_opt@" -i $runscript
      #sed -e "s@AAuse_radar_rv@$use_radar_rv@" -i $runscript
      #sed -e "s@AAuse_radar_rf@$use_radar_rf@" -i $runscript
      #sed -e "s@AAuse_radar_rqv@$use_radar_rqv@" -i $runscript
      #sed -e "s@AAuse_radar_rhv@$use_radar_rhv@" -i $runscript
      sed -e "s@AArf_noice@$rf_noice@" -i $runscript
      sed -e "s@AAuvlenscl@$uvlenscl@" -i $runscript
      sed -e "s@AAtlenscl@$tlenscl@" -i $runscript
      sed -e "s@AAqvlenscl@$qvlenscl@" -i $runscript
      sed -e "s@AAhydrolenscl@$hydrolenscl@" -i $runscript

      if( $da_type[$idastep] == "only_rv" ) then
        sed -e "s@AAuse_radar_rv@.true.@" -i $runscript
        sed -e "s@AAuse_radar_rf@.false.@" -i $runscript
        sed -e "s@AAuse_radar_rqv@.false.@" -i $runscript
        sed -e "s@AAuse_radar_rhv@.false.@" -i $runscript
        sed -e "s@AAnouterloop@$nouterloop@" -i $runscript
      endif

      if( $da_type[$idastep] == "only_rf" ) then
        sed -e "s@AAnouterloop@$nouterloop@" -i $runscript
        sed -e "s@AAuse_radar_rv@.false.@" -i $runscript
        sed -e "s@AAuse_radar_rf@$use_radar_rf@" -i $runscript
        sed -e "s@AAuse_radar_rqv@$use_radar_rqv@" -i $runscript
        sed -e "s@AAuse_radar_rhv@$use_radar_rhv@" -i $runscript
      endif
 
      if( $da_type[$idastep] == "both_rvrf" ) then
        sed -e "s@AAnouterloop@$nouterloop@" -i $runscript
        sed -e "s@AAuse_radar_rv@$use_radar_rv@" -i $runscript
        sed -e "s@AAuse_radar_rf@$use_radar_rf@" -i $runscript
        sed -e "s@AAuse_radar_rqv@$use_radar_rqv@" -i $runscript
        sed -e "s@AAuse_radar_rhv@$use_radar_rhv@" -i $runscript
      endif

      set chkwrd = "Done run_wrfda.csh for domain"

      if( ! -e $runscript.log ) then
      

       set jobid = `qsub $runscript`
       set jobid = `echo $jobid | awk -F "." '{printf "%s ",$1}'`
       echo "jobid is $jobid"
       set chkdone = 0
       while( $chkdone == 0 )

         set chkfile = wrfda.o$jobid
         if( -e $chkfile ) then
           set tem1 = `grep "$chkwrd" $chkfile `
           if( "$tem1" != "" ) then
             echo $tem1
             set chkdone = 1
           endif
         endif
         sleep 10
       end
       echo "Done" >& $runscript.log
      else
       echo "$runscript has completed"
      endif

      @ idastep++
    end

    if( $ccyy_e$mm_e$dd_e$hh_e$nn_e == $etime ) then
       set done = 1
       continue
    endif
#    echo "check da result"
#    exit
    ################################################################################
    ################################################################################
  
    set runscript = scscript_${expid}_$ccyy_e$mm_e$dd_e$hh_e${nn_e}.csh
    cp $scscript $runscript

#      sed -e "s@@@" -i $runscript
    sed -e "s@AADATE@$ccyy_e$mm_e$dd_e$hh_e$nn_e@" -i $runscript
    sed -e "s@AABIN_DIR@$BIN_DIR@" -i $runscript
    sed -e "s@AARC_DATDIR@$RC_DATDIR@" -i $runscript
    sed -e "s@AAWRF_SRCDIR@$WRF_SRCDIR@" -i $runscript
    sed -e "s@AARUN_BASEDIR@$RUN_BASEDIR@" -i $runscript
    sed -e "s@AAFCST_RANGE@$tint@" -i $runscript

    echo "Do forecast from $ccyy_e$mm_e$dd_e$hh_e$nn_e  "

#    echo "check wrf script"
#    exit
    if( ! -e $runscript.log  ) then
      set chkwrd = 'Done run_wrf_from_wrfda.csh'
      set jobid = `qsub $runscript`
      set jobid = `echo $jobid | awk -F "." '{printf "%s ",$1}'`
      echo "jobid is $jobid"
      set chkdone = 0
      while( $chkdone == 0 )
   
        set chkfile = wrf.o$jobid
        if( -e $chkfile ) then
          set tem1 = `grep "$chkwrd" $chkfile `
          if( "$tem1" != "" ) then      
            echo $tem1
            set chkdone = 1
          endif
        endif
        sleep 10

      end
      echo "Done" >& $runscript.log
    else
      echo "$runscript has completed"
    endif  

    set itime = $ccyy_e$mm_e$dd_e$hh_e$nn_e
#    
#    echo "check wrf result"
#    exit
   
#    cp $fcscript fcscript.csh
#    cp $scscript scscript.csh
#    cp $dascript dascript.csh

end


mv *${expid}*.log da*${expid}*.csh sc*${expid}*.csh $RUN_BASEDIR  #${expid}/

