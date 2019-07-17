#!/bin/csh

#set echo

set scscript = run_wrf_sub.csh  # subsequent cycles

set stime = 201806020000
set etime = 201806020100
set tint = 60   # min for cycle

set tint1 = 360  # min for  forecast

set END_DATE0 = 201806020600
set expnum = $1
if( $expnum == "" ) then
  echo "enter exp id looking like exp001"
  exit
endif
set IC_dir = /glade/scratch/wangs/02jun2018/wrfda_cycling/$expnum

# for wrf run

set BIN_DIR     = /glade/u/home/wangs/test/WRFDA/var/build
set WRF_SRCDIR  = /glade/u/home/wangs/WRFV3_test/
set RUN_BASEDIR = /glade/scratch/wangs/02jun2018/wrfda_cycling/$expnum

# for da dir 

set DA_SRCDIR    = /glade/u/home/wangs/test/WRFDA

# set a fake stime to creat a psudo directory of previous forecast cycle
#set stime = `${BIN_DIR}/da_advance_time.exe $stime -${tint}m -f ccyymmddhhnn`

set done = 0
@ itime = $stime
while( $done == 0 )

    set START_DATE = `${BIN_DIR}/da_advance_time.exe $itime 0 -w`
    set END_DATE =  `${BIN_DIR}/da_advance_time.exe $END_DATE0 0 -w` #$END_DATE0 #`${BIN_DIR}/da_advance_time.exe $itime ${tint1}m -w`

    set next_cyc = `${BIN_DIR}/da_advance_time.exe $itime ${tint}m -w`

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

    set ccyy_c = `echo $next_cyc | cut -c1-4`
    set mm_c   = `echo $next_cyc | cut -c6-7`
    set dd_c   = `echo $next_cyc | cut -c9-10`
    set hh_c   = `echo $next_cyc | cut -c12-13`
    set nn_c   = `echo $next_cyc | cut -c15-16`

    if( $ccyy_s$mm_s$dd_s$hh_s$nn_s == $etime ) then
       set done = 1
    endif

    ################################################################################
    ################################################################################
    
    set runscript = scscript_${expnum}_$ccyy_s$mm_s$dd_s$hh_s$nn_s.csh
    cp $scscript $runscript

    set RUN_BASEDIR1 = $RUN_BASEDIR/$ccyy_s$mm_s$dd_s$hh_s$nn_s/fcst/

#      sed -e "s@@@" -i $runscript
    sed -e "s@AADATE@$ccyy_s$mm_s$dd_s$hh_s$nn_s@" -i $runscript
    sed -e "s@AABIN_DIR@$BIN_DIR@" -i $runscript
    sed -e "s@AAWRF_SRCDIR@$WRF_SRCDIR@" -i $runscript
    sed -e "s@AARUN_BASEDIR@$RUN_BASEDIR1@" -i $runscript
    sed -e "s@AAFCST_RANGE@$tint1@" -i $runscript
    sed -e "s@AAIC_dir@$IC_dir@" -i $runscript
    sed -e "s@AAEND_DATE@$END_DATE@" -i $runscript

    echo "Do forecast from $ccyy_s$mm_s$dd_s$hh_s$nn_s  "

#    echo "check wrf script"
#    exit
    qsub $runscript
    #if( ! -e $runscript.log  ) then
    # set chkwrd = 'Done run_wrf.csh'
    # set jobid = `qsub $runscript`
    # set jobid = `echo $jobid | awk -F "." '{printf "%s ",$1}'`
    # echo "jobid is $jobid"
    # set chkdone = 0
    # while( $chkdone == 0 )
 
    #    set chkfile = wrfda_fcst.o$jobid
    #    if( -e $chkfile ) then
    #      set tem1 = `grep "$chkwrd" $chkfile `
    #      if( "$tem1" != "" ) then      
    #        echo $tem1
    #        set chkdone = 1
    #      endif
    #    endif
    #    sleep 10

    # end
    # echo "Done" >& $runscript.log
    #else
    # echo "$runscript has completed"
    #endif  

    set itime = $ccyy_c$mm_c$dd_c$hh_c$nn_c

    #mv $runscript.log $runscript ${RUN_BASEDIR}/  

end


