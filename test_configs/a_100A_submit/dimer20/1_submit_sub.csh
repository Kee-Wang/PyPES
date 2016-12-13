#!/bin/csh -f
            #This file split a file into nfile files, and creat standard molpro input file for each configuration.


            #-----------------------User define-----------------
            set nfile = 81 #No. configs
            set inputfile = 'dimer' #No suffix
            #---------------------------------------------------

            set n = 8
            set i = 1
            while ( $i <= $nfile )
              set tail = $i
              if ( $i < "10" ) then
                set tail = 00$i
              else if ($i < "100") then
                set tail = 0$i
              endif

              set k = `expr $i - 0`
              set j = `expr $k \* $n`
              set input = $inputfile$tail


            echo '*** co2h2o'        >> $input.temp.head
            echo 'memory,100,m'      >> $input.temp.head
            echo 'basis avtz '        >> $input.temp.head
            echo 'geomtype = xyz'    >> $input.temp.head
            echo 'geometry = {'      >> $input.temp.head


            #change here accroding to different input
            head -n $j $inputfile | tail -$n  >> $input.temp.config
            cat $input.temp.head  $input.temp.config >> $input
            rm $input.temp.head
            rm $input.temp.config

            echo ''                  >> $input
            echo '}'                 >> $input
            echo 'hf'                >> $input
            echo 'ccsd(t)-f12'       >> $input
            echo '---'               >> $input

            /usr/local/bin/molprop_2010_1_Linux_x86_64_i8 -n 4 $input
              set  i = `expr $i + 1`
            end
        