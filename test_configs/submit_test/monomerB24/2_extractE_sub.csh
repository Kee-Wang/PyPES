#!/bin/csh -f
               #When using this file, should be used for molpro outputfiles of the same batch. namexxx.out, and should be in same output format. Only the fist colum of $output will be extracted.
               #Need to set the number of output files need to be processed.
               #Need to set the number of atoms of the batch config
               #Need to set the order of atoms, i.e. where is the line of atom coordinates.
               #nfile=number of output files in this folder

               #please change following settings accordingly.

               #----------------------------User define-----------------
               set input = 'monomerB'
               set nfile = 15   #Num of file per node
               #--------------------------------------------------------


               set output = 'monomerB.abE'
               set natom = 5
               set lineAtom1 = 29
               set lineAtom2 = 30
               set lineAtom3 = 31
               set lineAtom4 = 32
               set lineAtom5 = 33
               set lineAtom6 = 34

               set i = 1
               while ( $i <= $nfile )
                 set tail = $i
                 if ( $i < "10" ) then
                   set tail = 00$i
                 else if ($i < "100") then
                   set tail = 0$i
                 endif

               #Make sure the line number of $output and coordinate are corret, using keyword1 and keyword2.
               #if ((head -n $keyword1 $input$tail.out | tail -1) == $linekey1 && (head -n $keyword2 $input$tail.out | tail -1) == $linekey2)


               #Number of atom
               echo $natom >> $output

               #Energy extracted from output files.

               tail -n 3 $input$tail.out |  head -1 | awk '{print $1}' >> $output



               ####What arrangement of atoms do you want?
               head -n $lineAtom1 $input$tail.out | tail -1 >> $output
               head -n $lineAtom2 $input$tail.out | tail -1 >> $output
               head -n $lineAtom3 $input$tail.out | tail -1 >> $output
               #head -n $lineAtom4 $input$tail.out | tail -1 >> $output
               #head -n $lineAtom5 $input$tail.out | tail -1 >> $output
               #head -n $lineAtom6 $input$tail.out | tail -1 >> $output
               #head -n $lineAtom7 $input$tail.out | tail -1 >> $output
               #head -n $lineAtom8 $input$tail.out | tail -1 >> $output
               #head -n $lineAtom9 $input$tail.out | tail -1 >> $output
               #head -n $lineAtom10 $input$tail.out | tail -1 >> $output

               set i = `expr $i + 1`

               #else
               #echo 'There is a bad file in file $i' >> badfilereport
               #endif
               end

                