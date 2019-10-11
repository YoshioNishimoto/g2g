C       VERSION HISTROTY
C 2015/02/16
C     - Bug in calculating Raman intensity is fixed
C 2015/02/13
C     - Fix half-width output in plot files
C 2015/02/03
C     - Bug in calculating ramanint is fixed
C 2015/01/07
C     - Add -go and -o options, to prevent from overwriting test.out
C       mistakenly
C     - Increased integer digits for some outputs
C     - Compatible with explicit Hessian contribution of atoms
C       but this may not be perfect
C 2015/06/01
C     - Some reading sections have been changed, because the format
C       of GAMESS output may be changed slightly.
C     - FMO Hessian/Raman
C 2016/06/03
C     - Resonance Raman spectrum with the short-time approximation
C     - Add -scal option.
C
      program test

      implicit double precision (a-h,o-z)

      character(80) file1,file2,line,filego
      character( 8) cdum
      character( 1) jc
      character(23) cfmt
      double precision, allocatable :: dipder(:,:),vec(:,:),freq(:),
     * reduced(:),irint(:),ramanact(:,:),depolar(:),zan(:),coord(:,:),
     * zmass(:),alphader(:,:),depu(:),ramanint(:,:),ramintout(:,:),
     * ramactout(:,:),irout(:),grad(:,:)
      integer, allocatable :: ind(:)
      double precision alpha(3,3)
      integer rottra(6)
      logical raman,normalize,alpder,aaa,high,excimg,plot,sq,loren,
     * excsus,skip,overwrite,fmo,rrs
C
C     ASSIGN INITIAL PARAMETERS
C
      hw = 1.0d+06
      temperature = 1.0d+06
      excimg = .true.
      iram = 1
      plot = .false.
      sq = .true.
      f0 = 19430.0d+00
      loren = .true.
      excsus = .false.
      overwrite = .false.
      filego = "test.out"
      fmo = .false.
      rrs = .false.
      scal = 1.0d+00
C
C     CHECK OPTIONS
C
      i = 1
      do
        cdum(1:8) = '        '
        call getarg(i,cdum)
        if (i.eq.1.and.cdum(1:1).ne.' '.and.cdum(1:3).ne.'-pl'
     *     .and.cdum(1:3).ne.'-h '.and.cdum(1:3).ne.'-H ') then
          write (*,*) "------- OPTIONS ------"
        else if (cdum(1:1).eq.' ') then
          exit
        end if
        i = i + 1
        if (cdum(1:3).eq."-hw") then
          call getarg(i,cdum)
          i = i + 1
          if (cdum(1:1).eq.' '.or.cdum(1:1).eq.'-') then
            write (*,*) "STRANGE ARGUMENT FOR -hw OPTION"
            stop
          end if
          read (cdum,*) hw
          write (*,'("  HALF-WIDTH  = ",f8.1," cm-1")') hw
        else if (cdum(1:2).eq."-t") then
          call getarg(i,cdum)
          i = i + 1
          if (cdum(1:1).eq.' '.or.cdum(1:1).eq.'-') then
            write (*,*) "STRANGE ARGUMENT FOR -t OPTION"
            stop
          end if
          read (cdum,*) temperature
          write (*,'("  TEMPERATURE = ",f8.2," K")') temperature
C       else if (cdum(1:2).eq."-l") then
C         call getarg(i,cdum)
C         i = i + 1
C         if (cdum(1:1).eq.' '.or.cdum(1:1).eq.'-') then
C           write (*,*) "STRANGE ARGUMENT FOR -l OPTION"
C           stop
C         end if
C         read (cdum,*) f0
C         write (*,'("  WAVE LENGTH = ",f8.2," nm")') f0
        else if (cdum(1:2).eq."-w") then
          call getarg(i,cdum)
          i = i + 1
          if (cdum(1:1).eq.' '.or.cdum(1:1).eq.'-') then
            write (*,*) "STRANGE ARGUMENT FOR -w OPTION"
            stop
          end if
          read (cdum,*) f0
          write (*,'("  WAVE NUMBER = ",f8.2," cm-1")') f0
        else if (cdum(1:3).eq."-go") then
          call getarg(i,cdum)
          i = i + 1
          if (cdum(1:1).eq.' ') then
            write (*,*) "STRANGE ARGUMENT FOR -go OPTION"
            stop
          end if
          filego(1:8) = '        '
          read (cdum,*) filego
          lf = len_trim(filego)
          overwrite = .true.
          write (*,'("  NAME OF GAUSSIAN-TYPE OUTPUT: ",A)')filego(1:lf)
        else if (cdum(1:5).eq."-scal") then
          call getarg(i,cdum)
          i = i + 1
          if (cdum(1:1).eq.' '.or.cdum(1:1).eq.'-') then
            write (*,*) "STRANGE ARGUMENT FOR -scal OPTION"
            stop
          end if
          read (cdum,*) scal
          write (*,'("  SCALING OF FREQUENCIES = ",f7.4)') scal
          if (scal.le.0.0d+00) then
            write (*,*) "STRANGE SCALING FACTOR?"
          end if
        else if (cdum(1:3).eq."-pl") then
          plot = .true.
        else if (cdum(1:2).eq."-h") then
          go to 500
        else if (cdum(1:2).eq."-H") then
          go to 400
        else if (cdum(1:2).eq."-a") then
          sq = .false.
        else if (cdum(1:2).eq."-e") then
          excsus = .true.
          write (*,*) " EXCLUDE SUSPICIOUS FROM RAMAN INTENSITY"
        else if (cdum(1:2).eq."-s") then
          sq = .true.
        else if (cdum(1:2).eq."-n") then
          iram = 1
          write (*,*) " USE NATURAL INCIDENT LIGHT"
        else if (cdum(1:2).eq."-p") then
          iram = 2
          write (*,*) " USE PLANE-POLARIZED LIGHT"
        else if (cdum(1:2).eq."-i") then
          excimg = .false.
          write (*,*) " INCLUDE IMAGINARY FREUQNEICES FOR PLOT"
        else if (cdum(1:2).eq."-l") then
          loren = .true.
          write (*,*) " USE LORENTZIAN FUNCTION"
        else if (cdum(1:2).eq."-g") then
          loren = .false.
          write (*,*) " USE GAUSSIAN FUNCTION"
        else if (cdum(1:2).eq."-o") then
          overwrite = .true.
          write (*,*) " OVERWRITE test.out ANYWAY"
        else if (cdum(1:4).eq."-fmo") then
          fmo = .true.
          write (*,*) " FMO OUTPUT"
        else if (cdum(1:4).eq."-rrs") then
          rrs = .true.
          write (*,*) " RESONANCE RAMAN SPECTRUM"
        else if (cdum(1:1).eq.' ') then
          exit
        else
          write (*,*) "STRANGE ARGUMENT NO",i-1
          write (*,'(" ARGUMENT = ", A8)') cdum
          stop
        end if
      end do
      if (plot) then
        call plot_dat(hw,excimg,sq,loren)
      end if
      if (i.ne.1) then
        write (*,*) "--- END OF OPTIONS ---"
      end if
C
C     START
C
      write (*,*) 
      write (*,*) " --------------------------------------"
      write (*,*) "  A Post-Processing Program for GAMESS"
      write (*,*) "       LAST UPDATED :: 2016/06/03"
      write (*,*) " --------------------------------------"
      write (*,*) 
      write (*,*) "NAME OF GAMESS OUTPUT FILE?"
      read (*,'(a80)') file1
      if (file1(1:2).eq."-h") go to 500
      write (*,*) 

      open (10,file=file1,status="old",iostat=is)
      if (is.ne.0) then
        write (*,*)
        write (*,*) "FAILED IN OPENING THE OUTPUT FILE"
        write (*,*)
        stop
      end if
      iatom = 0
      icoord = 0
      nrt = 6
      nat = 0
      do
        read(10,'(A)',end=100) line
        if (line(29:33).eq."RAMAN") then
          write (*,*)
          write (*,*) "  OUTPUT FILE INDICATES THAT"
          write (*,*) "  THIS IS A RAMAN CALCULATION"
          write (*,*)
          raman = .true.
          if (.not.fmo) go to 100
        else if (line(29:33).eq."HESSI") then
          write (*,*)
          write (*,*) "  OUTPUT FILE INDICATES THAT"
          write (*,*) "  THIS IS AN IR CALCULATION"
          write (*,*)
          raman = .false.
          if (.not.fmo) go to 100
        else if (line(19:44).eq."RECOGNIZED AS BEING LINEAR") then
          nrt = 5
          write (*,*)
          write (*,*) "INPUT MOLECULE SEEMS TO BE LINEAR"
          write (*,*) "ASSUME NUMBER OF ROT. AND TRANS. MODES IS 5"
        end if
        if (fmo) then
          !! valid only for DFTB?
          if (line(2:33).eq."Total number of basis functions:") then
            read (line(37:45),*) num
            write (*,*) "NUMBER OF ORBITALS = ", num
            l1 = num
            go to 100
          end if
          if (line(2:23).eq."Total number of atoms:") then
            read (line(37:45),*) nat
            write (*,*) "NUMBER OF ATOMS =    ", nat
C           write (*,*) "ALLOCATE SOME ARRAYS"
            ncomp = nat*3
            allocate (dipder(3,ncomp))
            allocate (vec(ncomp,ncomp))
C           write (*,*) "ALLOCATION FINISHED"
          end if
        else
          if (line(2:20).eq."NUMBER OF CARTESIAN") then
            read (line(50:52),*) num
            write (*,*) "NUMBER OF ORBITALS = ", num
            l1 = num
          end if
          if (line(2:22).eq."TOTAL NUMBER OF ATOMS") then
            read (line(49:52),*) nat
            write (*,*) "NUMBER OF ATOMS =    ", nat
C           write (*,*) "ALLOCATE SOME ARRAYS"
            ncomp = nat*3
            allocate (dipder(3,ncomp))
            allocate (vec(ncomp,ncomp))
C           write (*,*) "ALLOCATION FINISHED"
          end if
        end if
      end do

100   if (raman) then
        if (rrs) then
          write (*,*) "YOU CANNOT EVALUATE NON-RESONANCE AND RESONANCE "
     *               ,"RAMAN SIMULTANEOUSLY WITH THIS CODE"
          stop
        end if
        alpder = .true.
        write (*,*) "ALPHA POLARIZABILITY CAN BE OBTAINED FROM DAT FILE"
        write (*,*) "(DAT FILE OF GAMESS, WHICH CORRESPONDS TO .F07)"
        write (*,*) "DO YOU WANT TO DO SO? [Y/N]"
        read (*,'(a80)') line
        if (line(1:1).eq."y".or.line(1:1).eq."Y") then
          write (*,*)
          write (*,*) "GIVE THE NAME OF YOUR DAT FILE (CANCEL BY ""N"")"
          line(1:2) = "  "
          read (*,'(a80)') line
C         write (*,*) "line(1:2) = ", line(1:2)
          if ((line(1:1).eq."n".or.line(1:1).eq."N")
     *      .and.line(2:2).eq." ") then
            alpder = .false.
C           write (*,*) "do not read alpha derivative tensor"
            write (*,*)
            write (*,*) "ALPHA POLARIZABILITY WILL BE TAKEN FROM OUTPUT"
            write (*,*)
            io11 = 0
          else
            length = len_trim(line)
C           write (*,*) "trying to open: ", line(1:length)
            open (11,file=line,status="old",iostat=is)
            if (is.eq.0) then
              write (*,*)
              write (*,*) "THE DAT FILE IS SUCCESSFULLY OPENED"
              write (*,*)
            else
              write (*,*)
              write (*,*) "FAILED IN OPENING THE DAT FILE"
              write (*,*)
              stop
            end if
            alpder = .false.
            io11 = 1
          end if
        else
          write (*,*)
          write (*,*) "ALPHA POLARIZABILITY WILL BE TAKEN FROM OUTPUT"
          write (*,*)
        end if
      end if
      nrt0 = nrt
      if (nat.eq.0) write (*,*) "NUMBER OF ATOMS IS MISSING"
      if (num.eq.0) write (*,*) "NUMBER OF ORBITALS IS MISSING"
      if (nat.eq.0.or.num.eq.0) stop

      rewind (10)
C     write (*,*) "FINISHED READING VARIABLES"
      allocate (zan(nat))
      allocate (coord(3,nat))
      allocate (freq(ncomp))
      allocate (reduced(ncomp))
      allocate (irint(ncomp))
      allocate (zmass(nat))
      allocate (ind(3*nat))
      allocate (irout(5000))
      if (raman) then
C       write (*,*) "raman keyword was found"
        allocate (ramanact(2,ncomp))
        allocate (depolar(ncomp))
        allocate (alphader(6,ncomp))
        allocate (depu(ncomp))
        allocate (ramanint(2,ncomp))
        allocate (ramintout(2,5000))
        allocate (ramactout(2,5000))
C     else
C       write (*,*) "this is not RUNTYP=RAMAN"
      end if
      if (rrs) allocate (grad(3,nat))
C     write (*,*) "ALLOCATION FINISHED"
      if (fmo) coord = 0.0d+00

      nn = floor(dble(ncomp)/5.0d+00)
      mm = mod(ncomp,5)
C     write (jc,'(i1)') mm
C     write (*,*) "nn and mm = ",nn,mm
      do
        read(10,'(A)',end=300) line
        if (line(2:17).eq."ATOM      ATOMIC") then !! non-FMO
          read (10,*)
          do i = 1, nat
            read (10,*) cdum,zan(i),coord(1,i),coord(2,i),coord(3,i)
          end do
        end if
        if (line(11:34).eq."DIPOLE DERIVATIVE TENSOR") then
          write (*,*) "READING DIPOLE DERIVATIVE ..."
          read (10,*)
          read (10,*)
          read (10,*)
          do i = 1, ncomp
            read (10,'(21x,3(x,F14.9))') (dipder(j,i),j=1,3)
          end do
C         write (*,*) "finished."
        end if
        if (line(11:30).eq."ATOMIC WEIGHTS (AMU)") then
          write (*,*) "READING ATOMIC WEIGHTS ..."
          read (10,*)
          do i = 1, nat
            read (10,*),idum,cdum,zmass(i)
          end do
        end if
C       if(line(15:53).eq."ARE TAKEN AS ROTATIONS AND TRANSLATIONS")then
        if(line(21:59).eq."ARE TAKEN AS ROTATIONS AND TRANSLATIONS")then
          if (line(7:8).eq."**".or.line(12:13).eq."**") then
            write (*,*) "*** CANNOT GET THE INDICES OF ROTATION AND ",
     *                  "TRANSLATION ***"
            nstart = -1
            nlast = -1
          else
C           read (line(1:13),'(6x,i2,3x,i2)') nstart,nlast
            read (line(1:19),'(6x,i5,3x,i5)') nstart,nlast
            if (nstart.ne.1) then
              write (*,*)
              write (*,'(" *************************************")')
              write (*,'(" THERE MAY BE ",i2," IMAGINARY FREQUENCIES")')
     *        nstart-1
              write (*,'(" *************************************")')
              write (*,*)
              do i = 1, nstart-1
                ind(i) = -i
              end do
              do i = 1, nrt
C               rottra(nrt) = nstart + nrt-1
                rottra(nrt) = nstart + i-1
              end do
C             rottra(1) = nstart
C             rottra(2) = nstart+1
C             rottra(3) = nstart+2
C             rottra(4) = nstart+3
C             rottra(5) = nstart+4
C             rottra(6) = nstart+5
              do i = nlast+1, nat*3
                ind(i-nrt) = i
              end do
            else
              do i = 1, nat*3
                ind(i) = i + nrt
              end do
              do i = 1, nrt
                rottra(i) = i
              end do
            end if
          end if
        end if
        if (line(11:37).eq."ALPHA POLARIZABILITY TENSOR") then
          if (alpder) then
            read (10,*)
            read (10,*)
            read (10,*)
            read (10,'(12x, f15.5)') alpha(1,1)
            read (10,'(12x,2f15.5)') alpha(2,1),alpha(2,2)
            read (10,'(12x,3f15.5)') alpha(3,1),alpha(3,2),alpha(3,3)
            alpha(1,2) = alpha(2,1)
            alpha(1,3) = alpha(3,1)
            alpha(2,3) = alpha(3,2)
          end if
        end if
        if (line(11:41).eq."ALPHA POLARIZABILITY DERIVATIVE") then
          if (alpder) then
            read (10,*)
            read (10,*)
            read (10,*)
            do i = 1, ncomp
              read (10,'(17x,6f10.5)') (alphader(j,i),j=1,6)
            end do
          end if
        end if

        if ((raman.and.(line(6:22).eq."RAMAN INTENSITIES")).or.
     * ((.not.raman).and.(line(6:27).eq."REDUCED MASSES IN AMU."))) then
          write (*,*) "READING VIBRATIONAL ANALYSIS ..."
          do i = 1, nn
            read(10,*)      !! blank line
            read(10,'(A80)') line !! mode number
            if (line(6:60).eq.
     *    '*******************************************************')then
              write (*,*) "SKIP WARNING LINES..."
              read (10,*)
              read (10,*)
              read (10,*)
              read (10,*)
              read (10,*)
            end if
            do
              read(10,'(A)',end=300) line
              if (line(3:17).eq."     FREQUENCY:") then
                read (line(1:80),'(1x,16x,3x,5(f10.2,2x))')
     *                (freq(k),k=5*(i-1)+1,5*(i-1)+5)
C             else if (line(3:17).eq."      SYMMETRY:") then
C               read (line,'(1x,16x,3x,5(f20.2,2x))') line
              else if (line(3:17).eq."  REDUCED MASS:") then
                read (line(1:80),'(1x,16x,3x,5(f10.5,2x))')
     *                (reduced(k),k=5*(i-1)+1,5*(i-1)+5)
              else if (line(3:17).eq."  IR INTENSITY:") then
                read (line(1:80),'(1x,16x,3x,5(f10.5,2x))')
     *                (irint(k),k=5*(i-1)+1,5*(i-1)+5)
              else if (line(3:17).eq."RAMAN ACTIVITY:".and.alpder) then
C               read (line(1:80),'(1x,16x,3x,5(f12.3))')
                read (line(1:80),'(1x,16x,3x,5(f10.3,2x))')
     *                (ramanact(2,k),k=5*(i-1)+1,5*(i-1)+5)
              else if (line(3:17).eq."DEPOLARIZATION:") then
                read (line(1:80),'(1x,16x,3x,5(f10.3,2x))')
     *                (depolar(k),k=5*(i-1)+1,5*(i-1)+5)
              else if (line(3:17).eq."               ") then
                do j = 1, ncomp
                  read (10,'(20x,5f12.8)')
     *              (vec(j,k),k=5*(i-1)+1,5*(i-1)+5)
                end do
                read (10,*) !! blank
                read (10,*) !! trans. sayvetz x
                read (10,*) !! y
                read (10,*) !! z
                read (10,*) !! total
                read (10,*) !! blank
                read (10,*) !! rot. sayvetz x
                read (10,*) !! y
                read (10,*) !! z
                read (10,*) !! total
                exit
              end if
            end do
          end do
          if (mm.ne.0) then
            read (10,*) !! blank
            read (10,*) !! mode number
C           write (*,*) "read last line",mm
            do
              read(10,'(A)',end=300) line
              if (line(3:17).eq."     FREQUENCY:") then
                write (cfmt,'("(1x,16x,3x,",i1,"(f10.2,2x))")') mm
                read (line(1:80),fmt=cfmt)
     *                (freq(k),k=5*(i-1)+1,5*(i-1)+mm)
C             else if (line(3:17).eq."      SYMMETRY:") then
C               read (line,'(1x,16x,3x,5(f20.2,2x))') line
              else if (line(3:17).eq."  REDUCED MASS:") then
                write (cfmt,'("(1x,16x,3x,",i1,"(f10.5,2x))")') mm
                read (line(1:80),cfmt)
     *                (reduced(k),k=5*(i-1)+1,5*(i-1)+mm)
              else if (line(3:17).eq."  IR INTENSITY:") then
                write (cfmt,'("(1x,16x,3x,",i1,"(f10.5,2x))")') mm
                read (line(1:80),cfmt)
     *                (irint(k),k=5*(i-1)+1,5*(i-1)+mm)
              else if (line(3:17).eq."RAMAN ACTIVITY:".and.alpder) then
                write (cfmt,'("(1x,16x,3x,",i1,"(f12.3))")') mm
                read (line(1:80),cfmt)
     *                (ramanact(2,k),k=5*(i-1)+1,5*(i-1)+mm)
              else if (line(3:17).eq."DEPOLARIZATION:") then
                write (cfmt,'("(1x,16x,3x,",i1,"(f10.3,2x))")') mm
                read (line(1:80),cfmt)
     *                (depolar(k),k=5*(i-1)+1,5*(i-1)+mm)
              else if (line(3:17).eq."               ") then
                do j = 1, ncomp
                  if (mm.eq.1) then
                    read (10,'(20x,1f12.8)') vec(j,5*(i-1)+1)
                  else if (mm.eq.2) then
                    read (10,'(20x,2f12.8)')
     *                (vec(j,k),k=5*(i-1)+1,5*(i-1)+2)
                  else if (mm.eq.3) then
                    read (10,'(20x,3f12.8)')
     *                (vec(j,k),k=5*(i-1)+1,5*(i-1)+3)
                  else if (mm.eq.4) then
                    read (10,'(20x,4f12.8)')
     *                (vec(j,k),k=5*(i-1)+1,5*(i-1)+4)
                  end if
                end do
                exit
              end if
            end do
          end if
          do
            read(10,'(A)',end=300) line
            if (line(3:14).eq."MODE FREQ(CM") then
              if (raman.and.alpder) then
                do i = 1, ncomp
C                 read (10,*) idum,freq(i),cdum,reduced(i),irint(i),
C    *                        ramanact(2,i),depolar(i)
                  read (10,'(1x,i5,f12.3,4x,a4,1x,4f12.6)')
     *              idum,freq(i),cdum,reduced(i),irint(i),
     *              ramanact(2,i),depolar(i)
                end do
              else
                do i = 1, ncomp
C                 read (10,*) idum,freq(i),cdum,reduced(i),irint(i)
                  read (10,'(1x,i5,f12.3,4x,a4,1x,2f12.6)')
     *              idum,freq(i),cdum,reduced(i),irint(i)
                end do
              end if
              go to 200
            end if
          end do
        end if
      end do

200   continue
C     write (*,*) "dipder and vectors are read"
      close (10)

      if (.not.alpder.and.io11.eq.1) then
        write (*,*)
        write (*,*) "START READING .DAT FILE"
        do
          read (11,'(A)',end=300) line
          if (line(1:7).eq.' $ALPDR') then
C           write (*,*) "$ALPDR keyword found"
            write (*,*) "READING ALPHA DERIVATIVE ..."
            read (11,*) ! E
            read (11,*) ! alpha polarizability line
            read (11,'(6e13.6)') alpha(1,1),alpha(2,1),alpha(2,2),
     *        alpha(3,1),alpha(3,2),alpha(3,3)
            alpha(1,2) = alpha(2,1)
            alpha(1,3) = alpha(3,1)
            alpha(2,3) = alpha(3,2)
            alpha = alpha*(0.5291772108**3)
            read (11,*) ! alpha polarizability nuclear ...
            do i = 1, ncomp
              read (11,'(1p,6e13.6)') (alphader(j,i),j=1,6)
            end do
            exit
          end if
        end do
        close (11)
      end if
      write (*,*)
      write (*,*) "-------------------------------------------"
      write (*,*) "ALL VARIABLES AND PARAMETERS HAVE BEEN READ"
      write (*,*) "NOW, START POST-PROCESSING BY THIS PROGRAM "
      write (*,*) "-------------------------------------------"
      write (*,*)
      nzero = 0
      do i = 1, ncomp
        if (freq(i).le.1.0d+00.and.freq(i).ge.0.0d+00) nzero = nzero + 1
      end do
      if (nzero.gt.nrt) then
        write (*,'(" THERE ARE ",i3," NEAR-ZERO FREQUENCIES")') nzero
        write (*,'(" MAYBE -e OPTION IS USEFUL IN THIS CASE")')
        write (*,*)
        write (*,'(" BECAUSE OF PROJECTION OF EXPLICIT FREQUENCIES?")')
        read (*,'(a80)') line
        if (line(1:1).eq."y".or.line(1:1).eq."Y") then
          if (mod(nzero,3).ne.0) then
            write (*,*)
            write (*,*) "*** WARNING ***"
            write (*,*) "MAYBE THIS IS A LINEAR MOLECULE"
            write (*,*) "THIS MAY NOT WORK WELL FOR THAT CASE"
            write (*,*)
          end if
          fqprev = 1.0d+06
          nimg = 0.0d+00
          do i = 1, ncomp
            fqcurr = freq(i)
            if (fqprev.lt.fqcurr) then
              nimg = i
              exit
            end if
          end do
          natprj = (nzero-nrt)/3
          nstart = nzero-nrt+1
          nlast  = nzero
          write (*,*) "FOLLOWING VALUES ARE ONLY GUESS INFORMATION"
          write (*,'("   NUMBER OF IMAGINARY FREQUENCIES     :",i4)')
     *      nimg
          write (*,'("   NUMBER OF PROJECTED ATOMS           :",i4)')
     *      natprj
          write (*,'("   FIRST INDEX OF ROTATION/TRANSLATION :",i4)')
     *      nstart
          write (*,'("   LAST  INDEX OF ROTATION/TRANSLATION :",i4)')
     *      nlast
          write (*,*)
          do i = 1, nrt
            rottra(nrt) = nstart + nrt-1
          end do
          nrt = nrt + natprj*3
          do i = nlast+1, nat*3
            ind(i-nrt) = i
          end do
        else if (nstart.eq.-1.or.nlast.eq.-1) then
          fqprev = 1.0d+06
          nimg = 0.0d+00
          do i = 1, ncomp
            fqcurr = freq(i)
            if (fqprev.lt.fqcurr) then
              nimg = i
              exit
            end if
          end do
          nstart = nimg + 1
          nlast = nstart + nzero - 1
          write (*,*) "FOLLOWING VALUES ARE ONLY GUESS INFORMATION"
          write (*,'("   NUMBER OF IMAGINARY FREQUENCIES     :",i4)')
     *      nimg
          write (*,'("   FIRST INDEX OF ROTATION/TRANSLATION :",i4)')
     *      nstart
          write (*,'("   LAST  INDEX OF ROTATION/TRANSLATION :",i4)')
     *      nlast
          write (*,*)
          do i = 1, nimg
            ind(i) = -i
          end do
          do i = 1, nrt
C           rottra(nrt) = nstart + nrt-1
            rottra(nrt) = nlast - nrt + i
          end do
          do i = nlast+1, nat*3
C           ind(i-nrt) = i
            ind(i-nlast) = i
          end do
          nrt = nlast
          do i = 1, ncomp
            write (11,*) i,ind(i)
          end do
        end if
      end if
      if (nstart.eq.-1.or.nlast.eq.-1) then
        write (*,*) "somehow nstart or nlast is still negative..."
        stop
      end if
C
      if (rrs) then
        write (*,'(" RESONANCE RAMAN SPECTRUM WILL BE GENERATED USING
     * EXCITED-STATE GRADIENTS")')
        do i = 1, 80
          line(i:i) = ' '
        end do
        write (*,'(" NAME OF DAT FILE?")')
        read (*,'(a80)') line
        open (12,file=line,status="old",iostat=is)
        do
          read (12,'(a80)',end=300) line
          if (line(1:8).eq.' $GRAD  ') exit
        end do
C
        write (*,'(" $GRAD KEYWORD IS FOUND. START READING GRADIENT")')
        read (12,*)
        do i = 1, nat
          read (12,*,end=100) cdum,dum,(grad(j,i),j=1,3)
        end do
        close (12)
C       do i = 1, nat
C         write (*,'(i3,3f20.10)') i,(grad(j,i),j=1,3)
C       end do
      end if 
C
      if (scal.ne.1.0d+00) then
        call dscal(ncomp,scal,freq,1)
      end if
C
C     CALCULATE RAMAN ACTIVITY WITH NATURAL LIGHT (DEPOLAR (U))
C     NOTE THAT GAMESS SHOWS ACTIVITY OF PLANE-POLARIZED LIGHT BUT
C     MOLDEN SEEMS WANT TO USE THAT OF NATURAL INCIDENT LIGHT.
C
      if (hw.eq.1.0d+06) then
        hw = 2.0d+01 !! half-width in cm-1
        write (*,'(" HALF-WIDTH  =",f8.1," cm-1 (BUILT-IN)")') hw
      else
        write (*,'(" HALF-WIDTH  =",f8.1," cm-1 (USER DEFINED)")') hw
      end if
      if (raman) then
      if (temperature.eq.1.0d+06) then
        temperature = 2.9815d+02
        write (*,'(" TEMPERATURE =",f8.2," K    (BUILT-IN)")')
     *    temperature
      else
        write (*,'(" TEMPERATURE =",f8.2," K    (USER DEFINED)")')
     *    temperature
      end if
      if (f0.eq.19430.0d+00) then
        write (*,'(" WAVE NUMBER =",f8.2," cm-1 (BUILT-IN)")')
     *    f0
      else
        write (*,'(" WAVE NUMBER =",f8.2," cm-1 (USER DEFINED)")')
     *    f0
      end if
      end if
      write (*,*)
      if (.not.loren) then
        sigma = hw/(2.0d+00*sqrt(2.0d+00*log(2.0d+00)))
      end if
      !! See the following link for half-width of Lorentzian function
      !! http://mathworld.wolfram.com/LorentzianFunction.html
      hw = hw*0.5d+00 
      pi = 3.14159265359d+00
      if (raman) then
        units = 0.52917724924**2
        call dscal(6*ncomp,units,alphader,1)
        ki0 = 1.0d+00
C       f0 = 19430.0d+00
        h = 6.62606957d-34
        kb = 1.3806488d-23
        slight = 2.99792458d+08
        const = 4.799215664d-11*slight/temperature !! = hc/(kb*temperature)
        const = const*1.0d+02 !! convert frequency (cm-1) to m-1
        rammax1 = 0.0d+00
        rammax2 = 0.0d+00
        open (15,file="alpder_proj.dat")
        write (15,'("# number of mode, frequency (cm-1), xx, yy, zz,"
     *    ," xy, yz, zx(xz)")') 
        do i = 1, ncomp
          axx = ddot(ncomp,vec(1,i),1,alphader(1,1),6)
          axy = ddot(ncomp,vec(1,i),1,alphader(2,1),6)
          ayy = ddot(ncomp,vec(1,i),1,alphader(3,1),6)
          axz = ddot(ncomp,vec(1,i),1,alphader(4,1),6)
          ayz = ddot(ncomp,vec(1,i),1,alphader(5,1),6)
          azz = ddot(ncomp,vec(1,i),1,alphader(6,1),6)
          write (15,'(I4,F10.2,6(X,F20.10))') i,freq(i),axx,ayy,azz,
     *      axy,ayz,axz
          aprim2 = ((axx+ayy+azz)/3.0d+00)**2
          gprim2 = ((axx-ayy)**2 + (ayy-azz)**2 + (azz-axx)**2
     *           + 6.0d+00*(axy*axy + axz*axz + ayz*ayz))/2.0d+00
C         natural light
          rint = 4.5d+01*aprim2 + 1.3d+01*gprim2
          ramanact(1,i) = rint
          denom = 4.5d+01*aprim2 + 7.0d+00*gprim2
          depu(i) = 0.0d+00
          if (denom.gt.0.0d+00) depu(i) = 6.0d+00*gprim2 / denom
          ramanint(1,i) = ki0*(f0-freq(i))**4/
     *      (1.0d+00*freq(i)*(1.0d+00-exp(-freq(i)*const)))*rint
C         plane-polarized light (compatible with GAMESS and GAUSSIAN)
          rint = 4.5d+01*aprim2 + 7.0d+00*gprim2
          if (.not.alpder) ramanact(2,i) = rint
          ramanint(2,i) = ki0*(f0-freq(i))**4/
     *      (1.0d+00*freq(i)*(1.0d+00-exp(-freq(i)*const)))*rint
          if (freq(i).eq.0.0d+00) then
            ramanint(1,i) = 0.0d+00
            ramanint(2,i) = 0.0d+00
          end if
        end do
        close (15)
C
C       NORMALIZE RAMAN INTENSITY
C       NORMALIZATION IS IMPORTANT BECAUSE ABSOLUTE RAMAN INTENSITY
C       CALCULATED ABOVE DOES NOT INCLUDE CONVERSION FACTORS AND SO ON.
C       IGNORE ROTATION AND TRANSLATION MODES
C
        rammax1 = 0.0d+00
        rammax2 = 0.0d+00
        skip = .false.
        if (excsus) then
          fqprev = 1.0d+06
          skip = .true.
        end if
        do i = 1, ncomp-nrt
          nn = ind(i)
          if (nn.lt.0) then
            if (excimg) cycle
            nn = -nn
          end if
          if (excsus.and.skip) then
            fqcurr = freq(nn)
            if (fqprev.lt.fqcurr.and.fqcurr.gt.1.0d+00) then
              skip = .false.
              nskip = nn-1
            end if
            fqprev = freq(nn)
          end if
          if (skip) cycle
          if (ramanint(1,nn).gt.rammax1) rammax1 = ramanint(1,nn)
          if (ramanint(2,nn).gt.rammax2) rammax2 = ramanint(2,nn)
        end do
        if (excsus) then
          write (*,'(x,i2," MODES ARE SKIPPED FOR INTENSITY",
     *      " CALCULATION")') nskip
        end if
        rammax1 = 1.0d+00 / rammax1
        rammax2 = 1.0d+00 / rammax2
        call dscal(ncomp,rammax1,ramanint(1,1),2)
        call dscal(ncomp,rammax2,ramanint(2,1),2)
        open (12,file="raman_freq.dat")
        write (12,'("# number of mode, frequency (cm-1), act_natural,"
     *    ," act_polarized, int_natural, int-polarized")')
        write (12,'("# raman intensity is normalized")')
        if (nstart.ne.1) then
          do i = 1, nstart-1
            write (12,'(i4,f10.2,4(x,f20.10))') i,freq(i),ramanact(1,i),
     *        ramanact(2,i),ramanint(1,i),ramanint(2,i)
          end do
          do i = nlast+1, ncomp
            write (12,'(i4,f10.2,4(x,f20.10))') i,freq(i),ramanact(1,i),
     *        ramanact(2,i),ramanint(1,i),ramanint(2,i)
          end do
        else
          do i = 1, ncomp
            write (12,'(i4,f10.2,4(x,f20.10))') i,freq(i),ramanact(1,i),
     *        ramanact(2,i),ramanint(1,i),ramanint(2,i)
          end do
        end if
        close (12)
C
C       DAMP RAMAN INTENSITY
C       REFERENCE OF GAUSSIAN FUNCTION
C       http://mathworld.wolfram.com/GaussianFunction.html
C
        ramintout = 0.0d+00
        ramactout = 0.0d+00
        do i = 1, ncomp-nrt
          nn = ind(i)
          if (nn.lt.0) then
            if (excimg) cycle
            nn = -nn
          end if
          !! USE LORENTZIAN FUNCTION?
          if (loren) then
          if (ramanint(1,nn).gt.1.0d-10) then
            if (excsus.and.nn.le.nskip) then
              continue
            else 
              do indf = 1, 5000
                ramintout(1,indf) = ramintout(1,indf) + 1.0d+00*hw/
     *          (((dble(indf)-freq(nn))**2 + hw**2)*pi) * ramanint(1,nn)
              end do
            end if
          end if
          if (ramanint(2,nn).gt.1.0d-10) then
            if (excsus.and.nn.le.nskip) then
              continue
            else
              do indf = 1, 5000
                ramintout(2,indf) = ramintout(2,indf) + 1.0d+00*hw/
     *          (((dble(indf)-freq(nn))**2 + hw**2)*pi) * ramanint(2,nn)
              end do
            end if
          end if
          if (ramanact(1,nn).gt.1.0d-10) then
            do indf = 1, 5000
              ramactout(1,indf) = ramactout(1,indf) + 1.0d+00*hw/
     *          (((dble(indf)-freq(nn))**2 + hw**2)*pi) * ramanact(1,nn)
            end do
          end if
          if (ramanact(2,nn).gt.1.0d-10) then
            do indf = 1, 5000
              ramactout(2,indf) = ramactout(2,indf) + 1.0d+00*hw/
     *          (((dble(indf)-freq(nn))**2 + hw**2)*pi) * ramanact(2,nn)
            end do
          end if
          else
          !! USE GAUSSIAN FUNCTION
          if (ramanint(1,nn).gt.1.0d-10) then
            do indf = 1, 5000
              ramintout(1,indf) = ramintout(1,indf) + ramanint(1,nn)
     *          *exp(-((dble(indf)-freq(nn))**2)/(2.0d+00*sigma*sigma))
            end do
          end if
          if (ramanint(2,nn).gt.1.0d-10) then
            do indf = 1, 5000
              ramintout(2,indf) = ramintout(2,indf) + ramanint(2,nn)
     *          *exp(-((dble(indf)-freq(nn))**2)/(2.0d+00*sigma*sigma))
            end do
          end if
          if (ramanact(1,nn).gt.1.0d-10) then
            do indf = 1, 5000
              ramactout(1,indf) = ramactout(1,indf) + ramanact(1,nn)
     *          *exp(-((dble(indf)-freq(nn))**2)/(2.0d+00*sigma*sigma))
            end do
          end if
          if (ramanact(2,nn).gt.1.0d-10) then
            do indf = 1, 5000
              ramactout(2,indf) = ramactout(2,indf) + ramanact(2,nn)
     *          *exp(-((dble(indf)-freq(nn))**2)/(2.0d+00*sigma*sigma))
            end do
          end if
          end if
        end do
        if (loren) then
          call dscal(2*5000,hw*pi,ramintout,1)
          call dscal(2*5000,hw*pi,ramactout,1)
        else
C        call dscal(2*5000,1.0d+00/(sigma*sqrt(2.0d+00*pi)),ramintout,1)
C        call dscal(2*5000,1.0d+00/(sigma*sqrt(2.0d+00*pi)),ramactout,1)
        end if
        open (10,file="raman_int.dat",status="replace")
        open (11,file="raman_act.dat",status="replace")
        hwd = hw*2.0d+00
        if (loren) then
        write (10,'("# Lorentzian function with hw =",f8.1," cm-1")')hwd
        write (11,'("# Lorentzian function with hw =",f8.1," cm-1")')hwd
        else
        write (10,'("# Gaussian   function with hw =",f8.1," cm-1")')hwd
        write (11,'("# Gaussian   function with hw =",f8.1," cm-1")')hwd
        end if
        if (rammax1.lt.1.0d-10.or.rammax2.lt.1.0d-10) aaa=.true.
        do indf = 1, 5000
          write (10,'(i5,2(x,f20.10))')
     *      indf,ramintout(1,indf),ramintout(2,indf)
          if (aaa) then
            write (11,'(i5,2(x,e20.10))')
     *        indf,ramactout(1,indf),ramactout(2,indf)
          else
            write (11,'(i5,2(x,f20.10))')
     *        indf,ramactout(1,indf),ramactout(2,indf)
          end if
        end do
        close (10)
        close (11)
      end if
C
C     DAMP IR INTENSITY
C
      irout = 0.0d+00
      pi = 3.14159265359d+00
      do i = 1, ncomp-nrt
        nn = ind(i)
        if (nn.lt.0) then
          if (excimg) cycle
          nn = -nn
        end if
        !! USE LORENTZIAN FUNCTION?
        if (irint(nn).gt.1.0d-10) then
          if (loren) then
            do indf = 1, 5000
              irout(indf) = irout(indf) + 1.0d+00*hw/
     *          (((dble(indf)-freq(nn))**2 + hw**2)*pi) * irint(nn)
            end do
          else
            do indf = 1, 5000
              irout(indf) = irout(indf) + irint(nn)
     *          *exp(-((dble(indf)-freq(nn))**2)/(2.0d+00*sigma*sigma))
            end do
          end if
        end if
      end do
      if (loren) then
        call dscal(5000,hw*pi,irout,1)
      else
C       call dscal(5000,1.0d+00/(sigma*sqrt(2.0d+00*pi)),irout,1)
      end if
      open (12,file="ir_int.dat",status="replace")
      if (loren) then
        write (12,'("# Lorentzian function with hw =",f8.1," cm-1")')hw
      else
        write (12,'("# Gaussian   function with hw =",f8.1," cm-1")')hw
      end if
      do indf = 1, 5000
        write (12,'(i5,(x,f20.10))') indf,irout(indf)
      end do
      close (12)
C
      normalize = .true.
      normalize = .false.
      if (normalize) then
        do i = 1, ncomp
          tmp = 0.0d+00
          do jat = 1, nat
            do jk = 1, 3
              j = (jat-1)*3 + jk
              tmp = tmp + vec(j,i)**2!!*zmass(jat)
            end do
          end do
          fac = 1.0d+00/sqrt(tmp)
          do jat = 1, nat
            do jk = 1, 3
              j = (jat-1)*3 + jk
C             vec(j,i) = vec(j,i)*sqrt(zmass(jat))*fac
              vec(j,i) = vec(j,i)*fac
            end do
          end do
C         write (*,'(I3,x,F20.10)') i,fac
        end do
      end if
  
C     do i = 1, nn
C       do j = 1, ncomp
C         write (12,'(20x,5f12.8)') (vec(j,k),k=5*(i-1)+1,5*(i-1)+5)
C       end do
C       write (12,*)
C     end do
C     do j = 1, ncomp
C       if (mm.eq.1) then
C         write (12,'(20x,1f12.8)') vec(j,5*(i-1)+1)
C       else if (mm.eq.2) then
C         write (12,'(20x,2f12.8)') (vec(j,k),k=5*(i-1)+1,5*(i-1)+2)
C       else if (mm.eq.3) then
C         write (12,'(20x,3f12.8)') (vec(j,k),k=5*(i-1)+1,5*(i-1)+3)
C       else if (mm.eq.4) then
C         write (12,'(20x,4f12.8)') (vec(j,k),k=5*(i-1)+1,5*(i-1)+4)
C       end if
C     end do

      open (15,file="dipder_proj.dat")
      write (15,'("# number of mode, frequency (cm-1), x, y, z,
     *  intensity")')
      do i = 1, ncomp
        ddx = ddot(ncomp,vec(1,i),1,dipder(1,1),3)
        ddy = ddot(ncomp,vec(1,i),1,dipder(2,1),3)
        ddz = ddot(ncomp,vec(1,i),1,dipder(3,1),3)
        tmp = ddx*ddx + ddy*ddy + ddz*ddz
        if (i.lt.nstart) then
        write (15,'(I4,F10.2,4(X,F20.10))') i,-freq(i),ddx,ddy,ddz,tmp
        else
        write (15,'(I4,F10.2,4(X,F20.10))') i,freq(i),ddx,ddy,ddz,tmp
        end if
      end do
      close (15)
C
C     CONVERT THE UNIT OF COORDINATE FROM BOHR TO ANGSTROM
C
      call dscal(3*nat,0.529177211d+00,coord,1)
C
C     CONVERT THE UNIT OF IR INTENSITY TO KM/MOLE
C
      call dscal(3*nat,4.2255d+01,irint,1)
C
C     WRITE ALL INFORMATION TO A GAUSSIAN FORMAT FILE
C
      open (10,file=filego,status="new",iostat=is)
      if (is.eq.10.and..not.overwrite) then
        close (10)
        write (*,*) "MAYBE test.out ALREADY EXISTS?"
        write (*,*) "IN ORDER NOT TO ACCIDENTALLY OVERWRITE test.out, ",
     *              "PROGRAM IS MOMENTARILY SUSPENDED"
        write (*,*)
        write (*,*) "GIVE A DIFFERENT NAME, OR TYPE test.out AGAIN"
        read (*,*) filego
        open (10,file=filego)
      else if (is.ne.0.and..not.overwrite) then
        write (*,*) "SOMETHING WENT WRONG WHILE OPENING test.out"
      end if
C
      write (10,'(" Entering Gaussian System, Link 0=g09")')
      write (10,'(" ******************************************")')
      write (10,'(" Gaussian 09:  AM64L-G09RevC.01 23-Sep-2011")')
      write (10,*)
      write (10,'(" GradGradGradGradGradGradGradGradGradGradGradGradGrad
     *GradGradGradGradGrad")')
      write (10,'(" GradGradGradGradGradGradGradGradGradGradGradGradGrad
     *GradGradGradGradGrad")')
      write (10,*)
      write (10,'("                         Standard orientation:")')
      write (10,'(" ----------------------------------------------------
     *-----------------")')
      write (10,'(" Center     Atomic      Atomic             Coordinate
     *s (Angstroms)")')
      write (10,'(" Number     Number       Type             X          
     * Y           Z    ")')
      write (10,'(" ----------------------------------------------------
     *-----------------")')
      ndum = 0
      do i = 1, nat
        nzan = zan(i)
        write (10,'(3x,i4,8x,i3,11x,i1,4x,3f12.6)') i,nzan,ndum,
     *    (coord(j,i),j=1,3)
      end do
      write (10,'(" ----------------------------------------------------
     *-----------------")')
      write (10,'(i6," basis functions,",i6," primitive gaussians",i6,
     * " cartesian basis functions")') l1,l1,l1
      write (10,'(i6," alpha electrons",3x,i6," beta electrons")') l1,l1
      write (10,'("       nuclear repulsion energy ",f20.10,
     * " Hartrees.")') 0.0d+00
      write (10,'(" Initial guess orbital symmetries:")')
      write (10,'("       Occupied  (A)")')
      write (10,'("       Virtual   (A)")')
      write (10,'(" The electronic state of the initial guess is 1-A."
     * )')
      write (10,'(" SCF Done:  E(RHF) =",f16.10,"     A.U. after   ",i2,
     *  " cycles")') 0.0d+00,0
      write (10,'("             Convg  =    0.0000D-00",
     *"             -V/T =  2.0000")')
C
      write (10,*)
C
C     START TO WRITE FREQUENCY INFORMATION
C     NOTE THAT FIRST ROTATION AND TRANSLATION HAS TO BE REMOVED
C
      write (10,'(" Harmonic frequencies (cm**-1), IR intensities (KM/Mo
     *le), Raman scattering")')
      write (10,'(" activities (A**4/AMU), depolarization ratios for pla
     *ne and unpolarized")')
      write (10,'(" incident light, reduced masses (AMU), force constant
     *s (mDyne/A),")')
      write (10,'(" and normal coordinates:")')
      nn = floor(dble(ncomp-nrt)/3.0d+00)
      mm = mod(ncomp-nrt,3)
      do i = 1, nn
        nt1 = (i-1)*3 + 1
        nt2 = nt1 + 1
        nt3 = nt2 + 1
        n1 = ind(nt1)
        if (n1.lt.0) then
          n1 = abs(n1)
          freq(n1) = -freq(n1)
        end if
        n2 = ind(nt2)
        if (n2.lt.0) then
          n2 = abs(n2)
          freq(n2) = -freq(n2)
        end if
        n3 = ind(nt3)
        if (n3.lt.0) then
          n3 = abs(n3)
          freq(n3) = -freq(n3)
        end if
        write (10,'(18x,i4,19x,i4,19x,i4)') nt1,nt2,nt3
        write (10,'(21x,a,22x,a,22x,a)') "A","A","A"
        write (10,'(" Frequencies -- ",f10.4,13x,f10.4,13x,f10.4)')
     *    freq(n1),freq(n2),freq(n3)
        write (10,'(" Red. masses -- ",f10.4,13x,f10.4,13x,f10.4)')
     *    reduced(n1),reduced(n2),reduced(n3)
        write (10,'(" Frc consts  -- ",f10.4,13x,f10.4,13x,f10.4)')
     *    1.0d+00,1.0d+00,1.0d+00
        write (10,'(" IR Inten    -- ",f10.4,13x,f10.4,13x,f10.4)')
     *    irint(n1),irint(n2),irint(n3)
        if (raman) then
        write (10,'(" Raman Activ --" ,f11.4,12x,f11.4,12x,f11.4)')
     *    ramanact(iram,n1),ramanact(iram,n2),ramanact(iram,n3)
        write (10,'(" Depolar (P) -- ",f10.4,13x,f10.4,13x,f10.4)')
     *    depolar(n1),depolar(n2),depolar(n3)
        write (10,'(" Depolar (U) -- ",f10.4,13x,f10.4,13x,f10.4)')
     *    depu(n1),depu(n2),depu(n3)
        end if
        write (10,'("  Atom  AN      X      Y      Z        X      
     *Y      Z        X      Y      Z")')
        do j = 1, nat
          ndum = zan(j)
          nn1 = (j-1)*3 + 1
          nn2 = nn1 + 1
          nn3 = nn2 + 1
          write (10,'(2x,i4,x,i3,3(4x,f5.2,2x,f5.2,2x,f5.2))')
     *      j,ndum,vec(nn1,n1),vec(nn2,n1),vec(nn3,n1),vec(nn1,n2),
     *      vec(nn2,n2),vec(nn3,n2),vec(nn1,n3),vec(nn2,n3),vec(nn3,n3)
        end do
      end do
      if (mm.eq.1) then
        nt1 = (i-1)*3 + 1
        n1 = ind(nt1)
        if (n1.lt.0) then
          n1 = abs(n1)
          freq(n1) = -freq(n1)
        end if
        write (10,'(18x,i4)') nt1
        write (10,'(21x,a)') "A"
        write (10,'(" Frequencies -- ",f10.4)') freq(n1)
        write (10,'(" Red. masses -- ",f10.4)') reduced(n1)
        write (10,'(" Frc consts  -- ",f10.4)') 1.0d+00
        write (10,'(" IR Inten    -- ",f10.4)') irint(n1)
        if (raman) then
        write (10,'(" Raman Activ --" ,f11.4)') ramanact(iram,n1)
        write (10,'(" Depolar (P) -- ",f10.4)') depolar(n1)
        write (10,'(" Depolar (U) -- ",f10.4)') depu(n1)
        end if
        write (10,'("  Atom  AN      X      Y      Z")')
        do j = 1, nat
          ndum = zan(j)
          nn1 = (j-1)*3 + 1
          nn2 = nn1 + 1
          nn3 = nn2 + 1
          write (10,'(2x,i4,x,i3,4x,f5.2,2x,f5.2,2x,f5.2)')
     *      j,ndum,vec(nn1,n1),vec(nn2,n1),vec(nn3,n1)
        end do
      else if (mm.eq.2) then
        nt1 = (i-1)*3 + 1
        nt2 = nt1 + 1
        n1 = ind(nt1)
        if (n1.lt.0) then
          n1 = abs(n1)
          freq(n1) = -freq(n1)
        end if
        n2 = ind(nt2)
        if (n2.lt.0) then
          n2 = abs(n2)
          freq(n2) = -freq(n2)
        end if
        write (10,'(18x,i4,19x,i4)') nt1,nt2
        write (10,'(21x,a,22x,a)') "A","A"
        write (10,'(" Frequencies -- ",f10.4,13x,f10.4)')
     *    freq(n1),freq(n2)
        write (10,'(" Red. masses -- ",f10.4,13x,f10.4)')
     *    reduced(n1),reduced(n2)
        write (10,'(" Frc consts  -- ",f10.4,13x,f10.4)')
     *    1.0d+00,1.0d+00
        write (10,'(" IR Inten    -- ",f10.4,13x,f10.4)')
     *    irint(n1),irint(n2)
        if (raman) then
        write (10,'(" Raman Activ --" ,f11.4,12x,f11.4)')
     *    ramanact(iram,n1),ramanact(iram,n2)
        write (10,'(" Depolar (P) -- ",f10.4,13x,f10.4)')
     *    depolar(n1),depolar(n2)
        write (10,'(" Depolar (U) -- ",f10.4,13x,f10.4)')
     *    depu(n1),depu(n2)
        end if
        write (10,'("  Atom  AN      X      Y      Z        X      
     *Y      Z")')
        do j = 1, nat
          ndum = zan(j)
          nn1 = (j-1)*3 + 1
          nn2 = nn1 + 1
          nn3 = nn2 + 1
          write (10,'(2x,i4,x,i3,2(4x,f5.2,2x,f5.2,2x,f5.2))')
     *      j,ndum,vec(nn1,n1),vec(nn2,n1),vec(nn3,n1),vec(nn1,n2),
     *      vec(nn2,n2),vec(nn3,n2)
        end do
      end if
C
      write (10,*)
      write (10,'(" Temperature   298.150 Kelvin.  Pressure   1.00000
     * Atm.")')
      write (10,*)
      write (10,'(" GradGradGradGradGradGradGradGradGradGradGradGradGrad
     *GradGradGradGradGrad")')
      write (10,'(" Step number   1 out of a maximum of    2")')
      write (10,'(" GradGradGradGradGradGradGradGradGradGradGradGradGrad
     *GradGradGradGradGrad")')
      write (10,*)
      write (10,'(" Normal termination of Gaussian 09 at")')
C
      close (10)
C
C     CHECK THE QUALITY OF HESSIAN FROM ROTATION AND TRANSLATION MODES
C
      write (*,*) "FOLLOWING MODES ARE ASSIGNED AS ROTATION AND
     * TRANSLATION"
      write (*,*)
      tmpfreq = 0.0d+00
      high = .false.
      aaa  = .false.
      do i = nstart, nlast
        write (*,'(" MODE NO.",i3," ...",f7.3," cm-1")') i,freq(i)
        tmpfreq = tmpfreq + freq(i)
        if (freq(i).gt.1.0d-01) aaa  = .true.
        if (freq(i).gt.1.0d+00) high = .true.
      end do
      tmpfreq = tmpfreq/dble(nrt0) !! 6.0d+00
      write (*,*)
      write (*,'(" AVERAGE OF ",i1," FREQUENCIES ...",F7.3," cm-1")')
     * nrt0,tmpfreq
      if (high) then
      write (*,*)
      write (*,'(" AT LEAST ONE FREQUENCY IS HIGHER THAN 1.0 cm-1")')
      write (*,'(" PLEASE CHECK FREQUENCIES MANUALLY")')
      else if (aaa) then
      write (*,'(" ACCEPTABLE ... ALL FREQUENCIES ARE LESS THAN
     * 1.0 cm-1")')
      else
      write (*,'(" OK ... ALL FREUQNEICES ARE LESS THAN 0.1 cm-1")')
      end if
C
      if (rrs) then
        irint=0.0d+00
        do i = 1, ncomp
          irint(i) = ddot(ncomp,vec(1,i),1,grad,1)
          irint(i) = irint(i)*irint(i)
        end do
C
        rammax1 = 0.0d+00
        skip = .false.
        if (excsus) then
          fqprev = 1.0d+06
          skip = .true.
        end if
        do i = 1, ncomp-nrt
          nn = ind(i)
          if (nn.lt.0) then
            if (excimg) cycle
            nn = -nn
          end if
          if (excsus.and.skip) then
            fqcurr = freq(nn)
            if (fqprev.lt.fqcurr.and.fqcurr.gt.1.0d+00) then
              skip = .false.
              nskip = nn-1
            end if
            fqprev = freq(nn)
          end if
          if (skip) cycle
          if (irint(nn).gt.rammax1) rammax1 = irint(nn)
        end do
        if (excsus) then
          write (*,'(x,i2," MODES ARE SKIPPED FOR INTENSITY",
     *      " CALCULATION")') nskip
        end if
        rammax1 = 1.0d+00 / rammax1
        call dscal(ncomp,rammax1,irint,1)
C
        open (12,file="raman_freq_rrs.dat")
        write (12,'("# number of mode, frequency (cm-1), resonance-"
     *    ,"raman intensity")')
        write (12,'("# raman intensity is normalized")')
        if (nstart.ne.1) then
          do i = 1, nstart-1
            write (12,'(i4,f10.2,x,f20.10)') i,freq(i),irint(i)
          end do
          do i = nlast+1, ncomp
            write (12,'(i4,f10.2,x,f20.10)') i,freq(i),irint(i)
          end do
        else
          do i = 1, ncomp
            write (12,'(i4,f10.2,x,f20.10)') i,freq(i),irint(i)
          end do
        end if
        close (12)
C
        irout = 0.0d+00
        pi = 3.14159265359d+00
        do i = 1, ncomp-nrt
          nn = ind(i)
          if (nn.lt.0) then
            if (excimg) cycle
            nn = -nn
          end if
          !! USE LORENTZIAN FUNCTION?
          if (irint(nn).gt.1.0d-10) then
            if (loren) then
              do indf = 1, 5000
                irout(indf) = irout(indf) + 1.0d+00*hw/
     *            (((dble(indf)-freq(nn))**2 + hw**2)*pi) * irint(nn)
              end do
            else
              do indf = 1, 5000
                irout(indf) = irout(indf) + irint(nn)
     *           *exp(-((dble(indf)-freq(nn))**2)/(2.0d+00*sigma*sigma))
              end do
            end if
          end if
        end do
        if (loren) then
          call dscal(5000,hw*pi,irout,1)
        else
C         call dscal(5000,1.0d+00/(sigma*sqrt(2.0d+00*pi)),irout,1)
        end if
        open (12,file="raman_int_rrs.dat",status="replace")
        if (loren) then
         write (12,'("# Lorentzian function with hw =",f8.1," cm-1")')hw
        else
         write (12,'("# Gaussian   function with hw =",f8.1," cm-1")')hw
        end if
        do indf = 1, 5000
          write (12,'(i5,(x,f20.10))') indf,irout(indf)
        end do
        close (12)
      end if
C
C     write (*,*) "finished everything"
C     write (*,*) "deallocate everything"
 
      deallocate (dipder,vec)
      deallocate (zan,coord,freq,reduced,irint,zmass,ind,irout)
      if (raman) deallocate (ramanact,depolar,alphader,depu,ramanint,
     *  ramintout,ramactout)
      if (rrs) deallocate (grad)

      write (*,*)
      write (*,*) "  ------------------"
      write (*,*) "  NORMAL TERMINATION"
      write (*,*) "  ------------------"
      write (*,*)
      write (*,*) "  FOLLOWING FILES ARE GENERATED"
C     write (*,*) "  1. test.out"
      nc = len_trim(filego)
      write (*,*) "  1. ",filego(1:nc)
      write (*,*) "  2. ir_int.dat"
      write (*,*) "  3. dipder_proj.dat"
      if (raman) then
      write (*,*) "  4. raman_freq.dat"
      write (*,*) "  5. raman_act.dat"
      write (*,*) "  6. raman_int.dat"
      write (*,*) "  7. alpder_proj.dat"
      else if (rrs) then
      write (*,*) "  4. raman_freq_rrs.dat"
      write (*,*) "  5. raman_int_rrs.dat"
      end if
      write (*,*)
      write (*,*) "  IF YOU NEED BRIEF DETAILS OF OUTPUT FILES, PLEASE
     * TRY WITH ""-H"" OPTION"
      if (nstart.ne.1.and.excimg) then
      write (*,*)
      write (*,*) "  CONTRIBUTIONS OF IMAGINARY MODES TO PLOT FILES ARE
     * REMOVED"
      write (*,*) "  IF YOU WANT TO INCLUDE THEM, PLEASE
     * TRY WITH ""-i"" OPTION"
      end if
      write (*,*)

      stop

  400 continue
      write (*,*)
      write (*,*)"IR RELATING OUTPUT FILES (3 FILES)"
      write (*,*)"  1. test.out (default name)"
      write (*,*)"  - Molden or GaussView can read this file, and"
      write (*,*)"    visualize gemometry and vibrational vectors."
      write (*,*)"    File name can be altered with -go option"
      write (*,*)
      write (*,*)"  2. ir_int.dat"
      write (*,*)"  - IR intensity for drawing a plot"
      write (*,*)"    1st column: frequency (cm-1)"
      write (*,*)"    2nd column: IR intensity (DEBYE^2/AMU-ANGSTROM^2)"
      write (*,*)
      write (*,*)"  3. dipder_proj.dat"
      write (*,*)"  - Projected dipole moment derivatives to each
     * vibration mode"
      write (*,*)"    1st column: index of vibrational vector"
      write (*,*)"    2nd column: frequency (cm-1)"
      write (*,*)"    3rd column: d(mu_x)/dQ"
      write (*,*)"    4th column: d(mu_y)/dQ"
      write (*,*)"    5th column: d(mu_z)/dQ"
      write (*,*)"    6th column: IR intensity (DEBYE^2/AMU-ANGSTROM^2)"
      write (*,*)
      write (*,*)"RAMAN RELATING OUTPUT FILES (4 FILES)"
      write (*,*)"  4. raman_freq.dat"
      write (*,*)"  - Raman activity and intensity of each vibration
     * mode"
      write (*,*)"    1st column: index of vibrational vector"
      write (*,*)"    2nd column: frequency (cm-1)"
      write (*,*)"    3rd column: Raman activity  of natural light"
      write (*,*)"    4th column: Raman activity  of plane-polarized
     * light"
      write (*,*) "               (consistent with GAMESS)"
      write (*,*)"    5th column: Raman intensity of natural light"
      write (*,*)"    6th column: Raman intensity of plane-polarized
     * light"
      write (*,*)"  - Unit of Raman activity is ANGSTROM^4/AMU"
      write (*,*)"  - Note that intensity is normalized"
      write (*,*)
      write (*,*)"  5. raman_act.dat"
      write (*,*)"  - Raman activity for drawing a plot"
      write (*,*)"    1st column: frequency (cm-1)"
      write (*,*)"    2nd column: Raman activity  of natural light"
      write (*,*)"    3rd column: Raman activity  of plane-polarized
     * light"
      write (*,*)
      write (*,*)"  6. raman_int.dat"
      write (*,*)"  - Raman activity for drawing a plot"
      write (*,*)"    1st column: frequency (cm-1)"
      write (*,*)"    2nd column: Raman intensity of natural light"
      write (*,*)"    3rd column: Raman intensity of plane-polarized
     * light"
      write (*,*)
      write (*,*)"  7. alpder_proj.dat"
      write (*,*)"  - Projected alpha polarizability derivatives to each
     * vibration mode"
      write (*,*)"    1st column: index of vibrational vector"
      write (*,*)"    2nd column: frequency (cm-1)"
      write (*,*)"    3rd column: d(alpha_xx)/dQ"
      write (*,*)"    4th column: d(alpha_yy)/dQ"
      write (*,*)"    5th column: d(alpha_zz)/dQ"
      write (*,*)"    6th column: d(alpha_xy)/dQ"
      write (*,*)"    7th column: d(alpha_yz)/dQ"
      write (*,*)"    8th column: d(alpha_xz)/dQ"
      write (*,*)
      write (*,*)"OUTPUT FILES WITH -pl OPTION (1 OR 2 FILES)"
      write (*,*)"  8. dipder_plt.dat"
      write (*,*)"  - Projected dipole moment derivatives for plot"
      write (*,*)"    1st column: index of vibrational vector"
      write (*,*)"    2nd column: d(mu_x)/dQ"
      write (*,*)"    3rd column: d(mu_y)/dQ"
      write (*,*)"    4th column: d(mu_z)/dQ"
      write (*,*)"    5th column: IR intensity (DEBYE^2/AMU-ANGSTROM^2)"
      write (*,*)
      write (*,*)"  9. alpder_plt.dat"
      write (*,*)"  - Projected alpha polarizability derivatives for
     * plot"
      write (*,*)"    1st column: index of vibrational vector"
      write (*,*)"    2nd column: d(alpha_xx)/dQ"
      write (*,*)"    3rd column: d(alpha_yy)/dQ"
      write (*,*)"    4th column: d(alpha_zz)/dQ"
      write (*,*)"    5th column: d(alpha_xy)/dQ"
      write (*,*)"    6th column: d(alpha_yz)/dQ"
      write (*,*)"    7th column: d(alpha_xz)/dQ"
      write (*,*)

      stop

 500  continue
      write (*,*)
      write (*,*)"OPTIONS"
      write (*,*)"  -a"
      write (*,*)"     Use absolute values for -pl"
      write (*,*)"  -e"
      write (*,*)"     Exclude suspicious frequencies from Raman"
      write (*,*)"     intensity calculation.  Suspicious here means"
      write (*,*)"     frequencies higher than 1.0 cm-1 but ""below"""
      write (*,*)"     rotational or translational mode"
      write (*,*)"  -g"
      write (*,*)"     Use Gaussian function instead of Lorentzian"
      write (*,*)"     function for plot files.  Note that the top of"
      write (*,*)"     the Gaussian used here is 1, so the coefficient"
      write (*,*)"     (1/(sigma sqrt(2pi))) is omitted"
      write (*,*)"  -go"
      write (*,*)"     Name of Gaussian-type output."
      write (*,*)"     Default name is test.out (= -go test.out)"
      write (*,*)"  -h"
      write (*,*)"     Show this help"
      write (*,*)"  -H"
      write (*,*)"     Show brief explanations of output files"
      write (*,*)"  -i"
      write (*,*)"     Include imaginary frequencies in calculating"
      write (*,*)"     plot files"
      write (*,*)"  -l"
      write (*,*)"     Use Lorentzian function for plot files.  This is"
      write (*,*)"     the default choice.  Note that the top of the"
      write (*,*)"     Lorentzian function used here is 1, so the"
      write (*,*)"     coefficient (1/pi) is omitted"
      write (*,*)"  -n"
      write (*,*)"     Use Raman activity and intensity of natural"
      write (*,*)"     indient light (default)"
      write (*,*)"  -o"
      write (*,*)"     Overwrite test.out anyway"
      write (*,*)"  -p"
      write (*,*)"     Use Raman activity and intensity of"
      write (*,*)"     plane-polarized light."
      write (*,*)"  -s"
      write (*,*)"     Use squared values for -pl (default)"
      write (*,*)"  -t [value]"
      write (*,*)"     Change temperature in calculating Raman"
      write (*,*)"     intensity.  The default value is 298.15 K."
      write (*,*)"  -w"
      write (*,*)"     Wave number (cm-1) of incident light used when"
      write (*,*)"     Raman intensity is calculated.  The default"
      write (*,*)"     value is 19,430 cm-1."
      write (*,*)"  -hw [value]"
      write (*,*)"     Change half-width for generating plot files.  By"
      write (*,*)"     default, built-in half-width is set to 20 cm-1."
      write (*,*)"     This option works for -pl too."
      write (*,*)"  -pl"
      write (*,*)"     Process dipder_proj.dat and/or alpder_proj.dat"
      write (*,*)"     to generate a file suitable for plotting"
      write (*,*)"     projected dipole or alpha polarizability"
      write (*,*)"     derivatives.  Half-width can be changed by -hw"
      write (*,*)"     option, and the quantity calculated may be"
      write (*,*)"     changed with -a and -s options (I don't know"
      write (*,*)"     which is better)."
      write (*,*)"  -fmo"
      write (*,*)"     An option needed for processing FMO output"
      write (*,*)"     files"
      write (*,*)"  -rrs"
      write (*,*)"     Calculate Resonance-Raman spectra (intensity)"
      write (*,*)"     using ground-state second derivative (Hessian)"
      write (*,*)"     and excited-state first derivative (gradient)."
      write (*,*)"     A .dat file that contains such a gradient is"
      write (*,*)"     needed."
      write (*,*)"  -scal [value]"
      write (*,*)"     Frequencies will be scaled with [value]."
      write (*,*)"     For instance, -scal 0.89"
      write (*,*)

      stop

300   write (*,*) "ERROR TERMINATION"
      stop

      end program test
C
C-----------------------------------------------------------------------
C*MODULE BLAS1   *DECK DDOT
      DOUBLE PRECISION FUNCTION DDOT(N,DX,INCX,DY,INCY)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION DX(*),DY(*)
C
C     FORMS THE DOT PRODUCT OF TWO VECTORS.
C           DOT = DX(I) * DY(I)
C     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE.
C     JACK DONGARRA, LINPACK, 3/11/78.
C
      DDOT = 0.0D+00
      DTEMP = 0.0D+00
      IF(N.LE.0) RETURN
      IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20
C
C        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
C          NOT EQUAL TO 1
C
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DTEMP = DTEMP + DX(IX)*DY(IY)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      DDOT = DTEMP
      RETURN
C
C        CODE FOR BOTH INCREMENTS EQUAL TO 1
C
C
C        CLEAN-UP LOOP
C
   20 M = MOD(N,5)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        DTEMP = DTEMP + DX(I)*DY(I)
   30 CONTINUE
      IF( N .LT. 5 ) GO TO 60
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
        DTEMP = DTEMP + DX(I)*DY(I) + DX(I + 1)*DY(I + 1) +
     *   DX(I + 2)*DY(I + 2) + DX(I + 3)*DY(I + 3) + DX(I + 4)*DY(I + 4)
   50 CONTINUE
   60 DDOT = DTEMP
      RETURN
      END
C
C-----------------------------------------------------------------------
C*MODULE BLAS1   *DECK DSCAL
      SUBROUTINE DSCAL(N,DA,DX,INCX)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION DX(*)
C
C     SCALES A VECTOR BY A CONSTANT.
C           DX(I) = DA * DX(I)
C     USES UNROLLED LOOPS FOR INCREMENT EQUAL TO ONE.
C     JACK DONGARRA, LINPACK, 3/11/78.
C
      IF(N.LE.0) RETURN
      IF(INCX.EQ.1)GO TO 20
C
C        CODE FOR INCREMENT NOT EQUAL TO 1
C
      NINCX = N*INCX
      DO 10 I = 1,NINCX,INCX
        DX(I) = DA*DX(I)
   10 CONTINUE
      RETURN
C
C        CODE FOR INCREMENT EQUAL TO 1
C
C
C        CLEAN-UP LOOP
C
   20 M = MOD(N,5)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        DX(I) = DA*DX(I)
   30 CONTINUE
      IF( N .LT. 5 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
        DX(I) = DA*DX(I)
        DX(I + 1) = DA*DX(I + 1)
        DX(I + 2) = DA*DX(I + 2)
        DX(I + 3) = DA*DX(I + 3)
        DX(I + 4) = DA*DX(I + 4)
   50 CONTINUE
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      subroutine plot_dat(hw,excimg,sq,loren)
C
      implicit double precision (a-h,o-z)
C
      double precision, allocatable :: freq(:),ir(:,:),irout(:,:),
     *  raman(:,:),ramout(:,:)
      character(80) :: file1
      logical :: excimg,sq,loren,ok,dipole,alpha,square
C
      square = sq
C
      write (*,*)
      write (*,*)"WITH ""-pl"" OPTION, dipder_proj.dat AND/OR"
      write (*,*)"alpder_proj.dat ARE CONVERTED TO A FORMAT SUITABLE"
      write (*,*)"FOR PLOTTING"
      write (*,*)
      write (*,*)"GIVE ONE OF THREE NUMBERS"
      write (*,*)"  1: USE dipder_proj.dat"
      write (*,*)"  2: USE alpder_proj.dat"
      write (*,*)"  3: USE dipder_proj.dat AND alpder_proj.dat"
      write (*,*)
      read (*,'(a80)') file1
      write (*,*)
      ok = .false.
      dipole = .false.
      alpha = .false.
      if (file1(1:1).eq.'1') then
        nn = 1
        dipole = .true.
        ok = .true.
      else if (file1(1:1).eq.'2') then
        nn = 2
        alpha = .true.
        ok = .true.
      else if (file1(1:1).eq.'3') then
        nn = 3
        dipole = .true.
        alpha = .true.
        ok = .true.
      end if
C
      if (.not.ok) then
        write (*,*) "INPUT NUMBER IS NOT RECOGNIZED"
        stop
      end if
C
      ok = .true.
      if (dipole) then
        open (10,file="dipder_proj.dat",status="old",iostat=is)
        if (is.ne.0) then
          write (*,*) "FAILURE IN OPENING dipder_proj.dat"
          ok = .false.
        end if
      end if
      if (alpha) then
        open (11,file="alpder_proj.dat",status="old",iostat=is)
        if (is.ne.0) then
          write (*,*) "FAILURE IN OPENING alpder_proj.dat"
          ok = .false.
        end if
      end if
      if (.not.ok) then
        write (*,*)
        stop
      end if
C
      if (hw.eq.1.0d+06) then
        hw = 2.0d+01 !! half-width in cm-1
        write (*,'(" HALF-WIDTH  =",f8.1," cm-1 (BUILT-IN)")') hw
      else
        write (*,'(" HALF-WIDTH  =",f8.1," cm-1 (USER DEFINED)")') hw
      end if
      if (.not.loren) then
        sigma = hw/(2.0d+00*sqrt(2.0d+00*log(2.0d+00)))
      end if
      hw = hw * 0.5d+00
C
C     START FOR diper_proj.dat
C
      if (dipole) then
C
C       GET THE NUMBER OF MODES
C
        write (*,*)
        write (*,*) "READING dipder_proj.dat ..."
        write (*,*)
        nvec = 0
        nline = 0
        do
          read (10,*,end=10) file1
          nline = nline + 1
          if (file1(1:1).ne.'#') nvec = nvec + 1
        end do
   10   continue
        write (*,'("  NUMBER OF MODES =",I4)') nvec
        write (*,*)
C
C       ALLLOCATE FOR dipder_proj.dat
C
        allocate (freq(nvec))
        allocate (ir(4,nvec))
        allocate (irout(4,5000))
C
C       READ dipder_proj.dat
C       ir(1,*) :: ddx
C       ir(2,*) :: ddy
C       ir(3,*) :: ddz
C       ir(4,*) :: intensity
C
        rewind (10)
        do i = 1, nline - nvec
          read (10,*)
        end do
        nstart = 0
        do i = 1, nvec
          read (10,*) idum,freq(i),(ir(j,i),j=1,4)
          if (freq(i).lt.0.0d+00) nstart = i
        end do
        close (10)
        nstart = nstart+1
C
        irout = 0.0d+00
        pi = 3.14159265359d+00
        do i = 1, nvec
          if (i.le.nstart+5.and.i.gt.nstart-1) cycle
          ft = freq(i)
          if (ft.lt.0.0d+00) then
            if (excimg) cycle
            ft = -freq(i)
          end if
          !! USE LORENTZIAN FUNCTION?
          if (loren) then
          do j = 1, 4
            if (abs(ir(j,i)).gt.1.0d-10) then
              if (square.and.j.le.3) then
                do indf = 1, 5000
                  irout(j,indf) = irout(j,indf) + 1.0d+00*hw/
     *                (((dble(indf)-ft)**2 + hw**2)*pi)*ir(j,i)*ir(j,i)
                end do
              else
                do indf = 1, 5000
                  irout(j,indf) = irout(j,indf) + 1.0d+00*hw/
     *                (((dble(indf)-ft)**2 + hw**2)*pi)*abs(ir(j,i))
                end do
              end if
            end if
          end do
          else
          !! USE GAUSSIAN FUNCTION
          do j = 1, 4
            if (abs(ir(j,i)).gt.1.0d-10) then
              if (square.and.j.le.3) then
                do indf = 1, 5000
                  irout(j,indf) = irout(j,indf) + ir(j,i)*ir(j,i)
     *          *exp(-((dble(indf)-freq(nn))**2)/(2.0d+00*sigma*sigma))
                end do
              else
                do indf = 1, 5000
                  irout(j,indf) = irout(j,indf) + abs(ir(j,i))
     *          *exp(-((dble(indf)-freq(nn))**2)/(2.0d+00*sigma*sigma))
                end do
              end if
            end if
          end do
          end if
        end do
        if (loren) then
          call dscal(5000*4,hw*pi,irout,1)
        else
C        call dscal(5000*4,1.0d+00/(sigma*sqrt(2.0d+00*pi)),irout,1)
        end if
        open (12,file="dipder_plt.dat",status="replace")
        write (12,'("# wavenumber, x, y, z, intensity")')
        hwd = hw*2.0d+00
        if (loren) then
         write(12,'("# Lorentzian function with hw =",f8.1," cm-1")')hwd
        else
         write(12,'("# Gaussian   function with hw =",f8.1," cm-1")')hwd
        end if
        do indf = 1, 5000
          write (12,'(i5,4(x,f20.10))') indf,(irout(j,indf),j=1,4)
        end do
        close (12)
C
        deallocate (freq,ir,irout)
        write (*,*) "dipder_plt.dat IS GENERATED"
        if (square) then
        write (*,*) "SQUARED VALUES ARE PLOTTED FOR X, Y, AND Z"
        else
        write (*,*) "ABSOLUTE VALUES ARE PLOTTED FOR X, Y, AND Z"
        END IF
      end if
C
C     START FOR alpder_proj.dat
C
      if (alpha) then
C
C       GET THE NUMBER OF MODES
C
        write (*,*)
        write (*,*) "READING alpder_proj.dat ..."
        write (*,*)
        nvec = 0
        nline = 0
        do
          read (11,*,end=20) file1
          nline = nline + 1
          if (file1(1:1).ne.'#') nvec = nvec + 1
        end do
   20   continue
        write (*,'("  NUMBER OF MODES =",I4)') nvec
        write (*,*)
C
C       ALLLOCATE FOR dipder_proj.dat
C
        allocate (freq(nvec))
        allocate (raman(6,nvec))
        allocate (ramout(6,5000))
C
C       READ dipder_proj.dat
C       raman(1,*) :: xx
C       raman(2,*) :: yy
C       raman(3,*) :: zz
C       raman(4,*) :: xy
C       raman(5,*) :: yz
C       raman(6,*) :: xz
C
        rewind (11)
        do i = 1, nline - nvec
          read (11,*)
        end do
        nstart = 0
        do i = 1, nvec
          read (11,*) idum,freq(i),(raman(j,i),j=1,6)
          if (freq(i).lt.0.0d+00) nstart = i
        end do
        close (11)
        nstart = nstart+1
C
        ramout = 0.0d+00
        pi = 3.14159265359d+00
        do i = 1, nvec
          if (i.le.nstart+5.and.i.gt.nstart-1) cycle
          ft = freq(i)
          if (ft.lt.0.0d+00) then
            if (excimg) cycle
            ft = -freq(i)
          end if
          !! USE LORENTZIAN FUNCTION?
          if (loren) then
          do j = 1, 6
            if (abs(raman(j,i)).gt.1.0d-10) then
              if (square) then
                do indf = 1, 5000
                  ramout(j,indf) = ramout(j,indf) + 1.0d+00*hw/
     *                (((dble(indf)-ft)**2 + hw**2)*pi)
     *                *raman(j,i)*raman(j,i)
                end do
              else
                do indf = 1, 5000
                  ramout(j,indf) = ramout(j,indf) + 1.0d+00*hw/
     *                (((dble(indf)-ft)**2 + hw**2)*pi)*abs(raman(j,i))
                end do
              end if
            end if
          end do
          else
          !! USE GAUSSIAN FUNCTION
          do j = 1, 6
            if (abs(raman(j,i)).gt.1.0d-10) then
              if (square) then
                do indf = 1, 5000
                  ramout(j,indf) = ramout(j,indf)+ raman(j,i)*raman(j,i)
     *          *exp(-((dble(indf)-freq(nn))**2)/(2.0d+00*sigma*sigma))
                end do
              else
                do indf = 1, 5000
                  ramout(j,indf) = ramout(j,indf) + abs(raman(j,i))
     *          *exp(-((dble(indf)-freq(nn))**2)/(2.0d+00*sigma*sigma))
                end do
              end if
            end if
          end do
          end if
        end do
        if (loren) then
          call dscal(5000*6,hw*pi,ramout,1)
        else
C        call dscal(5000*6,1.0d+00/(sigma*sqrt(2.0d+00*pi)),ramout,1)
        end if
        open (12,file="alpder_plt.dat",status="replace")
        write (12,'("# wavenumber, xx, yy, zz, xy, yz, xz")')
        if (loren) then
         write (12,'("# Lorentzian function with hw =",f8.1," cm-1")')hw
        else
         write (12,'("# Gaussian   function with hw =",f8.1," cm-1")')hw
        end if
        do indf = 1, 5000
          write (12,'(i5,6(x,f20.10))') indf,(ramout(j,indf),j=1,6)
        end do
        close (12)
C
        deallocate (freq,raman,ramout)
        write (*,*) "alpder_proj.dat IS GENERATED"
        if (square) then
        write (*,*) "SQUARED VALUES ARE PLOTTED FOR XX, YY, ZZ, XY, YZ,
     * XZ"
        else
        write (*,*) "ABSOLUTE VALUES ARE PLOTTED FOR XX, YY, ZZ, XY, YZ,
     * XZ"
        END IF
      end if
C
      write (*,*)
      write (*,*) "  ------------------"
      write (*,*) "  NORMAL TERMINATION"
      write (*,*) "  ------------------"
      write (*,*)
      if (dipole.and.alpha) then
      write (*,*) "  FOLLOWING FILES ARE GENERATED"
      else
      write (*,*) "  FOLLOWING FILE IS GENERATED"
      end if
      if (dipole) then
      write (*,*) "  8. dipder_plt.dat"
      end if
      if (alpha) then
      write (*,*) "  9. alpder_plt.dat"
      end if
      write (*,*)
C
      stop
C
      end subroutine plot_dat
