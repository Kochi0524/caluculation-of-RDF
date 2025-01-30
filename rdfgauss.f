      program rdfgauss
      implicit none
c----- Variable decreation 
      integer natm,natmo, mx, my, mz, l,own,pair
      integer ia, ja, k, elemIndex, ix, iy, iz, i, ii
      integer, allocatable::zz(:), elemType(:), elemCounts(:)
      integer zztemp(10)
      real*8 hunit, dd, omega, r_max
      real*8 rho, dx, dy, dz, dis, sigma, gauss, dt, w, nfactor
      real*8 hm(3, 3)
      real*8, allocatable:: sa(:, :), ucoType(:)
      real*8, allocatable :: uco(:,:), co(:,:), rdf(:)
      real*8, parameter ::pi =  3.141592d0  
      character(len=30) inputfile,outputfile, selectFile

c-----file read
      inputfile = selectFile()
      open(13, file= inputfile, status='old')
      read(13, *)
      read(13, *)
      read(13, *) hm(1, 1), hm(2, 1), hm(3, 1)
      read(13, *) hm(1, 2), hm(2, 2), hm(3, 2)
      read(13, *) hm(1, 3), hm(2, 3), hm(3, 3)
      read(13, *)
      read(13, *) natmo

      allocate(sa(3, natmo), zz(natmo),
     &         elemType(natmo), co(3,natmo))
      do ia =  1, natmo
        read(13, *) zz(ia), sa(1, ia), sa(2, ia), sa(3, ia)
      end do
      close(13)

c-----count. the number of type of atom
      elemType(1) = 1
      zzTemp(:) = 0
      zzTemp(1) = zz(1)
      
      do ia =  2, natmo
        do ii =  1, 10
          if (zz(ia).eq.zzTemp(ii)) then
            elemType(ia) = ii
          exit
          else if (zzTemp(ii).eq.0) then
            zzTemp(ii) = zz(ia)
            elemType(ia) = ii
            exit
          end if
        end do
      end do

c-----count. the number of each atoms
      allocate(elemCounts(elemType(natmo)))
      elemCounts(:) = 0

      do elemIndex =  1, elemType(natmo)
        do ia =  1, natmo
          if (elemType(ia).eq.elemIndex) then
            elemCounts(elemIndex) = elemCounts(elemIndex) + 1
          end if
        end do
      end

c-----select own atom and pair atom number
      write(6, *) "Input own atom , pair atom"
      write(6, *) "Enter a less than or equal to ", elemType(natmo)
      read(5, *) own, pair


c-----define co
      do ia =  1, natmo
        do i =  1, 3
          co(i, ia) = sa(i, ia)
        end do
      end do

c-----setting coordinates of pairs
      r_max = 10.0d0
      mx = aint(r_max / hm(1,1)) + 1
      my = aint(r_max / hm(2,2)) + 1
      mz = aint(r_max / hm(3,3)) + 1

      natm = natmo*((2*mx + 1)*(2*my + 1)*(2*mz + 1))
      allocate(uco(3,natm), ucoType(natm))

      uco(:,:) = 0.0d0
      ucoType(:) = 0
      ja = 0

      do ix =  -mx, mx
      do iy =  -my, my
      do iz =  -mz, mz
        do ia =  1, natmo
          ja = ja + 1
          ucoType(ja) = elemType(ia)
          uco(1,ja) = co(1, ia)+hm(1,1)*dble(ix)
     &                         +hm(1,2)*dble(iy)
     &                         +hm(1,3)*dble(iz)
          uco(2,ja) = co(2, ia)+hm(2,1)*dble(ix)
     &                         +hm(2,2)*dble(iy)
     &                         +hm(2,3)*dble(iz)
          uco(3,ja) = co(3 ,ia)+hm(3,1)*dble(ix)
     &                         +hm(3,2)*dble(iy)
     &                         +hm(3,3)*dble(iz)
        end do
      end do        
      end do
      end do


c------calc.distance from own atom to pair atom
      omega = hm(1,1)*hm(2,2)*hm(3,3)+hm(1,2)*hm(3,2)*hm(1,3)
     &          +hm(1,3)*hm(2,1)*hm(3,2)-hm(1,1)*hm(2,3)*hm(3,2)
     &          -hm(1,2)*hm(2,1)*hm(3,3)-hm(1,3)*hm(2,2)*hm(3,1)  
      rho = dble(elemCounts(own)) / omega

      l = 100
      sigma = 0.2d0
      w = r_max / dble(l)
      allocate(rdf(l))

      do ia =  1, natmo
        if (elemType(ia).ne.own) cycle
        do ja =  1, natm
          if (ucoType(ja).ne.pair) cycle
          dx = uco(1,ja) - co(1,ia)
          dy = uco(2,ja) - co(2,ia)
          dz = uco(3,ja) - co(3,ia)
          dis = sqrt(dx**2 + dy**2 + dz**2)

          if (dis.lt.1.0d-8.or.dis.gt.l*1.5d0) cycle
          nfactor = 4.0d0*pi*rho*dble(elemCounts(pair))*(dis)**2
c         w:mesh width, (l-1)*w = r_max
          do k =  1, l
            dt = (k-1)*w
            gauss = dexp(-0.5d0* ((dt-dis)/sigma)**2) 
     &        / (sigma * dsqrt(2.0d0*pi))
            rdf(k) = rdf(k) + (gauss / nfactor)
          end do
        end do
      end do


c-----file write
      outputfile = selectFile()
      open(11, file=filePath, status='replace')
      do k =  1, l
        write(11, *) w*(k-1), rdf(k)
      end do
      close(11)

      end program rdfgauss

c-----choose filepath 
c-----data/ディレクトリに入力ファイルが存在するものとする
      function selectFile()
        implicit none
        character(len=30) filename, selectFile
        write(*, *) "Input file name"
        read(*, *) filename

        selectFile = "data/"//filename
        return 
      end function selectFile