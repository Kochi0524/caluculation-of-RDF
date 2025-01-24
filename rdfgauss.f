      program rdfgauss
      implicit none
c----- Variable decreation 
      integer natm,natmo, mx, my, mz, l
      integer ia, ja, k, elemIndex, ix, iy, iz, i
      integer, allocatable::zz(:), elemType(:), elemCounts(:)
      real*8 hunit, dd, omega, r_max
      real*8 rho, dx, dy, dz, dis, sigma, gauss, dt, w, nfactor
      real*8 hm(3, 3)
      real*8, allocatable:: ta(:, :), sa(:, :)
      real*8, allocatable :: uco(:,:), co(:,:), rdf(:)
      real*8, parameter ::pi =  3.141592d0 
      

c-----file read
      read(*, *)
      read(*, *)
      read(*, *) hm(1, 1), hm(2, 1), hm(3, 1)
      read(*, *) hm(1, 2), hm(2, 2), hm(3, 2)
      read(*, *) hm(1, 3), hm(2, 3), hm(3, 3)
      read(*, *)
      read(*, *) natmo

      allocate(sa(3, natmo), zz(natmo),
     &         ta(3, natmo), elemType(natmo), co(3,natmo))
      do ia =  1, natmo
        read(*, *) zz(ia), sa(1, ia), sa(2, ia), sa(3, ia)
      end do

c-----count. the number of type of atom
      elemType(1) = 1

      do ia =  2, natmo
        if (zz(ia).eq.zz(ia - 1) ) then  
          elemType(ia) = elemType(ia - 1)  
        else if (zz(ia).ne.zz(ia - 1)) then   
          elemType(ia) = elemType(ia - 1) + 1  
        end if  
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
      end do

c-----define co
      do ia =  1, natmo
        do i =  1, 3
          co(i, ia) = sa(i, ia)
        end do
      end do


c-----setting pairs
      r_max = 10.0d0
      mx = aint(r_max / hm(1,1)) + 1
      my = aint(r_max / hm(2,2)) + 1
      mz = aint(r_max / hm(3,3)) + 1

      natm = natmo*((2*mx + 1)*(2*my + 1)*(2*mz + 1))
      allocate(uco(3,natm))

      uco(:,:) = 0.0d0
      ja = 0

      do ix =  -mx, mx
      do iy =  -my, my
      do iz =  -mz, mz
c    ia max is natm or natmo ?
        do ia =  1, natmo
          ja = ja + 1
          uco(1,ja) = co(1,ia)+hm(1,1)*dble(ix)
     &                        +hm(1,2)*dble(iy)
     &                        +hm(1,3)*dble(iz)
          uco(2,ja) = co(1,ia)+hm(2,1)*dble(ix)
     &                        +hm(2,2)*dble(iy)
     &                        +hm(2,3)*dble(iz)
          uco(3,ja) = co(1,ia)+hm(3,1)*dble(ix)
     &                        +hm(3,2)*dble(iy)
     &                        +hm(3,3)*dble(iz)
        end do
      end do        
      end do
      end do


c------calc.distance from own atom to pair atom
      omega = hm(1,1)*hm(2,2)*hm(3,3)+hm(1,2)*hm(3,2)*hm(1,3)
     &          +hm(1,3)*hm(2,1)*hm(3,2)-hm(1,1)*hm(2,3)*hm(3,2)
     &          -hm(1,2)*hm(2,1)*hm(3,3)-hm(1,3)*hm(2,2)*hm(3,1)  
      rho = dble(natmo) / omega
c     l, sigma, w:temp
      l = 100
      sigma = 0.5d0
      w = r_max / dble(l)
      allocate(rdf(l))

      do ia =  1, natmo
        do ja =  1, natm
          dx = uco(1,ja) - co(1,ia)
          dy = uco(2,ja) - co(2,ia)
          dz = uco(3,ja) - co(3,ia)
          dis = sqrt(dx**2 + dy**2 + dz**2)

          if (dis.lt.1.0d-8.or.dis.gt.l*1.5d0) cycle
          nfactor = 4.0d0*pi*rho*dble(natmo)*(dis)**2
c         w:mesh width, (l-1)*w = r_max
          do k =  1, l
            dt = (k-1)*w
            gauss = dexp(-0.5d0* ((dt-dis)/sigma)**2) 
     &        / (sigma * dsqrt(2.0d0*pi))
            rdf(k) = rdf(k) + (gauss / nfactor)

c            if (ja.eq.natm.and.ia.eq.2) then
c              write(*, *) "gauss",gauss / nfactor, "dt", dt
c            end if
          end do
        end do
      end do


c-----TODO:file write
      open(11, file='rdf.dat', status='replace')
      do k =  1, l
        write(11, *) w*(k-1), rdf(k)
c        write(*, *) w*(k-1), rdf(k)
      end do
      close(11)

      end program rdfgauss