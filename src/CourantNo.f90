subroutine CourantNo
!
! Purpose: 
!   Calculate and output the mean and maximum Courant Numbers.
!
  use types
  use parameters, only: CoNum,meanCoNum, CoNumFixValue, CoNumFix, timestep, ltransient, itime, time
  use geometry, only: numCells, numInnerFaces, owner, neighbour, numBoundaries, bctype, nfaces, startFace, Vol
  use sparse_matrix, only: res
  use variables, only: flmass,den

  implicit none

  integer :: i, ijp, ijn, inp, ib, iface
  real(dp):: TotalVol,dt

  if (ltransient) then
  CoNum = 0.0_dp
  meanCoNum = 0.0_dp

  res = 0.0_dp

  !
  ! Suface sum of magnitude (i.e. absolute value) of mass flux phi, over inner faces only
  !

  ! Inner faces                                               
  do i=1,numInnerFaces                                                                                                               
    ijp = owner(i)
    ijn = neighbour(i)
    res(ijp) = res(ijp)+abs(flmass(i))
    res(ijn) = res(ijn)+abs(flmass(i))                                                                                                               
  enddo     

  ! Boundary faces
  do ib=1,numBoundaries

    if ( bctype(ib) == 'inlet' ) then
      ! Inlet boundaries (mass fluxes prescribed in routine 'bcin')
      do i=1,nfaces(ib)
        iface = startFace(ib) + i
        ijp = owner(iface)
        res(ijp)=res(ijp)+ abs(flmass(iface))
      enddo

    elseif ( bctype(ib) == 'outlet' ) then
      do i=1,nfaces(ib)
        iface = startFace(ib) + i
        ijp = owner(iface)
        res(ijp)=res(ijp)+abs(flmass(iface))
      enddo
    endif

  enddo


  ! Accumulate by looping trough cells
  TotalVol = 0.0_dp
 
  do inp=1,numCells

    CoNum = max( CoNum , res(inp)/(den(inp)*Vol(inp)) )

    meanCoNum = meanCoNum + res(inp)/den(inp)
    
    TotalVol = TotalVol + Vol(inp)

  enddo

  res = 0.0_dp

  CoNum = 0.5*CoNum*timestep
  meanCoNum = 0.5*meanCoNum/TotalVol*timestep

  !// If we keep the value of Courant Number fixed
  if( CoNumFix ) then
      dt = timestep
      timestep = CoNumFixValue/CoNum * timestep

      CoNum = CoNum * timestep/dt
      meanCoNum  = meanCoNum * timestep/dt
  endif

  time = time + timestep

  write(6,'(a)') ' '
  write(6,'(a,i7,2(a,es10.3))') "* Timestep No.",itime," Dt: ",timestep," Time = ",time
  write(6,'(2(a,es10.3))') "* Courant Number Average: ", meanCoNum," Max  = ", CoNum
  write(6,'(a)') ' '


endif

end subroutine
