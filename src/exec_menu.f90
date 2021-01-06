subroutine Exec_Menu(Iter)

    use,intrinsic::iso_c_binding
    use ogpf
    use variables,only:plt,plotarray
    use utils, only: formatSTR
    use parameters,only: lcal,nphi,resor,sor,time,itime,timestep,CoNum,meanCoNum
    
    integer,intent(in)::Iter
    
    interface
    function KBHIT()bind(C,name='_kbhit')
    integer::KBHIT
    end function
    function GETCH()bind(C,name='getche')
    integer::GETCH
    end function
    end interface
    
    integer::kbstat = 0
    character(1)::cmdstr
    integer::ch,i,maxsize
    character(len=:), allocatable :: varnames
    character(len=10) :: len10str
    character(len=7), dimension(nphi) ::  cs = &
  (/'U      ','V      ','W      ','Mass   ','TE      ','ED','T  ','VIS   ','VART','CON   ','EPOT   ' /)
    
    type(gpf)::gp
    real,dimension(:),allocatable::x
    
    kbstat=KBHIT()
    !sense keyboard input.
    if(kbstat.NE.0)then
        ch=getch()
        if(ch.eq.iachar('m'))then
            write(*,*)''
            WRITE(*,'(2x,a)') '------------------------------------------------------------------------------'
            write(*,'(2x,a)') 'Menu Invoked By Keyboard, Available Options Are 0~1:'
            ch=getch()
            select case(ch)
            case(ichar('0'))
                write(*,'(1x,a,a,$)')'Plot Residual By Gnuplot.'
                maxsize=floor(iter/float(2**plt%Level))
                allocate(x(maxsize))
                do i=1,maxsize
                    x(i)=(i-1)*2**plt%Level+1.0
                end do
                call gp%title('Normalization Residual')
                call gp%xlabel('Iters')
                call gp%ylabel('Residual')
                call gp%options('set key top right')
                call gp%options('set autoscale fix')
                call gp%options('set format y "%0.1e"')
                call gp%options('set grid lc "black" lw 2')
                call gp%options('set border linewidth 2')
                
                call gp%semilogy(x,plt%arraytmp(1:maxsize,:),lspec='with linespoints lt 1 title "U_x";&
            & with lp lt 2 lw 1 ps 1 title "U_y";&
            & with lp lt 3 lw 1 ps 1 title "U_z";&
            & with lp lt 4 lw 1 ps 1 title "Mass"')
                
                deallocate(x)
               
            case (ichar('1'))
                write(*,'(1x,a,a,$)')'Print variable informations.'
                write(*,'(2x,a,i7,2(a,es10.3))') "Timestep No.",itime," Dt: ",timestep," Time = ",time
                write(*,'(2x,a,es10.3,a,es10.3)') "Courant Number Average: ", meanCoNum," Max  = ", CoNum
                write(*,'(2x,a)')'______________________________________________________________________________'
                varnames = ''
                do i =1,NPHI
                    if(LCAL(I))then
                        write(len10str,'(a1,a9)')'|',trim(adjustl(cs(i)))
                        varnames = varnames//len10str
                    end if
                end do                
                write(*,'(2x,a,a)')'Var.Name',varnames
                varnames =''
                do i=1,NPHI
                    if(LCAL(i))then
                        write(len10str,'(a1,E9.2)')'|',sor(i)
                        varnames = varnames//len10str
                    end if
                end do
                write(*,'(2x,a,a)')'Tar.Res.',varnames
                varnames =''
                do i=1,NPHI
                    if(LCAL(i))then
                        write(len10str,'(a1,E9.2)')'|',resor(i)
                        varnames = varnames//len10str
                    end if
                end do
                write(*,'(2x,a,a)')'Cur.Res.',varnames       
                
            case default
                write(*,'(1x,a,a,$)')'Unexpected Input.'
            end select
                
            WRITE(*,'(2x,a)') '------------------------------------------------------------------------------'
        end if
        kbstat=0
    endif
    
end subroutine Exec_Menu