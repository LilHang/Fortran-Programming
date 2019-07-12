program grid_formatting
    implicit none
    integer:: num_node, num_ele, num_edge, ij, numbp, nobs
    integer:: i, j, k, l, num, kode, region_num, boundary_num, boundary_label, mat, lay, nmat, nlay, nin, size_geometry = 0
    integer, allocatable:: boundary_node(:), each_boundary_num(:), boundary_orientation(:), each_boundary_node(:), temp(:),&
                           geometry_information(:)
    logical:: is_initial_condition_q, prescribed_rd_rate, water_uptake, spatial_variability
    real:: par(10,10), x, z, h, conc, qn, beta, axz, bxz, dxz, angle, aniz1, aniz2
    real, allocatable:: soil_profile(:), initial_profile(:), initial_condition_water(:), initial_condition_c(:), width(:),&
                        coordinate(:,:), width_information(:)
    real, external:: fh, distance

    open(40, file='swms_2d.in\geometry information.txt', status='unknown')
    open(50, file='swms_2d.in\width information.txt', status='unknown')
    write(*,*) "open xxxx.msh gird file......"
    open(unit=10, file='swms_2d.in\nodal_drain.msh') ! 打开网格文件
    read(10, *) num_node, num_ele, num_edge ! 读入节点数、单元数，边数
    allocate(coordinate(2,num_node))

    write(*,"(/'touch a grid.in file......')")
    open(unit=20, file='swms_2d.in\grid.in', status='unknown')! 创建输入网格文件
    write(20,"('*** BLOCK H: NODAL INFORMATION **************************************************')")
    write(20,"('      NumNP     NumEl       IJ      NumBP     NObs')") ! 写入表头
    write(*,"(/'input parameter value of IJ, NumBP, Nobs respectively: ')")
    read(*,*) ij, numbp, nobs
    write(20,"(i10, i10, i10, i10, i8)") num_node, num_ele, ij, numbp, nobs ! 写入对应参数

    write(20,*) "   n  Code            x                 z                 h            Conc             &
    &Q           M        B         Axz       Bxz       Dxz" ! 写入表头
    write(*,"(/'is initial codintion soil water content ?'/&
                &'whether condsidering prescribed recharge/discharge rate, water uptake,&
                &spatial variability of hydraulic parameters ?'/&
                &'(t: yes, f: no, enter respectively)')")
    read(*,*) is_initial_condition_q, prescribed_rd_rate, water_uptake, spatial_variability

    call read_material_information(nmat, nlay, par)
    allocate(soil_profile(nmat + 1))
    write(*, "(/'input depth of soil interface (must be a decreasing sequence to zero !): ')")
    read(*,*) (soil_profile(i), i = 1, nmat + 1)

    write(*, "(/'input interface number of initial condition interface: ')")
    read(*,*) nin
    allocate(initial_condition_water(nin))
    allocate(initial_condition_c(nin))
    allocate(initial_profile(nin + 1))
    write(*, "('input depth of initial condition interface (must be a decreasing sequence to zero !): ')")
    read(*,*) (initial_profile(i), i = 1, nin + 1)
    write(*, "('input  initial condition (water) of interface: ')")
    read(*,*) (initial_condition_water(i), i = 1, nin)
    write(*, "('input  initial condition (concentration) of interface: ')")
    read(*,*) (initial_condition_c(i), i = 1, nin)

    do num = 1, num_node
        read(10, *) x, z, kode
        coordinate(1,num) = x
        coordinate(2,num) = z
        if(kode == 7 .or. kode == 8 .or. kode == 9) kode = 0
        i = 1
        do while (i <= nmat)
            if (z <= soil_profile(i) .and. z >= soil_profile(i + 1)) then
                mat = i
                exit
            else
                i = i + 1
            end if
        end do
        j = 1
        do while (j <= nin)
            if (z <= initial_profile(j) .and. z >= initial_profile(j + 1)) then
                h = initial_condition_water(j)
                conc = initial_condition_c(j)
                if (is_initial_condition_q) h= fh(h, par(1, mat))
                exit
            else
                j = j + 1
            end if
        end do

        if (prescribed_rd_rate) then
            call prescribed_rd()
        else
            qn = 0
        end if

        if (water_uptake) then
            call prescribed_rd()
        else
            beta = 0
        end if

        if (spatial_variability) then
            call variability()
        else
            axz = 1
            bxz = 1
            dxz = 1
        end if
        write(20, "(2i5,2f20.11,f15.4,2e15.2e3,i7,f12.4,3f10.2)") num, kode, x, z, h, conc, qn, mat, beta, axz, bxz, dxz
    end do
    write(*,"('nodal information done !')")

    write(20,"('*** BLOCK I: ELEMENT INFORMATION ************************************************')")
    write(20,"('   e    i    j    k    l   Angle  Aniz1    Aniz2    LayNum')")
    angle = 0
    aniz1 = 1
    aniz2 = 1
    do  num = 1, num_ele
        read(10,*) i, j, k, region_num
        if (region_num == 0) region_num = 1
        if (nlay == 1) then
            lay = nlay
        else
            lay = region_num
        end if
        write(20,"(i4,4i5,3f8.2,i7)") num, i, j, k, k, angle, aniz1, aniz2, lay
    end do
    write(*,"(/'element information done !')")

    write(20, "('*** BLOCK J: BOUNDARY GEOMETRY INFORMATION **************************************')")
    write(20, "('Node number array:')")
    write(*,"(/'input number of boundary type(0 not included): ')")
    read(*,*) boundary_num
    allocate (each_boundary_num(boundary_num))
    allocate(boundary_orientation(boundary_num))
    write(*,"('input number of nodes for each boundary type (edge * 2): ')")
    read(*,*) (each_boundary_num(i), i = 1, boundary_num)
    allocate(boundary_node(sum(each_boundary_num)))
    write(*,"('input boundary orientation for each boundary type(1: left, 2: right, 3: circle): ')")
    read(*,*) (boundary_orientation(i), i = 1, boundary_num)

    k = 1
    do num = 1, num_edge
        read(10,*) i, j, boundary_label
        if (boundary_label == 7 .or. boundary_label == 8 .or. boundary_label == 9) boundary_label = 0
        if (boundary_label /= 0) then
            boundary_node(k) = i
            boundary_node(k + 1) = j
            k = k + 2
        end if
    end do
    do k = 1, boundary_num
        allocate(temp(each_boundary_num(k)))
        if (k == 1) then
            temp = boundary_node(1:each_boundary_num(k))
        else
            num = sum(each_boundary_num(1:k-1))
            temp = boundary_node(num + 1 : num + each_boundary_num(k))
        end if

        if (boundary_orientation(k) == 1) then
            allocate(each_boundary_node(each_boundary_num(k) / 2 + 1))
            each_boundary_node(1) = temp(1)
            each_boundary_node(2) = temp(2)
            if (each_boundary_num(k) > 2) then
                i = 3
                do j = 4, each_boundary_num(k), 2
                    each_boundary_node(i) = temp(j)
                    i = i + 1
                end do
            end if
        else if (boundary_orientation(k) == 2) then
            allocate(each_boundary_node(each_boundary_num(k) / 2 + 1))
            each_boundary_node(1) = temp(each_boundary_num(k) - 1)
            each_boundary_node(2) = temp(each_boundary_num(k))
            if (each_boundary_num(k) > 2) then
                i = 3
                do j = each_boundary_num(k) - 2, 2, -2
                    each_boundary_node(i) = temp(j)
                    i = i + 1
                end do
            end if
        else if (boundary_orientation(k) == 3) then
            allocate(each_boundary_node(each_boundary_num(k) / 2 + 1))
            each_boundary_node(1) = temp(1)
            each_boundary_node(2) = temp(2)
            if (each_boundary_num(k) > 2) then
                i = 3
                do j = 4, each_boundary_num(k), 2
                    each_boundary_node(i) = temp(j)
                    i = i + 1
                end do
            end if
        else
            continue
        end if

        deallocate(temp)
        num = size(each_boundary_node)
        allocate(width(num))
        if (boundary_orientation(k) == 3) then
            width(1) = 0.5 * distance(coordinate(1,each_boundary_node(1)), coordinate(1, each_boundary_node(2))) + &
                       0.5 * distance(coordinate(1,each_boundary_node(1)),coordinate(1, each_boundary_node(num - 1)))
            do i = 2, num - 1
                width(i) = 0.5 * distance(coordinate(1, each_boundary_node(i)), coordinate(1, each_boundary_node(i - 1))) +&
                           0.5 * distance(coordinate(1, each_boundary_node(i)), coordinate(1, each_boundary_node(i + 1)))
            end do
        else
             width(1) = 0.5 * distance(coordinate(1,each_boundary_node(1)), coordinate(1, each_boundary_node(2)))
             width(num) =  0.5 * distance(coordinate(1,each_boundary_node(num)), coordinate(1, each_boundary_node(num - 1)))
            do i = 2, num - 1
                width(i) = 0.5 * distance(coordinate(1, each_boundary_node(i)), coordinate(1, each_boundary_node(i - 1))) +&
                           0.5 * distance(coordinate(1, each_boundary_node(i)), coordinate(1, each_boundary_node(i + 1)))
            end do
        end if
        if (boundary_orientation(k) == 3) num = num - 1
        size_geometry = size_geometry + num
        write(40, "(i4)") (each_boundary_node(l), l=1, num)
        write(50,"(f8.4)") (width(l), l=1, num)
        deallocate(each_boundary_node)
        deallocate(width)
    end do

    rewind(40)
    rewind(50)
    allocate(geometry_information(size_geometry))
    allocate(width_information(size_geometry))
    do i = 1, size_geometry
        read(40,*) geometry_information(i)
        read(50,*) width_information(i)
    end do

    write(20, "(10i6)") geometry_information
    write(20,"('width array')")
    write(20,"(10f12.4)") width_information
    write(20, "('length: ')")
    write(20, "('nobs: ')")
    write(20, "('*** END OF INPUT FILE GRID.IN *************************************************')")
    close(20)
    close(10)
    close(40)
    close(50)

    write(*, "('FINISH FORMATTING !')")
    read(*,*)



end program grid_formatting

subroutine read_material_information(nmat, nlay, par)

    implicit none
    integer:: i, m
    integer:: nmat, nlay, npar
    real:: par(10, 10), h1, hn

    open(30, file='SWMS_2D.IN\Selector.in', status='old')
    do i = 1,13
        read(30,*)
    end do

    read(30, *) nmat, nlay, h1, hn, npar
    read(30,*)
    do m = 1, nmat
        read(30, *) (par(i,m), i = 1, npar)
    end do
    close(30)
    return

end subroutine read_material_information

subroutine prescribed_rd()
    implicit none
    continue
end subroutine prescribed_rd

subroutine uptake()
    implicit none
    continue
end subroutine uptake

subroutine variability()
    implicit none
    continue
end subroutine variability

Real Function fh(qe, par)

Implicit Double Precision (A-H, O-Z)
Double Precision n, m
Real qe, par(9)

qr = par(1)
qs = par(2)
qa = par(3)
qm = par(4)
alfa = par(5)
n = par(6)

m = 1.D0 - 1.D0/n
hmin = -1.D300**(1.D0/n)/max(alfa, 1.D0) ! 确定一个最小负压即最大土壤吸力
qeem = (1.D0+(-alfa*hmin)**n)**(-m)      ! 根据hmin确定一个最小有效饱和度
qee = dmin1(dmax1(qe*(qs-qa)/(qm-qa),qeem), .999999999999999D0) ! min为了限制有效饱和度小于“1”，max为了限制有效饱和度大于qeem(最小有效饱和度)
fh = sngl(max(-1.D0/alfa*(qee**(-1.D0/m)-1.D0)**(1.D0/n),-1.D37)) ! 根据vg模型计算h并转换为单精度，max是为了限制土壤吸力小于1d37
Return

End Function fh

real function distance(a, b)
implicit none
    real:: a(2), b(2), x, y
    x = (a(1) - b(1)) ** 2
    y = (a(2) - b(2)) ** 2
    distance = sqrt(x + y)
    return
end function distance
