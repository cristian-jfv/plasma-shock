submodule (scheme) foamm_scheme
    use model, only: gamma, calc_primitives, calc_cons_quantities, &
                    alfven_speed, fast_speed, slow_speed, sound_speed, f
    use scheme, only: print_matrix
    implicit none

    contains



        module function flux(u, dx, dt) result(ans)
            real, dimension(:,:), intent(in) :: u
            real, intent(in):: dx, dt
            real, target, dimension(size(u,1), size(u,2)-1) :: avg_p
            real, dimension(size(u,1), size(u,2)-1) :: ans, p_star, u_star
            real, dimension(size(u,1)) :: dummy_p_star
            real, dimension(size(u,1), size(u,2)) :: p, f_complete
            real, dimension(size(u,2)-1) :: a, cf, cs, va
            real, dimension(7, size(u,2)-1) :: eigenvalues
            real, dimension(7,8) :: l_eigen, r_eigen
            real, dimension(8,8) :: a_mat
            integer :: p_cols, r, c
            real, pointer, dimension(:) :: pg, rho, bx
            real, pointer, dimension(:,:) :: b
            p_cols = size(p,2)

            ! Obtain the primitive variables from the conserved quantities
            call calc_primitives(u=u, p=p)
            call print_matrix(p, "p")
            ! Average from the neighbouring cells
            avg_p = 0.5*(p(:,1:p_cols-1) + p(:,2:p_cols))

            ! Calculate eigenvalues
            !        1   2   3   4   5   6   7   8
            ! P = [rho, vx, vy, vz, Bx, By, Bz, pg]
            eigenvalues = calc_eigenvalues(p=avg_p)
            call print_matrix(avg_p, "avg_p")
            call print_matrix(eigenvalues, "eigenvalues")
            rho => avg_p(1,:)
            bx => avg_p(5,:)
            pg => avg_p(8,:)
            b => avg_p(5:7,:)
            
            a = sound_speed(pg, rho)
            va = alfven_speed(bx, rho)
            cs = slow_speed(a, va, b=b, rho=rho, whole_b=1)
            cf = fast_speed(a, va, b=b, rho=rho, whole_b=1)
            
            ! Calculate p_star
            do c = 1, p_cols-1
                !print *, "c=", c, "***************************************************************************"

                ! Proceed colwise
                ! Calculate the left and right eigenvectors
                l_eigen = calc_left_eigenvectors(p=avg_p(:,c), a=a(c), cf=cf(c), cs=cs(c))
                !call print_matrix(l_eigen, "left eigenvectors")
                r_eigen = calc_right_eigenvectors(p=avg_p(:,c), a=a(c), cf=cf(c), cs=cs(c))
                !call print_matrix(r_eigen, "right eigenvectors")

                !print *, "pr - pl"
                !print *, p(:,c+1) - p(:,c)

                ! use p_star = p_l + sum(over negative eigenvalues)
                p_star(:,c) = 0.0
                do r = 1, 7
                    if (isnan(r_eigen(r,1)) .or. isnan(l_eigen(r,1))) then
                        !print *, "EIGENVECTOR IS NAN"
                        a_mat = calc_a_mat(avg_p(:,c), a=a(c))
                        !print *, "avg_p(:,c)"
                        !print *, avg_p(:,c)
                        !call print_matrix(a_mat, "A")
                        !print *, "pg=", pg(c)
                    else
                        !print *, " "
                        !print *, "r=", r
                        a_mat = calc_a_mat(avg_p(:,c), a=a(c))
                        !call print_matrix(a_mat, "A")
                        !print *, "check right eigenvector"
                        !call check_r_eigenvector(eigenvalues(r,c), a_mat, r_eigen(r,:))
                        !print *, "check left eigenvector"
                        !call check_l_eigenvector(eigenvalues(r,c), a_mat, l_eigen(r,:))

                        p_star(:,c) = p_star(:,c) + abs(eigenvalues(r,c))*r_eigen(r,:)*dot_product(l_eigen(r,:), p(:,c+1) - p(:,c))
                    end if
                end do
                !call print_matrix(calc_dudp(avg_p(:,c)), "du/dp")
                p_star(:,c) = matmul(calc_dudp(avg_p(:,c)), p_star(:,c))
                !print *, "p_star"
                !print *, p_star(:,c)

            end do

            call print_matrix(p_star, "p_star")

            call f(u, f_complete)
            ! call calc_cons_quantities(p_star, u_star)
            ans = 0.5*( f_complete(:,1:p_cols-1) + f_complete(:,2:p_cols) - p_star)

            !call f(u_star, ans)
        end function flux

        subroutine check_r_eigenvector(lambda, a_mat, r_eigen)
            real, dimension(8,8), intent(in) :: a_mat
            real, dimension(8), intent(in) :: r_eigen
            real, intent(in) :: lambda
            real, dimension(8) :: residual, dummy_a, dummy_b

            
            !print *, "eigenvetor: ", r_eigen
            
            dummy_a = lambda*r_eigen
            print *, "eigenvalue mult", dummy_a
            dummy_b = matmul(a_mat, r_eigen)
            print *, "matmul: ", dummy_b

            residual = dummy_a - dummy_b
            print *, "residual: ", residual

        end subroutine check_r_eigenvector

        subroutine check_l_eigenvector(lambda, a_mat, l_eigen)
            real, dimension(8,8), intent(in) :: a_mat
            real, dimension(8), intent(in) :: l_eigen
            real, intent(in) :: lambda
            real, dimension(8) :: residual, dummy_a, dummy_b

            !print *, "eigenvector: ", l_eigen
            
            dummy_a = lambda*l_eigen
            print *, "eigenvalue mult", dummy_a
            dummy_b = matmul(l_eigen, a_mat)
            print *, "matmul: ", dummy_b

            residual = dummy_a - dummy_b
            print *, "residual: ", residual

        end subroutine check_l_eigenvector


        function calc_a_mat(p, a) result(ans)
            real, target, dimension(8), intent(in) :: p
            real, intent(in) :: a
            real, dimension(8,8) :: ans
            real, pointer :: rho, vx, vy, vz, bx, by, bz, pg
            real :: v
            
            rho => p(1)
            vx => p(2)
            vy => p(3)
            vz => p(4)
            bx => p(5)
            by => p(6)
            bz => p(7)
            pg => p(8)
            
            v = sqrt(vx**2 + vy**2 + vz**2)
            
            ans(:,:) = 0.0
            !      p = [ rho,      vx,  vy,  bz,  bx,     by,     bz,    pg ]
            ans(1,:) = [  vx,      rho, 0.0, 0.0, 0.0,    0.0,    0.0,   0.0 ]
            ans(2,:) = [ 0.0,       vx, 0.0, 0.0, 0.0, by/rho, bz/rho, 1/rho ]
            ans(3,:) = [ 0.0,      0.0,  vx, 0.0, 0.0,-bx/rho,    0.0,   0.0 ]
            ans(4,:) = [ 0.0,      0.0, 0.0,  vx, 0.0,    0.0,-bx/rho,   0.0 ]
            ans(5,:) = [ 0.0,      0.0, 0.0, 0.0, 0.0,    0.0,    0.0,   0.0 ]
            ans(6,:) = [ 0.0,       by, -bx, 0.0, 0.0,     vx,    0.0,   0.0 ]
            ans(7,:) = [ 0.0,       bz, 0.0, -bx, 0.0,    0.0,     vx,   0.0 ]
            ans(8,:) = [ 0.0, rho*a**2, 0.0, 0.0, 0.0,    0.0,    0.0,    vx ]

        end function calc_a_mat

        function calc_dudp(p) result(ans)
            real, target, dimension(:), intent(in) :: p
            real, dimension(size(p,1), size(p,1)) :: ans
            real, pointer :: rho, vx, vy, vz, bx, by, bz, pg
            real :: v
            
            rho => p(1)
            vx => p(2)
            vy => p(3)
            vz => p(4)
            bx => p(5)
            by => p(6)
            bz => p(7)
            pg => p(8)
            
            v = sqrt(vx**2 + vy**2 + vz**2)

            ans(:,:) = 0.0
            
            ans(:,1) = (/ 1.0, vx, vy, vz, 0.0, 0.0, 0.0, 0.5*v**2 /)
            ans(2,2) = rho
            ans(3,3) = rho
            ans(4,4) = rho
            ans(5,5) = 1.0
            ans(6,6) = 1.0
            ans(7,7) = 1.0
            ans(8,:) = (/ 0.5*v**2, rho*vx, rho*vy, rho*vz, bx, by, bz, 1.0/(gamma - 1.0) /)
            
            
        end function calc_dudp

        function calc_left_eigenvectors(p, a, cf, cs) result(ans)
            real, target, dimension(8), intent(in) :: p
            real, intent(in) :: a, cf, cs
            real, dimension(7,8) :: ans
            real, pointer :: rho
            real, pointer, dimension(:) :: b
            rho => p(1)
            b => p(5:7)

            ! Each row has an eigenvector
            ans(1,:) = l1(rho, b, c=cf)
            ans(2,:) = l2(rho, b)
            ans(3,:) = l3(rho, b, c=cs)
            ans(4,:) = l4(a)
            ans(5,:) = l5(rho, b, c=cs)
            ans(6,:) = l6(rho, b)
            ans(7,:) = l7(rho, b, c=cf)

        end function calc_left_eigenvectors

        function calc_right_eigenvectors(p, a, cf, cs) result(ans)
            real, target, dimension(8), intent(in) :: p
            real, intent(in) :: a, cf, cs
            real, dimension(7,8) :: ans
            real, pointer, dimension(:) :: b
            real, pointer :: rho
            rho => p(1)
            b => p(5:7)
            
            ans(1,:) = r1(rho, b, c=cf, a=a)
            ans(2,:) = r2(rho, b)
            ans(3,:) = r3(rho, b, c=cs, a=a)
            ans(4,:) = r4()
            ans(5,:) = r5(rho, b, c=cs, a=a)
            ans(6,:) = r6(rho, b)
            ans(7,:) = r7(rho, b, c=cf, a=a)

        end function calc_right_eigenvectors

        function l1(rho, b, c) result(ans)
            real, target, dimension(3), intent(in) :: b
            real, intent(in) :: rho, c
            real, dimension(8) :: ans
            real, pointer :: bx, by, bz
            real :: delta
            ! P = [rho, vx, vy, vz, Bx, By, Bz, pg]

            bx => b(1)
            by => b(2)
            bz => b(3)
            delta = rho*c**2 - bx**2

            ans = (/ 0.0, -c, c*bx*by/delta, c*bx*bz/delta, 0.0, c**2*by/delta, c**2*bz/delta, 1.0/rho/)
            ans = ans/sqrt(dot_product(ans, ans))
            
        end function l1

        function l7(rho, b, c) result(ans)
            real, target, dimension(3), intent(in) :: b
            real, intent(in) :: rho, c
            real, dimension(8) :: ans
            real, pointer :: bx, by, bz
            real :: delta
            ! P = [rho, vx, vy, vz, Bx, By, Bz, pg]

            bx => b(1)
            by => b(2)
            bz => b(3)
            delta = rho*c**2 - bx**2

            ans = (/ 0.0, c, -c*bx*by/delta, -c*bx*bz/delta, 0.0, c**2*by/delta, c**2*bz/delta, 1.0/rho/)
            ans = ans/sqrt(dot_product(ans, ans))
        end function l7

        function csign(input) result(output)
            real, intent(in) :: input
            real :: output

            if (input > 0) then
                output = 1.0
            else if (input < 0) then
                output = -1.0
            else
                output = 0.0
            end if

        end function csign

        function l2(rho, b) result(ans)
            real, target, dimension(3), intent(in) :: b
            real, intent(in) :: rho
            real, dimension(8) :: ans
            real, pointer :: bx, by, bz
            real :: s
            
            bx => b(1)
            by => b(2)
            bz => b(3)
            s = csign(bx)

            ans = (/ 0.0, 0.0, s*bz, -s*by, 0.0, bz/sqrt(rho), -by/sqrt(rho), 0.0 /)
            
            if (dot_product(ans, ans) > 0) then
                ans = ans/sqrt(2.0*(by**2 + bz**2))
            end if

        end function l2

        function l6(rho, b) result(ans)
            real, target, dimension(3), intent(in) :: b
            real, intent(in) :: rho
            real, dimension(8) :: ans
            real, pointer :: bx, by, bz
            real :: s
            
            bx => b(1)
            by => b(2)
            bz => b(3)
            s = csign(bx)

            ans = (/ 0.0, 0.0, s*bz, -s*by, 0.0, -bz/sqrt(rho), by/sqrt(rho), 0.0 /)
            if (dot_product(ans, ans) > 0) then
                ans = ans/sqrt(2.0*(by**2 + bz**2))
            end if

        end function l6

        function l3(rho, b, c) result(ans)
            real, intent(in) :: rho, c
            real, target, dimension(3), intent(in) :: b
            real, dimension(8) :: ans
            real, pointer :: bx, by, bz
            real :: delta

            bx => b(1)
            by => b(2)
            bz => b(3)
            if (by**2 + bz**2 > 0.0) then
                delta = rho*c**2 - bx**2
            else
                delta = 1
            end if

            ans = (/ 0.0, -c, c*bx*by/delta, c*bx*bz/delta, 0.0, c**2*by/delta, c**2*bz/delta, 1.0/rho /)
            ans = ans/sqrt(dot_product(ans, ans))
        end function l3

        function l5(rho, b, c) result(ans)
            real, intent(in) :: rho, c
            real, target, dimension(3), intent(in) :: b
            real, dimension(8) :: ans
            real, pointer :: bx, by, bz
            real :: delta

            bx => b(1)
            by => b(2)
            bz => b(3)
            if (by**2 + bz**2 > 0.0) then
                delta = rho*c**2 - bx**2
            else
                delta = 1
            end if

            ans = (/ 0.0, c, -c*bx*by/delta, -c*bx*bz/delta, 0.0, c**2*by/delta, c**2*bz/delta, 1.0/rho /)
            ans = ans/sqrt(dot_product(ans, ans))
        end function l5

        function l4(a) result(ans)
            real, intent(in) :: a
            real, dimension(8) :: ans
            !        1  2  3  4  5  6  7  8
            ans = (/ 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0/a**2 /)
            ans = ans/sqrt(dot_product(ans, ans))
        end function l4

        function r1(rho, b, c, a) result(ans)
            real, intent(in) :: rho, c, a
            real, target, dimension(3), intent(in) :: b
            real, dimension(8) :: ans
            real, pointer :: bx, by, bz
            real :: delta

            bx => b(1)
            by => b(2)
            bz => b(3)
            delta = rho*c**2 - bx**2

            ans = (/ rho, -c, c*bx*by/delta, c*bx*bz/delta, 0.0, rho*c**2*by/delta, rho*c**2*bz/delta, rho*a**2 /)
            ans = ans/sqrt(dot_product(ans, ans))
        end function r1

        function r7(rho, b, c, a) result(ans)
            real, intent(in) :: rho, c, a
            real, target, dimension(3), intent(in) :: b
            real, dimension(8) :: ans
            real, pointer :: bx, by, bz
            real :: delta

            bx => b(1)
            by => b(2)
            bz => b(3)
            delta = rho*c**2 - bx**2

            ans = (/ rho, c, -c*bx*by/delta, -c*bx*bz/delta, 0.0, rho*c**2*by/delta, rho*c**2*bz/delta, rho*a**2 /)
            ans = ans/sqrt(dot_product(ans, ans))
        end function r7

        function r2(rho, b) result(ans)
            real, intent(in) :: rho
            real, target, dimension(3), intent(in) :: b
            real, dimension(8) :: ans
            real, pointer :: bx, by, bz
            real :: s

            bx => b(1)
            by => b(2)
            bz => b(3)

            s = csign(bx)
            
            ans = (/ 0.0, 0.0, s*bz, -s*by, 0.0, bz**sqrt(rho), -by*sqrt(rho), 0.0 /)
            if (dot_product(ans, ans) > 0) then
                ans = ans/sqrt(2.0*(by**2 + bz**2))
            end if

        end function r2

        function r6(rho, b) result(ans)
            real, intent(in) :: rho
            real, target, dimension(3), intent(in) :: b
            real, dimension(8) :: ans
            real, pointer :: bx, by, bz
            real :: s

            bx => b(1)
            by => b(2)
            bz => b(3)

            s = csign(bx)
            
            ans = (/ 0.0, 0.0, s*bz, -s*by, 0.0, -bz**sqrt(rho), by*sqrt(rho), 0.0 /)
            if (dot_product(ans, ans) > 0) then
                ans = ans/sqrt(2.0*(by**2 + bz**2))
            end if

        end function r6

        function r3(rho, b, c, a) result(ans)
            real, intent(in) :: rho, c, a
            real, target, dimension(3) :: b
            real, dimension(8) :: ans
            real, pointer :: bx, by, bz
            real :: delta

            bx => b(1)
            by => b(2)
            bz => b(3)

            if (by**2 + bz**2 > 0.0) then
                delta = rho*c**2 - bx**2
            else
                
                delta = 1.0
            end if
            !print *, "a=", a, "c=", c, "bx=", bx, "; by=", by, "; bz=", bz, "; delta=", delta, "rho=", rho

            ans = (/ rho, -c, c*bx*by/delta, c*bx*bz/delta, 0.0, rho*c**2*by/delta, rho*c**2*bz/delta, rho*a**2 /)
            !print *, "r3 before normalization", ans
            ans = ans/sqrt(dot_product(ans, ans))
        end function r3

        function r5(rho, b, c, a) result(ans)
            real, intent(in) :: rho, c, a
            real, target, dimension(3) :: b
            real, dimension(8) :: ans
            real, pointer :: bx, by, bz
            real :: delta

            bx => b(1)
            by => b(2)
            bz => b(3)
            if (by**2 + bz**2 > 0.0) then
                delta = rho*c**2 - bx**2
            else
                delta = 1
            end if

            !print *, "bx=", bx, "; by=", by, "; bz=", bz, "; delta=", delta

            ans = (/ rho, c, -c*bx*by/delta, c*bx*bz/delta, 0.0, rho*c**2*by/delta, rho*c**2*bz/delta, rho*a**2 /)
            !print *, "r5 before normalization", ans
            ans = ans/sqrt(dot_product(ans, ans))
        end function r5

        function r4() result(ans)
            real, dimension(8) :: ans
            ans = (/ 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
        end function r4

        function calc_eigenvalues(p) result(ans)
            ! Calculate the eigenvalues at every position given by the matrix p
            ! Asumes p is the averaged matrix at the interfaces of the cells
            real, target, dimension(:,:), intent(in) :: p
            real, dimension(size(p,2)) :: v, cs, ca, cf, a
            real, pointer, dimension(:,:) :: b
            real, pointer, dimension(:) :: pg, rho, bx
            real, dimension(7,size(p,2)) :: ans
            !       1    2   3   4   5   6   7   8
            ! P = [rho, vx, vy, vz, Bx, By, Bz, pg]

            rho => p(1,:)
            bx => p(5,:)
            pg => p(8,:)
            b => p(5:7,:)

            v = sqrt(sum(p(2:4,:)**2, dim=1))
            !v = p(2,:)

            ! a: sound speed
            a = sound_speed(pg, rho)
            ! ca: Alfven speed
            ca = alfven_speed(bx=bx, rho=rho)
            ! cf: Fast magnetosonic speed
            cf = fast_speed(a=a, va=ca, b=b, rho=rho, whole_b=1)
            ! cs: Slow magnetosonic speed
            cs = slow_speed(a=a, va=ca, b=b, rho=rho, whole_b=1)

            if (size(p,2) < 7) then
                !print *, "a=", a
                !print *, "ca=", ca
                !print *, "cf=", cf
                !print *, "cs=", cs
            end if

            ! lambda_1 = v - cf
            ans(1,:) = v - cf

            ! lambda 2: v - ca
            ans(2,:) = v - ca

            ! lambda 3: v - cs
            ans(3,:) = v - cs

            ! lambda 4: v
            ans(4,:) = v
            
            ! lambda 5: v + cs
            ans(5,:) = v + cs

            ! lambda 6: v + ca
            ans(6,:) = v + ca
        
            ! lambda 7: v + cf
            ans(7,:) = v + cf
    
        end function calc_eigenvalues


end submodule foamm_scheme 
