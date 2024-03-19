! Created by cf on 31/07/23.

module model

    private
    public calc_primitives, calc_cons_quantities, init_cons_quantities, gamma, f, &
            sound_speed, alfven_speed, fast_speed, slow_speed, pt, pg

    real :: gamma

    contains

        subroutine f(u, flux)
            real, target, intent(in) :: u(:,:)
            real, intent(out) :: flux(:,:)
            real, dimension(:), pointer :: rho, e, Bx, By, Bz
            real, pointer :: B(:,:)
            real, dimension(size(u,2)) :: vx, vy, vz, pt_val
            real, dimension(3, size(u,2)) :: v

            rho => u(1,:)
            vx = u(2,:)/rho
            vy = u(3,:)/rho
            vz = u(4,:)/rho

            v(1,:) = vx
            v(2,:) = vy
            v(3,:) = vz

            Bx => u(5,:)
            By => u(6,:)
            Bz => u(7,:)
            B => u(5:7,:)
            e => u(8,:)
            !print *, "after pointer assignation in f"

            pt_val = pt(pg(e, B, vx, vy, vz, rho), B)
            !print *, "after pt calculation in f"

            ! rho*vx
            flux(1,:) = u(2,:)
            ! rho*vx^2 + pt
            flux(2,:) = rho*vx*vx + pt_val - Bx*Bx
            ! rho*vx*vy - Bx*By
            flux(3,:) = rho*vx*vy - Bx*By
            ! rho*vx*vz - Bx*Bz
            flux(4,:) = rho*vx*vz - Bx*Bz
            ! 0
            flux(5,:) = 0
            ! By*vx - Bx*vy
            flux(6,:) = By*vx - Bx*vy
            ! Bz*vx - Bx*vz
            flux(7,:) = Bz*vx - Bx*vz
            ! vx*(e+pt) - Bx*(B@v)
            flux(8,:) = vx*(e + pt_val) - Bx*(sum(B*v, 1))
            !print *, "end of f subroutine"

        end subroutine f

        function pt(pg, B) result(ans)
            real, intent(in) :: pg(:)
            real, intent(in) :: B(:,:)
            real :: ans(size(pg))

            ans = pg + 0.5*sum(B*B, 1)

        end function pt

        function pg(e, B, vx, vy, vz, rho) result(ans)
            real, intent(in) :: e(:), rho(:), B(:,:)
            real, dimension(:), intent(in) :: vx, vy, vz
            real :: ans(size(e))

            ans = (gamma - 1.0)*(e - 0.5*sum(B*B, 1) - 0.5*rho*(vx*vx + vy*vy + vz*vz))

        end function pg

        function calc_e(B, v, pg, rho) result(ans)
            real, dimension(:,:), intent(in) :: B, v
            real, dimension(:), intent(in) :: pg, rho
            real :: ans(size(rho))

            ans = 0.5*sum(B*B, 1) + 0.5*rho*sum(v*v, 1) + pg/(gamma - 1.0)

        end function calc_e

        subroutine calc_primitives(u, p)
            real, target, intent(in) :: u(:,:)
            real, intent(out) :: p(:,:)
            real, pointer :: rho(:)
            integer :: rows, cols, i
            ! U = [rho, rho*vx, rho*vy, rho*vz, Bx, By, Bz, e]
            ! P = [rho, vx, vy, vz, Bx, By, Bz, pg]

            rho => u(1,:)

            p(1,:) = u(1,:)
            p(2,:) = u(2,:)/rho
            p(3,:) = u(3,:)/rho
            p(4,:) = u(4,:)/rho
            p(5:7,:) = u(5:7,:)
            p(8,:) = pg(e=u(8,:), B=u(5:7,:), vx=p(2,:), vy=p(3,:), vz=p(4,:), rho=rho)

            cols = size(u,2)
            if (cols < 1) then
                print *, "TEST CALCULATE PRIMITIVES"
                rows = size(u,1)
                do i = 1, rows
                    write(*,*) p(i,:)
                end do
                print *, "----------------------------------------"
            end if
        end subroutine calc_primitives

        subroutine calc_cons_quantities(p, u)
            real, target, intent(in) :: p(:,:)
            real, intent(out) :: u(:,:)
            real, pointer :: rho(:)
            ! U = [rho, rho*vx, rho*vy, rho*vz, Bx, By, Bz, e]
            ! P = [rho, vx, vy, vz, Bx, By, Bz, pg]

            rho => p(1,:)
            u(1,:) = rho
            u(2,:) = p(2,:)*rho
            u(3,:) = p(3,:)*rho
            u(4,:) = p(4,:)*rho
            u(5:7,:) = p(5:7,:)
            u(8,:) = calc_e(B=p(5:7,:), v=p(2:4,:), pg=p(8,:), rho=rho)

        end subroutine calc_cons_quantities

        subroutine init_cons_quantities(p, u)
            real, intent(in) :: p(:)
            real, intent(out) :: u(:,:)
            real :: e
            integer :: i
            ! U = [rho, rho*vx, rho*vy, rho*vz, Bx, By, Bz, e]
            ! P = [rho, vx, vy, vz, Bx, By, Bz, pg]

            print *, "p(1)=", p(1)
            print *, "p(2)=", p(2)

            u(1,:) = p(1)
            u(2,:) = p(2)*p(1)
            u(3,:) = p(3)*p(1)
            u(4,:) = p(4)*p(1)
            !u(2:4,:) = transpose(reshape(p(2:4)*p(1), (/size(u, 2), 3/)))
            !u(5:7,:) = reshape(p(5:7), (/3, 1/))
            u(5,:) = p(5)
            u(6,:) = p(6)
            u(7,:) = p(7)

            e = 0.5*dot_product(p(5:7), p(5:7)) + 0.5*p(1)*dot_product(p(2:4), p(2:4)) + p(8)/(gamma - 1.0)
            u(8,:) = e

            if (size(u, 2) < 7) then
                do i = 1,8
                    write(*,*) u(i,:)
                end do
            end if
        end subroutine init_cons_quantities

        function sound_speed(pg, rho) result(ans)
            real, dimension(:), intent(in) :: pg, rho
            real, dimension(size(pg)) :: ans

            ans = sqrt(gamma*pg/rho)

        end function sound_speed

        function alfven_speed(Bx, rho) result(ans)
            real, dimension(:), intent(in) :: Bx, rho
            real, dimension(size(rho)) :: ans

            ans = abs(Bx)/sqrt(rho)

        end function alfven_speed

        function fast_speed(a, va, b, rho, whole_b) result(ans)
            real, dimension(:), intent(in) :: a, va, rho
            real, dimension(3, size(a)), intent(in) :: b
            integer, intent(in) :: whole_b
            real, dimension(size(a)) :: ans, a2, va2, b2

            a2 = a**2
            va2 = va**2
            b2 = sum(b*b, dim=1)
            
            if (whole_b == 1) then
                ans = sqrt(0.5*(a2 + b2/rho + sqrt((a2 + b2/rho)**2 - 4*va2*a2)))
            else
                !ans = sqrt(0.5*(a2 + va2 + sqrt((a2 + va2)**2 - 4.0*va2*a2)))
                ans = 0.5*sqrt(a2 + va2 - sqrt((a2 + va2)**2 - 4*va2*a2))
            end if

        end function fast_speed

        function slow_speed(a, va, b, rho, whole_b) result(ans)
            real, dimension(:), intent(in) :: a, va, rho
            real, dimension(3, size(a)), intent(in) :: b
            integer, intent(in) :: whole_b
            real, dimension(size(a)) :: ans, a2, va2, b2

            a2 = a**2
            va2 = va**2
            b2 = sum(b*b, dim=1)

            if (whole_b == 1) then
                ans = sqrt(0.5*(a2 + b2/rho - sqrt((a2 + b2/rho)**2 - 4*va2*a2)))
            else
                !ans = sqrt(0.5*(a2 + va2 - sqrt((a2 + va2)**2 - 4.0*va2*a2)))
                ans = 0.5*sqrt(a2 + va2 - sqrt((a2 + va2)**2 - 4*va2*a2))
            end if

        end function slow_speed


end module model
