% Navier-Stokes solver,
% adapted for module SG212
% Computational Fluid Dynamics
%
% Depends on avg.m and DD.m
%


%------------------------------------

clear

% Parameters for lid-driven cavity flow
% The Richardson number is zero, i.e. passive scalar.

Pr = 0.71;     % Prandtl number
Re = 25;       % Reynolds number

dt = 0.001;    % Time step
Tf = 50;       % Final time
Lx = 1;        % Width of box
Ly = 1;        % Height of box
Nx = 20;       % Number of cells in x
Ny = 20;       % Number of cells in y
ig = 200;      % Number of iterations between output

% Boundary and initial conditions
Utop = 1.;
Ubottom = 0.;
Tbottom = 1.; Ttop = 0.;

%------------------------------------

% Number of iterations
Nit = Tf / dt;
% Spatial grid: Locations of cell corners
x = linspace(0,Lx,Nx+1);
y = linspace(0,Ly,Ny+1);

% Grid spacing
dx = Lx / Nx;
dy = Ly / Ny;
% Boundary conditions:
uN = x*0+Utop;          vN = zeros(1,Nx);
uS = x*0+Ubottom;       vS = zeros(1,Nx);
uW = zeros(1,Ny);       vW = zeros(1,Ny+1);
uE = zeros(1,Ny);       vE = zeros(1,Ny+1);
tN = ones(1,Nx+2)*Ttop; tS = ones(1,Nx+2)*Tbottom;
% Initial conditions
U = zeros(Nx-1,Ny); V = zeros(Nx,Ny-1);
% Linear profile for T
T = [];
for i=1:Nx
    T = [T;Tbottom+((avg(y,2)/Ly)*(Ttop-Tbottom))]; %#ok<*AGROW>
end
% Time series
tser = [];
Tser = [];

%------------------------------------

% Compute system matrices for pressure
% First set homogeneous Neumann condition all around
% Laplace operator on cell centres: Fxx + Fyy
Lp = kron(speye(Ny),DD(Nx,dx)) + kron(DD(Ny,dy),speye(Nx));
% Set one Dirichlet value to fix pressure at that point
Lp(1,:) = 0; Lp(1,1) = 1;

%------------------------------------

% Progress bar
fprintf(...
    '[         |         |         |         |         ]\n')

%------------------------------------

% Main loop over iterations

for k = 1:Nit

    % Include all boundary points for u and v (linear extrapolation
    % for ghost cells) into extended array (Ue,Ve)
    Ue = [uW;U;uE];
    Ve = [vS' V vN'];
    Ue = [(2*uS')-Ue(:,1) Ue (2*uN')-Ue(:,end)];
    Ve = [vW-Ve(1,:);Ve;vE-Ve(end,:)];

    % Averaged (Ua,Va) of u and v on corners
    Ua = avg(Ue,2);
    Va = avg(Ve,1);

    % Construct individual parts of nonlinear terms
    dUVdx = 1/dx * diff(Ua.*Va,1,1);
    dUVdy = 1/dy * diff(Ua.*Va,1,2);
    % For the remaining terms we need the velocity values at the pressure
    % nodes:
    Up = avg(avg(Ua,1),2); Vp = avg(avg(Va,1),2);
    dU2dx = 1/dx * diff(Up.*Up,1,1);
    dV2dy = 1/dy * diff(Vp.*Vp,1,2);

    % Treat viscosity explicitly
    viscu = diff(Ue(:,2:end-1),2,1)/dx^2 + diff(Ue(2:end-1,:),2,2)/dy^2;
    viscv = diff(Ve(:,2:end-1),2,1)/dx^2 + diff(Ve(2:end-1,:),2,2)/dy^2;

    % Compose final nonlinear term + explicit viscous terms
    U = U + ((dt/Re)*viscu) - (dt*(dU2dx+dUVdy(2:end-1,:)));
    V = V + ((dt/Re)*viscv) - (dt*(dUVdx(:,2:end-1)+dV2dy));

    % Pressure correction, Dirichlet P=0 at (1,1)
    rhs = (diff([uW;U;uE],1,1)/dx + diff([vS' V vN'],1,2)/dy)/dt;
    rhs = reshape(rhs,Nx*Ny,1);
    rhs(1) = 0;
    P = Lp \ rhs;
    P = reshape(P,Nx,Ny);

    % Apply pressure correction
    U = U - (dt*diff(P,1,1)/dx);
    V = V - (dt*diff(P,1,2)/dy);

    % Temperature equation
    Te = [T(1,:);T;T(end,:)];
    Te = [(2*tS')-Te(:,1) Te (2*tN')-Te(:,end)];
    Tu = avg(Te(:,2:end-1),1) .* Ue(:,2:end-1);
    Tv = avg(Te(2:end-1,:),2) .* Ve(2:end-1,:);
    H = -diff(Tu,1,1)/dx - diff(Tv,1,2)/dy + (diff(Te(:,2:end-1),2,1)/(dx^2)+diff(Te(2:end-1,:),2,2)/(dy^2))/(Pr*Re);
    T = T + dt*H;

    %------------------------------------

    % Progress bar
    if floor(51*k/Nit)>floor(51*(k-1)/Nit), fprintf('.'), end

    % Plot solution if needed
    if k==1||floor(k/ig)==k/ig

        % Compute divergence on cell centres
        if (1==1)
            div = diff([uW;U;uE])/dx + diff([vS' V vN'],1,2)/dy;

            figure(1); clf; hold on;
            contourf(avg(x,2),avg(y,2),div'); colorbar
            axis equal; axis([0 Lx 0 Ly]);
            % title(sprintf('Divergence at t=%g',k*dt))
            drawnow
        end

        % Compute velocity on cell corners
        % The velocity values have changed, so we need to reintroduce the
        % boundary values to get the extended grid.
        Ue = [uW;U;uE];
        Ve = [vS' V vN'];
        Ue = [(2*uS')-Ue(:,1) Ue (2*uN')-Ue(:,end)];
        Ve = [vW-Ve(1,:);Ve;vE-Ve(end,:)];
        Ua = avg(Ue,2);
        Va = avg(Ve,1);
        Len = sqrt(Ua.^2+Va.^2+eps);

        figure(2); clf; hold on;
        % contourf(avg(x,2),avg(y,2),P'); colorbar
        contourf(x,y,sqrt(Ua.^2+Va.^2)',20,'k-'); colorbar
        quiver(x,y,(Ua./Len)',(Va./Len)',.4,'k-')
        axis equal; axis([0 Lx 0 Ly]);
        % title(sprintf('u at t=%g',k*dt))
        drawnow

        % Compute temperature on cell corners
        % The temperature is only defined at cell centres, so we need to
        % reintroduce the boundary values
        Te = [T(1,:);T;T(end,:)];
        Te = [(2*tS')-Te(:,1) Te (2*tN')-Te(:,end)];
        Ta = avg(avg(Te,1),2);

        figure(3); clf; hold on;
        contourf(x,y,Ta',20,'k-'); colorbar
        quiver(x,y,(Ua./Len)',(Va./Len)',.4,'k-')
        axis equal; axis([0 Lx 0 Ly]);
        % title(sprintf('T at t=%g',k*dt))
        drawnow

        % Time history
        if (1==1)
            figure(5); hold on;
            tser = [tser k*dt];
            Tser = [Tser Ue(ceil((Nx+1)/2),ceil((Ny+1)/2))];
            plot(tser,abs(Tser))
            % title(sprintf('Probe signal at x=%g, y=%g',...
            %       x(ceil((Nx+1)/2)),y(ceil((Ny+1)/2))))
            set(gca,'yscale','log')
            xlabel('time t'); ylabel('u(t)');
        end
    end
end
fprintf('\n')
