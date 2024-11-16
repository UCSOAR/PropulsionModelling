% Composite cylinder stresses and strains for decoupled load cases

% TUBE GEOMETRY
R = 100;
L = 200;

% Layer 1
h1 = 0.62;
phi1 = 36.*pi./180;
E11 = 44;
E21 = 9.4;
G121 = 4;
v211 = 0.26;
v121 = 0.26;
G131 = 4;
G231 = 3;

% Layer 2
h2 = 0.60;
phi2 = pi./2;
E12 = 44;
E22 = 9.4;
G122 = 4;
v212 = 0.26;
v122 = 0.26;
G132 = 4;
G232 = 3;

% CONDITIONS
h = [h1;h2];
phi = [phi1;phi2];
E = [E11,E21; E12,E22];
G = [G121,G131,G231; G122,G132,G232];
v = [v121,v211; v122,v212];



%%%%%%%%%%%%%%%%%%%%%%%%%%% TESTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
z = h1;
P = 40;

statesAxial = [];
statesTorsion = [];
statesPressure = [];
for x = 0:1:L/2
    statesAxial = [statesAxial; Axial(x,z,P,R,h,phi,E,G,v)];
    statesTorsion = [statesTorsion; Torsion(x,z,P,R,h,phi,E,G,v)];
    %statesPressure = [statesPressure; Pressure(x,z,P,R,h,phi,E,G,v)];
end

linePlotStress(L,z,P,R,h,phi,E,G,v,1)

% Plot stress colored cylinder
stressAxial(L,z,P,R,h,phi,E,G,v,1)

% plot the deformed cylinder
strainAxial(L,z,P,R,h,phi,E,G,v,1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%% PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Emod = Emodstiffness(E,v,G)
    % modified modulus
    E1mod = E(:,1)./(1-v(:,1).*v(:,2));
    E2mod = E(:,2)./(1-v(:,1).*v(:,2));
    E12 = E1mod.*v(:,1) +2.*G(:,1);

    Emod = [E1mod,E2mod,E12];
end


function A = Astiffness(E,v,G,phi,Emod)
    % directions
    c = cos(phi);
    s = sin(phi);

    % stiffness coefficients
    A11 = Emod(:,1).*(c.^4) + Emod(:,2).*(s.^4) + 2.*Emod(:,3).*(c.^2).*(s.^2);
    A12 = Emod(:,1).*v(:,1) + (Emod(:,1)+Emod(:,2)-2.*Emod(:,3)).*(c.^2).*(s.^2);
    A21 = A12;
    A14 = (Emod(:,1).*(c.^2)-Emod(:,2).*(s.^2)-Emod(:,3).*((c.^2)-(s.^2))).*c.*s;
    A41 = A14;
    A22 = Emod(:,1).*(s.^4)+Emod(:,2).*(c.^4)+2.*Emod(:,3).*(c.^2).*(s.^2);
    A24 = (Emod(:,1).*(s.^2)-Emod(:,2).*(c.^2)+Emod(:,3).*((c.^2)-(s.^2))).*c.*s;
    A42 = A24;
    A44 = (Emod(:,1)+Emod(:,2)-2.*Emod(:,1).*v(:,1)).*(c.^2).*(s.^2) + G(:,1).*((c.^2)-(s.^2)).^2;
    A55 = G(:,2).*(c.^2)+G(:,3).*(s.^2);
    A56 = (G(:,2)-G(:,3)).*c.*s;
    A65 = A56;
    A66 = G(:,2).*(s.^2)+G(:,3).*(c.^2);

    A = [A11,A12,A21,A14,A41,A22,A24,A42,A44,A55,A56,A65,A66];

end


function I = Istiffness(A,h)
    % define layers
    t = 0;
    I0 = [];
    I1 = [];
    I2 = [];

    % Thin layer 0 I coefficients
    for i = 1:length(h)
        I0 = [I0; A(i,:).*h(i)];
    end
    I0 = sum(I0,1);
    
    % Thin layer 1 I coefficients
    for i = 1:length(h)
        t = [t,t(i)+h(i)];
        I1 = [I1; 0.5.*A(i,:).*h(i).*(t(i+1)-t(i))];
    end
    I1 = sum(I1,1);

    % Thin layer 2 I coefficients
    for i = 1:length(h)
        t = [t,t(i)+h(i)];
        I2 = [I2; 1./3.*A(i,:).*h(i).*(t(i+1).^2 + t(i+1).*t(i) + t(i).^2)];
    end
    I2 = sum(I2,1);

    I = [I0;I1;I2];
end
    

function BCD = BCDstiffness(I,e)
    % B stiffness
    B = I(1,:);

    % C stiffness
    C = I(2,:)-e.*I(1,:);

    % D stiffness
    D = I(3,:)-2.*e.*I(2,:)+(e.^2).*I(1,:);

    BCD = [B;C;D];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%% LOADS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function state = Axial(x,z,P,R,h,phi,E,G,v)

    % Stiffnesses
    Emod = Emodstiffness(E,v,G);
    A = Astiffness(E,v,G,phi,Emod);
    I = Istiffness(A,h);
    BCD = BCDstiffness(I,0);
    B = BCD(1,:);
    C = BCD(2,:);
    D = BCD(3,:);

%     phi1 = 36.*pi./180;
%     h2 = 0.60;
%     phi2 = pi./2;
%     R = 100;
%     
%     E11 = 44;
%     E21 = 9.4;
%     G12 = 4;
%     v211 = 0.26;
%     v121 = 0.26;
%     G13 = 4;
%     G23 = 3;
%     
%     A551 = 3.7;
%     A552 = 3;

    D(1,1) = I(3,1)-((I(2,1).^2)./I(1,1));
    htotal = sum(h);
    Sx = (htotal.^2)./sum(h./A(:,10));

%     D11 = 16.9;
%     Sx = 4.05;
    
%     A111 = 25;
%     A121 = 10;
%     A221 = 14.1;
%     A441 = 11.5;
%     
%     A112 = 9.5;
%     A122 = 2.5;
%     A222 = 44.7;
%     A442 = 4;
%     
%     
%     I110 = 21.2;
%     I120 = 7.7;
%     I220 = 35.6;
%     I440 = 9.5;
%     
%     I111 = 10.1;
%     I121 = 3.3;
%     I221 = 27.4;
%     I441 = 4.4;
%     
%     I112 = 21.7;
%     I122 = 5.9;
%     I222 = 94;
%     
%     
%     B11 = 21.2;
%     B12 = 7.7;
%     B22 = 35.6;
%     B21 = B12;
%     C11 = 10.1;
%     C12 = 3.3;
%     C22 = 27.4;
%     D12 = 5.9;
%     D22 = 94;
%     
%     
%     B12mod = 7.7;
%     B22mod = 35.3;
%     C12mod = 3.2;
%     C22mod = 26.5;

    B12mod = B(:,2) - C(:,2)./R;
    B22mod = B(:,6) - C(:,6)./R;
    C12mod = C(:,2) - D(:,2)./R;
    C22mod = C(:,6) - D(:,6)./R;
    BB = B(:,1).*B22mod - B12mod.*B(:,3);

%     E1mod1 = E11./(1-v121.*v211);
%     E2mod1 = E21./(1-v121.*v211);

    alphaSQ = 1./(2.*R).*(C(:,3)./D(:,1)+BB./(B(:,1).*Sx.*R));
    betaSQ = (BB./(D(:,1).*B(:,1).*(R.^2))).^(1./2);

    r = (1./2.*(betaSQ-alphaSQ)).^(1./2);
    t = (1./2.*(betaSQ+alphaSQ)).^(1./2);

%     r = 7.9./R;
%     t = 8.75./R;

    ex0 = -P.*B22mod./(2.*pi.*R.*BB).*(1+exp(-r.*x).*(0.11.*sin(t.*x)-0.052.*cos(t.*x)));
    ey0 = P.*B(:,3)./(2.*pi.*R.*BB).*(1+exp(-r.*x).*(0.51.*sin(t.*x)-0.24.*cos(t.*x)));
    gammaxy0 = 0;
    %thetax = P.*B(:,3)./(2.*pi.*R.*BB).*exp(-r.*x).*(6.3.*cos(t.*x)-2.3.*sin(t.*x));

    syms xvar
    thetaxPrime = diff(P.*B(:,3)./(2.*pi.*R.*BB).*exp(-r.*xvar).*(6.3.*cos(t.*xvar)-2.3.*sin(t.*xvar)));
    thetaPrime = subs(thetaxPrime, xvar, x);
    thetaPrime = eval(thetaPrime);

    kx = [];
    ky = [];
    kxy = [];
    ex = [];
    ey = [];
    gammaxy = [];
    stressx = [];
    stressy = [];
    shearxy = [];
    e1 = [];
    e2 = [];
    gamma12 = [];
    stress1 = [];
    stress2 = [];
    shear12 = [];

    for layer = 1:length(z)

        kx = [kx; thetaPrime];
        ky = [ky; -ey0./R];
        kxy = [kxy; 0];

        ex = [ex;ex0 + z(layer).*kx];
        ey = [ey;ey0 + z(layer).*ky];
        gammaxy = [gammaxy;gammaxy0 + z(layer).*kxy];

        
        stressx = [stressx; A(layer,1).*ex(layer) + A(layer,2).*ey(layer)];
        stressy = [stressy; 0];
        %shearxy = [shearxy; Gxy(layer).*gammaxy(layer)];

        e1 = [e1; ex(layer).*((cos(phi(layer))).^2) + ey(layer).*((sin(phi(layer))).^2) + gammaxy(layer).*sin(phi(layer)).*cos(phi(layer))];
        e2 = [e2; ex(layer).*((sin(phi(layer))).^2) + ey(layer).*((cos(phi(layer))).^2) - gammaxy(layer).*sin(phi(layer)).*cos(phi(layer))];
        %gamma12 = [gamma12; 2.*(ey(layer)-ex(layer)).*sin(phi(layer)).*cos(phi(layer)) + gammaxy(layer).*cos(2.*phi(layer))];

        stress1 = [stress1; Emod(layer,1).*(e1(layer)+v(layer,1).*e2(layer))];
        stress2 = [stress2; Emod(layer,2).*(e2(layer)+v(layer,2).*e1(layer))];
        %shear12 = [shear12; G(layer,1).*gamma12(layer)];
    end

    shearxz = P.*B22mod./(2.*pi.*BB.*R.^2).*exp(-r.*x).*((23.75.*cos(t.*x)-9.75.*sin(t.*x)).*z+(6.3.*cos(t.*x)+24.9.*sin(t.*x)).*(z.^2));
    stressz = -0.068.*P./(htotal.*R.^2).*(z+exp(-r.*x).*((0.18.*cos(t.*x)-0.0725.*sin(t.*x)).*z-(0.12.*cos(t.*x)+0.059.*sin(t.*x)).*(z.^2)+(0.05.*cos(t.*x)-0.076.*sin(t.*x).*(z.^3))));
    
    state = [x,kx,ky,ex,ey,stressx,stressy,e1,e2,stress1,stress2,stressz,shearxz];

end


function state = Torsion(x,z,P,R,h,phi,E,G,v)
    % Stiffnesses
    Emod = Emodstiffness(E,v,G);
    A = Astiffness(E,v,G,phi,Emod);
    I = Istiffness(A,h);
    BCD = BCDstiffness(I,I(2,9)/I(1,9));
    B = BCD(1,:);
    C = BCD(2,:);
    D = BCD(3,:);

    D(1,1) = I(3,1)-((I(2,1).^2)./I(1,1));
    htotal = sum(h);
    Sx = (htotal.^2)./sum(h./A(:,10));

    B12mod = B(:,2) - C(:,2)./R;
    B22mod = B(:,6) - C(:,6)./R;
    C12mod = C(:,2) - D(:,2)./R;
    C22mod = C(:,6) - D(:,6)./R;
    BB = B(:,1).*B22mod - B12mod.*B(:,3);

    gammaxy0 = P/(2*pi*(R^2)*B(:,9));

    kx = [];
    ky = [];
    kxy = [];
    ex = [];
    ey = [];
    gammaxy = [];
    stressx = [];
    stressy = [];
    shearxy = [];
    e1 = [];
    e2 = [];
    gamma12 = [];
    stress1 = [];
    stress2 = [];
    shear12 = [];

    for layer = 1:length(z)

        kx = [kx; 0];
        ky = [ky; 0];
        kxy = [kxy; 0];

        ex = [ex;0 + z(layer).*kx];
        ey = [ey;0 + z(layer).*ky];
        gammaxy = [gammaxy;gammaxy0 + z(layer).*kxy];

        
        stressx = [stressx; 0];
        stressy = [stressy; 0];
        shearxy = [shearxy; ex(layer).*A(:,5) + ey(layer).*A(:,8) + gammaxy(layer).*A(:,9)];

        e1 = [e1; ex(layer).*((cos(phi(layer))).^2) + ey(layer).*((sin(phi(layer))).^2) + gammaxy(layer).*sin(phi(layer)).*cos(phi(layer))];
        e2 = [e2; ex(layer).*((sin(phi(layer))).^2) + ey(layer).*((cos(phi(layer))).^2) - gammaxy(layer).*sin(phi(layer)).*cos(phi(layer))];
        gamma12 = [gamma12; 2.*(ey(layer)-ex(layer)).*sin(phi(layer)).*cos(phi(layer)) + gammaxy(layer).*cos(2.*phi(layer))];

        stress1 = [stress1; Emod(layer,1).*(e1(layer)+v(layer,1).*e2(layer))];
        stress2 = [stress2; Emod(layer,2).*(e2(layer)+v(layer,2).*e1(layer))];
        shear12 = [shear12; G(layer,1).*gamma12(layer)];
    end
    
    state = [x,gammaxy,gamma12,shear12];


end


function state = Pressure(x,z,P,R,h,phi,E,G,v)
    % Stiffnesses
    Emod = Emodstiffness(E,v,G);
    A = Astiffness(E,v,G,phi,Emod);
    I = Istiffness(A,h);
    BCD = BCDstiffness(I,0);
    B = BCD(1,:);
    C = BCD(2,:);
    D = BCD(3,:);

    BB = B(:,1).*B(:,6) - (B(:,2).^2);

    kx = [];
    ky = [];
    kxy = [];
    ex = [];
    ey = [];
    gammaxy = [];
    stressx = [];
    stressy = [];
    shearxy = [];
    e1 = [];
    e2 = [];
    gamma12 = [];
    stress1 = [];
    stress2 = [];
    shear12 = [];

    syms ex0 ey0
    eqn1 = (1/2*P*R - B(:,2).*(ey0))./B(:,1) == ex0;
    eqn2 = (P*R-B(:,2).*ex)./B(:,6) == ey0;
    soln = solve([eqn1,eqn2],[ex0,ey0]);

    ex0 = double(soln.ex0);
    ey0 = double(soln.ey0);
    gammaxy0 = 0;

    for layer = 1:length(z)

        kx = [kx; 0];
        ky = [ky; 0];
        kxy = [kxy; 0];

        ex = [ex;ex0 + z(layer).*kx];
        ey = [ey;ey0 + z(layer).*ky];
        gammaxy = [gammaxy;gammaxy0 + z(layer).*kxy];

        
        stressx = [stressx; A(layer,1).*ex(layer) + A(layer,2).*ey(layer)+A(layer,4).*gammaxy(layer)];
        stressy = [stressy; A(layer,3).*ex(layer) + A(layer,6).*ey(layer)+A(layer,7).*gammaxy(layer)];
        shearxy = [shearxy; A(layer,5).*ex(layer) + A(layer,8).*ey(layer)+A(layer,9).*gammaxy(layer)];

        e1 = [e1; ex(layer).*((cos(phi(layer))).^2) + ey(layer).*((sin(phi(layer))).^2) + gammaxy(layer).*sin(phi(layer)).*cos(phi(layer))];
        e2 = [e2; ex(layer).*((sin(phi(layer))).^2) + ey(layer).*((cos(phi(layer))).^2) - gammaxy(layer).*sin(phi(layer)).*cos(phi(layer))];
        gamma12 = [gamma12; 2.*(ey(layer)-ex(layer)).*sin(phi(layer)).*cos(phi(layer)) + gammaxy(layer).*cos(2.*phi(layer))];

%         stress1 = [stress1; Emod(layer,1).*(e1(layer)+v(layer,1).*e2(layer))];
%         stress2 = [stress2; Emod(layer,2).*(e2(layer)+v(layer,2).*e1(layer))];
%         shear12 = [shear12; G(layer,1).*gamma12(layer)];
        check = P*R/(2*(35.6*21.2-7.7^2))*Emod(layer,2)*((35.6-2*7.7)*((cos(phi(layer))^2)+v(layer,1)*(sin(phi(layer))^2))+(2*21.2-7.7)*((sin(phi(layer))^2)+v(layer,1)*(cos(phi(layer))^2)));
        stress1 = [stress1; P*R/(2*BB)*Emod(layer,2)*((B(layer,6)-2*B(layer,2))*((cos(phi(layer))^2)+v(layer,1)*(sin(phi(layer))^2))+(2*B(layer,1)-B(layer,2))*((sin(phi(layer))^2)+v(layer,1)*(cos(phi(layer))^2)))];
        stress2 = [stress2; P*R/(2*BB)*Emod(layer,2)*((B(layer,6)-2*B(layer,2))*((sin(phi(layer))^2)+v(layer,2)*(cos(phi(layer))^2))+(2*B(layer,1)-B(layer,2))*((cos(phi(layer))^2)+v(layer,2)*(sin(phi(layer))^2)))];
        shear12 = [shear12; P*R/(2*BB)*G(layer,1)*(2*B(layer,1)+B(layer,2)-B(layer,6))*sin(2*phi(layer))];

    end
    
    state = [x,kx,ky,ex,ey,stressx,stressy,e1,e2,stress1,stress2];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%% PLOTTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stressAxial(L,z,P,R,h,phi,E,G,v,layernumber)
    kx = [];
    ky = [];
    %kxy = [];
    ex = [];
    ey = [];
    %gammaxy = [];
    stressx = [];
    stressy = [];
    %shearxy = [];
    e1 = [];
    e2 = [];
    %gamma12 = [];
    stress1 = [];
    stress2 = [];
    %shear12 = [];
    stressz = [];
    shearxz = [];

    for x = 0:1:L/2
        state = Axial(x,z,P,R,h,phi,E,G,v);

        kx = [kx;state(layernumber,2)];
        ky = [ky;state(layernumber,3)];
        ex = [ex;state(layernumber,4)];
        ey = [ey;state(layernumber,5)];
        stressx = [stressx;state(layernumber,6)];
        stressy = [stressy;state(layernumber,7)];
        e1 = [e1;state(layernumber,8)];
        e2 = [e2;state(layernumber,9)];
        stress1 = [stress1;state(layernumber,10)];
        stress2 = [stress2;state(layernumber,11)];
        stressz = [stressz;state(layernumber,12)];
        shearxz = [shearxz;state(layernumber,13)];
    end

    vonMises = (0.5*((stress1-stress2).^2+(stress2-stressz).^2+(stress1-stressz).^2)).^(1/2);

    %plot(0:1:L/2,vonMises.')


    nTheta = 100;  % Number of angular divisions
    nZ = length(vonMises);  % Number of axial divisions based on the stress data size
    
    % Create the grid for the cylinder surface
    theta = linspace(0, 2*pi, nTheta);  % Angular positions
    z1 = linspace(0, L/2, nZ);           % Axial positions (half the length)
    z2 = linspace(L/2, L, nZ);
    z = [z1, z2];
    [Theta, Z] = meshgrid(theta, z);    % Create grid for angular and axial positions
    
    % Convert to Cartesian coordinates
    X = R * cos(Theta);  % X-coordinates of the cylinder surface
    Y = R * sin(Theta);  % Y-coordinates of the cylinder surface
    
    % Replicate stress values for the full cylinder
    StressCyl = repmat(vonMises(:), 1, nTheta);  % Repeat stress along the angular direction
    StressCylFull = [StressCyl; repmat(flip(vonMises(:)), 1, nTheta)];
    
    % Create a full stress matrix by mirroring the half
    %StressCylFull = [StressCyl, flip(StressCyl, 2)];  % Mirror the stress data
    
    % Mirror the X and Y coordinates to match the full cylinder
    %XFull = [X, -X];  % Mirror X coordinates
    %YFull = [Y, Y];   % Repeat Y coordinates
    %ZFull = [Z, Z + L/2];  % Extend Z coordinates to cover the full length
    
    % Plot the cylinder with mirrored stress color scale
    figure;
    surf(X, Y, Z, StressCylFull);  % Plot cylinder with stress mapped to color
    shading interp;                   % Smooth color interpolation
    colormap(jet);                    % Set colormap (customize if needed)
    c = colorbar;                         % Display color scale
    c.Label.String = 'Von Mises Stress (MPa)';
    caxis([min(vonMises) max(vonMises)]);  % Set color axis limits based on stress values
    
    % Enhance plot appearance
    xlabel('X (mm)');
    ylabel('Y (mm)');
    zlabel('Z (mm)');
    title('Composite Tube Stress Distribution');
    axis equal;                       % Ensure equal scaling of axes
    set(gca, 'FontSize', 12, 'FontWeight', 'bold');  % Improve axis font style

end


function strainAxial(L,z,P,R,h,phi,E,G,v,layernumber)
    kx = [];
    ky = [];
    %kxy = [];
    ex = [];
    ey = [];
    %gammaxy = [];
    stressx = [];
    stressy = [];
    %shearxy = [];
    e1 = [];
    e2 = [];
    %gamma12 = [];
    stress1 = [];
    stress2 = [];
    %shear12 = [];
    stressz = [];
    shearxz = [];

    x = 0:1:L/2;
    for x = 0:1:L/2
        state = Axial(x,z,P,R,h,phi,E,G,v);

        kx = [kx;state(layernumber,2)];
        ky = [ky;state(layernumber,3)];
        ex = [ex;state(layernumber,4)];
        ey = [ey;state(layernumber,5)];
        stressx = [stressx;state(layernumber,6)];
        stressy = [stressy;state(layernumber,7)];
        e1 = [e1;state(layernumber,8)];
        e2 = [e2;state(layernumber,9)];
        stress1 = [stress1;state(layernumber,10)];
        stress2 = [stress2;state(layernumber,11)];
        stressz = [stressz;state(layernumber,12)];
        shearxz = [shearxz;state(layernumber,13)];
    end
    
    vonMises = (0.5*((stress1-stress2).^2+(stress2-stressz).^2+(stress1-stressz).^2)).^(1/2);

    % new dimensions after deformation
    dy = ey.*R;
    x = 0:1:L/2;
    xdef = x(2);
    for i=2:length(ex)
        %dx = ex(i-1).*xdef(i-1);
        dx = ex(i-1).*(x(i)-x(i-1));
        xdef = [xdef xdef(i-1)+dx+(x(i)-x(i-1))];
    end
    xdef2 = xdef(end);
    ex2 = flip(ex);
    for i=2:length(ex)
        %dx = ex2(i-1).*xdef2(i-1);
        dx = ex2(i-1).*(x(i)-x(i-1));
        xdef2 = [xdef2 xdef2(i-1)+dx+(x(i)-x(i-1))];
    end
    ydef = linspace(R,R,length(dy))+100*dy.';
    ydef = [ydef,flip(ydef)];

    %xdef = [x + (ex.').*L, x + x(end) + flip(ex.').*L];
    size(xdef)
    size(ydef)

    nTheta = 100;  % Number of angular divisions
    nZ = length(x);  % Number of axial divisions based on the stress data size
    
    % Create the grid for the cylinder surface
    theta = linspace(0, 2*pi, nTheta);  % Angular positions
    z1 = xdef;        % Axial positions (half the length)
    z2 = xdef2;
    z = [z1, z2];
    %z = xdef;
    [Theta, Z] = meshgrid(theta, z);    % Create grid for angular and axial positions
    

    % Convert to Cartesian coordinates
    X = [];
    Y = [];
    for i = 1:length(z)
        X = [X;((ydef(i)).*cos(Theta(i,:)))];  % X-coordinates of the cylinder surface
        Y = [Y;((ydef(i)).*sin(Theta(i,:)))];  % Y-coordinates of the cylinder surface
    end
    
    % Replicate stress values for the full cylinder
    StressCyl = repmat(vonMises(:), 1, nTheta);  % Repeat stress along the angular direction
    StressCylFull = [StressCyl; repmat(flip(vonMises(:)), 1, nTheta)];
    
    % Plot the cylinder with mirrored stress color scale
    figure;
    surf(X, Y, Z, StressCylFull);  % Plot cylinder with stress mapped to color
    shading interp;                   % Smooth color interpolation
    colormap(jet);                    % Set colormap (customize if needed)
    c = colorbar;                         % Display color scale
    c.Label.String = 'Von Mises Stress (MPa)';
    caxis([min(vonMises) max(vonMises)]);  % Set color axis limits based on stress values
    
    % Enhance plot appearance
    xlabel('X (mm)');
    ylabel('Y (mm)');
    zlabel('Z (mm)');
    title('Composite Tube Deformation');
    axis equal;                       % Ensure equal scaling of axes
    set(gca, 'FontSize', 12, 'FontWeight', 'bold');  % Improve axis font style
end

function linePlotStress(L,z,P,R,h,phi,E,G,v,layernumber)
    kx = [];
    ky = [];
    %kxy = [];
    ex = [];
    ey = [];
    %gammaxy = [];
    stressx = [];
    stressy = [];
    %shearxy = [];
    e1 = [];
    e2 = [];
    %gamma12 = [];
    stress1 = [];
    stress2 = [];
    %shear12 = [];
    stressz = [];
    shearxz = [];

    for x = 0:1:L/2
        stateAxial = Axial(x,z,P,R,h,phi,E,G,v);

        kx = [kx;stateAxial(layernumber,2)];
        ky = [ky;stateAxial(layernumber,3)];
        ex = [ex;stateAxial(layernumber,4)];
        ey = [ey;stateAxial(layernumber,5)];
        stressx = [stressx;stateAxial(layernumber,6)];
        stressy = [stressy;stateAxial(layernumber,7)];
        e1 = [e1;stateAxial(layernumber,8)];
        e2 = [e2;stateAxial(layernumber,9)];
        stress1 = [stress1;stateAxial(layernumber,10)];
        stress2 = [stress2;stateAxial(layernumber,11)];
        stressz = [stressz;stateAxial(layernumber,12)];
        shearxz = [shearxz;stateAxial(layernumber,13)];
    end

    vonMisesAxial = (0.5*((stress1-stress2).^2+(stress2-stressz).^2+(stress1-stressz).^2)).^(1/2);


    shear12 = [];
    gammaxy = [];
    gamma12 = [];
    for x = 0:1:L/2
        stateTorsion = Torsion(x,z,P,R,h,phi,E,G,v);

        gammaxy = [gammaxy;stateTorsion(layernumber,2)];
        gamma12 = [gamma12;stateTorsion(layernumber,3)];
        shear12 = [shear12;stateTorsion(layernumber,4)];
    end

    vonMisesTorsion = shear12;

%     for x = 0:1:L/2
%         statePressure = Pressure(x,z,P,R,h,phi,E,G,v);
% 
%         kx = [kx;statePressure(layernumber,2)];
%         ky = [ky;statePressure(layernumber,3)];
%         ex = [ex;statePressure(layernumber,4)];
%         ey = [ey;statePressure(layernumber,5)];
%         stressx = [stressx;statePressure(layernumber,6)];
%         stressy = [stressy;statePressure(layernumber,7)];
%         e1 = [e1;statePressure(layernumber,8)];
%         e2 = [e2;statePressure(layernumber,9)];
%         stress1 = [stress1;statePressure(layernumber,10)];
%         stress2 = [stress2;statePressure(layernumber,11)];
%         stressz = [stressz;statePressure(layernumber,12)];
%         shearxz = [shearxz;statePressure(layernumber,13)];
%     end
% 
%     vonMisesPressure = (0.5*((stress1-stress2).^2+(stress2).^2+(stress1).^2)).^(1/2);


    nTheta = 100;  % Number of angular divisions
    nZ = length(vonMisesAxial);  % Number of axial divisions based on the stress data size
    
    % Create the grid for the cylinder surface
    theta = linspace(0, 2*pi, nTheta);  % Angular positions
    z1 = linspace(0, L/2, nZ);           % Axial positions (half the length)
    z2 = linspace(L/2, L, nZ);
    z = [z1, z2];

    % Replicate stress values for the full cylinder
    StressAxialFull = [vonMisesAxial; flip(vonMisesAxial(:))];
    StressTorsionFull = [vonMisesTorsion; flip(vonMisesTorsion(:))];
    %StressPressureFull = [StressPressure; flip(vonMisesPressure(:))];
    
    hold on
    plot(z, StressAxialFull, 'r')
    %plot(z, StressTorsionFull, 'm')
    %plot(z, StressPressureFull, 'o')
      
    % Enhance plot appearance
    xlabel('X (mm)');
    ylabel('Stress (MPa)');
    title('Composite Tube Stress vs Axial Coordinate');
    legend('Axial', 'Torsion', 'Pressure')
    hold off

end