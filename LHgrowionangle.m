height = READdate(:,1);
Neh = READdate(:,2);
arf1h = READdate(:,3);
arf2h = READdate(:,4);
arf3h = READdate(:,5);
arf4h = READdate(:,6);
arf5h = READdate(:,7);
arf6h = READdate(:,8);
Teh = READdate(:,9);
O = READdate(:,10);
N2 = READdate(:,11);
O2 = READdate(:,12);
Tih = READdate(:,13);
e = 1.6022 .* 10.^(-19);
me = 9.10956 .* 10.^(-31);
kb = 1.38 .* 10.^(-23);
mii = 16 .* me .* 1837;
E0 = 8.854 .* 10.^(-12);
Wpeh = (Neh .* e.^2 ./ (E0 .* me)).^(0.5);
Wpih = Wpeh ./ (1837 .* 16).^0.5;
Wce = 2 .* pi .* 1.187 .* 10.^6;
Wci = 2 .* pi .* 41.2;
W0 = 2 .* pi .* 95000;
beta = 4.2 .* 10.^(-6);
c = 3 .* 10.^8;
uc = 4 .* 10.^3;
vt = 0.1 .* 10.^3;
mi = 16;
mr = 18;
WLH = (Wpih.^2 .* Wce.^2 ./ (Wpeh.^2 + Wce.^2)).^0.5;
ven = 9.32 .* 10.^(-12) .* N2 .* (1 - 3.44 .* 10.^(-5) .* Teh) + 1.21 .* 10.^(-10) .* O2 .* (1 + 2.15 .* 10.^(-12) .* Teh .^0.5) .* Teh + 5.49 .* 10.^(-10) .* O .* Teh.^(0.5);
vei = 54 .* Neh ./ Teh.^1.5 ./ 10.^6;
veh = ven + vei;
rLHh = 0.5 .* veh .* (1 - WLH.^2 ./ W0.^2 + WLH.^2 ./ (Wce .* Wci));

angle = [-pi : pi./ 200 : pi]' - atan(1 ./ 4);%角度只用改这
R = [0.4 .* 10.^3 : 0.4 .* 10.^3 : 600 .* 10.^3];
elevation = 0;
[xx,yy,zz] = sph2cart(elevation,angle,R);

xx_y = xx(:);
zz_y = zz(:);

t = 59.9;
sx = [-600 .* 10.^3 : 400 : 600 .* 10.^3];
De1 = Wpeh.^2 ./ Wce.^2;
De2 = Wpeh.^2 ./ W0.^2;
De2g = 2 .* Wpeh.^2 ./ W0.^3;

kCE = 3.2 .* 10.^(-15);
kDR = 3.6 .* 10.^(-13);
E0 = 8.854 .* 10.^(-12);
e = 1.6022 .* 10.^(-19);
me = 9.10956 .* 10.^(-31);
mring = 18 .* me .* 1837;
angle_0 = - atan(1 ./ 4);
detangle = pi ./ 6;
vs = 3.2 .* 10.^3;
vm = 0.1 .* 10.^3;
ND0 = 3.94 .* 10.^24;
te = 60;
angle = [-pi : pi./ 200 : pi]' .* ones(401,1500);
R = [0.4 .* 10.^3 : 0.4 .* 10.^3 : 600 .* 10.^3] .* ones(401,1500);
Nrz1 = ones(201,1001);
Nrz2 = ones(201,1001);
Nez = ones(201,1001);

% Nr_10km = ones(201,400);
% Nr_20km = ones(201,400);
% Nr_50km = ones(201,400);
% Nr_80km = ones(201,400);

% v2z = 1000 .* ones(401,3001);
% v1z = (-2 .* (kb .* Tih ./ mii .* gradient(Neh .* arf1h) ./ (Neh .* arf1h) + kb .* Teh ./ me .* gradient(Neh) ./ Neh)).^0.5 .* ones(1,3001);
ELH = 10.^(-6) .* ones(201,1001);
 parfor hang = 1 : 201
    Ne = Neh(100:300,:) .* ones(1,1001);
    Wpe = (Ne .* e.^2 ./ (E0 .* me)).^(0.5);
    De1 = Wpe.^2 ./ Wce.^2;
    De2 = Wpe.^2 ./ W0.^2;
    De2g = 2 .* Wpe.^2 ./ W0.^3;
    Wpi = Wpe ./ (1837 .* 16).^0.5;
    Wpe = Wpeh(hang,:);
    Wpi = Wpih(hang,:);
    arf = 0;
    arf0 = 0;
    arf1 = arf1h(100:300,:) .* ones(1,1001);
    arf2 = arf2h(hang+99,:);
    arf3 = arf3h(hang+99,:);
    arf4 = arf4h(hang+99,:);
    arf5 = arf5h(hang+99,:);
    arf6 = arf6h(hang+99,:);
    Te = Teh(100:300,:);
    Ti = Tih(100:300,:);
    rLH = rLHh(hang+99,:);
    Nr1 = ones(201,1001);
    Nr2 = ones(201,1001);
    % Nhangr2 = ones(201,1001);
    h = [-100 : 100]' .* 1500 .* ones(1,1001);
    % v2z1 = 1000 .* ones(201,1001);
    % v2z2 = 1000 .* ones(201,1001);
    % [Negx,Negz] = gradient(Ne);
    % [NOpgx,NOpgz] = gradient(Ne .* arf1);
    % v1z = real((-2 .* (kb .* Ti ./ mii .* NOpgz ./ (Ne .* arf1) + kb .* Te ./ me .* Negz ./ Ne)).^0.5 ./500 .* ones(1,1001));
    % v1z = zeros(201,1001);
for tt = 0 : 1 : 50
    fi_angle = (erf((angle + angle_0) ./ detangle) - erf((angle - angle_0) ./ detangle)) ./ (exp(-detangle.^2 ./ 4) .* cos(angle_0) - erf(angle_0 ./ detangle));
    
    NH2O = ND0 .* fi_angle .* (erf((R - (60 + 0.1 .* tt) .* vs) .* vs ./ ((2).^0.5 .* R .* vm)) - erf((R - (60 + 0.1 .* tt - min(60 + 0.1 .* tt,te)) .* vs) .* vs ./ ((2).^0.5 .* R .*vm))) ./ (8 .* pi .* R.^2 .* vs) .* 4 ./ 7;

    sx = [-200 .* 10.^3 + 4000 .* (0.1 .* tt) : 400 : 200 .* 10.^3 + 4000 .* (0.1 .* tt)] .* ones(201,1); %方向为负
    NH2O_y = NH2O(:);
%     xx_yy = xx_y - 710 .* tt;%方向为负
    F = scatteredInterpolant(xx_y,zz_y,NH2O_y);
    NH2Obf = F(sx,h);
    Nr1 = [Nr1(:,2:1001) ones(201,1)];
    Nr2 = [Nr2(:,2:1001) ones(201,1)];
    Ne = [Ne(:,2:1001) Neh(100:300,:)];
%     Nr1 = [ones(201,1) Nr1(:,1:1000)];
%     Nr2 = [ones(201,1) Nr2(:,1:1000)];
%     Ne = [Neh(100:300,:) Ne(:,1:1000)];
    % [v1zgx,v1zgz] = gradient(v1z.^2);
    % [v2zgx1,v2zgz1] = gradient(v2z1.^2);
    % [v2zgx2,v2zgz2] = gradient(v2z2.^2);
    % [Negx,Negz] = gradient(Ne);
    % [NOpgx,NOpgz] = gradient(Ne .* arf1);
    % [Nrgx1,Nrgz1] = gradient(Nr1);
    % [Nrgx2,Nrgz2] = gradient(Nr2);
    
    % v1z = v1z + 0.1 .* ((-0.5 .* v1zgz - kb .* Ti ./ mii .* NOpgz ./ (Ne .* arf1) ./ 250000 - kb .* Te ./ me .* Negz ./ Ne ./ 250000) ./ 1500 - 9.7 .* 10.^(-14) .* Nr1 .* (v1z - v2z1) - 9.7 .* 10.^(-14) .* Nr2 .* (v1z - v2z2));
    % v1z = v1z + 0.1 .* (-9.7 .* 10.^(-14) .* (Nr1 + Nr2) .* (v1z - v2z));
    % v2z = v2z + 0.1 .* ((-0.5 .* v2zgz - 100.^2 .* Nrgz ./ (Nr1+Nr2) - kb .* Te ./ me .* Negz ./ Ne) ./ 1500 - 9.7 .* 10.^(-14) .* (Ne .* arf1) .* (v2z - v1z));
    % v2z1 = v2z1 + 0.1 .* ((-0.5 .* v2zgz1 - 100.^2 .* Nrgz1 ./ Nr1 ./ 250000 - kb .* Te ./ me .* Negz ./ Ne ./ 250000) ./ 1500 - 9.7 .* 10.^(-14) .* (Ne .* arf1) .* (v2z1 - v1z));
    % v2z2 = v2z2 + 0.1 .* ((-0.5 .* v2zgz2 - 100.^2 .* Nrgz2 ./ Nr2 ./ 250000 - kb .* Te ./ me .* Negz ./ Ne ./ 250000) ./ 1500 - 9.7 .* 10.^(-14) .* (Ne .* arf1) .* (v2z2 - v1z));
%     Nr_y = Nr(:);
%     F = scatteredInterpolant(xx_y,zz_y,Nr_y);
%     Nr = F(sx,h);
    Nrmax = kCE .* arf1 .* Ne .* NH2Obf ./ (kDR .* Ne + kCE .* NH2Obf);
    % [av1zgx,av1zgz] = gradient(arf1 .* v1z);
    % [Nrv21zgx,Nrv21zgz] = gradient(Nr1 .* v2z1);
    % [Nrv22zgx,Nrv22zgz] = gradient(Nr2 .* v2z2);
    % Nrv21zgz = diff(Nr1);
    % Nrv22zgz = diff(Nr2);
    % detNr1 = kCE .* arf1 .* Ne .* NH2Obf .* 0.1 - kDR .* Nr1 .* Ne .* 0.1 - Nrv21zgz ./ 1500 .* 100;
    % detNr2 = -kDR .* Nr2 .* Ne .* 0.1 - Nrv22zgz ./ 1500 .* 100;

    %det1
    detNr1 = kCE .* arf1 .* Ne .* NH2Obf .* 0.1 - kDR .* Nr1 .* Ne .* 0.1 - Nr1 ./ 1500 .* 100 + [Nr1(2:201,:);0.95.*Nr1(201,:)] ./ 1500 .* 100;
    detNr2 = -kDR .* Nr2 .* Ne .* 0.1 - Nr2 ./ 1500 .* 100 + [Nr2(2:201,:);0.95.*Nr2(201,:)] ./ 1500 .* 100;
   
    %det2
    % [Nrv21zgx,Nrv21zgz] = gradient(Nr1);
    % [Nrv22zgx,Nrv22zgz] = gradient(Nr2);
    % detNr1 = kCE .* arf1 .* Ne .* NH2Obf .* 0.1 - kDR .* Nr1 .* Ne .* 0.1 + Nrv21zgz ./ 1500 .* 100;
    % detNr2 = -kDR .* Nr2 .* Ne .* 0.1 + Nrv22zgz ./ 1500 .* 100;
   
    % detNr1 = kCE .* arf1 .* Ne .* NH2Obf .* 0.1 - kDR .* Nr1 .* Ne .* 0.1;
    % detNr2 =  -kDR .* Nr2 .* Ne .* 0.1;
    % detarf1 = -kCE .* arf1 .* NH2Obf .* 0.1 - av1zgz ./ 1500 .* 0.1;
    % detarf1 = -kCE .* arf1 .* NH2Obf .* 0.1;

    detNe = -kDR .* (Nr1 + Nr2) .* Ne .* 0.1;
    % detNr1(isnan(detNr1)) = 0;
    % detNr2(isnan(detNr2)) = 0;
    % detNe(isnan(detNe)) = 0;

    % detarf1(isnan(detarf1)) = 0;
    
    Nr1 = Nr1 + detNr1;
    Nr2 = Nr2 + detNr2;
    % Nr1(isnan(Nr1)) = 1;
    % Nr2(isnan(Nr2)) = 1;
    % Nr1(Nr1<1) = 1;
    % Nr2(Nr2<1) = 1;
    Ne = Ne + detNe;
    Ne(isnan(Ne)) = Neh(hang,:);
    % arf1 = arf1 + detarf1;
    % arf1(isnan(arf1)) = arf1h(hang,:);
    % arf1(arf1 < 0) = 0;

    Wpe = (Ne .* e.^2 ./ (E0 .* me)).^(0.5);
    De1 = Wpe.^2 ./ Wce.^2;
    De2 = Wpe.^2 ./ W0.^2;
    De2g = 2 .* Wpe.^2 ./ W0.^3;
    Wpi = Wpe ./ (1837 .* 16).^0.5;
    % Nr1(Nr1<1) = 1;
    % Nr2(Nr2<1) = 1;
%     Nr_y = Nr(:);
%     F = scatteredInterpolant(xx_y,zz_y,Nr_y);
%     Nrzbf = F(sx,h);

    for sxi = 1 : 1001

        if Nr1(hang,sxi) > Nrmax(hang,sxi)
            Nr1(hang,sxi) = Nrmax(hang,sxi);
        end
        % if Nr2(hang,sxi) > Nrmax(hang,sxi)
        %     Nr2(hang,sxi) = Nrmax(hang,sxi);
        % end
        if Nr1(hang,sxi) < 1
            Nr1(hang,sxi) = 1;
        end
        if Nr2(hang,sxi) < 1
            Nr2(hang,sxi) = 1;
        end
        % Nr1(isnan(Nr1)) = 1;
        arf0 = Nr2(hang,sxi) ./ Ne(hang,sxi);
        arf = Nr1(hang,sxi) ./ Ne(hang,sxi);
        Di1 = (8./9) .* arf0 .* Wpi(hang,sxi).^2 ./ W0.^2 + (arf1(hang,sxi) - arf - arf0) .* Wpi(hang,sxi).^2 ./ W0.^2 + (8./7) .* arf2 .* Wpi(hang,sxi).^2 ./ W0.^2 + 8 .* arf3 .* Wpi(hang,sxi).^2 ./ W0.^2 + 4 .* arf4 .* Wpi(hang,sxi).^2 ./ W0.^2 + 0.5 .* arf5 .* Wpi(hang,sxi).^2 ./ W0.^2 + 8./15 .* arf6 .* Wpi(hang,sxi).^2 ./ W0.^2;
        Di2 = arf .* (mi ./ mr) .* Wpi(hang,sxi).^2 .* uc.^1.5 ./ (4 .* W0.^2 .* vt.^1.5);
        Di1g = 2 .* ((8./9) .* arf0 .* Wpi(hang,sxi).^2 ./ W0.^3 + (arf1(hang,sxi) - arf - arf0) .* Wpi(hang,sxi).^2 ./ W0.^3 + (8./7) .* arf2 .* Wpi(hang,sxi).^2 ./ W0.^3 + 8 .* arf3 .* Wpi(hang,sxi).^2 ./ W0.^3 + 4 .* arf4 .* Wpi(hang,sxi).^2 ./ W0.^3 + 0.5 .* arf5 .* Wpi(hang,sxi).^2 ./ W0.^3 + 8./15 .* arf6 .* Wpi(hang,sxi).^2 ./ W0.^3);
        Di2g = 2 .* arf .* (mi ./ mr) .* Wpi(hang,sxi).^2 .* uc.^1.5 ./ (4 .* W0.^3 .* vt.^1.5);
            
        aa = De2g(hang,sxi) .* (De1(hang,sxi) + De2(hang,sxi));
        bb = -(De2g(hang,sxi) + Di1g + Di2g) .* (De1(hang,sxi) + De2(hang,sxi)) + De2g(hang,sxi) .* (1 - Di1 - Di2 - De2(hang,sxi));
        cc = -(Di1g + Di2g + De2g(hang,sxi)) .* (1 - Di1 - Di2 - De2(hang,sxi)) + Di2 .* Di2g;

        sintheta2jj = (-bb - (bb.^2 - 4 .* aa .* cc).^0.5) ./ (2 .* aa);
        sintheta = (sintheta2jj).^0.5;

        %r = -Di2 ./ ((1 - sintheta.^2) .* De2g + Di1g + Di2g);

        % detWLH0 = 1./4 .*  E0 .* W0 .* (De2g(hang,:) .* (1 - sintheta.^2) + Di1g + Di2g) .* ELH(hang,sxi).^2 .* (exp((Di2 ./ ((1 - sintheta.^2) .* De2g(hang,:) + Di1g + Di2g) - rLH) .* 2 .* 0.1) - 1);
        detWLH0 = 1./4 .*  E0 .* W0 .* (De2g(hang,sxi) .* (1 - sintheta.^2) + Di1g + Di2g) .* ELH(hang,sxi).^2 .* (exp((Di2 ./ ((1 - sintheta.^2) .* De2g(hang,sxi) + Di1g + Di2g) - rLH) .* 2 .* 0.1) - exp(-rLH .* 0.2));
        WrS = 0.5 .* mring .* uc.^2 .* Nr1(hang,sxi) .* 2.4 .* 10.^(-4);
%         WrS = 0.5 .* mring .* uc.^2 .* 10.^8 .* 4.8 .* 10.^(-6);
        if WrS >= detWLH0
            ELH(hang,sxi) = ELH(hang,sxi) .* exp((Di2 ./ ((1 - sintheta.^2) .* De2g(hang,sxi) + Di1g + Di2g) - rLH) .* 0.1);
            if ELH(hang,sxi) < 10.^(-6)
                ELH(hang,sxi) = 10.^(-6);
            end
                Nr2(hang,sxi) = Nr2(hang,sxi) + Nr1(hang,sxi) .* detWLH0 ./ WrS;
                Nr1(hang,sxi) = Nr1(hang,sxi) .* (1 - detWLH0 ./ WrS);
        else
            ELH(hang,sxi) = ((1./4 .* E0 .* W0 .* (De2g(hang,sxi) .* (1 - sintheta.^2) + Di1g + Di2g) .* ELH(hang,sxi).^2 .* exp(-rLH .* 0.2) + WrS) ./ (1./4 .* E0 .* W0 .* (De2g(hang,sxi) .* (1 - sintheta.^2) + Di1g + Di2g))).^0.5;
            Nr2(hang,sxi) = Nr2(hang,sxi) + Nr1(hang,sxi);
            Nr1(hang,sxi) = 1;
        end
       
    end
%     Nr_y = Nr(:); 
    % sx_y = sx(:);
    % h_y = [(-201 + hang) .* 1500 .* ones(1,3001)]';
    % Nr_y = Nr(:);
    % G = scatteredInterpolant(sx_y,h_y,Nr_y);
    % sx = [-600 .* 10.^3 + 4000 .* (60.1 + 0.1 .* tt - 60) : 400 : 600 .* 10.^3 + 4000 .* (60.1 + 0.1 .* tt - 60)]; 
    % Nrzbf = G(sxF,h);

    % Nr_10km(hang,tt) = Nr1(hang,526-tt) + Nr2(hang,526-tt);
    % Nr_20km(hang,tt) = Nr1(hang,551-tt) + Nr2(hang,551-tt)
    % Nr_50km(hang,tt) = Nr1(hang,626-tt) + Nr2(hang,626-tt)
    % Nr_80km(hang,tt) = Nr1(hang,701-tt) + Nr2(hang,701-tt)
end
    % Nr1(isnan(Nr1)) = 1;
    % Nr2(isnan(Nr2)) = 1;
    Nrz1(hang,:) = Nr1(hang,:);
    Nrz2(hang,:) = Nr2(hang,:);
    Nez(hang,:) = Ne(hang,:);
 
 end
 ELH(isnan(ELH)) = 10.^(-6);
 % Nrz1(isnan(Nrz1)) = 1;     


 % Nrz1(Nrz1<1) = 1;
 % Nrz2(isnan(Nrz2)) = 1;
 % Nrz2(Nrz2<1) = 1;

sx = [-200 .* 10.^3 + 4000 .* (90 - 60) : 400 : 200 .* 10.^3 + 4000 .* (90 - 60)]; 
pcolor(sx,height(100:300,:), ELH);
shading flat

ELHR2Dtu = 20 .* log10(ELH ./ 10.^(-6));
pcolor(sx,height(100:300,:), real(ELHR2Dtu));
shading flat

pcolor(sx,height(100:300,:), log10(Nrz1+Nrz2));
shading flat

ttt = [1 : 400];
pcolor(ttt,height(100:300,:), log10(Nr_50km));
shading flat


save ELH95K-60-65-fuhe-shuyun-angle-700km-diff-new-gradient-24-4.mat ELH
save Nr195K-60-65-fuhe-shuyun-angle-700km-new-gradient-24-4.mat Nrz1
save Nr295K-60-65-fuhe-shuyun-angle-700km-new-gradient-24-4.mat Nrz2
save Nez95K-60-65-fuhe-shuyun-angle-700km-new-gradient-24-4.mat Nez


max(max(ELHR2Dtu))
% ELH(184:201,:) = 10.^(-6);
%cz
h = [-100 : 100]' .* 1500 .* ones(201,1001);
sx = [-200 .* 10.^3 + 4000 .* (80 - 60) : 400 : 200 .* 10.^3 + 4000 .* (80 - 60)] .* ones(201,1001); 
hcz = [-1500 : 1500]' .* 100 .* ones(3001,4001);
sxcz = [-200 .* 10.^3 + 4000 .* (80 - 60) : 100 : 200 .* 10.^3 + 4000 .* (80 - 60)] .* ones(3001,4001); 
h = h(:);
sx = sx(:);
ELH = ELH(:);

F = scatteredInterpolant(sx,h,ELH);
ELHcz = F(sxcz,hcz);
ELHcz(ELHcz<10.^(-6)) = 10.^(-6);

pcolor(sxcz,hcz, ELHcz);
shading flat

ELHR2Dtu = 20 .* log10(ELHcz ./ 10.^(-6));
pcolor(sxcz,hcz, ELHR2Dtu);
shading flat

%25.2K
height = [READdate(:,1)]';
height2 = [400 : 0.3 : 1000];
READdate2 = interp1(height,READdate,height2);
Neh = READdate2(:,2);
arf1h = READdate2(:,3);
arf2h = READdate2(:,4);
arf3h = READdate2(:,5);
arf4h = READdate2(:,6);
arf5h = READdate2(:,7);
arf6h = READdate2(:,8);
Teh = READdate2(:,9);
O = READdate2(:,10);
N2 = READdate2(:,11);
O2 = READdate2(:,12);
Tih = READdate2(:,13);
e = 1.6022 .* 10.^(-19);
me = 9.10956 .* 10.^(-31);
kb = 1.38 .* 10.^(-23);
mii = 16 .* me .* 1837;
E0 = 8.854 .* 10.^(-12);
Wpeh = (Neh .* e.^2 ./ (E0 .* me)).^(0.5);
Wpih = Wpeh ./ (1837 .* 16).^0.5;
Wce = 2 .* pi .* 1.187 .* 10.^6;
Wci = 2 .* pi .* 41.2;
W0 = 2 .* pi .* 428000;
beta = 4.2 .* 10.^(-6);
c = 3 .* 10.^8;
uc = 4 .* 10.^3;
vt = 0.1 .* 10.^3;
mi = 16;
mr = 18;
WLH = (Wpih.^2 .* Wce.^2 ./ (Wpeh.^2 + Wce.^2)).^0.5;
ven = 9.32 .* 10.^(-12) .* N2 .* (1 - 3.44 .* 10.^(-5) .* Teh) + 1.21 .* 10.^(-10) .* O2 .* (1 + 2.15 .* 10.^(-12) .* Teh .^0.5) .* Teh + 5.49 .* 10.^(-10) .* O .* Teh.^(0.5);
vei = 54 .* Neh ./ Teh.^1.5 ./ 10.^6;
veh = ven + vei;
rLHh = 0.5 .* veh .* (1 - WLH.^2 ./ W0.^2 + WLH.^2 ./ (Wce .* Wci));

angle = [-pi : pi./ 200 : pi]' - atan(1 ./ 4);
R = [0.4 .* 10.^3 : 0.4 .* 10.^3 : 600 .* 10.^3];
elevation = 0;
[xx,yy,zz] = sph2cart(elevation,angle,R);

xx_y = xx(:);
zz_y = zz(:);

t = 59.9;
sx = [-600 .* 10.^3 : 400 : 600 .* 10.^3];
De1 = Wpeh.^2 ./ Wce.^2;
De2 = Wpeh.^2 ./ W0.^2;
De2g = 2 .* Wpeh.^2 ./ W0.^3;

kCE = 3.2 .* 10.^(-15);
kDR = 3.6 .* 10.^(-13);
E0 = 8.854 .* 10.^(-12);
e = 1.6022 .* 10.^(-19);
me = 9.10956 .* 10.^(-31);
mring = 18 .* me .* 1837;
angle_0 = - atan(1 ./ 4);
detangle = pi ./ 6;
vs = 3.2 .* 10.^3;
vm = 0.1 .* 10.^3;
ND0 = 3.94 .* 10.^24;
te = 60;
angle = [-pi : pi./ 200 : pi]' .* ones(401,1500);
R = [0.4 .* 10.^3 : 0.4 .* 10.^3 : 600 .* 10.^3] .* ones(401,1500);
Nrz1 = ones(201,1001);
Nrz2 = ones(201,1001);
Nez = ones(201,1001);
% v2z = 1000 .* ones(401,3001);
% v1z = (-2 .* (kb .* Tih ./ mii .* gradient(Neh .* arf1h) ./ (Neh .* arf1h) + kb .* Teh ./ me .* gradient(Neh) ./ Neh)).^0.5 .* ones(1,3001);
ELH = 10.^(-6) .* ones(201,1001);
 parfor hang = 1 : 201
    Ne = Neh(900:1100,:) .* ones(1,1001);
    Wpe = (Ne .* e.^2 ./ (E0 .* me)).^(0.5);
    De1 = Wpe.^2 ./ Wce.^2;
    De2 = Wpe.^2 ./ W0.^2;
    De2g = 2 .* Wpe.^2 ./ W0.^3;
    Wpi = Wpe ./ (1837 .* 16).^0.5;
    Wpe = Wpeh(hang,:);
    Wpi = Wpih(hang,:);
    arf = 0;
    arf0 = 0;
    arf1 = arf1h(900:1100,:) .* ones(1,1001);
    arf2 = arf2h(hang+899,:);
    arf3 = arf3h(hang+899,:);
    arf4 = arf4h(hang+899,:);
    arf5 = arf5h(hang+899,:);
    arf6 = arf6h(hang+899,:);
    Te = Teh(900:1100,:);
    Ti = Tih(900:1100,:);
    rLH = rLHh(hang+899,:);
    Nr1 = ones(201,1001);
    Nr2 = ones(201,1001);
    % Nhangr2 = ones(201,1001);
    h = [-100 : 100]' .* 300 .* ones(1,1001);
    % v2z1 = 1000 .* ones(201,1001);
    % v2z2 = 1000 .* ones(201,1001);
    % [Negx,Negz] = gradient(Ne);
    % [NOpgx,NOpgz] = gradient(Ne .* arf1);
    % v1z = real((-2 .* (kb .* Ti ./ mii .* NOpgz ./ (Ne .* arf1) + kb .* Te ./ me .* Negz ./ Ne)).^0.5 ./500 .* ones(1,1001));
    % v1z = zeros(201,1001);
for tt = 1 : 1 : 500
    fi_angle = (erf((angle + angle_0) ./ detangle) - erf((angle - angle_0) ./ detangle)) ./ (exp(-detangle.^2 ./ 4) .* cos(angle_0) - erf(angle_0 ./ detangle));
    
    NH2O = ND0 .* fi_angle .* (erf((R - (60 + 0.02 .* tt) .* vs) .* vs ./ ((2).^0.5 .* R .* vm)) - erf((R - (60 + 0.02 .* tt - min(60 + 0.02 .* tt,te)) .* vs) .* vs ./ ((2).^0.5 .* R .*vm))) ./ (8 .* pi .* R.^2 .* vs) .* 4 ./ 7;

    sx = [-40 .* 10.^3 + 4000 .* (65 + 0.02 .* tt - 65) : 80 : 40 .* 10.^3 + 4000 .* (65 + 0.02 .* tt - 65)] .* ones(201,1); 
    NH2O_y = NH2O(:);
    F = scatteredInterpolant(xx_y,zz_y,NH2O_y);
    NH2Obf = F(sx,h);
    Nr1 = [Nr1(:,2:1001) ones(201,1)];
    Nr2 = [Nr2(:,2:1001) ones(201,1)];
    Ne = [Ne(:,2:1001) Neh(900:1100,:)];
    % [v1zgx,v1zgz] = gradient(v1z.^2);
    % [v2zgx1,v2zgz1] = gradient(v2z1.^2);
    % [v2zgx2,v2zgz2] = gradient(v2z2.^2);
    % [Negx,Negz] = gradient(Ne);
    % [NOpgx,NOpgz] = gradient(Ne .* arf1);
    % [Nrgx1,Nrgz1] = gradient(Nr1);
    % [Nrgx2,Nrgz2] = gradient(Nr2);
    
    % v1z = v1z + 0.1 .* ((-0.5 .* v1zgz - kb .* Ti ./ mii .* NOpgz ./ (Ne .* arf1) ./ 250000 - kb .* Te ./ me .* Negz ./ Ne ./ 250000) ./ 1500 - 9.7 .* 10.^(-14) .* Nr1 .* (v1z - v2z1) - 9.7 .* 10.^(-14) .* Nr2 .* (v1z - v2z2));
    % v1z = v1z + 0.1 .* (-9.7 .* 10.^(-14) .* (Nr1 + Nr2) .* (v1z - v2z));
    % v2z = v2z + 0.1 .* ((-0.5 .* v2zgz - 100.^2 .* Nrgz ./ (Nr1+Nr2) - kb .* Te ./ me .* Negz ./ Ne) ./ 1500 - 9.7 .* 10.^(-14) .* (Ne .* arf1) .* (v2z - v1z));
    % v2z1 = v2z1 + 0.1 .* ((-0.5 .* v2zgz1 - 100.^2 .* Nrgz1 ./ Nr1 ./ 250000 - kb .* Te ./ me .* Negz ./ Ne ./ 250000) ./ 1500 - 9.7 .* 10.^(-14) .* (Ne .* arf1) .* (v2z1 - v1z));
    % v2z2 = v2z2 + 0.1 .* ((-0.5 .* v2zgz2 - 100.^2 .* Nrgz2 ./ Nr2 ./ 250000 - kb .* Te ./ me .* Negz ./ Ne ./ 250000) ./ 1500 - 9.7 .* 10.^(-14) .* (Ne .* arf1) .* (v2z2 - v1z));
%     Nr_y = Nr(:);
%     F = scatteredInterpolant(xx_y,zz_y,Nr_y);
%     Nr = F(sx,h);
    Nrmax = kCE .* arf1 .* Ne .* NH2Obf ./ (kDR .* Ne + kCE .* NH2Obf);
    % [av1zgx,av1zgz] = gradient(arf1 .* v1z);
    % [Nrv21zgx,Nrv21zgz] = gradient(Nr1 .* v2z1);
    % [Nrv22zgx,Nrv22zgz] = gradient(Nr2 .* v2z2);
    % Nrv21zgz = diff(Nr1);
    % Nrv22zgz = diff(Nr2);
    % detNr1 = kCE .* arf1 .* Ne .* NH2Obf .* 0.1 - kDR .* Nr1 .* Ne .* 0.1 - Nrv21zgz ./ 1500 .* 100;
    % detNr2 = -kDR .* Nr2 .* Ne .* 0.1 - Nrv22zgz ./ 1500 .* 100;

    %det1
    detNr1 = kCE .* arf1 .* Ne .* NH2Obf .* 0.02 - kDR .* Nr1 .* Ne .* 0.02 - Nr1 ./ 300 .* 20 + [Nr1(2:201,:);0.95.*Nr1(201,:)] ./ 300 .* 20;
    detNr2 = -kDR .* Nr2 .* Ne .* 0.02 - Nr2 ./ 300 .* 20 + [Nr2(2:201,:);0.95.*Nr2(201,:)] ./ 300 .* 20;
   
    %det2
    % [Nrv21zgx,Nrv21zgz] = gradient(Nr1);
    % [Nrv22zgx,Nrv22zgz] = gradient(Nr2);
    % detNr1 = kCE .* arf1 .* Ne .* NH2Obf .* 0.1 - kDR .* Nr1 .* Ne .* 0.1 + Nrv21zgz ./ 1500 .* 100;
    % detNr2 = -kDR .* Nr2 .* Ne .* 0.1 + Nrv22zgz ./ 1500 .* 100;
   
    % detNr1 = kCE .* arf1 .* Ne .* NH2Obf .* 0.1 - kDR .* Nr1 .* Ne .* 0.1;
    % detNr2 =  -kDR .* Nr2 .* Ne .* 0.1;
    % detarf1 = -kCE .* arf1 .* NH2Obf .* 0.1 - av1zgz ./ 1500 .* 0.1;
    % detarf1 = -kCE .* arf1 .* NH2Obf .* 0.1;

    detNe = -kDR .* (Nr1 + Nr2) .* Ne .* 0.02;
    % detNr1(isnan(detNr1)) = 0;
    % detNr2(isnan(detNr2)) = 0;
    % detNe(isnan(detNe)) = 0;

    % detarf1(isnan(detarf1)) = 0;
    
    Nr1 = Nr1 + detNr1;
    Nr2 = Nr2 + detNr2;
    % Nr1(isnan(Nr1)) = 1;
    % Nr2(isnan(Nr2)) = 1;
    % Nr1(Nr1<1) = 1;
    % Nr2(Nr2<1) = 1;
    Ne = Ne + detNe;
    Ne(isnan(Ne)) = Neh(hang,:);
    % arf1 = arf1 + detarf1;
    % arf1(isnan(arf1)) = arf1h(hang,:);
    % arf1(arf1 < 0) = 0;

    Wpe = (Ne .* e.^2 ./ (E0 .* me)).^(0.5);
    De1 = Wpe.^2 ./ Wce.^2;
    De2 = Wpe.^2 ./ W0.^2;
    De2g = 2 .* Wpe.^2 ./ W0.^3;
    Wpi = Wpe ./ (1837 .* 16).^0.5;
    % Nr1(Nr1<1) = 1;
    % Nr2(Nr2<1) = 1;
%     Nr_y = Nr(:);
%     F = scatteredInterpolant(xx_y,zz_y,Nr_y);
%     Nrzbf = F(sx,h);

    for sxi = 1 : 1001

        if Nr1(hang,sxi) > Nrmax(hang,sxi)
            Nr1(hang,sxi) = Nrmax(hang,sxi);
        end
        % if Nr2(hang,sxi) > Nrmax(hang,sxi)
        %     Nr2(hang,sxi) = Nrmax(hang,sxi);
        % end
        if Nr1(hang,sxi) < 1
            Nr1(hang,sxi) = 1;
        end
        if Nr2(hang,sxi) < 1
            Nr2(hang,sxi) = 1;
        end
        % Nr1(isnan(Nr1)) = 1;
        arf0 = Nr2(hang,sxi) ./ Ne(hang,sxi);
        arf = Nr1(hang,sxi) ./ Ne(hang,sxi);
        Di1 = (8./9) .* arf0 .* Wpi(hang,sxi).^2 ./ W0.^2 + (arf1(hang,sxi) - arf - arf0) .* Wpi(hang,sxi).^2 ./ W0.^2 + (8./7) .* arf2 .* Wpi(hang,sxi).^2 ./ W0.^2 + 8 .* arf3 .* Wpi(hang,sxi).^2 ./ W0.^2 + 4 .* arf4 .* Wpi(hang,sxi).^2 ./ W0.^2 + 0.5 .* arf5 .* Wpi(hang,sxi).^2 ./ W0.^2 + 8./15 .* arf6 .* Wpi(hang,sxi).^2 ./ W0.^2;
        Di2 = arf .* (mi ./ mr) .* Wpi(hang,sxi).^2 .* uc.^1.5 ./ (4 .* W0.^2 .* vt.^1.5);
        Di1g = 2 .* ((8./9) .* arf0 .* Wpi(hang,sxi).^2 ./ W0.^3 + (arf1(hang,sxi) - arf - arf0) .* Wpi(hang,sxi).^2 ./ W0.^3 + (8./7) .* arf2 .* Wpi(hang,sxi).^2 ./ W0.^3 + 8 .* arf3 .* Wpi(hang,sxi).^2 ./ W0.^3 + 4 .* arf4 .* Wpi(hang,sxi).^2 ./ W0.^3 + 0.5 .* arf5 .* Wpi(hang,sxi).^2 ./ W0.^3 + 8./15 .* arf6 .* Wpi(hang,sxi).^2 ./ W0.^3);
        Di2g = 2 .* arf .* (mi ./ mr) .* Wpi(hang,sxi).^2 .* uc.^1.5 ./ (4 .* W0.^3 .* vt.^1.5);
            
        aa = De2g(hang,sxi) .* (De1(hang,sxi) + De2(hang,sxi));
        bb = -(De2g(hang,sxi) + Di1g + Di2g) .* (De1(hang,sxi) + De2(hang,sxi)) + De2g(hang,sxi) .* (1 - Di1 - Di2 - De2(hang,sxi));
        cc = -(Di1g + Di2g + De2g(hang,sxi)) .* (1 - Di1 - Di2 - De2(hang,sxi)) + Di2 .* Di2g;

        sintheta2jj = (-bb - (bb.^2 - 4 .* aa .* cc).^0.5) ./ (2 .* aa);
        sintheta = (sintheta2jj).^0.5;

        %r = -Di2 ./ ((1 - sintheta.^2) .* De2g + Di1g + Di2g);

        % detWLH0 = 1./4 .*  E0 .* W0 .* (De2g(hang,:) .* (1 - sintheta.^2) + Di1g + Di2g) .* ELH(hang,sxi).^2 .* (exp((Di2 ./ ((1 - sintheta.^2) .* De2g(hang,:) + Di1g + Di2g) - rLH) .* 2 .* 0.1) - 1);
        detWLH0 = 1./4 .*  E0 .* W0 .* (De2g(hang,sxi) .* (1 - sintheta.^2) + Di1g + Di2g) .* ELH(hang,sxi).^2 .* (exp((Di2 ./ ((1 - sintheta.^2) .* De2g(hang,sxi) + Di1g + Di2g) - rLH) .* 2 .* 0.02) - exp(-rLH .* 0.04));
        WrS = 0.5 .* mring .* uc.^2 .* Nr1(hang,sxi) .* 2.4 .* 10.^(-4);
%         WrS = 0.5 .* mring .* uc.^2 .* 10.^8 .* 4.8 .* 10.^(-6);
        if WrS >= detWLH0
            ELH(hang,sxi) = ELH(hang,sxi) .* exp((Di2 ./ ((1 - sintheta.^2) .* De2g(hang,sxi) + Di1g + Di2g) - rLH) .* 0.02);
            if ELH(hang,sxi) < 10.^(-6)
                ELH(hang,sxi) = 10.^(-6);
            end
                Nr2(hang,sxi) = Nr2(hang,sxi) + Nr1(hang,sxi) .* detWLH0 ./ WrS;
                Nr1(hang,sxi) = Nr1(hang,sxi) .* (1 - detWLH0 ./ WrS);
        else
            ELH(hang,sxi) = ((1./4 .* E0 .* W0 .* (De2g(hang,sxi) .* (1 - sintheta.^2) + Di1g + Di2g) .* ELH(hang,sxi).^2 .* exp(-rLH .* 0.04) + WrS) ./ (1./4 .* E0 .* W0 .* (De2g(hang,sxi) .* (1 - sintheta.^2) + Di1g + Di2g))).^0.5;
            Nr2(hang,sxi) = Nr2(hang,sxi) + Nr1(hang,sxi);
            Nr1(hang,sxi) = 1;
        end
       
    end
%     Nr_y = Nr(:); 
    % sx_y = sx(:);
    % h_y = [(-201 + hang) .* 1500 .* ones(1,3001)]';
    % Nr_y = Nr(:);
    % G = scatteredInterpolant(sx_y,h_y,Nr_y);
    % sx = [-600 .* 10.^3 + 4000 .* (60.1 + 0.1 .* tt - 60) : 400 : 600 .* 10.^3 + 4000 .* (60.1 + 0.1 .* tt - 60)]; 
    % Nrzbf = G(sxF,h);
end
    % Nr1(isnan(Nr1)) = 1;
    % Nr2(isnan(Nr2)) = 1;
    Nrz1(hang,:) = Nr1(hang,:);
    Nrz2(hang,:) = Nr2(hang,:);
    Nez(hang,:) = Ne(hang,:);
 end
 ELH(isnan(ELH)) = 10.^(-6);
 % Nrz1(isnan(Nrz1)) = 1;     


 % Nrz1(Nrz1<1) = 1;
 % Nrz2(isnan(Nrz2)) = 1;
 % Nrz2(Nrz2<1) = 1;

sx = [-40 .* 10.^3 + 4000 .* (100 - 60) : 80 : 40 .* 10.^3 + 4000 .* (100 - 60)]; 
pcolor(sx,[height2(:,900:1100)]', ELH);
shading flat

ELHR2Dtu = 20 .* log10(ELH ./ 10.^(-6));
pcolor(sx,[height2(:,900:1100)]', real(ELHR2Dtu));
shading flat



pcolor(sx,height(100:300,:), log10(Nrz1+Nrz2));
shading flat


save ELH428K-60-70-fuhe-shuyun-angle-700km-diff-new-gradient-24-4.mat ELH
save Nr1428K-60-70-fuhe-shuyun-angle-700km-new-gradient-24-4.mat Nrz1
save Nr2428K-60-70-fuhe-shuyun-angle-700km-new-gradient-24-4.mat Nrz2
save Nez428K-60-70-fuhe-shuyun-angle-700km-new-gradient-24-4.mat Nez

max(max(ELHR2Dtu))


%<60s
height = [READdate(:,1)]';
height2 = [180 : 0.2 : 780];
READdate2 = interp1(height,READdate,height2);
Neh = READdate2(:,2);
arf1h = READdate2(:,3);
arf2h = READdate2(:,4);
arf3h = READdate2(:,5);
arf4h = READdate2(:,6);
arf5h = READdate2(:,7);
arf6h = READdate2(:,8);
Teh = READdate2(:,9);
O = READdate2(:,10);
N2 = READdate2(:,11);
O2 = READdate2(:,12);
Tih = READdate2(:,13);
e = 1.6022 .* 10.^(-19);
me = 9.10956 .* 10.^(-31);
kb = 1.38 .* 10.^(-23);
mii = 16 .* me .* 1837;
E0 = 8.854 .* 10.^(-12);
Wpeh = (Neh .* e.^2 ./ (E0 .* me)).^(0.5);
Wpih = Wpeh ./ (1837 .* 16).^0.5;
Wce = 2 .* pi .* 1.187 .* 10.^6;
Wci = 2 .* pi .* 41.2;
W0 = 2 .* pi .* 30000;
beta = 4.2 .* 10.^(-6);
c = 3 .* 10.^8;
uc = 4 .* 10.^3;
vt = 0.1 .* 10.^3;
mi = 16;
mr = 18;
WLH = (Wpih.^2 .* Wce.^2 ./ (Wpeh.^2 + Wce.^2)).^0.5;
ven = 9.32 .* 10.^(-12) .* N2 .* (1 - 3.44 .* 10.^(-5) .* Teh) + 1.21 .* 10.^(-10) .* O2 .* (1 + 2.15 .* 10.^(-12) .* Teh .^0.5) .* Teh + 5.49 .* 10.^(-10) .* O .* Teh.^(0.5);
vei = 54 .* Neh ./ Teh.^1.5 ./ 10.^6;
veh = ven + vei;
rLHh = 0.5 .* veh .* (1 - WLH.^2 ./ W0.^2 + WLH.^2 ./ (Wce .* Wci));

angle = [-pi : pi./ 200 : pi]' - atan(1 ./ 4);
R = [0.4 .* 10.^3 : 0.4 .* 10.^3 : 600 .* 10.^3];
elevation = 0;
[xx,yy,zz] = sph2cart(elevation,angle,R);

xx_y = xx(:);
zz_y = zz(:);

t = 59.9;
sx = [-600 .* 10.^3 : 400 : 600 .* 10.^3];
De1 = Wpeh.^2 ./ Wce.^2;
De2 = Wpeh.^2 ./ W0.^2;
De2g = 2 .* Wpeh.^2 ./ W0.^3;

kCE = 3.2 .* 10.^(-15);
kDR = 3.6 .* 10.^(-13);
E0 = 8.854 .* 10.^(-12);
e = 1.6022 .* 10.^(-19);
me = 9.10956 .* 10.^(-31);
mring = 18 .* me .* 1837;
angle_0 = - atan(1 ./ 4);
detangle = pi ./ 6;
vs = 3.2 .* 10.^3;
vm = 0.1 .* 10.^3;
ND0 = 3.94 .* 10.^24;
te = 60;
angle = [-pi : pi./ 200 : pi]' .* ones(401,1500);
R = [0.4 .* 10.^3 : 0.4 .* 10.^3 : 600 .* 10.^3] .* ones(401,1500);
Nrz1 = ones(201,1001);
Nrz2 = ones(201,1001);
Nez = ones(201,1001);
% v2z = 1000 .* ones(401,3001);
% v1z = (-2 .* (kb .* Tih ./ mii .* gradient(Neh .* arf1h) ./ (Neh .* arf1h) + kb .* Teh ./ me .* gradient(Neh) ./ Neh)).^0.5 .* ones(1,3001);
ELH = 10.^(-6) .* ones(201,1001);
 parfor hang = 1 : 201
    Ne = Neh(900:1100,:) .* ones(1,1001);
    Wpe = (Ne .* e.^2 ./ (E0 .* me)).^(0.5);
    De1 = Wpe.^2 ./ Wce.^2;
    De2 = Wpe.^2 ./ W0.^2;
    De2g = 2 .* Wpe.^2 ./ W0.^3;
    Wpi = Wpe ./ (1837 .* 16).^0.5;
    Wpe = Wpeh(hang,:);
    Wpi = Wpih(hang,:);
    arf = 0;
    arf0 = 0;
    arf1 = arf1h(1400:1600,:) .* ones(1,1001);
    arf2 = arf2h(hang+1399,:);
    arf3 = arf3h(hang+1399,:);
    arf4 = arf4h(hang+1399,:);
    arf5 = arf5h(hang+1399,:);
    arf6 = arf6h(hang+1399,:);
    Te = Teh(1400:1600,:);
    Ti = Tih(1400:1600,:);
    rLH = rLHh(hang+1399,:);
    Nr1 = ones(201,1001);
    Nr2 = ones(201,1001);
    % Nhangr2 = ones(201,1001);
    h = [-100 : 100]' .* 200 .* ones(1,1001);
    % v2z1 = 1000 .* ones(201,1001);
    % v2z2 = 1000 .* ones(201,1001);
    % [Negx,Negz] = gradient(Ne);
    % [NOpgx,NOpgz] = gradient(Ne .* arf1);
    % v1z = real((-2 .* (kb .* Ti ./ mii .* NOpgz ./ (Ne .* arf1) + kb .* Te ./ me .* Negz ./ Ne)).^0.5 ./500 .* ones(1,1001));
    % v1z = zeros(201,1001);
for tt = 0 : 1 : 200
    fi_angle = (erf((angle + angle_0) ./ detangle) - erf((angle - angle_0) ./ detangle)) ./ (exp(-detangle.^2 ./ 4) .* cos(angle_0) - erf(angle_0 ./ detangle));
    
    NH2O = ND0 .* fi_angle .* (erf((R - (28 + 0.01 .* tt) .* vs) .* vs ./ ((2).^0.5 .* R .* vm)) - erf((R - (28 + 0.01 .* tt - min(28 + 0.01 .* tt,te)) .* vs) .* vs ./ ((2).^0.5 .* R .*vm))) ./ (8 .* pi .* R.^2 .* vs) .* 4 ./ 7;

    sx = [-20 .* 10.^3 + 4000 .* (0.01 .* tt) : 40 : 20 .* 10.^3 + 4000 .* (0.01 .* tt)] .* ones(201,1); 
    NH2O_y = NH2O(:);
    F = scatteredInterpolant(xx_y,zz_y,NH2O_y);
    NH2Obf = F(sx,h);
    Nr1 = [Nr1(:,2:1001) ones(201,1)];
    Nr2 = [Nr2(:,2:1001) ones(201,1)];
    Ne = [Ne(:,2:1001) Neh(1400:1600,:)];
    
    % [v1zgx,v1zgz] = gradient(v1z.^2);
    % [v2zgx1,v2zgz1] = gradient(v2z1.^2);
    % [v2zgx2,v2zgz2] = gradient(v2z2.^2);
    % [Negx,Negz] = gradient(Ne);
    % [NOpgx,NOpgz] = gradient(Ne .* arf1);
    % [Nrgx1,Nrgz1] = gradient(Nr1);
    % [Nrgx2,Nrgz2] = gradient(Nr2);
    
    % v1z = v1z + 0.1 .* ((-0.5 .* v1zgz - kb .* Ti ./ mii .* NOpgz ./ (Ne .* arf1) ./ 250000 - kb .* Te ./ me .* Negz ./ Ne ./ 250000) ./ 1500 - 9.7 .* 10.^(-14) .* Nr1 .* (v1z - v2z1) - 9.7 .* 10.^(-14) .* Nr2 .* (v1z - v2z2));
    % v1z = v1z + 0.1 .* (-9.7 .* 10.^(-14) .* (Nr1 + Nr2) .* (v1z - v2z));
    % v2z = v2z + 0.1 .* ((-0.5 .* v2zgz - 100.^2 .* Nrgz ./ (Nr1+Nr2) - kb .* Te ./ me .* Negz ./ Ne) ./ 1500 - 9.7 .* 10.^(-14) .* (Ne .* arf1) .* (v2z - v1z));
    % v2z1 = v2z1 + 0.1 .* ((-0.5 .* v2zgz1 - 100.^2 .* Nrgz1 ./ Nr1 ./ 250000 - kb .* Te ./ me .* Negz ./ Ne ./ 250000) ./ 1500 - 9.7 .* 10.^(-14) .* (Ne .* arf1) .* (v2z1 - v1z));
    % v2z2 = v2z2 + 0.1 .* ((-0.5 .* v2zgz2 - 100.^2 .* Nrgz2 ./ Nr2 ./ 250000 - kb .* Te ./ me .* Negz ./ Ne ./ 250000) ./ 1500 - 9.7 .* 10.^(-14) .* (Ne .* arf1) .* (v2z2 - v1z));
%     Nr_y = Nr(:);
%     F = scatteredInterpolant(xx_y,zz_y,Nr_y);
%     Nr = F(sx,h);
    Nrmax = kCE .* arf1 .* Ne .* NH2Obf ./ (kDR .* Ne + kCE .* NH2Obf);
    % [av1zgx,av1zgz] = gradient(arf1 .* v1z);
    % [Nrv21zgx,Nrv21zgz] = gradient(Nr1 .* v2z1);
    % [Nrv22zgx,Nrv22zgz] = gradient(Nr2 .* v2z2);
    % Nrv21zgz = diff(Nr1);
    % Nrv22zgz = diff(Nr2);
    % detNr1 = kCE .* arf1 .* Ne .* NH2Obf .* 0.1 - kDR .* Nr1 .* Ne .* 0.1 - Nrv21zgz ./ 1500 .* 100;
    % detNr2 = -kDR .* Nr2 .* Ne .* 0.1 - Nrv22zgz ./ 1500 .* 100;

    %det1
    detNr1 = kCE .* arf1 .* Ne .* NH2Obf .* 0.02 - kDR .* Nr1 .* Ne .* 0.02 - Nr1 ./ 300 .* 20 + [Nr1(2:201,:);0.95.*Nr1(201,:)] ./ 300 .* 20;
    detNr2 = -kDR .* Nr2 .* Ne .* 0.02 - Nr2 ./ 300 .* 20 + [Nr2(2:201,:);0.95.*Nr2(201,:)] ./ 300 .* 20;
   
    %det2
    % [Nrv21zgx,Nrv21zgz] = gradient(Nr1);
    % [Nrv22zgx,Nrv22zgz] = gradient(Nr2);
    % detNr1 = kCE .* arf1 .* Ne .* NH2Obf .* 0.1 - kDR .* Nr1 .* Ne .* 0.1 + Nrv21zgz ./ 1500 .* 100;
    % detNr2 = -kDR .* Nr2 .* Ne .* 0.1 + Nrv22zgz ./ 1500 .* 100;
   
    % detNr1 = kCE .* arf1 .* Ne .* NH2Obf .* 0.1 - kDR .* Nr1 .* Ne .* 0.1;
    % detNr2 =  -kDR .* Nr2 .* Ne .* 0.1;
    % detarf1 = -kCE .* arf1 .* NH2Obf .* 0.1 - av1zgz ./ 1500 .* 0.1;
    % detarf1 = -kCE .* arf1 .* NH2Obf .* 0.1;

    detNe = -kDR .* (Nr1 + Nr2) .* Ne .* 0.02;
    % detNr1(isnan(detNr1)) = 0;
    % detNr2(isnan(detNr2)) = 0;
    % detNe(isnan(detNe)) = 0;

    % detarf1(isnan(detarf1)) = 0;
    
    Nr1 = Nr1 + detNr1;
    Nr2 = Nr2 + detNr2;
    % Nr1(isnan(Nr1)) = 1;
    % Nr2(isnan(Nr2)) = 1;
    % Nr1(Nr1<1) = 1;
    % Nr2(Nr2<1) = 1;
    Ne = Ne + detNe;
    Ne(isnan(Ne)) = Neh(hang,:);
    arf1 = arf1 - detNr1 ./ Ne - detNr2 ./ Ne;
    % arf1 = arf1 + detarf1;
    % arf1(isnan(arf1)) = arf1h(hang,:);
    % arf1(arf1 < 0) = 0;

    Wpe = (Ne .* e.^2 ./ (E0 .* me)).^(0.5);
    De1 = Wpe.^2 ./ Wce.^2;
    De2 = Wpe.^2 ./ W0.^2;
    De2g = 2 .* Wpe.^2 ./ W0.^3;
    Wpi = Wpe ./ (1837 .* 16).^0.5;
    % Nr1(Nr1<1) = 1;
    % Nr2(Nr2<1) = 1;
%     Nr_y = Nr(:);
%     F = scatteredInterpolant(xx_y,zz_y,Nr_y);
%     Nrzbf = F(sx,h);

    for sxi = 1 : 1001

        if Nr1(hang,sxi) > Nrmax(hang,sxi)
            Nr1(hang,sxi) = Nrmax(hang,sxi);
        end
        % if Nr2(hang,sxi) > Nrmax(hang,sxi)
        %     Nr2(hang,sxi) = Nrmax(hang,sxi);
        % end
        if Nr1(hang,sxi) < 1
            Nr1(hang,sxi) = 1;
        end
        if Nr2(hang,sxi) < 1
            Nr2(hang,sxi) = 1;
        end
        % Nr1(isnan(Nr1)) = 1;
        arf0 = Nr2(hang,sxi) ./ Ne(hang,sxi);
        arf = Nr1(hang,sxi) ./ Ne(hang,sxi);
        Di1 = (8./9) .* arf0 .* Wpi(hang,sxi).^2 ./ W0.^2 + (arf1(hang,sxi)) .* Wpi(hang,sxi).^2 ./ W0.^2 + (8./7) .* arf2 .* Wpi(hang,sxi).^2 ./ W0.^2 + 8 .* arf3 .* Wpi(hang,sxi).^2 ./ W0.^2 + 4 .* arf4 .* Wpi(hang,sxi).^2 ./ W0.^2 + 0.5 .* arf5 .* Wpi(hang,sxi).^2 ./ W0.^2 + 8./15 .* arf6 .* Wpi(hang,sxi).^2 ./ W0.^2;
        Di2 = arf .* (mi ./ mr) .* Wpi(hang,sxi).^2 .* uc.^1.5 ./ (4 .* W0.^2 .* vt.^1.5);
        Di1g = 2 .* ((8./9) .* arf0 .* Wpi(hang,sxi).^2 ./ W0.^3 + (arf1(haang,sxi)) .* Wpi(hang,sxi).^2 ./ W0.^3 + (8./7) .* arf2 .* Wpi(hang,sxi).^2 ./ W0.^3 + 8 .* arf3 .* Wpi(hang,sxi).^2 ./ W0.^3 + 4 .* arf4 .* Wpi(hang,sxi).^2 ./ W0.^3 + 0.5 .* arf5 .* Wpi(hang,sxi).^2 ./ W0.^3 + 8./15 .* arf6 .* Wpi(hang,sxi).^2 ./ W0.^3);
        Di2g = 2 .* arf .* (mi ./ mr) .* Wpi(hang,sxi).^2 .* uc.^1.5 ./ (4 .* W0.^3 .* vt.^1.5);
            
        aa = De2g(hang,sxi) .* (De1(hang,sxi) + De2(hang,sxi));
        bb = -(De2g(hang,sxi) + Di1g + Di2g) .* (De1(hang,sxi) + De2(hang,sxi)) + De2g(hang,sxi) .* (1 - Di1 - Di2 - De2(hang,sxi));
        cc = -(Di1g + Di2g + De2g(hang,sxi)) .* (1 - Di1 - Di2 - De2(hang,sxi)) + Di2 .* Di2g;

        sintheta2jj = (-bb - (bb.^2 - 4 .* aa .* cc).^0.5) ./ (2 .* aa);
        sintheta = (sintheta2jj).^0.5;

        %r = -Di2 ./ ((1 - sintheta.^2) .* De2g + Di1g + Di2g);

        % detWLH0 = 1./4 .*  E0 .* W0 .* (De2g(hang,:) .* (1 - sintheta.^2) + Di1g + Di2g) .* ELH(hang,sxi).^2 .* (exp((Di2 ./ ((1 - sintheta.^2) .* De2g(hang,:) + Di1g + Di2g) - rLH) .* 2 .* 0.1) - 1);
        detWLH0 = 1./4 .*  E0 .* W0 .* (De2g(hang,sxi) .* (1 - sintheta.^2) + Di1g + Di2g) .* ELH(hang,sxi).^2 .* (exp((Di2 ./ ((1 - sintheta.^2) .* De2g(hang,sxi) + Di1g + Di2g) - rLH) .* 2 .* 0.01) - exp(-rLH .* 0.02));
        WrS = 0.5 .* mring .* uc.^2 .* Nr1(hang,sxi) .* 2.4 .* 10.^(-4);
%         WrS = 0.5 .* mring .* uc.^2 .* 10.^8 .* 4.8 .* 10.^(-6);
        if WrS >= detWLH0
            ELH(hang,sxi) = ELH(hang,sxi) .* exp((Di2 ./ ((1 - sintheta.^2) .* De2g(hang,sxi) + Di1g + Di2g) - rLH) .* 0.01);
            if ELH(hang,sxi) < 10.^(-6)
                ELH(hang,sxi) = 10.^(-6);
            end
                Nr2(hang,sxi) = Nr2(hang,sxi) + Nr1(hang,sxi) .* detWLH0 ./ WrS;
                Nr1(hang,sxi) = Nr1(hang,sxi) .* (1 - detWLH0 ./ WrS);
        else
            ELH(hang,sxi) = ((1./4 .* E0 .* W0 .* (De2g(hang,sxi) .* (1 - sintheta.^2) + Di1g + Di2g) .* ELH(hang,sxi).^2 .* exp(-rLH .* 0.02) + WrS) ./ (1./4 .* E0 .* W0 .* (De2g(hang,sxi) .* (1 - sintheta.^2) + Di1g + Di2g))).^0.5;
            Nr2(hang,sxi) = Nr2(hang,sxi) + Nr1(hang,sxi);
            Nr1(hang,sxi) = 1;
        end
       
    end
%     Nr_y = Nr(:); 
    % sx_y = sx(:);
    % h_y = [(-201 + hang) .* 1500 .* ones(1,3001)]';
    % Nr_y = Nr(:);
    % G = scatteredInterpolant(sx_y,h_y,Nr_y);
    % sx = [-600 .* 10.^3 + 4000 .* (60.1 + 0.1 .* tt - 60) : 400 : 600 .* 10.^3 + 4000 .* (60.1 + 0.1 .* tt - 60)]; 
    % Nrzbf = G(sxF,h);
end
    % Nr1(isnan(Nr1)) = 1;
    % Nr2(isnan(Nr2)) = 1;
    Nrz1(hang,:) = Nr1(hang,:);
    Nrz2(hang,:) = Nr2(hang,:);
    Nez(hang,:) = Ne(hang,:);
 end
 ELH(isnan(ELH)) = 10.^(-6);
 % Nrz1(isnan(Nrz1)) = 1;     


 % Nrz1(Nrz1<1) = 1;
 % Nrz2(isnan(Nrz2)) = 1;
 % Nrz2(Nrz2<1) = 1;

sx = [-40 .* 10.^3 + 4000 .* (100 - 60) : 80 : 40 .* 10.^3 + 4000 .* (100 - 60)]; 
pcolor(sx,[height2(:,1400:1600)]', ELH);
shading flat

ELHR2Dtu = 20 .* log10(ELH ./ 10.^(-6));
pcolor(sx,[height2(:,1400:1600)]', real(ELHR2Dtu));
shading flat



pcolor(sx,height(100:300,:), log10(Nrz1+Nrz2));
shading flat


save ELH30K-28-30-fuhe-shuyun-angle-480km-diff-new-gradient-24-4.mat ELH
save Nr130K-28-30-fuhe-shuyun-angle-480km-new-gradient-24-4.mat Nrz1
save Nr230K-28-30-fuhe-shuyun-angle-480km-new-gradient-24-4.mat Nrz2
save Nez30K-28-30-fuhe-shuyun-angle-480km-new-gradient-24-4.mat Nez

max(max(ELHR2Dtu))