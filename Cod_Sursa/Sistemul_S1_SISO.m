clear; clc; close all;

set(groot,'defaultTextInterpreter','none');
set(groot,'defaultAxesTickLabelInterpreter','none');
set(groot,'defaultLegendInterpreter','none');

rng(42,'twister');

k1 = 1.00;
k2 = 0.15;
b0 = 1.00;

n = 2;
m = 1;

f = @(x) [x(2);
         -k1*x(1) - k2*x(1)^3];

g = @(x) [0;
          b0];

x_echilibru = zeros(n,1);

u_echilibru = (k1*x_echilibru(1) + k2*x_echilibru(1)^3)/b0;

tau_exact   = @(x) [x(1) - x_echilibru(1);
                    x(2) - x_echilibru(2)];

delta_exact = @(x) -k1*x(1) - k2*x(1)^3;
gamma_exact = @(x) b0;

echilibru_reziduu = norm(f(x_echilibru) + g(x_echilibru)*u_echilibru);

r1 = 2;
grad_relativ_complet = (r1 == n);

Ac = [0 1;
      0 0];

Bc = [0;
      1];

Q_lqr = diag([20 2]);
R_lqr = 0.5;
K_feedback = -lqr(Ac,Bc,Q_lqr,R_lqr);

poli_liniarizat = eig(Ac + Bc*K_feedback);
rang_controlabilitate = rank(matrice_controlabilitate(Ac,Bc));

Z = @(x) [1;
          x(1);
          x(2);
          x(1)^2;
          x(1)*x(2);
          x(2)^2;
          x(1)^3;
          x(2)^3];

Y = @(x) [1;
          x(1);
          x(2);
          x(1)^2;
          x(1)*x(2);
          x(2)^2;
          x(1)^3;
          x(2)^3;
          x(1)^4];

W = @(x) [1;
          x(1);
          x(2);
          x(1)^2;
          x(2)^2];

dimZ = numel(Z(x_echilibru));
dimY = numel(Y(x_echilibru));
dimW = numel(W(x_echilibru));

numeZ = {'1','x1','x2','x1^2','x1*x2','x2^2','x1^3','x2^3'};
numeY = {'1','x1','x2','x1^2','x1*x2','x2^2','x1^3','x2^3','x1^4'};
numeW = {'1','x1','x2','x1^2','x2^2'};

nr_coef = n*dimZ + m*dimY + m*dimW;

T_ideal = zeros(n,dimZ);
T_ideal(1,1) = -x_echilibru(1);
T_ideal(1,2) =  1;
T_ideal(2,1) = -x_echilibru(2);
T_ideal(2,3) =  1;

N_ideal = zeros(m,dimY);
N_ideal(1,2) = -k1;
N_ideal(1,7) = -k2;

M_ideal = zeros(m,dimW);
M_ideal(1,1) = b0;

scala_ref = 1/gamma_exact(x_echilibru);
T_ref = scala_ref*T_ideal;
N_ref = scala_ref*N_ideal;
M_ref = scala_ref*M_ideal;
coef_ref = [T_ref(:); N_ref(:); M_ref(:)];

tau_ref = @(x) T_ref*Z(x);
delta_ref = @(x) N_ref*Y(x);
gamma_ref = @(x) M_ref*W(x);

dt = 0.05;
numar_traiectorii = 32;
numar_pasi = 175;
amplitudine_intrare = 1.20;

numar_traiectorii_validare = 8;
numar_pasi_validare = 150;

snr_x_db = 25;

fereastra_filtrare = 13;
ordin_polinom = 3;
margine = floor(fereastra_filtrare/2);
pas_integral = 5;

lambda_grid = logspace(-10,0,65);

rho_jacobian_echilibru = 2e-1;

prag_minJtau = 5e-2;
prag_minGamma = 5e-2;
prag_reziduu_validare = 1e-1;

regiune_test = [x_echilibru(1)-0.80, x_echilibru(1)+0.80;
                x_echilibru(2)-0.80, x_echilibru(2)+0.80];

x0_cl = x_echilibru + [0.45; -0.10];
timp_cl = linspace(0,10,1001);
u_saturatie = 20;

optiuni_ode = odeset('RelTol',1e-6,'AbsTol',1e-8,'MaxStep',0.02,'Events',@eveniment_oprire);

print_titlu('verificari initiale');

fprintf('S1: sistem SISO masa-arc cu rigiditate cubica.\n');
fprintf('Punct de echilibru ales: xe = [%.3f; %.3f]\n',x_echilibru(1),x_echilibru(2));
fprintf('Comanda de echilibru folosita la generarea datelor: u_echilibru = %.4e\n',u_echilibru);
fprintf('Reziduu echilibru model exact ||f(xe)+g(xe)ue|| = %.4e\n',echilibru_reziduu);
fprintf('Grad relativ presupus: r=%d, n=%d -> %s\n',r1,n,text_ok(grad_relativ_complet));
fprintf('Rang controlabilitate(Ac,Bc)=%d din %d -> %s\n',rang_controlabilitate,n,text_ok(rang_controlabilitate==n));
fprintf('Poli forma liniarizata: ');
disp(poli_liniarizat.');
fprintf('Normalizari: tau(xe)=0, gamma(xe)=1, J_tau(xe) aproximativ I.\n');
fprintf('Reziduurile sunt normalizate ca eroare medie pe ecuatie din F(D).\n');
fprintf('Numar coeficienti necunoscuti = %d.\n',nr_coef);

print_titlu('colectare date de baza');

rng(42,'twister');
[x_clean,u_data,xdot_clean] = colecteaza_date_simulare(dt,numar_traiectorii,numar_pasi, ...
    f,g,n,m,amplitudine_intrare,x_echilibru,u_echilibru);

rng(4242,'twister');
[x_val,u_val,xdot_val] = colecteaza_date_simulare(dt,numar_traiectorii_validare,numar_pasi_validare, ...
    f,g,n,m,amplitudine_intrare,x_echilibru,u_echilibru);

putere_medie_x = mean(sum((x_clean - x_echilibru).^2,1));
sigma_x = sqrt(putere_medie_x/(n*10^(snr_x_db/10)));

rng(123,'twister');
x_noisy = x_clean + sigma_x*randn(size(x_clean));

rng(321,'twister');
x_val_noisy = x_val + sigma_x*randn(size(x_val));

fprintf('Puncte identificare: %d\n',size(x_clean,2));
fprintf('Puncte validare separata: %d\n',size(x_val,2));
fprintf('SNR x pentru cazul zgomotos = %.1f dB, sigma zgomot x = %.4e\n',snr_x_db,sigma_x);
fprintf('u este comanda aplicata si salvata, deci este folosita fara zgomot in ambele cazuri.\n');

print_titlu('diagnostic pe date curate exacte');

FD_clean = construieste_FD_derivativ(x_clean,u_data,xdot_clean,Z,Y,W,Ac,Bc,n,m,dimZ,dimY,dimW);
FD_val_clean = construieste_FD_derivativ(x_val,u_val,xdot_val,Z,Y,W,Ac,Bc,n,m,dimZ,dimY,dimW);

rez_ref_clean = calc_reziduu_relativ(FD_clean,coef_ref);
rez_ref_val = calc_reziduu_relativ(FD_val_clean,coef_ref);

fprintf('FD curat exact: [%d x %d]\n',size(FD_clean,1),size(FD_clean,2));
fprintf('Reziduu solutie exacta scalata pe FD curat ID  = %.4e\n',rez_ref_clean);
fprintf('Reziduu solutie exacta scalata pe FD curat VAL = %.4e\n',rez_ref_val);
fprintf('Acest test foloseste xdot exact din simulare si verifica consistenta ecuatiilor.\n');

[Ceq,beq,Cj,j_target] = construieste_conditii_normalizare(x_echilibru,Z,W,n,m,dimZ,dimY,dimW);

print_titlu('identificare comparativa');

cazuri = repmat(init_caz_identificare(),2,1);

cazuri(1) = ruleaza_caz_depersis_curat('Fara zgomot - de baza',Inf,0, ...
    x_clean,x_val,xdot_clean,xdot_val,u_data,u_val,FD_clean,FD_val_clean, ...
    dt,numar_traiectorii,numar_pasi,numar_traiectorii_validare,numar_pasi_validare, ...
    fereastra_filtrare,ordin_polinom,margine,pas_integral, ...
    Z,Y,W,Ac,Bc,n,m,dimZ,dimY,dimW,coef_ref,lambda_grid,Ceq,beq,Cj,j_target, ...
    rho_jacobian_echilibru,tau_ref,delta_ref,gamma_ref,K_feedback,regiune_test, ...
    x_echilibru,prag_minJtau,prag_minGamma,prag_reziduu_validare, ...
    f,g,x0_cl,timp_cl,optiuni_ode,u_saturatie);

cazuri(2) = ruleaza_caz_zgomot('Cu zgomot',snr_x_db,sigma_x, ...
    x_noisy,x_val_noisy,x_clean,x_val,xdot_clean,xdot_val,u_data,u_val, ...
    dt,numar_traiectorii,numar_pasi,numar_traiectorii_validare,numar_pasi_validare, ...
    fereastra_filtrare,ordin_polinom,margine,pas_integral, ...
    Z,Y,W,Ac,Bc,n,m,dimZ,dimY,dimW,coef_ref,lambda_grid,Ceq,beq,Cj,j_target, ...
    rho_jacobian_echilibru,tau_ref,delta_ref,gamma_ref,K_feedback,regiune_test, ...
    x_echilibru,prag_minJtau,prag_minGamma,prag_reziduu_validare, ...
    f,g,x0_cl,timp_cl,optiuni_ode,u_saturatie);

print_titlu('comparatie finala');

for ic = 1:numel(cazuri)
    c = cazuri(ic);
    fprintf('\nCaz: %s\n',c.nume);
    fprintf('  sigma x = %.4e\n',c.sigma_x);
    fprintf('  RMSE x masurat = %.4e\n',c.rmse.x_noisy);
    fprintf('  RMSE xdot filtrat = %.4e\n',c.rmse.xdot_filt);
    fprintf('  reziduu ref integral = %.4e\n',c.rez_ref_int);
    if isnan(c.idx_final)
        fprintf('  metoda finala = nicio metoda acceptata\n');
    else
        rf = c.rezultate(c.idx_final);
        sf = c.sim.metode(c.idx_final);
        fprintf('  metoda finala = %s\n',c.metoda_finala);
        fprintf('  rez VAL = %.4e, minJ = %.4e, min gamma = %.4e\n',rf.rez_val,rf.minJtau,rf.minGamma);
        fprintf('  closed-loop ||x(T)-xe|| = %.4e, max|u| = %.4e\n',sf.err_finala,sf.u_max);
    end
end

raport = struct();
raport.n = n;
raport.m = m;
raport.k1 = k1;
raport.k2 = k2;
raport.b0 = b0;
raport.gamma0_exact = gamma_exact(x_echilibru);
raport.x_echilibru = x_echilibru;
raport.u_echilibru = u_echilibru;
raport.echilibru_reziduu = echilibru_reziduu;
raport.Ac = Ac;
raport.Bc = Bc;
raport.K_feedback = K_feedback;
raport.poli_liniarizat = poli_liniarizat;
raport.rang_controlabilitate = rang_controlabilitate;
raport.grad_relativ_complet = grad_relativ_complet;
raport.r1 = r1;
raport.Z = Z;
raport.Y = Y;
raport.W = W;
raport.dimZ = dimZ;
raport.dimY = dimY;
raport.dimW = dimW;
raport.numeZ = numeZ;
raport.numeY = numeY;
raport.numeW = numeW;
raport.nr_coef = nr_coef;
raport.T_ref = T_ref;
raport.N_ref = N_ref;
raport.M_ref = M_ref;
raport.coef_ref = coef_ref;
raport.tau_ref = tau_ref;
raport.delta_ref = delta_ref;
raport.gamma_ref = gamma_ref;
raport.dt = dt;
raport.numar_traiectorii = numar_traiectorii;
raport.numar_pasi = numar_pasi;
raport.numar_traiectorii_validare = numar_traiectorii_validare;
raport.numar_pasi_validare = numar_pasi_validare;
raport.snr_x_db = snr_x_db;
raport.sigma_x = sigma_x;
raport.x_clean = x_clean;
raport.x_noisy = x_noisy;
raport.xdot_clean = xdot_clean;
raport.u_data = u_data;
raport.FD_clean = FD_clean;
raport.FD_val_clean = FD_val_clean;
raport.rez_ref_clean = rez_ref_clean;
raport.rez_ref_val = rez_ref_val;
raport.Ceq = Ceq;
raport.beq = beq;
raport.Cj = Cj;
raport.j_target = j_target;
raport.rho_jacobian_echilibru = rho_jacobian_echilibru;
raport.prag_minJtau = prag_minJtau;
raport.prag_minGamma = prag_minGamma;
raport.prag_reziduu_validare = prag_reziduu_validare;
raport.timp_cl = timp_cl;
raport.cazuri = cazuri;

try
    deschide_dashboard_S1(raport);
catch ME
    warning('Dashboard-ul nu a putut fi deschis: %s',ME.message);
    afiseaza_eroare_dashboard(ME);
end

function [x_data,u_data,xdot_data] = colecteaza_date_simulare(dt,ntraj,npasi,f,g,n,m,amp,x_eq,u_eq)
    x_data = zeros(n,ntraj*npasi);
    u_data = zeros(m,ntraj*npasi);
    xdot_data = zeros(n,ntraj*npasi);
    idx0 = 0;

    for it = 1:ntraj
        x = x_eq + [0.90; 0.70].*(rand(n,1)-0.5);

        for k = 1:npasi
            idx = idx0 + k;

            if mod(it,3) == 0
                u = u_eq + amp*(2*rand(m,1)-1);
            elseif mod(it,3) == 1
                u = u_eq + 0.75*amp*sin(0.15*k + 0.10*it) + 0.20*amp*(2*rand(m,1)-1);
            else
                u = u_eq + 0.55*amp*cos(0.11*k - 0.05*it) + 0.25*amp*(2*rand(m,1)-1);
            end

            u = max(min(u,u_eq+amp),u_eq-amp);

            xdot = f(x) + g(x)*u;
            x_data(:,idx) = x;
            u_data(:,idx) = u;
            xdot_data(:,idx) = xdot;

            [~,xs] = ode45(@(~,xx) f(xx) + g(xx)*u,[0 dt],x);
            x = xs(end,:)';

            if norm(x-x_eq) > 3.5
                x = x_eq + [0.90; 0.70].*(rand(n,1)-0.5);
            end
        end

        idx0 = idx0 + npasi;
    end
end

function [x_est,xdot_est] = estimeaza_brut(x_noisy,dt,ntraj,npasi)
    x_est = zeros(size(x_noisy));
    xdot_est = zeros(size(x_noisy));

    for it = 1:ntraj
        idx = (it-1)*npasi + (1:npasi);
        for s = 1:size(x_noisy,1)
            semnal = x_noisy(s,idx);
            x_est(s,idx) = semnal;
            xdot_est(s,idx) = gradient(semnal,dt);
        end
    end
end

function [x_est,xdot_est,rez_punct] = estimeaza_polinom_local(x_noisy,dt,ntraj,npasi,fereastra,ordin)
    if mod(fereastra,2) == 0
        fereastra = fereastra + 1;
    end

    jumatate = floor(fereastra/2);
    ordin = max(2,min(ordin,fereastra-2));

    x_est = zeros(size(x_noisy));
    xdot_est = zeros(size(x_noisy));
    rez_punct = zeros(1,size(x_noisy,2));
    n = size(x_noisy,1);

    for it = 1:ntraj
        idx_tr = (it-1)*npasi + (1:npasi);
        for s = 1:n
            semnal = x_noisy(s,idx_tr);
            for k = 1:npasi
                k1 = max(1,k-jumatate);
                k2 = min(npasi,k+jumatate);
                idx_loc = k1:k2;
                tloc = ((idx_loc-k)')*dt;
                yloc = semnal(idx_loc)';
                ordin_loc = min(ordin,numel(idx_loc)-1);

                p = polyfit(tloc,yloc,ordin_loc);
                x_est(s,idx_tr(k)) = polyval(p,0);
                dp = polyder(p);
                xdot_est(s,idx_tr(k)) = polyval(dp,0);

                yfit = polyval(p,tloc);
                rez = norm(yfit-yloc)/max(sqrt(numel(idx_loc)),1);
                rez_punct(idx_tr(k)) = rez_punct(idx_tr(k)) + rez/n;
            end
        end
    end
end

function idx = construieste_index_interior(ntraj,npasi,margine)
    idx = [];
    for it = 1:ntraj
        ii = (it-1)*npasi + (1:npasi);
        idx = [idx ii(margine+1:npasi-margine)];
    end
end

function FD = construieste_FD_derivativ(x,u,xdot,Z,Y,W,Ac,Bc,n,m,dimZ,dimY,dimW)
    L = size(x,2);
    nr_coef = n*dimZ + m*dimY + m*dimW;
    FD = zeros(n*L,nr_coef);

    for i = 1:L
        xi = x(:,i);
        ui = u(:,i);
        xdi = xdot(:,i);

        Zi = Z(xi);
        Yi = Y(xi);
        Wi = W(xi);

        dZ = jacobian_numeric(Z,xi,dimZ,n);

        blocT = kron(Zi',Ac) - kron((dZ*xdi)',eye(n));
        blocN = kron(Yi',Bc);
        blocM = kron((Wi*ui)',Bc);

        FD((i-1)*n+(1:n),:) = [blocT blocN blocM];
    end
end

function FD = construieste_FD_integral(x,u,Z,Y,W,Ac,Bc,n,m,dimZ,dimY,dimW,dt,ntraj,npasi,pas_int,margine)
    nr_max = ntraj*max(npasi - 2*margine - pas_int,1);
    nr_coef = n*dimZ + m*dimY + m*dimW;
    FD = zeros(n*nr_max,nr_coef);
    lin = 0;

    for it = 1:ntraj
        idx_tr = (it-1)*npasi + (1:npasi);
        k_start = margine + 1;
        k_stop = npasi - margine - pas_int;

        for k = k_start:k_stop
            idx = idx_tr(k:k+pas_int);
            tloc = (0:pas_int)*dt;

            Zval = zeros(dimZ,numel(idx));
            Yval = zeros(dimY,numel(idx));
            Wuval = zeros(dimW,numel(idx));

            for q = 1:numel(idx)
                xq = x(:,idx(q));
                uq = u(:,idx(q));
                Zval(:,q) = Z(xq);
                Yval(:,q) = Y(xq);
                Wuval(:,q) = W(xq)*uq;
            end

            deltaZ = Zval(:,end)-Zval(:,1);
            intZ = trapz(tloc,Zval,2);
            intY = trapz(tloc,Yval,2);
            intWu = trapz(tloc,Wuval,2);

            blocT = kron(deltaZ',eye(n)) - kron(intZ',Ac);
            blocN = -kron(intY',Bc);
            blocM = -kron(intWu',Bc);

            lin = lin + 1;
            FD((lin-1)*n+(1:n),:) = [blocT blocN blocM];
        end
    end

    FD = FD(1:lin*n,:);
end

function J = jacobian_numeric(fun,x,dim,n)
    h = 1e-6;
    J = zeros(dim,n);

    for k = 1:n
        e = zeros(n,1);
        e(k) = h;
        J(:,k) = (fun(x+e)-fun(x-e))/(2*h);
    end
end

function [Ceq,beq,Cj,j_target] = construieste_conditii_normalizare(x_eq,Z,W,n,m,dimZ,dimY,dimW)
    nr_coef = n*dimZ + m*dimY + m*dimW;
    offsetM = n*dimZ + m*dimY;

    Z_eq = Z(x_eq);
    Ctau = zeros(n,nr_coef);
    Ctau(:,1:n*dimZ) = kron(Z_eq',eye(n));

    W_eq = W(x_eq);
    Cgamma = zeros(1,nr_coef);
    Cgamma(offsetM + (1:dimW)) = W_eq(:)';

    Ceq = [Cgamma; Ctau];
    beq = [1; zeros(n,1)];

    dZ_eq = jacobian_numeric(Z,x_eq,dimZ,n);
    Cj = zeros(n*n,nr_coef);
    Cj(:,1:n*dimZ) = kron(dZ_eq',eye(n));
    I = eye(n);
    j_target = I(:);
end

function rezultat = identifica_depersis_svd(nume,FD,FD_val,Ceq,beq,Cj,j_target,rhoJ,n,m,dimZ,dimY,dimW,Z,Y,W,tau_ref,delta_ref,gamma_ref,Kfb,regiune,x_eq,pragJ,pragG,pragRez)

    [FDs,scale] = scaleaza_coloane_cu_scara(FD);
    Ceq_z = Ceq ./ scale(:)';
    Cj_z = Cj ./ scale(:)';

    rezultat = init_rezultat();
    rezultat.nume = nume;
    rezultat.lambda = 0;

    if isempty(FDs)
        rezultat.status = 'RESPINS';
        return;
    end

    [~,S,V] = svd(FDs,'econ');
    valori_singulare = diag(S);

    if isempty(valori_singulare)
        rezultat.status = 'RESPINS';
        return;
    end
    smax = max(valori_singulare);
    if smax < eps
        smax = 1;
    end
    valori_relative = valori_singulare/smax;
    prag_spatiu_nul = 1e-8;
    idx_null = find(valori_relative <= prag_spatiu_nul);

    if isempty(idx_null)
        rezultat.status = 'RESPINS_SPATIU_NUL';
        rezultat.rez_train = valori_relative(end);
        rezultat.rez_val = NaN;
        rezultat.observatie = sprintf('Nu s-a gasit spatiu nul numeric: sigma_min/sigma_max = %.4e.',valori_relative(end));
        return;
    end

    baza_null = V(:,idx_null);
    [a,ok_norm,eroare_norm] = rezolva_in_spatiul_nul(baza_null,Ceq_z,beq,Cj_z,j_target,rhoJ);

    if ~ok_norm
        rezultat.status = 'RESPINS_NORMALIZARE';
        rezultat.observatie = sprintf('Metoda De Persis directa pe date zgomotoase nu produce o solutie admisibila. Eroare normalizare = %.4e.',eroare_norm);
        return;
    end

    z = baza_null*a;
    v = z ./ scale(:);

    if any(~isfinite(v))
        rezultat.status = 'RESPINS';
        rezultat.observatie = 'Solutia SVD contine valori numerice nefinite.';
        return;
    end

    [T,N,M] = extrage_TNM(v,n,m,dimZ,dimY,dimW);
    model = creeaza_model(T,N,M,Z,Y,W);

    rez_train = calc_reziduu_relativ(FD,v);
    rez_val = calc_reziduu_relativ(FD_val,v);
    [minJ,minG,u_rms_reg] = evalueaza_admisibilitate_si_comanda(model,Kfb,regiune,n);
    metric = evalueaza_model_pe_grila(model,tau_ref,delta_ref,gamma_ref,Kfb,regiune);

    tau_ech_norm = norm(model.tau(x_eq));
    gamma_ech = model.gamma(x_eq);
    errJ = norm(Cj*v - j_target);

    rezultat.T = T;
    rezultat.N = N;
    rezultat.M = M;
    rezultat.coef = v;
    rezultat.model = model;
    rezultat.rez_train = rez_train;
    rezultat.rez_val = rez_val;
    rezultat.minJtau = minJ;
    rezultat.minGamma = minG;
    rezultat.u_rms_reg = u_rms_reg;
    rezultat.errJ_ech = errJ;
    rezultat.tau_ech_norm = tau_ech_norm;
    rezultat.gamma_ech = gamma_ech;
    rezultat.coef_norm = norm(v);
    rezultat.err_tau = metric.err_tau;
    rezultat.err_delta = metric.err_delta;
    rezultat.err_gamma = metric.err_gamma;
    rezultat.err_u = metric.err_u;
    rezultat.score = rez_val + 0.10*rez_train;

    if rezultat.minJtau >= pragJ && rezultat.minGamma >= pragG && rezultat.rez_val <= pragRez
        rezultat.status = 'OK';
    elseif rezultat.rez_val > pragRez
        rezultat.status = 'RESPINS_REZIDUU';
    else
        rezultat.status = 'RESPINS_ADMIS';
    end
end

function [a,ok_norm,eroare_norm] = rezolva_in_spatiul_nul(baza_null,Ceq_z,beq,Cj_z,j_target,~)

    k = size(baza_null,2);
    Aeq = Ceq_z*baza_null;
    Aj = Cj_z*baza_null;

    a = NaN(k,1);
    ok_norm = false;
    eroare_norm = inf;

    if k == 0 || isempty(Aeq)
        return;
    end

    a0 = pinv_svd(Aeq)*beq;
    eroare_norm = norm(Aeq*a0 - beq);

    if eroare_norm > 1e-7
        return;
    end

    Nliber = null(Aeq,'r');

    if isempty(Nliber)
        a = a0;
    else
        Ared = Aj*Nliber;
        bred = j_target - Aj*a0;
        y = pinv_svd(Ared)*bred;
        a = a0 + Nliber*y;
    end

    eroare_norm = norm(Aeq*a - beq);
    ok_norm = isfinite(eroare_norm) && eroare_norm <= 1e-7 && all(isfinite(a));
end

function X = pinv_svd(A)
    if isempty(A)
        X = A';
        return;
    end

    [U,S,V] = svd(A,'econ');
    s = diag(S);

    if isempty(s)
        X = zeros(size(A,2),size(A,1));
        return;
    end

    tol = max(size(A))*eps(max(s));
    sinv = zeros(size(s));
    idx = s > tol;
    sinv(idx) = 1./s(idx);
    X = V*diag(sinv)*U';
end

function rezultat = identifica_tikhonov(nume,FD,FD_val,lambda_grid,Ceq,beq,Cj,j_target,rhoJ,n,m,dimZ,dimY,dimW,Z,Y,W,tau_ref,delta_ref,gamma_ref,Kfb,regiune,x_eq,pragJ,pragG,pragRez)
    [FDs,scale] = scaleaza_coloane_cu_scara(FD);
    Ceq_z = Ceq ./ scale(:)';
    Cj_z = Cj ./ scale(:)';

    p = size(FDs,2);
    rezultat = init_rezultat();
    rezultat.nume = nume;

    best_score = inf;

    for il = 1:numel(lambda_grid)
        lambda = lambda_grid(il);

        H = (FDs'*FDs)/max(size(FDs,1),1) + lambda*eye(p) + rhoJ*(Cj_z'*Cj_z);
        rhs0 = rhoJ*(Cj_z'*j_target);

        KKT = [H Ceq_z'; Ceq_z zeros(size(Ceq_z,1))];
        rhs = [rhs0; beq];

        sol = pinv_svd(KKT)*rhs;

        z = sol(1:p);
        v = z ./ scale(:);

        if any(~isfinite(v))
            continue;
        end

        [T,N,M] = extrage_TNM(v,n,m,dimZ,dimY,dimW);
        model = creeaza_model(T,N,M,Z,Y,W);

        rez_train = calc_reziduu_relativ(FD,v);
        rez_val = calc_reziduu_relativ(FD_val,v);
        [minJ,minG,u_rms_reg] = evalueaza_admisibilitate_si_comanda(model,Kfb,regiune,n);
        metric = evalueaza_model_pe_grila(model,tau_ref,delta_ref,gamma_ref,Kfb,regiune);

        tau_ech_norm = norm(model.tau(x_eq));
        gamma_ech = model.gamma(x_eq);
        errJ = norm(Cj*v - j_target);

        penal = 0;
        penal = penal + 20*max(0,pragJ-minJ)^2;
        penal = penal + 20*max(0,pragG-minG)^2;
        penal = penal + 0.02*errJ^2;
        penal = penal + 0.01*u_rms_reg;

        score = rez_val + 0.10*rez_train + penal;

        if score < best_score
            best_score = score;
            rezultat.T = T;
            rezultat.N = N;
            rezultat.M = M;
            rezultat.coef = v;
            rezultat.model = model;
            rezultat.lambda = lambda;
            rezultat.rez_train = rez_train;
            rezultat.rez_val = rez_val;
            rezultat.minJtau = minJ;
            rezultat.minGamma = minG;
            rezultat.u_rms_reg = u_rms_reg;
            rezultat.errJ_ech = errJ;
            rezultat.tau_ech_norm = tau_ech_norm;
            rezultat.gamma_ech = gamma_ech;
            rezultat.coef_norm = norm(v);
            rezultat.err_tau = metric.err_tau;
            rezultat.err_delta = metric.err_delta;
            rezultat.err_gamma = metric.err_gamma;
            rezultat.err_u = metric.err_u;
            rezultat.score = score;
        end
    end

    if isempty(rezultat.coef)
        rezultat.status = 'RESPINS_ADMIS';
        return;
    end

    if rezultat.minJtau >= pragJ && rezultat.minGamma >= pragG && rezultat.rez_val <= pragRez
        rezultat.status = 'OK';
    elseif rezultat.rez_val > pragRez
        rezultat.status = 'RESPINS_REZIDUU';
    else
        rezultat.status = 'RESPINS_ADMIS';
    end
end

function r = init_rezultat()
    r = struct();
    r.nume = '';
    r.T = [];
    r.N = [];
    r.M = [];
    r.coef = [];
    r.model = [];
    r.lambda = NaN;
    r.rez_train = inf;
    r.rez_val = inf;
    r.minJtau = 0;
    r.minGamma = 0;
    r.u_rms_reg = inf;
    r.errJ_ech = inf;
    r.tau_ech_norm = inf;
    r.gamma_ech = NaN;
    r.coef_norm = inf;
    r.err_tau = NaN;
    r.err_delta = NaN;
    r.err_gamma = NaN;
    r.err_u = NaN;
    r.score = inf;
    r.status = 'RESPINS';
    r.observatie = '';
    r.sim = [];
end

function [T,N,M] = extrage_TNM(v,n,m,dimZ,dimY,dimW)
    iT = n*dimZ;
    iN = iT + m*dimY;
    T = reshape(v(1:iT),[n,dimZ]);
    N = reshape(v(iT+1:iN),[m,dimY]);
    M = reshape(v(iN+1:iN+m*dimW),[m,dimW]);
end

function model = creeaza_model(T,N,M,Z,Y,W)
    model.T = T;
    model.N = N;
    model.M = M;
    model.tau = @(x) T*Z(x);
    model.delta = @(x) N*Y(x);
    model.gamma = @(x) M*W(x);
end

function metric = evalueaza_model_pe_grila(model,tau_ref,delta_ref,gamma_ref,Kfb,regiune)
    x1v = linspace(regiune(1,1),regiune(1,2),25);
    x2v = linspace(regiune(2,1),regiune(2,2),25);

    e_tau = [];
    e_delta = [];
    e_gamma = [];
    e_u = [];

    for x1 = x1v
        for x2 = x2v
            x = [x1;x2];

            tau_h = model.tau(x);
            tau_e = tau_ref(x);
            delta_h = model.delta(x);
            delta_e = delta_ref(x);
            gamma_h = model.gamma(x);
            gamma_e = gamma_ref(x);

            e_tau(end+1,1) = norm(tau_h-tau_e)/max(norm(tau_e),1e-3);
            e_delta(end+1,1) = norm(delta_h-delta_e)/max(norm(delta_e),0.1);
            e_gamma(end+1,1) = abs(gamma_h-gamma_e)/max(abs(gamma_e),0.1);

            u_e = control_fl(x,tau_ref,delta_ref,gamma_ref,Kfb,20);
            u_h = control_fl(x,model.tau,model.delta,model.gamma,Kfb,20);
            e_u(end+1,1) = abs(u_h-u_e);
        end
    end

    metric.err_tau = mean(e_tau)*100;
    metric.err_delta = mean(e_delta)*100;
    metric.err_gamma = mean(e_gamma)*100;
    metric.err_u = mean(e_u);
end

function [minJ,minG,u_rms] = evalueaza_admisibilitate_si_comanda(model,Kfb,regiune,n)
    x1v = linspace(regiune(1,1),regiune(1,2),21);
    x2v = linspace(regiune(2,1),regiune(2,2),21);

    minJ = inf;
    minG = inf;
    uvals = [];
    h = 1e-6;

    for x1 = x1v
        for x2 = x2v
            x = [x1;x2];

            J = zeros(n,n);
            for k = 1:n
                e = zeros(n,1);
                e(k) = h;
                J(:,k) = (model.tau(x+e)-model.tau(x-e))/(2*h);
            end

            sJ = svd(J);
            minJ = min(minJ,sJ(end));
            minG = min(minG,abs(model.gamma(x)));

            uvals(end+1,1) = control_fl(x,model.tau,model.delta,model.gamma,Kfb,20);
        end
    end

    u_rms = sqrt(mean(uvals.^2));
end

function u = control_fl(x,tau,delta,gamma,Kfb,u_sat)
    G = gamma(x);
    rhs = Kfb*tau(x) - delta(x);

    if abs(G) < 1e-8
        G = sign(G + (G==0))*1e-8;
    end

    u = rhs/G;
    u = max(min(u,u_sat),-u_sat);
end

function s = simuleaza_un_caz(f,g,tau,delta,gamma,Kfb,x0,x_eq,tspan,opts,u_sat,prag_final,prag_u)
    s = init_simulare(numel(x0),1);
    lege = @(x) control_fl(x,tau,delta,gamma,Kfb,u_sat);

    try
        [tout,xout] = ode45(@(~,x) f(x)+g(x)*lege(x),tspan,x0,opts);
        X = xout';
        U = zeros(1,numel(tout));

        for k = 1:numel(tout)
            U(k) = lege(X(:,k));
        end

        err = vecnorm(X - x_eq);
        norma_u = abs(U);

        s.simulat = true;
        s.T = tout';
        s.X = X;
        s.U = U;
        s.err_traj = err;
        s.err_finala = err(end);
        s.u_max = max(norma_u);
        s.u_rms = sqrt(mean(norma_u.^2));
        s.ok = s.err_finala <= prag_final && s.u_max <= prag_u;
    catch
        s.simulat = false;
        s.ok = false;
    end
end

function s = init_simulare(n,m)
    s = struct();
    s.simulat = false;
    s.ok = false;
    s.T = [];
    s.X = NaN(n,0);
    s.U = NaN(m,0);
    s.err_traj = [];
    s.err_finala = NaN;
    s.u_max = NaN;
    s.u_rms = NaN;
end

function idx_final = alege_metoda_finala(rezultate,sim)

    idx_final = selecteaza_metoda(rezultate,sim,false);

    if isnan(idx_final)
        idx_final = selecteaza_metoda(rezultate,sim,true);
    end
end

function idx = selecteaza_metoda(rezultate,sim,include_test_direct)
    idx = NaN;
    best = inf;

    for i = 1:numel(rezultate)
        if ~include_test_direct && contains(rezultate(i).nume,'direct pe zgomot')
            continue;
        end

        if ~strcmp(rezultate(i).status,'OK')
            continue;
        end

        if isempty(sim.metode(i).T) || ~sim.metode(i).simulat || ~sim.metode(i).ok
            continue;
        end

        scor = rezultate(i).rez_val + 0.5*sim.metode(i).err_finala + 0.02*sim.metode(i).u_rms;

        if scor < best
            best = scor;
            idx = i;
        end
    end
end

function caz = init_caz_identificare()
    caz = struct();
    caz.nume = "";
    caz.snr_x_db = NaN;
    caz.sigma_x = NaN;
    caz.x_masurat = [];
    caz.x_val_masurat = [];
    caz.x_raw = [];
    caz.xdot_raw = [];
    caz.x_filt = [];
    caz.xdot_filt = [];
    caz.x_val_raw = [];
    caz.xdot_val_raw = [];
    caz.x_val_filt = [];
    caz.xdot_val_filt = [];
    caz.idx_interior = [];
    caz.idx_val_interior = [];
    caz.rmse = struct('x_noisy',NaN,'x_filt',NaN,'xdot_raw',NaN,'xdot_filt',NaN);
    caz.FD_raw = [];
    caz.FD_filt = [];
    caz.FD_int = [];
    caz.FD_val_raw = [];
    caz.FD_val_filt = [];
    caz.FD_val_int = [];
    caz.rez_ref_int = NaN;
    caz.rez_ref_val_int = NaN;
    caz.sv_int = [];
    caz.rezultate = [];
    caz.sim = [];
    caz.idx_final = NaN;
    caz.metoda_finala = 'nicio metoda acceptata';
    caz.prag_eroare_finala = NaN;
    caz.prag_umax = NaN;
end

function caz = ruleaza_caz_depersis_curat(nume_caz,snr_val,sigma_x,x_clean,x_val,xdot_clean,xdot_val,u_data,u_val,FD_clean,FD_val_clean, ...
    dt,ntraj,npasi,ntraj_val,npasi_val,fereastra,ordin,margine,pas_integral, ...
    Z,Y,W,Ac,Bc,n,m,dimZ,dimY,dimW,coef_ref,lambda_grid,Ceq,beq,Cj,j_target,rhoJ, ...
    tau_ref,delta_ref,gamma_ref,Kfb,regiune,x_eq,pragJ,pragG,pragRez,f,g,x0_cl,tspan,opts,u_sat)

    print_titlu(['caz ' char(nume_caz)]);

    caz = init_caz_identificare();
    caz.nume = string(nume_caz);
    caz.snr_x_db = snr_val;
    caz.sigma_x = sigma_x;
    caz.x_masurat = x_clean;
    caz.x_val_masurat = x_val;

    [caz.x_raw,caz.xdot_raw] = estimeaza_brut(x_clean,dt,ntraj,npasi);
    caz.x_filt = x_clean;
    caz.xdot_filt = xdot_clean;
    caz.idx_interior = 1:size(x_clean,2);

    [caz.x_val_raw,caz.xdot_val_raw] = estimeaza_brut(x_val,dt,ntraj_val,npasi_val);
    caz.x_val_filt = x_val;
    caz.xdot_val_filt = xdot_val;
    caz.idx_val_interior = 1:size(x_val,2);

    caz.rmse.x_noisy = 0;
    caz.rmse.x_filt = 0;
    caz.rmse.xdot_raw = calc_rmse(caz.xdot_raw,xdot_clean);
    caz.rmse.xdot_filt = 0;

    caz.FD_raw = FD_clean;
    caz.FD_filt = FD_clean;
    caz.FD_val_raw = FD_val_clean;
    caz.FD_val_filt = FD_val_clean;

    caz.FD_int = construieste_FD_integral(x_clean,u_data,Z,Y,W,Ac,Bc,n,m,dimZ,dimY,dimW,dt,ntraj,npasi,pas_integral,margine);
    caz.FD_val_int = construieste_FD_integral(x_val,u_val,Z,Y,W,Ac,Bc,n,m,dimZ,dimY,dimW,dt,ntraj_val,npasi_val,pas_integral,margine);

    caz.sv_int = svd(scaleaza_coloane(caz.FD_int),'econ');
    caz.rez_ref_int = calc_reziduu_relativ(caz.FD_int,coef_ref);
    caz.rez_ref_val_int = calc_reziduu_relativ(caz.FD_val_int,coef_ref);

    fprintf('SNR x = %s, sigma x = %.4e\n',snr_text(snr_val),sigma_x);
    fprintf('Metoda fara zgomot: din De Persis prin spatiul nul SVD al lui F(D).\n');
    fprintf('RMSE x masurat = %.4e\n',caz.rmse.x_noisy);
    fprintf('RMSE xdot brut, diagnostic = %.4e\n',caz.rmse.xdot_raw);
    fprintf('FD De Persis: [%d x %d]\n',size(caz.FD_filt,1),size(caz.FD_filt,2));
    fprintf('FD integral curat:  [%d x %d]\n',size(caz.FD_int,1),size(caz.FD_int,2));
    fprintf('Reziduu solutie de referinta pe FD De Persis = %.4e\n',calc_reziduu_relativ(caz.FD_filt,coef_ref));
    fprintf('Reziduu solutie exacta scalata pe FD integral  = %.4e\n',caz.rez_ref_int);

    caz.rezultate = repmat(init_rezultat(),1,1);

    caz.rezultate(1) = identifica_depersis_svd('spatiu nul SVD', ...
        caz.FD_filt,caz.FD_val_filt,Ceq,beq,Cj,j_target,rhoJ, ...
        n,m,dimZ,dimY,dimW,Z,Y,W,tau_ref,delta_ref,gamma_ref,Kfb, ...
        regiune,x_eq,pragJ,pragG,pragRez);

    for i = 1:numel(caz.rezultate)
        fprintf('\nMetoda: %s\n',caz.rezultate(i).nume);
        if contains(caz.rezultate(i).nume,'De Persis')
            fprintf('  rezolvare                = spatiu nul SVD, fara Tikhonov\n');
        else
            fprintf('  lambda ales              = %.4e\n',caz.rezultate(i).lambda);
        end
        if isempty(caz.rezultate(i).coef)
            fprintf('  observatie               = %s\n',caz.rezultate(i).observatie);
            fprintf('  status                   = %s\n',caz.rezultate(i).status);
            continue;
        end
        fprintf('  reziduu identificare     = %.4e\n',caz.rezultate(i).rez_train);
        fprintf('  reziduu validare         = %.4e\n',caz.rezultate(i).rez_val);
        fprintf('  min sigma(J_tau)         = %.4e\n',caz.rezultate(i).minJtau);
        fprintf('  min |gamma|              = %.4e\n',caz.rezultate(i).minGamma);
        fprintf('  ||tau(xe)||             = %.4e\n',caz.rezultate(i).tau_ech_norm);
        fprintf('  |gamma(xe)-1|           = %.4e\n',abs(caz.rezultate(i).gamma_ech-1));
        fprintf('  err J_tau(xe)           = %.4e\n',caz.rezultate(i).errJ_ech);
        fprintf('  status                   = %s\n',caz.rezultate(i).status);
    end

    caz.prag_eroare_finala = 0.03;
    caz.prag_umax = 18;

    caz.sim = struct();
    caz.sim.exact = simuleaza_un_caz(f,g,tau_ref,delta_ref,gamma_ref,Kfb,x0_cl,x_eq,tspan,opts,u_sat,caz.prag_eroare_finala,caz.prag_umax);
    caz.sim.metode = repmat(init_simulare(n,m),numel(caz.rezultate),1);

    fprintf('\nReferinta exacta ||x(T)-xe||=%.4e, max|u|=%.4e, RMS|u|=%.4e\n', ...
        caz.sim.exact.err_finala,caz.sim.exact.u_max,caz.sim.exact.u_rms);

    for i = 1:numel(caz.rezultate)
        if ~strcmp(caz.rezultate(i).status,'OK')
            fprintf('%-24s metoda respinsa numeric: %s\n',caz.rezultate(i).nume,status_afisat(caz.rezultate(i).status));
            continue;
        end

        caz.sim.metode(i) = simuleaza_un_caz(f,g,caz.rezultate(i).model.tau,caz.rezultate(i).model.delta,caz.rezultate(i).model.gamma, ...
            Kfb,x0_cl,x_eq,tspan,opts,u_sat,caz.prag_eroare_finala,caz.prag_umax);

        fprintf('%-24s ||x(T)-xe||=%.4e, max|u|=%.4e, RMS|u|=%.4e, closed-loop=%s\n', ...
            caz.rezultate(i).nume,caz.sim.metode(i).err_finala,caz.sim.metode(i).u_max,caz.sim.metode(i).u_rms,text_ok(caz.sim.metode(i).ok));

        if ~caz.sim.metode(i).ok
            caz.rezultate(i).status = 'RESPINS_CL';
        end
    end

    caz.idx_final = alege_metoda_finala(caz.rezultate,caz.sim);
    if isnan(caz.idx_final)
        caz.metoda_finala = 'nicio metoda acceptata';
        fprintf('\nNu exista metoda finala acceptata pentru cazul %s.\n',nume_caz);
    else
        caz.metoda_finala = caz.rezultate(caz.idx_final).nume;
        fprintf('\nMetoda finala pentru cazul %s: %s\n',nume_caz,caz.metoda_finala);
    end
end

function caz = ruleaza_caz_zgomot(nume_caz,snr_val,sigma_x,x_masurat,x_val_masurat,x_clean,x_val,xdot_clean,xdot_val,u_data,u_val, ...
    dt,ntraj,npasi,ntraj_val,npasi_val,fereastra,ordin,margine,pas_integral, ...
    Z,Y,W,Ac,Bc,n,m,dimZ,dimY,dimW,coef_ref,lambda_grid,Ceq,beq,Cj,j_target,rhoJ, ...
    tau_ref,delta_ref,gamma_ref,Kfb,regiune,x_eq,pragJ,pragG,pragRez,f,g,x0_cl,tspan,opts,u_sat)

    print_titlu(['caz ' char(nume_caz)]);

    caz = init_caz_identificare();
    caz.nume = string(nume_caz);
    caz.snr_x_db = snr_val;
    caz.sigma_x = sigma_x;
    caz.x_masurat = x_masurat;
    caz.x_val_masurat = x_val_masurat;

    [caz.x_raw,caz.xdot_raw] = estimeaza_brut(x_masurat,dt,ntraj,npasi);
    [caz.x_filt,caz.xdot_filt,~] = estimeaza_polinom_local(x_masurat,dt,ntraj,npasi,fereastra,ordin);
    caz.idx_interior = construieste_index_interior(ntraj,npasi,margine);

    [caz.x_val_raw,caz.xdot_val_raw] = estimeaza_brut(x_val_masurat,dt,ntraj_val,npasi_val);
    [caz.x_val_filt,caz.xdot_val_filt,~] = estimeaza_polinom_local(x_val_masurat,dt,ntraj_val,npasi_val,fereastra,ordin);
    caz.idx_val_interior = construieste_index_interior(ntraj_val,npasi_val,margine);

    caz.rmse.x_noisy = calc_rmse(x_masurat,x_clean);
    caz.rmse.x_filt = calc_rmse(caz.x_filt(:,caz.idx_interior),x_clean(:,caz.idx_interior));
    caz.rmse.xdot_raw = calc_rmse(caz.xdot_raw(:,caz.idx_interior),xdot_clean(:,caz.idx_interior));
    caz.rmse.xdot_filt = calc_rmse(caz.xdot_filt(:,caz.idx_interior),xdot_clean(:,caz.idx_interior));

    caz.FD_raw = construieste_FD_derivativ(caz.x_raw,u_data,caz.xdot_raw,Z,Y,W,Ac,Bc,n,m,dimZ,dimY,dimW);
    caz.FD_filt = construieste_FD_derivativ(caz.x_filt(:,caz.idx_interior),u_data(:,caz.idx_interior),caz.xdot_filt(:,caz.idx_interior),Z,Y,W,Ac,Bc,n,m,dimZ,dimY,dimW);
    caz.FD_int = construieste_FD_integral(caz.x_filt,u_data,Z,Y,W,Ac,Bc,n,m,dimZ,dimY,dimW,dt,ntraj,npasi,pas_integral,margine);

    caz.FD_val_raw = construieste_FD_derivativ(caz.x_val_raw,u_val,caz.xdot_val_raw,Z,Y,W,Ac,Bc,n,m,dimZ,dimY,dimW);
    caz.FD_val_filt = construieste_FD_derivativ(caz.x_val_filt(:,caz.idx_val_interior),u_val(:,caz.idx_val_interior),caz.xdot_val_filt(:,caz.idx_val_interior),Z,Y,W,Ac,Bc,n,m,dimZ,dimY,dimW);
    caz.FD_val_int = construieste_FD_integral(caz.x_val_filt,u_val,Z,Y,W,Ac,Bc,n,m,dimZ,dimY,dimW,dt,ntraj_val,npasi_val,pas_integral,margine);

    caz.sv_int = svd(scaleaza_coloane(caz.FD_int),'econ');
    caz.rez_ref_int = calc_reziduu_relativ(caz.FD_int,coef_ref);
    caz.rez_ref_val_int = calc_reziduu_relativ(caz.FD_val_int,coef_ref);

    fprintf('SNR x = %s, sigma x = %.4e\n',snr_text(snr_val),sigma_x);
    fprintf('RMSE x masurat = %.4e\n',caz.rmse.x_noisy);
    fprintf('RMSE x filtrat = %.4e\n',caz.rmse.x_filt);
    fprintf('RMSE xdot brut = %.4e\n',caz.rmse.xdot_raw);
    fprintf('RMSE xdot filtrat = %.4e\n',caz.rmse.xdot_filt);
    fprintf('FD brut, pentru De Persis direct pe zgomot: [%d x %d]\n',size(caz.FD_raw,1),size(caz.FD_raw,2));
    fprintf('FD filtrat, pentru forma derivativa + regularizare Tikhonov:   [%d x %d]\n',size(caz.FD_filt,1),size(caz.FD_filt,2));
    fprintf('FD integral, pentru forma integrala + regularizare Tikhonov:   [%d x %d]\n',size(caz.FD_int,1),size(caz.FD_int,2));
    fprintf('Reziduu solutie exacta scalata pe FD integral = %.4e\n',caz.rez_ref_int);

    caz.rezultate = repmat(init_rezultat(),3,1);

    caz.rezultate(1) = identifica_depersis_svd('Direct pe zgomot', ...
        caz.FD_raw,caz.FD_val_raw,Ceq,beq,Cj,j_target,rhoJ, ...
        n,m,dimZ,dimY,dimW,Z,Y,W,tau_ref,delta_ref,gamma_ref,Kfb, ...
        regiune,x_eq,pragJ,pragG,pragRez);

    caz.rezultate(2) = identifica_tikhonov('Forma derivativa filtrata + reg. Tikhonov', ...
        caz.FD_filt,caz.FD_val_filt,lambda_grid,Ceq,beq,Cj,j_target,rhoJ, ...
        n,m,dimZ,dimY,dimW,Z,Y,W,tau_ref,delta_ref,gamma_ref,Kfb, ...
        regiune,x_eq,pragJ,pragG,pragRez);

    caz.rezultate(3) = identifica_tikhonov('Forma integrala + reg. Tikhonov', ...
        caz.FD_int,caz.FD_val_int,lambda_grid,Ceq,beq,Cj,j_target,rhoJ, ...
        n,m,dimZ,dimY,dimW,Z,Y,W,tau_ref,delta_ref,gamma_ref,Kfb, ...
        regiune,x_eq,pragJ,pragG,pragRez);

    for i = 1:numel(caz.rezultate)
        fprintf('\nMetoda: %s\n',caz.rezultate(i).nume);
        if contains(caz.rezultate(i).nume,'De Persis')
            fprintf('  rezolvare                = spatiu nul SVD, fara Tikhonov\n');
        else
            fprintf('  lambda ales              = %.4e\n',caz.rezultate(i).lambda);
        end
        if isempty(caz.rezultate(i).coef)
            fprintf('  observatie               = %s\n',caz.rezultate(i).observatie);
            fprintf('  status                   = %s\n',caz.rezultate(i).status);
            continue;
        end
        fprintf('  reziduu identificare     = %.4e\n',caz.rezultate(i).rez_train);
        fprintf('  reziduu validare         = %.4e\n',caz.rezultate(i).rez_val);
        fprintf('  min sigma(J_tau)         = %.4e\n',caz.rezultate(i).minJtau);
        fprintf('  min |gamma|              = %.4e\n',caz.rezultate(i).minGamma);
        fprintf('  ||tau(xe)||             = %.4e\n',caz.rezultate(i).tau_ech_norm);
        fprintf('  |gamma(xe)-1|           = %.4e\n',abs(caz.rezultate(i).gamma_ech-1));
        fprintf('  err J_tau(xe)           = %.4e\n',caz.rezultate(i).errJ_ech);
        fprintf('  status                   = %s\n',caz.rezultate(i).status);
    end

    caz.prag_eroare_finala = max(0.03,5*sigma_x);
    caz.prag_umax = 18;

    caz.sim = struct();
    caz.sim.exact = simuleaza_un_caz(f,g,tau_ref,delta_ref,gamma_ref,Kfb,x0_cl,x_eq,tspan,opts,u_sat,caz.prag_eroare_finala,caz.prag_umax);
    caz.sim.metode = repmat(init_simulare(n,m),numel(caz.rezultate),1);

    fprintf('\nReferinta exacta ||x(T)-xe||=%.4e, max|u|=%.4e, RMS|u|=%.4e\n', ...
        caz.sim.exact.err_finala,caz.sim.exact.u_max,caz.sim.exact.u_rms);

    for i = 1:numel(caz.rezultate)
        if ~strcmp(caz.rezultate(i).status,'OK')
            fprintf('%-20s metoda respinsa numeric: %s\n',caz.rezultate(i).nume,status_afisat(caz.rezultate(i).status));
            continue;
        end

        caz.sim.metode(i) = simuleaza_un_caz(f,g,caz.rezultate(i).model.tau,caz.rezultate(i).model.delta,caz.rezultate(i).model.gamma, ...
            Kfb,x0_cl,x_eq,tspan,opts,u_sat,caz.prag_eroare_finala,caz.prag_umax);

        fprintf('%-20s ||x(T)-xe||=%.4e, max|u|=%.4e, RMS|u|=%.4e, closed-loop=%s\n', ...
            caz.rezultate(i).nume,caz.sim.metode(i).err_finala,caz.sim.metode(i).u_max,caz.sim.metode(i).u_rms,text_ok(caz.sim.metode(i).ok));

        if ~caz.sim.metode(i).ok
            caz.rezultate(i).status = 'RESPINS_CL';
        end
    end

    caz.idx_final = alege_metoda_finala(caz.rezultate,caz.sim);
    if isnan(caz.idx_final)
        caz.metoda_finala = 'nicio metoda acceptata';
        fprintf('\nNu exista metoda finala acceptata pentru cazul %s.\n',nume_caz);
    else
        caz.metoda_finala = caz.rezultate(caz.idx_final).nume;
        fprintf('\nMetoda finala pentru cazul %s: %s\n',nume_caz,caz.metoda_finala);
    end
end

function txt = snr_text(snr_val)
    if isinf(snr_val)
        txt = 'Inf (fara zgomot)';
    else
        txt = sprintf('%.1f dB',snr_val);
    end
end

function [As,scale] = scaleaza_coloane_cu_scara(A)
    scale = vecnorm(A,2,1);
    scale(scale < 1e-12) = 1;
    As = A ./ scale;
end

function Fn = scaleaza_coloane(A)
    sc = vecnorm(A,2,1);
    sc(sc < 1e-12) = 1;
    Fn = A ./ sc;
end

function rez = calc_reziduu_relativ(FD,v)
    if isempty(FD)
        rez = NaN;
    else
        rez = norm(FD*v)/max(sqrt(size(FD,1))*norm(v),1e-12);
    end
end

function e = calc_rmse(A,B)
    d = A(:)-B(:);
    e = sqrt(mean(d.^2));
end

function C = matrice_controlabilitate(A,B)
    n = size(A,1);
    C = B;
    Ap = eye(n);

    for k = 1:n-1
        Ap = Ap*A;
        C = [C Ap*B];
    end
end

function [value,isterminal,direction] = eveniment_oprire(~,x)
    value = norm(x) - 40;
    isterminal = 1;
    direction = 1;
end

function txt = text_ok(cond)
    if cond
        txt = 'OK';
    else
        txt = 'ATENTIE';
    end
end

function print_titlu(txt)
    fprintf('\n--------------------------------------------------------------------\n');
    fprintf('%s\n',upper(txt));
    fprintf('--------------------------------------------------------------------\n');
end

function afiseaza_eroare_dashboard(ME)
    try
        for k = 1:numel(ME.stack)
            fprintf('  eroare dashboard: %s, linia %d\n',ME.stack(k).name,ME.stack(k).line);
        end
    catch
    end
end

function deschide_dashboard_S1(r)

    t_final  = tabel_comparatie_finala(r);
    t_date   = tabel_date_comparatie(r);
    t_fd     = tabel_fd_comparatie(r);
    t_metode = tabel_metode_comparatie(r);
    t_closed = tabel_closed_comparatie(r);

    ecran = get(groot,'ScreenSize');
    latime = min(1750,max(1360,ecran(3)-60));
    inaltime = min(980,max(840,ecran(4)-80));

    fig = uifigure('Name','Sistemul S1', ...
                   'Position',[30 30 latime inaltime], ...
                   'Color',[0.945 0.952 0.965]);
    try
        fig.WindowState = 'maximized';
    catch
    end

    main = uigridlayout(fig,[2 1]);
    main.RowHeight = {90,'1x'};
    main.Padding = [14 10 14 12];
    main.RowSpacing = 8;
    main.BackgroundColor = [0.940 0.948 0.965];

    top = uigridlayout(main,[1 6]);
    top.ColumnWidth = {160,150,150,170,600,'1x'};
    top.Padding = [12 10 12 10];
    top.ColumnSpacing = 12;
    top.BackgroundColor = [0.930 0.940 0.960];

    text_status_ui(top,'Sistem','SISO S1',true);
    text_status_ui(top,'Grad relativ',sprintf('r=%d, n=%d',r.r1,r.n),r.grad_relativ_complet);
    text_status_ui(top,'Zgomot stari',sprintf('SNR %.0f dB',r.snr_x_db),true);
    text_status_ui(top,'Coeficienti',sprintf('%d necunoscute',r.nr_coef),true);

    metoda_curat = metoda_finala_scurta(r.cazuri(1));
    metoda_zgomot = metoda_finala_scurta(r.cazuri(2));
    text_status_ui(top,'Metode finale',"Spațiul nul SVD (date neperturbate) | Forma Integrală + Tikhonov (date zgomotoase)", ...
        ~isnan(r.cazuri(1).idx_final) && ~isnan(r.cazuri(2).idx_final));

    tabs = uitabgroup(main);
    tab1 = uitab(tabs,'Title','1. Rezumat');
    tab2 = uitab(tabs,'Title','2. Date');
    tab3 = uitab(tabs,'Title','3. Matrici F(D)');
    tab4 = uitab(tabs,'Title','4. Identificare');
    tab5 = uitab(tabs,'Title','5. Transformare');
    tab6 = uitab(tabs,'Title','6. Bucla inchisa - validare');

    construieste_tab_rezumat(tab1,r,t_final,t_metode);
    construieste_tab_date(tab2,r,t_date);
    construieste_tab_fd(tab3,r,t_fd);
    construieste_tab_identificare(tab4,r,t_metode);
    construieste_tab_transformare(tab5,r);
    construieste_tab_bucla_inchisa(tab6,r,t_closed);

    tabs.SelectedTab = tab1;
    drawnow;
end

function construieste_tab_rezumat(tab,r,t_final,t_metode)
    lay = uigridlayout(tab,[1 2]);
    lay.ColumnWidth = {640,'1x'};
    lay.Padding = [14 12 14 14];
    lay.ColumnSpacing = 12;
    lay.BackgroundColor = [0.950 0.955 0.965];

    p1 = creeaza_panou(lay,'1.1 Sinteza experimentului');
    g1 = uigridlayout(p1,[3 1]);
   g1.RowHeight = {82,180,'1x'};
    g1.Padding = [10 8 10 10];
    g1.RowSpacing = 8;

txt = "Sistemul este evaluat comparativ in doua regimuri de lucru: cazul ideal, fara zgomot, si cazul practic, cu zgomot pe starile masurate. " + ...
      "Semnalul de intrare u este considerat cunoscut, deoarece reprezinta comanda aplicata procesului si salvata in setul de date. "  + ...
      "Pentru fiecare regim este selectata metoda de identificare care respecta conditiile de validare si conduce la stabilizare in bucla inchisa." + newline + ...
      "Sinteza evidentiaza metodele finale, reziduurile de validare si indicatorii principali folositi in comparatie.";
    adauga_text(g1,txt,12);

    adauga_tabel_rezumat_final(g1,t_final);
    grafic_sinteza_finala(uiaxes(g1),r);

    p2 = creeaza_panou(lay,'1.2 Grafice principale');
    g2 = uigridlayout(p2,[2 2]);
    g2.Padding = [12 10 12 12];
    g2.RowSpacing = 10;
    g2.ColumnSpacing = 10;

    ax1 = uiaxes(g2);
ax1.Layout.Row = 1;
ax1.Layout.Column = 1;
grafic_stari_bucla_inchisa(ax1,r,1);

ax2 = uiaxes(g2);
ax2.Layout.Row = 1;
ax2.Layout.Column = 2;
grafic_stari_bucla_inchisa(ax2,r,2);

ax3 = uiaxes(g2);
ax3.Layout.Row = 2;
ax3.Layout.Column = 1;
grafic_norma_bucla_inchisa(ax3,r);

ax4 = uiaxes(g2);
ax4.Layout.Row = 2;
ax4.Layout.Column = 2;
grafic_rmse_date(ax4,r);
end

function construieste_tab_date(tab,r,t_date)
    lay = uigridlayout(tab,[1 2]);
    lay.ColumnWidth = {620,'1x'};
    lay.Padding = [16 14 16 16];
    lay.ColumnSpacing = 14;
    lay.BackgroundColor = [0.950 0.955 0.965];

    p1 = creeaza_panou(lay,'2.1 Calitatea datelor');

    g1 = uigridlayout(p1,[3 1]);
    g1.RowHeight = {105,250,'1x'};
    g1.Padding = [12 10 12 12];
    g1.RowSpacing = 9;

    txt = "Datele de identificare sunt obtinute prin simularea sistemului neliniar S1 in jurul punctului de echilibru." + newline + ...
          "Pentru cazul ideal se utilizeaza starile curate si derivata exacta furnizata de modelul de simulare." + newline + ...
          "Pentru cazul zgomotos, starile masurate sunt afectate de perturbatii, sunt filtrate, iar derivata este estimata numeric." + newline + ...
          "Aceasta comparatie evidentiaza sensibilitatea metodelor derivative si avantajul formularii integrale.";

    adauga_text(g1,txt,13);

    adauga_tabel_date(g1,t_date);

    grafic_rmse_date(uiaxes(g1),r);

    p2 = creeaza_panou(lay,'2.2 Evolutia starilor, derivatelor estimate si a comenzii');

    g2 = uigridlayout(p2,[1 1]);
    g2.Padding = [10 8 10 10];

    group = uitabgroup(g2);

    for ic = 1:numel(r.cazuri)
        for s = 1:r.n
            caz_txt = lower(char(nume_caz_scurt(r.cazuri(ic).nume)));

            tb = uitab(group,'Title',sprintf('Starea x%d %s',s,caz_txt));

            gl = uigridlayout(tb,[2 1]);
            gl.RowHeight = {'1x','1x'};
            gl.Padding = [10 8 10 10];
            gl.RowSpacing = 10;

            grafic_stare_date(uiaxes(gl),r,ic,s);
            grafic_derivata_date(uiaxes(gl),r,ic,s);
        end
    end

    tb_u = uitab(group,'Title','Comanda u aplicata');

    gl_u = uigridlayout(tb_u,[1 1]);
    gl_u.Padding = [10 8 10 10];

    grafic_intrare(uiaxes(gl_u),r);
end

function construieste_tab_fd(tab,r,t_fd)
    lay = uigridlayout(tab,[1 2]);
    lay.ColumnWidth = {640,'1x'};
    lay.Padding = [16 14 16 16];
    lay.ColumnSpacing = 14;
    lay.BackgroundColor = [0.950 0.955 0.965];

    p1 = creeaza_panou(lay,'3.1 Constructia matricilor F(D)');

g1 = uigridlayout(p1,[3 1]);
g1.RowHeight = {90,305,'1x'};
g1.Padding = [12 8 12 10];
g1.RowSpacing = 7;

txt = "Matricea F(D) sintetizeaza ecuatiile de identificare si le exprima liniar in coeficientii necunoscuti ai matricilor T, N si M." + newline + ...
      "Pentru date curate, F(D) foloseste starile si derivata exacta, oferind un reper de consistenta numerica." + newline + ...
      "Pentru date zgomotoase, sunt comparate forma bruta, forma derivativa filtrata si forma integrala cu regularizare Tikhonov." + newline + ...
      "Reziduul de referinta pe date curate este " + sprintf('%.4e',r.rez_ref_clean) + ".";

adauga_text(g1,txt,12);

adauga_tabel(g1,t_fd, ...
    ["Caz","Matrice","Date","Rol","Randuri","Coloane"], ...
    ["Caz","Matrice","Date","Rol","Randuri","Coloane"], ...
    {105,112,96,'1x',82,82},31,10.4,10.9);

grafic_dimensiuni_fd(uiaxes(g1),r);

p2 = creeaza_panou(lay,'3.2 Diagnostic numeric al matricilor de identificare');

g2 = uigridlayout(p2,[2 2]);
g2.Padding = [12 10 12 12];
g2.RowSpacing = 12;
g2.ColumnSpacing = 12;

grafic_valori_singulare(uiaxes(g2),r.cazuri(1).FD_int, ...
    'Valori singulare ale matricei F(D) integrale - fara zgomot');

grafic_valori_singulare(uiaxes(g2),r.cazuri(2).FD_int, ...
    'Valori singulare ale matricei F(D) integrale - cu zgomot');

ax = uiaxes(g2);
ax.Layout.Row = 2;
ax.Layout.Column = [1 2];

grafic_reziduu_integral(ax,r);
end

function construieste_tab_identificare(tab,r,t_metode)
    lay = uigridlayout(tab,[1 2]);
    lay.ColumnWidth = {760,'1x'};
    lay.Padding = [16 14 16 16];
    lay.ColumnSpacing = 14;
    lay.BackgroundColor = [0.950 0.955 0.965];

    p1 = creeaza_panou(lay,'4.1 Sinteza metodelor de identificare');

    g1 = uigridlayout(p1,[3 1]);
   g1.RowHeight = {82,330,'1x'};
    g1.Padding = [12 10 12 12];
    g1.RowSpacing = 9;

    txt = "Sunt comparate metodele utilizate pentru identificarea functiilor tau(x), delta(x) si gamma(x)." + newline + ...
          "Metoda SVD aplicata direct pe date zgomotoase este pastrata ca test de sensibilitate, in timp ce variantele regularizate sunt analizate din perspectiva reziduurilor, admisibilitatii si validarii in bucla inchisa." + newline + ...
          "Pentru cazul zgomotos, metoda finala este selectata numai daca respecta criteriile numerice si produce un raspuns stabil in bucla inchisa.";

    adauga_text(g1,txt,12);

    adauga_taburi_identificare(g1,t_metode);

    adauga_text(g1,text_identificare_sinteza(r),13.0);

    p2 = creeaza_panou(lay,'4.2 Indicatori de performanta ai metodelor');

    g2 = uigridlayout(p2,[1 1]);
    g2.Padding = [10 8 10 10];

    group = uitabgroup(g2);

    tb1 = uitab(group,'Title','Reziduuri');
    gl1 = uigridlayout(tb1,[1 1]);
    gl1.Padding = [10 8 10 10];

    grafic_reziduuri_identificare_validare(uiaxes(gl1),r);

    tb2 = uitab(group,'Title','Admisibilitate');
    gl2 = uigridlayout(tb2,[2 1]);
    gl2.Padding = [10 8 10 10];
    gl2.RowSpacing = 10;

    grafic_prag_metode(uiaxes(gl2),r,'minJtau',r.prag_minJtau, ...
        'Conditia de difeomorfism: min sigma(J_tau)');

    grafic_prag_metode(uiaxes(gl2),r,'minGamma',r.prag_minGamma, ...
        'Conditia de decuplare: min |gamma(x)|');

    tb3 = uitab(group,'Title','Erori fata de referinta');
    gl3 = uigridlayout(tb3,[2 2]);
    gl3.Padding = [10 8 10 10];
    gl3.RowSpacing = 10;
    gl3.ColumnSpacing = 10;

    grafic_metrici_metode(uiaxes(gl3),r,'err_tau', ...
        'Eroarea relativa a transformarii tau(x)','[%]');

    grafic_metrici_metode(uiaxes(gl3),r,'err_delta', ...
        'Eroarea relativa a termenului delta(x)','[%]');

    grafic_metrici_metode(uiaxes(gl3),r,'err_gamma', ...
        'Eroarea relativa a termenului gamma(x)','[%]');

    grafic_metrici_metode(uiaxes(gl3),r,'err_u', ...
        'Eroarea medie a comenzii de control','u');
end

function construieste_tab_transformare(tab,r)
    lay = uigridlayout(tab,[1 2]);
    lay.ColumnWidth = {640,'1x'};
    lay.Padding = [16 14 16 16];
    lay.ColumnSpacing = 14;
    lay.BackgroundColor = [0.950 0.955 0.965];

    p1 = creeaza_panou(lay,'5.1 Normalizari si verificari de admisibilitate');
    g1 = uigridlayout(p1,[4 1]);
    g1.RowHeight = {88,120,'1x','1x'};
    g1.Padding = [14 11 14 13];
    g1.RowSpacing = 9;

   txt = "Normalizarile impuse fixeaza originea, scara si orientarea locala a transformarii identificate." + newline + ...
      "Conditia tau(xe)=0 fixeaza originea in noile coordonate, iar gamma(xe)=1 elimina ambiguitatea de scala a solutiei." + newline + ...
      "Orientarea locala este controlata prin apropierea jacobianului J_tau(xe) de matricea identitate." + newline + ...
      "Admisibilitatea este verificata prin valoarea singulara minima a lui J_tau si prin valoarea minima a termenului de decuplare |gamma(x)|.";
    adauga_text(g1,txt,12.3);

    adauga_tabel(g1,tabel_transformare_finala(r), ...
        ["Caz","Metoda","Tau_xe","Gamma_xe","Err_Jxe","Min_Jtau","Min_gamma"], ...
        ["Caz","Metoda","tau(xe)","gamma(xe)","err Jxe","min Jtau","min gamma"], ...
        {82,125,70,78,78,78,'1x'},38,10,10.5);

    axNorm = uiaxes(g1);
    grafic_normalizari_transformare(axNorm,r);

    axAdm = uiaxes(g1);
    grafic_admisibilitate_transformare(axAdm,r);

    p2 = creeaza_panou(lay,'5.2 Compararea coeficientilor identificati cu referinta exacta');
    g2 = uigridlayout(p2,[1 1]);
    g2.Padding = [10 8 10 10];

    group_caz = uitabgroup(g2);
    for ic = 1:numel(r.cazuri)
        c = r.cazuri(ic);
        tb_caz = uitab(group_caz,'Title',char(nume_caz_scurt(c.nume)));
        inner = uitabgroup(tb_caz);

        [Tsel,Nsel,Msel,nume] = obtine_model_final_caz(c,r);
        nume_txt = char(strjoin(string(nume)," "));

        tbT = uitab(inner,'Title','Transformare T');
        axT = uiaxes(uigridlayout(tbT,[1 1]));
        grafic_eroare_matrice(axT,r.T_ref,Tsel,sprintf('T - %s',nume_txt),r.numeZ,{'tau1','tau2'});

        tbN = uitab(inner,'Title','Termen neliniar N');
        axN = uiaxes(uigridlayout(tbN,[1 1]));
        grafic_eroare_matrice(axN,r.N_ref,Nsel,sprintf('N - %s',nume_txt),r.numeY,{'delta'});

        tbM = uitab(inner,'Title','Decuplare M');
        axM = uiaxes(uigridlayout(tbM,[1 1]));
        grafic_eroare_matrice(axM,r.M_ref,Msel,sprintf('M - %s',nume_txt),r.numeW,{'gamma'});
    end
end

function grafic_normalizari_transformare(ax,r)
    [labels,tau0,gamma0,errJ] = date_transformare_finale(r);
    cla(ax);

    if isempty(labels)
        text(ax,0.5,0.5,'Nu exista solutii finale pentru acest grafic','HorizontalAlignment','center');
        axis(ax,'off');
        return;
    end

    vals = [max(tau0,1e-14), max(gamma0,1e-14), max(errJ,1e-14)];
    bar(ax,vals,'grouped');
    set(ax,'YScale','log');
    ax.XTick = 1:numel(labels);
    ax.XTickLabel = cellstr(labels);
    title(ax,'Erori ale conditiilor de normalizare');
    xlabel(ax,'Caz analizat');
    ylabel(ax,'Valoare pe scala logaritmica');
    legend(ax,{'||tau(xe)||','|gamma(xe)-1|','||J_tau(xe)-I||'},'Location','best');
    grid(ax,'on');
    stil_axe(ax);
end

function grafic_admisibilitate_transformare(ax,r)
    [labels,~,~,~,minJ,minG] = date_transformare_finale(r);
    cla(ax);

    if isempty(labels)
        text(ax,0.5,0.5,'Nu exista solutii finale pentru acest grafic','HorizontalAlignment','center');
        axis(ax,'off');
        return;
    end

   vals = [max(minJ,1e-14), max(minG,1e-14)];
bar(ax,vals,'grouped');
hold(ax,'on');

prag = r.prag_minJtau;

yline(ax,prag,'--k','prag', ...
    'LineWidth',1.1, ...
    'LabelHorizontalAlignment','right', ...
    'LabelVerticalAlignment','bottom');

set(ax,'YScale','log');

ax.XTick = 1:numel(labels);
ax.XTickLabel = cellstr(labels);

title(ax,'Indicatori de admisibilitate pe regiunea de test');
xlabel(ax,'Caz analizat');
ylabel(ax,'Valoare minima');

legend(ax,{'min sigma(J_tau)','min |gamma|','prag'},'Location','best');

grid(ax,'on');
stil_axe(ax);
end

function [labels,tau0,gamma0,errJ,minJ,minG] = date_transformare_finale(r)
    labels = strings(0,1);
    tau0 = [];
    gamma0 = [];
    errJ = [];
    minJ = [];
    minG = [];

    for ic = 1:numel(r.cazuri)
        c = r.cazuri(ic);
        if isnan(c.idx_final)
            continue;
        end

        rr = c.rezultate(c.idx_final);
        labels(end+1,1) = nume_caz_scurt(c.nume);
        tau0(end+1,1) = abs(rr.tau_ech_norm);
        gamma0(end+1,1) = abs(rr.gamma_ech - 1);
        errJ(end+1,1) = abs(rr.errJ_ech);
        minJ(end+1,1) = rr.minJtau;
        minG(end+1,1) = rr.minGamma;
    end
end

function grafic_sinteza_finala(ax,r)
    labels = strings(0,1);
    rez = [];
    err = [];

    for ic = 1:numel(r.cazuri)
        c = r.cazuri(ic);
        if isnan(c.idx_final)
            continue;
        end
        rr = c.rezultate(c.idx_final);
        ss = c.sim.metode(c.idx_final);
        labels(end+1,1) = nume_caz_scurt(c.nume);
        rez(end+1,1) = max(rr.rez_val,1e-14);
        err(end+1,1) = max(ss.err_finala,1e-14);
    end

    cla(ax);
    if isempty(labels)
        text(ax,0.5,0.5,'Nu exista metode finale selectate','HorizontalAlignment','center');
        axis(ax,'off');
        return;
    end

    bar(ax,[rez err],'grouped');
    set(ax,'YScale','log');
    ax.XTick = 1:numel(labels);
    ax.XTickLabel = cellstr(labels);
    title(ax,'Validarea metodelor finale prin reziduu si eroare in bucla inchisa');
    xlabel(ax,'Caz analizat');
    ylabel(ax,'Valoare pe scala logaritmica');
    legend(ax,{'Reziduu validare','Eroare finala in bucla inchisa'},'Location','best');
    grid(ax,'on');
    stil_axe(ax);
end

function grafic_metrici_bucla_inchisa(ax,r)
    labels = strings(0,1);
    err = [];
    umax = [];
    urms = [];

    s0 = r.cazuri(1).sim.exact;
    if s0.simulat
        labels(end+1,1) = "Referinta";
        err(end+1,1) = max(s0.err_finala,1e-14);
        umax(end+1,1) = max(s0.u_max,1e-14);
        urms(end+1,1) = max(s0.u_rms,1e-14);
    end

    for ic = 1:numel(r.cazuri)
        c = r.cazuri(ic);
        if isnan(c.idx_final)
            continue;
        end
        si = c.sim.metode(c.idx_final);
        if ~isempty(si.T) && si.simulat
            labels(end+1,1) = nume_caz_scurt(c.nume);
            err(end+1,1) = max(si.err_finala,1e-14);
            umax(end+1,1) = max(si.u_max,1e-14);
            urms(end+1,1) = max(si.u_rms,1e-14);
        end
    end

    cla(ax);
    if isempty(labels)
        text(ax,0.5,0.5,'Nu exista simulari in bucla inchisa','HorizontalAlignment','center');
        axis(ax,'off');
        return;
    end

    bar(ax,[err umax urms],'grouped');
    set(ax,'YScale','log');
    ax.XTick = 1:numel(labels);
    ax.XTickLabel = cellstr(labels);
    title(ax,'Indicatori ai raspunsului in bucla inchisa');
    xlabel(ax,'Caz analizat');
    ylabel(ax,'Valoare pe scala logaritmica');
    legend(ax,{'Eroare finala','Max |u|','RMS |u|'},'Location','best');
    grid(ax,'on');
    stil_axe(ax);
end

function txt = text_identificare_sinteza(r)
    buc = strings(0,1);

    c1 = r.cazuri(1);
    c2 = r.cazuri(2);

    if ~isnan(c1.idx_final)
        rr1 = c1.rezultate(c1.idx_final);

        buc(end+1) = sprintf([ ...
            'Pentru cazul ideal, fara zgomot, metoda bazata pe spatiul nul SVD furnizeaza o solutie admisibila. ', ...
            'Pe setul de validare se obtine un reziduu de %.3e, cu min sigma(J_tau)=%.3e si min |gamma(x)|=%.3e.'], ...
            rr1.rez_val, rr1.minJtau, rr1.minGamma);
    else
        buc(end+1) = "Pentru cazul ideal, fara zgomot, nu a fost selectata o metoda finala admisibila.";
    end

    if ~isnan(c2.idx_final)
        rr2 = c2.rezultate(c2.idx_final);

        buc(end+1) = sprintf([ ...
            'Pentru cazul cu zgomot pe stari, metoda selectata este forma integrala cu regularizare Tikhonov. ', ...
            'Aceasta respecta conditiile de admisibilitate si ramane stabila in bucla inchisa, avand reziduu de validare %.3e, ', ...
            'min sigma(J_tau)=%.3e si min |gamma(x)|=%.3e.'], ...
            rr2.rez_val, rr2.minJtau, rr2.minGamma);
    else
        buc(end+1) = "Pentru cazul cu zgomot pe stari, nu a fost selectata o metoda finala admisibila.";
    end

    txt = strjoin(buc,newline);
end
function construieste_tab_bucla_inchisa(tab,r,t_closed)
    lay = uigridlayout(tab,[1 2]);
    lay.ColumnWidth = {650,'1x'};
    lay.Padding = [14 12 14 14];
    lay.ColumnSpacing = 12;
    lay.BackgroundColor = [0.950 0.955 0.965];

    p1 = creeaza_panou(lay,'6.1 Validare in bucla inchisa');
    g1 = uigridlayout(p1,[4 1]);
    g1.RowHeight = {90,175,'1x',92};
    g1.Padding = [12 10 12 12];
    g1.RowSpacing = 9;

   txt = "Validarea in bucla inchisa verifica daca regulatorul obtinut din date stabilizeaza sistemul neliniar." + newline + ...
      "Modelul exact este folosit doar pentru simularea procesului, iar legea de control este construita din functiile identificate." + newline + ...
      "Comparatia urmareste raspunsul starilor, eroarea finala fata de echilibru si nivelul comenzii aplicate." + newline + ...
      "Astfel se verifica daca solutia identificata este utilizabila nu doar numeric, ci si in regim de control.";
    adauga_text(g1,txt,12.5);

    adauga_tabel(g1,t_closed, ...
        ["Caz","Metoda","Err_finala","Max_u","RMS_u","Status"], ...
        ["Caz","Metoda","||x(T)-xe||","max |u|","RMS |u|","Status"], ...
        {110,150,90,80,80,'1x'},35,10.4,10.9);

    grafic_metrici_bucla_inchisa(uiaxes(g1),r);
    adauga_text(g1,text_interpretare_bucla_inchisa(r),13.5);

    p2 = creeaza_panou(lay,'6.2 Raspunsul sistemului controlat');
    g2 = uigridlayout(p2,[1 1]);
    g2.Padding = [10 8 10 10];

    group = uitabgroup(g2);

    tb1 = uitab(group,'Title','Stari');
    gl1 = uigridlayout(tb1,[2 1]);
    gl1.RowHeight = {'1x','1x'};
    gl1.Padding = [10 8 10 10];
    gl1.RowSpacing = 10;
    grafic_stari_bucla_inchisa(uiaxes(gl1),r,1);
    grafic_stari_bucla_inchisa(uiaxes(gl1),r,2);

    tb2 = uitab(group,'Title','Norma erorii');
    grafic_norma_bucla_inchisa(uiaxes(uigridlayout(tb2,[1 1])),r);

    tb3 = uitab(group,'Title','Comanda');
    grafic_comanda_bucla_inchisa(uiaxes(uigridlayout(tb3,[1 1])),r);

    tb4 = uitab(group,'Title','Plan de faza');
    grafic_plan_faza(uiaxes(uigridlayout(tb4,[1 1])),r);
end

function t = tabel_comparatie_finala(r)
    nc = numel(r.cazuri);
    Caz = strings(nc,1);
    SNR = strings(nc,1);
    Sigma_x = strings(nc,1);
    Metoda_finala = strings(nc,1);
    Rez_VAL = strings(nc,1);
    Min_Jtau = strings(nc,1);
    Min_gamma = strings(nc,1);
    Err_finala = strings(nc,1);
    Max_u = strings(nc,1);
    Status = strings(nc,1);

    for i = 1:nc
        c = r.cazuri(i);
        Caz(i) = nume_caz_scurt(c.nume);
        SNR(i) = string(snr_text(c.snr_x_db));
        Sigma_x(i) = format_numar_tabel(c.sigma_x);
        Metoda_finala(i) = metoda_finala_scurta(c);

        if ~isnan(c.idx_final)
            rr = c.rezultate(c.idx_final);
            ss = c.sim.metode(c.idx_final);
            Rez_VAL(i) = format_numar_tabel(rr.rez_val);
            Min_Jtau(i) = format_numar_tabel(rr.minJtau);
            Min_gamma(i) = format_numar_tabel(rr.minGamma);
            Err_finala(i) = format_numar_tabel(ss.err_finala);
            Max_u(i) = format_numar_tabel(ss.u_max);
            Status(i) = status_afisat(rr.status);
        else
            Rez_VAL(i) = "-";
            Min_Jtau(i) = "-";
            Min_gamma(i) = "-";
            Err_finala(i) = "-";
            Max_u(i) = "-";
            Status(i) = "fara metoda";
        end
    end

    t = table(Caz,SNR,Sigma_x,Metoda_finala,Rez_VAL,Min_Jtau,Min_gamma,Err_finala,Max_u,Status);
end

function t = tabel_date_comparatie(r)
    nc = numel(r.cazuri);
    Caz = strings(nc,1);
    SNR = strings(nc,1);
    Sigma_x = strings(nc,1);
    RMSE_x_masurat = strings(nc,1);
    RMSE_x_filtrat = strings(nc,1);
    RMSE_xdot_brut = strings(nc,1);
    RMSE_xdot_folosit = strings(nc,1);
    Observatie = strings(nc,1);

    for i = 1:nc
        c = r.cazuri(i);
        Caz(i) = nume_caz_scurt(c.nume);
        SNR(i) = string(snr_text(c.snr_x_db));
        Sigma_x(i) = format_numar_tabel(c.sigma_x);
        RMSE_x_masurat(i) = format_numar_tabel(c.rmse.x_noisy);
        RMSE_x_filtrat(i) = format_numar_tabel(c.rmse.x_filt);
        RMSE_xdot_brut(i) = format_numar_tabel(c.rmse.xdot_raw);
        RMSE_xdot_folosit(i) = format_numar_tabel(c.rmse.xdot_filt);

        if isinf(c.snr_x_db)
            Observatie(i) = "xdot folosit = exact din simulare";
        else
            Observatie(i) = "xdot folosit = estimat dupa filtrare";
        end
    end

    t = table(Caz,SNR,Sigma_x,RMSE_x_masurat,RMSE_x_filtrat,RMSE_xdot_brut,RMSE_xdot_folosit,Observatie);
end

function t = tabel_fd_comparatie(r)
    Caz = strings(0,1);
    Matrice = strings(0,1);
    Date = strings(0,1);
    Rol = strings(0,1);
    Randuri = strings(0,1);
    Coloane = strings(0,1);

    for ic = 1:numel(r.cazuri)
        c = r.cazuri(ic);

               if ic == 1
            [Caz,Matrice,Date,Rol,Randuri,Coloane] = adauga_linie_fd_simplu( ...
                Caz,Matrice,Date,Rol,Randuri,Coloane, ...
                "Fara zgomot","FD baza","identificare", ...
                "constructie cu x, u si xdot exacte",c.FD_filt);

            [Caz,Matrice,Date,Rol,Randuri,Coloane] = adauga_linie_fd_simplu( ...
                Caz,Matrice,Date,Rol,Randuri,Coloane, ...
                "Fara zgomot","FD baza","validare", ...
                "verificare pe set separat",c.FD_val_filt);
        else
            [Caz,Matrice,Date,Rol,Randuri,Coloane] = adauga_linie_fd_simplu( ...
                Caz,Matrice,Date,Rol,Randuri,Coloane, ...
                "Cu zgomot","FD brut","identificare", ...
                "test SVD pe date zgomotoase",c.FD_raw);

            [Caz,Matrice,Date,Rol,Randuri,Coloane] = adauga_linie_fd_simplu( ...
                Caz,Matrice,Date,Rol,Randuri,Coloane, ...
                "Cu zgomot","FD filtrat","identificare", ...
                "forma derivativa filtrata + Tikhonov",c.FD_filt);

            [Caz,Matrice,Date,Rol,Randuri,Coloane] = adauga_linie_fd_simplu( ...
                Caz,Matrice,Date,Rol,Randuri,Coloane, ...
                "Cu zgomot","FD integral","identificare", ...
                "forma integrala + Tikhonov",c.FD_int);

            [Caz,Matrice,Date,Rol,Randuri,Coloane] = adauga_linie_fd_simplu( ...
                Caz,Matrice,Date,Rol,Randuri,Coloane, ...
                "Cu zgomot","FD brut","validare", ...
                "validare SVD directa",c.FD_val_raw);

            [Caz,Matrice,Date,Rol,Randuri,Coloane] = adauga_linie_fd_simplu( ...
                Caz,Matrice,Date,Rol,Randuri,Coloane, ...
                "Cu zgomot","FD filtrat","validare", ...
                "validare forma derivativa",c.FD_val_filt);

            [Caz,Matrice,Date,Rol,Randuri,Coloane] = adauga_linie_fd_simplu( ...
                Caz,Matrice,Date,Rol,Randuri,Coloane, ...
                "Cu zgomot","FD integral","validare", ...
                "validare forma integrala",c.FD_val_int);
               end
    end

    t = table(Caz,Matrice,Date,Rol,Randuri,Coloane);
end

function [Caz,Matrice,Date,Rol,Randuri,Coloane] = adauga_linie_fd_simplu( ...
    Caz,Matrice,Date,Rol,Randuri,Coloane, ...
    caz_txt,matrice_txt,date_txt,rol_txt,F)

    Caz(end+1,1) = string(caz_txt);
    Matrice(end+1,1) = string(matrice_txt);
    Date(end+1,1) = string(date_txt);
    Rol(end+1,1) = string(rol_txt);

    if isempty(F)
        Randuri(end+1,1) = "-";
        Coloane(end+1,1) = "-";
    else
        Randuri(end+1,1) = string(size(F,1));
        Coloane(end+1,1) = string(size(F,2));
    end
end

function t = tabel_metode_comparatie(r)
    Caz = strings(0,1);
    Metoda = strings(0,1);
    Rol = strings(0,1);
    Lambda = strings(0,1);
    Rez_ID = strings(0,1);
    Rez_VAL = strings(0,1);
    Min_Jtau = strings(0,1);
    Min_gamma = strings(0,1);
    Err_Jxe = strings(0,1);
    Err_finala = strings(0,1);
    Max_u = strings(0,1);
    Observatie = strings(0,1);
    Status = strings(0,1);

    for ic = 1:numel(r.cazuri)
        c = r.cazuri(ic);
        for i = 1:numel(c.rezultate)
            rr = c.rezultate(i);
            Caz(end+1,1) = nume_caz_scurt(c.nume);
            Metoda(end+1,1) = nume_metoda_scurt(rr.nume);
            Rol(end+1,1) = rol_metoda_scurt(c.nume,rr.nume,rr.status);
            Lambda(end+1,1) = format_lambda_tabel(rr.lambda);

            if isfield(rr,'rez_train') && isscalar(rr.rez_train) && isfinite(rr.rez_train)
                Rez_ID(end+1,1) = format_numar_tabel(rr.rez_train);
            else
                Rez_ID(end+1,1) = "-";
            end

            if isfield(rr,'rez_val') && isscalar(rr.rez_val) && isfinite(rr.rez_val)
                Rez_VAL(end+1,1) = format_numar_tabel(rr.rez_val);
            else
                Rez_VAL(end+1,1) = "-";
            end

            if isfield(rr,'minJtau') && isscalar(rr.minJtau) && isfinite(rr.minJtau) && rr.minJtau > 0
                Min_Jtau(end+1,1) = format_numar_tabel(rr.minJtau);
            else
                Min_Jtau(end+1,1) = "-";
            end

            if isfield(rr,'minGamma') && isscalar(rr.minGamma) && isfinite(rr.minGamma) && rr.minGamma > 0
                Min_gamma(end+1,1) = format_numar_tabel(rr.minGamma);
            else
                Min_gamma(end+1,1) = "-";
            end

            if isfield(rr,'errJ_ech') && isscalar(rr.errJ_ech) && isfinite(rr.errJ_ech)
                Err_Jxe(end+1,1) = format_numar_tabel(rr.errJ_ech);
            else
                Err_Jxe(end+1,1) = "-";
            end

            si = c.sim.metode(i);
            if ~isempty(si.T) && si.simulat
                Err_finala(end+1,1) = format_numar_tabel(si.err_finala);
                Max_u(end+1,1) = format_numar_tabel(si.u_max);
            else
                Err_finala(end+1,1) = "-";
                Max_u(end+1,1) = "-";
            end

            Observatie(end+1,1) = observatie_metoda_tabel(rr);
            Status(end+1,1) = status_afisat(rr.status);
        end
    end

    t = table(Caz,Metoda,Rol,Lambda,Rez_ID,Rez_VAL,Min_Jtau,Min_gamma,Err_Jxe,Err_finala,Max_u,Observatie,Status);
end

function txt = observatie_metoda_tabel(rr)
    status = string(rr.status);

    if status == "OK"
        txt = "solutie acceptata";
    elseif status == "RESPINS_SPATIU_NUL"
        txt = "spatiu nul numeric neclar";
    elseif status == "RESPINS_NORMALIZARE"
        txt = "normalizari nerespectate";
    elseif status == "RESPINS_REZIDUU"
        txt = "reziduu de validare peste prag";
    elseif status == "RESPINS_ADMIS"
        txt = "conditii de admisibilitate neindeplinite";
    elseif status == "RESPINS_CL"
        txt = "respins la validarea in bucla inchisa";
    elseif isfield(rr,'observatie') && strlength(string(rr.observatie)) > 0
        txt = string(rr.observatie);
    else
        txt = "metoda respinsa numeric";
    end
end

function t = tabel_transformare_finala(r)
    Caz = strings(0,1);
    Metoda = strings(0,1);
    Tau_xe = strings(0,1);
    Gamma_xe = strings(0,1);
    Err_Jxe = strings(0,1);
    Min_Jtau = strings(0,1);
    Min_gamma = strings(0,1);

    for ic = 1:numel(r.cazuri)
        c = r.cazuri(ic);
        if isnan(c.idx_final)
            continue;
        end
        rr = c.rezultate(c.idx_final);
        Caz(end+1,1) = nume_caz_scurt(c.nume);
        Metoda(end+1,1) = nume_metoda_scurt(rr.nume);
        Tau_xe(end+1,1) = format_numar_tabel(rr.tau_ech_norm);
        Gamma_xe(end+1,1) = format_numar_tabel(rr.gamma_ech);
        Err_Jxe(end+1,1) = format_numar_tabel(rr.errJ_ech);
        Min_Jtau(end+1,1) = format_numar_tabel(rr.minJtau);
        Min_gamma(end+1,1) = format_numar_tabel(rr.minGamma);
    end

    t = table(Caz,Metoda,Tau_xe,Gamma_xe,Err_Jxe,Min_Jtau,Min_gamma);
end

function t = tabel_closed_comparatie(r)
    Caz = strings(0,1);
    Metoda = strings(0,1);
    Err_finala = strings(0,1);
    Max_u = strings(0,1);
    RMS_u = strings(0,1);
    Status = strings(0,1);

    s0 = r.cazuri(1).sim.exact;
    Caz(end+1,1) = "Referinta exacta";
    Metoda(end+1,1) = "exact";
    Err_finala(end+1,1) = format_numar_tabel(s0.err_finala);
    Max_u(end+1,1) = format_numar_tabel(s0.u_max);
    RMS_u(end+1,1) = format_numar_tabel(s0.u_rms);
    Status(end+1,1) = "simulat";

    for ic = 1:numel(r.cazuri)
        c = r.cazuri(ic);
        for i = 1:numel(c.rezultate)
            rr = c.rezultate(i);
            Caz(end+1,1) = nume_caz_scurt(c.nume);
            Metoda(end+1,1) = nume_metoda_scurt(rr.nume);
            si = c.sim.metode(i);
            if ~isempty(si.T) && si.simulat
                Err_finala(end+1,1) = format_numar_tabel(si.err_finala);
                Max_u(end+1,1) = format_numar_tabel(si.u_max);
                RMS_u(end+1,1) = format_numar_tabel(si.u_rms);
            else
                Err_finala(end+1,1) = "-";
                Max_u(end+1,1) = "-";
                RMS_u(end+1,1) = "-";
            end
            Status(end+1,1) = status_afisat(rr.status);
        end
    end

    t = table(Caz,Metoda,Err_finala,Max_u,RMS_u,Status);
end

function [labels,vals,status] = colecteaza_metrici(r,field)
    labels = strings(0,1);
    vals = [];
    status = strings(0,1);

    for ic = 1:numel(r.cazuri)
        c = r.cazuri(ic);
        for i = 1:numel(c.rezultate)
            rr = c.rezultate(i);

            if ~afiseaza_metoda_in_grafic(rr)
                continue;
            end
            if ~isfield(rr,field)
                continue;
            end
            val = rr.(field);
            if isempty(val) || ~isscalar(val) || ~isfinite(val)
                continue;
            end

            labels(end+1,1) = eticheta_metoda_grafic(c.nume,rr.nume);
            vals(end+1,1) = val;
            status(end+1,1) = string(rr.status);
        end
    end
end

function ok = afiseaza_metoda_in_grafic(rr)
    status = string(rr.status);
    if status == "RESPINS_NORMALIZARE" || status == "RESPINS_SPATIU_NUL" || status == "RESPINS_SVD"
        ok = false;
        return;
    end
    ok = isfinite(rr.rez_train) && isfinite(rr.rez_val);
end

function grafic_metrici_metode(ax,r,field,titlu,ylab)
    [labels,vals,status] = colecteaza_metrici(r,field);
    cla(ax);

    if isempty(vals)
        text(ax,0.5,0.5,'Nu exista valori numerice pentru acest grafic','HorizontalAlignment','center');
        axis(ax,'off');
        return;
    end

    b = bar(ax,max(vals,1e-14));
    b.FaceColor = 'flat';

    for i = 1:numel(vals)
        if status(i) == "OK"
            b.CData(i,:) = [0.15 0.40 0.70];
        else
            b.CData(i,:) = [0.70 0.25 0.20];
        end
    end

    vv = vals(isfinite(vals) & vals > 0);
    if ~isempty(vv) && max(vv)/max(min(vv),1e-14) > 30
        set(ax,'YScale','log');
    end

    ax.XTick = 1:numel(labels);
    ax.XTickLabel = cellstr(labels);
    xtickangle(ax,12);
    title(ax,titlu);
    ylabel(ax,ylab);
    grid(ax,'on');
    stil_axe(ax);
end

function grafic_prag_metode(ax,r,field,prag,titlu)
    [labels,vals,~] = colecteaza_metrici(r,field);
    cla(ax);

    if isempty(vals)
        text(ax,0.5,0.5,'Nu exista valori numerice pentru acest grafic','HorizontalAlignment','center');
        axis(ax,'off');
        return;
    end

    hold(ax,'on');
    b = bar(ax,max(vals,1e-14));
    b.FaceColor = 'flat';

    for i = 1:numel(vals)
        if vals(i) >= prag
            b.CData(i,:) = [0.15 0.40 0.70];
        else
            b.CData(i,:) = [0.70 0.25 0.20];
        end
    end

    yline(ax,prag,'--k','prag','LineWidth',1.1);
    set(ax,'YScale','log');
    ax.XTick = 1:numel(labels);
    ax.XTickLabel = cellstr(labels);
    xtickangle(ax,12);
    title(ax,titlu);
    ylabel(ax,'Valoare');
    grid(ax,'on');
    stil_axe(ax);
end

function grafic_rmse_date(ax,r)
    c = r.cazuri(2);

    valori = [c.rmse.x_noisy, c.rmse.x_filt, c.rmse.xdot_raw, c.rmse.xdot_filt];
    etichete = {'x masurat','x filtrat','xdot brut','xdot filtrat'};

    cla(ax);
    b = bar(ax,max(valori,1e-14));
    b.FaceColor = [0.15 0.40 0.70];
    set(ax,'YScale','log');
    ax.XTick = 1:numel(etichete);
    ax.XTickLabel = etichete;
    xtickangle(ax,12);
    title(ax,'Influenta zgomotului asupra starilor si derivatelor estimate');
    ylabel(ax,'RMSE');
    grid(ax,'on');
    stil_axe(ax);
end

function grafic_stare_date(ax,r,ic,s)
    c = r.cazuri(ic);

    idx = 1:min(r.numar_pasi,size(r.x_clean,2));
    t = (0:numel(idx)-1)*r.dt;

    caz_txt = lower(char(nume_caz_scurt(c.nume)));

    cla(ax);
    hold(ax,'on');

    plot(ax,t,r.x_clean(s,idx), ...
        'k-', ...
        'LineWidth',1.4, ...
        'DisplayName','curata');

    plot(ax,t,c.x_masurat(s,idx), ...
        '.', ...
        'MarkerSize',6, ...
        'DisplayName','masurata');

    plot(ax,t,c.x_filt(s,idx), ...
        '--', ...
        'LineWidth',1.2, ...
        'DisplayName','filtrata');

    yline(ax,r.x_echilibru(s), ...
        ':', ...
        'xe', ...
        'HandleVisibility','off');

    title(ax,sprintf('Starea x%d in datele de identificare - %s',s,caz_txt));
    xlabel(ax,'Timp [s]');
    ylabel(ax,sprintf('Starea x%d',s));

    xlim(ax,[t(1) t(end)]);

    grid(ax,'on');
    legend(ax,'Location','best');
    stil_axe(ax);
end

function grafic_derivata_date(ax,r,ic,s)
    c = r.cazuri(ic);

    idx = 1:min(r.numar_pasi,size(r.x_clean,2));
    t = (0:numel(idx)-1)*r.dt;

    caz_txt = lower(char(nume_caz_scurt(c.nume)));

    cla(ax);
    hold(ax,'on');

    plot(ax,t,r.xdot_clean(s,idx), ...
        'k-', ...
        'LineWidth',1.4, ...
        'DisplayName','exacta');

    plot(ax,t,c.xdot_raw(s,idx), ...
        ':', ...
        'LineWidth',1.0, ...
        'DisplayName','bruta');

    plot(ax,t,c.xdot_filt(s,idx), ...
        '--', ...
        'LineWidth',1.2, ...
        'DisplayName','filtrata');

    title(ax,sprintf('Derivata dx%d/dt estimata - %s',s,caz_txt));
    xlabel(ax,'Timp [s]');
    ylabel(ax,sprintf('Derivata dx%d/dt',s));

    xlim(ax,[t(1) t(end)]);

    grid(ax,'on');
    legend(ax,'Location','best');
    stil_axe(ax);
end

function grafic_dimensiuni_fd(ax,r)
    vals = zeros(numel(r.cazuri),3);
    labels = strings(numel(r.cazuri),1);

    for i = 1:numel(r.cazuri)
        c = r.cazuri(i);
        labels(i) = nume_caz_scurt(c.nume);
        vals(i,:) = [size(c.FD_raw,1), size(c.FD_filt,1), size(c.FD_int,1)];
    end

    cla(ax);
    hold(ax,'on');
    bar(ax,vals,'grouped');
    yline(ax,r.nr_coef,'--k',sprintf('coef = %d',r.nr_coef),'LineWidth',1.1);
    ax.XTick = 1:numel(labels);
    ax.XTickLabel = cellstr(labels);
    legend(ax,{'FD brut','FD filtrat','FD integral','coef'},'Location','best');
    title(ax,'Numarul de ecuatii disponibile in matricea F(D)');
    ylabel(ax,'Numar de ecuatii');
    grid(ax,'on');
    stil_axe(ax);
end

function grafic_reziduu_integral(ax,r)
    vals = zeros(numel(r.cazuri),2);
    labels = strings(numel(r.cazuri),1);

    for i = 1:numel(r.cazuri)
        c = r.cazuri(i);
        labels(i) = nume_caz_scurt(c.nume);
        vals(i,:) = [c.rez_ref_int c.rez_ref_val_int];
    end

    cla(ax);
    bar(ax,max(vals,1e-14),'grouped');
    set(ax,'YScale','log');
    ax.XTick = 1:numel(labels);
    ax.XTickLabel = cellstr(labels);
    legend(ax,{'identificare','validare'},'Location','best');
    title(ax,'Consistenta formei integrale pe date de identificare si validare');
    ylabel(ax,'Reziduu normalizat');
    grid(ax,'on');
    stil_axe(ax);
end

function grafic_norma_bucla_inchisa(ax,r)
    cla(ax);
    hold(ax,'on');

    s0 = r.cazuri(1).sim.exact;
    if s0.simulat
        semilogy(ax,s0.T,max(s0.err_traj,1e-14),'k-','LineWidth',1.8,'DisplayName','referinta exacta');
    end

    for ic = 1:numel(r.cazuri)
        c = r.cazuri(ic);
        if isnan(c.idx_final)
            continue;
        end
        si = c.sim.metode(c.idx_final);
        if ~isempty(si.T) && si.simulat
            semilogy(ax,si.T,max(si.err_traj,1e-14),'-','LineWidth',1.4, ...
                'DisplayName',[char(nume_caz_scurt(c.nume)) ' - ' char(nume_metoda_scurt(c.metoda_finala))]);
        end
    end

    title(ax,'Convergenta erorii de stare in bucla inchisa');
    xlabel(ax,'Timp [s]');
    ylabel(ax,'Norma erorii de stare');
    grid(ax,'on');
    legend(ax,'Location','best');
    stil_axe(ax);
end

function grafic_stari_bucla_inchisa(ax,r,s)
    cla(ax);
    hold(ax,'on');

    s0 = r.cazuri(1).sim.exact;
    if s0.simulat
        plot(ax,s0.T,s0.X(s,:),'k-','LineWidth',1.8,'DisplayName','referinta exacta');
    end

    for ic = 1:numel(r.cazuri)
        c = r.cazuri(ic);
        if isnan(c.idx_final)
            continue;
        end
        si = c.sim.metode(c.idx_final);
        if ~isempty(si.T) && si.simulat
            plot(ax,si.T,si.X(s,:),'-','LineWidth',1.35, ...
                'DisplayName',[char(nume_caz_scurt(c.nume)) ' - ' char(nume_metoda_scurt(c.metoda_finala))]);
        end
    end

    yline(ax,r.x_echilibru(s),':','xe','HandleVisibility','off');
    title(ax,sprintf('Evolutia starii x%d in bucla inchisa',s));
    xlabel(ax,'Timp [s]');
    ylabel(ax,sprintf('Stare x%d',s));
    grid(ax,'on');
    legend(ax,'Location','best');
    stil_axe(ax);
end

function grafic_comanda_bucla_inchisa(ax,r)
    cla(ax);
    hold(ax,'on');

    s0 = r.cazuri(1).sim.exact;
    if s0.simulat
        plot(ax,s0.T,s0.U,'k-','LineWidth',1.8,'DisplayName','referinta exacta');
    end

    for ic = 1:numel(r.cazuri)
        c = r.cazuri(ic);
        if isnan(c.idx_final)
            continue;
        end
        si = c.sim.metode(c.idx_final);
        if ~isempty(si.T) && si.simulat
            plot(ax,si.T,si.U,'-','LineWidth',1.35, ...
                'DisplayName',[char(nume_caz_scurt(c.nume)) ' - ' char(nume_metoda_scurt(c.metoda_finala))]);
        end
    end

    yline(ax,r.u_echilibru,'--','u_echilibru','HandleVisibility','off');
    title(ax,'Semnalul de comanda aplicat in bucla inchisa');
    xlabel(ax,'Timp [s]');
    ylabel(ax,'Comanda u');
    grid(ax,'on');
    legend(ax,'Location','best');
    stil_axe(ax);
end

function grafic_plan_faza(ax,r)
    cla(ax);
    hold(ax,'on');

    s0 = r.cazuri(1).sim.exact;
    if s0.simulat
        plot(ax,s0.X(1,:),s0.X(2,:),'k-','LineWidth',1.8,'DisplayName','referinta exacta');
    end

    for ic = 1:numel(r.cazuri)
        c = r.cazuri(ic);
        if isnan(c.idx_final)
            continue;
        end
        si = c.sim.metode(c.idx_final);
        if ~isempty(si.T) && si.simulat
            plot(ax,si.X(1,:),si.X(2,:),'-','LineWidth',1.35, ...
                'DisplayName',[char(nume_caz_scurt(c.nume)) ' - ' char(nume_metoda_scurt(c.metoda_finala))]);
        end
    end

    plot(ax,r.x_echilibru(1),r.x_echilibru(2),'rx','LineWidth',1.8,'MarkerSize',8,'HandleVisibility','off');
    title(ax,'Traiectoria sistemului in planul de faza x1-x2');
    xlabel(ax,'Starea x1');
    ylabel(ax,'Starea x2');
    grid(ax,'on');
    legend(ax,'Location','best');
    stil_axe(ax);
end

function [T,N,M,nume] = obtine_model_final_caz(c,r)
    if isnan(c.idx_final)
        T = r.T_ref;
        N = r.N_ref;
        M = r.M_ref;
        nume = "referinta";
    else
        T = c.rezultate(c.idx_final).T;
        N = c.rezultate(c.idx_final).N;
        M = c.rezultate(c.idx_final).M;
        nume = nume_metoda_scurt(c.rezultate(c.idx_final).nume);
    end
end

function txt = text_interpretare_bucla_inchisa(r)
    buc = strings(0,1);

    c1 = r.cazuri(1);
    c2 = r.cazuri(2);

    if ~isnan(c1.idx_final)
        s1 = c1.sim.metode(c1.idx_final);
        buc(end+1) = sprintf(['Pentru cazul ideal, fara zgomot, metoda selectata reproduce raspunsul de referinta si stabilizeaza sistemul. ', ...
            'Eroarea finala este %.3e, iar valoarea maxima a comenzii este %.3e.'], ...
            s1.err_finala, s1.u_max);
    end

    if ~isnan(c2.idx_final)
        s2 = c2.sim.metode(c2.idx_final);
        buc(end+1) = sprintf(['Pentru cazul cu zgomot pe stari, forma integrala cu regularizare Tikhonov mentine convergenta catre echilibru. ', ...
            'Eroarea finala este %.3e, iar comanda ramane limitata, cu max|u|=%.3e.'], ...
            s2.err_finala, s2.u_max);
    end

    txt = strjoin(buc,newline);
end

function rmet = rol_metoda(nume_metoda,nume_caz)
    rmet = rol_metoda_scurt(nume_caz,nume_metoda,"");
end

function s = status_afisat(status_cod)
    switch string(status_cod)
        case "OK"
            s = "Acceptat";
        case "RESPINS_REZIDUU"
            s = "respins: reziduu mare";
        case "RESPINS_ADMIS"
            s = "respins: gamma/Jtau";
        case "RESPINS_CL"
            s = "respins: bucla inchisa";
        case "RESPINS_SPATIU_NUL"
            s = "respins: fara spatiu nul clar";
        case "RESPINS_NORMALIZARE"
            s = "respins: solutie SVD neadmisibila";
        case "RESPINS_SVD"
            s = "respins: SVD nevalid";
        otherwise
            s = string(status_cod);
    end
end

function grafic_intrare(ax,r)
    idx = 1:min(3*r.numar_pasi,size(r.u_data,2));
    t = (0:numel(idx)-1)*r.dt;

    cla(ax);
    plot(ax,t,r.u_data(1,idx),'LineWidth',1.2);
    hold(ax,'on');

    for it = 1:2
        xline(ax,it*r.numar_pasi*r.dt,':','HandleVisibility','off');
    end

    yline(ax,r.u_echilibru,'--','u_echilibru');
    title(ax,'Semnalul de intrare folosit la generarea datelor');
    xlabel(ax,'Timp [s]');
    ylabel(ax,'Comanda u');
    grid(ax,'on');
    stil_axe(ax);
end

function grafic_valori_singulare(ax,F,titlu)
    s = svd(scaleaza_coloane(F),'econ');

    cla(ax);
    semilogy(ax,max(s,1e-14),'o-','LineWidth',1.1,'MarkerSize',4);
    title(ax,titlu);
    xlabel(ax,'Index valoare singulara');
    ylabel(ax,'Valoare singulara');
    grid(ax,'on');
    stil_axe(ax);
end

function grafic_eroare_matrice(ax,Aref,Aest,titlu,xlabels,ylabels)
    D = Aest - Aref;

    cla(ax);
    imagesc(ax,D);
    cb = colorbar(ax);
    try
        cb.Label.String = 'Diferenta fata de referinta';
    catch
    end

    lim = max(abs(D(:)));
    if ~isfinite(lim) || lim < 1e-12
        lim = 1;
    end
    caxis(ax,[-lim lim]);

    titlu_txt = char(strjoin(string(titlu)," "));
    title(ax,sprintf('Abatere coeficienti %s',titlu_txt));
    xlabel(ax,'Termeni din dictionar');
    ylabel(ax,'Rand coeficient');

    ax.XTick = 1:size(D,2);
    if numel(xlabels) == size(D,2)
        ax.XTickLabel = xlabels;
    end

    ax.YTick = 1:size(D,1);
    if numel(ylabels) == size(D,1)
        ax.YTickLabel = ylabels;
    end

    xtickangle(ax,45);
    try
        ax.FontSize = 9.5;
    catch
    end
    stil_axe(ax);
end

function txt = nume_caz_scurt(nume_caz)
    caz = string(nume_caz);
    if contains(caz,"Fara zgomot","IgnoreCase",true)
        txt = "Fara zgomot";
    elseif contains(caz,"Cu zgomot","IgnoreCase",true)
        txt = "Cu zgomot";
    else
        txt = caz;
    end
end

function txt = nume_metoda_scurt(nume_metoda)
    metoda = string(nume_metoda);
    if contains(metoda,"spatiu nul","IgnoreCase",true)
        txt = "Metoda de baza-SVD";
    elseif contains(metoda,"direct","IgnoreCase",true)
        txt = "Metoda de baza directa";
    elseif contains(metoda,"derivativa","IgnoreCase",true)
        txt = "Derivativ + Tikhonov";
    elseif contains(metoda,"integrala","IgnoreCase",true)
        txt = "Integral + Tikhonov";
    else
        txt = metoda;
    end
end

function txt = eticheta_metoda_grafic(nume_caz,nume_metoda)
    caz = string(nume_caz);
    metoda = string(nume_metoda);
    if contains(caz,"Fara zgomot","IgnoreCase",true)
        txt = "Date curate: SVD";
    elseif contains(metoda,"derivativa","IgnoreCase",true)
        txt = "Zgomot: derivativ + Tik.";
    elseif contains(metoda,"integrala","IgnoreCase",true)
        txt = "Zgomot: integral + Tik.";
    elseif contains(metoda,"direct","IgnoreCase",true)
        txt = "Zgomot: SVD direct";
    else
        txt = nume_caz_scurt(caz) + ": " + nume_metoda_scurt(metoda);
    end
end

function txt = metoda_finala_scurta(c)
    if isnan(c.idx_final)
        txt = "fara metoda";
    else
        txt = nume_metoda_scurt(c.metoda_finala);
    end
end

function txt = rol_metoda_scurt(nume_caz,nume_metoda,status)
    caz = string(nume_caz);
    metoda = string(nume_metoda);
    status = string(status);

    if contains(caz,"Fara zgomot","IgnoreCase",true)
        txt = "caz ideal";
    elseif contains(metoda,"direct","IgnoreCase",true)
        txt = "test sensibilitate";
    elseif contains(metoda,"derivativa","IgnoreCase",true)
        txt = "varianta intermediara";
    elseif contains(metoda,"integrala","IgnoreCase",true)
        txt = "varianta propusa";
    elseif status ~= "OK"
        txt = "metoda respinsa";
    else
        txt = "-";
    end
end

function ok = are_solutie_pentru_tabel(rr)
    status = string(rr.status);
    if status == "RESPINS_NORMALIZARE" || status == "RESPINS_SPATIU_NUL" || status == "RESPINS_SVD"
        ok = false;
        return;
    end
    ok = isfinite(rr.rez_train) && isfinite(rr.rez_val);
end

function txt = format_lambda_tabel(x)
    if isempty(x) || ~isscalar(x) || ~isfinite(x) || abs(x) < 1e-14
        txt = "0";
    else
        txt = string(sprintf('%.4e',x));
    end
end

function txt = format_numar_tabel(x)
    if isempty(x) || ~isscalar(x) || ~isfinite(x)
        txt = "-";
    elseif abs(x) < 1e-12
        txt = "0";
    else
        txt = string(sprintf('%.4e',x));
    end
end

function adauga_tabel_date(parent,t_date)

    fontHeaderDate = 12.5;
    fontBodyDate   = 12.0;

    if height(t_date) < 2
        adauga_tabel(parent,t_date, ...
            ["Caz","SNR","Sigma_x","RMSE_x_masurat","RMSE_x_filtrat","RMSE_xdot_brut","RMSE_xdot_folosit","Observatie"], ...
            ["Caz","SNR","sigma_x","RMSE x","RMSE x filt.","RMSE xdot","RMSE xdot filt.","Observatie"], ...
            {75,70,75,82,82,82,90,'1x'},34,11.0,11.5);
        return;
    end

    box = uipanel(parent);
    box.BackgroundColor = [1 1 1];
    box.BorderType = 'line';

    indicatori = ["SNR"; ...
                  "sigma_x"; ...
                  "RMSE x masurat"; ...
                  "RMSE x filtrat"; ...
                  "RMSE xdot brut"; ...
                  "RMSE xdot folosit"; ...
                  "Observatie"];

    vals1 = [string(t_date.SNR(1)); ...
             string(t_date.Sigma_x(1)); ...
             string(t_date.RMSE_x_masurat(1)); ...
             string(t_date.RMSE_x_filtrat(1)); ...
             string(t_date.RMSE_xdot_brut(1)); ...
             string(t_date.RMSE_xdot_folosit(1)); ...
             string(t_date.Observatie(1))];

    vals2 = [string(t_date.SNR(2)); ...
             string(t_date.Sigma_x(2)); ...
             string(t_date.RMSE_x_masurat(2)); ...
             string(t_date.RMSE_x_filtrat(2)); ...
             string(t_date.RMSE_xdot_brut(2)); ...
             string(t_date.RMSE_xdot_folosit(2)); ...
             string(t_date.Observatie(2))];

    g = uigridlayout(box,[numel(indicatori)+1 3]);

    g.Padding = [4 3 4 3];
    g.RowSpacing = 1;
    g.ColumnSpacing = 1;

    g.RowHeight = [{30}, repmat({29},1,numel(indicatori)-1), {48}];

    g.ColumnWidth = {185,'1x','1x'};

    g.BackgroundColor = [0.82 0.84 0.88];

    headers = ["Indicator", string(t_date.Caz(1)), string(t_date.Caz(2))];

    for j = 1:3
        adauga_celula_tabel(g,1,j,headers(j),true,[0.90 0.93 0.98],fontHeaderDate);
    end

    for i = 1:numel(indicatori)
        bg = [1 1 1];

        if mod(i,2) == 0
            bg = [0.965 0.970 0.982];
        end

        adauga_celula_tabel(g,i+1,1,indicatori(i),true,bg,fontBodyDate);
        adauga_celula_tabel(g,i+1,2,vals1(i),false,bg,fontBodyDate);
        adauga_celula_tabel(g,i+1,3,vals2(i),false,bg,fontBodyDate);
    end
end

function adauga_taburi_identificare(parent,t_metode)
    group = uitabgroup(parent);
tb1 = uitab(group,'Title','Rol si status');
g1 = uigridlayout(tb1,[1 1]);
g1.Padding = [8 8 8 8];

adauga_tabel(g1,t_metode, ...
    ["Caz","Metoda","Rol","Lambda","Observatie","Status"], ...
    ["Caz","Metoda","Rol","lambda","Observatie","Status"], ...
    {95,210,150,80,210,'1x'},58,10.2,10.8);

    tb2 = uitab(group,'Title','Reziduuri');
    g2 = uigridlayout(tb2,[1 1]);
    g2.Padding = [8 8 8 8];

    adauga_tabel(g2,t_metode, ...
        ["Caz","Metoda","Rez_ID","Rez_VAL"], ...
        ["Caz","Metoda","Rez. identificare","Rez. validare"], ...
        {110,'1x',155,155},48,10.6,11);

    tb3 = uitab(group,'Title','Admisibilitate');
    g3 = uigridlayout(tb3,[1 1]);
    g3.Padding = [8 8 8 8];

    adauga_tabel(g3,t_metode, ...
        ["Caz","Metoda","Min_Jtau","Min_gamma","Err_Jxe"], ...
        ["Caz","Metoda","min sigma(J_tau)","min |gamma|","err J_tau(xe)"], ...
        {110,'1x',155,145,145},48,10.6,11);

    tb4 = uitab(group,'Title','Bucla inchisa');
    g4 = uigridlayout(tb4,[1 1]);
    g4.Padding = [8 8 8 8];

    adauga_tabel(g4,t_metode, ...
        ["Caz","Metoda","Err_finala","Max_u","Status"], ...
        ["Caz","Metoda","||x(T)-xe||","max |u|","Status"], ...
        {110,'1x',150,130,170},48,10.6,11);
end

function adauga_tabel(parent,t,cols,headers,widths,rowHeight,fontBody,fontHeader)
    box = uipanel(parent);
    box.BackgroundColor = [1 1 1];
    box.BorderType = 'line';

    nr = height(t);
    nc = numel(cols);

    g = uigridlayout(box,[nr+1 nc]);
    g.Padding = [5 4 5 4];
    g.RowSpacing = 1;
    g.ColumnSpacing = 1;

    g.RowHeight = [{28}, repmat({rowHeight},1,nr)];
    g.ColumnWidth = widths;

    g.BackgroundColor = [1 1 1];

    for j = 1:nc
        adauga_celula_tabel(g,1,j,string(headers(j)),true,[0.90 0.93 0.98],fontHeader);
    end

    for i = 1:nr
        bg = [1 1 1];

        if mod(i,2) == 0
            bg = [0.965 0.970 0.982];
        end

        for j = 1:nc
            nume_col = char(cols(j));
            val = string(t.(nume_col)(i));

            boldCell = (j == 1 || j == 2);

            adauga_celula_tabel(g,i+1,j,val,boldCell,bg,fontBody);
        end
    end
end

function p = creeaza_panou(parent,titlu)
    p = uipanel(parent,'Title',titlu);
    p.BackgroundColor = [1 1 1];
    p.BorderType = 'line';
    p.FontName = 'Arial';
    p.FontWeight = 'bold';
    p.FontSize = 12;
end

function text_status_ui(parent,titlu,valoare,ok)
    card = uipanel(parent);
    card.BackgroundColor = [1 1 1];
    card.BorderType = 'line';

    gl = uigridlayout(card,[2 1]);
    gl.RowHeight = {24,'1x'};
    gl.Padding = [12 7 12 8];
    gl.RowSpacing = 2;
    gl.BackgroundColor = [1 1 1];

    l1 = uilabel(gl);
    l1.Text = string(titlu);
    l1.FontName = 'Arial';
    l1.FontSize = 11.2;
    l1.FontColor = [0.34 0.36 0.43];
    l1.HorizontalAlignment = 'left';
    l1.VerticalAlignment = 'bottom';

    l2 = uilabel(gl);
    l2.Text = string(valoare);
    l2.FontName = 'Arial';
    l2.FontSize = 13.0;
    l2.FontWeight = 'bold';
    l2.WordWrap = 'on';
    l2.Interpreter = 'none';
    l2.HorizontalAlignment = 'left';
    l2.VerticalAlignment = 'top';

    if ok
        l2.FontColor = [0.02 0.22 0.46];
    else
        l2.FontColor = [0.65 0.12 0.10];
    end
end

function h = adauga_text(parent,txt,fs)
    h = uitextarea(parent);
    txt = replace(string(txt),"\n",newline);
    h.Value = cellstr(splitlines(txt));
    h.Editable = 'off';
    h.FontName = 'Arial';
    h.FontSize = fs;
    h.BackgroundColor = [1 1 1];
end

function stil_tabel(uit)
    try
        uit.FontName = 'Arial';
        uit.FontSize = 11;
        uit.RowStriping = 'on';
        uit.ColumnSortable = true;
    catch
    end
end

function stil_axe(ax)
    try
        ax.FontName = 'Arial';
        ax.FontSize = 10.5;
        ax.Box = 'on';
        ax.GridAlpha = 0.18;
    catch
    end
end
function adauga_tabel_rezumat_final(parent,t_final)
    box = uipanel(parent);
    box.BackgroundColor = [1 1 1];
    box.BorderType = 'line';

    g = uigridlayout(box,[3 6]);
    g.Padding = [4 3 4 3];
    g.RowSpacing = 1;
    g.ColumnSpacing = 1;

    g.RowHeight = {30,74,74};
    g.ColumnWidth = {78,68,70,118,150,'1x'};

    g.BackgroundColor = [0.82 0.84 0.88];

    fontHeaderFinal = 12.5;
    fontBodyFinal   = 11.4;

    headers = ["Caz","SNR","sigma_x","Metoda finala","Validare","Bucla inchisa"];

    for j = 1:numel(headers)
        adauga_celula_tabel(g,1,j,headers(j),true,[0.90 0.93 0.98],fontHeaderFinal);
    end

    idx_fara = find(contains(string(t_final.Caz),"Fara","IgnoreCase",true),1);
    idx_zgom = find(contains(string(t_final.Caz),"Cu","IgnoreCase",true),1);

    adauga_linie_rezumat_final(g,2,t_final,idx_fara,fontBodyFinal);
    adauga_linie_rezumat_final(g,3,t_final,idx_zgom,fontBodyFinal);
end

function adauga_linie_rezumat_final(g,row,t,idx,fontBodyFinal)
    if isempty(idx)
        valori = ["-","-","-","-","-","-"];
    else
        valori = strings(1,6);

        valori(1) = string(t.Caz(idx));
        valori(2) = string(t.SNR(idx));
        valori(3) = string(t.Sigma_x(idx));
        valori(4) = string(t.Metoda_finala(idx));

        valori(5) = "Rez. val. = " + string(t.Rez_VAL(idx)) + newline + ...
                    "min Jtau = "  + string(t.Min_Jtau(idx)) + newline + ...
                    "min gamma = " + string(t.Min_gamma(idx));

        valori(6) = "||x(T)-xe|| = " + string(t.Err_finala(idx)) + newline + ...
                    "max |u| = "     + string(t.Max_u(idx));
    end

    for j = 1:6

        boldCell = (j == 1);

        adauga_celula_tabel(g,row,j,valori(j),boldCell,[1 1 1],fontBodyFinal);
    end
end

function adauga_tabel_rezumat_metode(parent,t_metode)
    box = uipanel(parent);
    box.BackgroundColor = [1 1 1];
    box.BorderType = 'line';

    nr = height(t_metode);

    g = uigridlayout(box,[nr+1 6]);
    g.Padding = [4 3 4 3];
    g.RowSpacing = 1;
    g.ColumnSpacing = 1;
   g.RowHeight = [{24}, repmat({54},1,nr)];
g.ColumnWidth = {108,78,104,104,66,'1x'};
    g.BackgroundColor = [0.82 0.84 0.88];

    headers = ["Metoda","Rol","Reziduuri","Admisib.","CL","Observatie"];
    for j = 1:numel(headers)
        adauga_celula_tabel(g,1,j,headers(j),true,[0.90 0.93 0.98],9.4);
    end

    for i = 1:nr
        metoda_txt = string(t_metode.Caz(i)) + newline + string(t_metode.Metoda(i));

        rol_txt = string(t_metode.Rol(i)) + newline + ...
                  "lambda " + string(t_metode.Lambda(i));

        rez_txt = "ID " + string(t_metode.Rez_ID(i)) + newline + ...
                  "VAL " + string(t_metode.Rez_VAL(i));

        adm_txt = "minJ " + string(t_metode.Min_Jtau(i)) + newline + ...
                  "minG " + string(t_metode.Min_gamma(i)) + newline + ...
                  "Jxe " + string(t_metode.Err_Jxe(i));

        cl_txt = "xT " + string(t_metode.Err_finala(i)) + newline + ...
                 "u " + string(t_metode.Max_u(i));

        obs_txt = string(t_metode.Observatie(i)) + newline + string(t_metode.Status(i));

        bg = [1 1 1];
        if mod(i,2) == 0
            bg = [0.965 0.970 0.982];
        end

      adauga_celula_tabel(g,i+1,1,metoda_txt,false,bg,8.7);
adauga_celula_tabel(g,i+1,2,rol_txt,false,bg,8.7);
adauga_celula_tabel(g,i+1,3,rez_txt,false,bg,8.7);
adauga_celula_tabel(g,i+1,4,adm_txt,false,bg,8.7);
adauga_celula_tabel(g,i+1,5,cl_txt,false,bg,8.7);
adauga_celula_tabel(g,i+1,6,obs_txt,false,bg,8.5);
    end
end

function h = adauga_celula_tabel(parent,row,col,txt,bold,bg,fs)
    h = uilabel(parent);
    h.Layout.Row = row;
    h.Layout.Column = col;
    h.Text = string(txt);
    h.WordWrap = 'on';
    h.FontName = 'Arial';
    h.FontSize = fs;
    h.BackgroundColor = bg;
    h.FontColor = [0.05 0.06 0.08];
    h.HorizontalAlignment = 'left';
    h.VerticalAlignment = 'center';

    if bold
        h.FontWeight = 'bold';
    else
        h.FontWeight = 'normal';
    end
end

function grafic_reziduuri_identificare_validare(ax,r)
    labels = strings(0,1);
    rez_id = [];
    rez_val = [];

    for ic = 1:numel(r.cazuri)
        c = r.cazuri(ic);

        for i = 1:numel(c.rezultate)
            rr = c.rezultate(i);

            are_id  = isfield(rr,'rez_train') && isscalar(rr.rez_train) && isfinite(rr.rez_train);
            are_val = isfield(rr,'rez_val')   && isscalar(rr.rez_val)   && isfinite(rr.rez_val);

            if ~(are_id || are_val)
                continue;
            end

            labels(end+1,1) = nume_caz_scurt(c.nume) + ": " + nume_metoda_scurt(rr.nume);

            if are_id
                rez_id(end+1,1) = max(rr.rez_train,1e-14);
            else
                rez_id(end+1,1) = NaN;
            end

            if are_val
                rez_val(end+1,1) = max(rr.rez_val,1e-14);
            else
                rez_val(end+1,1) = NaN;
            end
        end
    end

    cla(ax);

    if isempty(labels)
        text(ax,0.5,0.5,'Nu exista reziduuri disponibile pentru afisare', ...
            'HorizontalAlignment','center');
        axis(ax,'off');
        return;
    end

    vals = [rez_id rez_val];

    bar(ax,vals,'grouped');

    set(ax,'YScale','log');

    ax.XTick = 1:numel(labels);
    ax.XTickLabel = cellstr(labels);
    ax.XTickLabelRotation = 15;

    title(ax,'Compararea reziduurilor pe identificare si validare');
    xlabel(ax,'Metoda analizata');
    ylabel(ax,'Reziduu normalizat');

    legend(ax,{'Identificare','Validare'},'Location','best');

    grid(ax,'on');
    stil_axe(ax);
end
