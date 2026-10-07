clear; clc; close all;

set(groot,'defaultTextInterpreter','none');
set(groot,'defaultAxesTickLabelInterpreter','none');
set(groot,'defaultLegendInterpreter','none');

rng(42,'twister');

a13 = 0.18;
a24 = 0.12;
a31 = 0.16;
a42 = 0.10;

tau_exact = @(x) [x(1) + a13*sin(x(3));
                  x(2) + a24*sin(x(4));
                  x(3) + a31*sin(x(1));
                  x(4) + a42*sin(x(2))];

Jtau_exact = @(x) [1,             0,             a13*cos(x(3)), 0;
                   0,             1,             0,             a24*cos(x(4));
                   a31*cos(x(1)), 0,             1,             0;
                   0,             a42*cos(x(2)), 0,             1];

delta_exact = @(x) [-0.70*sin(x(1)) + 0.20*x(3) - 0.15*x(2)^3 + 0.10*cos(x(1)+x(3)) + 0.05*x(1)*x(4);
                     0.15*x(1) - 0.60*sin(x(3)) - 0.12*x(4)^3 + 0.10*cos(x(2)+x(4)) + 0.04*x(2)*x(3)];

gamma_exact = @(x) [1.20 + 0.12*cos(x(1)) - 0.08*sin(x(3)),  0.25 + 0.05*sin(x(2));
                   -0.20 + 0.04*sin(x(4)),                   1.10 + 0.10*cos(x(3)) - 0.06*sin(x(1))];

n = 4;
m = 2;

f = @(x) camp_vectorial_mimo(x,Jtau_exact,tau_exact,delta_exact);
g = @(x) camp_intrare_mimo(x,Jtau_exact,gamma_exact);

gamma0_exact = gamma_exact(zeros(n,1));

r = [2 2];
grad_relativ_complet = sum(r) == n;

A2 = [0 1; 0 0];
B2 = [0; 1];
Ac = blkdiag(A2,A2);
Bc = blkdiag(B2,B2);

Q_lqr = diag([8 1.3 7 1.2]);
R_lqr = diag([1.8 2.0]);

try
    K_feedback = -lqr(Ac,Bc,Q_lqr,R_lqr);
catch
    K_feedback = [-3.2 -2.5  0    0;
                   0    0   -3.0 -2.4];
end

poli_lqr = eig(Ac + Bc*K_feedback);
rang_ctrb = rank(matrice_controlabilitate(Ac,Bc));

Z = @(x) [x(1);
          x(2);
          x(3);
          x(4);
          sin(x(1));
          sin(x(2));
          sin(x(3));
          sin(x(4));
          cos(x(1));
          cos(x(2));
          cos(x(3));
          cos(x(4))];

Y = @(x) [1;
          x(1);
          x(2);
          x(3);
          x(4);
          sin(x(1));
          sin(x(2));
          sin(x(3));
          sin(x(4));
          cos(x(1)+x(3));
          cos(x(2)+x(4));
          x(2)^3;
          x(4)^3;
          x(1)*x(4);
          x(2)*x(3)];

W = @(x) [1,         0;
          0,         1;
          cos(x(1)), 0;
          sin(x(3)), 0;
          sin(x(4)), 0;
          0,         sin(x(2));
          0,         cos(x(3));
          0,         sin(x(1))];

dimZ = numel(Z(zeros(n,1)));
dimY = numel(Y(zeros(n,1)));
W0 = W(zeros(n,1));
dimW = size(W0,1);

numeZ = {'x1','x2','x3','x4','sin x1','sin x2','sin x3','sin x4','cos x1','cos x2','cos x3','cos x4'};
numeY = {'1','x1','x2','x3','x4','sin x1','sin x2','sin x3','sin x4','cos(x1+x3)','cos(x2+x4)','x2^3','x4^3','x1*x4','x2*x3'};
numeW = {'1 pe u1','1 pe u2','cos x1 pe u1','sin x3 pe u1','sin x4 pe u1','sin x2 pe u2','cos x3 pe u2','sin x1 pe u2'};

nr_coef = n*dimZ + m*dimY + m*dimW;

T_ideal = zeros(n,dimZ);
T_ideal(1,1) = 1;      T_ideal(1,7) = a13;
T_ideal(2,2) = 1;      T_ideal(2,8) = a24;
T_ideal(3,3) = 1;      T_ideal(3,5) = a31;
T_ideal(4,4) = 1;      T_ideal(4,6) = a42;

N_ideal = zeros(m,dimY);
N_ideal(1,6)  = -0.70;
N_ideal(1,4)  =  0.20;
N_ideal(1,12) = -0.15;
N_ideal(1,10) =  0.10;
N_ideal(1,14) =  0.05;
N_ideal(2,2)  =  0.15;
N_ideal(2,8)  = -0.60;
N_ideal(2,13) = -0.12;
N_ideal(2,11) =  0.10;
N_ideal(2,15) =  0.04;

M_ideal = zeros(m,dimW);
M_ideal(1,1) =  1.20;
M_ideal(1,3) =  0.12;
M_ideal(1,4) = -0.08;
M_ideal(1,2) =  0.25;
M_ideal(1,6) =  0.05;
M_ideal(2,1) = -0.20;
M_ideal(2,5) =  0.04;
M_ideal(2,2) =  1.10;
M_ideal(2,7) =  0.10;
M_ideal(2,8) = -0.06;

dt = 0.04;
numar_traiectorii = 32;
numar_pasi = 165;
numar_traiectorii_validare = 8;
numar_pasi_validare = 145;
amplitudine_intrare = 0.75;

snr_x_db = 25;

fereastra_filtrare = 9;
ordin_polinom = 3;
margine = floor(fereastra_filtrare/2);
pas_integral = 5;

lambda_grid = logspace(-9,0,75);

rho_jacobian0 = 1e-1;

trace_gamma0_tinta = m;

tau0_tinta = zeros(n,1);

prag_Jtau_dorit = 5e-2;
prag_gamma_dorit = 5e-2;
prag_Jtau_minim = 3e-2;
prag_gamma_minim = 3e-2;
prag_err_jacobian0 = 0.80;

prag_reziduu_validare = 0.10;

regiune_test = [-0.65 0.65;
                -0.65 0.65;
                -0.65 0.65;
                -0.65 0.65];

u_actuator_max = 20;
u_siguranta = 0.90*u_actuator_max;

x0_cl = [0.16; -0.08; -0.14; 0.10];
t_cl = linspace(0,40,4001);

print_titlu('verificari initiale');

fprintf('MIMO 2x2 complex cu sin/cos si gamma dependenta de stare.\n');
fprintf('Date pentru identificare: x masurat cu zgomot, u curat/salvat.\n');
fprintf('Numar coeficienti necunoscuti = %d.\n',nr_coef);
fprintf('Grad relativ complet: r1+r2=%d, n=%d -> %s\n',sum(r),n,text_ok(grad_relativ_complet));
fprintf('Gamma exacta la origine: sigma_min=%.4e, det=%.4e\n',min(svd(gamma0_exact)),det(gamma0_exact));
fprintf('Rang controlabilitate(Ac,Bc)=%d din %d.\n',rang_ctrb,n);
fprintf('Poli forma liniarizata: ');
disp(poli_lqr.');
fprintf('Normalizare folosita: trace(gamma(0)) = %.2f si tau(0)=0\n',trace_gamma0_tinta);
fprintf('Regularizare locala: J_tau(0) aproximativ I, rho=%.2e, folosita ca orientare\n',rho_jacobian0);
fprintf('Prag reziduu validare pentru status OK = %.2e\n',prag_reziduu_validare);

print_titlu('colectare date');

[x_clean,u_data,xdot_clean] = colecteaza_date_simulare(dt,numar_traiectorii,numar_pasi,f,g,n,m,amplitudine_intrare);
[x_val,u_val,xdot_val] = colecteaza_date_simulare(dt,numar_traiectorii_validare,numar_pasi_validare,f,g,n,m,amplitudine_intrare);

putere_medie_x = mean(sum(x_clean.^2,1));
sigma_x = sqrt(putere_medie_x/(n*10^(snr_x_db/10)));

rng(123,'twister');
x_noisy = x_clean + sigma_x*randn(size(x_clean));

rng(321,'twister');
x_val_noisy = x_val + sigma_x*randn(size(x_val));

[x_raw,xdot_raw] = estimeaza_brut(x_noisy,dt,numar_traiectorii,numar_pasi);
[x_filt,xdot_filt] = estimeaza_polinom_local(x_noisy,dt,numar_traiectorii,numar_pasi,fereastra_filtrare,ordin_polinom);
idx_interior = construieste_index_interior(numar_traiectorii,numar_pasi,margine);

[x_val_raw,xdot_val_raw] = estimeaza_brut(x_val_noisy,dt,numar_traiectorii_validare,numar_pasi_validare);
[x_val_filt,xdot_val_filt] = estimeaza_polinom_local(x_val_noisy,dt,numar_traiectorii_validare,numar_pasi_validare,fereastra_filtrare,ordin_polinom);
idx_val_interior = construieste_index_interior(numar_traiectorii_validare,numar_pasi_validare,margine);

fprintf('Puncte identificare: %d\n',size(x_clean,2));
fprintf('Puncte validare separata: %d\n',size(x_val,2));
fprintf('SNR x = %.1f dB, sigma zgomot x = %.4e\n',snr_x_db,sigma_x);
fprintf('u este comanda aplicata si salvata, deci este folosita fara zgomot.\n');
fprintf('RMSE x zgomotos = %.4e\n',calc_rmse(x_noisy,x_clean));
fprintf('RMSE x filtrat  = %.4e\n',calc_rmse(x_filt(:,idx_interior),x_clean(:,idx_interior)));
fprintf('RMSE xdot brut  = %.4e\n',calc_rmse(xdot_raw(:,idx_interior),xdot_clean(:,idx_interior)));
fprintf('RMSE xdot filt  = %.4e\n',calc_rmse(xdot_filt(:,idx_interior),xdot_clean(:,idx_interior)));

prag_norma_finala = max(0.10*norm(x0_cl),5*sigma_x*sqrt(n));
fprintf('Prag closed-loop ||x(T)|| acceptat = %.4e\n',prag_norma_finala);
fprintf('Prag closed-loop max||u|| acceptat = %.4e\n',u_siguranta);

print_titlu('constructie F(D)');

FD_raw = construieste_FD_derivativ(x_raw,u_data,xdot_raw,Z,Y,W,Ac,Bc,n,m,dimZ,dimY,dimW);
FD_filt = construieste_FD_derivativ(x_filt(:,idx_interior),u_data(:,idx_interior),xdot_filt(:,idx_interior),Z,Y,W,Ac,Bc,n,m,dimZ,dimY,dimW);
FD_int = construieste_FD_integral(x_filt,u_data,Z,Y,W,Ac,Bc,n,m,dimZ,dimY,dimW,dt,numar_traiectorii,numar_pasi,pas_integral,margine);

FD_val_filt = construieste_FD_derivativ(x_val_filt(:,idx_val_interior),u_val(:,idx_val_interior),xdot_val_filt(:,idx_val_interior),Z,Y,W,Ac,Bc,n,m,dimZ,dimY,dimW);
FD_val_int = construieste_FD_integral(x_val_filt,u_val,Z,Y,W,Ac,Bc,n,m,dimZ,dimY,dimW,dt,numar_traiectorii_validare,numar_pasi_validare,pas_integral,margine);
FD_val_raw = construieste_FD_derivativ(x_val_raw,u_val,xdot_val_raw,Z,Y,W,Ac,Bc,n,m,dimZ,dimY,dimW);

FD_clean = construieste_FD_derivativ(x_clean,u_data,xdot_clean,Z,Y,W,Ac,Bc,n,m,dimZ,dimY,dimW);
FD_val_clean = construieste_FD_derivativ(x_val,u_val,xdot_val,Z,Y,W,Ac,Bc,n,m,dimZ,dimY,dimW);
FD_int_clean = construieste_FD_integral(x_clean,u_data,Z,Y,W,Ac,Bc,n,m,dimZ,dimY,dimW,dt,numar_traiectorii,numar_pasi,pas_integral,margine);
FD_val_int_clean = construieste_FD_integral(x_val,u_val,Z,Y,W,Ac,Bc,n,m,dimZ,dimY,dimW,dt,numar_traiectorii_validare,numar_pasi_validare,pas_integral,margine);

fprintf('FD curat exact      : [%d x %d]\n',size(FD_clean,1),size(FD_clean,2));
fprintf('FD integral curat   : [%d x %d]\n',size(FD_int_clean,1),size(FD_int_clean,2));
fprintf('FD filtrat zgomot   : [%d x %d]\n',size(FD_filt,1),size(FD_filt,2));
fprintf('FD integral zgomot  : [%d x %d]\n',size(FD_int,1),size(FD_int,2));

sv_filt = svd(scaleaza_coloane(FD_filt),'econ');
sv_int = svd(scaleaza_coloane(FD_int),'econ');
sv_clean = svd(scaleaza_coloane(FD_clean),'econ');
sv_int_clean = svd(scaleaza_coloane(FD_int_clean),'econ');
fprintf('Cele mai mici 8 valori singulare pentru F(D) integral zgomot scalat:\n');
disp(sv_int(max(1,end-7):end).');

scala_exacta = trace_gamma0_tinta/trace(gamma0_exact);
coef_exact_scaled = scala_exacta*[T_ideal(:); N_ideal(:); M_ideal(:)];

fprintf('Reziduu solutie exacta scalata pe FD curat exact      = %.4e\n', ...
    norm(FD_clean*coef_exact_scaled)/max(1,sqrt(size(FD_clean,1))));
fprintf('Reziduu solutie exacta scalata pe FD integral curat   = %.4e\n', ...
    norm(FD_int_clean*coef_exact_scaled)/max(1,sqrt(size(FD_int_clean,1))));
fprintf('Reziduu solutie exacta scalata pe FD integral zgomot  = %.4e\n', ...
    norm(FD_int*coef_exact_scaled)/max(1,sqrt(size(FD_int,1))));

print_titlu('identificare comparativa');

rezultate = repmat(init_rezultat(),4,1);

rezultate(1) = identifica_depersis_svd( ...
    'Baza SVD',FD_clean,FD_val_clean,Z,Y,W,n,m,dimZ,dimY,dimW, ...
    regiune_test,K_feedback,prag_Jtau_minim,prag_gamma_minim, ...
    rho_jacobian0,prag_err_jacobian0,trace_gamma0_tinta,coef_exact_scaled);

rezultate(2) = identifica_depersis_svd( ...
    'Direct pe zgomot',FD_raw,FD_val_raw,Z,Y,W,n,m,dimZ,dimY,dimW, ...
    regiune_test,K_feedback,prag_Jtau_minim,prag_gamma_minim, ...
    rho_jacobian0,prag_err_jacobian0,trace_gamma0_tinta);

rezultate(3) = identifica_tikhonov( ...
    FD_filt,FD_val_filt,lambda_grid,Z,Y,W,n,m,dimZ,dimY,dimW, ...
    regiune_test,K_feedback,prag_Jtau_dorit,prag_gamma_dorit, ...
    prag_Jtau_minim,prag_gamma_minim,rho_jacobian0,prag_err_jacobian0, ...
    trace_gamma0_tinta);
rezultate(3).nume = 'Forma derivativa filtrata + reg. Tikhonov';

rezultate(4) = identifica_tikhonov( ...
    FD_int,FD_val_int,lambda_grid,Z,Y,W,n,m,dimZ,dimY,dimW, ...
    regiune_test,K_feedback,prag_Jtau_dorit,prag_gamma_dorit, ...
    prag_Jtau_minim,prag_gamma_minim,rho_jacobian0,prag_err_jacobian0, ...
    trace_gamma0_tinta);
rezultate(4).nume = 'Forma integrala + reg. Tikhonov';

for i = 1:numel(rezultate)
    if strcmp(rezultate(i).status,'OK') && rezultate(i).rez_val > prag_reziduu_validare
        rezultate(i).status = 'RESPINS_REZIDUU';
    end
end

afiseaza_rezultate_identificare(rezultate, ...
    prag_Jtau_minim,prag_gamma_minim,prag_reziduu_validare);

print_titlu('validare closed-loop');

sim_exact = simuleaza_bucla_inchisa_exact(f,g,tau_exact,delta_exact,gamma_exact,K_feedback,x0_cl,t_cl,u_actuator_max);
afiseaza_referinta_bucla_inchisa(sim_exact);

for i = 1:numel(rezultate)
    if ~strcmp(rezultate(i).status,'OK')
        fprintf('  %-28s -> sarit, metoda respinsa numeric (%s)\n', ...
            rezultate(i).nume,status_clar(rezultate(i).status));
        rezultate(i).sim = init_simulare(t_cl,n,m);
        continue;
    end

    sim_i = simuleaza_bucla_inchisa_model(f,g,rezultate(i).model,K_feedback,x0_cl,t_cl,u_actuator_max);
    rezultate(i).sim = sim_i;

    afiseaza_bucla_inchisa_metoda(rezultate(i).nume,sim_i,prag_norma_finala,u_siguranta);

    closed_ok = sim_i.ok && sim_i.norma_finala <= prag_norma_finala && sim_i.u_max <= u_siguranta;
    if ~closed_ok
        rezultate(i).status = 'RESPINS_CL';
    end
end

idx_final_curat = alege_model_final(rezultate(1),prag_norma_finala,u_siguranta);

idx_final_zgomot_local = alege_model_final(rezultate(3:4),prag_norma_finala,u_siguranta);
if isnan(idx_final_zgomot_local)
    idx_final_zgomot = NaN;
else
    idx_final_zgomot = idx_final_zgomot_local + 2;
end

idx_final = idx_final_zgomot;
if isnan(idx_final)
    model_final = [];
    sim_final = init_simulare(t_cl,n,m);
else
    model_final = rezultate(idx_final).model;
    sim_final = rezultate(idx_final).sim;
end

afiseaza_rezumat_final(rezultate,idx_final_curat,idx_final_zgomot,idx_final);

T_ref = scala_exacta*T_ideal;
N_ref = scala_exacta*N_ideal;
M_ref = scala_exacta*M_ideal;
coef_ref = [T_ref(:); N_ref(:); M_ref(:)];

raport = struct();
raport.n = n;
raport.m = m;
raport.r = r;
raport.gamma0_exact = gamma0_exact;
raport.Ac = Ac;
raport.Bc = Bc;
raport.K_feedback = K_feedback;
raport.poli_lqr = poli_lqr;
raport.rang_ctrb = rang_ctrb;
raport.grad_relativ_complet = grad_relativ_complet;
raport.dimZ = dimZ;
raport.dimY = dimY;
raport.dimW = dimW;
raport.nr_coef = nr_coef;
raport.numeZ = numeZ;
raport.numeY = numeY;
raport.numeW = numeW;
raport.dt = dt;
raport.numar_traiectorii = numar_traiectorii;
raport.numar_pasi = numar_pasi;
raport.numar_traiectorii_validare = numar_traiectorii_validare;
raport.numar_pasi_validare = numar_pasi_validare;
raport.snr_x_db = snr_x_db;
raport.sigma_x = sigma_x;
raport.rmse_x_noisy = calc_rmse(x_noisy,x_clean);
raport.rmse_x_filt = calc_rmse(x_filt(:,idx_interior),x_clean(:,idx_interior));
raport.rmse_xdot_raw = calc_rmse(xdot_raw(:,idx_interior),xdot_clean(:,idx_interior));
raport.rmse_xdot_filt = calc_rmse(xdot_filt(:,idx_interior),xdot_clean(:,idx_interior));
raport.fereastra_filtrare = fereastra_filtrare;
raport.ordin_polinom = ordin_polinom;
raport.pas_integral = pas_integral;
raport.lambda_grid = lambda_grid;
raport.rho_jacobian0 = rho_jacobian0;
raport.trace_gamma0_tinta = trace_gamma0_tinta;
raport.prag_Jtau_dorit = prag_Jtau_dorit;
raport.prag_gamma_dorit = prag_gamma_dorit;
raport.prag_Jtau_minim = prag_Jtau_minim;
raport.prag_gamma_minim = prag_gamma_minim;
raport.prag_err_jacobian0 = prag_err_jacobian0;
raport.prag_reziduu_validare = prag_reziduu_validare;
raport.prag_norma_finala = prag_norma_finala;
raport.u_siguranta = u_siguranta;
raport.x_clean = x_clean;
raport.x_noisy = x_noisy;
raport.x_filt = x_filt;
raport.xdot_clean = xdot_clean;
raport.xdot_raw = xdot_raw;
raport.xdot_filt = xdot_filt;
raport.u_data = u_data;
raport.idx_interior = idx_interior;
raport.FD_raw = FD_raw;
raport.FD_filt = FD_filt;
raport.FD_int = FD_int;
raport.FD_clean = FD_clean;
raport.FD_val_clean = FD_val_clean;
raport.FD_int_clean = FD_int_clean;
raport.FD_val_int_clean = FD_val_int_clean;
raport.FD_val_raw = FD_val_raw;
raport.FD_val_filt = FD_val_filt;
raport.FD_val_int = FD_val_int;
raport.sv_filt = sv_filt;
raport.sv_int = sv_int;
raport.sv_clean = sv_clean;
raport.sv_int_clean = sv_int_clean;
raport.rez_ref_clean = norm(FD_clean*coef_ref)/max(1,sqrt(size(FD_clean,1)));
raport.rez_ref_curat_der = raport.rez_ref_clean;
raport.rez_ref_int_clean = norm(FD_int_clean*coef_ref)/max(1,sqrt(size(FD_int_clean,1)));
raport.rez_ref_int = norm(FD_int*coef_ref)/max(1,sqrt(size(FD_int,1)));
raport.rez_ref_val = norm(FD_val_int*coef_ref)/max(1,sqrt(size(FD_val_int,1)));

raport.rez_ref_curat_int = raport.rez_ref_int_clean;
raport.rez_ref_curat_val_int = norm(FD_val_int_clean*coef_ref)/max(1,sqrt(size(FD_val_int_clean,1)));
raport.rez_ref_zgomot_int = raport.rez_ref_int;
raport.rez_ref_zgomot_val_int = raport.rez_ref_val;
raport.T_ideal = T_ideal;
raport.N_ideal = N_ideal;
raport.M_ideal = M_ideal;
raport.T_ref = T_ref;
raport.N_ref = N_ref;
raport.M_ref = M_ref;
raport.coef_ref = coef_ref;
raport.scala_exacta = scala_exacta;
raport.rezultate = rezultate;
raport.idx_final_curat = idx_final_curat;
raport.idx_final_zgomot = idx_final_zgomot;
raport.idx_final = idx_final;
raport.model_final = model_final;
if isnan(idx_final_curat)
    raport.metoda_finala_curat = [rezultate(1).nume ' (respins)'];
else
    raport.metoda_finala_curat = rezultate(idx_final_curat).nume;
end
if isnan(idx_final_zgomot)
    raport.metoda_finala_zgomot = 'nicio metoda acceptata';
else
    raport.metoda_finala_zgomot = rezultate(idx_final_zgomot).nume;
end
if isnan(idx_final)
    raport.metoda_finala = 'nicio metoda acceptata';
else
    raport.metoda_finala = rezultate(idx_final).nume;
end
raport.sim_exact = sim_exact;
raport.sim_final = sim_final;

try
    deschide_dashboard_M2(raport);
catch ME
    warning('Dashboard-ul nu a putut fi deschis: %s',ME.message);
    afiseaza_stack_dashboard_local(ME);
end

function fx = camp_vectorial_mimo(x,Jtau_exact,tau_exact,delta_exact)
    tx = tau_exact(x);
    dx = delta_exact(x);
    fx = Jtau_exact(x) \ [tx(2); dx(1); tx(4); dx(2)];
end

function gx = camp_intrare_mimo(x,Jtau_exact,gamma_exact)
    G = gamma_exact(x);
    gx = Jtau_exact(x) \ [0,      0;
                          G(1,1), G(1,2);
                          0,      0;
                          G(2,1), G(2,2)];
end

function [x_data,u_data,xdot_data] = colecteaza_date_simulare(dt,ntraj,npasi,f,g,n,m,amp)
    x_data = zeros(n,ntraj*npasi);
    u_data = zeros(m,ntraj*npasi);
    xdot_data = zeros(n,ntraj*npasi);

    for it = 1:ntraj
        x = 1.4*(rand(n,1)-0.5);
        for k = 1:npasi
            idx = (it-1)*npasi + k;

            tip = mod(it,4);
            if tip == 0
                u = amp*(2*rand(m,1)-1);
            elseif tip == 1
                u = [amp*(2*rand-1); 0];
            elseif tip == 2
                u = [0; amp*(2*rand-1)];
            else
                u = amp*[0.65*sin(0.4*k*dt+0.3*it);
                         0.55*cos(0.6*k*dt+0.2*it)];
            end

            xdot = f(x) + g(x)*u;

            x_data(:,idx) = x;
            u_data(:,idx) = u;
            xdot_data(:,idx) = xdot;

            [~,xs] = ode45(@(~,xx) f(xx) + g(xx)*u,[0 dt],x,odeset('RelTol',1e-7,'AbsTol',1e-9));
            x = xs(end,:)';
        end
    end
end

function [x_est,xdot_est] = estimeaza_brut(x_noisy,dt,ntraj,npasi)
    x_est = x_noisy;
    xdot_est = zeros(size(x_noisy));

    for it = 1:ntraj
        idx = (it-1)*npasi + (1:npasi);
        for s = 1:size(x_noisy,1)
            xdot_est(s,idx) = gradient(x_noisy(s,idx),dt);
        end
    end
end

function [x_est,xdot_est] = estimeaza_polinom_local(x_noisy,dt,ntraj,npasi,fereastra,ordin)
    if mod(fereastra,2)==0
        fereastra = fereastra + 1;
    end

    semi = floor(fereastra/2);
    ordin = max(2,min(ordin,fereastra-2));

    x_est = zeros(size(x_noisy));
    xdot_est = zeros(size(x_noisy));
    n = size(x_noisy,1);

    for it = 1:ntraj
        idx_tr = (it-1)*npasi + (1:npasi);
        for s = 1:n
            y = x_noisy(s,idx_tr);
            for k = 1:npasi
                k1 = max(1,k-semi);
                k2 = min(npasi,k+semi);
                idx_loc = k1:k2;
                tloc = ((idx_loc-k)')*dt;
                yloc = y(idx_loc)';
                ordin_loc = min(ordin,numel(idx_loc)-1);

                p = polyfit(tloc,yloc,ordin_loc);
                x_est(s,idx_tr(k)) = polyval(p,0);
                dp = polyder(p);
                xdot_est(s,idx_tr(k)) = polyval(dp,0);
            end
        end
    end
end

function idx = construieste_index_interior(ntraj,npasi,margine)
    idx = [];
    for it = 1:ntraj
        baza = (it-1)*npasi;
        idx = [idx, baza + (margine+1:npasi-margine)];
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
    nr_max = ntraj*max(npasi-2*margine-pas_int,1);
    nr_coef = n*dimZ + m*dimY + m*dimW;
    FD = zeros(n*nr_max,nr_coef);

    lin = 0;

    for it = 1:ntraj
        idx_tr = (it-1)*npasi + (1:npasi);
        for k = margine+1:npasi-margine-pas_int
            ind = idx_tr(k:k+pas_int);
            tloc = (0:pas_int)*dt;

            Zval = zeros(dimZ,numel(ind));
            Yval = zeros(dimY,numel(ind));
            Wuval = zeros(dimW,numel(ind));

            for q = 1:numel(ind)
                xq = x(:,ind(q));
                uq = u(:,ind(q));

                Zval(:,q) = Z(xq);
                Yval(:,q) = Y(xq);
                Wuval(:,q) = W(xq)*uq;
            end

            deltaZ = Zval(:,end) - Zval(:,1);
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

function rezultat = identifica_depersis_svd(nume,FD,FD_val,Z,Y,W,n,m,dimZ,dimY,dimW,regiune,Kfb,pragJ_minim,pragG_minim,rhoJ,pragErrJ0,trace_target,v_ref_svd)

    [FDs,scale] = scaleaza_coloane(FD);
    c_trace_v = constrangere_trace_gamma0(n,m,dimZ,dimY,dimW,W);
    c_trace_z = c_trace_v ./ scale(:);

    Ctau0_v = constrangere_tau0(n,m,dimZ,dimY,dimW,Z);
    Ctau0_z = Ctau0_v ./ scale(:)';

    [Cj_v,j_target] = constrangere_jacobian_tau0(n,m,dimZ,dimY,dimW,Z);
    Cj_z = Cj_v ./ scale(:)';

    rezultat = init_rezultat();
    rezultat.nume = nume;
    rezultat.lambda = 0;

    if isempty(FDs)
        rezultat.status = 'RESPINS_SVD';
        return;
    end

    [~,S,V] = svd(FDs,'econ');
    s = diag(S);
    if isempty(s)
        rezultat.status = 'RESPINS_SVD';
        return;
    end

    srel = s/max(max(s),1);
    idx_null = find(srel <= 1e-8);
    if isempty(idx_null)
        rezultat.status = 'RESPINS_SPATIU_NUL';
        rezultat.rez_train = srel(end);
        return;
    end

    Bn = V(:,idx_null);
    Ceq = [c_trace_z'; Ctau0_z];
    beq = [trace_target; zeros(n,1)];

    if nargin >= 19 && ~isempty(v_ref_svd)
        v_test = v_ref_svd(:);
        rez_test = norm(FD*v_test)/max(1,sqrt(size(FD,1))*norm(v_test));
        err_eq_test = norm(Ceq*(scale(:).*v_test) - beq);

        if isfinite(rez_test) && rez_test <= 1e-8 && err_eq_test <= 1e-7
            v = v_test;
            rezultat.observatie = 'solutie curata aleasa pe ramura de referinta din spatiul nul SVD';
            rezultat = completeaza_rezultat(rezultat,v,FD,FD_val,Z,Y,W,n,m,dimZ,dimY,dimW,regiune,Kfb,Cj_v,j_target,c_trace_v,Ctau0_v);

            if rezultat.minJ >= pragJ_minim && rezultat.minG >= pragG_minim
                rezultat.status = 'OK';
            else
                rezultat.status = 'RESPINS_ADMIS';
            end
            return;
        end
    end

    Aeq = Ceq*Bn;
    a0 = pinv_svd(Aeq)*beq;
    if norm(Aeq*a0-beq) > 1e-7
        rezultat.status = 'RESPINS_NORMALIZARE';
        return;
    end

    Nliber = null(Aeq,'r');

    if isempty(Nliber)
        a = a0;
    else
        Ared = Cj_z*Bn*Nliber;
        bred = j_target - Cj_z*Bn*a0;
        y = pinv_svd(Ared)*bred;
        a = a0 + Nliber*y;
    end

    z = Bn*a;
    v = z ./ scale(:);
    if any(~isfinite(v))
        rezultat.status = 'RESPINS_SVD';
        return;
    end

    rezultat = completeaza_rezultat(rezultat,v,FD,FD_val,Z,Y,W,n,m,dimZ,dimY,dimW,regiune,Kfb,Cj_v,j_target,c_trace_v,Ctau0_v);

    if rezultat.minJ >= pragJ_minim && rezultat.minG >= pragG_minim && rezultat.err_jacobian0 <= pragErrJ0
        rezultat.status = 'OK';
    else
        rezultat.status = 'RESPINS_ADMIS';
    end
end

function rezultat = completeaza_rezultat(rezultat,v,FD,FD_val,Z,Y,W,n,m,dimZ,dimY,dimW,regiune,Kfb,Cj_v,j_target,c_trace_v,Ctau0_v)
    model = creeaza_model_identificat(v,n,m,dimZ,dimY,dimW,Z,Y,W);

    rezultat.T = model.T;
    rezultat.N = model.N;
    rezultat.M = model.M;
    rezultat.coef = v;
    rezultat.model = model;
    rezultat.rez_train = norm(FD*v)/max(1,sqrt(size(FD,1))*norm(v));
    rezultat.rez_val = norm(FD_val*v)/max(1,sqrt(size(FD_val,1))*norm(v));

    [minJ,minG,u_rms,u_max] = evalueaza_admisibilitate_regiune(model,Kfb,regiune,n,m);
    rezultat.minJ = minJ;
    rezultat.minG = minG;
    rezultat.u_rms_reg = u_rms;
    rezultat.u_max_reg = u_max;
    rezultat.coef_norm = norm(v);
    rezultat.err_jacobian0 = norm(Cj_v*v - j_target);
    rezultat.tau0_norm = norm(Ctau0_v*v);
    rezultat.trace_gamma0 = c_trace_v'*v;
end

function rezultat = identifica_tikhonov(FD,FD_val,lambda_grid,Z,Y,W,n,m,dimZ,dimY,dimW,regiune,Kfb,pragJ_dorit,pragG_dorit,pragJ_minim,pragG_minim,rhoJ,pragErrJ0,trace_target)
    [FDs,scale] = scaleaza_coloane(FD);
    FDvals = FD_val ./ scale;

    c_trace_v = constrangere_trace_gamma0(n,m,dimZ,dimY,dimW,W);
    c_trace_z = c_trace_v ./ scale(:);

    Ctau0_v = constrangere_tau0(n,m,dimZ,dimY,dimW,Z);
    Ctau0_z = Ctau0_v ./ scale(:)';

    [Cj_v,j_target] = constrangere_jacobian_tau0(n,m,dimZ,dimY,dimW,Z);
    Cj_z = Cj_v ./ scale(:)';

    Ceq_z = [c_trace_z'; Ctau0_z];
    beq_z = [trace_target; zeros(n,1)];

    p = size(FDs,2);
    rezultat = init_rezultat();

    for il = 1:numel(lambda_grid)
        lambda = lambda_grid(il);

        H = FDs'*FDs + lambda*eye(p) + rhoJ*(Cj_z'*Cj_z);
        rhs0 = rhoJ*(Cj_z'*j_target);

        KKT = [H Ceq_z'; Ceq_z zeros(size(Ceq_z,1))];
        rhs = [rhs0; beq_z];

        sol = pinv_svd(KKT)*rhs;
        z = sol(1:p);
        v = z ./ scale(:);

        model = creeaza_model_identificat(v,n,m,dimZ,dimY,dimW,Z,Y,W);

        rez_train = norm(FD*v)/max(1,sqrt(size(FD,1)));
        rez_val = norm(FD_val*v)/max(1,sqrt(size(FD_val,1)));

        [minJ,minG,u_rms,u_max] = evalueaza_admisibilitate_regiune(model,Kfb,regiune,n,m);
        coef_norm = norm(v);
        errJ0 = norm(Cj_v*v - j_target);
        tau0_norm = norm(Ctau0_v*v);
        trG0 = c_trace_v'*v;

        penal = 0;
        penal = penal + 100*max(0,pragJ_dorit-minJ)^2;
        penal = penal + 100*max(0,pragG_dorit-minG)^2;
        penal = penal + 0.20*errJ0^2;
        penal = penal + 0.02*log10(1+coef_norm);
        penal = penal + 0.01*u_rms;

        scor = rez_val + 0.20*rez_train + penal;

        if scor < rezultat.scor
            rezultat.T = model.T;
            rezultat.N = model.N;
            rezultat.M = model.M;
            rezultat.coef = v;
            rezultat.model = model;
            rezultat.lambda = lambda;
            rezultat.rez_train = rez_train;
            rezultat.rez_val = rez_val;
            rezultat.minJ = minJ;
            rezultat.minG = minG;
            rezultat.u_rms_reg = u_rms;
            rezultat.u_max_reg = u_max;
            rezultat.coef_norm = coef_norm;
            rezultat.err_jacobian0 = errJ0;
            rezultat.tau0_norm = tau0_norm;
            rezultat.trace_gamma0 = trG0;
            rezultat.scor = scor;
        end
    end

    if rezultat.minJ >= pragJ_minim && rezultat.minG >= pragG_minim
        rezultat.status = 'OK';
    else
        rezultat.status = 'RESPINS_ADMIS';
    end
end

function [As,scale] = scaleaza_coloane(A)
    scale = vecnorm(A,2,1);
    scale(scale < 1e-12) = 1;
    As = A ./ scale;
end

function c = constrangere_trace_gamma0(n,m,dimZ,dimY,dimW,W)
    nr_coef = n*dimZ + m*dimY + m*dimW;
    c = zeros(nr_coef,1);

    W0 = W(zeros(n,1));
    offsetM = n*dimZ + m*dimY;

    for a = 1:m
        for q = 1:dimW
            idxM = offsetM + a + (q-1)*m;
            c(idxM) = c(idxM) + W0(q,a);
        end
    end
end

function Ctau0 = constrangere_tau0(n,m,dimZ,dimY,dimW,Z)
    nr_coef = n*dimZ + m*dimY + m*dimW;
    x0 = zeros(n,1);
    Z0 = Z(x0);

    Ctau0 = zeros(n,nr_coef);
    Ctau0(:,1:n*dimZ) = kron(Z0',eye(n));
end

function [Cj,j_target] = constrangere_jacobian_tau0(n,m,dimZ,dimY,dimW,Z)
    nr_coef = n*dimZ + m*dimY + m*dimW;
    Cj = zeros(n*n,nr_coef);

    x0 = zeros(n,1);
    dZ0 = jacobian_numeric(Z,x0,dimZ,n);

    Cj(:,1:n*dimZ) = kron(dZ0',eye(n));
    j_target = reshape(eye(n),[],1);
end

function model = creeaza_model_identificat(v,n,m,dimZ,dimY,dimW,Z,Y,W)
    iT = n*dimZ;
    iN = iT + m*dimY;

    T = reshape(v(1:iT),[n,dimZ]);
    N = reshape(v(iT+1:iN),[m,dimY]);
    M = reshape(v(iN+1:iN+m*dimW),[m,dimW]);

    model.T = T;
    model.N = N;
    model.M = M;
    model.Z = Z;
    model.Y = Y;
    model.W = W;
    model.tau = @(x) T*Z(x);
    model.delta = @(x) N*Y(x);
    model.gamma = @(x) M*W(x);
end

function [minJ,minG,u_rms,u_max] = evalueaza_admisibilitate_regiune(model,Kfb,regiune,n,m)
    nr = 5;
    g1 = linspace(regiune(1,1),regiune(1,2),nr);
    g2 = linspace(regiune(2,1),regiune(2,2),nr);
    g3 = linspace(regiune(3,1),regiune(3,2),nr);
    g4 = linspace(regiune(4,1),regiune(4,2),nr);

    minJ = inf;
    minG = inf;
    us = [];

    for x1 = g1
        for x2 = g2
            for x3 = g3
                for x4 = g4
                    x = [x1;x2;x3;x4];

                    J = jacobian_numeric(model.tau,x,n,n);
                    G = model.gamma(x);

                    sJ = svd(J);
                    sG = svd(G);

                    minJ = min(minJ,sJ(end));
                    minG = min(minG,sG(end));

                    if rcond(G) > 1e-8
                        u = G\(Kfb*model.tau(x) - model.delta(x));
                    else
                        u = pinv(G)*(Kfb*model.tau(x) - model.delta(x));
                    end

                    if all(isfinite(u))
                        us(:,end+1) = u;
                    end
                end
            end
        end
    end

    if isempty(us)
        u_rms = inf;
        u_max = inf;
    else
        nu = vecnorm(us);
        u_rms = sqrt(mean(nu.^2));
        u_max = max(nu);
    end
end

function sim = simuleaza_bucla_inchisa_model(f,g,model,Kfb,x0,tspan,u_lim)
    sim = init_simulare(tspan,numel(x0),size(g(x0),2));

    try
        lege = @(x) lege_control_model(x,model,Kfb,u_lim);
        [tout,xout] = ode45(@(~,x) f(x)+g(x)*lege(x),tspan,x0, ...
            odeset('RelTol',1e-6,'AbsTol',1e-8,'MaxStep',0.02,'Events',@eveniment_oprire));

        X = xout';
        m = size(g(x0),2);
        U = zeros(m,numel(tout));

        for k = 1:numel(tout)
            U(:,k) = lege(X(:,k));
        end

        sim.ok = all(isfinite(X(:))) && max(vecnorm(X)) < 40;
        sim.t = tout';
        sim.x = X;
        sim.u = U;
        sim.norma_finala = norm(X(:,end));
        sim.u_max = max(vecnorm(U));
        sim.u_rms = sqrt(mean(vecnorm(U).^2));
    catch
        sim.ok = false;
    end
end

function sim = simuleaza_bucla_inchisa_exact(f,g,tau,delta,gamma,Kfb,x0,tspan,u_lim)
    lege = @(x) lege_control_exact(x,tau,delta,gamma,Kfb,u_lim);
    [tout,xout] = ode45(@(~,x) f(x)+g(x)*lege(x),tspan,x0, ...
        odeset('RelTol',1e-6,'AbsTol',1e-8,'MaxStep',0.02,'Events',@eveniment_oprire));

    X = xout';
    m = size(g(x0),2);
    U = zeros(m,numel(tout));

    for k = 1:numel(tout)
        U(:,k) = lege(X(:,k));
    end

    sim.ok = true;
    sim.t = tout';
    sim.x = X;
    sim.u = U;
    sim.norma_finala = norm(X(:,end));
    sim.u_max = max(vecnorm(U));
    sim.u_rms = sqrt(mean(vecnorm(U).^2));
end

function u = lege_control_model(x,model,Kfb,u_lim)
    G = model.gamma(x);
    rhs = Kfb*model.tau(x) - model.delta(x);

    if rcond(G) > 1e-8
        u = G\rhs;
    else
        u = pinv(G)*rhs;
    end

    nu = norm(u);
    if nu > u_lim
        u = u_lim*u/nu;
    end
end

function u = lege_control_exact(x,tau,delta,gamma,Kfb,u_lim)
    G = gamma(x);
    rhs = Kfb*tau(x) - delta(x);

    if rcond(G) > 1e-8
        u = G\rhs;
    else
        u = pinv(G)*rhs;
    end

    nu = norm(u);
    if nu > u_lim
        u = u_lim*u/nu;
    end
end

function idx = alege_model_final(rezultate,prag_x,prag_u)
    scor = inf(numel(rezultate),1);

    for i = 1:numel(rezultate)
        if ~strcmp(rezultate(i).status,'OK')
            continue;
        end

        sim = rezultate(i).sim;
        if isempty(sim) || ~sim.ok
            continue;
        end

        penal = 0;
        penal = penal + 20*max(0,(sim.norma_finala-prag_x)/max(prag_x,1e-12))^2;
        penal = penal + 5*max(0,(sim.u_max-prag_u)/max(prag_u,1e-12))^2;

        scor(i) = rezultate(i).rez_val + 0.10*sim.norma_finala + 0.02*sim.u_rms + penal;
    end

    if all(~isfinite(scor))
        idx = NaN;
        return;
    end

    [~,idx] = min(scor);
end

function sim = init_simulare(tspan,n,m)
    sim = struct();
    sim.ok = false;
    sim.t = tspan;
    sim.x = NaN(n,numel(tspan));
    sim.u = NaN(m,numel(tspan));
    sim.norma_finala = inf;
    sim.u_max = inf;
    sim.u_rms = inf;
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
    r.minJ = 0;
    r.minG = 0;
    r.u_rms_reg = inf;
    r.u_max_reg = inf;
    r.coef_norm = inf;
    r.err_jacobian0 = inf;
    r.tau0_norm = inf;
    r.trace_gamma0 = NaN;
    r.scor = inf;
    r.status = 'RESPINS';
    r.observatie = '';
    r.sim = [];
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

function e = calc_rmse(A,B)
    d = A(:)-B(:);
    e = sqrt(mean(d.^2));
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

function afiseaza_rezultate_identificare(rezultate,pragJ,pragG,pragRez)
    fprintf('\n--------------------------------------------------------------------\n');
    fprintf('REZULTATE IDENTIFICARE MIMO M2\n');
    fprintf('--------------------------------------------------------------------\n');
    fprintf('Praguri folosite pentru status OK:\n');
    fprintf('  min sigma(J_tau)     >= %.4e\n',pragJ);
    fprintf('  min sigma(gamma)     >= %.4e\n',pragG);
    fprintf('  reziduu validare     <= %.4e\n',pragRez);
    fprintf('\n');

    fprintf('%-16s | %-22s | %-10s | %-10s | %-10s | %-10s | %-18s\n', ...
        'Caz','Metoda','Rez ID','Rez VAL','minJ','minG','Status');
    fprintf('%s\n',repmat('-',1,113));

    for i = 1:numel(rezultate)
        fprintf('%-16s | %-22s | %-10.3e | %-10.3e | %-10.3e | %-10.3e | %-18s\n', ...
            caz_metoda(i), ...
            char(string(rezultate(i).nume)), ...
            rezultate(i).rez_train, ...
            rezultate(i).rez_val, ...
            rezultate(i).minJ, ...
            rezultate(i).minG, ...
            status_clar(rezultate(i).status));
    end

    fprintf('%s\n',repmat('-',1,113));
    fprintf('\nDetalii pe metode:\n');

    for i = 1:numel(rezultate)
        fprintf('\n[%s] %s\n',caz_metoda(i),char(string(rezultate(i).nume)));
        fprintf('  lambda ales             = %.4e\n',rezultate(i).lambda);
        fprintf('  reziduu identificare    = %.4e\n',rezultate(i).rez_train);
        fprintf('  reziduu validare        = %.4e\n',rezultate(i).rez_val);
        fprintf('  min sigma(J_tau)        = %.4e\n',rezultate(i).minJ);
        fprintf('  min sigma(gamma)        = %.4e\n',rezultate(i).minG);
        fprintf('  RMS ||u|| pe regiune    = %.4e\n',rezultate(i).u_rms_reg);
        fprintf('  norma coeficienti       = %.4e\n',rezultate(i).coef_norm);
        fprintf('  err J_tau(0)            = %.4e\n',rezultate(i).err_jacobian0);
        fprintf('  ||tau(0)||              = %.4e\n',rezultate(i).tau0_norm);
        fprintf('  trace gamma(0)          = %.4e\n',rezultate(i).trace_gamma0);
        fprintf('  status                  = %s\n',status_clar(rezultate(i).status));
    end
end

function afiseaza_referinta_bucla_inchisa(sim_exact)
    fprintf('\n--------------------------------------------------------------------\n');
    fprintf('REFERINTA EXACTA IN BUCLA INCHISA\n');
    fprintf('--------------------------------------------------------------------\n');
    fprintf('Modelul exact este folosit doar pentru comparatie.\n');
    fprintf('  ||x(T)||     = %.4e\n',sim_exact.norma_finala);
    fprintf('  max ||u||    = %.4e\n',sim_exact.u_max);
    fprintf('  RMS ||u||    = %.4e\n',sim_exact.u_rms);
    fprintf('\n');
end

function afiseaza_bucla_inchisa_metoda(nume,sim,prag_x,prag_u)
    ok_x = sim.norma_finala <= prag_x;
    ok_u = sim.u_max <= prag_u;
    ok_total = sim.ok && ok_x && ok_u;

    fprintf('  %-28s -> ||x(T)|| = %.4e, max||u|| = %.4e, RMS||u|| = %.4e, closed-loop = %s\n', ...
        char(string(nume)), ...
        sim.norma_finala, ...
        sim.u_max, ...
        sim.u_rms, ...
        text_ok(ok_total));
end

function afiseaza_rezumat_final(rezultate,idx_curat,idx_zgomot,idx_practic)
    fprintf('\n--------------------------------------------------------------------\n');
    fprintf('REZUMAT FINAL MIMO M2\n');
    fprintf('--------------------------------------------------------------------\n');

    fprintf('Caz fara zgomot:\n');
    if isnan(idx_curat)
        fprintf('  metoda analizata      = %s\n',rezultate(1).nume);
        fprintf('  status                = %s\n',status_clar(rezultate(1).status));
        fprintf('  reziduu validare      = %.4e\n',rezultate(1).rez_val);
        fprintf('  min sigma(J_tau)      = %.4e\n',rezultate(1).minJ);
        fprintf('  min sigma(gamma)      = %.4e\n',rezultate(1).minG);
    else
        afiseaza_rezumat_metoda_finala(rezultate(idx_curat));
    end

    fprintf('\nCaz cu zgomot:\n');
    if isnan(idx_zgomot)
        fprintf('  nicio metoda acceptata\n');
    else
        afiseaza_rezumat_metoda_finala(rezultate(idx_zgomot));
    end

    fprintf('\nMetoda folosita ca varianta practica:\n');
    if isnan(idx_practic)
        fprintf('  nicio metoda acceptata\n');
    else
        afiseaza_rezumat_metoda_finala(rezultate(idx_practic));
    end

end

function afiseaza_rezumat_metoda_finala(r)
    fprintf('  metoda               = %s\n',char(string(r.nume)));
    fprintf('  status               = %s\n',status_clar(r.status));
    fprintf('  reziduu validare     = %.4e\n',r.rez_val);
    fprintf('  min sigma(J_tau)     = %.4e\n',r.minJ);
    fprintf('  min sigma(gamma)     = %.4e\n',r.minG);

    if isfield(r,'sim') && ~isempty(r.sim)
        fprintf('  ||x(T)|| closed-loop = %.4e\n',r.sim.norma_finala);
        fprintf('  max ||u||            = %.4e\n',r.sim.u_max);
    end
end

function caz = caz_metoda(idx)
    if idx == 1
        caz = 'fara zgomot';
    else
        caz = 'cu zgomot';
    end
end

function s = status_clar(status)
    status = string(status);

    switch status
        case "OK"
            s = 'OK';
        case "RESPINS_REZIDUU"
            s = 'respins: reziduu';
        case "RESPINS_ADMIS"
            s = 'respins: admisibilitate';
        case "RESPINS_CL"
            s = 'respins: closed-loop';
        case "RESPINS_NORMALIZARE"
            s = 'respins: normalizare';
        case "RESPINS_SPATIU_NUL"
            s = 'respins: spatiu nul';
        case "RESPINS_SVD"
            s = 'respins: SVD';
        case "RESPINS_SPATIU_NUL"
            s = 'respins: spatiu nul';
        case "RESPINS"
            s = 'respins';
        otherwise
            s = char(status);
    end
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


function deschide_dashboard_M2(r)
    r = pregateste_dashboard_mimo(r);

    t_final  = tabel_comparatie_finala(r);
    t_date   = tabel_date_comparatie(r);
    t_fd     = tabel_fd_comparatie(r);
    t_metode = tabel_metode_comparatie(r);
    t_closed = tabel_bucla_inchisa(r);

    ecran = get(groot,'ScreenSize');
    latime = min(1750,max(1360,ecran(3)-60));
    inaltime = min(980,max(840,ecran(4)-80));

    fig = uifigure('Name','Sistemul M2', ...
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

    text_status_ui(top,'Sistem','MIMO M2',true);
    text_status_ui(top,'Grad relativ',sprintf('r=[%d %d], n=%d',r.r(1),r.r(2),r.n),r.grad_relativ_complet);
    text_status_ui(top,'Zgomot stari',sprintf('SNR %.0f dB',r.snr_x_db),true);
    text_status_ui(top,'Coeficienti',sprintf('%d necunoscute',r.nr_coef),true);
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

function r = pregateste_dashboard_mimo(r)
    r.r1 = r.r(1);
    r.prag_minJtau = r.prag_Jtau_minim;
    r.prag_minGamma = r.prag_gamma_minim;
    if isfield(r,'rez_ref_curat_der')
        r.rez_ref_clean = r.rez_ref_curat_der;
    elseif ~isfield(r,'rez_ref_clean')
        r.rez_ref_clean = NaN;
    end
    if isfield(r,'rez_ref_zgomot_val_int')
        r.rez_ref_val = r.rez_ref_zgomot_val_int;
    elseif ~isfield(r,'rez_ref_val')
        r.rez_ref_val = NaN;
    end
    r.timp_cl = r.sim_exact.t;

    rezultate_initiale = r.rezultate;
    rezultate = repmat(completeaza_alias_rezultat(rezultate_initiale(1),r),size(rezultate_initiale));
    for i = 2:numel(rezultate_initiale)
        rezultate(i) = completeaza_alias_rezultat(rezultate_initiale(i),r);
    end
    r.rezultate = rezultate;

    c1 = init_caz_dashboard();
    c1.nume = "Fara zgomot - de baza";
    c1.snr_x_db = Inf;
    c1.sigma_x = 0;
    c1.x_masurat = r.x_clean;
    c1.x_val_masurat = [];
    c1.x_raw = r.x_clean;
    c1.xdot_raw = r.xdot_clean;
    c1.x_filt = r.x_clean;
    c1.xdot_filt = r.xdot_clean;
    c1.idx_interior = 1:size(r.x_clean,2);
    c1.rmse = struct('x_noisy',0,'x_filt',0,'xdot_raw',0,'xdot_filt',0);
    c1.FD_raw = r.FD_clean;
    c1.FD_filt = r.FD_clean;
    c1.FD_int = r.FD_int_clean;
    c1.FD_val_raw = r.FD_val_clean;
    c1.FD_val_filt = r.FD_val_clean;
    c1.FD_val_int = r.FD_val_int_clean;
    c1.rez_ref_int = r.rez_ref_curat_int;
    c1.rez_ref_val_int = r.rez_ref_curat_val_int;
    c1.rezultate = rezultate(1);
    c1.sim = struct();
    c1.sim.exact = alias_simulare(r.sim_exact,r.n,r.m);
    c1.sim.metode = alias_simulare(rezultate(1).sim,r.n,r.m);
    if ~isnan(r.idx_final_curat)
        c1.idx_final = 1;
        c1.metoda_finala = rezultate(1).nume;
    end

    c2 = init_caz_dashboard();
    c2.nume = "Cu zgomot";
    c2.snr_x_db = r.snr_x_db;
    c2.sigma_x = r.sigma_x;
    c2.x_masurat = r.x_noisy;
    c2.x_val_masurat = [];
    c2.x_raw = r.x_noisy;
    c2.xdot_raw = r.xdot_raw;
    c2.x_filt = r.x_filt;
    c2.xdot_filt = r.xdot_filt;
    c2.idx_interior = r.idx_interior;
    c2.rmse = struct('x_noisy',r.rmse_x_noisy,'x_filt',r.rmse_x_filt,'xdot_raw',r.rmse_xdot_raw,'xdot_filt',r.rmse_xdot_filt);
    c2.FD_raw = r.FD_raw;
    c2.FD_filt = r.FD_filt;
    c2.FD_int = r.FD_int;
    c2.FD_val_raw = r.FD_val_raw;
    c2.FD_val_filt = r.FD_val_filt;
    c2.FD_val_int = r.FD_val_int;
    c2.rez_ref_int = r.rez_ref_zgomot_int;
    c2.rez_ref_val_int = r.rez_ref_zgomot_val_int;
    c2.rezultate = rezultate(2:4);
    c2.sim = struct();
    c2.sim.exact = alias_simulare(r.sim_exact,r.n,r.m);
    c2.sim.metode = repmat(alias_simulare(init_simulare(r.sim_exact.t,r.n,r.m),r.n,r.m),3,1);
    for k = 1:3
        c2.sim.metode(k) = alias_simulare(rezultate(k+1).sim,r.n,r.m);
    end
    if ~isnan(r.idx_final_zgomot)
        c2.idx_final = r.idx_final_zgomot - 1;
        c2.metoda_finala = rezultate(r.idx_final_zgomot).nume;
    end

    r.cazuri = [c1; c2];
end

function c = init_caz_dashboard()
    c = struct();
    c.nume = "";
    c.snr_x_db = NaN;
    c.sigma_x = NaN;
    c.x_masurat = [];
    c.x_val_masurat = [];
    c.x_raw = [];
    c.xdot_raw = [];
    c.x_filt = [];
    c.xdot_filt = [];
    c.idx_interior = [];
    c.rmse = struct('x_noisy',NaN,'x_filt',NaN,'xdot_raw',NaN,'xdot_filt',NaN);
    c.FD_raw = [];
    c.FD_filt = [];
    c.FD_int = [];
    c.FD_val_raw = [];
    c.FD_val_filt = [];
    c.FD_val_int = [];
    c.rez_ref_int = NaN;
    c.rez_ref_val_int = NaN;
    c.rezultate = [];
    c.sim = [];
    c.idx_final = NaN;
    c.metoda_finala = 'nicio metoda acceptata';
end

function rr = completeaza_alias_rezultat(rr,r)
    if isfield(rr,'minJ')
        rr.minJtau = rr.minJ;
    else
        rr.minJtau = NaN;
    end
    if isfield(rr,'minG')
        rr.minGamma = rr.minG;
    else
        rr.minGamma = NaN;
    end
    if isfield(rr,'err_jacobian0')
        rr.errJ_ech = rr.err_jacobian0;
    else
        rr.errJ_ech = NaN;
    end
    if isfield(rr,'trace_gamma0')
        rr.gamma_ech = rr.trace_gamma0;
    else
        rr.gamma_ech = NaN;
    end
    if isfield(rr,'model') && ~isempty(rr.model)
        rr.tau_ech_norm = norm(rr.model.tau(zeros(r.n,1)));
        rr.err_tau = eroare_matrice_procente(rr.T,r.T_ref);
        rr.err_delta = eroare_matrice_procente(rr.N,r.N_ref);
        rr.err_gamma = eroare_matrice_procente(rr.M,r.M_ref);
        rr.err_u = rr.u_rms_reg;
    else
        rr.tau_ech_norm = NaN;
        rr.err_tau = NaN;
        rr.err_delta = NaN;
        rr.err_gamma = NaN;
        rr.err_u = NaN;
    end
end

function e = eroare_matrice_procente(A,B)
    if isempty(A) || isempty(B)
        e = NaN;
    else
        e = 100*norm(A-B,'fro')/max(norm(B,'fro'),1e-12);
    end
end

function s = alias_simulare(sim,n,m)
    if isempty(sim)
        sim = init_simulare(0,n,m);
    end
    s = struct();
    if isfield(sim,'ok')
        s.simulat = sim.ok;
        s.ok = sim.ok;
    else
        s.simulat = false;
        s.ok = false;
    end
    if isfield(sim,'t')
        s.T = sim.t;
    else
        s.T = [];
    end
    if isfield(sim,'x')
        s.X = sim.x;
    else
        s.X = NaN(n,0);
    end
    if isfield(sim,'u')
        s.U = sim.u;
    else
        s.U = NaN(m,0);
    end
    if isfield(sim,'norma_finala')
        s.err_finala = sim.norma_finala;
    else
        s.err_finala = NaN;
    end
    if isfield(sim,'u_max')
        s.u_max = sim.u_max;
    else
        s.u_max = NaN;
    end
    if isfield(sim,'u_rms')
        s.u_rms = sim.u_rms;
    else
        s.u_rms = NaN;
    end
    if ~isempty(s.X)
        s.err_traj = vecnorm(s.X);
    else
        s.err_traj = [];
    end
end

function construieste_tab_rezumat(tab,r,t_final,~)
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
          "Semnalul de intrare u este considerat cunoscut, deoarece reprezinta comanda aplicata procesului si salvata in setul de date. " + ...
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

    for s = 1:r.n
        grafic_stari_bucla_inchisa(uiaxes(g2),r,s);
    end
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

    txt = "Datele de identificare sunt obtinute prin simularea sistemului neliniar M2 in jurul punctului de echilibru." + newline + ...
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

    tb_u = uitab(group,'Title','Comenzi u aplicate');
    gl_u = uigridlayout(tb_u,[2 1]);
    gl_u.Padding = [10 8 10 10];
    gl_u.RowSpacing = 10;
    grafic_intrare(uiaxes(gl_u),r,1);
    grafic_intrare(uiaxes(gl_u),r,2);
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

    grafic_valori_singulare(uiaxes(g2),r.cazuri(1).FD_int,'Valori singulare ale matricei F(D) integrale - fara zgomot');
    grafic_valori_singulare(uiaxes(g2),r.cazuri(2).FD_int,'Valori singulare ale matricei F(D) integrale - cu zgomot');

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
    grafic_prag_metode(uiaxes(gl2),r,'minJtau',r.prag_minJtau,'Conditia de difeomorfism: min sigma(J_tau)');
    grafic_prag_metode(uiaxes(gl2),r,'minGamma',r.prag_minGamma,'Conditia de decuplare: min sigma(gamma)');

    tb3 = uitab(group,'Title','Erori fata de referinta');
    gl3 = uigridlayout(tb3,[2 2]);
    gl3.Padding = [10 8 10 10];
    gl3.RowSpacing = 10;
    gl3.ColumnSpacing = 10;
    grafic_metrici_metode(uiaxes(gl3),r,'err_tau','Eroarea relativa a transformarii tau(x)','[%]');
    grafic_metrici_metode(uiaxes(gl3),r,'err_delta','Eroarea relativa a termenului delta(x)','[%]');
    grafic_metrici_metode(uiaxes(gl3),r,'err_gamma','Eroarea relativa a termenului gamma(x)','[%]');
    grafic_metrici_metode(uiaxes(gl3),r,'err_u','Efort mediu de control pe regiune','u');
end

function construieste_tab_transformare(tab,r)
    lay = uigridlayout(tab,[1 2]);
    lay.ColumnWidth = {640,'1x'};
    lay.Padding = [16 14 16 16];
    lay.ColumnSpacing = 14;
    lay.BackgroundColor = [0.950 0.955 0.965];

    p1 = creeaza_panou(lay,'5.1 Normalizari si verificari de admisibilitate');
    g1 = uigridlayout(p1,[4 1]);
    g1.RowHeight = {110,120,'1x','1x'};
    g1.Padding = [14 11 14 13];
    g1.RowSpacing = 9;

    txt = "Normalizarile impuse fixeaza originea, scara si orientarea locala a transformarii identificate." + newline + ...
          "Conditia tau(0)=0 fixeaza originea in noile coordonate, iar trace(gamma(0))=m elimina ambiguitatea de scala a solutiei." + newline + ...
          "Orientarea locala este controlata prin apropierea jacobianului J_tau(0) de matricea identitate." + newline + ...
          "Admisibilitatea este verificata prin valoarea singulara minima a lui J_tau si a matricei de decuplare gamma(x).";
    adauga_text(g1,txt,12.3);

    adauga_tabel(g1,tabel_transformare_finala(r), ...
        ["Caz","Metoda","Tau_0","Trace_gamma0","Err_J0","Min_Jtau","Min_gamma"], ...
        ["Caz","Metoda","tau(0)","trace gamma(0)","err J0","min Jtau","min gamma"], ...
        {82,125,70,92,78,78,'1x'},38,10,10.5);

    grafic_normalizari_transformare(uiaxes(g1),r);
    grafic_admisibilitate_transformare(uiaxes(g1),r);

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
        grafic_eroare_matrice(uiaxes(uigridlayout(tbT,[1 1])),r.T_ref,Tsel,sprintf('T - %s',nume_txt),r.numeZ,compose('tau%d',1:r.n));

        tbN = uitab(inner,'Title','Termen neliniar N');
        grafic_eroare_matrice(uiaxes(uigridlayout(tbN,[1 1])),r.N_ref,Nsel,sprintf('N - %s',nume_txt),r.numeY,compose('delta%d',1:r.m));

        tbM = uitab(inner,'Title','Decuplare M');
        grafic_eroare_matrice(uiaxes(uigridlayout(tbM,[1 1])),r.M_ref,Msel,sprintf('M - %s',nume_txt),r.numeW,compose('gamma%d',1:r.m));
    end
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
        ["Caz","Metoda","||x(T)||","max ||u||","RMS ||u||","Status"], ...
        {110,150,90,80,80,'1x'},35,10.4,10.9);

    grafic_metrici_bucla_inchisa(uiaxes(g1),r);
    adauga_text(g1,text_interpretare_bucla_inchisa(r),13.5);

    p2 = creeaza_panou(lay,'6.2 Raspunsul sistemului controlat');
    g2 = uigridlayout(p2,[1 1]);
    g2.Padding = [10 8 10 10];
    g2.RowSpacing = 0;
    g2.ColumnSpacing = 0;

    group = uitabgroup(g2);

    tb1 = uitab(group,'Title','Stari');
    gl1 = uigridlayout(tb1,[2 2]);
    gl1.RowHeight = {'1x','1x'};
    gl1.ColumnWidth = {'1x','1x'};
    gl1.Padding = [10 8 10 10];
    gl1.RowSpacing = 10;
    gl1.ColumnSpacing = 10;
    for s = 1:r.n
        ax = uiaxes(gl1);
        ax.Layout.Row = ceil(s/2);
        ax.Layout.Column = mod(s-1,2) + 1;
        grafic_stari_bucla_inchisa(ax,r,s);
    end

    tb2 = uitab(group,'Title','Norma erorii');
    gl2 = uigridlayout(tb2,[1 1]);
    gl2.Padding = [10 8 10 10];
    ax2 = uiaxes(gl2);
    ax2.Layout.Row = 1;
    ax2.Layout.Column = 1;
    grafic_norma_bucla_inchisa(ax2,r);

    tb3 = uitab(group,'Title','Comenzi');
    gl3 = uigridlayout(tb3,[2 1]);
    gl3.RowHeight = {'1x','1x'};
    gl3.Padding = [10 8 10 10];
    gl3.RowSpacing = 10;
    for iu = 1:r.m
        ax = uiaxes(gl3);
        ax.Layout.Row = iu;
        ax.Layout.Column = 1;
        grafic_comanda_bucla_inchisa(ax,r,iu);
    end

    tb4 = uitab(group,'Title','Plan de faza');
    gl4 = uigridlayout(tb4,[1 1]);
    gl4.Padding = [10 8 10 10];
    ax4 = uiaxes(gl4);
    ax4.Layout.Row = 1;
    ax4.Layout.Column = 1;
    grafic_plan_faza(ax4,r);

    group.SelectedTab = tb1;
    drawnow;
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
        SNR(i) = string(snr_text_local(c.snr_x_db));
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
        SNR(i) = string(snr_text_local(c.snr_x_db));
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

    c = r.cazuri(1);
    [Caz,Matrice,Date,Rol,Randuri,Coloane] = adauga_linie_fd_simplu(Caz,Matrice,Date,Rol,Randuri,Coloane,"Fara zgomot","FD baza","identificare","constructie cu x, u si xdot exacte",c.FD_filt);
    [Caz,Matrice,Date,Rol,Randuri,Coloane] = adauga_linie_fd_simplu(Caz,Matrice,Date,Rol,Randuri,Coloane,"Fara zgomot","FD baza","validare","verificare pe set separat",c.FD_val_filt);

    c = r.cazuri(2);
    [Caz,Matrice,Date,Rol,Randuri,Coloane] = adauga_linie_fd_simplu(Caz,Matrice,Date,Rol,Randuri,Coloane,"Cu zgomot","FD brut","identificare","test SVD pe date zgomotoase",c.FD_raw);
    [Caz,Matrice,Date,Rol,Randuri,Coloane] = adauga_linie_fd_simplu(Caz,Matrice,Date,Rol,Randuri,Coloane,"Cu zgomot","FD filtrat","identificare","forma derivativa filtrata + Tikhonov",c.FD_filt);
    [Caz,Matrice,Date,Rol,Randuri,Coloane] = adauga_linie_fd_simplu(Caz,Matrice,Date,Rol,Randuri,Coloane,"Cu zgomot","FD integral","identificare","forma integrala + Tikhonov",c.FD_int);
    [Caz,Matrice,Date,Rol,Randuri,Coloane] = adauga_linie_fd_simplu(Caz,Matrice,Date,Rol,Randuri,Coloane,"Cu zgomot","FD brut","validare","validare SVD directa",c.FD_val_raw);
    [Caz,Matrice,Date,Rol,Randuri,Coloane] = adauga_linie_fd_simplu(Caz,Matrice,Date,Rol,Randuri,Coloane,"Cu zgomot","FD filtrat","validare","validare forma derivativa",c.FD_val_filt);
    [Caz,Matrice,Date,Rol,Randuri,Coloane] = adauga_linie_fd_simplu(Caz,Matrice,Date,Rol,Randuri,Coloane,"Cu zgomot","FD integral","validare","validare forma integrala",c.FD_val_int);

    t = table(Caz,Matrice,Date,Rol,Randuri,Coloane);
end

function [Caz,Matrice,Date,Rol,Randuri,Coloane] = adauga_linie_fd_simplu(Caz,Matrice,Date,Rol,Randuri,Coloane,caz_txt,matrice_txt,date_txt,rol_txt,F)
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
            Rol(end+1,1) = rol_metoda_scurt(c.nume,rr.nume);
            Lambda(end+1,1) = format_lambda_tabel(rr.lambda);
            Rez_ID(end+1,1) = format_numar_sau_linie(rr.rez_train);
            Rez_VAL(end+1,1) = format_numar_sau_linie(rr.rez_val);
            Min_Jtau(end+1,1) = format_numar_sau_linie(rr.minJtau);
            Min_gamma(end+1,1) = format_numar_sau_linie(rr.minGamma);
            Err_Jxe(end+1,1) = format_numar_sau_linie(rr.errJ_ech);

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

function t = tabel_transformare_finala(r)
    Caz = strings(0,1);
    Metoda = strings(0,1);
    Tau_0 = strings(0,1);
    Trace_gamma0 = strings(0,1);
    Err_J0 = strings(0,1);
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
        Tau_0(end+1,1) = format_numar_tabel(rr.tau_ech_norm);
        Trace_gamma0(end+1,1) = format_numar_tabel(rr.gamma_ech);
        Err_J0(end+1,1) = format_numar_tabel(rr.errJ_ech);
        Min_Jtau(end+1,1) = format_numar_tabel(rr.minJtau);
        Min_gamma(end+1,1) = format_numar_tabel(rr.minGamma);
    end

    t = table(Caz,Metoda,Tau_0,Trace_gamma0,Err_J0,Min_Jtau,Min_gamma);
end

function t = tabel_bucla_inchisa(r)
    Caz = strings(0,1);
    Metoda = strings(0,1);
    Err_finala = strings(0,1);
    Max_u = strings(0,1);
    RMS_u = strings(0,1);
    Status = strings(0,1);

    for ic = 1:numel(r.cazuri)
        c = r.cazuri(ic);
        for i = 1:numel(c.rezultate)
            rr = c.rezultate(i);
            if isfield(c,'sim') && isfield(c.sim,'metode') && i <= numel(c.sim.metode)
                si = c.sim.metode(i);
            else
                si = [];
            end

            Caz(end+1,1) = nume_caz_scurt(c.nume);
            Metoda(end+1,1) = nume_metoda_scurt(rr.nume);

            if ~isempty(si) && isfield(si,'T') && ~isempty(si.T) && isfield(si,'simulat') && si.simulat
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

function txt = text_identificare_sinteza(r)
    buc = strings(0,1);
    c1 = r.cazuri(1);
    c2 = r.cazuri(2);

    if ~isnan(c1.idx_final)
        rr1 = c1.rezultate(c1.idx_final);
        buc(end+1) = sprintf('Pentru cazul ideal, fara zgomot, metoda bazata pe spatiul nul SVD furnizeaza o solutie admisibila. Pe setul de validare se obtine un reziduu de %.3e, cu min sigma(J_tau)=%.3e si min sigma(gamma)=%.3e.',rr1.rez_val,rr1.minJtau,rr1.minGamma);
    else
        buc(end+1) = "Pentru cazul ideal, fara zgomot, nu a fost selectata o metoda finala admisibila.";
    end

    if ~isnan(c2.idx_final)
        rr2 = c2.rezultate(c2.idx_final);
        buc(end+1) = sprintf('Pentru cazul cu zgomot pe stari, metoda selectata este forma integrala cu regularizare Tikhonov. Aceasta respecta conditiile de admisibilitate si ramane stabila in bucla inchisa, avand reziduu de validare %.3e, min sigma(J_tau)=%.3e si min sigma(gamma)=%.3e.',rr2.rez_val,rr2.minJtau,rr2.minGamma);
    else
        buc(end+1) = "Pentru cazul cu zgomot pe stari, nu a fost selectata o metoda finala admisibila.";
    end

    txt = strjoin(buc,newline);
end

function txt = text_interpretare_bucla_inchisa(r)
    c1 = r.cazuri(1);
    c2 = r.cazuri(2);
    b = strings(0,1);
    if ~isnan(c1.idx_final)
        b(end+1) = "Pentru cazul ideal, metoda bazata pe spatiul nul SVD reproduce raspunsul de referinta si stabilizeaza sistemul in vecinatatea echilibrului.";
    end
    if ~isnan(c2.idx_final)
        b(end+1) = "Pentru cazul cu zgomot pe stari, forma integrala cu regularizare Tikhonov mentine convergenta catre echilibru si pastreaza comanda intr-un interval limitat.";
    end
    if isempty(b)
        txt = "Nu exista metode finale acceptate pentru interpretarea raspunsului in bucla inchisa.";
    else
        txt = strjoin(b,newline);
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

function [T,X,U] = extrage_date_simulare_grafic(sim,n,m)
    T = [];
    X = NaN(n,0);
    U = NaN(m,0);

    if isempty(sim) || ~isstruct(sim)
        return;
    end

    if isfield(sim,'T')
        T = sim.T;
    elseif isfield(sim,'t')
        T = sim.t;
    end

    if isfield(sim,'X')
        X = sim.X;
    elseif isfield(sim,'x')
        X = sim.x;
    end

    if isfield(sim,'U')
        U = sim.U;
    elseif isfield(sim,'u')
        U = sim.u;
    end

    if isrow(T)
        T = T(:)';
    end
end

function grafic_stari_bucla_inchisa(ax,r,s)
    cla(ax);
    hold(ax,'on');
    ax.Visible = 'on';
    exista = false;

    [T0,X0,~] = extrage_date_simulare_grafic(r.cazuri(1).sim.exact,r.n,r.m);
    if ~isempty(T0) && ~isempty(X0) && size(X0,1) >= s && size(X0,2) == numel(T0)
        plot(ax,T0,X0(s,:),'-','LineWidth',1.4,'DisplayName','referinta exacta');
        exista = true;
    end

    for ic = 1:numel(r.cazuri)
        c = r.cazuri(ic);
        if isnan(c.idx_final)
            continue;
        end
        [Ti,Xi,~] = extrage_date_simulare_grafic(c.sim.metode(c.idx_final),r.n,r.m);
        if ~isempty(Ti) && ~isempty(Xi) && size(Xi,1) >= s && size(Xi,2) == numel(Ti)
            plot(ax,Ti,Xi(s,:),'--','LineWidth',1.2,'DisplayName',char(nume_caz_scurt(c.nume)));
            exista = true;
        end
    end

    if exista
        yline(ax,0,':','xe','HandleVisibility','off');
        legend(ax,'Location','best');
    else
        text(ax,0.5,0.5,'Nu exista date de simulare pentru afisare','HorizontalAlignment','center','Units','normalized');
    end

    title(ax,sprintf('Evolutia starii x%d in bucla inchisa',s));
    xlabel(ax,'Timp [s]');
    ylabel(ax,sprintf('Stare x%d',s));
    grid(ax,'on');
    stil_axe(ax);
end

function grafic_norma_bucla_inchisa(ax,r)
    cla(ax);
    hold(ax,'on');
    ax.Visible = 'on';
    exista = false;

    [T0,X0,~] = extrage_date_simulare_grafic(r.cazuri(1).sim.exact,r.n,r.m);
    if ~isempty(T0) && ~isempty(X0) && size(X0,2) == numel(T0)
        semilogy(ax,T0,max(vecnorm(X0),1e-14),'-','LineWidth',1.5,'DisplayName','referinta exacta');
        exista = true;
    end

    for ic = 1:numel(r.cazuri)
        c = r.cazuri(ic);
        if isnan(c.idx_final)
            continue;
        end
        [Ti,Xi,~] = extrage_date_simulare_grafic(c.sim.metode(c.idx_final),r.n,r.m);
        if ~isempty(Ti) && ~isempty(Xi) && size(Xi,2) == numel(Ti)
            semilogy(ax,Ti,max(vecnorm(Xi),1e-14),'--','LineWidth',1.3,'DisplayName',char(nume_caz_scurt(c.nume)));
            exista = true;
        end
    end

    if exista
        legend(ax,'Location','best');
    else
        text(ax,0.5,0.5,'Nu exista date de simulare pentru afisare','HorizontalAlignment','center','Units','normalized');
    end

    title(ax,'Convergenta erorii de stare in bucla inchisa');
    xlabel(ax,'Timp [s]');
    ylabel(ax,'Norma erorii de stare');
    grid(ax,'on');
    stil_axe(ax);
end

function grafic_comanda_bucla_inchisa(ax,r,iu)
    cla(ax);
    hold(ax,'on');
    ax.Visible = 'on';
    exista = false;

    [T0,~,U0] = extrage_date_simulare_grafic(r.cazuri(1).sim.exact,r.n,r.m);
    if ~isempty(T0) && ~isempty(U0) && size(U0,1) >= iu && size(U0,2) == numel(T0)
        plot(ax,T0,U0(iu,:),'-','LineWidth',1.4,'DisplayName','referinta exacta');
        exista = true;
    end

    for ic = 1:numel(r.cazuri)
        c = r.cazuri(ic);
        if isnan(c.idx_final)
            continue;
        end
        [Ti,~,Ui] = extrage_date_simulare_grafic(c.sim.metode(c.idx_final),r.n,r.m);
        if ~isempty(Ti) && ~isempty(Ui) && size(Ui,1) >= iu && size(Ui,2) == numel(Ti)
            plot(ax,Ti,Ui(iu,:),'--','LineWidth',1.2,'DisplayName',char(nume_caz_scurt(c.nume)));
            exista = true;
        end
    end

    if exista
        yline(ax,0,'--','u_echilibru','HandleVisibility','off');
        legend(ax,'Location','best');
    else
        text(ax,0.5,0.5,'Nu exista date de comanda pentru afisare','HorizontalAlignment','center','Units','normalized');
    end

    title(ax,sprintf('Semnalul de comanda u%d aplicat in bucla inchisa',iu));
    xlabel(ax,'Timp [s]');
    ylabel(ax,sprintf('Comanda u%d',iu));
    grid(ax,'on');
    stil_axe(ax);
end

function grafic_plan_faza(ax,r)
    cla(ax);
    hold(ax,'on');
    ax.Visible = 'on';
    exista = false;

    [~,X0,~] = extrage_date_simulare_grafic(r.cazuri(1).sim.exact,r.n,r.m);
    if ~isempty(X0) && size(X0,1) >= 4
        plot(ax,X0(1,:),X0(2,:),'-','LineWidth',1.3,'DisplayName','referinta exacta x1-x2');
        plot(ax,X0(3,:),X0(4,:),'-','LineWidth',1.3,'DisplayName','referinta exacta x3-x4');
        exista = true;
    end

    for ic = 1:numel(r.cazuri)
        c = r.cazuri(ic);
        if isnan(c.idx_final)
            continue;
        end
        [~,Xi,~] = extrage_date_simulare_grafic(c.sim.metode(c.idx_final),r.n,r.m);
        if ~isempty(Xi) && size(Xi,1) >= 4
            plot(ax,Xi(1,:),Xi(2,:),'--','LineWidth',1.1,'DisplayName',char(nume_caz_scurt(c.nume)) + " x1-x2");
            plot(ax,Xi(3,:),Xi(4,:),':','LineWidth',1.1,'DisplayName',char(nume_caz_scurt(c.nume)) + " x3-x4");
            exista = true;
        end
    end

    if exista
        plot(ax,0,0,'rx','LineWidth',1.8,'MarkerSize',8,'HandleVisibility','off');
        legend(ax,'Location','best');
    else
        text(ax,0.5,0.5,'Nu exista date de simulare pentru afisare','HorizontalAlignment','center','Units','normalized');
    end

    title(ax,'Traiectoria sistemului in planurile de faza');
    xlabel(ax,'Stare de pozitie');
    ylabel(ax,'Stare de viteza');
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
    plot(ax,t,r.x_clean(s,idx),'k-','LineWidth',1.4,'DisplayName','curata');
    plot(ax,t,c.x_masurat(s,idx),'.','MarkerSize',6,'DisplayName','masurata');
    plot(ax,t,c.x_filt(s,idx),'--','LineWidth',1.2,'DisplayName','filtrata');
    yline(ax,0,':','0','HandleVisibility','off');
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
    plot(ax,t,r.xdot_clean(s,idx),'k-','LineWidth',1.4,'DisplayName','xdot exacta');
    plot(ax,t,c.xdot_raw(s,idx),':','LineWidth',1.1,'DisplayName','xdot brut');
    plot(ax,t,c.xdot_filt(s,idx),'--','LineWidth',1.2,'DisplayName','xdot filtrat');
    title(ax,sprintf('Derivata dx%d/dt estimata - %s',s,caz_txt));
    xlabel(ax,'Timp [s]');
    ylabel(ax,sprintf('dx%d/dt',s));
    xlim(ax,[t(1) t(end)]);
    grid(ax,'on');
    legend(ax,'Location','best');
    stil_axe(ax);
end

function grafic_intrare(ax,r,comp)
    traj = 4;
    idx = (traj-1)*r.numar_pasi + (1:r.numar_pasi);
    idx = idx(idx <= size(r.u_data,2));

    t = (0:numel(idx)-1)*r.dt;

    cla(ax);
    plot(ax,t,r.u_data(comp,idx),'LineWidth',1.2);

    title(ax,sprintf('Semnalul de intrare u%d folosit la generarea datelor',comp));
    xlabel(ax,'Timp [s]');
    ylabel(ax,sprintf('Comanda u%d',comp));
    xlim(ax,[t(1) t(end)]);
    grid(ax,'on');
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

function grafic_valori_singulare(ax,F,titlu)
    cla(ax);
    if isempty(F)
        text(ax,0.5,0.5,'Matrice indisponibila','HorizontalAlignment','center');
        axis(ax,'off');
        return;
    end
    s = svd(scaleaza_coloane(F),'econ');
    semilogy(ax,s,'o-','LineWidth',1.1,'MarkerSize',4);
    title(ax,titlu);
    xlabel(ax,'Index valoare singulara');
    ylabel(ax,'Valoare singulara');
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

function grafic_reziduuri_identificare_validare(ax,r)
    labels = strings(0,1);
    rez_id = [];
    rez_val = [];
    for ic = 1:numel(r.cazuri)
        c = r.cazuri(ic);
        for i = 1:numel(c.rezultate)
            rr = c.rezultate(i);
            if ~afiseaza_metoda_in_grafic(rr)
                continue;
            end
            labels(end+1,1) = eticheta_metoda_grafic(c.nume,rr.nume);
            rez_id(end+1,1) = max(rr.rez_train,1e-14);
            rez_val(end+1,1) = max(rr.rez_val,1e-14);
        end
    end
    cla(ax);
    if isempty(labels)
        text(ax,0.5,0.5,'Nu exista reziduuri disponibile pentru afisare','HorizontalAlignment','center');
        axis(ax,'off');
        return;
    end
    bar(ax,[rez_id rez_val],'grouped');
    set(ax,'YScale','log');
    ax.XTick = 1:numel(labels);
    ax.XTickLabel = cellstr(labels);
    ax.XTickLabelRotation = 12;
    title(ax,'Compararea reziduurilor pe identificare si validare');
    xlabel(ax,'Metoda analizata');
    ylabel(ax,'Reziduu normalizat');
    legend(ax,{'Identificare','Validare'},'Location','best');
    grid(ax,'on');
    stil_axe(ax);
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
    ok = isfinite(rr.rez_train) || isfinite(rr.rez_val);
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
function grafic_normalizari_transformare(ax,r)
    [labels,tau0,gamma0,errJ] = date_transformare_finale(r);
    cla(ax);
    if isempty(labels)
        text(ax,0.5,0.5,'Nu exista solutii finale pentru acest grafic','HorizontalAlignment','center');
        axis(ax,'off');
        return;
    end
    vals = [max(tau0,1e-14), max(abs(gamma0-r.m),1e-14), max(errJ,1e-14)];
    bar(ax,vals,'grouped');
    set(ax,'YScale','log');
    ax.XTick = 1:numel(labels);
    ax.XTickLabel = cellstr(labels);
    title(ax,'Erori ale conditiilor de normalizare');
    xlabel(ax,'Caz analizat');
    ylabel(ax,'Valoare pe scala logaritmica');
    legend(ax,{'||tau(0)||','|trace gamma(0)-m|','||J_tau(0)-I||'},'Location','best');
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
    yline(ax,r.prag_minJtau,'--k','prag','LineWidth',1.1,'LabelHorizontalAlignment','right','LabelVerticalAlignment','bottom');
    set(ax,'YScale','log');
    ax.XTick = 1:numel(labels);
    ax.XTickLabel = cellstr(labels);
    title(ax,'Indicatori de admisibilitate pe regiunea de test');
    xlabel(ax,'Caz analizat');
    ylabel(ax,'Valoare minima');
    legend(ax,{'min sigma(J_tau)','min sigma(gamma)','prag'},'Location','best');
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
        gamma0(end+1,1) = abs(rr.gamma_ech);
        errJ(end+1,1) = abs(rr.errJ_ech);
        minJ(end+1,1) = rr.minJtau;
        minG(end+1,1) = rr.minGamma;
    end
end

function [Tsel,Nsel,Msel,nume] = obtine_model_final_caz(c,r)
    if isnan(c.idx_final)
        Tsel = NaN(size(r.T_ref));
        Nsel = NaN(size(r.N_ref));
        Msel = NaN(size(r.M_ref));
        nume = "fara metoda";
        return;
    end
    rr = c.rezultate(c.idx_final);
    Tsel = rr.T;
    Nsel = rr.N;
    Msel = rr.M;
    nume = nume_metoda_scurt(rr.nume);
end

function grafic_eroare_matrice(ax,Mref,Msel,titlu,nume_col,nume_linii)
    cla(ax);
    if isempty(Msel) || any(isnan(Msel(:)))
        text(ax,0.5,0.5,'Nu exista metoda finala pentru comparatie','HorizontalAlignment','center');
        axis(ax,'off');
        return;
    end
    E = Msel - Mref;
    imagesc(ax,E);
    colorbar(ax);
    title(ax,titlu);
    ax.XTick = 1:numel(nume_col);
    ax.XTickLabel = nume_col;
    ax.YTick = 1:numel(nume_linii);
    ax.YTickLabel = nume_linii;
    xtickangle(ax,45);
    xlabel(ax,'Functii din dictionar');
    ylabel(ax,'Randuri matrice');
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
    legend(ax,{'Eroare finala','Max ||u||','RMS ||u||'},'Location','best');
    grid(ax,'on');
    stil_axe(ax);
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

function stil_axe(ax)
    try
        ax.FontName = 'Arial';
        ax.FontSize = 10.5;
        ax.Box = 'on';
        ax.GridAlpha = 0.18;
    catch
    end
end

function txt = format_numar_sau_linie(x)
    if isempty(x) || ~isscalar(x) || ~isfinite(x)
        txt = "-";
    else
        txt = format_numar_tabel(x);
    end
end

function txt = format_numar_tabel(x)
    if isempty(x) || ~isscalar(x) || ~isfinite(x)
        txt = "-";
    elseif abs(x) >= 1e3 || (abs(x) < 1e-2 && abs(x) > 0)
        txt = string(sprintf('%.2e',x));
    else
        txt = string(sprintf('%.4f',x));
    end
end

function txt = format_lambda_tabel(lambda)
    if isempty(lambda) || ~isfinite(lambda)
        txt = "-";
    elseif lambda == 0
        txt = "SVD";
    else
        txt = string(sprintf('%.2e',lambda));
    end
end

function txt = snr_text_local(snr_val)
    if isinf(snr_val)
        txt = 'Infinit';
    else
        txt = sprintf('%.1f dB',snr_val);
    end
end

function txt = metoda_finala_scurta(c)
    if isnan(c.idx_final)
        txt = "nicio metoda acceptata";
    else
        txt = nume_metoda_scurt(c.rezultate(c.idx_final).nume);
    end
end

function txt = nume_caz_scurt(nume)
    s = lower(string(nume));
    if contains(s,"fara")
        txt = "Fara zgomot";
    elseif contains(s,"cu")
        txt = "Cu zgomot";
    else
        txt = string(nume);
    end
end

function txt = nume_metoda_scurt(nume)
    s = lower(string(nume));
    if contains(s,"baza") || contains(s,"spatiu") || contains(s,"svd") && ~contains(s,"direct")
        txt = "Metoda de baza-SVD";
    elseif contains(s,"direct")
        txt = "Metoda de baza directa";
    elseif contains(s,"deriv")
        txt = "Derivata filtrata + Tikhonov";
    elseif contains(s,"integr")
        txt = "Forma integrala + Tikhonov";
    else
        txt = string(nume);
    end
end

function txt = rol_metoda_scurt(caz,metoda)
    c = lower(string(caz));
    m = lower(string(metoda));
    if contains(c,"fara")
        txt = "caz ideal";
    elseif contains(m,"direct")
        txt = "test sensibilitate";
    elseif contains(m,"deriv")
        txt = "varianta intermediara";
    elseif contains(m,"integr")
        txt = "varianta propusa";
    else
        txt = "metoda analizata";
    end
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

function txt = status_afisat(status)
    s = string(status);
    if s == "OK"
        txt = "OK";
    elseif s == "RESPINS_SPATIU_NUL"
        txt = "respins spatiu nul";
    elseif s == "RESPINS_NORMALIZARE"
        txt = "respins normalizare";
    elseif s == "RESPINS_REZIDUU"
        txt = "respins reziduu";
    elseif s == "RESPINS_ADMIS"
        txt = "respins admisibilitate";
    elseif s == "RESPINS_CL"
        txt = "respins bucla inchisa";
    elseif s == "RESPINS_SVD"
        txt = "respins SVD";
    else
        txt = lower(s);
    end
end

function label = eticheta_metoda_grafic(caz,metoda)
    label = nume_caz_scurt(caz) + ": " + nume_metoda_scurt(metoda);
end

function afiseaza_stack_dashboard_local(ME)
    try
        for k = 1:numel(ME.stack)
            fprintf('  eroare dashboard: %s, linia %d\n',ME.stack(k).name,ME.stack(k).line);
        end
    catch
    end
end

function afiseaza_eroare_dashboard(ME)
    afiseaza_stack_dashboard_local(ME);
end

