% compressorModel
% 11-stage axial compressor model compatible with burnerModel and turbineModel
%
% x(1) = overall compressor pressure ratio pi_c

%%
% 1) Inputs and constants

pi_c = x(1);

gamma = getOr(data,'gamma_03',1.4);
cp    = getOr(data,'Cp_03',1000);       % J/kg-K
eta_c = getOr(data,'eta_c',0.90);

T0 = getOr(data,'T0',220);              % K
P0 = getOr(data,'P0',23000);            % Pa
V0 = getOr(data,'V0',250);              % m/s

mdot3 = getOr(data,'mdot3',getOr(data,'mdot_core',[]));

% FIXED number of stages (11 for project requirement)
N_stages = 11;

%% 
% 2) Compressor inlet (station 2)

R = cp*(gamma-1)/gamma;
a0 = sqrt(gamma*R*T0);
M0 = V0/a0;

Tt2 = T0*(1 + (gamma-1)/2*M0^2);
Pt2 = P0*(1 + (gamma-1)/2*M0^2)^(gamma/(gamma-1));

%% 
% 3) Per-stage pressure ratio

pi_stage = pi_c^(1/N_stages);

%%
% 4) Allocate arrays (stage-by-stage)

Tt = zeros(N_stages+1,1);
Pt = zeros(N_stages+1,1);
W_stage = zeros(N_stages,1);

Tt(1) = Tt2;
Pt(1) = Pt2;

%%
% 5) 11-stage compression loop

for i = 1:N_stages

    % Isentropic temperature ratio for this stage
    tau_is = pi_stage^((gamma-1)/gamma);

    % Actual temperature rise (efficiency included)
    dTt = (Tt(i)/eta_c)*(tau_is - 1);

    % Update totals
    Tt(i+1) = Tt(i) + dTt;
    Pt(i+1) = Pt(i) * pi_stage;

    % Stage specific work
    W_stage(i) = cp * dTt;
end

%% 
% 6) Compressor exit (station 3)

Tt3 = Tt(end);
Pt3 = Pt(end);

% Total specific work (sum of stages)
W_comp = sum(W_stage);

% Shaft power required
P_comp = [];
if ~isempty(mdot3)
    P_comp = mdot3 * W_comp;
end

%% 
% 7) Package outputs (for turbine & burner)

compressor.Tt2 = Tt2;
compressor.Pt2 = Pt2;
compressor.Tt3 = Tt3;
compressor.Pt3 = Pt3;

compressor.pi_c = pi_c;
compressor.N_stages = N_stages;

compressor.W_stage = W_stage;      % J/kg per stage
compressor.W_comp  = W_comp;       % J/kg total
compressor.mdot3   = mdot3;

if ~isempty(P_comp)
    compressor.P_required_shaft = P_comp;
end

%% 
% 8) Stage-by-stage results table (not necessary)

Stage = (1:N_stages).';
Tt_in  = Tt(1:end-1);
Tt_out = Tt(2:end);
Pt_in  = Pt(1:end-1);
Pt_out = Pt(2:end);
pi_s   = Pt_out ./ Pt_in;

compressor.stageTable = table( ...
    Stage, ...
    Tt_in, Tt_out, ...
    Pt_in, Pt_out, ...
    pi_s, ...
    W_stage, ...
    'VariableNames', ...
    {'Stage','Tt_in_K','Tt_out_K','Pt_in_Pa','Pt_out_Pa','pi_stage','W_stage_Jpkg'} );

compressor.success = true;

end

%%
% Helper function
% 
function v = getOr(s, field, defaultVal)
    if isstruct(s) && isfield(s,field) && ~isempty(s.(field))
        v = s.(field);
    else
        v = defaultVal;
    end
end
