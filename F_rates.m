function F = F_rates(m_vec, DynamicsParams)
% F(i,j) is the rate that a node will change from state i to state j given that it has m_vec(1) neighbours in state 1, m_vec(2) neighbours in state 2, ..., m_vec(n) neighbours in state n. DynamicsParams is a vector of model parameters (e.g., DynamicsParams = [beta mu] in the case of SIS).


% % SI/SIS
% beta = DynamicsParams(1); % Infection rate
% mu   = DynamicsParams(2); % Recovery rate (set = 0 for SI)
% F = [0 beta*m_vec(2);
%     mu 0];

% % Cooperative SIS (Simple)
% beta   = DynamicsParams(1); % Infection rate
% lambda = DynamicsParams(2); % Accentuation parameter
% F = [0 beta*(m_vec(2) + m_vec(4)) beta*(m_vec(3) + m_vec(4)) 0;
%     1 0 0 beta*lambda*(m_vec(3) + m_vec(4));
%     1 0 0 beta*lambda*(m_vec(2) + m_vec(4));
%     0 1 1 0];

% Fredrickson Andersen model
T = DynamicsParams(1);  % temperature
f = DynamicsParams(2);  % facilitation parameter
spin_down = m_vec(1) + m_vec(3);
if(spin_down >= f)
    F = [0 0 0 1; 0 0 exp(-1/T) 0; 0 0 0 1; 0 0 exp(-1/T) 0];
else
    F = zeros(4);
end


% % % Cooperative SIS (Sanz et al. "Dynamics of Interacting Diseases",Phys Rev X, 2014)
% % Baseline parameters
% lambda1 = DynamicsParams(1); % Baseline infectiousness disease 1
% lambda2 = DynamicsParams(2); % Baseline infectiousness disease 2
% mu1 = DynamicsParams(3);  % Baseline recovery rate disease 1
% mu2 = DynamicsParams(4);  % Baseline recovery rate disease 2
% % Variational parameters (See Table 1 in Sanz et al.)
% beta1a = DynamicsParams(5);
% beta2a = DynamicsParams(6);
% beta1b = DynamicsParams(7);
% beta2b = DynamicsParams(8);
% eta1 = DynamicsParams(9);
% eta2 = DynamicsParams(10);
% F = [0  lambda1*(m_vec(2)+beta1b*m_vec(4)) lambda2*(m_vec(3)+beta2b*m_vec(4)) 0;
%     mu1 0         0        beta2a*lambda2*(m_vec(3)+beta2b*m_vec(4));
%     mu2 0         0        beta1a*lambda1*(m_vec(2)+beta1b*m_vec(4));
%     0   eta2*mu2  eta1*mu1 0];

% % Bass diffusion model
% c = DynamicsParams(1);
% d = DynamicsParams(2);
% F = [0 d*m_vec(2) + c;
%     0 0];

% % Majority Vote
% Q = DynamicsParams(1);
% if(m_vec(1) > m_vec(2))
%     F = [0 Q;
%         1-Q 0];
% elseif(m_vec(1) == m_vec(2))
%     F = [0 0.5;
%         0.5 0];
% else
%     F = [0 1-Q;
%         Q 0];
% end

% % % Ising Glauber
% J = DynamicsParams(1);
% T = DynamicsParams(2);
% F = [0 1/(1+exp(2*J*(m_vec(1)-m_vec(2))/T));
%     1/(1+exp(2*J*(m_vec(2)-m_vec(1))/T)) 0];


