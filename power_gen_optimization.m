%% CE 295 - Energy Systems and Control
%   Optimal Economic Dispatch in Distribution Feeders with Renewables
%   Prof. Moura
%   Last updated: February 20, 2018

% TSANG_NATHAN_HW3.m

clear; close all;
fs = 15;    % Font Size for plots

%% 13 Node IEEE Test Feeder Parameters

%%% Node (aka Bus) Data
% l_j^P: Active power consumption [MW]
l_P = [0; 0.2; 0; 0.4; 0.17; 0.23; 1.155; 0; 0.17; 0.843; 0; 0.17; 0.128];

% l_j^Q: Reactive power consumption [MVAr]
l_Q = [0; 0.116; 0; 0.29; 0.125; 0.132; 0.66; 0; 0.151; 0.462; 0; 0.08; 0.086];

% l_j^S: Apparent power consumption [MVA]
l_S = sqrt(l_P.^2 + l_Q.^2);

% s_j,max: Maximal generating power [MW]
s_max = [5; 0; 0; 3; 0; 0; 0; 0; 0; 3; 0; 0; 0];

% c_j: Marginal generation cost [USD/MW]
c = [100; 0; 0; 150; 0; 0; 0; 0; 0; 50; 0; 0; 0];

% V_min, V_max: Minimum and maximum nodal voltages [V]
v_min = 0.95;
v_max = 1.05;

%%% Edge (aka Line) Data
% r_ij: Resistance [p.u.]
r = xlsread('HW3_Data.xlsx', 'Line-Data', 'B6:N18');

% x_ij: Reactance [p.u.]
x = xlsread('HW3_Data.xlsx', 'Line-Data', 'B24:N36');

% I_max_ij: Maximal line current [p.u.]
I_max = xlsread('HW3_Data.xlsx', 'Line-Data', 'B42:N54');

% A_ij: Adjacency matrix; A_ij = 1 if i is parent of j
A = xlsread('HW3_Data.xlsx', 'Line-Data', 'B60:N72');


%%% Set Data (add +1 everywhere for Matlab indexing)
% \rho(j): Parent node of node j
rho = [0; 0; 1; 2; 1; 4; 1; 6; 6; 8; 6; 10; 10]+1;

% Node index set - CREATED BY NATE
j = [0 1 2 3 4 5 6 7 8 9 10 11 12];

%% Problem 1

% Plot active and reactive power consumption
figure(1);
bar(j, [l_P,l_Q]);
legend('Real power (MW)','Reactive power (MVAR)')
xlabel('Node index number')
ylabel('Power consumption')
set(gca,'FontSize',fs)

%% Problem 2

% Assumptions:
%   - Disregard the entire network diagram
%   - Balance supply and demand, without any network considerations
%   - Goal is to minimize generation costs, given by c^T s

% Solve with CVX
cvx_begin
    variables s(13,1) p(13,1) q(13,1)   % declare your optimization variables here
    minimize(c'*s)                      % objective function here
    subject to                          % constraints 
        % Balance power generation with power consumption
        sum(p) == sum(l_P);
        sum(q) == sum(l_Q);
        
        % Loop over each node
        for jj = 1:13
            
            % Non-negative power generation
            q(jj) >= 0;
            p(jj) >= 0;
            
            % Compute apparent power from active & reactive power
            norm([p(jj) q(jj)],2) <= s(jj);
            
        end
        
        % Apparent Power Limits
        s <= s_max;
        
cvx_end

% Output Results
fprintf(1,'------------------- PROBLEM 2 --------------------\n');
fprintf(1,'--------------------------------------------------\n');
fprintf(1,'Minimum Generating Cost : %4.2f USD\n',cvx_optval);
fprintf(1,'\n');
fprintf(1,'Node 0 [Grid]  Gen Power : p_0 = %1.3f MW | q_0 = %1.3f MW | s_0 = %1.3f MW\n',p(1),q(1),s(1));
fprintf(1,'Node 3 [Gas]   Gen Power : p_3 = %1.3f MW | q_3 = %1.3f MW | s_3 = %1.3f MW\n',p(4),q(4),s(4));
fprintf(1,'Node 9 [Solar] Gen Power : p_9 = %1.3f MW | q_9 = %1.3f MW | s_9 = %1.3f MW\n',p(10),q(10),s(10));
fprintf(1,'\n');
fprintf(1,'Total active power   : %1.3f MW   consumed | %1.3f MW   generated\n',sum(l_P),sum(p));
fprintf(1,'Total reactive power : %1.3f MVAr consumed | %1.3f MVAr generated\n',sum(l_Q),sum(q));
fprintf(1,'Total apparent power : %1.3f MVA  consumed | %1.3f MVA  generated\n',sum(l_S),sum(s));

%% Problem 3

% Assumptions:
%   - Disregard L_ij, the squared magnitude of complex line current
%   - Disregard nodal voltage equation
%   - Disregard nodal voltage limits
%   - Disregard maximum line current
%   - Goal is to minimize generation costs, given by c^T s

% Solve with CVX
cvx_begin
    variables p(13,1) q(13,1) s(13,1) P(13,13) Q(13,13)
    dual variable mu_s
    minimize(c' * s)
       subject to
    
        % Boundary condition for power line flows
        P( 1 , 1 ) == 0;
        Q( 1 , 1 ) == 0;
        
        % Loop over each node
        for jj = 1:13
            
            % Parent node, i = \rho(j)
            ii = rho(jj);    
            
            % Line Power Flows
            P(ii,jj) == l_P(jj) - p(jj) + sum(A(jj,:).* P(jj,:));
            Q(ii,jj) == l_Q(jj) - q(jj) + sum(A(jj,:).* Q(jj,:));
            
            % Compute apparent power from active & reactive power
            norm([p(jj) q(jj)],2) <= s(jj);        
            q(jj) >= 0;
            p(jj) >= 0;
            
        end
        
        % Apparent Power Limits
        s <= s_max : mu_s;
        
cvx_end

% Output Results
fprintf(1,'------------------- PROBLEM 3 --------------------\n');
fprintf(1,'--------------------------------------------------\n');
fprintf(1,'Minimum Generating Cost : %4.2f USD\n',cvx_optval);
fprintf(1,'\n');
fprintf(1,'Node 0 [Grid]  Gen Power : p_0 = %1.3f MW | q_0 = %1.3f MW | s_0 = %1.3f MW || mu_s0 = %3.0f USD/MW\n',p(1),q(1),s(1),mu_s(1));
fprintf(1,'Node 3 [Gas]   Gen Power : p_3 = %1.3f MW | q_3 = %1.3f MW | s_3 = %1.3f MW || mu_s3 = %3.0f USD/MW\n',p(4),q(4),s(4),mu_s(4));
fprintf(1,'Node 9 [Solar] Gen Power : p_9 = %1.3f MW | q_9 = %1.3f MW | s_9 = %1.3f MW || mu_s9 = %3.0f USD/MW\n',p(10),q(10),s(10),mu_s(10));
fprintf(1,'\n');
fprintf(1,'Total active power   : %1.3f MW   consumed | %1.3f MW   generated\n',sum(l_P),sum(p));
fprintf(1,'Total reactive power : %1.3f MVAr consumed | %1.3f MVAr generated\n',sum(l_Q),sum(q));
fprintf(1,'Total apparent power : %1.3f MVA  consumed | %1.3f MVA  generated\n',sum(l_S),sum(s));


%% Problem 4

% Assumptions:
%   - Add back all previously disregarded terms and constraints
%   - Relax squared line current equation into inequality
%   - Goal is to minimize generation costs, given by c^T s

% Solve with CVX
cvx_begin
    variables p(13,1) q(13,1) s(13,1) P(13,13) Q(13,13) L(13,13) V(13,1)
    dual variables mu_s mu_L mu_vmin mu_vmax
    minimize(c'*s)
    subject to
    
        % Boundary condition for power line flows
        P( 1 , 1 ) == 0;
        Q( 1 , 1 ) == 0;
        
        % Boundary condition for squared line current
        L( 1 , 1 ) == 0;
        
        % Fix node 0 voltage to be 1 "per unit" (p.u.)
        V(1) == 1;
        
        % Loop over each node
        for jj = 1:13
            
            % Parent node, i = \rho(j)
            ii = rho(jj);
            
            % Line Power Flows
            P(ii,jj) == l_P(jj) - p(jj) + sum(A(jj,:).* P(jj,:)) + r(ii,jj)*L(ii,jj);
            Q(ii,jj) == l_Q(jj) - q(jj) + sum(A(jj,:).* Q(jj,:)) + x(ii,jj)*L(ii,jj);
            
            % Nodal voltage
            V(jj) == V(ii)+(r(ii,jj)^2 + x(ii,jj)^2) * L(ii,jj) - 2*(r(ii,jj) * P(ii,jj) + x(ii,jj) * Q(ii,jj));
            
            % Squared current magnitude on lines
            L(ii,jj) >= quad_over_lin(P(ii,jj),V(jj))+ quad_over_lin(Q(ii,jj),V(jj));
            
            % Compute apparent power from active & reactive power
            norm([p(jj) q(jj)],2) <= s(jj);        
            q(jj) >= 0;
            p(jj) >= 0;
            
        end
        
        % Squared line current limits
        L <= I_max.^2 : mu_L;
            
        % Nodal voltage limits
        V <= v_max.^2 : mu_vmax;
        V >= v_min.^2 : mu_vmin;
        
        % Apparent Power Limits
        s <= s_max : mu_s;
        
cvx_end

% Output Results
fprintf(1,'------------------- PROBLEM 4 --------------------\n');
fprintf(1,'--------------------------------------------------\n');
fprintf(1,'Minimum Generating Cost : %4.2f USD\n',cvx_optval);
fprintf(1,'\n');
fprintf(1,'Node 0 [Grid]  Gen Power : p_0 = %1.3f MW | q_0 = %1.3f MW | s_0 = %1.3f MW || mu_s0 = %3.0f USD/MW | mu_vmax0 = %3.2f p.u.| mu_vmin0 = %3.2f p.u./MW\n',p(1),q(1),s(1),mu_s(1),mu_vmax(1),mu_vmin(1));
fprintf(1,'Node 3 [Gas]   Gen Power : p_3 = %1.3f MW | q_3 = %1.3f MW | s_3 = %1.3f MW || mu_s3 = %3.0f USD/MW | mu_vmax3 = %3.2f p.u.| mu_vmin3 = %3.2f p.u./MW\n',p(4),q(4),s(4),mu_s(4),mu_vmax(4),mu_vmin(4));
fprintf(1,'Node 9 [Solar] Gen Power : p_9 = %1.3f MW | q_9 = %1.3f MW | s_9 = %1.3f MW || mu_s9 = %3.0f USD/MW | mu_vmax9 = %3.2f p.u.| mu_vmin9 = %3.2f p.u./MW\n',p(10),q(10),s(10),mu_s(10),mu_vmax(10),mu_vmin(10));
fprintf(1,'\n');
fprintf(1,'Total active power   : %1.3f MW   consumed | %1.3f MW   generated\n',sum(l_P),sum(p));
fprintf(1,'Total reactive power : %1.3f MVAr consumed | %1.3f MVAr generated\n',sum(l_Q),sum(q));
fprintf(1,'Total apparent power : %1.3f MVA  consumed | %1.3f MVA  generated\n',sum(l_S),sum(s));
fprintf(1,'\n');
for jj = 1:13
    fprintf(1,'Node %2.0f Voltage : %1.3f p.u.\n',jj,sqrt(V(jj)));
end


%% Problem 5

% Assumptions:
%   - Assume solar generator at node 9 has uncertain power capacity
%   - Goal is to minimize generation costs, given by c^T s, in face of uncertainty

% Solve with CVX

% define new variables
a_bar = [-1.25; -1.25; 1];
E = diag([0.25 0.25 0]);
b = 0;

cvx_begin
    variables p(13,1) q(13,1) s(13,1) P(13,13) Q(13,13) L(13,13) V(13,1) y(3,1)
    dual variables mu_s mu_L mu_vmin mu_vmax
    minimize(c'*s)
    subject to
    
        % Boundary condition for power line flows
        P( 1 , 1 ) == 0;
        Q( 1 , 1 ) == 0;
        
        % Boundary condition for squared line current
        L( 1 , 1 ) == 0;
        
        % Fix node 0 voltage to be 1 "per unit" (p.u.)
        V(1) == 1;
        
        % Loop over each node
        for jj = 1:13
            % Parent node, i = \rho(j)
            ii = rho(jj);
            if jj == 10
                a_bar'*y + norm([E'*y],2) <= b;
                y(1:2,1) >= 0;
                y(1:2,1) <= 1;
                s(jj) == y(3);  
                
                % Line Power Flows
                P(ii,jj) == l_P(jj) - p(jj) + sum(A(jj,:).* P(jj,:)) + r(ii,jj)*L(ii,jj);
                Q(ii,jj) == l_Q(jj) - q(jj) + sum(A(jj,:).* Q(jj,:)) + x(ii,jj)*L(ii,jj);

                % Nodal voltage
                V(jj) == V(ii)+(r(ii,jj)^2 + x(ii,jj)^2) * L(ii,jj) - 2*(r(ii,jj) * P(ii,jj) + x(ii,jj) * Q(ii,jj));

                % Squared current magnitude on lines
                L(ii,jj) >= quad_over_lin(P(ii,jj),V(jj))+ quad_over_lin(Q(ii,jj),V(jj));

                % Compute apparent power from active & reactive power
                norm([p(jj) q(jj)],2) <= s(jj);        
                q(jj) >= 0;
                p(jj) >= 0;

            else 

            % Line Power Flows
            P(ii,jj) == l_P(jj) - p(jj) + sum(A(jj,:).* P(jj,:)) + r(ii,jj)*L(ii,jj);
            Q(ii,jj) == l_Q(jj) - q(jj) + sum(A(jj,:).* Q(jj,:)) + x(ii,jj)*L(ii,jj);

            % Nodal voltage
            V(jj) == V(ii)+(r(ii,jj)^2 + x(ii,jj)^2) * L(ii,jj) - 2*(r(ii,jj) * P(ii,jj) + x(ii,jj) * Q(ii,jj));

            % Squared current magnitude on lines
            L(ii,jj) >= quad_over_lin(P(ii,jj),V(jj))+ quad_over_lin(Q(ii,jj),V(jj));

            % Compute apparent power from active & reactive power
            norm([p(jj) q(jj)],2) <= s(jj);        
            q(jj) >= 0;
            p(jj) >= 0;
            end
            
        end
        
        % Squared line current limits
        L <= I_max.^2 : mu_L;
            
        % Nodal voltage limits
        V <= v_max.^2 : mu_vmax;
        V >= v_min.^2 : mu_vmin;
        
        % Apparent Power Limits
        s <= s_max : mu_s;
    
cvx_end

% Output Results
fprintf(1,'------------------- PROBLEM 5 --------------------\n');
fprintf(1,'--------------------------------------------------\n');
fprintf(1,'Minimum Generating Cost : %4.2f USD\n',cvx_optval);
fprintf(1,'\n');
fprintf(1,'Node 0 [Grid]  Gen Power : p_0 = %1.3f MW | q_0 = %1.3f MW | s_0 = %1.3f MW || mu_s0 = %3.0f USD/MW | mu_vmax0 = %3.2f p.u.| mu_vmin0 = %3.2f p.u./MW\n',p(1),q(1),s(1),mu_s(1),mu_vmax(1),mu_vmin(1));
fprintf(1,'Node 3 [Gas]   Gen Power : p_3 = %1.3f MW | q_3 = %1.3f MW | s_3 = %1.3f MW || mu_s3 = %3.0f USD/MW | mu_vmax3 = %3.2f p.u.| mu_vmin3 = %3.2f p.u./MW\n',p(4),q(4),s(4),mu_s(4),mu_vmax(4),mu_vmin(4));
fprintf(1,'Node 9 [Solar] Gen Power : p_9 = %1.3f MW | q_9 = %1.3f MW | s_9 = %1.3f MW || mu_s9 = %3.0f USD/MW | mu_vmax9 = %3.2f p.u.| mu_vmin9 = %3.2f p.u./MW\n',p(10),q(10),s(10),mu_s(10),mu_vmax(10),mu_vmin(10));
fprintf(1,'\n');
fprintf(1,'Total active power   : %1.3f MW   consumed | %1.3f MW   generated\n',sum(l_P),sum(p));
fprintf(1,'Total reactive power : %1.3f MVAr consumed | %1.3f MVAr generated\n',sum(l_Q),sum(q));
fprintf(1,'Total apparent power : %1.3f MVA  consumed | %1.3f MVA  generated\n',sum(l_S),sum(s));
fprintf(1,'\n');
for jj = 1:13
    fprintf(1,'Node %2.0f Voltage : %1.3f p.u.\n',jj,sqrt(V(jj)));
end
