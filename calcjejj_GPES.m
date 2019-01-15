%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% This code is for the Development of Generalized Potenial%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%for Carbon Clusters %%%%%%%%%%%%%%%%%%%%%%%%%

function [JE,JtJ,normJE]=calcjejj_GPES(net_fc,...
    net_fr_CC,...
    net_ftheta_CCC,...
    net_fphi_CCCC,...
    rs, clust_size, Q, type, Ex)

% There will be 3 Jacobian matrices, one each for fc, fr, and ftheta

%%%% Jacobian matrix for fr
%J_fr will have dimensions of Q x ( 3N + 1 )
% First layer W (N x 1), b (N x 1), Second Layer W (1 x N ), b (1x1)

W1CC = net_fr_CC.IW{1,1};
W2CC = net_fr_CC.LW{2,1};
B1CC = net_fr_CC.b{1};
B2CC = net_fr_CC.b{2};


S1 = 25;
S2=1;
Total2BodyNumParam=3*S1+1;

W1CCC = net_ftheta_CCC.IW{1,1};
W2CCC = net_ftheta_CCC.LW{2,1};
B1CCC = net_ftheta_CCC.b{1};
B2CCC = net_ftheta_CCC.b{2};

S1 = 25;
Total3BodyNumParam=5*S1+1;

W1CCCC = net_fphi_CCCC.IW{1,1};
W2CCCC = net_fphi_CCCC.LW{2,1};
B1CCCC = net_fphi_CCCC.b{1};
B2CCCC = net_fphi_CCCC.b{2};

S1 = 5;
Total4BodyNumParam=8*S1+1;


% %%%%%%%
for iQ = 1:Q
    %     iQ
    jac(1:Total2BodyNumParam)=0;
    ciQ = 0;
    ciQCC = 0;
    ciQCCC=0;
    ciQCCCC=0;

    for i=1:1:clust_size(iQ)
        for j=i+1:1:clust_size(iQ)

            ciQCC = ciQCC+1;
            r_fr_CC(ciQCC) = rs(i,j,iQ); %#ok<AGROW>
            for k=j+1:1:clust_size(iQ)

                ciQCCC=ciQCCC+1;
                r_ftheta_CCC(1,ciQCCC)=rs(i,j,iQ); %#ok<AGROW>
                r_ftheta_CCC(2,ciQCCC)=rs(i,k,iQ); %#ok<AGROW>
                r_ftheta_CCC(3,ciQCCC)=rs(j,k,iQ); %#ok<AGROW>
                for l=k+1:1:clust_size(iQ)
                    ciQCCCC=ciQCCCC+1;
                    r_fphi_CCCC(1,ciQCCCC)=rs(i,j,iQ); %#ok<AGROW>
                    r_fphi_CCCC(2,ciQCCCC)=rs(i,k,iQ); %#ok<AGROW>
                    r_fphi_CCCC(3,ciQCCCC)=rs(i,l,iQ); %#ok<AGROW>
                    r_fphi_CCCC(4,ciQCCCC)=rs(j,k,iQ); %#ok<AGROW>
                    r_fphi_CCCC(5,ciQCCCC)=rs(j,l,iQ); %#ok<AGROW>
                    r_fphi_CCCC(6,ciQCCCC)=rs(k,l,iQ); %#ok<AGROW>
                end
            end

        end
    end


    % IMPORTANT NOTE: The size of 'ones' vector depends on the overall
    % combination within the cluster. Rightnow, it is dynamic only for this
    % case



    if ciQCC>0
        A1CC = logsig(W1CC*r_fr_CC + B1CC*ones(1,(ciQCC)*1));
        A2CC = W2CC*A1CC + B2CC*ones(1,(ciQCC)*1);
        % FIND JACOBIAN
        A1CC = kron(A1CC,ones(1,S2));
        D2CC = nnmdlin(A2CC);

        D1CC = nnmdlog(A1CC,D2CC,W2CC);
        jac1 = nnlmarq(kron(r_fr_CC,ones(1,S2)),D1CC);
        jac2 = nnlmarq(A1CC,D2CC);

        jac_fr_CC = [jac1,D1CC',jac2,D2CC']; % Make jacbian for current iQ configuration
        for j= 1: Total2BodyNumParam
            jac(1,j)= 0; %#ok<AGROW>
            for i = 1:(ciQCC)
                jac(1,j) = jac(1,j)+jac_fr_CC(i,j);  %#ok<AGROW> % sum jacobian for current iQ configuration along all rows to get 1 row

            end
        end
    end

    for j= 1:Total2BodyNumParam
        J_fr(iQ,j)= jac(1,j); %#ok<AGROW>
    end

    if ciQCCC>0
        A1CCC = logsig(W1CCC*r_ftheta_CCC + B1CCC*ones(1,(ciQCCC)*1));
        A2CCC = W2CCC*A1CCC+B2CCC*ones(1,(ciQCCC)*1);

        % FIND JACOBIAN
        A1CCC = kron(A1CCC,ones(1,S2));
        D2CCC = nnmdlin(A2CCC);

        D1CCC = nnmdlog(A1CCC,D2CCC,W2CCC);
        jac1 = nnlmarq(kron(r_ftheta_CCC,ones(1,S2)),D1CCC);
        jac2 = nnlmarq(A1CCC,D2CCC);

        jac_ftheta_CCC = [jac1,D1CCC',jac2,D2CCC'];
        %%%%
        for j= 1: Total3BodyNumParam
            jac(1,j)= 0;
            for i = 1:(ciQCCC)
                jac(1,j) = jac(1,j)+jac_ftheta_CCC(i,j);  % sum jacobian for current iQ configuration along all rows to get 1 row

            end
        end
    end


    for j= 1:Total3BodyNumParam
        J_ftheta(iQ,j)= jac(1,j);
    end
    
    if ciQCCCC>0
        A1CCCC = logsig(W1CCCC*r_fphi_CCCC + B1CCCC*ones(1,(ciQCCCC)*1));
        A2CCCC = W2CCCC*A1CCCC+B2CCCC*ones(1,(ciQCCCC)*1);

        % FIND JACOBIAN
        A1CCCC = kron(A1CCCC,ones(1,S2));
        D2CCCC = nnmdlin(A2CCCC);

        D1CCCC = nnmdlog(A1CCCC,D2CCCC,W2CCCC);
        jac1 = nnlmarq(kron(r_fphi_CCCC,ones(1,S2)),D1CCCC);
        jac2 = nnlmarq(A1CCCC,D2CCCC);

        jac_fphi_CCCC = [jac1,D1CCCC',jac2,D2CCCC'];
        %%%%
        for j= 1: Total4BodyNumParam
            jac(1,j)= 0;
            for i = 1:(ciQCCCC)
                jac(1,j) = jac(1,j)+jac_fphi_CCCC(i,j);  % sum jacobian for current iQ configuration along all rows to get 1 row

            end
        end
    end
    
    for j= 1:Total4BodyNumParam
        J_fphi(iQ,j)= jac(1,j);
    end

end
J=horzcat(J_fr,J_ftheta, J_fphi);
JE=J'*Ex;
JtJ=J'*J;
normJE=sqrt(JE'*JE);



i=1; %dummy line for breakpoint




