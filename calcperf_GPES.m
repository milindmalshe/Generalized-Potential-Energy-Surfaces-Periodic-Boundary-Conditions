%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% This code is for the Development of Generalized Potenial%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%for CArbon clusters with PBC %%%%%%%%%%%%%%%%

function [perf,Ex,Vhat]=calcperf_GPES(net_fc, ...
    net_fr_CC, ...
    net_ftheta_CCC,...
    rs, V, Q, clust_size, type)

perf=0.0;

for iQ=1:1:Q

    %converting rs to my format
    for i=1:1:clust_size(iQ)
        for j=i+1:1:clust_size(iQ)
            r(i,j)=rs(i,j,iQ); %#ok<AGROW>
        end
    end


    %calculating potential energy for iQ configuration

    frCC=0.0;
    fthetaCCC=0.0;

    fc=1;%sim(net_fc,rij);
    for i=1:1:clust_size(iQ)
        for j=i+1:1:clust_size(iQ)

            rij=r(i,j);


            frCC=frCC+purelin(net_fr_CC.lw{2,1}*logsig(net_fr_CC.IW{1,1}*rij+net_fr_CC.b{1,1})+net_fr_CC.b{2,1});%sim(net_fr_CC,rij);

            for k=j+1:1:clust_size(iQ)

                rik=rs(i,k,iQ);
                rjk=rs(j,k,iQ);

                fthetaCCC=fthetaCCC+purelin(net_ftheta_CCC.lw{2,1}*logsig(net_ftheta_CCC.IW{1,1}*[rij;rik;rjk]+net_ftheta_CCC.b{1,1})+net_ftheta_CCC.b{2,1});%sim(net_ftheta_CCC,[rij; rik ;rjk]);

            end

        end
    end


    Vhat_iQ=fc*(frCC + fthetaCCC);

    Vhat(iQ)=Vhat_iQ; %#ok<AGROW>
    Ex(iQ,1)=V(iQ)-Vhat(iQ);  %#ok<AGROW>
    %SSE
    perf=perf+Ex(iQ,1)*Ex(iQ,1);



end


i=1;%dummy line

