% function [Fx, Fy, Fz]=calc_force_GPES(net_fc, ...
%                                       net_fr_CC, ...
%                                       net_ftheta_CCC,...
%                                       rs, V, Q, clust_size, type)


xq=[0.01140371
    1.78140365
    1.78140339
    2.66640154
    0.0113991
    0.89639689
    0.89639689
    2.66639666

    ]    ;
yq=[0.01140356
    1.78140347
    0.01139922
    0.89639666
    1.78140353
    2.66640133
    0.89639696
    2.6663965


    ]    ;
zq=[0.01140356
    0.01139863
    1.78140347
    0.89639673
    1.78140339
    0.89639666
    2.66640136
    2.66639597


    ]    ;
W1CC=net_fr_CC.IW{1,1};
b1CC=net_fr_CC.b{1,1};
W2CC=net_fr_CC.lw{2,1};
b2CC=net_fr_CC.b{2,1};

W1CCC=net_ftheta_CCC.IW{1,1};
b1CCC=net_ftheta_CCC.b{1,1};
W2CCC=net_ftheta_CCC.lw{2,1};
b2CCC=net_ftheta_CCC.b{2,1};

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
    
    dV_dx(1:clust_size(iQ))=0.0;
    dV_dy(1:clust_size(iQ))=0.0;
    dV_dz(1:clust_size(iQ))=0.0;
    
    dV_drij=0.0;
    dV_drik=0.0;
    dV_drjk=0.0;

    fc=1;%sim(net_fc,rij);
    for i=1:1:clust_size(iQ)
        for j=i+1:1:clust_size(iQ)

            rij=r(i,j);


            frCC=frCC+sim(net_fr_CC,rij);
            %derivative calculations
            n1CC=W1CC*rij+b1CC;
            a1CC=logsig(n1CC);
            %             frCCtemp=W2CC*a1CC+b2CC;
            da1_dn1CC=logsig('dn',n1CC,a1CC);
            dV_drij=W2CC*(W1CC.*da1_dn1CC);
            dV_dx(i)=dV_dx(i)-dV_drij*(xq(i)-xq(j))/rij;
            dV_dy(i)=dV_dy(i)-dV_drij*(yq(i)-yq(j))/rij;
            dV_dz(i)=dV_dz(i)-dV_drij*(zq(i)-zq(j))/rij;

            dV_dx(j)=dV_dx(j)+dV_drij*(xq(i)-xq(j))/rij;
            dV_dy(j)=dV_dy(j)+dV_drij*(yq(i)-yq(j))/rij;
            dV_dz(j)=dV_dz(j)+dV_drij*(zq(i)-zq(j))/rij;


            for k=j+1:1:clust_size(iQ)

                rik=r(i,k);
                rjk=r(j,k);
                rv=[rij; rik ;rjk];
                fthetaCCC=fthetaCCC+sim(net_ftheta_CCC,rv);
                %derivative calculations
                n1CCC=W1CCC*rv+b1CCC;
                a1CCC=logsig(n1CCC);
                da1_dn1CCC=logsig('dn',n1CCC,a1CCC);
                dV_drij=W2CCC*(W1CCC(:,1).*da1_dn1CCC);
                dV_dx(i)=dV_dx(i)-dV_drij*(xq(i)-xq(j))/rij;
                dV_dy(i)=dV_dy(i)-dV_drij*(yq(i)-yq(j))/rij;
                dV_dz(i)=dV_dz(i)-dV_drij*(zq(i)-zq(j))/rij;

                dV_dx(j)=dV_dx(j)+dV_drij*(xq(i)-xq(j))/rij;
                dV_dy(j)=dV_dy(j)+dV_drij*(yq(i)-yq(j))/rij;
                dV_dz(j)=dV_dz(j)+dV_drij*(zq(i)-zq(j))/rij;

                dV_drik=W2CCC*(W1CCC(:,2).*da1_dn1CCC);

                dV_dx(i)=dV_dx(i)-dV_drik*(xq(i)-xq(k))/rik;
                dV_dy(i)=dV_dy(i)-dV_drik*(yq(i)-yq(k))/rik;
                dV_dz(i)=dV_dz(i)-dV_drik*(zq(i)-zq(k))/rik;

                dV_dx(k)=dV_dx(k)+dV_drik*(xq(i)-xq(k))/rik;
                dV_dy(k)=dV_dy(k)+dV_drik*(yq(i)-yq(k))/rik;
                dV_dz(k)=dV_dz(k)+dV_drik*(zq(i)-zq(k))/rik;

                dV_drjk=W2CCC*(W1CCC(:,3).*da1_dn1CCC);

                dV_dx(j)=dV_dx(j)-dV_drjk*(xq(j)-xq(k))/rjk;
                dV_dy(j)=dV_dy(j)-dV_drjk*(yq(j)-yq(k))/rjk;
                dV_dz(j)=dV_dz(j)-dV_drjk*(zq(j)-zq(k))/rjk;

                dV_dx(k)=dV_dx(k)+dV_drjk*(xq(j)-xq(k))/rjk;
                dV_dy(k)=dV_dy(k)+dV_drjk*(yq(j)-yq(k))/rjk;
                dV_dz(k)=dV_dz(k)+dV_drjk*(zq(j)-zq(k))/rjk;

            end

        end
    end



    Vhat_iQ=fc*(frCC + fthetaCCC);

    Vhat(iQ)=Vhat_iQ; %#ok<AGROW>
    Ex(iQ,1)=V(iQ)-Vhat(iQ);  %#ok<AGROW>
    %SSE
    perft=Ex(iQ,1)*Ex(iQ,1);



end

i=1;
