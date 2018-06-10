% Xu 2007 (c)
% The Elmira Conflict for Status QUO (new)
%______________________________________
% parameters 

delta = 1;  % episode, delta = |UA|
n = 3; % the number of decision maker
m = 9; % the number of state

J = zeros(m,m,n);
Jp = zeros(m,m,n);
Pp = zeros(m,m,n);
Pne = zeros(m,m,n);
Phne = zeros(m,m,n);

payoff = zeros(m,n);

Mh = zeros(m,m,n);
Mhp = zeros(m,m,n);


M = zeros(m,m,n,n,delta);
Mp = zeros(m,m,n,n,delta);

MQ = zeros(m,m,n,n,delta);
MQp = zeros(m,m,n,n,delta);
MQt = zeros(m,m,n,n);
MQtp = zeros(m,m,n,n);
MQt1=zeros(m,m,n);
MQt1p=zeros(m,m,n);
                           
                           
Mnash = zeros(n,m);
Mgmr = zeros(m,m,n);
Msmr = zeros(m,m,n);
Mseq = zeros(m,m,n);

E = ones(m,m); Ii = eye(m,m);

%___________________________________________
% set initial value
J(1,2,1) = 1; J(3,4,1) = 1; J(5,6,1) = 1; J(7,8,1) = 1; 

payoff(:,1) = [4,3,8,7,5,2,9,6,1]';

J(1,3,2) = 1; J(1,9,2) = 1; J(2,4,2) = 1; J(2,9,2) = 1;
J(3,9,2) = 1; J(4,9,2) = 1; J(5,7,2) = 1; J(5,9,2) = 1;
J(6,8,2) = 1; J(6,9,2) = 1; J(7,9,2) = 1; J(8,9,2) = 1;

payoff(:,2) = [9,2,4,8,6,1,3,7,5]';

J(1,5,3) = 1; J(5,1,3) = 1; J(2,6,3) = 1; J(6,2,3) = 1; 
J(7,3,3) = 1; J(3,7,3) = 1; J(4,8,3) = 1; J(8,4,3) = 1;

payoff(:,3) = [6,2,8,3,7,4,9,5,1]';



%______________________________________________
% calculate Pp

for i = 1:n
    for j = 1:m
        for k = 1:m
            Pp(j,k,i) = sign(0.5*(abs(payoff(k,i)-payoff(j,i)) + payoff(k,i)-payoff(j,i)));
        end
    end
end

%______________________________________________
% calculate Pne

for i = 1:n
    Pne(:,:,i) = ones(m,m) - Pp(:,:,i) - eye(m,m);
end

%______________________________________________
% calculate Jp
for i = 1:n
    Jp(:,:,i) = J(:,:,i).*Pp(:,:,i);
end
Jp = sign(Jp);

%______________________________________________
% calculate M, Mp
for i = 1:n
    for j = 1:n
        if j ~= i
            M(:,:,i,j,1) = J(:,:,j);     % initial
            Mp(:,:,i,j,1) = Jp(:,:,j);
        else
            M(:,:,i,j,1) = zeros(m,m);    
            Mp(:,:,i,j,1) = zeros(m,m);
        end
    end
end

for k = 2:delta
    for i = 1:n
        
        for j = 1:n
            
            %if j ~= i
                tempM = zeros(m,m); 
                tempMp = zeros(m,m);

                for t = 1:n
                    %if (t ~= i) && (t ~= j)
                    if t ~= j
                        tempM = tempM + M(:,:,i,t,k-1);  
                        tempMp = tempMp + Mp(:,:,i,t,k-1);
                    end
                end 

                M(:,:,i,j,k) = sign(J(:,:,j)*tempM);
                Mp(:,:,i,j,k) = sign(Jp(:,:,j)*tempMp);
            %end
            
        end
        
    end
end

%_____________________________________________
% calculate MQ, MQp
for i = 1:n
    for j = 1:n
%         if j ~= i
            MQ(:,:,i,j,1) = J(:,:,j);     % initial
            MQp(:,:,i,j,1) = Jp(:,:,j);
%         else
%             MQ(:,:,i,j,1) = zeros(m,m);    
%             MQp(:,:,i,j,1) = zeros(m,m);
%         end
    end
end

for k = 2:delta
    for i = 1:n
        
        for j = 1:n
            
%             if j ~= i
                tempMQ = zeros(m,m); 
                tempMQp = zeros(m,m);

                for t = 1:n
                     %if (t ~= i) && (t ~= j)
                    if t ~= j
                        tempMQ = tempMQ + MQ(:,:,i,t,k-1);  
                        tempMQp = tempMQp + MQp(:,:,i,t,k-1);
                    end
%                 end 

                MQ(:,:,i,j,k) = sign(tempMQ*J(:,:,j));
                MQp(:,:,i,j,k) = sign(tempMQp*Jp(:,:,j));
            end
            
        end
        
    end
end

% calculate MQh, MQhp

        
            for i=1:n
                for j=1:n
                    for k = 1:delta
                    MQt(:,:,i,j) = MQt(:,:,i,j)+MQ(:,:,i,j,k);
                    MQtp(:,:,i,j) = MQtp(:,:,i,j)+ MQp(:,:,i,j,k) ;
                    end
                end
           end
        
% calculate MQt1, MQt1p

        
            
                for j=1:n
                    for i=1:n   
                        for k = 1:delta
                           MQt1(:,:,j) = MQt1(:,:,j)+MQt(:,:,i,j);
                           MQt1p(:,:,j) = MQt1p(:,:,j)+ MQtp(:,:,i,j) ;
                        end
                    end
                    MQt1(:,:,j) = sign(MQt1(:,:,j));
                    MQt1p(:,:,j) = sign(MQt1p(:,:,j));
                end

    
                      
    