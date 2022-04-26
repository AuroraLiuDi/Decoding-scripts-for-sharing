%%Title: code for replotting the figures that appeared in the paper
%%Author: Liu Di 201911061110
%%Important information:

%%1)The code is so long that it takes more than 10 minutes to run as a
%%whole. I strongly suggest the reader to run the code part by part,
%%in order to prevent getting too bored while waiting.

%%2)In Figure 1 part, you can get 2 images, one depicting the tuning curve,
%%and the other replicating the Fig1 in artical

%%3)For more information, please check the docx file in the folder
%% Figure1: 
% Preparation:simulating the neural response
x=(-pi:0.01:pi)*180/pi;%all the angles
fi=zeros(4,length(-pi:0.01:pi));%initialization
a=-0.14;%set a to original value
for i=1:4
    thetai=-135+90*(i-1);%favorite angle
    j=1;
    for theta=-pi:0.01:pi
        fi(i,j)=(cos(theta-thetai*pi/180)-a)/(1-a);%tuning curve
        %if you want to have a sense of what the simulation data looks
        %like, please run the codes below:
%         noise_sigma=0.1;
%         fi(i,j)=fi(i,j)+noise_sigma*randn(1);%the variability increased proportionaly to the firing rate
%         if fi(i,j)>=0
%             fi(i,j)=fi(i,j);
%         else
%             fi(i,j)=0;
%         end
        j=j+1;
    end
    x=(-pi:0.01:pi)*180/pi;
    y=fi(i,:);
    plot(x,y);
    hold on
    title("tuning curve");
    axis([-180,180,0,1]);
end

% 1.1Maximum likelyhood (method 2: correct!)
clear;
a=-0.14;
noise_sigma=0.1;%noise sigma=0.1*max(fi)=0.1
%simulation
for rep=1:1:100
    p=1;
    for theta=-pi:0.01:pi
        %simulation
        for i=1:4
            thetai=-3*pi/4+pi*(i-1)/2;
            vi(i)=(cos(theta-thetai)-a)/(1-a);
            if vi(i)>=0
                vi(i)=vi(i);
            else
                vi(i)=0;
            end
            ri(i)=vi(i)+noise_sigma*randn(1);
            if ri(i)>=0
                ri(i)=ri(i);
            else
                ri(i)=0;
            end
            %according to the likelyhood function, the likely hood can achieve the
            %highest value when:
            theta_est=theta-pi/2:0.01:theta+pi/2;
            for j=1:1:length(theta_est)
                fi(i)=(cos(theta_est(j)-thetai)-a)/(1-a);
                if fi(i)<0 fi(i)=0;
                end
                likelyhood_max(i,j)=-(ri(i)-fi(i))*sin(theta_est(j)-thetai)/(noise_sigma^2*(1-a));
            end
        end
        likelyhood_max=sum(likelyhood_max);
        [x,marker]=find(abs(likelyhood_max)==min(abs(likelyhood_max)));
        decoding_theta=theta_est(marker);
        decoding_error(rep,p)=abs(decoding_theta-theta)*180/pi;%caulculate the error
        p=p+1;
    end
end
figure (2);
subplot(2,3,1);
x=(-pi:0.01:pi)*180/pi;
y=mean(decoding_error);
plot(x,y);
title("Maximum likelyhood");
axis([-180,180,0,8]);
xlabel("theta in degrees");
ylabel("error in degrees");

% 1.2Bayesian (easy version...)
clear;
a=-0.14;
noise_sigma=0.1;
for repeat=1:1:100
    j=1;
    for theta=-pi:0.01:pi
        %simulation
        for i=1:4
            thetai=-3*pi/4+pi*(i-1)/2;
            vi(i)=(cos(theta-thetai)-a)/(1-a);
            if vi(i)<0 vi(i)=0;end

            for N=1:1:100
                ri(i,N)=vi(i)+noise_sigma*randn(1);
                if ri(i,N)<0 ri(i,N)=0;end
            end
            %"a easier version of Bayesian"---calculate the average value
            %of V in P(r|V) distribution:
            error_theta(i)=mean(ri(i,:));
%             for delta=1:1:1000
%                 cumulator(delta)=ri(i,delta)*numel(find(ri(i,:)==ri(i,delta)))/1000;%the average of distribution...(this is so awkward, and is not Bayesian at all)
%                 %Prv=normpdf(ri(i,delta),vi(i),noise_sigma);
%                 %cumulator(delta)=1*vi(i)*Prv;
%             end
%             error_theta(i)=sum(cumulator);
            %decode_theta(i)=acos(error_theta(i)*(1-a)+a)+thetai*pi/180;
        end
        %error(repeat,j)= mean(decode_theta-theta);
        error(repeat,j)=acos(dot(error_theta,vi)/(norm(error_theta)*norm(vi)))*180/pi;
        j=j+1;
    end
end
x=(-pi:0.01:pi)*180/pi;
y=sqrt(mean(error));
subplot(2,3,2);
plot(x,y);
axis([-180,180,0,2]);
title("Bayesian");
xlabel("theta in degrees");
ylabel("error in degrees");

% 1.3 least square
clear;
a=-0.14;
noise_sigma=0.1;
for rep=1:1:100
    p=1;
    for theta=-pi:0.01:pi
        %simulation
        for i=1:4
            thetai=-3*pi/4+pi*(i-1)/2;
            vi(i)=(cos(theta-thetai)-a)/(1-a);
            if vi(i)>=0
                vi(i)=vi(i);
            else
                vi(i)=0;
            end
            ri(i)=vi(i)+noise_sigma*randn(1);
            if ri(i)>=0
                ri(i)=ri(i);
            else
                ri(i)=0;
            end
            %according to the square function, the sum can achieve the
            %least value when:
            theta_est=theta-pi/2:0.01:theta+pi/2;
            for j=1:1:length(theta_est)
                fi(i)=(cos(theta_est(j)-thetai)-a)/(1-a);
                if fi(i)<0 fi(i)=0;
                end
                least_square(i,j)=2*(ri(i)-fi(i))*sin(theta_est(j)-thetai)/(noise_sigma*(1-a));
            end
        end
        least_square=sum(least_square);
        [x,marker]=find(abs(least_square)==min(abs(least_square)));
        decoding_theta=theta_est(marker);
        decoding_error(rep,p)=abs(decoding_theta-theta)*180/pi;%calculate the error
        p=p+1;
    end
end
x=(-pi:0.01:pi)*180/pi;
y=mean(decoding_error);
subplot(2,3,3);
plot(x,y);
axis([-180, 180, 0, 8]);
title("least square");
xlabel("theta in degrees");
ylabel("error in degrees");

%1.4vector
clear;
a=-0.14;
optimal=(cos(0)-a)/(1-a);%optimal=1
%establish the C matrix
c(1,:)=[-1,-1]/norm([-1,-1]);
c(2,:)=[1,-1]/norm([1,-1]);
c(3,:)=[1,1]/norm([1,1]);
c(4,:)=[-1,1]/norm([-1,1]);
angle_original=-pi:0.01:pi;
v_input=[cos(angle_original);sin(angle_original)];
noise_sigma=0.1;
for rep=1:1:2000
    for theta=1:1:length(angle_original)
        for i=1:4
            thetai=-3*pi/4+pi*(i-1)/2;
            vi(i)=(dot(v_input(:,theta),c(i,:))-a)/(1-a);
            if vi(i)<0 vi(i)=0;
            end
            %for col=1:1:100
                ri(i)=vi(i)+noise_sigma*randn(1);
                if ri(i)<0 ri(i)=0;
                end
                %popular_decoder_inner(col,:)=ri(i)*c(i,:);
            %end
            %popular_decoder_outsider(i,:)=mean(popular_decoder_inner);
            popular_decoder_outsider(i,:)=ri(i)*c(i,:);
        end
        popular_decoder(theta,:)=sum(popular_decoder_outsider);
        error(rep,theta)=acos(dot(popular_decoder(theta,:)',v_input(:,theta))/(norm(popular_decoder(theta,:))*norm(v_input(:,theta))))*180/pi;
    end
end
subplot(2,3,4);
x=(-pi:0.01:pi)*180/pi;
y=nanmean(error);
plot(x,y);
title("vector");
axis([-180 180 0 8])
xlabel("theta in degrees");
ylabel("error in degrees");

%OLE
clear;
a=-0.14;
noise_sigma=0.1;
angles=-pi:0.01:pi;
v_input=[cos(angles);sin(angles)];
for rep=1:1:2000
for theta=1:1:length(angles)
    for i=1:4
        thetai=-135+90*(i-1);
        vi(i,theta)=(cos(angles(theta)-thetai*pi/180)-a)/(1-a);
        if vi(i,theta)<0 vi(i,theta)=0;
        end
        ri(i,theta)=vi(i,theta)+noise_sigma*randn(1);
        if ri(i,theta)<0 ri(i,theta)=0;
        end
    end
end
for j=1:4
    for theta=1:1:length(angles)
        Lj(:,theta)=v_input(:,theta)*vi(j,theta);
        for q=1:4
                Qij_each_q(q,theta)=vi(q,theta)*vi(j,theta);
        end
    end
    L(j,:)=sum(Lj');
    for q=1:4
        if j~=q Qij(j,q)=sum(Qij_each_q(q,:));
        else Qij(j,q)=noise_sigma^2+sum(Qij_each_q(q,:));
        end
    end
end
Qij=inv(Qij);
for i=1:4
    for j=1:4
        Di_each_i(j,:)=L(j,:)*Qij(i,j);
    end
    Di(i,:)=sum(Di_each_i);
end
for decode_theta=1:1:length(angles)
    for i=1:4
        vest_each_i(i,:)=ri(i,decode_theta)*Di(i,:);
    end
    vest(:,decode_theta)=sum(vest_each_i)';
    error(rep,decode_theta)=acos(dot(vest(:,decode_theta),v_input(:,decode_theta))/(norm(vest(:,decode_theta))*norm(v_input(:,decode_theta))))*180/pi;
end
end
subplot(2,3,6);
x=angles*180/pi;
plot(x,mean(error));
axis([-180,180,0,8]);
title("OLE");
xlabel("theta in degrees");
ylabel("error in degrees");
%% Figure2 
%Figure A
clear;
mark=1;
for a=-0.5:0.01:0.5
    optimal=(cos(0)-a)/(1-a);
    c(1,:)=[-1,-1]/norm([-1,-1]);
    c(2,:)=[1,-1]/norm([1,-1]);
    c(3,:)=[1,1]/norm([1,1]);
    c(4,:)=[-1,1]/norm([-1,1]);
    for rep=1:1:100
        for i=1:4
            thetai=-135+90*(i-1);
            j=1;
            for theta=-pi:0.01:pi
                vi(i,j)=(cos(theta-thetai*pi/180)-a)/(1-a);
                if vi(i,j)>=0
                    vi(i,j)=vi(i,j);
                else
                    vi(i,j)=0;
                end
                noise_sigma(i,j)=0.1;
                ri(i,j)=vi(i,j)+noise_sigma(i,j)*randn(1);%the variability increased proportionaly to the firing rate
                if ri(i,j)>=0
                    ri(i,j)=ri(i,j);
                else
                    ri(i,j)=0;
                end
                j=j+1;
            end
        end
        vi=(-pi:0.01:pi)*180/pi;
        for theta=1:1:length(ri)
            for i=1:4
                vest(i,:)=ri(i,theta)/optimal.*c(i,:);
            end
            error_theta=sum(vest);
            error(rep,theta)=(abs(acos(error_theta(1)/norm(error_theta))*180/pi)-abs(vi(theta)))^2;
        end
    end
    mean_error(mark)=mean(sqrt(mean(error)));
    max_error(mark)=max(sqrt(mean(error)));
    mark=mark+1;
end
a=-0.5:0.01:0.5;
figure(3)
subplot(1,2,1);
plot(a,mean_error);
hold on
plot(a,max_error);
hold off
legend('mean','max');
axis square
title("Panel A");
% Figure 2 B
a=-0.14;
thetai=0;
j=1;
for theta=-pi:0.01:pi
    fi(j)=(cos(theta-0)-a)/(1-a);
    if fi(j)>=0
        fi(j)=fi(j);
    else
        fi(j)=0;
    end
    j=j+1;
end
f_theta_a1=fi;
a=0;
j=1;
for theta=-pi:0.01:pi
    fi(j)=(cos(theta-0)-a)/(1-a);
    if fi(j)>=0
        fi(j)=fi(j);
    else
        fi(j)=0;
    end
    j=j+1;
end
f_theta_a2=fi;
x=(-pi:0.01:pi)*180/pi;
subplot(1,2,2);
plot(x,f_theta_a1);
hold on
plot(x,f_theta_a2);
hold off
legend ('a1=-0.14','a2=0');
axis square
title("Panel B");

%% Figure 3 A
clear;
a=-0.14;
optimal=(cos(0)-a)/(1-a);
angle_original=-pi:0.01:pi;
v_input=[cos(angle_original);sin(angle_original)];
%randomly choose 5 angles
random_num=[20,70,315,375,400];
ci=v_input(:,random_num);
noise_sigma=0.1
for rep=1:1:500
    for theta=1:1:length(angle_original)
        for i=1:5
            vi(i)=dot(v_input(:,theta),ci(:,i));
            if vi(i)<0 vi(i)=0;
            end
            for col=1:1:10
                ri(i)=vi(i)+noise_sigma*randn(1);
                if ri(i)<0 ri(i)=0;
                end
                popular_decoder_inner(col,:)=(ri(i)*ci(:,i))';
            end
            popular_decoder_outsider(i,:)=mean(popular_decoder_inner);
        end
        popular_decoder(theta,:)=sum(popular_decoder_outsider);
        error(rep,theta)=acos(dot(popular_decoder(theta,:)',v_input(:,theta))/(norm(popular_decoder(theta,:))*norm(v_input(:,theta))))*180/pi;
    end
end
x=(-pi:0.01:pi)*180/pi;
y=nanmean(error);
figure (4)
subplot(2,1,1);
plot(x,y);
title("vector");
axis([-180,180,0,50]);
xlabel("theta");
ylabel("error");
axis square

% Figure 3 B
clear;
a=-0.14;
noise_sigma=0.1;
angles=-pi:0.01:pi;
v_input=[cos(angles);sin(angles)];
%randomly choose 5 angles
random_num=[170,189,350,490,590];
ci=v_input(:,random_num);
for rep=1:1:500
for theta=1:1:length(angles)
    for i=1:5
        vi(i,theta)=dot(v_input(:,theta),ci(:,i));
        if vi(i,theta)<0 vi(i,theta)=0;
        end
        ri(i,theta)=vi(i,theta)+noise_sigma*randn(1);
        if ri(i,theta)<0 ri(i,theta)=0;
        end
    end
end
for j=1:5
    for theta=1:1:length(angles)
        Lj(:,theta)=v_input(:,theta)*vi(j,theta);
        for q=1:5
                Qij_each_q(q,theta)=vi(q,theta)*vi(j,theta);
        end
    end
    L(j,:)=sum(Lj');
    for q=1:5
        if j~=q Qij(j,q)=sum(Qij_each_q(q,:));
        else Qij(j,q)=noise_sigma^2+sum(Qij_each_q(q,:));
        end
    end
end
Qij=inv(Qij);
for i=1:5
    for j=1:5
        Di_each_i(j,:)=L(j,:)*Qij(i,j);
    end
    Di(i,:)=sum(Di_each_i);
end
for decode_theta=1:1:length(angles)
    for i=1:5
        vest_each_i(i,:)=ri(i,decode_theta)*Di(i,:);
    end
    vest(:,decode_theta)=sum(vest_each_i)';
    error(rep,decode_theta)=acos(dot(vest(:,decode_theta),v_input(:,decode_theta))/(norm(vest(:,decode_theta))*norm(v_input(:,decode_theta))))*180/pi;
end
end

subplot(2,1,2);
x=angles*180/pi;
plot(x,mean(error));
axis([-180,180,0,50]);
title("OLE");
xlabel("theta");
ylabel("error");
axis square

%% Figure 4: the half cosines
%least squares
clear;
a=-0.14;
angles=-pi:0.1:pi;
noise_sigma=0.1;
mark=1;
for n=1:1:100
    rand_number=randi(numel(angles),1,n);
    thetai=angles(rand_number);
    p=1;
    for theta=-pi:0.1:pi
        for i=1:n
            vi(i)=(cos(theta-thetai(i))-a)/(1-a);
            if vi(i)<0 vi(i)=0;end
            ri(i)=vi(i)+noise_sigma*randn(1);
            if ri(i)<0 ri(i)=0;end
            %according to the square function, the sum can achieve the
            %least value when:
            theta_est=theta-pi/2:0.1:theta+pi/2;
            for j=1:1:length(theta_est)
                fi(i)=(cos(theta_est(j)-thetai(i))-a)/(1-a);
                if fi(i)<0 fi(i)=0;
                end
                least_square(i,j)=2*(ri(i)-fi(i))*sin(theta_est(j)-thetai(i))/(noise_sigma*(1-a));
            end
        end
        least_square=sum(least_square);
        [x,marker]=find(abs(least_square)==min(abs(least_square)));
        if length(marker)>1 marker=marker(1);end
        decoding_theta=theta_est(marker);
        decoding_error(p)=abs(decoding_theta-theta)*180/pi;
        p=p+1;
    end
    error_LS(mark)=mean(decoding_error);
    mark=mark+1;
end
x=1:1:100;
figure(5)
subplot(1,2,2)
loglog(x,error_LS);
hold on


%OLE
clear;
a=-0.14;
noise_sigma=0.1;
angles=-pi:0.01:pi;
v_input=[cos(angles);sin(angles)];
mark=1;
for n=2:10:100
rand_number=randi(numel(angles),1,n);
ci=v_input(:,rand_number);
for theta=1:1:length(angles)
    for i=1:n
        vi(i,theta)=dot(v_input(:,theta),ci(:,i));
        if vi(i,theta)<0 vi(i,theta)=0;
        end
        ri(i,theta)=vi(i,theta)+noise_sigma*randn(1);
        if ri(i,theta)<0 ri(i,theta)=0;
        end
    end
end
for j=1:n
    for theta=1:1:length(angles)
        Lj(:,theta)=v_input(:,theta)*vi(j,theta);
        for q=1:n
                Qij_each_q(q,theta)=vi(q,theta)*vi(j,theta);
        end
    end
    L(j,:)=sum(Lj');
    for q=1:n
        if j~=q Qij(j,q)=sum(Qij_each_q(q,:));
        else Qij(j,q)=noise_sigma^2+sum(Qij_each_q(q,:));
        end
    end
end
Qij=inv(Qij);
for i=1:n
    for j=1:n
        Di_each_i(j,:)=L(j,:)*Qij(i,j);
    end
    Di(i,:)=sum(Di_each_i);
end
for decode_theta=1:1:length(angles)
    for i=1:n
        vest_each_i(i,:)=ri(i,decode_theta)*Di(i,:);
    end
    vest(:,decode_theta)=sum(vest_each_i)';
    error(decode_theta)=acos(dot(vest(:,decode_theta),v_input(:,decode_theta))/(norm(vest(:,decode_theta))*norm(v_input(:,decode_theta))))*180/pi;
end
error_OLE(mark)=mean(error(decode_theta));
mark=mark+1;
end
x=1:10:100;
loglog(x,error_OLE);
hold on
%vector
clear;
a=-0.14;
optimal=1;
angles=-pi:0.1:pi;
v_input=[cos(angles);sin(angles)];
noise_sigma=0.1
mark=1;
for n=2:10:1000
    rand_number=randi(numel(angles),1,n);
    c=v_input(:,rand_number);
    for theta=1:1:length(angles)
        for i=1:n
            thetai=-3*pi/4+pi*(i-1)/2;
            vi(i)=(dot(v_input(:,theta),c(:,i))-a)/(1-a);
            if vi(i)<0 vi(i)=0;
            end
                ri(i)=vi(i)+noise_sigma*randn(1);
                if ri(i)<0 ri(i)=0;
                end
            popular_decoder_outsider(i,:)=ri(i)*c(:,i);
        end
        popular_decoder(theta,:)=sum(popular_decoder_outsider);
        error(theta)=acos(dot(popular_decoder(theta,:)',v_input(:,theta))/(norm(popular_decoder(theta,:))*norm(v_input(:,theta))))*180/pi;
    end
    error_vec(mark)=mean(error);
    mark=mark+1;
end
x=1:10:1000;
loglog(x,error_vec);
legend('LS','OLE','Vector');
xlabel('number of neurons');
ylabel('error');
title("half cosine");
hold off
% Figure 4 : full cosine
%least squares
clear;
a=-0.14;
angles=-pi:0.1:pi;
noise_sigma=0.1;
mark=1;
for n=1:1:80
    rand_number=randi(numel(angles),1,n);
    thetai=angles(rand_number);
    p=1;
    for theta=-pi:0.1:pi
        for i=1:n
            vi(i)=(cos(theta-thetai(i))-a)/(1-a);
            ri(i)=vi(i)+noise_sigma*randn(1);
            %according to the square function, the sum can achieve the
            %least value when:
            theta_est=theta-pi/2:0.1:theta+pi/2;
            for j=1:1:length(theta_est)
                fi(i)=(cos(theta_est(j)-thetai(i))-a)/(1-a);
                least_square(i,j)=2*(ri(i)-fi(i))*sin(theta_est(j)-thetai(i))/(noise_sigma*(1-a));
            end
        end
        least_square=sum(least_square);
        [x,marker]=find(abs(least_square)==min(abs(least_square)));
        if length(marker)>1 marker=marker(1);end
        decoding_theta=theta_est(marker);
        decoding_error(p)=abs(decoding_theta-theta)*180/pi;
        p=p+1;
    end
    error_LS(mark)=mean(decoding_error);
    mark=mark+1;
end
x=1:1:80;
subplot(1,2,1);
loglog(x,error_LS);
axis([0,1000,1,120]);
hold on

%OLE
clear;
a=-0.14;
noise_sigma=0.1;
angles=-pi:0.05:pi;
v_input=[cos(angles);sin(angles)];
mark=1;
for n=2:1:50
rand_number=randi(numel(angles),1,n);
ci=v_input(:,rand_number);
for theta=1:1:length(angles)
    for i=1:n
        vi(i,theta)=dot(v_input(:,theta),ci(:,i));
        ri(i,theta)=vi(i,theta)+noise_sigma*randn(1);
    end
end
for j=1:n
    for theta=1:1:length(angles)
        Lj(:,theta)=v_input(:,theta)*vi(j,theta);
        for q=1:n
                Qij_each_q(q,theta)=vi(q,theta)*vi(j,theta);
        end
    end
    L(j,:)=sum(Lj');
    for q=1:n
        if j~=q Qij(j,q)=sum(Qij_each_q(q,:));
        else Qij(j,q)=noise_sigma^2+sum(Qij_each_q(q,:));
        end
    end
end
Qij=inv(Qij);
for i=1:n
    for j=1:n
        Di_each_i(j,:)=L(j,:)*Qij(i,j);
    end
    Di(i,:)=sum(Di_each_i);
end
for decode_theta=1:1:length(angles)
    for i=1:n
        vest_each_i(i,:)=ri(i,decode_theta)*Di(i,:);
    end
    vest(:,decode_theta)=sum(vest_each_i)';
    error(decode_theta)=acos(dot(vest(:,decode_theta),v_input(:,decode_theta))/(norm(vest(:,decode_theta))*norm(v_input(:,decode_theta))))*180/pi;
end
error_OLE(mark)=mean(error(decode_theta));
mark=mark+1;
end
x=2:1:50;
loglog(x,error_OLE);
axis([0,1000,1,120]);
hold on

%vector
clear;
a=-0.14;
optimal=1;
angles=-pi:0.01:pi;
v_input=[cos(angles);sin(angles)];
noise_sigma=0.1;
mark=1;
for n=2:10:1000
    rand_number=randi(numel(angles),1,n);
    c=v_input(:,rand_number);
    for theta=1:1:length(angles)
        for i=1:n
            thetai=-3*pi/4+pi*(i-1)/2;
            vi(i)=(dot(v_input(:,theta),c(:,i))-a)/(1-a);
                ri(i)=vi(i)+noise_sigma*randn(1);
            popular_decoder_outsider(i,:)=ri(i)*c(:,i);
        end
        popular_decoder(theta,:)=sum(popular_decoder_outsider);
        error(theta)=acos(dot(popular_decoder(theta,:)',v_input(:,theta))/(norm(popular_decoder(theta,:))*norm(v_input(:,theta))))*180/pi;
    end
    error_vec(mark)=mean(error);
    mark=mark+1;
end
x=1:10:1000;
loglog(x,error_vec);
axis([0,1000,1,120]);
legend('LS','OLE','Vector');
xlabel('number of neurons');
ylabel('error');
title("full cosine");
hold off

%% figure 5
%least squares
clear;
a=-0.14;
angles=-pi:0.1:pi;
angles_favourite=[-pi:0.1:0,1:0.1:pi];
noise_sigma=0.1;
mark=1;
for n=1:1:80
    rand_number=randi(numel(angles_favourite),1,n);
    thetai=angles_favourite(rand_number);
    p=1;
    for theta=-pi:0.1:pi
        for i=1:n
            vi(i)=(cos(theta-thetai(i))-a)/(1-a);
            ri(i)=vi(i)+noise_sigma*randn(1);
            %according to the square function, the sum can achieve the
            %least value when:
            theta_est=theta-pi/2:0.1:theta+pi/2;
            for j=1:1:length(theta_est)
                fi(i)=(cos(theta_est(j)-thetai(i))-a)/(1-a);
                least_square(i,j)=2*(ri(i)-fi(i))*sin(theta_est(j)-thetai(i))/(noise_sigma*(1-a));
            end
        end
        least_square=sum(least_square);
        [x,marker]=find(abs(least_square)==min(abs(least_square)));
        if length(marker)>1 marker=marker(1);end
        decoding_theta=theta_est(marker);
        decoding_error(p)=abs(decoding_theta-theta)*180/pi;
        p=p+1;
    end
    error_LS(mark)=mean(decoding_error);
    mark=mark+1;
end
x=1:1:80;
figure(6);
loglog(x,error_LS);
axis([0,1000,1,120]);
hold on

%OLE
clear;
a=-0.14;
noise_sigma=0.1;
angles=-pi:0.01:pi;
angles_favourite=[-pi:0.1:0,1:0.1:pi];
v_input=[cos(angles);sin(angles)];
v_input_favourite=[cos(angles_favourite);sin(angles_favourite)];
mark=1;
for n=2:1:50
rand_number=randi(numel(angles_favourite),1,n);
ci=v_input_favourite(:,rand_number);
for theta=1:1:length(angles)
    for i=1:n
        vi(i,theta)=dot(v_input(:,theta),ci(:,i));
        ri(i,theta)=vi(i,theta)+noise_sigma*randn(1);
    end
end
for j=1:n
    for theta=1:1:length(angles)
        Lj(:,theta)=v_input(:,theta)*vi(j,theta);
        for q=1:n
                Qij_each_q(q,theta)=vi(q,theta)*vi(j,theta);
        end
    end
    L(j,:)=sum(Lj');
    for q=1:n
        if j~=q Qij(j,q)=sum(Qij_each_q(q,:));
        else Qij(j,q)=noise_sigma^2+sum(Qij_each_q(q,:));
        end
    end
end
Qij=inv(Qij);
for i=1:n
    for j=1:n
        Di_each_i(j,:)=L(j,:)*Qij(i,j);
    end
    Di(i,:)=sum(Di_each_i);
end
for decode_theta=1:1:length(angles)
    for i=1:n
        vest_each_i(i,:)=ri(i,decode_theta)*Di(i,:);
    end
    vest(:,decode_theta)=sum(vest_each_i)';
    error(decode_theta)=acos(dot(vest(:,decode_theta),v_input(:,decode_theta))/(norm(vest(:,decode_theta))*norm(v_input(:,decode_theta))))*180/pi;
end
error_OLE(mark)=mean(error(decode_theta));
mark=mark+1;
end
x=2:1:50;
loglog(x,error_OLE);
axis([0,1000,1,120]);
hold on

%vector
clear;
a=-0.14;
optimal=1;
angles=-pi:0.1:pi;
v_input=[cos(angles);sin(angles)];
angles_favourite=[-pi:0.1:0,1:0.1:pi];
v_input_favourite=[cos(angles_favourite);sin(angles_favourite)];
noise_sigma=0.1
mark=1;
for n=2:10:1000
    rand_number=randi(numel(angles_favourite),1,n);
    c=v_input_favourite(:,rand_number);
    for theta=1:1:length(angles)
        for i=1:n
            thetai=-3*pi/4+pi*(i-1)/2;
            vi(i)=(dot(v_input(:,theta),c(:,i))-a)/(1-a);
                ri(i)=vi(i)+noise_sigma*randn(1);
            popular_decoder_outsider(i,:)=ri(i)*c(:,i);
        end
        popular_decoder(theta,:)=sum(popular_decoder_outsider);
        error(theta)=acos(dot(popular_decoder(theta,:)',v_input(:,theta))/(norm(popular_decoder(theta,:))*norm(v_input(:,theta))))*180/pi;
    end
    error_vec(mark)=mean(error);
    mark=mark+1;
end
x=1:10:1000;
loglog(x,error_vec);
axis([0,1000,1,120]);
legend('LS','OLE','Vector');
xlabel('number of neurons');
ylabel('error');
title("biased full cosine");
hold off

%% Appendix：1.1maximum likelyhood (method 1：Not correct for some reason)
%This part takes down some funny mistakes that I made along the way. Run to
%have fun~

% clear;
% a=-0.14;
% for repeat=1:1:10
%     j=1;
%     for theta=-pi:0.01:pi
%         for i=1:4
%             thetai=-135+90*(i-1);
%             vi(i)=(cos(theta-thetai*pi/180)-a)/(1-a);
%             if vi(i)>=0
%                 vi(i)=vi(i);
%             else
%                 vi(i)=0;
%             end
%             noise_sigma=0.1;
%             for N=1:1:1000
%                 ri(i,N)=vi(i)+noise_sigma*randn(1);
%                 if ri(i,N)>=0
%                     ri(i,N)=ri(i,N);
%                 else
%                     ri(i,N)=0;
%                 end
%             end
%             p=mle('norm',ri(i,:));%maximize the p(vest|r)
%             error_theta(i)=p(1);
%             %decode_theta(i)=acos(error_theta(i)*(1-a)+a)+thetai*pi/180;
%         end
%         %error(repeat,j)= mean(decode_theta-theta)*180/pi;
%         error(repeat,j)=acos(dot(error_theta,vi)/(norm(error_theta)*norm(vi)))*180/pi;
%         j=j+1;
%     end
% end
% figure(7);
% subplot(1,2,1);
% x=(-pi:0.01:pi)*180/pi;
% y=mean(error);
% plot(x,y);
% title("ML failure");
% 
% %% 1.2 Bayesian (hard version,I don't know why, but there is no error here...so weird)
% clear
% a=-0.14;
% angle_original=-pi:0.01:pi;
% v_input=[cos(angle_original);sin(angle_original)];
% noise_sigma=0.1;
% for i=1:4
%     thetai=-3*pi/4+pi*(i-1)/2;
%     for theta=1:1:length(angle_original)
%         vi(i,theta)=(cos(theta-thetai)-a)/(1-a);
%         if vi(i,theta)<0 vi(i,theta)=0;
%         end
%         ri(i,theta)=vi(i,theta)+noise_sigma*randn(1);
%         if ri(i,theta)<0 ri(i,theta)=0;
%         end
%     end
% end
% for i=1:4
%     for theta=1:1:length(angle_original)
%         additional_angle=theta-pi/4:0.01:theta+pi/4;
%         for decode_theta=1:1:length(additional_angle)
%             p_r_v(decode_theta)=normcdf(ri(i,theta),vi(i,decode_theta),noise_sigma);
%             if p_r_v(decode_theta)>0.5 p_r_v(decode_theta)=1-p_r_v(decode_theta);
%             end
%         decode_each_theta_inner(decode_theta)=angle_original(theta)*p_r_v(decode_theta);
%     end
%     decode_each_theta(i,theta)=sum(decode_each_theta_inner)/sum(p_r_v);
%     end
% end
%     decode_vest=mean(decode_each_theta);
%     for theta=1:1:length(angle_original)
%         error(theta)=abs(decode_vest(theta)-angle_original(theta))*180/pi;
%     end
%     subplot(1,2,2);
%     x=angle_original*180/pi;
%     plot(x,error);
%     title("strange Bayesian");