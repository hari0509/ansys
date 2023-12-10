%% Calculation of Mu & Epsilon from Data file
% Development Date : 18.04.2020
format short;

%% Insertion of Data files.
S11_amp = dlmread('/MATLAB Drive/sample_ansys/S Parameter Plot 11_mag.csv');
S11_phase = dlmread('/MATLAB Drive/sample_ansys/S Parameter Plot 11_phase.csv');
S21_amp = dlmread('/MATLAB Drive/sample_ansys/S Parameter Plot 21_mag.csv');
S21_phase = dlmread('/MATLAB Drive/sample_ansys/S Parameter Plot 21_phase.csv');
C = 3e8;                        % Velocity of light in free space in m/s
Freq = S11_amp(:,1);            % Frequency zone of Operation
lamda0 = C./Freq;                % Free Space Wavelength
lamdaC1 = .0457;
lamdaC2 = .0316;

L = 0.025;                       % Sample length in meter

% Starting of the calculation...
% Stage - 1 : Calculation of S11 & S21
S11 = S11_amp(:,2).*(cos(S11_phase(:,2)*pi/180)+1i*sin(S11_phase(:,2)*pi/180));
S21 = S21_amp(:,2).*(cos(S21_phase(:,2)*pi/180)+1i*sin(S21_phase(:,2)*pi/180));

% Stage - 2 : Calculation of various relevant parameters
X = ((S11.^2)-(S21.^2)+1.0)./(2.*S11);
Y1 = X + sqrt((X.^2)-1.0);
Y2 = X - sqrt((X.^2)-1.0);
for ii = 1:length(Y1)
    if abs(Y1(ii)) >= 1
        Y(ii) = Y2(ii);
%         disp(ii)
    else
        Y(ii) = Y1(ii);
    end
end
T = (S11+S21-Y')./(1.0-((S11+S21).*Y'));
Q = -((log(1./T))./((2*pi*L)).^2);
R = sqrt(Q);
%E = sqrt((1./(lamda0).^2)-(1/(lamdaC).^2));
%a = dlmread('A.txt');
%%b = dlmread('B.txt');
%q = a+1i.*b;
Result = [Q, R];
Filename = ['Cu_param dated - ', num2str(date), '.dat'];
fid1 = fopen(Filename,'w');
% Writing of the Variables at the end of the programme
            data1 = [(Freq'/1e9); real(Q)'; imag(Q)'; real(R)';imag(R)'];     
            fprintf(fid1,'\n%12.8f %12.8f %12.8f %12.8f %12.8f',data1);
            fclose(fid1);
            a = dlmread('/MATLAB Drive/sample_ansys/gamma_re.csv');
            b = dlmread('/MATLAB Drive/sample_ansys/gamma_im.csv');
            q = a+1i.*b;

% Stage - 3 : Final calculation of relative permittivity & permiability
for ii = 1:length(Freq)
    ff = Freq(ii);
    if ff <= 12e9
        lamdaC = lamdaC1;
    else
        lamdaC = lamdaC2;
    end
    I(ii) = sqrt((1./lamda0(ii).^2)-(1./lamdaC.^2));

ee = I';

Mu_eff = (q.*((1+Y')./(1-Y')))./(I(ii));
Eps_eff = (((Q+(1./(lamdaC).^2))).*((lamda0(ii)).^2))./(Mu_eff);
end
% Storing the generated values in a Data file
Filename = ['Cu_Relative Permittivity & Permiability file, dated-', num2str(date), '_Starting freq.-',num2str(Freq(1)/1e9),...
            ' GHz_End freq.-',num2str(Freq(end)/1e9), ' GHz_Sample length-',num2str(L),' m.dat']
% Opening of the File:
            fid1 = fopen(Filename,'w');
% $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $
            % Writing Variable Names:
            dataA1 = ['Frequency', ' Real_S11', ' Imag_S11', ' Real_S21', ' Imag_S21', ' Real_Energy', ' Imag_Energy',...
                       ' Real_Effec._Mu', ' Imag_Effec._Mu', ' Real_Effec.Eps', ' Imag_Effec.Eps'];      
            fprintf(fid1,'%s %s %s %s %s %s %s %s %s%s %s\n',dataA1);
% $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $
% Writing Units corresponding to the Variables:
            dataA2 = [' GHz', ' -', ' -', ' -', ' -', ' Joule', ' -', ' -', ' -', ' -', ' -']';      
            fprintf(fid1,'\n%s %s %s %s %s %s %s %s %s %s %s',dataA2);
% $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $ $
% Writing of the Variables at the end of the programme
            data1 = [(Freq'/1e9); real(S11)'; imag(S11)'; real(S21)'; imag(S21)'; real(ee)'; imag(ee)';...
                       real(Mu_eff)'; imag(Mu_eff)'; real(Eps_eff)'; imag(Eps_eff)'];      
            fprintf(fid1,'\n%12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f',data1);
            fclose(fid1);

figure(1)
subplot(2,1,1), plot(Freq,real(Mu_eff),'r')
title('Effective Permiability...');
% xlabel('Frequency');
ylabel('Real Mu');
subplot(2,1,2), plot(Freq,imag(Mu_eff),'b')
% title('Effective Permiability...');
xlabel('Frequency');
ylabel('Imag. Mu');

figure(2)
subplot(2,1,1), plot(Freq,real(Eps_eff),'r')
title('Effective Permittivity...');
% xlabel('Frequency');
ylabel('Real Eps');
subplot(2,1,2), plot(Freq,imag(Eps_eff),'b')
% title('Effective Permittivity...');
xlabel('Frequency');
ylabel('Imag. Eps');
