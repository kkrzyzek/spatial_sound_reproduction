%Spatial Sound Reproduction for 16 speakers surround system
%Auralization Lab

clear all
close all
clear sound
clc

%polozenie glosnikow we wsp sferycznych
fi = [0 45 90 135 180 -135 -90 -45 68 158 -113 -23 23 113 -158 -68];
teta = [0 0 0 0 0 0 0 0 -45 -45 -45 -45 45 45 45 45];

%tor1 - lokalizacja dzwieku we wsp sferycznych
%FI=[0 10 20 80 90 100 0 20 70 80 90 100]; %(-180,180)
%TETA=[0 0 20 30 40 50 0 -20 20 30 40 50]; %(-90,90)

%tor2 zrodla pozornego - spirala po sferze
c = 5.0;
t = 0:0.2:10*pi;

x1 = sin(t/(2*c)).*cos(t);
y1 = sin(t/(2*c)).*sin(t);
z1 = cos(t/(2*c));

L123 = [x1',y1',z1'];

figure(1)
plot3(x1,y1,z1) 
grid on
title('Tor dzwieku')
print('tor_dzwieku','-dpng')

%%%%%%%%%%%%%%%%%%%%%%%
%TRIANGULARYZACJA TORU%

INDEX3 = [];
INDEX = [];
TEMP = [];

%polozenia glosnikow we wsp kartezjanskich
x = cos(teta).*cos(fi);
y = cos(teta).*sin(fi);
z = sin(teta);

%polozenie pozornego zrodla we wsp kartezjanskich dla toru 1
%x1=cos(TETA).*cos(FI);
%y1=cos(TETA).*sin(FI);
%z1=sin(TETA);

for i = 1:length(t)
    
    for k = 1:length(fi)
        
        temp = sqrt((x1(i)-x(k))^2+(y1(i)-y(k))^2+(z1(i)-z(k))^2);
        TEMP(k,i) = temp;
        
    end
    
    [tmp,idx] = sort(TEMP(:,i));
    INDEX3 = [idx(1); idx(2); idx(3)];
    INDEX = [INDEX INDEX3];
    
end

%wczytanie nagrania MONO
[s, fs] = audioread('helicopter_mono.wav');

%skrocenie sygnalu
u = length(s)-mod(length(s), length(t));
s = s(1:u);

%tworzenie macierzy dzwiekow dla kazdego z glosnikow + korekcja wzmocnienia
S = [];
m = length(s)/length(t);

for p=1:length(x1)

    g1 = TEMP(INDEX(1,p),p);
    g2 = TEMP(INDEX(2,p),p);
    g3 = TEMP(INDEX(3,p),p);
    
    g1n = g1/(sqrt(g1^2+g2^2+g3^2));
    g2n = g2/(sqrt(g1^2+g2^2+g3^2));
    g3n = g3/(sqrt(g1^2+g2^2+g3^2));
    
    S(INDEX(1,p),(p-1)*m+1:m*((p-1)+1)) = s((p-1)*m+1:m*((p-1)+1))*g1n;
    S(INDEX(2,p),(p-1)*m+1:m*((p-1)+1)) = s((p-1)*m+1:m*((p-1)+1))*g2n;
    S(INDEX(3,p),(p-1)*m+1:m*((p-1)+1)) = s((p-1)*m+1:m*((p-1)+1))*g3n;

end

%zapis 16 nagran, dla kazdego z glosnikow
for f=1:length(fi)
    
    filename = sprintf('speaker%d.wav', f); 
    audiowrite(filename, S(f,:), fs);

end