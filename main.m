clc;clearvars;
% aynı scriptte function calismasi icin gerekli.
mainfun = localfunctions;

% format: BN SN dH Si
OLCULER = [
    4 1 43.156 0.65
    2 1 19.218 0.80
    2 3 33.524 1.00
    4 3 57.440 1.40
    4 2 23.962 1.5 % 23.962
    1 3 14.267 1.95];

% format: noktaNo yaklasikYukseklik
YAKLASIK = [
    1 123.829
    2 104.635
    3 138.113
    4 80.673];

d = 1; % nivelman aglari icin datum parametresi
% f = n(olcu sayisi) - u(bilinmeyen sayisi) + d(datum parametresi)
% f = size(OLCULER, 1) - size(YAKLASIK, 1) + d; 

% agirlik matrislerinin olusturulmasi
P = diag(1 ./ OLCULER(:, 4));

% V = Ax - l

% Katsayilar matrisi ve otelenmis olculer vektorunu elinizle bulabilirsiniz veya bunlari
% cozen algoritma yazabilirsiniz.
% katsayilar matrisi
A = [
    1 0 0 -1
    1 -1 0 0
    0 -1 1 0
    0 0 1 -1
    0 1 0 -1
    -1 0 1 0];

% l otelenmis olcu vektoru
l = [
    OLCULER(1, 3) - (YAKLASIK(1, 2) - YAKLASIK(4, 2))
    OLCULER(2, 3) - (YAKLASIK(1, 2) - YAKLASIK(2, 2))
    OLCULER(3, 3) - (YAKLASIK(3, 2) - YAKLASIK(2, 2))
    OLCULER(4, 3) - (YAKLASIK(3, 2) - YAKLASIK(4, 2))
    OLCULER(5, 3) - (YAKLASIK(2, 2) - YAKLASIK(4, 2))
    OLCULER(6, 3) - (YAKLASIK(3, 2) - YAKLASIK(1, 2))
    ]; % m

% l vektorunu mm yapıyoruz
l = round(l * 1e3) ;%mm

% eger uyusumsuz olcu varsa dongude bunu bize vericektir.
while true
    % olculer cikarildigi icin tekrar serbestlik derecesi hesaplanir.
    f = size(A, 1) - size(A, 2) + d;
    % N, n -> Qxx
    % d = 1
    % [
    %    1/sqrt(p)
    %    1/sqrt(p)
    %](uxd)
    % x = (A' * P * A)^+ * A' * P * l;

    %normal denklemlerin katsayilar matrisi
    N = A' * P * A ;% N^+ = Qxx
    %sabit terimler vektoru
    n = A' * P * l; 

    %dengeleme bilinmeyeni
    x = round(pinv(N)*n, 2);
    % bilinmeyenlerin kesin degeri
    H = round(YAKLASIK(:, 2) + x * 0.001, 4);

    %duzeltmeler.
    V = A * x - l;

    % --duyarlılık hesaplari--

    % birim olcunu ortalama hatasi
    m0 = sqrt((V' * P * V) / f);
    % olculerin ortalama hatasi 
    mi = m0 ./ sqrt(diag(P));
    %bilinmeyenlerin ortalama hatasi
    mx = m0 .* sqrt(diag(pinv(N)));
    %dengeli olculerin ters agirlik matrisi
    Qll = A * pinv(N) * A';
    %dengeli olculerin ortalama hatası
    ml = m0 * sqrt(diag(Qll)); 
    % duzeltmelerin ters agirlik matrisi
    Qvv = P^-1 - Qll;
    % duzeltmelerin ortalama hatasi
    mv =  m0 * sqrt(diag(Qvv));
    
    % student testi.
    [bl, index] = studentTest(V, P, Qvv, f); 
    if bl
        % eger uyusumsuz olcu yoksa dongu sonlanir.
        break
    else
        % eger uyusumsuz olcu varsa dongu devam eder.
        % bu matrislerdeki uyusumsuz olcuye denk gelen olculeri siliyoruz.
        V(index) = [];
        P = diag(P); P(index) = [];
        P = diag(P);
        l(index) = [];
        A(index, :) = [];
    end 
end



% @Function, uyusumsuz olculer testi.
% girdiler: düzeltme denklemleri(V), agirlik matrisi(P_), 
% Qvv(duzeltmelerin ters agirlik matrisi), f (serbestlik derecesi)
% 
% ciktilar: bl = boolean: uyusumsuz olcu yoksa true, varsa false doner
% index = hangi olcu uyusumsuz ise onun indexini verir.
function [bl, index] = studentTest(V, P_, Qvv, f)
s0 = sqrt((1 / (f - 1)) * ((V'*P_*V) - (V.^2) ./ (diag(Qvv))));
T = abs(V) ./ (s0.*sqrt(diag(Qvv)));
T_ = max(T);
bl = true;
index = 0;
% yuzde 5 lik two tailed.
alfa = 5 ;
% two tailed oldugundan -> (alfa / 100) / 2
t = tinv(1 - (alfa / 100)/2, f - 1);
    if T_ > t
        index = find(T_ == T); 
        fprintf('%d. olcu uyusumsuzdur\n', index);
        bl = false;
    else
        fprintf('uyusumsuz olcu yoktur.\n')
    end
end



